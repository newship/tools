// Copyright 2015 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// This file contains the Go wrapper for the constant-time, 64-bit assembly
// implementation of P256. The optimizations performed here are described in
// detail in:
// S.Gueron and V.Krasnov, "Fast prime field elliptic-curve cryptography with
//                          256-bit primes"
// https://link.springer.com/article/10.1007%2Fs13389-014-0090-x
// https://eprint.iacr.org/2013/816.pdf

// +build amd64,!field arm64,!field

package sm2

import (
	"errors"
	"io"
	"math/big"
	"sync"
)

type (
	p256Curve struct {
		*CurveParams
		pMinusN *big.Int
	}

	p256Point struct {
		xyz [12]uint64
	}
)

var (
	pSM2            = p256Curve{CurveParams: sm2Params}
	p256Precomputed *[43][32 * 8]uint64
	precomputeOnce  sync.Once
	n2minus         *big.Int
	bigOne          = new(big.Int).SetInt64(1)
	two             = new(big.Int).SetInt64(2)
	one             = []uint64{1, 0, 0, 0}
	p256MontOneA    = [4]uint64{1, 0xffffffff, 0, 0x100000000}
	p256MontOne     = p256MontOneA[:]
)
var errZeroParam = errors.New("zero parameter")
var errParam = errors.New("error parameter")
var errSqrt = errors.New("error sqrt")

// nRR 2^512 mod N
var nRR = []uint64{0x901192af7c114f20, 0x3464504ade6fa2fa, 0x620fc84c3affe0d4,
	0x1eb5e412a22b3d3b}

// p256Mul operates in a Montgomery domain with R = 2^256 mod p, where p is the
// underlying field of the curve. (See initP256 for the value.) Thus rr here is
// R×R mod p. See comment in Inverse about how this is used.
var rr = []uint64{0x0000000200000003, 0x00000002ffffffff, 0x0000000100000001,
	0x0000000400000002}

func initSM2() {
	// See FIPS 186-3, section D.2.3
	//pSM2.CurveParams = sm2Params
	n2minus = new(big.Int).Sub(sm2Params.N, two)
	pSM2.pMinusN = new(big.Int).Sub(pSM2.P, pSM2.N)
}

func (curve p256Curve) Params() *CurveParams {
	return curve.CurveParams
}

// Functions implemented in p256_asm_*64.s
// Add modulo sm2 Prime256
//go:noescape
func p256Add(res, in1, in2 []uint64)

// Sub modulo sm2 Prime256
//go:noescape
func p256Sub(res, in1, in2 []uint64)

// Functions implemented in p256_asm_*64.s
// Montgomery multiplication modulo P256
//go:noescape
func p256Mul(res, in1, in2 []uint64)

// Montgomery square modulo P256, repeated n times (n >= 1)
//go:noescape
func p256Sqr(res, in []uint64, n int)

// Montgomery multiplication by 1
//go:noescape
func p256FromMont(res, in []uint64)

// iff cond == 1  val <- -val
//go:noescape
func p256NegCond(val []uint64, cond int)

// if cond == 0 res <- b; else res <- a
//go:noescape
func p256MovCond(res, a, b []uint64, cond int)

// Endianness swap
//go:noescape
func p256BigToLittle(res []uint64, in []byte)

//go:noescape
func p256LittleToBig(res []byte, in []uint64)

// Constant time table access
//go:noescape
func p256Select(point, table []uint64, idx int)

//go:noescape
func p256SelectBase(point, table []uint64, idx int)

// Montgomery multiplication modulo Ord(G)
//go:noescape
func p256OrdMul(res, in1, in2 []uint64)

// Montgomery square modulo Ord(G), repeated n times
//go:noescape
func p256OrdSqr(res, in []uint64, n int)

// Point add with in2 being affine point
// If sign == 1 -> in2 = -in2
// If sel == 0 -> res = in1
// if zero == 0 -> res = in2
//go:noescape
func p256PointAddAffineAsm(res, in1, in2 []uint64, sign, sel, zero int)

// Point add. Returns one if the two input points were equal and zero
// otherwise. (Note that, due to the way that the equations work out, some
// representations of ∞ are considered equal to everything by this function.)
//go:noescape
func p256PointAddAsm(res, in1, in2 []uint64) int

// Point double
//go:noescape
func p256PointDoubleAsm(res, in []uint64)

func (curve p256Curve) Inverse(k *big.Int) *big.Int {
	if k.Sign() < 0 {
		// This should never happen.
		k = new(big.Int).Neg(k)
	}

	if k.Cmp(curve.N) >= 0 {
		// This should never happen.
		k = new(big.Int).Mod(k, curve.N)
	}
	return new(big.Int).ModInverse(k, curve.N)
}

// fromBig converts a *big.Int into a format used by this code.
func fromBig(out []uint64, big *big.Int) {
	for i := range out {
		out[i] = 0
	}
	outLen := len(out)

	for i, v := range big.Bits() {
		if i >= outLen {
			break
		}
		out[i] = uint64(v)
	}
}

// toBig convert u64 slice format to *big.Int
func toBig(x []uint64) *big.Int {
	var res [4]big.Word
	res[0] = big.Word(x[0])
	res[1] = big.Word(x[1])
	res[2] = big.Word(x[2])
	res[3] = big.Word(x[3])
	return new(big.Int).SetBits(res[:])
}

// p256GetScalar endian-swaps the big-endian scalar value from in and writes it
// to out. If the scalar is equal or greater than the order of the group, it's
// reduced modulo that order.
func p256GetScalar(out []uint64, in []byte) {
	if len(in) > 32 {
		in = in[:32]
	}
	n := new(big.Int).SetBytes(in)

	if n.Cmp(pSM2.N) >= 0 {
		n.Sub(n, pSM2.N)
	}
	fromBig(out, n)
}

func maybeReduceModP(in *big.Int) *big.Int {
	if in.Cmp(pSM2.P) < 0 {
		return in
	}
	//return new(big.Int).Mod(in, pSM2.P)
	return in.SetBits(in.Bits()[:4])
}

func (curve p256Curve) combinedMult(bigX, bigY *big.Int, baseScalar, scalar []byte) *p256Point {
	scalarReversed := make([]uint64, 4)
	var r1, r2 p256Point
	p256GetScalar(scalarReversed, baseScalar)
	r1IsInfinity := scalarIsZero(scalarReversed)
	r1.p256BaseMult(scalarReversed)

	p256GetScalar(scalarReversed, scalar)
	r2IsInfinity := scalarIsZero(scalarReversed)
	fromBig(r2.xyz[0:4], maybeReduceModP(bigX))
	fromBig(r2.xyz[4:8], maybeReduceModP(bigY))
	p256Mul(r2.xyz[0:4], r2.xyz[0:4], rr[:])
	p256Mul(r2.xyz[4:8], r2.xyz[4:8], rr[:])

	// This sets r2's Z value to 1, in the Montgomery domain.
	copy(r2.xyz[8:], p256MontOne)

	r2.p256ScalarMult(scalarReversed)

	var sum, double p256Point
	pointsEqual := p256PointAddAsm(sum.xyz[:], r1.xyz[:], r2.xyz[:])
	p256PointDoubleAsm(double.xyz[:], r1.xyz[:])
	sum.CopyConditional(&double, pointsEqual)
	sum.CopyConditional(&r1, r2IsInfinity)
	sum.CopyConditional(&r2, r1IsInfinity)

	return &sum
}

func (curve p256Curve) CombinedMult(bigX, bigY *big.Int, baseScalar, scalar []byte) (x, y *big.Int) {
	sum := curve.combinedMult(bigX, bigY, baseScalar, scalar)
	return sum.p256PointToAffine()
}

func (curve p256Curve) ScalarBaseMult(scalar []byte) (x, y *big.Int) {
	scalarReversed := make([]uint64, 4)
	p256GetScalar(scalarReversed, scalar)

	var r p256Point
	r.p256BaseMult(scalarReversed)
	return r.p256PointToAffine()
}

func (curve p256Curve) ScalarMult(bigX, bigY *big.Int, scalar []byte) (x, y *big.Int) {
	scalarReversed := make([]uint64, 4)
	p256GetScalar(scalarReversed, scalar)

	var r p256Point
	fromBig(r.xyz[0:4], maybeReduceModP(bigX))
	fromBig(r.xyz[4:8], maybeReduceModP(bigY))
	p256Mul(r.xyz[0:4], r.xyz[0:4], rr[:])
	p256Mul(r.xyz[4:8], r.xyz[4:8], rr[:])
	// This sets r2's Z value to 1, in the Montgomery domain.
	copy(r.xyz[8:], p256MontOne)

	r.p256ScalarMult(scalarReversed)
	return r.p256PointToAffine()
}

func (c p256Curve) Verify(r, s, msg, px, py *big.Int) bool {
	N := c.N
	if r.Sign() <= 0 || s.Sign() <= 0 {
		return false
	}
	if r.Cmp(N) >= 0 || s.Cmp(N) >= 0 {
		return false
	}
	t := new(big.Int).Add(r, s)
	if t.Cmp(N) >= 0 {
		t.Sub(t, N)
	}
	if t.Sign() == 0 {
		return false
	}
	pt := c.combinedMult(px, py, s.Bytes(), t.Bytes())
	x1 := new(big.Int).Sub(r, msg)
	if x1.Sign() < 0 {
		x1.Add(x1, N)
	}
	// zz = pt.z ^2
	var zz [4]uint64
	p256Sqr(zz[:], pt.xyz[8:], 1)
	if p256ProdEqual(pt.xyz[:4], zz[:], x1) {
		return true
	}
	if x1.Cmp(c.pMinusN) < 0 {
		x1.Add(x1, N)
		if p256ProdEqual(pt.xyz[:4], zz[:], x1) {
			return true
		}
	}
	return false
}

func (c p256Curve) Sign(rand io.Reader, msg, secret *big.Int) (r, s *big.Int,
	v uint, err error) {
	var kB [32]byte
	N := c.N
	dInv := new(big.Int).Add(secret, bigOne)
	if dInv.ModInverse(dInv, N) == nil {
		return nil, nil, 0, errZeroParam
	}
	var res, xp, yp [4]uint64
	fromBig(xp[:], dInv)
	p256OrdMul(xp[:], xp[:], nRR)
	for {
		rand.Read(kB[:])
		k := new(big.Int).SetBytes(kB[:])
		if k.Cmp(n2minus) >= 0 {
			k.Sub(k, n2minus)
		}
		k.Add(k, bigOne)
		x1, y1 := c.ScalarBaseMult(k.Bytes())
		r = new(big.Int).Add(x1, msg)
		if r.Cmp(N) >= 0 {
			r.Sub(r, N)
		}
		if r.Sign() == 0 {
			continue
		}
		v = y1.Bit(0)
		k.Add(k, r)
		if k.Cmp(N) >= 0 {
			k.Sub(k, N)
		}
		if k.Sign() == 0 {
			continue
		}
		fromBig(yp[:], k)
		p256OrdMul(yp[:], yp[:], nRR)
		p256OrdMul(res[:], xp[:], yp[:])
		p256OrdMul(res[:], res[:], one)
		s = toBig(res[:])
		if s.Sign() == 0 {
			continue
		}
		s.Sub(s, r)
		if s.Sign() < 0 {
			s.Add(s, N)
		}
		if s.Sign() != 0 {
			break
		}
	}
	return
}

func RecoverPoint(x1 *big.Int, v uint) (y1 *big.Int, err error) {
	var xp, t1 [4]uint64
	c := pSM2
	if x1.Sign() <= 0 || x1.Cmp(c.N) >= 0 {
		return nil, errParam
	}
	fromBig(xp[:], x1)
	p256Mul(xp[:], xp[:], rr)
	p256Sqr(t1[:], xp[:], 1)
	p256Mul(t1[:], t1[:], xp[:])
	// t1 = x1^3
	p256FromMont(t1[:], t1[:])
	fromBig(xp[:], x1)
	p256Sub(t1[:], t1[:], xp[:])
	p256Sub(t1[:], t1[:], xp[:])
	p256Sub(t1[:], t1[:], xp[:])
	fromBig(xp[:], c.B)
	p256Add(t1[:], t1[:], xp[:])
	copy(xp[:], t1[:])
	if !p256Sqrt(t1[:], xp[:]) {
		return nil, errSqrt
	}
	tt := toBig(t1[:])
	if (v ^ tt.Bit(0)) != 0 {
		tt.Sub(c.P, tt)
	}
	return tt, nil
}

func (c p256Curve) Recover(r, s, msg *big.Int, v uint) (pubX, pubY *big.Int, err error) {
	if r.Sign() <= 0 || s.Sign() <= 0 {
		return nil, nil, errParam
	}
	if r.Cmp(c.N) >= 0 || s.Cmp(c.N) >= 0 {
		return nil, nil, errParam
	}
	u1 := new(big.Int).Add(r, s)
	if u1.Cmp(c.N) >= 0 {
		u1.Sub(u1, c.N)
	}
	x1 := new(big.Int).Sub(r, msg)
	if u1.Sign() == 0 || x1.Sign() == 0 {
		return nil, nil, errParam
	} else if x1.Sign() < 0 {
		x1.Add(x1, c.N)
	}
	var y1 *big.Int
	if y1, err = RecoverPoint(x1, v); err != nil {
		return
	}
	if u1.ModInverse(u1, c.N) == nil {
		return nil, nil, errParam
	}
	var xp, yp, res [4]uint64
	fromBig(xp[:], u1)
	fromBig(yp[:], s)
	p256OrdMul(xp[:], xp[:], nRR) // to montgomery form
	p256OrdMul(yp[:], yp[:], nRR)
	p256OrdMul(res[:], xp[:], yp[:])
	p256OrdMul(res[:], res[:], one) // convert from montgomery
	u2 := toBig(res[:])
	u2.Sub(c.N, u2) // u2 = -u2 = - u1*s
	pubX, pubY = c.CombinedMult(x1, y1, u2.Bytes(), u1.Bytes())
	return
}

// uint64IsZero returns 1 if x is zero and zero otherwise.
func uint64IsZero(x uint64) int {
	x = ^x
	x &= x >> 32
	x &= x >> 16
	x &= x >> 8
	x &= x >> 4
	x &= x >> 2
	x &= x >> 1
	return int(x & 1)
}

// scalarIsZero returns 1 if scalar represents the zero value, and zero
// otherwise.
func scalarIsZero(scalar []uint64) int {
	return uint64IsZero(scalar[0] | scalar[1] | scalar[2] | scalar[3])
}

func (p *p256Point) p256PointToAffine() (x, y *big.Int) {
	zInv := make([]uint64, 4)
	zInvSq := make([]uint64, 4)
	p256Inverse(zInv, p.xyz[8:12])
	p256Sqr(zInvSq, zInv, 1)
	p256Mul(zInv, zInv, zInvSq)

	p256Mul(zInvSq, p.xyz[0:4], zInvSq)
	p256Mul(zInv, p.xyz[4:8], zInv)

	p256FromMont(zInvSq, zInvSq)
	p256FromMont(zInv, zInv)

	return toBig(zInvSq), toBig(zInv)
}

// return prod == xp * yp Mod prime p
//	yp = y convert to montgomery form
func p256ProdEqual(prod, xp []uint64, y *big.Int) bool {
	var yp [4]uint64
	fromBig(yp[:], y)
	p256Mul(yp[:], yp[:], rr)
	p256Mul(yp[:], xp, yp[:])
	return prod[0] == yp[0] && prod[1] == yp[1] && prod[2] == yp[2] &&
		prod[3] == yp[3]
}

// CopyConditional copies overwrites p with src if v == 1, and leaves p
// unchanged if v == 0.
func (p *p256Point) CopyConditional(src *p256Point, v int) {
	pMask := uint64(v) - 1
	srcMask := ^pMask

	for i, n := range p.xyz {
		p.xyz[i] = (n & pMask) | (src.xyz[i] & srcMask)
	}
}

// p256Inverse sets out to in^-1 mod p.
func p256Inverse(out, in []uint64) {
	var stack [5 * 4]uint64
	p2 := stack[4*0 : 4*0+4]
	p4 := stack[4*1 : 4*1+4]
	p8 := stack[4*2 : 4*2+4]
	p16 := stack[4*3 : 4*3+4]
	p32 := stack[4*4 : 4*4+4]
	//p3 := stack[4*5 : 4*5+4]

	p256Sqr(out, in, 1)
	p256Mul(p2, out, in) // 3*p

	//p256Sqr(out, p2, 1)
	//p256Mul(p3, out, in) // 7*p

	p256Sqr(out, p2, 2)
	p256Mul(p4, out, p2) // f*p

	p256Sqr(out, p4, 4)
	p256Mul(p8, out, p4) // ff*p

	p256Sqr(out, p8, 8)
	p256Mul(p16, out, p8) // ffff*p

	p256Sqr(out, p16, 16)
	p256Mul(p32, out, p16) // ffffffff*p

	p256Sqr(out, p16, 8)
	p256Mul(out, out, p8) // ffffff*p
	p256Sqr(out, out, 4)
	p256Mul(out, out, p4) // fffffff*p
	p256Sqr(out, out, 2)
	p256Mul(out, out, p2)
	p256Sqr(out, out, 1)
	p256Mul(out, out, in)
	p256Sqr(out, out, 1) // fffffffe*p

	p256Sqr(out, out, 32)
	p256Mul(out, out, p32) // fffffffeffffffff

	p256Sqr(out, out, 32)
	p256Mul(out, out, p32)
	p256Sqr(out, out, 32)
	p256Mul(out, out, p32) // ... ffffffffffffffff
	p256Sqr(out, out, 32)
	p256Mul(out, out, p32)
	p256Sqr(out, out, 32) // ... ffffffff00000000

	p256Sqr(out, out, 32)
	p256Mul(out, out, p32)

	p256Sqr(out, out, 16)
	p256Mul(out, out, p16)

	p256Sqr(out, out, 8)
	p256Mul(out, out, p8)

	p256Sqr(out, out, 4)
	p256Mul(out, out, p4)

	p256Sqr(out, out, 2)
	p256Mul(out, out, p2)

	p256Sqr(out, out, 2)
	p256Mul(out, out, in)
}

// p256ExpQuadP sets out to in^quadP mod p.
func p256ExpQuadP(out, in []uint64) {
	var stack [5 * 4]uint64
	p2 := stack[4*0 : 4*0+4]
	p4 := stack[4*1 : 4*1+4]
	p8 := stack[4*2 : 4*2+4]
	p16 := stack[4*3 : 4*3+4]
	p32 := stack[4*4 : 4*4+4]
	//p3 := stack[4*5 : 4*5+4]

	p256Sqr(out, in, 1)
	p256Mul(p2, out, in) // 3*p

	//p256Sqr(out, p2, 1)
	//p256Mul(p3, out, in) // 7*p

	p256Sqr(out, p2, 2)
	p256Mul(p4, out, p2) // f*p

	p256Sqr(out, p4, 4)
	p256Mul(p8, out, p4) // ff*p

	p256Sqr(out, p8, 8)
	p256Mul(p16, out, p8) // ffff*p

	p256Sqr(out, p16, 16)
	p256Mul(p32, out, p16) // ffffffff*p

	p256Sqr(out, p16, 8)
	p256Mul(out, out, p8) // ffffff*p
	p256Sqr(out, out, 4)
	p256Mul(out, out, p4) // fffffff*p
	p256Sqr(out, out, 2)
	p256Mul(out, out, p2)
	p256Sqr(out, out, 1)
	p256Mul(out, out, in)
	p256Sqr(out, out, 1) // fffffffe*p

	p256Sqr(out, out, 32)
	p256Mul(out, out, p32) // fffffffeffffffff

	p256Sqr(out, out, 32)
	p256Mul(out, out, p32)
	p256Sqr(out, out, 32)
	p256Mul(out, out, p32) // ... ffffffffffffffff
	p256Sqr(out, out, 32)
	p256Mul(out, out, p32)
	p256Sqr(out, out, 32) // ... ffffffff00000000

	p256Sqr(out, out, 32)
	p256Mul(out, out, p32)

	p256Sqr(out, out, 16)
	p256Mul(out, out, p16)

	p256Sqr(out, out, 8)
	p256Mul(out, out, p8)

	p256Sqr(out, out, 4)
	p256Mul(out, out, p4)

	p256Sqr(out, out, 2)
	p256Mul(out, out, p2)
}

// p256ExpQuadP sets out to in^-1 mod p.
func p256Sqrt(out, in []uint64) bool {
	var xp, a0, a1 [4]uint64
	p256Mul(xp[:], in, rr)
	p256ExpQuadP(a1[:], xp[:])
	p256Mul(out, a1[:], xp[:])
	p256Mul(a0[:], out, a1[:])
	p256FromMont(out, out)
	return a0 == p256MontOneA
	/*
		return a0[0] == p256MontOne[0] && a0[1] == p256MontOne[1] &&
			a0[2] == p256MontOne[2] && a0[3] == p256MontOne[3]
	*/
}

func (p *p256Point) p256StorePoint(r *[16 * 4 * 3]uint64, index int) {
	copy(r[index*12:], p.xyz[:])
}

func boothW5(in uint) (int, int) {
	var s uint = ^((in >> 5) - 1)
	var d uint = (1 << 6) - in - 1
	d = (d & s) | (in & (^s))
	d = (d >> 1) + (d & 1)
	return int(d), int(s & 1)
}

func boothW6(in uint) (int, int) {
	var s uint = ^((in >> 6) - 1)
	var d uint = (1 << 7) - in - 1
	d = (d & s) | (in & (^s))
	d = (d >> 1) + (d & 1)
	return int(d), int(s & 1)
}

func initTable() {
	p256Precomputed = new([43][32 * 8]uint64)

	// basePoint in montgomery form
	// filled later w/ Gx,Gy, MontOne
	basePoint := []uint64{
		0, 0, 0xd6a1ed99ac24c3c3, 0x91167a5ee1c13b05,
		0, 0, 0x8d4cfb066e2a48f8, 0x63cd65d481d735bd,
		0x0000000000000001, 0xffffffff, 0, 0x100000000,
	}
	// convert Gx, Gy to montgomery form
	fromBig(basePoint[:4], sm2Params.Gx)
	fromBig(basePoint[4:8], sm2Params.Gy)
	p256Mul(basePoint[:4], basePoint[:4], rr)
	p256Mul(basePoint[4:8], basePoint[4:8], rr)
	copy(basePoint[8:], p256MontOne)
	t1 := make([]uint64, 12)
	t2 := make([]uint64, 12)
	copy(t2, basePoint)

	zInv := make([]uint64, 4)
	zInvSq := make([]uint64, 4)
	for j := 0; j < 32; j++ {
		copy(t1, t2)
		for i := 0; i < 43; i++ {
			// The window size is 6 so we need to double 6 times.
			if i != 0 {
				for k := 0; k < 6; k++ {
					p256PointDoubleAsm(t1, t1)
				}
			}
			// Convert the point to affine form. (Its values are
			// still in Montgomery form however.)
			p256Inverse(zInv, t1[8:12])
			p256Sqr(zInvSq, zInv, 1)
			p256Mul(zInv, zInv, zInvSq)

			p256Mul(t1[:4], t1[:4], zInvSq)
			p256Mul(t1[4:8], t1[4:8], zInv)

			copy(t1[8:12], basePoint[8:12])
			// Update the table entry
			copy(p256Precomputed[i][j*8:], t1[:8])
		}
		if j == 0 {
			p256PointDoubleAsm(t2, basePoint)
		} else {
			p256PointAddAsm(t2, t2, basePoint)
		}
	}
}

func (p *p256Point) p256BaseMult(scalar []uint64) {
	precomputeOnce.Do(initTable)

	wvalue := (scalar[0] << 1) & 0x7f
	sel, sign := boothW6(uint(wvalue))
	p256SelectBase(p.xyz[0:8], p256Precomputed[0][0:], sel)
	p256NegCond(p.xyz[4:8], sign)

	// (This is one, in the Montgomery domain.)
	copy(p.xyz[8:], p256MontOne)

	var t0 p256Point
	// (This is one, in the Montgomery domain.)
	copy(t0.xyz[8:], p256MontOne)

	index := uint(5)
	zero := sel

	for i := 1; i < 43; i++ {
		if index < 192 {
			wvalue = ((scalar[index/64] >> (index % 64)) + (scalar[index/64+1] << (64 - (index % 64)))) & 0x7f
		} else {
			wvalue = (scalar[index/64] >> (index % 64)) & 0x7f
		}
		index += 6
		sel, sign = boothW6(uint(wvalue))
		p256SelectBase(t0.xyz[0:8], p256Precomputed[i][0:], sel)
		p256PointAddAffineAsm(p.xyz[0:12], p.xyz[0:12], t0.xyz[0:8], sign, sel, zero)
		zero |= sel
	}
}

func (p *p256Point) p256ScalarMult(scalar []uint64) {
	// precomp is a table of precomputed points that stores powers of p
	// from p^1 to p^16.
	var precomp [16 * 4 * 3]uint64
	var t0, t1, t2, t3 p256Point

	// Prepare the table
	p.p256StorePoint(&precomp, 0) // 1

	p256PointDoubleAsm(t0.xyz[:], p.xyz[:])
	p256PointDoubleAsm(t1.xyz[:], t0.xyz[:])
	p256PointDoubleAsm(t2.xyz[:], t1.xyz[:])
	p256PointDoubleAsm(t3.xyz[:], t2.xyz[:])
	t0.p256StorePoint(&precomp, 1)  // 2
	t1.p256StorePoint(&precomp, 3)  // 4
	t2.p256StorePoint(&precomp, 7)  // 8
	t3.p256StorePoint(&precomp, 15) // 16

	p256PointAddAsm(t0.xyz[:], t0.xyz[:], p.xyz[:])
	p256PointAddAsm(t1.xyz[:], t1.xyz[:], p.xyz[:])
	p256PointAddAsm(t2.xyz[:], t2.xyz[:], p.xyz[:])
	t0.p256StorePoint(&precomp, 2) // 3
	t1.p256StorePoint(&precomp, 4) // 5
	t2.p256StorePoint(&precomp, 8) // 9

	p256PointDoubleAsm(t0.xyz[:], t0.xyz[:])
	p256PointDoubleAsm(t1.xyz[:], t1.xyz[:])
	t0.p256StorePoint(&precomp, 5) // 6
	t1.p256StorePoint(&precomp, 9) // 10

	p256PointAddAsm(t2.xyz[:], t0.xyz[:], p.xyz[:])
	p256PointAddAsm(t1.xyz[:], t1.xyz[:], p.xyz[:])
	t2.p256StorePoint(&precomp, 6)  // 7
	t1.p256StorePoint(&precomp, 10) // 11

	p256PointDoubleAsm(t0.xyz[:], t0.xyz[:])
	p256PointDoubleAsm(t2.xyz[:], t2.xyz[:])
	t0.p256StorePoint(&precomp, 11) // 12
	t2.p256StorePoint(&precomp, 13) // 14

	p256PointAddAsm(t0.xyz[:], t0.xyz[:], p.xyz[:])
	p256PointAddAsm(t2.xyz[:], t2.xyz[:], p.xyz[:])
	t0.p256StorePoint(&precomp, 12) // 13
	t2.p256StorePoint(&precomp, 14) // 15

	// Start scanning the window from top bit
	index := uint(254)
	var sel, sign int

	wvalue := (scalar[index/64] >> (index % 64)) & 0x3f
	sel, _ = boothW5(uint(wvalue))

	p256Select(p.xyz[0:12], precomp[0:], sel)
	zero := sel

	for index > 4 {
		index -= 5
		p256PointDoubleAsm(p.xyz[:], p.xyz[:])
		p256PointDoubleAsm(p.xyz[:], p.xyz[:])
		p256PointDoubleAsm(p.xyz[:], p.xyz[:])
		p256PointDoubleAsm(p.xyz[:], p.xyz[:])
		p256PointDoubleAsm(p.xyz[:], p.xyz[:])

		if index < 192 {
			wvalue = ((scalar[index/64] >> (index % 64)) + (scalar[index/64+1] << (64 - (index % 64)))) & 0x3f
		} else {
			wvalue = (scalar[index/64] >> (index % 64)) & 0x3f
		}

		sel, sign = boothW5(uint(wvalue))

		p256Select(t0.xyz[0:], precomp[0:], sel)
		p256NegCond(t0.xyz[4:8], sign)
		p256PointAddAsm(t1.xyz[:], p.xyz[:], t0.xyz[:])
		p256MovCond(t1.xyz[0:12], t1.xyz[0:12], p.xyz[0:12], sel)
		p256MovCond(p.xyz[0:12], t1.xyz[0:12], t0.xyz[0:12], zero)
		zero |= sel
	}

	p256PointDoubleAsm(p.xyz[:], p.xyz[:])
	p256PointDoubleAsm(p.xyz[:], p.xyz[:])
	p256PointDoubleAsm(p.xyz[:], p.xyz[:])
	p256PointDoubleAsm(p.xyz[:], p.xyz[:])
	p256PointDoubleAsm(p.xyz[:], p.xyz[:])

	wvalue = (scalar[0] << 1) & 0x3f
	sel, sign = boothW5(uint(wvalue))

	p256Select(t0.xyz[0:], precomp[0:], sel)
	p256NegCond(t0.xyz[4:8], sign)
	p256PointAddAsm(t1.xyz[:], p.xyz[:], t0.xyz[:])
	p256MovCond(t1.xyz[0:12], t1.xyz[0:12], p.xyz[0:12], sel)
	p256MovCond(p.xyz[0:12], t1.xyz[0:12], t0.xyz[0:12], zero)
}
