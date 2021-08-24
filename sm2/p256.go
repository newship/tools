// Copyright 2013 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build !amd64,!arm64	field

package sm2

// This file contains a constant-time, 32-bit implementation of SM2 Prime 256.

import (
	"errors"
	"math/big"
	"unsafe"
)

const (
	p256Limbs     = 9
	p256LimbsProd = 18
)

var errParam = errors.New("error parameter")
var errSqrt = errors.New("error sqrt")

type p256Curve struct {
	*CurveParams
}

type feType [p256Limbs]uint32

var (
	// p256FeRR contains RR mod p - the Montgomery constant
	// (2**257)^2.
	p256FeRR feType
	pSM2     = p256Curve{CurveParams: sm2Params}
	bigOne   = new(big.Int).SetInt64(1)
)

func initSM2() {
	// p256One is the number 1 as a field element.
	tmp := new(big.Int).Lsh(bigOne, 257)
	tmp.Mod(tmp, sm2Params.P)
	p256FromBigRaw(&p256One, tmp)
	tmp.Mul(tmp, tmp)
	tmp.Mod(tmp, sm2Params.P)
	p256FromBigRaw(&p256FeRR, tmp)
	// p256P is the prime modulus as a field element.
	p256FromBigRaw(&p256P, sm2Params.P)
	tmp.Lsh(sm2Params.P, 1)
	// p2562P is the twice prime modulus as a field element.
	p256FromBigRaw(&p2562P, tmp)
	//p256FromBig(&p256One, bigOne)
	initTable()
}

func (curve p256Curve) Params() *CurveParams {
	return curve.CurveParams
}

// p256GetScalar endian-swaps the big-endian scalar value from in and writes it
// to out. If the scalar is equal or greater than the order of the group, it's
// reduced modulo that order.
func p256GetScalar(out *[32]byte, in []byte) {
	/*
		n := new(big.Int).SetBytes(in)
		var scalarBytes []byte

		if n.Cmp(sm2Params.N) >= 0 {
			n.Mod(n, sm2Params.N)
			scalarBytes = n.Bytes()
		} else {
			scalarBytes = in
		}
	*/
	var scalarBytes [32]byte
	copy(scalarBytes[:], in)
	for i := 0; i < 32; i++ {
		out[32-(1+i)] = scalarBytes[i]
	}
}

func RecoverPoint(x1 *big.Int, v uint) (y1 *big.Int, err error) {
	var xp, t1 feType
	c := pSM2
	if x1.Sign() <= 0 || x1.Cmp(c.N) >= 0 {
		return nil, errParam
	}
	p256FromBig(&xp, x1)
	p256Square(&t1, &xp)
	p256Mult(&t1, &t1, &xp)
	// t1 = x1^3
	p256Diff(&t1, &t1, &xp)
	p256Diff(&t1, &t1, &xp)
	p256Diff(&t1, &t1, &xp)
	p256FromBig(&xp, c.B)
	p256Sum(&t1, &t1, &xp)
	if !p256Sqrt(&xp, &t1) {
		return nil, errSqrt
	}
	tt := p256ToBig(&xp)
	if (v ^ tt.Bit(0)) != 0 {
		tt.Sub(c.P, tt)
	}
	return tt, nil
}

func (curve p256Curve) ScalarBaseMult(scalar []byte) (x, y *big.Int) {
	//return curve.ScalarMult(curve.Gx, curve.Gy, scalar)
	var scalarReversed [32]byte
	p256GetScalar(&scalarReversed, scalar)

	var x1, y1, z1 feType
	p256ScalarBaseMult(&x1, &y1, &z1, &scalarReversed)
	return p256ToAffine(&x1, &y1, &z1)
}

func (curve p256Curve) ScalarMult(bigX, bigY *big.Int, scalar []byte) (x, y *big.Int) {
	var scalarReversed [32]byte
	p256GetScalar(&scalarReversed, scalar)

	var px, py, x1, y1, z1 feType
	p256FromBig(&px, bigX)
	p256FromBig(&py, bigY)
	//z1 = p256One
	p256ScalarMult(&x1, &y1, &z1, &px, &py, &scalarReversed)
	return p256ToAffine(&x1, &y1, &z1)
}

// Field elements are represented as nine, unsigned 32-bit words.
//
// The value of an field element is:
//   x[0] + (x[1] * 2**29) + (x[2] * 2**57) + ... + (x[8] * 2**228)
//
// That is, each limb is alternately 29 or 28-bits wide in little-endian
// order.
//
// This means that a field element hits 2**257, rather than 2**256 as we would
// like. A 28, 29, ... pattern would cause us to hit 2**256, but that causes
// problems when multiplying as terms end up one bit short of a limb which
// would require much bit-shifting to correct.
//
// Finally, the values stored in a field element are in Montgomery form. So the
// value |y| is stored as (y*R) mod p, where p is the P-256 prime and R is
// 2**257.

const (
	bottom28Bits = 0xfffffff
	bottom29Bits = 0x1fffffff
)

var (
	// p256One is the number 1 as a field element.
	p256One  feType
	p256Zero = feType{0, 0, 0, 0, 0, 0, 0, 0, 0}
	// p256P is the prime modulus as a field element.
	p256P feType
	// p2562P is the twice prime modulus as a field element.
	p2562P feType
)

// p256Precomputed contains precomputed values to aid the calculation of scalar
// multiples of the base point, G. It's actually two, equal length, tables
// concatenated.
//
// The first table contains (x,y) field element pairs for 16 multiples of the
// base point, G.
//
//   Index  |  Index (binary) | Value
//       0  |           0000  | 0G (all zeros, omitted)
//       1  |           0001  | G
//       2  |           0010  | 2**64G
//       3  |           0011  | 2**64G + G
//       4  |           0100  | 2**128G
//       5  |           0101  | 2**128G + G
//       6  |           0110  | 2**128G + 2**64G
//       7  |           0111  | 2**128G + 2**64G + G
//       8  |           1000  | 2**192G
//       9  |           1001  | 2**192G + G
//      10  |           1010  | 2**192G + 2**64G
//      11  |           1011  | 2**192G + 2**64G + G
//      12  |           1100  | 2**192G + 2**128G
//      13  |           1101  | 2**192G + 2**128G + G
//      14  |           1110  | 2**192G + 2**128G + 2**64G
//      15  |           1111  | 2**192G + 2**128G + 2**64G + G
//
// The second table follows the same style, but the terms are 2**32G,
// 2**96G, 2**160G, 2**224G.
//
// This is ~2KB of data.
var p256Precomputed = [p256Limbs * 2 * 15 * 2]uint32{
	0x11522878, 0xe730d41, 0xdb60179, 0x4afe2ff, 0x12883add, 0xcaddd88, 0x119e7edc, 0xd4a6eab, 0x3120bee,
	0x1d2aac15, 0xf25357c, 0x19e45cdd, 0x5c721d0, 0x1992c5a5, 0xa237487, 0x154ba21, 0x14b10bb, 0xae3fe3,
}

func setPrecomputed(idx uint, x, y *feType) {
	if idx >= 2*15 {
		return
	}
	idx *= p256Limbs * 2
	//var xx, yy *feType
	xx := (*feType)(unsafe.Pointer(&p256Precomputed[idx]))
	idx += p256Limbs
	yy := (*feType)(unsafe.Pointer(&p256Precomputed[idx]))
	*xx = *x
	*yy = *y
}

func getPrecomputed(idx uint) (x, y *feType) {
	if idx >= 2*15 {
		return
	}
	idx *= p256Limbs * 2
	//var xx, yy *feType
	x = (*feType)(unsafe.Pointer(&p256Precomputed[idx]))
	idx += p256Limbs
	y = (*feType)(unsafe.Pointer(&p256Precomputed[idx]))
	return
}

func initTable() {
	var x, y, z feType
	var xx, yy feType
	p256FromBig(&x, sm2Params.Gx)
	p256FromBig(&y, sm2Params.Gy)
	setPrecomputed(0, &x, &y)
	for i := uint(1); i <= 8; i <<= 1 {
		x2, y2 := getPrecomputed(i - 1) // [0][i]
		x = *x2
		y = *y2
		z = p256One
		for j := 0; j < 32; j++ {
			p256PointDouble(&x, &y, &z, &x, &y, &z)
		}
		p256PointToAffine(&xx, &yy, &x, &y, &z)
		setPrecomputed(14+i, &xx, &yy) // [1][i]
		if i == 8 {
			break
		}
		x = xx
		y = yy
		z = p256One
		for j := 0; j < 32; j++ {
			p256PointDouble(&x, &y, &z, &x, &y, &z)
		}
		p256PointToAffine(&xx, &yy, &x, &y, &z)
		setPrecomputed(2*i-1, &xx, &yy) // [0][2*i]
	}
	for i := uint(0); i < 2; i++ {
		x1, y1 := getPrecomputed(i*15 + 1) // [i][2]
		x = *x1
		y = *y1
		x2, y2 := getPrecomputed(i*15 + 3) // [i][4]
		p256PointAddMixed(&x, &y, &z, x1, y1, &p256One, x2, y2)
		p256PointToAffine(&xx, &yy, &x, &y, &z)
		setPrecomputed(i*15+5, &xx, &yy)   // [i][6]
		x3, y3 := getPrecomputed(i*15 + 7) // [i][8]
		p256PointAddMixed(&x, &y, &z, x1, y1, &p256One, x3, y3)
		p256PointToAffine(&xx, &yy, &x, &y, &z)
		setPrecomputed(i*15+9, &xx, &yy) // [i][10]
		p256PointAddMixed(&x, &y, &z, x2, y2, &p256One, x3, y3)
		p256PointToAffine(&xx, &yy, &x, &y, &z)
		setPrecomputed(i*15+11, &xx, &yy) // [i][12]
		p256PointAddMixed(&x, &y, &z, &xx, &yy, &p256One, x1, y1)
		p256PointToAffine(&xx, &yy, &x, &y, &z)
		setPrecomputed(i*15+13, &xx, &yy) // [i][14]
		for j := uint(1); j < 8; j++ {
			x1, y1 := getPrecomputed(i*15 + j*2 - 1) // [i][2*j]
			x2, y2 := getPrecomputed(i * 15)         // [i][1]
			p256PointAddMixed(&x, &y, &z, x2, y2, &p256One, x1, y1)
			p256PointToAffine(&xx, &yy, &x, &y, &z)
			setPrecomputed(i*15+j*2, &xx, &yy) // [i][2*j]
		}
	}
}

// Field element operations:

// nonZeroToAllOnes returns:
//   0xffffffff for 0 < x <= 2**31
//   0 for x == 0 or x > 2**31.
func nonZeroToAllOnes(x uint32) uint32 {
	return ((x - 1) >> 31) - 1
}

// p256ReduceCarry adds a multiple of p in order to cancel |carry|,
// which is a term at 2**257.
//
// On entry: carry < 2**3, inout[0,2,...] < 2**29, inout[1,3,...] < 2**28.
// On exit: inout[0,2,..] < 2**30, inout[1,3,...] < 2**29.
func p256ReduceCarry(inout *feType, carry uint32) {
	carry_mask := nonZeroToAllOnes(carry)

	inout[0] += carry << 1
	// carry < 2**3 thus (carry << 8) < 2**11 and we added 2**28 in the
	// previous line therefore this doesn't underflow.
	// -2⁶⁵
	inout[2] += 0x20000000 & carry_mask
	inout[2] -= carry << 8
	inout[3] -= 1 & carry_mask
	// +2⁹⁷
	inout[3] += carry << 11
	// +2²²⁵
	inout[7] += carry << 25
	cc := inout[7] >> 28
	inout[7] &= bottom28Bits
	inout[8] += cc
}

// p256Sum sets out = in+in2.
//
// On entry, in[i]+in2[i] must not overflow a 32-bit word.
// On exit: out[0,2,...] < 2**30, out[1,3,...] < 2**29
func p256Sum(out, in, in2 *feType) {
	carry := uint32(0)
	for i := 0; ; i++ {
		out[i] = in[i] + in2[i]
		out[i] += carry
		carry = out[i] >> 29
		out[i] &= bottom29Bits

		i++
		if i == p256Limbs {
			break
		}

		out[i] = in[i] + in2[i]
		out[i] += carry
		carry = out[i] >> 28
		out[i] &= bottom28Bits
	}
	p256ReduceCarry(out, carry)
}

// p256NormalWeak sets out = in normal 29,28...
//
// On entry, in[i] must not overflow a 32-bit word.
// On exit: out[0,2,...] < 2**29, out[1,3,...] < 2**28
// On exit: out < 2²⁵⁶, < 2P
func p256NormalWeak(inout *feType) {
	carry := uint32(0)
	for i := 0; ; i++ {
		inout[i] += carry
		carry = inout[i] >> 29
		inout[i] &= bottom29Bits

		i++
		if i == p256Limbs {
			break
		}

		inout[i] += carry
		carry = inout[i] >> 28
		inout[i] &= bottom28Bits
	}

	carry <<= 1
	carry |= inout[8] >> 28
	inout[8] &= bottom28Bits
	// carry < 2³
	for carry != 0 {
		// try reduce carry over 2^256
		// +1
		inout[0] += carry
		cc := inout[0] >> 29
		inout[0] &= bottom29Bits
		inout[1] += cc
		cc = inout[1] >> 28
		inout[1] &= bottom28Bits
		// carry < 2**3 thus (carry << 8) < 2**11 and we added 2**28 in the
		// previous line therefore this doesn't underflow.
		// -2⁶⁴
		inout[2] += 0x20000000 + cc
		inout[2] -= carry << 7
		cc = inout[2] >> 29
		inout[2] &= bottom29Bits
		// +2⁹⁶
		inout[3] += (carry << 10)
		inout[3] += cc - 1
		cc = inout[3] >> 28
		inout[3] &= bottom28Bits
		inout[4] += cc
		cc = inout[4] >> 29
		inout[4] &= bottom29Bits
		inout[5] += cc
		cc = inout[5] >> 28
		inout[5] &= bottom28Bits
		inout[6] += cc
		cc = inout[6] >> 29
		inout[6] &= bottom29Bits
		// +2²²⁴
		inout[7] += carry << 24
		inout[7] += cc
		cc = inout[7] >> 28
		inout[7] &= bottom28Bits
		inout[8] += cc
		carry = inout[8] >> 28
		inout[8] &= bottom28Bits
	}
}

const (
	two31m3    = 1<<31 - 1<<3
	two30m2    = 1<<30 - 1<<2
	two31p10m2 = 1<<31 + 1<<10 - 1<<2
	two30m13m2 = 1<<30 - 1<<13 - 1<<2
	two31m2    = 1<<31 - 1<<2
	two30m27m2 = 1<<30 - 1<<27 - 1<<2
)

// p256Zero31 is 0 mod p.
var p256Zero31 = feType{two31m3, two30m2, two31p10m2, two30m13m2, two31m2, two30m2, two31m2, two30m27m2, two31m2}

// p256Diff sets out = in-in2.
//
// On entry: in[0,2,...] < 2**30, in[1,3,...] < 2**29 and
//           in2[0,2,...] < 2**30, in2[1,3,...] < 2**29.
// On exit: out[0,2,...] < 2**30, out[1,3,...] < 2**29.
func p256Diff(out, in, in2 *feType) {
	var carry uint32

	for i := 0; ; i++ {
		out[i] = in[i] - in2[i]
		out[i] += p256Zero31[i]
		out[i] += carry
		carry = out[i] >> 29
		out[i] &= bottom29Bits

		i++
		if i == p256Limbs {
			break
		}

		out[i] = in[i] - in2[i]
		out[i] += p256Zero31[i]
		out[i] += carry
		carry = out[i] >> 28
		out[i] &= bottom28Bits
	}

	p256ReduceCarry(out, carry)
}

// p256Mod sets out = in modulo prime P.
//
// On entry: in[0,2,...] < 2**29, in[1,3,...] < 2**28
// On exit: out[0,2,...] < 2**29, out[1,3,...] < 2**28.
func p256Mod(out, in *feType) {
	carry := uint32(0)

	for i := 0; ; i++ {
		out[i] = in[i] - p256P[i] - carry
		carry = out[i] >> 29
		if carry != 0 {
			carry = 1
		}
		out[i] &= bottom29Bits

		i++
		if i == p256Limbs {
			break
		}

		out[i] = in[i] - p256P[i] - carry
		carry = out[i] >> 28
		if carry != 0 {
			carry = 1
		}
		out[i] &= bottom28Bits
	}

	// if carry not zero, copy in
	xMask := nonZeroToAllOnes(carry)
	p256CopyConditional(out, in, xMask)
}

// p256Equal returns  in == in2.
//
// On entry: in[0,2,...] < 2**30, in[1,3,...] < 2**29 and
//           in2[0,2,...] < 2**30, in2[1,3,...] < 2**29.
// On exit: true/false
func p256Equal(in, in2 *feType) bool {
	if *in == *in2 {
		return true
	}
	var res feType
	p256Diff(&res, in, in2)
	p256NormalWeak(&res)
	if res == p256P {
		return true
	}
	return res == p256Zero
}

// p256Reduction sets out = in/R mod p where tmp contains 32-bit words with
// the same 29,28,... bit positions as an field element.
//
// The values in field elements are in Montgomery form: x*R mod p where R =
// 2**257.  We wish to divide by R in order for the result also to be
// converted from Montgomery form.
//
// On entry: in[0,2,...] < 2**30, in[1,3,...] < 2**29
// On exit: out[0,2,...] < 2**30, out[1,3,...] < 2**29
func p256Reduction(inout *[p256LimbsProd]uint32) {
	var x, xMask uint32
	// Montgomery elimination of terms:
	//
	// Since R is 2**257, we can divide by R with a bitwise shift if we can
	// ensure that the right-most 257 bits are all zero. We can make that true
	// by adding multiplies of p without affecting the value.
	//
	// So we eliminate limbs from right to left. Since the bottom 29 bits of p
	// are all ones, then by adding in[0]*p to in we'll make in[0] == 0.
	// We can do that for 8 further limbs and then right shift to eliminate the
	// extra factor of R.
	for i := 0; ; i += 2 {
		inout[i+1] += inout[i] >> 29
		x = inout[i] & bottom29Bits
		xMask = nonZeroToAllOnes(x)
		inout[i] = 0

		// The bounds calculations for this loop are tricky. Each iteration of
		// the loop eliminates two words by adding values to words to their
		// right.
		//
		// The following table contains the amounts added to each word (as an
		// offset from the value of i at the top of the loop). The amounts are
		// accounted for from the first and second half of the loop separately
		// and are written as, for example, 28 to mean a value <2**28.
		//
		// Word:                   3   4   5   6   7   8   9   10
		// Added in top half:     28  11      29  21  29  28
		//                                        28  29
		//                                            29
		// Added in bottom half:      29  10      28  21  28   28
		//                                            29
		//
		// The value that is currently offset 7 will be offset 5 for the next
		// iteration and then offset 3 for the iteration after that. Therefore
		// the total value added will be the values added at 7, 5 and 3.
		//
		// The following table accumulates these values. The sums at the bottom
		// are written as, for example, 29+28, to mean a value < 2**29+2**28.
		//
		// Word:                   3   4   5   6   7   8   9  10  11  12  13
		//                        28  11  10  29  21  29  28  28  28  28  28
		//                            29  28  11  28  29  28  29  28  29  28
		//                                    29  28  21  21  29  21  29  21
		//                                        10  29  28  21  28  21  28
		//                                        28  29  28  29  28  29  28
		//                                            11  10  29  10  29  10
		//                                            29  28  11  28  11
		//                                                    29      29
		//                        --------------------------------------------
		//                                                30+ 31+ 30+ 31+ 30+
		//                                                28+ 29+ 28+ 29+ 21+
		//                                                21+ 28+ 21+ 28+ 10
		//                                                10  21+ 10  21+
		//                                                    11      11
		//
		// So the greatest amount is added to inout[10] and inout[12]. If
		// inout[10/12] has an initial value of <2**29, then the maximum value
		// will be < 2**31 + 2**30 + 2**28 + 2**21 + 2**11, which is < 2**32,
		// as required.
		// 2⁶⁴
		inout[i+2] += (x << 7) & bottom29Bits
		inout[i+3] += (x >> 22)

		// -2⁹⁶
		inout[i+3] += 0x10000000 & xMask
		inout[i+4] += (0x20000000 - 1) & xMask
		inout[i+3] -= (x << 10) & bottom28Bits
		inout[i+4] -= (x >> 18)

		inout[i+5] += (0x10000000 - 1) & xMask
		inout[i+6] += (0x20000000 - 1) & xMask

		// At position 200, which is the starting bit position for word 7, we
		// have a factor of 0xf000000 = 2**28 - 2**24.
		// -2²²⁴
		inout[i+7] += (0x10000000 - 1) & xMask
		inout[i+7] -= (x << 24) & bottom28Bits
		// trickly add 2^28, (x-1)^256 + 2^256-2^224
		inout[i+8] += (0x10000000 - 1) & xMask
		inout[i+8] -= (x >> 4)
		// 2²⁵⁶
		inout[i+8] += ((x - 1) << 28) & xMask & bottom29Bits
		inout[i+9] += ((x - 1) >> 1) & xMask

		if i+1 == p256Limbs {
			break
		}
		inout[i+2] += inout[i+1] >> 28
		x = inout[i+1] & bottom28Bits
		xMask = nonZeroToAllOnes(x)
		inout[i+1] = 0

		// 2⁶⁴
		inout[i+3] += (x << 7) & bottom28Bits
		inout[i+4] += (x >> 21)

		// -2⁹⁶
		inout[i+4] += 0x20000000 & xMask
		inout[i+5] += (0x10000000 - 1) & xMask
		inout[i+4] -= (x << 11) & bottom29Bits
		inout[i+5] -= (x >> 18)

		inout[i+6] += (0x20000000 - 1) & xMask
		inout[i+7] += (0x10000000 - 1) & xMask

		// At position 199, which is the starting bit of the 8th word when
		// dealing with a context starting on an odd word, we have a factor of
		// 0x1e000000 = 2**29 - 2**25. Since we have not updated i, the 8th
		// word from i+1 is i+8.
		inout[i+8] += (0x20000000 - 1) & xMask
		inout[i+8] -= (x << 25) & bottom29Bits
		inout[i+9] += (0x10000000 - 1) & xMask
		inout[i+9] -= (x >> 4)
		// 2²⁵⁶
		inout[i+10] += (x - 1) & xMask
	}
}

func p256ProdHigh(out *feType, in *[p256LimbsProd]uint32) {
	// We merge the right shift with a carry chain. The words above 2**257 have
	// widths of 28,29,... which we need to correct when copying them down.
	carry := uint32(0)
	for i := 0; i < 8; i++ {
		// The maximum value of in[i + 9] occurs on the first iteration and
		// is < 2**30+2**29+2**28. Adding 2**29 (from in[i + 10]) is
		// therefore safe.
		out[i] = in[i+9]
		out[i] += carry
		out[i] += (in[i+10] << 28) & bottom29Bits
		carry = out[i] >> 29
		out[i] &= bottom29Bits

		i++
		out[i] = in[i+9] >> 1
		out[i] += carry
		carry = out[i] >> 28
		out[i] &= bottom28Bits
	}
	out[8] = in[17]
	out[8] += carry
	carry = out[8] >> 29
	out[8] &= bottom29Bits

	p256ReduceCarry(out, carry)
}

// p256ReduceDegree sets out = tmp/R mod p where tmp contains 64-bit words with
// the same 29,28,... bit positions as an field element.
//
// The values in field elements are in Montgomery form: x*R mod p where R =
// 2**257. Since we just multiplied two Montgomery values together, the result
// is x*y*R*R mod p. We wish to divide by R in order for the result also to be
// in Montgomery form.
//
// On entry: tmp[i] < 2**64
// On exit: out[0,2,...] < 2**30, out[1,3,...] < 2**29
func p256ReduceDegree(out *feType, tmp [17]uint64) {
	// The following table may be helpful when reading this code:
	//
	// Limb number:   0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10...
	// Width (bits):  29| 28| 29| 28| 29| 28| 29| 28| 29| 28| 29
	// Start bit:     0 | 29| 57| 86|114|143|171|200|228|257|285
	//   (odd phase): 0 | 28| 57| 85|114|142|171|199|228|256|285
	var tmp2 [p256LimbsProd]uint32
	var carry uint32

	// tmp contains 64-bit words with the same 29,28,29-bit positions as an
	// field element. So the top of an element of tmp might overlap with
	// another element two positions down. The following loop eliminates
	// this overlap.
	tmp2[0] = uint32(tmp[0]) & bottom29Bits

	tmp2[1] = uint32(tmp[0]>>29) & bottom28Bits
	tmp2[1] += uint32(tmp[1]) & bottom28Bits
	carry = tmp2[1] >> 28
	tmp2[1] &= bottom28Bits

	for i := 2; i < 17; i++ {
		tmp2[i] = uint32(tmp[i-2] >> 57)
		tmp2[i] += uint32(tmp[i-1]>>28) & bottom29Bits
		tmp2[i] += uint32(tmp[i]) & bottom29Bits
		tmp2[i] += carry
		carry = tmp2[i] >> 29
		tmp2[i] &= bottom29Bits

		i++
		if i == 17 {
			break
		}
		tmp2[i] = uint32(tmp[i-2] >> 57)
		tmp2[i] += uint32(tmp[i-1]>>29) & bottom28Bits
		tmp2[i] += uint32(tmp[i]) & bottom28Bits
		tmp2[i] += carry
		carry = tmp2[i] >> 28
		tmp2[i] &= bottom28Bits
	}

	tmp2[17] = uint32(tmp[15] >> 57)
	tmp2[17] += uint32(tmp[16] >> 29)
	tmp2[17] += carry
	p256Reduction(&tmp2)
	p256ProdHigh(out, &tmp2)
}

// p256Square sets out=in*in.
//
// On entry: in[0,2,...] < 2**30, in[1,3,...] < 2**29.
// On exit: out[0,2,...] < 2**30, out[1,3,...] < 2**29.
func p256Square(out, in *feType) {
	var tmp [17]uint64

	tmp[0] = uint64(in[0]) * uint64(in[0])
	tmp[1] = uint64(in[0]) * (uint64(in[1]) << 1)
	tmp[2] = uint64(in[0])*(uint64(in[2])<<1) +
		uint64(in[1])*(uint64(in[1])<<1)
	tmp[3] = uint64(in[0])*(uint64(in[3])<<1) +
		uint64(in[1])*(uint64(in[2])<<1)
	tmp[4] = uint64(in[0])*(uint64(in[4])<<1) +
		uint64(in[1])*(uint64(in[3])<<2) +
		uint64(in[2])*uint64(in[2])
	tmp[5] = uint64(in[0])*(uint64(in[5])<<1) +
		uint64(in[1])*(uint64(in[4])<<1) +
		uint64(in[2])*(uint64(in[3])<<1)
	tmp[6] = uint64(in[0])*(uint64(in[6])<<1) +
		uint64(in[1])*(uint64(in[5])<<2) +
		uint64(in[2])*(uint64(in[4])<<1) +
		uint64(in[3])*(uint64(in[3])<<1)
	tmp[7] = uint64(in[0])*(uint64(in[7])<<1) +
		uint64(in[1])*(uint64(in[6])<<1) +
		uint64(in[2])*(uint64(in[5])<<1) +
		uint64(in[3])*(uint64(in[4])<<1)
	// tmp[8] has the greatest value of 2**61 + 2**60 + 2**61 + 2**60 + 2**60,
	// which is < 2**64 as required.
	tmp[8] = uint64(in[0])*(uint64(in[8])<<1) +
		uint64(in[1])*(uint64(in[7])<<2) +
		uint64(in[2])*(uint64(in[6])<<1) +
		uint64(in[3])*(uint64(in[5])<<2) +
		uint64(in[4])*uint64(in[4])
	tmp[9] = uint64(in[1])*(uint64(in[8])<<1) +
		uint64(in[2])*(uint64(in[7])<<1) +
		uint64(in[3])*(uint64(in[6])<<1) +
		uint64(in[4])*(uint64(in[5])<<1)
	tmp[10] = uint64(in[2])*(uint64(in[8])<<1) +
		uint64(in[3])*(uint64(in[7])<<2) +
		uint64(in[4])*(uint64(in[6])<<1) +
		uint64(in[5])*(uint64(in[5])<<1)
	tmp[11] = uint64(in[3])*(uint64(in[8])<<1) +
		uint64(in[4])*(uint64(in[7])<<1) +
		uint64(in[5])*(uint64(in[6])<<1)
	tmp[12] = uint64(in[4])*(uint64(in[8])<<1) +
		uint64(in[5])*(uint64(in[7])<<2) +
		uint64(in[6])*uint64(in[6])
	tmp[13] = uint64(in[5])*(uint64(in[8])<<1) +
		uint64(in[6])*(uint64(in[7])<<1)
	tmp[14] = uint64(in[6])*(uint64(in[8])<<1) +
		uint64(in[7])*(uint64(in[7])<<1)
	tmp[15] = uint64(in[7]) * (uint64(in[8]) << 1)
	tmp[16] = uint64(in[8]) * uint64(in[8])

	p256ReduceDegree(out, tmp)
}

// p256Mult sets out=in*in2.
//
// On entry: in[0,2,...] < 2**30, in[1,3,...] < 2**29 and
//           in2[0,2,...] < 2**30, in2[1,3,...] < 2**29.
// On exit: out[0,2,...] < 2**30, out[1,3,...] < 2**29.
func p256Mult(out, in, in2 *feType) {
	var tmp [17]uint64

	tmp[0] = uint64(in[0]) * uint64(in2[0])
	tmp[1] = uint64(in[0])*(uint64(in2[1])<<0) +
		uint64(in[1])*(uint64(in2[0])<<0)
	tmp[2] = uint64(in[0])*(uint64(in2[2])<<0) +
		uint64(in[1])*(uint64(in2[1])<<1) +
		uint64(in[2])*(uint64(in2[0])<<0)
	tmp[3] = uint64(in[0])*(uint64(in2[3])<<0) +
		uint64(in[1])*(uint64(in2[2])<<0) +
		uint64(in[2])*(uint64(in2[1])<<0) +
		uint64(in[3])*(uint64(in2[0])<<0)
	tmp[4] = uint64(in[0])*(uint64(in2[4])<<0) +
		uint64(in[1])*(uint64(in2[3])<<1) +
		uint64(in[2])*(uint64(in2[2])<<0) +
		uint64(in[3])*(uint64(in2[1])<<1) +
		uint64(in[4])*(uint64(in2[0])<<0)
	tmp[5] = uint64(in[0])*(uint64(in2[5])<<0) +
		uint64(in[1])*(uint64(in2[4])<<0) +
		uint64(in[2])*(uint64(in2[3])<<0) +
		uint64(in[3])*(uint64(in2[2])<<0) +
		uint64(in[4])*(uint64(in2[1])<<0) +
		uint64(in[5])*(uint64(in2[0])<<0)
	tmp[6] = uint64(in[0])*(uint64(in2[6])<<0) +
		uint64(in[1])*(uint64(in2[5])<<1) +
		uint64(in[2])*(uint64(in2[4])<<0) +
		uint64(in[3])*(uint64(in2[3])<<1) +
		uint64(in[4])*(uint64(in2[2])<<0) +
		uint64(in[5])*(uint64(in2[1])<<1) +
		uint64(in[6])*(uint64(in2[0])<<0)
	tmp[7] = uint64(in[0])*(uint64(in2[7])<<0) +
		uint64(in[1])*(uint64(in2[6])<<0) +
		uint64(in[2])*(uint64(in2[5])<<0) +
		uint64(in[3])*(uint64(in2[4])<<0) +
		uint64(in[4])*(uint64(in2[3])<<0) +
		uint64(in[5])*(uint64(in2[2])<<0) +
		uint64(in[6])*(uint64(in2[1])<<0) +
		uint64(in[7])*(uint64(in2[0])<<0)
	// tmp[8] has the greatest value but doesn't overflow. See logic in
	// p256Square.
	tmp[8] = uint64(in[0])*(uint64(in2[8])<<0) +
		uint64(in[1])*(uint64(in2[7])<<1) +
		uint64(in[2])*(uint64(in2[6])<<0) +
		uint64(in[3])*(uint64(in2[5])<<1) +
		uint64(in[4])*(uint64(in2[4])<<0) +
		uint64(in[5])*(uint64(in2[3])<<1) +
		uint64(in[6])*(uint64(in2[2])<<0) +
		uint64(in[7])*(uint64(in2[1])<<1) +
		uint64(in[8])*(uint64(in2[0])<<0)
	tmp[9] = uint64(in[1])*(uint64(in2[8])<<0) +
		uint64(in[2])*(uint64(in2[7])<<0) +
		uint64(in[3])*(uint64(in2[6])<<0) +
		uint64(in[4])*(uint64(in2[5])<<0) +
		uint64(in[5])*(uint64(in2[4])<<0) +
		uint64(in[6])*(uint64(in2[3])<<0) +
		uint64(in[7])*(uint64(in2[2])<<0) +
		uint64(in[8])*(uint64(in2[1])<<0)
	tmp[10] = uint64(in[2])*(uint64(in2[8])<<0) +
		uint64(in[3])*(uint64(in2[7])<<1) +
		uint64(in[4])*(uint64(in2[6])<<0) +
		uint64(in[5])*(uint64(in2[5])<<1) +
		uint64(in[6])*(uint64(in2[4])<<0) +
		uint64(in[7])*(uint64(in2[3])<<1) +
		uint64(in[8])*(uint64(in2[2])<<0)
	tmp[11] = uint64(in[3])*(uint64(in2[8])<<0) +
		uint64(in[4])*(uint64(in2[7])<<0) +
		uint64(in[5])*(uint64(in2[6])<<0) +
		uint64(in[6])*(uint64(in2[5])<<0) +
		uint64(in[7])*(uint64(in2[4])<<0) +
		uint64(in[8])*(uint64(in2[3])<<0)
	tmp[12] = uint64(in[4])*(uint64(in2[8])<<0) +
		uint64(in[5])*(uint64(in2[7])<<1) +
		uint64(in[6])*(uint64(in2[6])<<0) +
		uint64(in[7])*(uint64(in2[5])<<1) +
		uint64(in[8])*(uint64(in2[4])<<0)
	tmp[13] = uint64(in[5])*(uint64(in2[8])<<0) +
		uint64(in[6])*(uint64(in2[7])<<0) +
		uint64(in[7])*(uint64(in2[6])<<0) +
		uint64(in[8])*(uint64(in2[5])<<0)
	tmp[14] = uint64(in[6])*(uint64(in2[8])<<0) +
		uint64(in[7])*(uint64(in2[7])<<1) +
		uint64(in[8])*(uint64(in2[6])<<0)
	tmp[15] = uint64(in[7])*(uint64(in2[8])<<0) +
		uint64(in[8])*(uint64(in2[7])<<0)
	tmp[16] = uint64(in[8]) * (uint64(in2[8]) << 0)

	p256ReduceDegree(out, tmp)
}

func p256Assign(out, in *feType) {
	*out = *in
}

func (out *feType) p256Sqr(in *feType, N int) {
	p256Square(out, in)
	for i := 1; i < N; i++ {
		p256Square(out, out)
	}
}

func precomputeExp(table *[5]feType, in *feType) {
	var ftmp feType
	// each e_I will hold |in|^{2^I - 1}
	e2 := &table[0]
	e4 := &table[1]
	e8 := &table[2]
	e16 := &table[3]
	e32 := &table[4]

	p256Square(&ftmp, in)    // 2^1
	p256Mult(e2, &ftmp, in)  // 2^2 - 2^0
	p256Square(&ftmp, e2)    // 2^3 - 2^1
	p256Square(&ftmp, &ftmp) // 2^4 - 2^2
	p256Mult(e4, &ftmp, e2)  // 2^4 - 2^0
	ftmp.p256Sqr(e4, 4)
	p256Mult(e8, &ftmp, e4) // 2^8 - 2^0
	ftmp.p256Sqr(e8, 8)
	p256Mult(e16, &ftmp, e8)  // 2^16 - 2^0
	ftmp.p256Sqr(e16, 16)     // 2^32 - 2^16
	p256Mult(e32, &ftmp, e16) // 2^32 - 2^0
}

// p256Invert calculates |out| = |in|^{-1}
//
// Based on Fermat's Little Theorem:
//   a^p = a (mod p)
//   a^{p-1} = 1 (mod p)
//   a^{p-2} = a^{-1} (mod p)
func p256Invert(out, in *feType) {
	var table [5]feType
	// each e_I will hold |in|^{2^I - 1}
	e2 := &table[0]
	e4 := &table[1]
	e8 := &table[2]
	e16 := &table[3]
	e32 := &table[4]
	precomputeExp(&table, in)

	// start calc in^p-2
	out.p256Sqr(e16, 8)
	p256Mult(out, out, e8) // 2^24 - 2^0
	out.p256Sqr(out, 4)
	p256Mult(out, out, e4) // 2^28 - 2^0
	p256Square(out, out)
	p256Square(out, out)
	p256Mult(out, out, e2) // 2^30 - 2^0
	p256Square(out, out)
	p256Mult(out, out, in) // 2^31 - 2^0
	p256Square(out, out)   // 2^32 - 2^1 = 2^32 -2^0 -2^0

	out.p256Sqr(out, 32)
	p256Mult(out, out, e32) // 2^64 - 2^32 - 2^0
	out.p256Sqr(out, 32)
	p256Mult(out, out, e32) // 2^96 - 2^64 - 2^0
	out.p256Sqr(out, 32)
	p256Mult(out, out, e32) // 2^128 - 2^96 - 2^0
	out.p256Sqr(out, 32)
	p256Mult(out, out, e32) // 2^160 - 2^128 - 2^0
	out.p256Sqr(out, 32)    // 2^192 - 2^160 - 2^32

	out.p256Sqr(out, 32)
	p256Mult(out, out, e32) // 2^224 - 2^192 - 2^64 + 2^32 - 1

	out.p256Sqr(out, 16)
	p256Mult(out, out, e16) // 2^240 - 2^208 - 2^80 + 2^48 - 1
	out.p256Sqr(out, 8)
	p256Mult(out, out, e8) // 2^248 - 2^216 - 2^88 + 2^56 - 1
	out.p256Sqr(out, 4)
	p256Mult(out, out, e4) // 2^252 - 2^220 - 2^92 + 2^60 - 1
	p256Square(out, out)
	p256Square(out, out)
	p256Mult(out, out, e2) // 2^254 - 2^222 - 2^94 + 2^62 - 1
	p256Square(out, out)
	p256Square(out, out)   // 2^256 - 2^224 - 2^96 + 2^64 - 4
	p256Mult(out, out, in) // 2^256 - 2^224 - 2^96 + 2^64 - 3
}

// p256ExpQuadP calculates |out| = |in|^{P/4}
//
func p256ExpQuadP(out, in *feType) {
	var table [5]feType
	// each e_I will hold |in|^{2^I - 1}
	e2 := &table[0]
	e4 := &table[1]
	e8 := &table[2]
	e16 := &table[3]
	e32 := &table[4]
	precomputeExp(&table, in)

	// start calc in^p-2
	out.p256Sqr(e16, 8)
	p256Mult(out, out, e8) // 2^24 - 2^0
	out.p256Sqr(out, 4)
	p256Mult(out, out, e4) // 2^28 - 2^0
	p256Square(out, out)
	p256Square(out, out)
	p256Mult(out, out, e2) // 2^30 - 2^0
	p256Square(out, out)
	p256Mult(out, out, in) // 2^31 - 2^0
	p256Square(out, out)   // 2^32 - 2^1 = 2^32 -2^0 -2^0

	out.p256Sqr(out, 32)
	p256Mult(out, out, e32) // 2^64 - 2^32 - 2^0
	out.p256Sqr(out, 32)
	p256Mult(out, out, e32) // 2^96 - 2^64 - 2^0
	out.p256Sqr(out, 32)
	p256Mult(out, out, e32) // 2^128 - 2^96 - 2^0
	out.p256Sqr(out, 32)
	p256Mult(out, out, e32) // 2^160 - 2^128 - 2^0
	out.p256Sqr(out, 32)    // 2^192 - 2^160 - 2^32

	out.p256Sqr(out, 32)
	p256Mult(out, out, e32) // 2^224 - 2^192 - 2^64 + 2^32 - 1

	out.p256Sqr(out, 16)
	p256Mult(out, out, e16) // 2^240 - 2^208 - 2^80 + 2^48 - 1
	out.p256Sqr(out, 8)
	p256Mult(out, out, e8) // 2^248 - 2^216 - 2^88 + 2^56 - 1
	out.p256Sqr(out, 4)
	p256Mult(out, out, e4) // 2^252 - 2^220 - 2^92 + 2^60 - 1
	p256Square(out, out)
	p256Square(out, out)
	p256Mult(out, out, e2) // 2^254 - 2^222 - 2^94 + 2^62 - 1
}

// p256Sqrt sets out=in^(1/2).
//
// On entry: in[0,2,...] < 2**30, in[1,3,...] < 2**29.
// On exit: out[0,2,...] < 2**30, out[1,3,...] < 2**29.
func p256Sqrt(out, in *feType) bool {
	var a0, a1 feType
	p256ExpQuadP(&a1, in)
	p256Mult(out, &a1, in)
	p256Mult(&a0, out, &a1)
	p256NormalWeak(&a0)
	return p256Equal(&a0, &p256One)
}

// p256Scalar3 sets out=3*out.
//
// On entry: out[0,2,...] < 2**30, out[1,3,...] < 2**29.
// On exit: out[0,2,...] < 2**30, out[1,3,...] < 2**29.
func p256Scalar3(out *feType) {
	var carry uint32

	for i := 0; ; i++ {
		out[i] *= 3
		out[i] += carry
		carry = out[i] >> 29
		out[i] &= bottom29Bits

		i++
		if i == p256Limbs {
			break
		}

		out[i] *= 3
		out[i] += carry
		carry = out[i] >> 28
		out[i] &= bottom28Bits
	}

	p256ReduceCarry(out, carry)
}

// p256Scalar4 sets out=4*out.
//
// On entry: out[0,2,...] < 2**30, out[1,3,...] < 2**29.
// On exit: out[0,2,...] < 2**30, out[1,3,...] < 2**29.
func p256Scalar4(out *feType) {
	var carry, nextCarry uint32

	for i := 0; ; i++ {
		nextCarry = out[i] >> 27
		out[i] <<= 2
		out[i] &= bottom29Bits
		out[i] += carry
		carry = nextCarry + (out[i] >> 29)
		out[i] &= bottom29Bits

		i++
		if i == p256Limbs {
			break
		}
		nextCarry = out[i] >> 26
		out[i] <<= 2
		out[i] &= bottom28Bits
		out[i] += carry
		carry = nextCarry + (out[i] >> 28)
		out[i] &= bottom28Bits
	}

	p256ReduceCarry(out, carry)
}

// p256Scalar8 sets out=8*out.
//
// On entry: out[0,2,...] < 2**30, out[1,3,...] < 2**29.
// On exit: out[0,2,...] < 2**30, out[1,3,...] < 2**29.
func p256Scalar8(out *feType) {
	var carry, nextCarry uint32

	for i := 0; ; i++ {
		nextCarry = out[i] >> 26
		out[i] <<= 3
		out[i] &= bottom29Bits
		out[i] += carry
		carry = nextCarry + (out[i] >> 29)
		out[i] &= bottom29Bits

		i++
		if i == p256Limbs {
			break
		}
		nextCarry = out[i] >> 25
		out[i] <<= 3
		out[i] &= bottom28Bits
		out[i] += carry
		carry = nextCarry + (out[i] >> 28)
		out[i] &= bottom28Bits
	}

	p256ReduceCarry(out, carry)
}

// Group operations:
//
// Elements of the elliptic curve group are represented in Jacobian
// coordinates: (x, y, z). An affine point (x', y') is x'=x/z**2, y'=y/z**3 in
// Jacobian form.

// p256PointDouble sets {xOut,yOut,zOut} = 2*{x,y,z}.
//
// See https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
func p256PointDouble(xOut, yOut, zOut, x, y, z *feType) {
	var delta, gamma, alpha, beta, tmp, tmp2 feType

	p256Square(&delta, z)
	p256Square(&gamma, y)
	p256Mult(&beta, x, &gamma)

	p256Sum(&tmp, x, &delta)
	p256Diff(&tmp2, x, &delta)
	p256Mult(&alpha, &tmp, &tmp2)
	p256Scalar3(&alpha)

	p256Sum(&tmp, y, z)
	p256Square(&tmp, &tmp)
	p256Diff(&tmp, &tmp, &gamma)
	p256Diff(zOut, &tmp, &delta)

	p256Scalar4(&beta)
	p256Square(xOut, &alpha)
	p256Diff(xOut, xOut, &beta)
	p256Diff(xOut, xOut, &beta)

	p256Diff(&tmp, &beta, xOut)
	p256Mult(&tmp, &alpha, &tmp)
	p256Square(&tmp2, &gamma)
	p256Scalar8(&tmp2)
	p256Diff(yOut, &tmp, &tmp2)
}

// p256PointAddMixed sets {xOut,yOut,zOut} = {x1,y1,z1} + {x2,y2,1}.
// (i.e. the second point is affine.)
//
// See https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
//
// Note that this function does not handle P+P, infinity+P nor P+infinity
// correctly.
func p256PointAddMixed(xOut, yOut, zOut, x1, y1, z1, x2, y2 *feType) {
	var z1z1, z1z1z1, s2, u2, h, i, j, r, rr, v, tmp feType

	p256Square(&z1z1, z1)
	p256Sum(&tmp, z1, z1)

	p256Mult(&u2, x2, &z1z1)
	p256Mult(&z1z1z1, z1, &z1z1)
	p256Mult(&s2, y2, &z1z1z1)
	p256Diff(&h, &u2, x1)
	p256Sum(&i, &h, &h)
	p256Square(&i, &i)
	p256Mult(&j, &h, &i)
	p256Diff(&r, &s2, y1)
	p256Sum(&r, &r, &r)
	p256Mult(&v, x1, &i)

	p256Mult(zOut, &tmp, &h)
	p256Square(&rr, &r)
	p256Diff(xOut, &rr, &j)
	p256Diff(xOut, xOut, &v)
	p256Diff(xOut, xOut, &v)

	p256Diff(&tmp, &v, xOut)
	p256Mult(yOut, &tmp, &r)
	p256Mult(&tmp, y1, &j)
	p256Diff(yOut, yOut, &tmp)
	p256Diff(yOut, yOut, &tmp)
}

// p256PointAdd sets {xOut,yOut,zOut} = {x1,y1,z1} + {x2,y2,z2}.
//
// See https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
//
// Note that this function does not handle P+P, infinity+P nor P+infinity
// correctly.
func p256PointAdd(xOut, yOut, zOut, x1, y1, z1, x2, y2, z2 *feType) {
	var z1z1, z1z1z1, z2z2, z2z2z2, s1, s2, u1, u2, h, i, j, r, rr, v, tmp feType

	p256Square(&z1z1, z1)
	p256Square(&z2z2, z2)
	p256Mult(&u1, x1, &z2z2)

	p256Sum(&tmp, z1, z2)
	p256Square(&tmp, &tmp)
	p256Diff(&tmp, &tmp, &z1z1)
	p256Diff(&tmp, &tmp, &z2z2)

	p256Mult(&z2z2z2, z2, &z2z2)
	p256Mult(&s1, y1, &z2z2z2)

	p256Mult(&u2, x2, &z1z1)
	p256Mult(&z1z1z1, z1, &z1z1)
	p256Mult(&s2, y2, &z1z1z1)
	p256Diff(&h, &u2, &u1)
	p256Sum(&i, &h, &h)
	p256Square(&i, &i)
	p256Mult(&j, &h, &i)
	p256Diff(&r, &s2, &s1)
	p256Sum(&r, &r, &r)
	p256Mult(&v, &u1, &i)

	p256Mult(zOut, &tmp, &h)
	p256Square(&rr, &r)
	p256Diff(xOut, &rr, &j)
	p256Diff(xOut, xOut, &v)
	p256Diff(xOut, xOut, &v)

	p256Diff(&tmp, &v, xOut)
	p256Mult(yOut, &tmp, &r)
	p256Mult(&tmp, &s1, &j)
	p256Diff(yOut, yOut, &tmp)
	p256Diff(yOut, yOut, &tmp)
}

// p256CopyConditional sets out=in if mask = 0xffffffff in constant time.
//
// On entry: mask is either 0 or 0xffffffff.
func p256CopyConditional(out, in *feType, mask uint32) {
	for i := 0; i < p256Limbs; i++ {
		tmp := mask & (in[i] ^ out[i])
		out[i] ^= tmp
	}
}

// p256SelectAffinePoint sets {out_x,out_y} to the index'th entry of table.
// On entry: index < 16, table[0] must be zero.
func p256SelectAffinePoint(xOut, yOut *feType, table []uint32, index uint32) {
	for i := range xOut {
		xOut[i] = 0
	}
	for i := range yOut {
		yOut[i] = 0
	}

	for i := uint32(1); i < 16; i++ {
		mask := i ^ index
		mask |= mask >> 2
		mask |= mask >> 1
		mask &= 1
		mask--
		for j := range xOut {
			xOut[j] |= table[0] & mask
			table = table[1:]
		}
		for j := range yOut {
			yOut[j] |= table[0] & mask
			table = table[1:]
		}
	}
}

// p256SelectJacobianPoint sets {out_x,out_y,out_z} to the index'th entry of
// table.
// On entry: index < 16, table[0] must be zero.
func p256SelectJacobianPoint(xOut, yOut, zOut *feType, table *[16][3]feType, index uint32) {
	for i := range xOut {
		xOut[i] = 0
	}
	for i := range yOut {
		yOut[i] = 0
	}
	for i := range zOut {
		zOut[i] = 0
	}

	// The implicit value at index 0 is all zero. We don't need to perform that
	// iteration of the loop because we already set out_* to zero.
	for i := uint32(1); i < 16; i++ {
		mask := i ^ index
		mask |= mask >> 2
		mask |= mask >> 1
		mask &= 1
		mask--
		for j := range xOut {
			xOut[j] |= table[i][0][j] & mask
		}
		for j := range yOut {
			yOut[j] |= table[i][1][j] & mask
		}
		for j := range zOut {
			zOut[j] |= table[i][2][j] & mask
		}
	}
}

// p256GetBit returns the bit'th bit of scalar.
func p256GetBit(scalar *[32]uint8, bit uint) uint32 {
	return uint32(((scalar[bit>>3]) >> (bit & 7)) & 1)
}

// p256ScalarBaseMult sets {xOut,yOut,zOut} = scalar*G where scalar is a
// little-endian number. Note that the value of scalar must be less than the
// order of the group.
func p256ScalarBaseMult(xOut, yOut, zOut *feType, scalar *[32]uint8) {
	nIsInfinityMask := ^uint32(0)
	var pIsNoninfiniteMask, mask, tableOffset uint32
	var px, py, tx, ty, tz feType

	for i := range xOut {
		xOut[i] = 0
	}
	for i := range yOut {
		yOut[i] = 0
	}
	for i := range zOut {
		zOut[i] = 0
	}

	// The loop adds bits at positions 0, 64, 128 and 192, followed by
	// positions 32,96,160 and 224 and does this 32 times.
	for i := uint(0); i < 32; i++ {
		if i != 0 {
			p256PointDouble(xOut, yOut, zOut, xOut, yOut, zOut)
		}
		tableOffset = 0
		for j := uint(0); j <= 32; j += 32 {
			bit0 := p256GetBit(scalar, 31-i+j)
			bit1 := p256GetBit(scalar, 95-i+j)
			bit2 := p256GetBit(scalar, 159-i+j)
			bit3 := p256GetBit(scalar, 223-i+j)
			index := bit0 | (bit1 << 1) | (bit2 << 2) | (bit3 << 3)

			p256SelectAffinePoint(&px, &py, p256Precomputed[tableOffset:], index)
			tableOffset += 30 * p256Limbs

			// Since scalar is less than the order of the group, we know that
			// {xOut,yOut,zOut} != {px,py,1}, unless both are zero, which we handle
			// below.
			p256PointAddMixed(&tx, &ty, &tz, xOut, yOut, zOut, &px, &py)
			// The result of pointAddMixed is incorrect if {xOut,yOut,zOut} is zero
			// (a.k.a.  the point at infinity). We handle that situation by
			// copying the point from the table.
			p256CopyConditional(xOut, &px, nIsInfinityMask)
			p256CopyConditional(yOut, &py, nIsInfinityMask)
			p256CopyConditional(zOut, &p256One, nIsInfinityMask)

			// Equally, the result is also wrong if the point from the table is
			// zero, which happens when the index is zero. We handle that by
			// only copying from {tx,ty,tz} to {xOut,yOut,zOut} if index != 0.
			pIsNoninfiniteMask = nonZeroToAllOnes(index)
			mask = pIsNoninfiniteMask & ^nIsInfinityMask
			p256CopyConditional(xOut, &tx, mask)
			p256CopyConditional(yOut, &ty, mask)
			p256CopyConditional(zOut, &tz, mask)
			// If p was not zero, then n is now non-zero.
			nIsInfinityMask &^= pIsNoninfiniteMask
		}
	}
}

// p256PointToAffine converts a Jacobian point to an affine point. If the input
// is the point at infinity then it returns (0, 0) in constant time.
func p256PointToAffine(xOut, yOut, x, y, z *feType) {
	var zInv, zInvSq feType

	p256Invert(&zInv, z)
	p256Square(&zInvSq, &zInv)
	p256Mult(xOut, x, &zInvSq)
	p256Mult(&zInv, &zInv, &zInvSq)
	p256Mult(yOut, y, &zInv)
}

// p256ToAffine returns a pair of *big.Int containing the affine representation
// of {x,y,z}.
func p256ToAffine(x, y, z *feType) (xOut, yOut *big.Int) {
	var xx, yy feType
	p256PointToAffine(&xx, &yy, x, y, z)
	return p256ToBig(&xx), p256ToBig(&yy)
}

// p256ScalarMult sets {xOut,yOut,zOut} = scalar*{x,y}.
func p256ScalarMult(xOut, yOut, zOut, x, y *feType, scalar *[32]uint8) {
	var px, py, pz, tx, ty, tz feType
	var precomp [16][3]feType
	var nIsInfinityMask, index, pIsNoninfiniteMask, mask uint32

	// We precompute 0,1,2,... times {x,y}.
	precomp[1][0] = *x
	precomp[1][1] = *y
	precomp[1][2] = p256One

	for i := 2; i < 16; i += 2 {
		p256PointDouble(&precomp[i][0], &precomp[i][1], &precomp[i][2], &precomp[i/2][0], &precomp[i/2][1], &precomp[i/2][2])
		p256PointAddMixed(&precomp[i+1][0], &precomp[i+1][1], &precomp[i+1][2], &precomp[i][0], &precomp[i][1], &precomp[i][2], x, y)
	}

	for i := range xOut {
		xOut[i] = 0
	}
	for i := range yOut {
		yOut[i] = 0
	}
	for i := range zOut {
		zOut[i] = 0
	}
	nIsInfinityMask = ^uint32(0)

	// We add in a window of four bits each iteration and do this 64 times.
	for i := 0; i < 64; i++ {
		if i != 0 {
			p256PointDouble(xOut, yOut, zOut, xOut, yOut, zOut)
			p256PointDouble(xOut, yOut, zOut, xOut, yOut, zOut)
			p256PointDouble(xOut, yOut, zOut, xOut, yOut, zOut)
			p256PointDouble(xOut, yOut, zOut, xOut, yOut, zOut)
		}

		index = uint32(scalar[31-i/2])
		if (i & 1) == 1 {
			index &= 15
		} else {
			index >>= 4
		}

		// See the comments in scalarBaseMult about handling infinities.
		p256SelectJacobianPoint(&px, &py, &pz, &precomp, index)
		p256PointAdd(&tx, &ty, &tz, xOut, yOut, zOut, &px, &py, &pz)
		p256CopyConditional(xOut, &px, nIsInfinityMask)
		p256CopyConditional(yOut, &py, nIsInfinityMask)
		p256CopyConditional(zOut, &pz, nIsInfinityMask)

		pIsNoninfiniteMask = nonZeroToAllOnes(index)
		mask = pIsNoninfiniteMask & ^nIsInfinityMask
		p256CopyConditional(xOut, &tx, mask)
		p256CopyConditional(yOut, &ty, mask)
		p256CopyConditional(zOut, &tz, mask)
		nIsInfinityMask &^= pIsNoninfiniteMask
	}
}

// p256FromBigRaw sets out = in.
func p256FromBigRaw(out *feType, in *big.Int) {
	var tmp [4]big.Word
	copy(tmp[:], in.Bits())

	out[0] = uint32(tmp[0]) & bottom29Bits
	out[1] = uint32(tmp[0]>>29) & bottom28Bits
	out[2] = uint32(tmp[0]>>57) | uint32(tmp[1]<<7)
	out[2] &= bottom29Bits
	out[3] = uint32(tmp[1]>>22) & bottom28Bits
	out[4] = uint32(tmp[1]>>50) | uint32(tmp[2]<<14)
	out[4] &= bottom29Bits
	out[5] = uint32(tmp[2]>>15) & bottom28Bits
	out[6] = uint32(tmp[2]>>43) | uint32(tmp[3]<<21)
	out[6] &= bottom29Bits
	out[7] = uint32(tmp[3]>>8) & bottom28Bits
	out[8] = uint32(tmp[3] >> 36)
}

// p256FromBig sets out = R*in.
func p256FromBig(out *feType, in *big.Int) {
	p256FromBigRaw(out, in)
	p256Mult(out, out, &p256FeRR)
	p256NormalWeak(out)
}

// p256ToBigRaw returns a *big.Int containing the value of in raw form.
func p256ToBigRaw(in *feType) *big.Int {
	var bits [5]big.Word
	bits[0] = big.Word(in[0])
	bits[0] |= big.Word(in[1]) << 29
	bits[0] |= big.Word(in[2]) << 57
	bits[1] = big.Word(in[2] >> 7)
	bits[1] |= big.Word(in[3]) << 22
	bits[1] |= big.Word(in[4]) << 50
	bits[2] = big.Word(in[4] >> 14)
	bits[2] |= big.Word(in[5]) << 15
	bits[2] |= big.Word(in[6]) << 43
	bits[3] = big.Word(in[6] >> 21)
	bits[3] |= big.Word(in[7]) << 8
	bits[3] |= big.Word(in[8]) << 36
	bits[4] = big.Word(in[8] >> 28)
	result := new(big.Int).SetBits(bits[:])
	return result
}

// p256ToBig returns a *big.Int containing the value of in montgomery form.
func p256ToBig(in *feType) *big.Int {
	var btmp [p256LimbsProd]uint32
	copy(btmp[:], in[:])
	p256Reduction(&btmp)
	var tmp feType
	var tmp1 feType
	p256ProdHigh(&tmp, &btmp)
	p256NormalWeak(&tmp)
	p256Mod(&tmp1, &tmp)
	return p256ToBigRaw(&tmp1)
}
