package ecc

// #cgo CXXFLAGS: -O3 -Wpedantic -I../include -Wno-uninitialized -std=c++11
// #include "ecc.h"
import "C"

// for gcc 4.8.5, __builtin_add/sub_overflow improve 10% -- 20%
// -fvar-tracking-assignments
// WITH_SM2_MULTP using shift instead of multiply, no improvement @x86_64/arm64
// cgo CXXFLAGS: -O2 -Wpedantic -Wall -std=gnu++11 -DWITH_SM2_MULTP

import (
	"crypto/elliptic"
	"errors"
	"io"
	"math/big"
	"unsafe"
)

// A Curve represents a short-form Weierstrass curve with a=-3.
// See https://www.hyperelliptic.org/EFD/g1p/auto-shortw.html
type Curve = elliptic.Curve

// CurveParams contains the parameters of an elliptic curve and also provides
// a generic, non-constant time implementation of Curve.
type CurveParams = elliptic.CurveParams

// eccCurve contains the parameters of an elliptic curve and also provides
// a generic, non-constant time implementation of Curve.
type eccCurve struct {
	*CurveParams
	hnd    C.CURVE_HND
	inited bool
}

var errParam = errors.New("error parameter")
var sm2c eccCurve
var bigOne *big.Int

func init() {
	if sm2c.inited {
		return
	}
	sm2c.inited = true
	bigOne = big.NewInt(1)
	cHnd := C.get_curve(C.ECC_CURVE_SM2)
	if cHnd == C.CURVE_HND(uintptr(0)) {
		return
	}
	sm2c.hnd = cHnd
	sm2c.CurveParams = getCurveParams(0)
}

func SM2C() *eccCurve {
	if sm2c.CurveParams == nil {
		return nil
	}
	return &sm2c
}

// Marshal converts a point into the uncompressed form specified in section 4.3.6 of ANSI X9.62.
func Marshal(curve Curve, x, y *big.Int) []byte {
	byteLen := (curve.Params().BitSize + 7) >> 3

	ret := make([]byte, 1+2*byteLen)
	ret[0] = 4 // uncompressed point

	xBytes := x.Bytes()
	copy(ret[1+byteLen-len(xBytes):], xBytes)
	yBytes := y.Bytes()
	copy(ret[1+2*byteLen-len(yBytes):], yBytes)
	return ret
}

// Unmarshal converts a point, serialized by Marshal, into an x, y pair.
// It is an error if the point is not in uncompressed form or is not on the curve.
// On error, x = nil.
func Unmarshal(curve Curve, data []byte) (x, y *big.Int) {
	//byteLen := (curve.Params().BitSize + 7) >> 3
	const byteLen = 32
	if len(data) == 1+byteLen {
		if data[0] != 0x2 && data[0] != 0x3 {
			return
		}
		p := curve.Params().P
		x = new(big.Int).SetBytes(data[1:])
		if x.Cmp(p) >= 0 {
			return nil, nil
		}
		if y1, err := RecoverPoint(x, uint(data[0]&1)); err != nil {
			return nil, nil
		} else {
			y = y1
		}
		return
	}
	if len(data) != 1+2*byteLen {
		return
	}
	if data[0] != 4 { // uncompressed form
		return
	}
	p := curve.Params().P
	x = new(big.Int).SetBytes(data[1 : 1+byteLen])
	y = new(big.Int).SetBytes(data[1+byteLen:])
	if x.Cmp(p) >= 0 || y.Cmp(p) >= 0 {
		return nil, nil
	}
	if !curve.IsOnCurve(x, y) {
		return nil, nil
	}
	return
}

func vliModMultMontP(x, y []big.Word) *big.Int {
	var r [4]big.Word
	C.mont_sm2_mod_mult_p((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&x[0])), (*C.u64)(unsafe.Pointer(&y[0])))
	return new(big.Int).SetBits(r[:4])
}

func vliModMultMontN(x, y []big.Word) *big.Int {
	var r [4]big.Word
	C.mont_sm2_mod_mult_n((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&x[0])), (*C.u64)(unsafe.Pointer(&y[0])))
	return new(big.Int).SetBits(r[:4])
}

func RecoverPoint(x1 *big.Int, v uint) (y1 *big.Int, err error) {
	cHnd := sm2c.hnd
	var r [4]big.Word
	x := fromWordSlice(x1.Bits())
	if C.point_recoverY((*C.u64)(unsafe.Pointer(&r[0])), (*C.u64)(unsafe.Pointer(x)),
		C.uint(v), cHnd) != 0 {
		y1 = new(big.Int).SetBits(r[:])
	} else {
		// set err
		return nil, errParam
	}
	return
}

func getCurveParams(curveId uint) *CurveParams {
	var cHnd C.CURVE_HND
	if curveId == 0 {
		cHnd = sm2c.hnd
	} else {
		cHnd = C.get_curve(C.uint(curveId))
	}
	if cHnd == C.CURVE_HND(uintptr(0)) {
		return nil
	}
	var p, n, b, gx, gy [4]big.Word
	C.get_curve_params((*C.u64)(unsafe.Pointer(&p[0])),
		(*C.u64)(unsafe.Pointer(&n[0])), (*C.u64)(unsafe.Pointer(&b[0])),
		(*C.u64)(unsafe.Pointer(&gx[0])), (*C.u64)(unsafe.Pointer(&gy[0])),
		cHnd)
	sm2Params := &CurveParams{Name: "SM2-c"}
	sm2Params.P = new(big.Int).SetBits(p[:])
	sm2Params.N = new(big.Int).SetBits(n[:])
	sm2Params.B = new(big.Int).SetBits(b[:])
	sm2Params.Gx = new(big.Int).SetBits(gx[:])
	sm2Params.Gy = new(big.Int).SetBits(gy[:])
	sm2Params.BitSize = 256
	return sm2Params
}

func (c eccCurve) newPoint(x, y, z *big.Int) *C.Point {
	var pt C.Point
	if z == nil {
		z = bigOne
	}
	pt.x = *fromWordSlice(x.Bits())
	pt.y = *fromWordSlice(y.Bits())
	pt.z = *fromWordSlice(z.Bits())
	return &pt
}

func (c eccCurve) Inverse(k *big.Int) *big.Int {
	bnInv := vliModInv(k.Bits(), c.Params().N.Bits())
	return new(big.Int).SetBits(bnInv)
}

func (c eccCurve) Add(x1, y1, x2, y2 *big.Int) (rx, ry *big.Int) {
	pt1 := c.newPoint(x1, y1, nil)
	pt2 := c.newPoint(x2, y2, nil)
	var pt C.Point
	C.point_add(&pt, pt1, pt2, c.hnd)
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	return
}

func (c eccCurve) addJacobian(x1, y1, z1, x2, y2, z2 *big.Int) (rx, ry, rz *big.Int) {
	pt1 := c.newPoint(x1, y1, z1)
	pt2 := c.newPoint(x2, y2, z2)
	var pt C.Point
	C.point_add_jacobian(&pt, pt1, pt2, c.hnd)
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	rz = new(big.Int).SetBits(toWordSlice(pt.z))
	return
}

func (c eccCurve) doubleJacobian(x, y, z *big.Int) (rx, ry, rz *big.Int) {
	pt1 := c.newPoint(x, y, z)
	var pt C.Point
	C.point_double_jacobian(&pt, pt1, c.hnd)
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	rz = new(big.Int).SetBits(toWordSlice(pt.z))
	return
}

func (c eccCurve) Double(x, y *big.Int) (rx, ry *big.Int) {
	pt1 := c.newPoint(x, y, nil)
	var pt C.Point
	C.point_double(&pt, pt1, c.hnd)
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	return
}

func (c eccCurve) ScalarBaseMult(k []byte) (rx, ry *big.Int) {
	var pt C.Point
	scal := new(big.Int).SetBytes(k)
	var ss [4]big.Word
	copy(ss[:], scal.Bits())
	C.point_cmult(&pt, nil, nil, (*C.u64)(unsafe.Pointer(&ss[0])), c.hnd, nil)
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	return
}

func (c eccCurve) ScalarMult(x, y *big.Int, k []byte) (rx, ry *big.Int) {
	pt1 := c.newPoint(x, y, nil)
	var pt C.Point
	scal := new(big.Int).SetBytes(k)
	var ss [4]big.Word
	copy(ss[:], scal.Bits())
	C.point_mult(&pt, pt1, (*C.u64)(unsafe.Pointer(&ss[0])), c.hnd, nil)
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	return
}

func (c eccCurve) cMult(x, y *big.Int, k, gk, sBuff []byte) (rx, ry *big.Int) {
	pt1 := c.newPoint(x, y, nil)
	var pt C.Point
	scal := new(big.Int).SetBytes(k)
	var ss [4]big.Word
	copy(ss[:], scal.Bits())
	scal.SetBytes(gk)
	var gs [4]big.Word
	copy(gs[:], scal.Bits())
	if sBuff != nil {
		if len(sBuff) < 2000 {
			sBuff = make([]byte, 2048)
		}
		C.point_cmult(&pt, pt1, (*C.u64)(unsafe.Pointer(&ss[0])),
			(*C.u64)(unsafe.Pointer(&gs[0])), c.hnd, unsafe.Pointer(&sBuff[0]))
	} else {
		C.point_cmult(&pt, pt1, (*C.u64)(unsafe.Pointer(&ss[0])),
			(*C.u64)(unsafe.Pointer(&gs[0])), c.hnd, nil)
	}
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	return
}

func (c eccCurve) CombinedMult(x, y *big.Int, k, gk []byte) (rx, ry *big.Int) {
	return c.cMult(x, y, k, gk, nil)
}

func (c eccCurve) affineFromJacobian(x, y, z *big.Int) (xOut, yOut *big.Int) {
	var xb, yb [4]big.Word
	pt := c.newPoint(x, y, z)
	C.affine_from_jacobian((*C.u64)(unsafe.Pointer(&xb[0])),
		(*C.u64)(unsafe.Pointer(&yb[0])), pt, c.hnd)
	xOut = new(big.Int).SetBits(xb[:])
	yOut = new(big.Int).SetBits(yb[:])
	return
}

func (c eccCurve) Verify(rB, sB, msgB, px, py *big.Int) bool {
	var r, s, msg [4]big.Word
	copy(r[:], rB.Bits())
	copy(s[:], sB.Bits())
	copy(msg[:], msgB.Bits())
	pt := c.newPoint(px, py, nil)
	return C.ecc_verify((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&s[0])),
		(*C.u64)(unsafe.Pointer(&msg[0])), pt, c.hnd) != 0
}

func (c eccCurve) Sign(rand io.Reader, msgB, secret *big.Int) (r, s *big.Int, v uint, err error) {
	var rw, sw, msg [4]big.Word
	copy(msg[:], msgB.Bits())
	pt := c.newPoint(bigOne, bigOne, secret)
	v = uint(C.ecc_sign((*C.u64)(unsafe.Pointer(&rw[0])),
		(*C.u64)(unsafe.Pointer(&sw[0])),
		(*C.u64)(unsafe.Pointer(&msg[0])), pt, c.hnd))
	r = new(big.Int).SetBits(rw[:])
	s = new(big.Int).SetBits(sw[:])
	return
}

func (c eccCurve) Recover(rB, sB, msgB *big.Int, v uint) (pubX, pubY *big.Int, err error) {
	var r, s, msg [4]big.Word
	copy(r[:], rB.Bits())
	copy(s[:], sB.Bits())
	copy(msg[:], msgB.Bits())
	var pt C.Point
	ret := C.ecc_recover((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&s[0])),
		(*C.u64)(unsafe.Pointer(&msg[0])), C.uint(v), &pt, c.hnd)
	if ret != 0 {
		pubX = new(big.Int).SetBits(toWordSlice(pt.x))
		pubY = new(big.Int).SetBits(toWordSlice(pt.y))
	} else {
		return nil, nil, errParam
	}
	return
}
