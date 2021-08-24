// Copyright 2013 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build !amd64,!arm64	field

package sm2

import (
	"crypto/rand"
	"gitee.com/jkuang/go-fastecdsa"
	"math/big"
	"testing"
)

func init() {
	SM2go()
}

func feMontMult(x, y *big.Int) *big.Int {
	var xp, yp, res feType
	p256FromBig(&xp, x)
	p256FromBig(&yp, y)
	p256Mult(&res, &xp, &yp)
	return p256ToBig(&res)
}

func feMontSqr(x *big.Int) *big.Int {
	var xp, res feType
	p256FromBig(&xp, x)
	p256Square(&res, &xp)
	return p256ToBig(&res)
}

func feSqrt(y *big.Int) (rt *big.Int) {
	var yp, res feType
	p256FromBig(&yp, y)
	if p256Sqrt(&res, &yp) {
		rt = p256ToBig(&res)
	}
	return
}

func bigExp(x *big.Int, y int64) *big.Int {
	return new(big.Int).Exp(x, big.NewInt(y), sm2Params.P)
}

func feMontRed(y *big.Int) *big.Int {
	var yp, rr feType
	var res [p256LimbsProd]uint32
	p256FromBigRaw(&yp, y)
	copy(res[:], yp[:])
	p256Reduction(&res)
	p256ProdHigh(&yp, &res)
	p256NormalWeak(&yp)
	p256Mod(&rr, &yp)
	return p256ToBigRaw(&rr)
}

func TestFeMontRed(t *testing.T) {
	montOne := feType{2, 0, (1 << 29) - (1 << 8), (1 << 11) - 1, 0, 0, 0, 1 << 25}
	if montOne != p256One {
		t.Logf("p256One diff montOne:\n%v\n%v", p256One, montOne)
		t.Fail()
	}
	bOne := p256ToBig(&p256One)
	if bOne.Cmp(bigOne) != 0 {
		t.Logf("p256One to Big not bigOne: len(%d) %s", bOne.BitLen(), bOne.Text(16))
		t.Fail()
	}
}

func TestExpQuad(t *testing.T) {
	quadP := new(big.Int).Rsh(sm2Params.P, 2)
	var tmp, res feType
	p256FromBig(&tmp, x1)
	p256ExpQuadP(&res, &tmp)
	eQuad := new(big.Int).Exp(x1, quadP, sm2Params.P)
	re := p256ToBig(&res)
	if eQuad.Cmp(re) != 0 {
		t.Logf("ExpQuad diff:\n%s vs\n%s", eQuad.Text(16), re.Text(16))
		t.Fail()
	}
}

func TestFeMontMul(t *testing.T) {
	prod := new(big.Int).Mul(x1, y1)
	m1 := new(big.Int).Mod(prod, sm2Params.P)
	m2 := feMontMult(x1, y1)
	if m1.Cmp(m2) != 0 {
		t.Logf("MontMulMod step 1 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
	prod.Mul(x1, x1)
	prod.Mul(prod, y1)
	m1.Mod(prod, sm2Params.P)
	m2 = feMontSqr(x1)
	m2 = feMontMult(m2, y1)
	if m1.Cmp(m2) != 0 {
		t.Logf("MontMulMod step 1a diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
	prod.Mul(x2, y2)
	m1.Mod(prod, sm2Params.P)
	m2 = feMontMult(x2, y2)
	if m1.Cmp(m2) != 0 {
		t.Logf("MontMulMod step2 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
}

func TestExpPrecomput(t *testing.T) {
	var table [5]feType
	var tmp feType
	//xx := x2 // works ok
	xx := x1
	p256FromBig(&tmp, xx)
	precomputeExp(&table, &tmp)
	tt := bigExp(xx, 3)
	res := p256ToBig(&table[0])
	if tt.Cmp(res) != 0 {
		t.Logf("Exp e2 diff:\n%s vs \n%s", tt.Text(16), res.Text(16))
		t.Fail()
	}
	tt = bigExp(xx, 15)
	res = p256ToBig(&table[1])
	if tt.Cmp(res) != 0 {
		t.Logf("Exp e4 diff:\n%s vs \n%s", tt.Text(16), res.Text(16))
		t.Fail()
	}
	tt = bigExp(xx, 255)
	res = p256ToBig(&table[2])
	if tt.Cmp(res) != 0 {
		t.Logf("Exp e8 diff:\n%s vs \n%s", tt.Text(16), res.Text(16))
		t.Fail()
	}
	tt = bigExp(xx, 0xffff)
	res = p256ToBig(&table[3])
	if tt.Cmp(res) != 0 {
		t.Logf("Exp e16 diff:\n%s vs \n%s", tt.Text(16), res.Text(16))
		t.Fail()
	}
	tt = bigExp(xx, 0xffffffff)
	res = p256ToBig(&table[4])
	if tt.Cmp(res) != 0 {
		t.Logf("Exp e32 diff:\n%s vs \n%s", tt.Text(16), res.Text(16))
		t.Fail()
	}
}

func TestFeMontSqr(t *testing.T) {
	prod := new(big.Int).Mul(x1, x1)
	m1 := new(big.Int).Mod(prod, sm2Params.P)
	m2 := feMontSqr(x1)
	if m1.Cmp(m2) != 0 {
		t.Logf("MontSqrMod step 1 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
	if m3 := feSqrt(m2); m3 == nil {
		t.Log("Can't get Sqrt")
		t.Fail()
	} else if m3.Cmp(x1) != 0 {
		t.Logf("ModSqrt step 1 diff:\n%s vs\n%s", x1.Text(16), m3.Text(16))
		t.Fail()
	}
	prod = new(big.Int).Mul(x2, x2)
	m1 = new(big.Int).Mod(prod, sm2Params.P)
	m2 = feMontSqr(x2)
	if m1.Cmp(m2) != 0 {
		t.Logf("MontSqrMod step2 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
	if m3 := feSqrt(m2); m3 == nil {
		t.Log("Can't get Sqrt")
		t.Fail()
	} else if m3.Cmp(x2) != 0 {
		m3.Sub(sm2Params.P, m3)
		if m3.Cmp(x2) != 0 {
			t.Logf("ModSqrt step 2 diff:\n%s vs\n%s", x1.Text(16), m3.Text(16))
			t.Fail()
		}
	}
}

func TestFeInverse(t *testing.T) {
	p := sm2Params.P
	RR := x1
	Rinv := new(big.Int).ModInverse(RR, p)
	if rMul := feMontMult(Rinv, RR); rMul.Cmp(bigOne) != 0 {
		t.Logf("RR * Rinv != 1, rMul: %s", rMul.Text(16))
		t.Fail()
	}
	var res, yy feType
	p256FromBig(&yy, RR)
	p256Invert(&res, &yy)
	RinvA := p256ToBig(&res)
	if rMul := feMontMult(RinvA, RR); rMul.Cmp(bigOne) != 0 {
		t.Logf("RR * RinvA != 1, step1 rMul: %s", rMul.Text(16))
		t.Fail()
	}
	if Rinv.Cmp(RinvA) != 0 {
		t.Logf("p256Inverse step1 diff:\n%s vs\n%s", Rinv.Text(16), RinvA.Text(16))
		t.Fail()
	}
	RR = x2
	Rinv.ModInverse(RR, p)
	if rMul := feMontMult(Rinv, RR); rMul.Cmp(bigOne) != 0 {
		t.Logf("RR * Rinv != 1, rMul: %s", rMul.Text(16))
		t.Fail()
	}
	p256FromBig(&yy, RR)
	p256Invert(&res, &yy)
	RinvA = p256ToBig(&res)
	if rMul := feMontMult(RinvA, RR); rMul.Cmp(bigOne) != 0 {
		t.Logf("RR * RinvA != 1, step1 rMul: %s", rMul.Text(16))
		t.Fail()
	}
	if Rinv.Cmp(RinvA) != 0 {
		t.Logf("p256Inverse step1 diff:\n%s vs\n%s", Rinv.Text(16), RinvA.Text(16))
		t.Fail()
	}
}

func TestPointRecover(t *testing.T) {
	c := pSM2
	px, py := c.ScalarBaseMult(d1.Bytes())
	v := py.Bit(0)
	if py2, err := RecoverPoint(px, v); err != nil {
		t.Log("Can't recover step1 pointY, error:", err)
		t.Fail()
	} else if py2.Cmp(py) != 0 {
		t.Logf("RecoverPoint step1 diff:\n%s vs\n%s", py.Text(16), py2.Text(16))
		t.Fail()
	}
	px, py = c.ScalarBaseMult(d2.Bytes())
	v = py.Bit(0)
	if py2, err := RecoverPoint(px, v); err != nil {
		t.Log("Can't recover step2 pointY, error:", err)
		t.Fail()
	} else if py2.Cmp(py) != 0 {
		t.Logf("RecoverPoint step2 diff:\n%s vs\n%s", py.Text(16), py2.Text(16))
		t.Fail()
	}
}

func BenchmarkFeInverse(b *testing.B) {
	priv, _ := fastecdsa.GenerateKey(P256(), rand.Reader)
	var res, yy feType
	p256FromBig(&yy, priv.PublicKey.X)

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			p256Invert(&res, &yy)
		}
	})
}

func BenchmarkFeMontModMul(b *testing.B) {
	//c := pSM2
	c := sm2g
	var xp, rrP, res feType
	p256FromBigRaw(&xp, x1)
	p256FromBigRaw(&rrP, c.rr)

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			p256Mult(&res, &xp, &rrP)
		}
	})
}

func BenchmarkFeMontSqr(b *testing.B) {
	var res feType
	p256FromBig(&res, x1)
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			p256Square(&res, &res)
		}
	})
}

func BenchmarkPointRecover(b *testing.B) {
	c := pSM2
	px, py := c.ScalarBaseMult(d1.Bytes())
	v := py.Bit(0)
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			RecoverPoint(px, v)
		}
	})
}

func BenchmarkFeECADD(b *testing.B) {
	curve := SM2()

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = curve.Add(x1, y1, x2, y2)
		}
	})
}

func BenchmarkFeECDBL(b *testing.B) {
	curve := SM2()

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = curve.Double(x1, y1)
		}
	})
}

func BenchmarkFeECMULT(b *testing.B) {
	Curve := SM2()
	goGx := Curve.Params().Gx
	goGy := Curve.Params().Gy

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = Curve.ScalarMult(goGx, goGy, d1.Bytes())
		}
	})
}

func BenchmarkFeECGMULT(b *testing.B) {
	Curve := SM2()

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = Curve.ScalarBaseMult(d1.Bytes())
		}
	})
}
