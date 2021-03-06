/*
 * Copyright(c) 2020 Jesse Kuang  <jkuang@21cn.com>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *  * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#pragma once
#ifndef __CURVE_IMPL_HPP__
#define __CURVE_IMPL_HPP__

#include <string.h>
#include "vli.hpp"
#include "vli_bn.hpp"
#include "ecc_impl.hpp"
#include "mont.hpp"

namespace ecc {

using namespace vli;

constexpr int W = 5;
constexpr int wSize = 1<<(W-1);
constexpr int BaseW = 6;	// 7 for 37 *64 pre computed Points, 6 for 43*32
constexpr int wBaseSize = 1<<(BaseW-1);


template<const uint N, typename curveT> forceinline
void pre_compute(const curveT& curve, point_t<N> pre_comp[wSize],
		const point_t<N>& p) noexcept
{
	pre_comp[0] = p;
	curve.point_double(pre_comp[1], p.x, p.y);
	for (int i = 2; i < wSize; ++i) {
		if ((i & 1) == 0) {
			curve.point_add(pre_comp[i], pre_comp[i-1], p.x, p.y);
		} else {
			curve.point_double(pre_comp[i], pre_comp[i >> 1]);
		}
	}
}

template<const uint N=4> forceinline
static constexpr int nwBaseNAF() { return (64*N -1 + BaseW) / BaseW; }


/**
 * struct ecc_curve - definition of elliptic curve
 *
 * @name:	Short name of the curve.
 * @g:		Generator point of the curve.
 * @p:		Prime number, if Barrett's reduction is used for this curve
 * @n:		Order of the curve group.
 * @a:		Curve parameter a.
 * @b:		Curve parameter b.
 */
template<const uint N=4, const bool A_is_n3=true>
class ecc_curve {
public:
	using felem_t = bignum<N>;
	typedef spoint_t<N>	gwNAF_t[wBaseSize];
	static constexpr int maxBaseNAF = nwBaseNAF<N>();
	ecc_curve(const char *_name, const u64 *_gx, const u64 *_gy, const u64 *p,
			const u64 *n, const u64 *a, const u64 *b, const u64* rrP,
			const u64* rrN, const u64 k0P, const u64 k0N):
		name(_name), gx(_gx), gy(_gy),
		_p(p), _n(n), _a(a), _b(b), rr_p(rrP), rr_n(rrN), k0_p(k0P), k0_n(k0N),
		_a_is_neg3(A_is_n3)
	{
		init();
		static_assert(N > 3, "curve only support 256Bits or more");
	}
	ecc_curve(ecc_curve &&) = default;
	static const ecc_curve* new_ecc_curve(const char *_name, const u64 *_gx,
			const u64 *_gy, const u64 *p, const u64 *n, const u64 *a,
			const u64 *b) noexcept
	{
		static_assert(N > 3, "curve only support 256Bits or more");
		bignum<N>	prr, nrr;
		u64		k0P, k0N;
		bignum<N>	prime(p);
		bignum<N>	n_prime(n);
		k0P = calcK0<N>(prime);
		calcRR<N>(prr, prime);
		k0N = calcK0<N>(n_prime);
		calcRR<N>(nrr, n_prime);
		return new ecc_curve(_name, _gx, _gy, p, n, a, b, prr.data(),
						nrr.data(), k0P, k0N);
	}
	static const ecc_curve* new_ecc_curve(const char *_name, const u8 *params)
			noexcept
	{
		u64		vp[N], va[N], vb[N], vx[N], vy[N], vn[N];
		vli_from_be64<N>(vp, params);
		vli_from_be64<N>(va, params + 8*N);
		vli_from_be64<N>(vb, params + 16*N);
		vli_from_be64<N>(vx, params + 24*N);
		vli_from_be64<N>(vy, params + 32*N);
		vli_from_be64<N>(vn, params + 40*N);
		return new_ecc_curve(_name, vx, vy, vp, vn, va, vb);
	}
	const uint ndigits() const { return N; }
	explicit operator bool() const noexcept
	{
		//return name != nullptr && name[0] != 0;
		return name != "";
	}
	void getP(u64 *v) const noexcept { _p.set(v); }
	void getN(u64 *v) const noexcept { _n.set(v); }
	void getA(u64 *v) const noexcept { _a.set(v); }
	void getB(u64 *v) const noexcept { _b.set(v); }
	void getGx(u64 *v) const noexcept { gx.set(v); }
	void getGy(u64 *v) const noexcept { gy.set(v); }
	const felem_t& getGx() const noexcept { return gx; }
	const felem_t& getGy() const noexcept { return gy; }
	const felem_t& montParamA() const noexcept { return _mont_a; }
	const felem_t& montParamB() const noexcept { return _mont_b; }
	const felem_t& paramP() const noexcept { return _p; }
	const felem_t& paramN() const noexcept { return _n; }
	const felem_t& paramA() const noexcept { return _a; }
	const felem_t& paramB() const noexcept { return _b; }
	const bool a_is_pminus3() const noexcept { return _a_is_neg3; }
	const bool a_is_zero() const noexcept { return _a_is_zero; }
	const felem_t& mont_one() const noexcept { return _mont_one; }
	const felem_t& mont_three() const noexcept { return _mont_three; }
	const felem_t& quadP() const noexcept { return _quadP; }
#ifndef	EXHAUSTIVE_TEST_ORDER
	const felem_t& P_minus_N() const noexcept { return _p_minus_n; }
#endif
#if	defined(WITH_ECDSA) && !defined(NO_HALF_N)
	const felem_t& halfN() const noexcept { return _half_n; }
#endif
	const size_t presBuffSize() const noexcept { return sizeof(presBuff); }
	void to_montgomery(felem_t& res, const u64 *x) const noexcept
	{
		felem_t   *xx = reinterpret_cast<felem_t *>(const_cast<u64 *>(x));
		res.mont_mult(*xx, rr_p, _p, k0_p);
	}
	void to_montgomery(felem_t& res, const felem_t& x) const noexcept
	{
		res.mont_mult(x, rr_p, _p, k0_p);
	}
	void to_montgomeryN(felem_t& res, const felem_t& x) const noexcept
	{
		res.mont_mult(x, rr_n, _n, k0_n);
	}
	void from_montgomery(felem_t& res, const felem_t& y) const noexcept
	{
		res.mont_reduction(y, _p, k0_p);
	}
	void from_montgomery(u64* result, const felem_t& y) const noexcept
	{
		felem_t   *res = reinterpret_cast<felem_t *>(result);
		res->mont_reduction(y, _p, k0_p);
	}
	void from_montgomeryN(felem_t& res, const felem_t& y) const noexcept
	{
		res.mont_reduction(y, _n, k0_n);
	}
	void mont_nmult(felem_t& res, const felem_t& left, const felem_t& right)
	const noexcept
	{
		res.mont_mult(left, right, _n, k0_n);
	}
	void mont_mmult(felem_t& res, const felem_t& left, const felem_t& right)
	const noexcept
	{
		res.mont_mult(left, right, _p, k0_p);
	}
	void mont_msqr(felem_t& res, const felem_t left, const uint nTimes=1)
	const noexcept
	{
		res.mont_sqr(left, _p, k0_p);
		for (uint i=1; i < nTimes; i++) res.mont_sqr(res, _p, k0_p);
	}
	// left,right less than p, result may large than p
	void mod_add(felem_t& res, const felem_t& left, const felem_t& right)
	const noexcept
	{
#ifndef	ommit
		if (res.add(left, right)) {
			res.sub_from(_p);
		}
#else
		res.mod_add(left, right, _p);
#endif
	}
	// left,right less than p, result may large than p
	void mod_add_to(felem_t& res, const felem_t& right) const noexcept
	{
#ifndef	ommit
		if (res.add_to(right)) {
			res.sub_from(_p);
		}
#else
		res.mod_add_to(right, _p);
#endif
	}
	// left,right less than p
	void mod_sub(felem_t& res, const felem_t& left, const felem_t& right)
	const noexcept
	{
		if (res.sub(left, right)) {
			res.add_to(_p);
		}
	}
	void mod_sub_from(felem_t& res, const felem_t& right) const noexcept
	{
		if (res.sub_from(right)) {
			res.add_to(_p);
		}
	}
	void mont_mult2(felem_t& res, const felem_t& left) const noexcept
	{
#ifdef	ommit
		this->mod_add(res, left, left);
#else
		if (res.lshift1(left) != 0) {
			res.sub_from(_p);
		}
#endif
	}
	void mont_mult4(felem_t& res) const noexcept
	{
		this->mont_mult2(res, res);
		this->mont_mult2(res, res);
	}
	void mont_mult8(felem_t& res) const noexcept
	{
		this->mont_mult2(res, res);
		this->mont_mult2(res, res);
		this->mont_mult2(res, res);
	}
	void mont_div2(felem_t& res) const noexcept
	{
		if ( likely(res.is_even()) ) {
			res.rshift1();
		} else {
			bool carry = res.add_to(this->_p);
			res.rshift1(carry);
		}
	}
	bool prod_equal(const felem_t& prod, const felem_t& x, const felem_t& yp)
	const noexcept
	{
		felem_t	xp;
		this->to_montgomery(xp, x);
		this->mont_mmult(xp, xp, yp);
		return prod == xp;
	}
	void
	mod_exp(felem_t& res, const felem_t &x1, const felem_t& y1) const noexcept
	{
		felem_t		tmp;
		to_montgomery(tmp, x1);
		ecc::mont_mod_exp<N>(*this, res, tmp, y1);
		from_montgomery(res, res);
	}
	// prime specially be 4 * k + 3
	bool mont_sqrt(felem_t& res, const felem_t& xp) const noexcept
	{
		felem_t	a1;
		//this->mont_mod_exp(a1, xp, this->quadP());
		ecc::mont_mod_exp<N>(*this, a1, xp, this->quadP());
		felem_t	a0;
		mont_mmult(res, a1, xp);
		mont_mmult(a0, res, a1);
		return (a0 == this->mont_one());
	}
	bool mod_sqrt(felem_t& res, const felem_t& x) const noexcept
	{
		felem_t	xp;
		to_montgomery(xp, x);
		auto ret = mont_sqrt(res, xp);
		from_montgomery(res, res);
		return ret;
	}
	bool is_on_curve(const felem_t& x1, const felem_t& y1) const noexcept
	{
		felem_t	xp, yp, tt, yy;
		to_montgomery(xp, x1);
		to_montgomery(yp, y1);
		mont_msqr(tt, xp);
		// mont_a also works
#if	__cplusplus >= 201703L
		if constexpr(A_is_n3)
#else
		if (this->_a_is_neg3)
#endif
		{
			if ( tt.sub_from(this->_mont_three) ) tt.add_to(this->_p);
		} else {
			if (! this->_a_is_zero) return false;
			//if ( tt.add_to(this->_mont_a) ) tt.sub_from(this->_p);
		}
		// tt = x^3 + ax
		mont_mmult(tt, tt, xp);
		// tt = x^3 + ax, to normal bignum
		tt.mod_add_to(this->_mont_b, this->_p);
		mont_msqr(yy, yp);
		return yy == tt;
	}
	// affined point in montgomery form
	void to_affined(point_t<N>& pt, const u64 *x1, const u64 *y1) const noexcept
	{
		to_montgomery(pt.x, x1);
		to_montgomery(pt.y, y1);
		pt.z = mont_one();
	}
	void apply_z(felem_t& x, felem_t& y, const point_t<N>& pt) const noexcept
	{
		felem_t	z, t;
		from_montgomery(z, pt.z);
		mod_inv<N>(z, z, this->_p);	// z = p.z^-1
		to_montgomery(z, z);
		this->mont_msqr(t, z);	// t = z^2
		this->mont_mmult(x, pt.x, t);	// x1 * z^2
		this->mont_mmult(t, t, z);	// t = z^3
		this->mont_mmult(y, pt.y, t);	// y1 * z^3
		// montgomery reduction
		from_montgomery(x, x);
		from_montgomery(y, y);
	}
	void apply_z_mont(point_t<N>& pt) const noexcept
	{
		if ( unlikely(pt.z.is_zero() || pt.z == mont_one()) ) return;
		felem_t	t;
		from_montgomery(pt.z, pt.z);
		mod_inv<N>(pt.z, pt.z, this->_p);
		to_montgomery(pt.z, pt.z);		// z = p.z^-1
		this->mont_msqr(t, pt.z);	// t1 = z^2
		this->mont_mmult(pt.x, pt.x, t);	// x1 * z^2
		this->mont_mmult(t, t, pt.z);	// t1 = z^3
		this->mont_mmult(pt.y, pt.y, t);	// y1 * z^3
		pt.z = this->mont_one();
	}
	bool point_eq(const point_t<N>& p, const point_t<N>& q) const noexcept
	{
		if (p == q) return true;
		if ( likely(p.z == mont_one()) ) {
			felem_t	z(q.z);
			felem_t	x2, y2;
			from_montgomery(z, z);
			mod_inv<N>(z, z, this->_p);
			to_montgomery(z, z);
			this->mont_msqr(z, z);	// z = z^2
			this->mont_mmult(x2, p.x, z);	// x1 * z^2
			this->mont_mmult(z, z, z);	// z = z^3
			this->mont_mmult(y2, p.y, z);	// y1 * z^3
			return q.x == x2 && q.y == y2;
		} else {
			felem_t	x1, x2, y1, y2;
			this->apply_z(x1, y1, p);
			this->apply_z(x2, y2, q);
			return x1 == x2 && y1 == y2;
		}
	}
	bool point_eq(const spoint_t<N>& p, const point_t<N>& q) const noexcept
	{
		if (q.z != this->mont_one()) return false;
		return p.x == q.x && p.y == q.y;
	}
	void point_neg(point_t<N>& q, const point_t<N>& p) const noexcept
	{
		if ( unlikely(&q != &p) ) {
			q.x = p.x;
			q.z = p.z;
		}
		q.y.sub(this->_p, p.y);
	}
	void scalar_mult(point_t<N>& q, const spoint_t<N>& p, const felem_t& scalar,
			void *scratchBuff=nullptr) const noexcept
	{
		point_t<N>	tmp;
		uint	nbits = N*64; //scalar.num_bits();
		q.clear();
		if ( unlikely(scalar.is_zero()) ) return;
		point_t<N>	*pres;
		if (scratchBuff == nullptr) {
			pres = (point_t<N> *)this->presBuff;
		} else {
			pres = (point_t<N> *)scratchBuff;
		}
		to_montgomery(tmp.x, p.x);
		to_montgomery(tmp.y, p.y);
		tmp.z = this->mont_one();
		pre_compute<N>(*this, pres, tmp);
		int	off = nbits % W;
		if (off == 0) off = W;
		--nbits;
		int	idx = nbits - off;
		bool skip = true;
		{
			uint	bits;
			uint	digit;
			bits = vli_get_bits<N, W+1>(scalar.data(), idx);
			recode_scalar_bits<W>(digit, bits);
			if (digit != 0) {
				--digit;
				// sign MUST BE zero
				q = pres[digit];
				skip = false;
			}
		}
		for (; idx >= 0; ) {
			idx -= W;
			if (!skip) {
				for (int j=0; j<W ; ++j) point_double(q, q);
			}
			uint	bits;
			uint	digit;
			bits = vli_get_bits<N, W+1>(scalar.data(), idx);
			auto sign = recode_scalar_bits<W>(digit, bits);
			if (digit == 0) continue;
			--digit;
#ifndef	NO_CONDITIONAL_COPY
			tmp = pres[digit];
			felem_t	ny;
			ny.sub(this->_p, tmp.y);
			ny.copy_conditional(tmp.y, sign-1);
			tmp.y = ny;
#else
			if ( sign ) point_neg(tmp, pres[digit]); else tmp = pres[digit];
#endif
			if (!skip) point_add(q, q, tmp); else {
				q = tmp;
				skip = false;
			}
		}
	}
	// scalar MUST less than N, result never be infinite
	void scalar_mult(point_t<N>& q, const point_t<N>& p, const felem_t& scalar,
			void *scratchBuff=nullptr) const noexcept
	{
		if ( unlikely(p.is_zero()) ) {
			q = p;
			return;
		}
		spoint_t<N>	pp(p.x, p.y);
		scalar_mult(q, pp, scalar, scratchBuff);
		// montgomery reduction
#ifdef	ommit
		if ( unlikely(q.z.is_zero()) ) {
			q.x.clear();
			q.y.clear();
			return;
		}
#else
		if ( unlikely(q.z.is_zero()) ) return;
#endif
		this->apply_z_mont(q);
		this->from_montgomery(q.x, q.x);
		this->from_montgomery(q.y, q.y);
		q.z = felem_t(1);
	}
	void point_double(point_t<N>& q, const point_t<N>& p) const noexcept
	{
#if	defined(WITH_DBL_2004hmv) && __cplusplus >= 201703L
		if constexpr(A_is_n3)
		{
			point_double3n_jacob(*this, q.x, q.y, q.z, p.x, p.y, p.z);
			return;
		} else
#endif
		point_double_jacob<A_is_n3>(*this, q.x, q.y, q.z, p.x, p.y, p.z);
	}
	void point_double(point_t<N>& q, const felem_t& x1, const felem_t& y1)
	const noexcept
	{
		point_doublez_jacob<A_is_n3>(*this, q.x, q.y, q.z, x1, y1);
	}
	void point_add(point_t<N>& q, const point_t<N>& p1, const point_t<N>& p2)
	const noexcept
	{
		point_add_jacob<A_is_n3>(*this, q.x, q.y, q.z, p1.x, p1.y, p1.z,
						p2.x, p2.y, p2.z);
	}
	void point_add(point_t<N>& q, const point_t<N>& p1, const felem_t& x2,
					const felem_t& y2) const noexcept
	{
		point_addz_jacob<A_is_n3>(*this, q.x, q.y, q.z, p1.x, p1.y, p1.z,
						x2, y2);
	}
	const gwNAF_t&	select_base_NAF(const uint iLvl = 0) const noexcept {
		if ( unlikely(iLvl >= nBaseNAF) ) {
			static gwNAF_t	dummy;
			return dummy;
		}
		return g_precomps[iLvl];
	}
	bool select_base_point(spoint_t<N>& pt, const uint idx, const uint iLvl = 0)
		const noexcept
	{
		if ( unlikely(idx >= wBaseSize) ) return false;
		if ( unlikely(iLvl >= nBaseNAF) ) return false;
		pt = g_precomps[iLvl][idx];
		return true;
	}
	void scalar_mult_base(point_t<N>& q, const felem_t& scalar) const noexcept
	{
		if ( unlikely(nBaseNAF == 0) ) return;
		if ( unlikely(scalar.is_zero()) ) return;
		scalar_mult_base_internal(q, scalar);
		// montgomery reduction
		if ( unlikely(q.z.is_zero()) ) {
			q.x.clear();
			q.y.clear();
			return;
		}
		this->apply_z_mont(q);
		this->from_montgomery(q.x, q.x);
		this->from_montgomery(q.y, q.y);
		q.z = felem_t(1);
	}
	// scalar may be zero, g_scalar may be zero
	void cmult(point_t<N>& q, const point_t<N>& p, const felem_t& scalar,
			const felem_t& g_scalar, void *scratchBuff=nullptr) const noexcept
	{
		if ( unlikely(nBaseNAF == 0) ) return;
		if ( unlikely(g_scalar.is_zero()) ) return;
		scalar_mult_base_internal(q, g_scalar);
		if ( likely(!scalar.is_zero()) ) {
			point_t<N>	tmp;
			spoint_t<N> pp(p.x, p.y);
			scalar_mult(tmp, pp, scalar, scratchBuff);
			point_add(q, q, tmp);
		}
		// no montgomery reduction
	}
	void combined_mult(point_t<N>& q, const point_t<N>& p,
			const felem_t& scalar, const felem_t& g_scalar,
			void *scratchBuff=nullptr) const noexcept
	{
		if ( unlikely(nBaseNAF == 0) ) return;
		if ( unlikely(g_scalar.is_zero()) ) return;
		this->cmult(q, p, scalar, g_scalar, scratchBuff);
		// montgomery reduction
		if ( unlikely(q.z.is_zero()) ) {
			q.x.clear();
			q.y.clear();
			return;
		}
		this->apply_z_mont(q);
		this->from_montgomery(q.x, q.x);
		this->from_montgomery(q.y, q.y);
		q.z = felem_t(1);
	}
protected:
	bool initTable() noexcept
	{
		point_t<N>	G;
		static_assert(N == 4, "only 256 bits curve supported");
#if	__cplusplus >= 201703L
		if constexpr(BaseW == 6) {
			static_assert(nwBaseNAF<N>() == 43, "only 256 bits curve supported");
			static_assert(maxBaseNAF == 43, "only 256 bits curve supported");
		}
#endif
		static_assert((N * 64 % BaseW) != 0, "curve Bits MUUST not multiples of BaseW");
		if (nBaseNAF != 0) return true;
		nBaseNAF = maxBaseNAF;
		to_montgomery(G.x, this->getGx());
		to_montgomery(G.y, this->getGy());
		G.z = this->mont_one();

		point_t<N>	t1, t2;
		t2 = G;

		for (int j = 0; j < wBaseSize; ++j) {
			t1 = t2;
			for (int i = 0; i < maxBaseNAF; ++i) {
				// The window size is 6 so we need to double 6 times.
				// baseW is the window size
				if (i != 0) {
					this->point_double(t1, t1.x, t1.y);
					for (int k = 1; k < BaseW; ++k) {
						this->point_double(t1, t1);
					}
				}
				// Convert the point to affine form. (Its values are
				// still in Montgomery form however.)
				apply_z_mont(t1);

				// Update the table entry
				g_precomps[i][j].x = t1.x;
				g_precomps[i][j].y = t1.y;
			}
			if (j == 0) {
				this->point_double(t2, G.x, G.y);
			} else {
				this->point_add(t2, t2, G.x, G.y);
			}
		}
		return true;
	}
	void scalar_mult_base_internal(point_t<N>& q, const felem_t& scalar)
	const noexcept
	{
		if ( unlikely(nBaseNAF == 0) ) return;
		q.clear();
		if ( unlikely(scalar.is_zero()) ) return;
		bool skip = true;
		int	idx = -1;
		for (int iLvl = 0; iLvl < maxBaseNAF ; ++iLvl)
		{
			spoint_t<N>	tmp;
			uint	digit;
			uint	bits = vli_get_bits<N, BaseW+1>(scalar.data(), idx);
			auto sign = recode_scalar_bits<BaseW>(digit, bits);
			idx += BaseW;
			if (digit == 0) continue;
			--digit;
#ifndef	NO_CONDITIONAL_COPY
			if ( unlikely(select_base_point(tmp, digit, iLvl)) ) {
				// assert, should panic
			}
			felem_t	ny;
			ny.sub(this->_p, tmp.y);
			ny.copy_conditional(tmp.y, sign-1);
			tmp.y = ny;
#else
			if ( sign )
				point_neg(tmp, g_precomps[iLvl][digit]);
			else
				tmp = g_precomps[iLvl][digit];
#endif
			if (!skip) point_add(q, q, tmp.x, tmp.y); else {
				q.x = tmp.x;
				q.y = tmp.y;
				q.z = this->mont_one();
				skip = false;
			}
		}
	}
	bool init() noexcept
	{
		if (_inited) return _inited;
		if (unlikely(k0_p == 0))
		{
			// no calc k0 and rr after instantiation
			return false;
		}
		felem_t	t1;
		t1.clear();
		_mont_one.sub(t1, _p);
#ifndef	EXHAUSTIVE_TEST_ORDER
		_p_minus_n.sub(_p, _n);
#endif
		mont_mult2(_mont_three, _mont_one);
		mod_add_to(_mont_three, _mont_one);
		_quadP = _p;
		_quadP.rshift1();
		_quadP.rshift1();
#if	defined(WITH_ECDSA) && !defined(NO_HALF_N)
		_half_n = _n;
		_half_n.rshift1();
#endif
		if (unlikely( !_a_is_neg3 )) {
			if (likely(_a.is_zero())) _a_is_zero = true;
		}
		to_montgomery(_mont_a, this->_a);
		to_montgomery(_mont_b, this->_b);
		// should verify calc K0 and RR
		_inited = true;
		return this->initTable();
		//return true;
	}
	const std::string name;
	const felem_t gx;
	const felem_t gy;
	const felem_t _p;
	const felem_t _n;
	const felem_t _a;
	const felem_t _b;
	const felem_t rr_p = {};
	const felem_t rr_n = {};
	gwNAF_t	g_precomps[maxBaseNAF];
	felem_t _mont_one;
	felem_t _mont_three;
	felem_t _mont_a;
	felem_t _mont_b;
	felem_t	_quadP;
#ifndef	EXHAUSTIVE_TEST_ORDER
	felem_t	_p_minus_n;
#endif
#if	defined(WITH_ECDSA) && !defined(NO_HALF_N)
	felem_t _half_n;
#endif
	const u64	k0_p = 0;
	const u64	k0_n = 0;
	const uint _ndigits = N;
	const bool _a_is_neg3 = false;
	uint	nBaseNAF=0;
	bool _a_is_zero = false;
	bool _inited = false;
	point_t<N>  presBuff[wSize];
};


// curve256 for SM2
class alignas(64) curve256 : public ecc_curve<4,true> {
public:
	using felem_t = bignum<4>;
	curve256(const char *_name, const u64 *_gx, const u64 *_gy, const u64 *p,
			const u64 *n, const u64 *a, const u64 *b) :
		ecc_curve<4>(_name, _gx, _gy, p, n, a, b, sm2_p_rr, sm2_n_rr,
					sm2_p_k0, sm2_n_k0)
	{
		static_assert(sm2_p_k0 == 1, "MUST be sm2");
	}
	const felem_t& mont_one() const noexcept { return this->_mont_one; }
	void to_montgomery(felem_t& res, const u64 *x) const noexcept
	{
		u64   *resp = reinterpret_cast<u64 *>(&res);
		sm2p_mult(resp, x, this->rr_p.data());
	}
	void to_montgomery(felem_t& res, const felem_t& x) const noexcept
	{
		u64   *resp = reinterpret_cast<u64 *>(&res);
		sm2p_mult(resp, x.data(), this->rr_p.data());
	}
	void from_montgomery(felem_t& res, const felem_t& y) const noexcept
	{
		u64   *resp = reinterpret_cast<u64 *>(&res);
		sm2p_reduction(resp, y.data());
	}
	void from_montgomery(u64* result, const felem_t& y) const noexcept
	{
		sm2p_reduction(result, y.data());
	}
	void apply_z_mont(point_t<4>& pt) const noexcept
	{
		if (pt.z.is_zero() || pt.z == mont_one()) return;
		felem_t	t, z;
		from_montgomery(z, pt.z);
		mod_inv<4>(z, z, this->_p);
		to_montgomery(z, z);		// z = p.z^-1
		this->mont_msqr(t, z);	// t = z^2
		this->mont_mmult(pt.x, pt.x, t);	// x1 * z^2
		this->mont_mmult(t, t, z);	// t = z^3
		this->mont_mmult(pt.y, pt.y, t);	// y1 * z^3
		pt.z = this->mont_one();
	}
	void mont_mmult(felem_t& res, const felem_t& left, const felem_t& right)
	const noexcept
	{
		u64   *resp = reinterpret_cast<u64 *>(&res);
		sm2p_mult(resp, left.data(), right.data());
	}
	void mont_msqr(felem_t& res, const felem_t left, const uint nTimes=1)
	const noexcept
	{
		u64   *resp = reinterpret_cast<u64 *>(&res);
		sm2p_sqrN(resp, left.data());
		for (uint i=1; i < nTimes; i++) sm2p_sqrN(resp, resp);
	}
	// left,right less than p, result may large than p
	void mod_add(felem_t& res, const felem_t& left, const felem_t& right)
	const noexcept
	{
		auto carry = res.add(left, right);
		u64   *resp = reinterpret_cast<u64 *>(&res);
		sm2p_mod(resp, res.data(), carry);
	}
	// left,right less than p, result may large than p
	void mod_add_to(felem_t& res, const felem_t& right) const noexcept
	{
		auto carry = res.add_to(right);
		u64   *resp = reinterpret_cast<u64 *>(&res);
		sm2p_mod(resp, res.data(), carry);
	}
	void mont_mult2(felem_t& res, const felem_t& left) const noexcept
	{
#ifdef	ommit
		this->mod_add(res, left, left);
#else
		if (res.lshift1(left) != 0) {
			//this->carry_reduce(res, 1);
			res.sub_from(this->_p);
		}
#endif
	}
	void mont_mult4(felem_t& res) const noexcept
	{
		u64	carry;
		if ((carry=res.lshift(2)) != 0) {
			this->carry_reduce(res, carry);
		}
	}
	void mont_mult8(felem_t& res) const noexcept
	{
		u64	carry;
		if ((carry=res.lshift(3)) != 0) {
			this->carry_reduce(res, carry);
		}
	}
	void mont_div2(felem_t& res) const noexcept
	{
		if ( likely(res.is_even()) ) {
			res.rshift1();
		} else {
			bool carry = res.add_to(this->_p);
			res.rshift1(carry);
		}
	}
	bool prod_equal(const felem_t& prod, const felem_t& x, const felem_t& yp)
	const noexcept
	{
		felem_t	xp;
		this->to_montgomery(xp, x);
		this->mont_mmult(xp, xp, yp);
		return prod == xp;
	}
	void
	mod_exp(felem_t& res, const felem_t &x1, const felem_t& y1) const noexcept
	{
		felem_t		tmp;
		to_montgomery(tmp, x1);
		ecc::mont_mod_exp<4>(*this, res, tmp, y1);
		from_montgomery(res, res);
	}
	bool mont_sqrt(felem_t& res, const felem_t& xp) const noexcept
	{
		felem_t	a1;
		//ecc::mont_mod_exp<4>(*this, a1, xp, this->quadP());
		this->mont_mod_exp_quadP(a1, xp);
		felem_t	a0;
		mont_mmult(res, a1, xp);
		mont_mmult(a0, res, a1);
		return (a0 == this->mont_one());
	}
	bool mod_sqrt(felem_t& res, const felem_t& x) const noexcept
	{
		felem_t	xp;
		to_montgomery(xp, x);
		auto ret = mont_sqrt(res, xp);
		from_montgomery(res, res);
		return ret;
	}
	// affined point in montgomery form
	void to_affined(point_t<4>& pt, const u64 *x1, const u64 *y1) const noexcept
	{
		to_montgomery(pt.x, x1);
		to_montgomery(pt.y, y1);
		pt.z = mont_one();
	}
	void apply_z(felem_t& x, felem_t& y, const point_t<4>& pt) const noexcept
	{
		felem_t	z, t;
		from_montgomery(z, pt.z);
		mod_inv<4>(z, z, this->_p);
		to_montgomery(z, z);	// z = pt.z^-1
		this->mont_msqr(t, z);	// t = z^2
		this->mont_mmult(x, pt.x, t);	// x1 * z^2
		this->mont_mmult(t, t, z);	// t = z^3
		this->mont_mmult(y, pt.y, t);	// y1 * z^3
		// montgomery reduction
		from_montgomery(x, x);
		from_montgomery(y, y);
	}
	void scalar_mult(point_t<4>& q, const spoint_t<4>& p, const felem_t& scalar,
			void *scratchBuff=nullptr) const noexcept
	{
		point_t<4>	tmp;
		uint	nbits = 4*64;
		q.clear();
		if ( unlikely(scalar.is_zero()) ) return;
		point_t<4>	*pres;
		if (scratchBuff == nullptr) {
			pres = (point_t<4> *)this->presBuff;
		} else {
			pres = (point_t<4> *)scratchBuff;
		}
		to_montgomery(tmp.x, p.x);
		to_montgomery(tmp.y, p.y);
		tmp.z = this->mont_one();
		pre_compute<4>(*this, pres, tmp);
		int	off = nbits % W;
		if (off == 0) off = W;
		--nbits;
		int	idx = nbits - off;
		bool skip = true;
		{
			uint	bits;
			uint	digit;
			bits = vli_get_bits<4, W+1>(scalar.data(), idx);
			recode_scalar_bits<W>(digit, bits);
			if (digit != 0) {
				--digit;
				// sign MUST BE zero
				q = pres[digit];
				skip = false;
			}
		}
		for (; idx >= 0; ) {
			idx -= W;
			if (!skip) {
				for (int j=0; j<W ; ++j) point_double(q, q);
			}
			uint	bits;
			uint	digit;
			bits = vli_get_bits<4, W+1>(scalar.data(), idx);
			auto sign = recode_scalar_bits<W>(digit, bits);
			if (digit == 0) continue;
			--digit;
#ifndef	NO_CONDITIONAL_COPY
			tmp = pres[digit];
			felem_t	ny;
			ny.sub(this->_p, tmp.y);
			ny.copy_conditional(tmp.y, sign-1);
			tmp.y = ny;
#else
			if ( sign ) point_neg(tmp, pres[digit]); else tmp = pres[digit];
#endif
			if (!skip) point_add(q, q, tmp); else {
				q = tmp;
				skip = false;
			}
		}
	}
	void scalar_mult(point_t<4>& q, const point_t<4>& p, const felem_t& scalar,
			void *scratchBuff=nullptr) const noexcept
	{
		if ( unlikely(p.is_zero()) ) {
			q = p;
			return;
		}
		spoint_t<4>	pp(p.x, p.y);
		scalar_mult(q, pp, scalar, scratchBuff);
		// montgomery reduction
		if ( unlikely(q.z.is_zero()) ) return;
		this->apply_z_mont(q);
		this->from_montgomery(q.x, q.x);
		this->from_montgomery(q.y, q.y);
		q.z = felem_t(1);
	}
	void point_double(point_t<4>& q, const point_t<4>& p) const noexcept
	{
		if ( p.z == this->mont_one() )
			point_doublez_jacob<true>(*this, q.x, q.y, q.z, p.x, p.y);
		else
#if	defined(WITH_DBL_2004hmv) && __cplusplus >= 201703L
			point_double3n_jacob(*this, q.x, q.y, q.z, p.x, p.y, p.z);
#else
			point_double_jacob<true>(*this, q.x, q.y, q.z, p.x, p.y, p.z);
#endif
	}
	void point_double(point_t<4>& q, const felem_t& x1, const felem_t& y1)
	const noexcept
	{
		point_doublez_jacob<true>(*this, q.x, q.y, q.z, x1, y1);
	}
	void point_add(point_t<4>& q, const point_t<4>& p1, const point_t<4>& p2)
	const noexcept
	{
		point_add_jacob<true>(*this, q.x, q.y, q.z, p1.x, p1.y, p1.z,
						p2.x, p2.y, p2.z);
	}
	void point_add(point_t<4>& q, const point_t<4>& p1, const felem_t& x2,
			const felem_t& y2) const noexcept
	{
		point_addz_jacob<true>(*this, q.x, q.y, q.z, p1.x, p1.y, p1.z, x2, y2);
	}
	void scalar_mult_base(point_t<4>& q, const felem_t& scalar) const noexcept
	{
		if ( unlikely(this->nBaseNAF == 0) ) return;
		if ( unlikely(scalar.is_zero()) ) return;
		scalar_mult_base_internal(q, scalar);
		// montgomery reduction
		if ( unlikely(q.z.is_zero()) ) return;
		this->apply_z_mont(q);
		this->from_montgomery(q.x, q.x);
		this->from_montgomery(q.y, q.y);
		q.z = felem_t(1);
	}
	// scalar may be zero, g_scalar may be zero
	void cmult(point_t<4>& q, const point_t<4>& p, const felem_t& scalar,
			const felem_t& g_scalar, void *scratchBuff=nullptr) const noexcept
	{
		if ( unlikely(this->nBaseNAF == 0) ) return;
		if ( unlikely(g_scalar.is_zero()) ) return;
		scalar_mult_base_internal(q, g_scalar);
		if ( likely(!scalar.is_zero()) ) {
			point_t<4>	tmp;
			spoint_t<4> pp(p.x, p.y);
			scalar_mult(tmp, pp, scalar, scratchBuff);
			point_add(q, q, tmp);
		}
		// no montgomery reduction
	}
	void combined_mult(point_t<4>& q, const point_t<4>& p,
			const felem_t& scalar, const felem_t& g_scalar,
			void *scratchBuff=nullptr) const noexcept
	{
		if ( unlikely(this->nBaseNAF == 0) ) return;
		if ( unlikely(g_scalar.is_zero()) ) 
		{
			scalar_mult(q, p, scalar);
			return;
		}
		this->cmult(q, p, scalar, g_scalar, scratchBuff);
		// montgomery reduction
		if ( unlikely(q.z.is_zero()) ) {
			q.x.clear();
			q.y.clear();
			return;
		}
		this->apply_z_mont(q);
		this->from_montgomery(q.x, q.x);
		this->from_montgomery(q.y, q.y);
		q.z = felem_t(1);
	}
private:
	void carry_reduce(felem_t& res, const u64 carry) const noexcept
	{
		//static_assert(Pk0 == 1, "MUST be sm2");
		// carry < 2^32
		u64		u = carry & ((1L<<32) -1);
		u64		cc[4];
		cc[0] = u;
		cc[1] = (u << 32) - u;
		cc[2] = 0;
		cc[3] = u << 32;
		if (res.add_to(cc)) res.sub_from(this->_p);
	}
	// x^quadP modulo prime, xp is x in motgomery form
	void mont_mod_exp_quadP(felem_t& res, const felem_t& xp)
	const noexcept
	{
		res = xp;

		felem_t	bn_tbl[5];	// _1, _11, _101, _111, _1111
		felem_t	tmp, x8, x16, x32;
		// precompute
		bn_tbl[0] = xp;
		mont_msqr(tmp, xp);
		mont_mmult(bn_tbl[1], tmp, xp);
		mont_mmult(bn_tbl[2], tmp, bn_tbl[1]);	// _101
		mont_mmult(bn_tbl[3], tmp, bn_tbl[2]);	// _111
		mont_msqr(tmp, bn_tbl[2]);	// _1010
		mont_mmult(bn_tbl[4], bn_tbl[2], tmp);	// _1111
		mont_msqr(tmp, bn_tbl[4], 4);
		mont_mmult(x8, tmp, bn_tbl[4]);
		mont_msqr(tmp, x8, 8);
		mont_mmult(x16, tmp, x8);
		mont_msqr(tmp, x16, 16);
		mont_mmult(x32, tmp, x16);

		mont_msqr(res, x16, 8);
		mont_mmult(res, res, x8);
		mont_msqr(res, res, 4);
		mont_mmult(res, res, bn_tbl[4]);	// _1111
		mont_msqr(res, res, 3);
		mont_mmult(res, res, bn_tbl[3]);	// _111
		mont_msqr(res, res);
		mont_msqr(res, res, 32);
		mont_mmult(res, res, x32);
		mont_msqr(res, res, 32);
		mont_mmult(res, res, x32);
		mont_msqr(res, res, 32);
		mont_mmult(res, res, x32);
		mont_msqr(res, res, 32);
		mont_mmult(res, res, x32);
		mont_msqr(res, res, 32);
		mont_msqr(res, res, 32);
		mont_mmult(res, res, x32);
		mont_msqr(res, res, 16);
		mont_mmult(res, res, x16);
		mont_msqr(res, res, 8);
		mont_mmult(res, res, x8);
		mont_msqr(res, res, 4);
		mont_mmult(res, res, bn_tbl[4]);	// _1111
		mont_msqr(res, res, 2);
		mont_mmult(res, res, bn_tbl[1]);	// _11
	}
	void scalar_mult_base_internal(point_t<4>& q, const felem_t& scalar)
	const noexcept
	{
		if ( unlikely(this->nBaseNAF == 0) ) return;
		q.clear();
		if ( unlikely(scalar.is_zero()) ) return;
		bool skip = true;
		int	idx = -1;
		for (int iLvl = 0; iLvl < nwBaseNAF<4>() ; ++iLvl)
		{
			spoint_t<4>	tmp;
			uint	digit;
			uint	bits = vli_get_bits<4, BaseW+1>(scalar.data(), idx);
			auto sign = recode_scalar_bits<BaseW>(digit, bits);
			idx += BaseW;
			if (digit == 0) continue;
			--digit;
#ifndef	NO_CONDITIONAL_COPY
			if ( unlikely(this->select_base_point(tmp, digit, iLvl)) ) {
				// assert, should panic
			}
			felem_t	ny;
			ny.sub(this->_p, tmp.y);
			ny.copy_conditional(tmp.y, sign-1);
			tmp.y = ny;
#else
			if ( sign )
				point_neg(tmp, g_precomps[iLvl][digit]);
			else
				tmp = g_precomps[iLvl][digit];
#endif
			if (!skip) point_add(q, q, tmp.x, tmp.y); else {
				q.x = tmp.x;
				q.y = tmp.y;
				q.z = this->mont_one();
				skip = false;
			}
		}
	}
};


struct slice_t {
	u64	*data;
	int64_t	len;
	int64_t	cap;
	bool isZero() {
		if (len <= 0) return true;
		for (int i=0;i<len;++i) {
			if (data[i] != 0) return false;
		}
		return true;
	}
	explicit operator bool () { return len != 0; }
	template <const uint N>slice_t(u64 vd[N]) noexcept : data(vd),
			len(N), cap(N)
	{
		vli_clear<N>(data);
	}
	void normal() {
		if (len <= 0) return;
		int	i;
		for(i=len-1;i>=0;i--) {
			if (data[i] != 0) break;
		}
		len = i+1;
	}
};


}	// namespace ecc

#endif	//	__CURVE_IMPL_HPP__
