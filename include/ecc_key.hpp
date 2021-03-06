/*
 * Copyright (c) 2013, Kenneth MacKay
 * All rights reserved.
 * Copyright(c) 2020 Jesse Kuang  <jkuang@21cn.com>
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
#ifndef	__ECC_KEY_H__
#define __ECC_KEY_H__

#include "cdefs.h"
#include <time.h>
#include <unistd.h>
#include <random>
#ifdef	__x86_64__
#include <x86intrin.h>
#endif
#include "ecc_impl.hpp"


namespace ecc {

using namespace vli;


template<const uint N=4>
class	bn_random {
public:
	bn_random(const bn_random&) = delete;
	static bn_random&	Instance() noexcept {
		static	bn_random	bn_random_;
		static	bool	_inited = false;
		if (!_inited) {
			u64		seeds[4];
#ifdef	WITH_GETENTROPY
			if (getentropy((void *)seeds, sizeof(seeds)) < 0)
#endif
			{
				seeds[0] = clock();
				struct timespec tp;
				clock_gettime(CLOCK_REALTIME, &tp);
				seeds[1] = tp.tv_sec << 24 | (tp.tv_nsec & 0xffffff);
				clock_gettime(CLOCK_MONOTONIC, &tp);
				seeds[2] = tp.tv_sec << 24 | (tp.tv_nsec & 0xffffff);
#ifdef	__x86_64__
				unsigned int aux;
				__rdtscp(&aux);
				seeds[3] = aux;
#else
				seeds[3] = clock();
				clock_gettime(CLOCK_REALTIME, &tp);
				seeds[3] += tp.tv_nsec;
#endif
			}
			for (int i=0; i<4; ++i) bn_random_._rd[i].seed(seeds[i]);
			_inited = true;
		}
		return bn_random_;
	}
	bignum<N> get_random() noexcept {
		u64		dd[N];
		for (uint i=0;i<N;++i) dd[i] = _rd[i&3]();
		return bignum<N>(dd);
	}
private:
	bn_random() {
		_rd[0].seed(clock());
	}
	std::mt19937_64	_rd[4];
};


/*
 * ECC private keys are generated using the method of extra random bits,
 * equivalent to that described in FIPS 186-4, Appendix B.4.1.
 *
 * d = (c mod(n???1))	where c is a string of random bits, 64 bits longer
 *						than requested, if d == 0 then d++
 * 0 <= c mod(n-1) <= n-2  and implies that
 * 1 <= d <= n-2
 *
 * This method generates a private key uniformly distributed in the range
 * [1, n-2].
 */
// definitions for Private/Public key
template<const uint N=4>
class	private_key {
public:
	using felem_t = bignum<N>;
	using public_key = spoint_t<N>;
	private_key() = default;
	template<typename curveT>
	private_key(const curveT& curve)
	{
		auto&	rd = bn_random<N>::Instance();
		felem_t secret;
		do {
			secret = rd.get_random();
		} while ( !init(curve, secret) );
		calc_pa(curve);
	}
	template<typename curveT>
	private_key(const curveT& curve, const felem_t& secret,
				const public_key& pk) : _pa(pk)
	{
		init(curve, secret);
	}
	template<typename curveT>
	private_key(const curveT& curve, const felem_t& secret)
	{
		init(curve, secret);
	}
	explicit operator bool() const noexcept { return _inited; }
	const felem_t&		D() const noexcept { return *(&_d); }
#ifdef	WITH_MONT_D
	const felem_t&		mont_d() const noexcept { return *(&_mont_d); }
#endif
	// Di() return (1 + dA)^-1 in montgomery form modN
	const felem_t&		Di() const noexcept { return *(&_dInv); }
	const public_key&	PubKey() const noexcept { return *(&_pa); }
protected:
	template<typename curveT>
	bool init(const curveT& curve, const felem_t& secret) noexcept {
		_inited = false;
		bignum<N>	n1;
		n1.usub(curve.paramN(), 1);
		// mod N-1
		if (secret >= n1) _d.sub(secret, n1); else _d = secret;
		if (_d.is_zero()) _d.uadd_to(1);
		// _d in [1 .. N-2]
		_dInv = _d;
		_dInv.uadd_to(1);
		if (_dInv == curve.paramN()) return false;
		// _dInv = (1 + d)^-1
		mod_inv<N>(_dInv, _dInv, curve.paramN());
		curve.to_montgomeryN(_dInv, _dInv);
#ifdef	WITH_MONT_D
		curve.to_montgomeryN(_mont_d, _d);
#endif
		_inited = true;
		return true;
	}
	template<typename curveT>
	void calc_pa(const curveT& curve) noexcept {
		if ( unlikely(!_inited) ) return;
		point_t<N>	pt;
		curve.scalar_mult_base(pt, _d);
		_pa.x= pt.x;
		_pa.y = pt.y;
	}
	bool		_inited = false;
	felem_t		_d;			// secret
#ifdef	WITH_MONT_D
	felem_t		_mont_d;	// secret in montgomery form
#endif
	felem_t		_dInv;		// (1+d)^-1  for SM2 sign
	public_key	_pa;
};


template<const uint N=4, typename curveT>
forceinline static
void gen_keypair(const curveT& curve, bignum<N>& secret, bignum<N>& pubX,
		bignum<N>& pubY) noexcept
{
	auto&	rd = bn_random<N>::Instance();
	bignum<N>	n1;
	n1.usub(curve.paramN(), 2);
	// mod N-2, add_to 1
	secret = rd.get_random();
	if ( unlikely(secret >= n1) ) secret.sub_from(n1);
	secret.uadd_to(1);
	point_t<N>	pt;
	curve.scalar_mult_base(pt, secret);
	pubX = pt.x;
	pubY = pt.y;
}

/*
 * SM2 sign/verify
 * 			e = HASH (msg)
 * 			Pa ...	Public key, Px/Py
 * 			dA ...	private key
 * 	sign
 * 			k = random   256bits   Label getk
 * 			x1, y1 = k * G
 * 			r = e + x1
 *			if r == 0 modN  goto getk
 *			s = (1 + dA)^-1 * (k - r * dA)
 *			s = (1 + dA)^-1 * (k + r) - r
 *			if s == 0 modN goto getk
 *			if r + s == 0 modN goto getk
 *			return r, s, y1_is_odd
 *
 * 	verify
 * 			t = r + s modN
 * 			if t == 0 fail
 * 			x2, y2 = s*G + t*Pa
 * 			if x2 + e == r modN success else fail
 *
 *	recover public key
 *			Psx = r - e
 *			Psy = get_Py ( Psx, Psy_is_odd )
 *			Ps ...	Psx/Psy
 *			u1 = (r + s)^-1
 *			u2 = s * u1
 *			Pa = (r+s)^1 * ( k*G - s*G) = (r+s)^-1 * Ps - s*(r+s)^-1 * G
 *			Pa = u1 * Ps - u2 * G
 */

// sign return ecInd,  0  for Py(of sign) even, 1 odd
template<const uint N=4, typename curveT>
forceinline static
int ec_sign(const curveT& curve, bignum<N>& r, bignum<N>& s,
		const private_key<N>& priv, const bignum<N>& msg) noexcept
{
	bignum<N>	k, x1;	// y1 reuse tmp
	bignum<N>	tmp;	// r+s
	int		ret=0;
	do {
		do {
			gen_keypair<N>(curve, k, x1, tmp);
			r.mod_add(msg, x1, curve.paramN());
		} while (r.is_zero());
		ret = tmp.is_odd();
#ifdef	WITH_MONT_D
		bignum<N>	rp;	// r+s
		curve.to_montgomeryN(k, k);
		curve.to_montgomeryN(rp, r);
		// tmp = r * dA
		curve.mont_nmult(tmp, rp, priv.mont_d());
		// k = k - r * dA
		if (k.sub_from(tmp)) k.add_to(curve.paramN());
		// tmp = (1 + dA)^-1 * (k - r * dA)
		curve.mont_nmult(tmp, priv.Di(), k);
		curve.from_montgomeryN(s, tmp);
#else
		tmp.mod_add(k, r, curve.paramN());
		curve.to_montgomeryN(tmp, tmp);
		curve.mont_nmult(tmp, tmp, priv.Di());
		curve.from_montgomeryN(tmp, tmp);
		s.mod_sub(tmp, r, curve.paramN());
#endif
		if (s.is_zero()) continue;
		// SM2 sign s may above halfN
		tmp.mod_add(r, s, curve.paramN());
	} while (tmp.is_zero());
	return ret;
}

// verify return bool, true for success(ok)
template<const uint N=4, typename curveT>
forceinline static
bool ec_verify(const curveT& curve, const bignum<N>& r, const bignum<N>& s,
		const spoint_t<N>&	pub, const bignum<N>& msg) noexcept
{
	bignum<N>	t;
	if ( unlikely(r.is_zero()) ) return false;
	if ( unlikely(s.is_zero()) ) return false;
	if ( unlikely(curve.paramN() < r) ) return false;
	if ( unlikely(curve.paramN() < s) ) return false;
	t.mod_add(r, s, curve.paramN());
	if ( unlikely(t.is_zero()) ) return false;
	point_t<N>	q, p(pub.x, pub.y);
	// q = s*G + t * Pub
	// x2, y2 = s*G + t * Pub
#ifdef	EXHAUSTIVE_TEST_ORDER
	curve.combined_mult(q, p, t, s);
	if ( unlikely(q.x.is_zero()) ) return false;
	// t = x2 + msg modN
	t.mod_add(q.x, msg, curve.paramN());
	if ( unlikely(t.is_zero()) ) return false;
	return r == t;
#else
	curve.cmult(q, p, t, s);
	if ( unlikely(q.z.is_zero()) ) return false;
	// t = r - msg
	t.mod_sub(r, msg, curve.paramN());
	bignum<N>	zz;
	curve.mont_msqr(zz, q.z);
	if (curve.prod_equal(q.x, t, zz)) return true;
	if (t < curve.P_minus_N()) {
		t.add_to(curve.paramN());
		if (curve.prod_equal(q.x, t, zz)) return true;
	}
	return false;
#endif
}

// vecover public key, from r/s/v/msg, v for pubY is odd
// return true for successfully recover public key
template<const uint N=4, typename curveT>
forceinline static
bool ec_recover(const curveT& curve, spoint_t<N>&  pub, const bignum<N>& r,
				const bignum<N>& s, const int v, const bignum<N>& msg) noexcept
{
	point_t<N>	p;
	// p.x = r - msg
	p.x.mod_sub(r, msg, curve.paramN());
	if ( unlikely(!pointY_recover(curve, p.y, p.x, v)) ) return false;
	bignum<N>	u1, u2;
	u1.mod_add(r, s, curve.paramN());
	if ( unlikely(u1.is_zero()) ) return false;
	// u1 = (r + s)^-1
	mod_inv<N>(u1, u1, curve.paramN());
	bignum<N>	u1p;
	curve.to_montgomeryN(u2, s);
	curve.to_montgomeryN(u1p, u1);
	curve.mont_nmult(u2, u2, u1p);
	// u2 = s * u1
	curve.from_montgomeryN(u2, u2);
	// 0 < u2 < N 
	if ( unlikely(u2.is_zero()) ) return false;
	// u2 = -u2 = - s* u1
	u2.sub(curve.paramN(), u2);
	// Pa = u1 * Ps - u2 * G
	point_t<N>	q;
	curve.combined_mult(q, p, u1, u2);
	pub.x = q.x;
	pub.y = q.y;
	return true;
}

}	// namespace ecc

#endif	//	__ECC_KEY_H__
