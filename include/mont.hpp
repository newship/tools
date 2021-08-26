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
#ifndef __MONT_HPP__
#define __MONT_HPP__

#include "vli.hpp"
#include "curve_const.hpp"		// include sm2_p, sm2_p_rr
#include "mont_asm.h"

// SM2 prime optimize
// p is 2^256 - 2^224 - 2^96 + 2^64 -1
// p is 2²⁵⁶ - 2²²⁴ - 2⁹⁶ + 2⁶⁴ -1
forceinline static void vli_sm2_multP(u64 *result, const u64 u) noexcept
{
	u64	t_low, t_high;
	t_low = u << 32;	// ^192, ^96
	t_high = u >> 32;
	u64		carry = 0;
	result[0] = u64_subc(0, u, carry);
	result[1] = u64_subc(u, t_low, carry);
	result[2] = u64_subc(0, t_high, carry);
	result[3] = u64_subc(0, t_low, carry);
	result[4] = u64_subc(u, t_high, carry);
	//result[5] = carry;
}


// secp256k1 prime optimize
// p is 2^256 - 2^32 - 0x3d1 = 2^256 - 0x1000003d1
// secp256k1 p: 2²⁵⁶ - 2³² - 2¹⁰ + 2⁵ + 2⁴ - 1
forceinline static void vli_btc_multP(u64 *result, const u64 u) noexcept
{
	u64		carry = 0;
	uint128_t	pd;
	pd.mul_64_64(0x1000003d1, u);
	result[0] = u64_subc(0, pd.m_low(), carry);
	result[1] = u64_subc(0, pd.m_high(), carry);
	result[2] = u64_subc(0, 0, carry);
	result[3] = u64_subc(0, 0, carry);
	result[4] = u64_subc(u, 0, carry);
}

/* result > mod (result = mod + remainder), so subtract mod to
 * get remainder.
 */
static forceinline void
sm2p_mod(u64 *res, const u64 *left, const bool carry) noexcept
{
#if	defined(__x86_64__) || defined(__aarch64__)
	sm2p_mod_asm(res, left, carry);
#else
	vli_mod<4>(res, left, sm2_p, carry);
#endif
}


forceinline static void
sm2p_reductionStep(u64& r0, u64& r1, u64& r2, u64& r3, u64& carry)
{
	u64	u = r0;
	r0 = carry;		// rshift1w
	u64 cc = 0;
	r1 = u64_addc(r1, u, cc);
	r2 = u64_addcz(r2, cc);
	r3 = u64_addcz(r3, cc);
	r0 = u64_addc(r0, u, cc);
	carry = cc;
	u64 t_low = u << 32;	// ^192
	u64 t_high = u >> 32;
	cc = 0;
	r1 = u64_subc(r1, t_low, cc);
	r2 = u64_subc(r2, t_high, cc);
	r3 = u64_subc(r3, t_low, cc);
	r0 = u64_subc(r0, t_high, cc);
	carry -= cc;
}

forceinline static void sm2p_reductionN(u64 *result,
		const u64 *y, const bool isProd=false) noexcept
{
	u64	r0, r1, r2, r3;
	vli4_load(y, r0, r1, r2, r3);
	u64 carry = 0;
	sm2p_reductionStep(r0, r1, r2, r3, carry);
	sm2p_reductionStep(r1, r2, r3, r0, carry);
	sm2p_reductionStep(r2, r3, r0, r1, carry);
	sm2p_reductionStep(r3, r0, r1, r2, carry);

	// add high 256 bits
	if ( unlikely(isProd) )
	{
		u64	cc=0;
		r0 = u64_addc(r0, y[4], cc);
		r1 = u64_addc(r1, y[5], cc);
		r2 = u64_addc(r2, y[6], cc);
		r3 = u64_addc(r3, y[7], cc);
		carry += cc;
	}
	// sm2p_mod
	{
		u64	cc=0;
		u64 s0 = u64_subc(r0, sm2_p[0], cc);
		u64 s1 = u64_subc(r1, sm2_p[1], cc);
		u64 s2 = u64_subc(r2, sm2_p[2], cc);
		u64 s3 = u64_subc(r3, sm2_p[3], cc);
		u64_subcz(carry, cc);
		if (cc != 0) vli4_save(result, r0, r1, r2, r3); else
			vli4_save(result, s0, s1, s2, s3);
	}
}


forceinline static void sm2p_reduction(u64 *result, const u64 *y) noexcept
{
#if	defined(__x86_64__) || defined(__aarch64__)
	u64	r0=y[0], r1=y[1], r2=y[2] ,r3=y[3];
	sm2p_reductionS0(r0, r1, r2, r3);
	sm2p_reductionS2(result, r0, r1, r2, r3);
#else
	sm2p_reductionN(result, y);
#endif
}

forceinline static void sm2p_reductionProd(u64 *result, const u64 *y) noexcept
{
	//u64	r0, r1, r2 ,r3;
	u64	r0=y[0], r1=y[1], r2=y[2] ,r3=y[3];
	sm2p_reductionS0(r0, r1, r2, r3);
	// S1, with C/C++ u64_addc
#if	defined(__aarch64__)
	register u64 res0 asm("x4") = r0;
	register u64 res1 asm("x5") = r1;
	register u64 res2 asm("x6") = r2;
	register u64 res3 asm("x7") = r3;
	register u64 r4 asm("x9") = y[4];
	register u64 r5 asm("x10") = y[5];
	register u64 r6 asm("x11") = y[6];
	register u64 r7 asm("x12") = y[7];
	register u64 carry asm("x0") = 0;
	asm volatile(
	// Add bits [511:256] of the product y
		"ADDS	x4, x9, x4\n"
		"ADCS	x5, x10, x5\n"
		"ADCS	x6, x11, x6\n"
		"ADCS	x7, x12, x7\n"
		"ADC	x0, xzr, xzr\n"	// x0 = carry
		: "+r" (carry), "+r" (res0), "+r" (res1), "+r" (res2), "+r" (res3)
		: "r" (r4), "r" (r5), "r" (r6), "r" (r7)
		: "cc");
	sm2p_reductionS2(result, res0, res1, res2, res3, carry);
#else
	u64	carry=0;
	{
		u64     cc=0;
		r0 = u64_addc(r0, y[4], cc);
		r1 = u64_addc(r1, y[5], cc);
		r2 = u64_addc(r2, y[6], cc);
		r3 = u64_addc(r3, y[7], cc);
		carry += cc;
	}
	sm2p_reductionS2(result, r0, r1, r2, r3, carry);
#endif
}

template<const uint N> forceinline
static void vli_mont_reduction(u64 *result, const u64 *y,
		const u64 *prime, const u64 k0) noexcept
{
	u64	s[N + 1];
	u64	r[N];
	vli_set<N>(r, y);
	for (uint i=0; i < N; i++) {
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		s[N] += vli_add_to<N>(r, s);
		vli_rshift1w<N>(r, s[N]);
	}
	vli_mod<N>(result, r, prime);
}

template<const uint N> forceinline
static void vli_mont_mult(u64 *result, const u64 *x,
		const u64 *y, const u64 *prime, const u64 k0) noexcept
{
	u64	s[N + 1];
	u64	r[N + 1];
	vli_clear<N + 1>(r);
	for (uint i=0; i < N;i++) {
		vli_umult2<N>(s, x, y[i]);
		vli_add_to<N + 1>(r, s);
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		u = vli_add_to<N + 1>(r, s);
		vli_rshift1w<N + 1>(r, u);
	}
	vli_mod<N>(result, r, prime, r[N] != 0);
}

template<const uint N> forceinline static void
vli_mont_sqr(u64 *result, const u64 *x, const u64 *prime, const u64 k0) noexcept
{
	u64	s[N + 1];
	u64	r[N + 1];
	vli_clear<N + 1>(r);
	for (uint i=0; i < N;i++) {
		vli_umult2<N>(s, x, x[i]);
		vli_add_to<N + 1>(r, s);
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		u = vli_add_to<N + 1>(r, s);
		vli_rshift1w<N + 1>(r, u);
	}
	vli_mod<N>(result, r, prime, r[N] != 0);
}


template<const uint N, const u64 k0> forceinline
static void mont_reduction(u64 *result, const u64 *y, const u64 *prime) noexcept
{
	u64	s[N + 1];
	u64	r[N];
	vli_set<N>(r, y);
	for (uint i=0; i < N; i++) {
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		s[N] += vli_add_to<N>(r, s);
		vli_rshift1w<N>(r, s[N]);
	}
	vli_mod<N>(result, r, prime);
}

template<const uint N, const u64 k0> forceinline static void
mont_mult(u64 *result, const u64 *x, const u64 *y, const u64 *prime) noexcept
{
	u64	s[N + 1];
	u64	r[N + 1];
	vli_clear<N + 1>(r);
	for (uint i=0; i < N;i++) {
		vli_umult2<N>(s, x, y[i]);
		vli_add_to<N + 1>(r, s);
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		u = vli_add_to<N + 1>(r, s);
		vli_rshift1w<N + 1>(r, u);
	}
	vli_mod<N>(result, r, prime, r[N] != 0);
}


// secp256k1 prime optimize
// p is 2^256 - 2^32 - 0x3d1 = 2^256 - 0x1000003d1
/* result > mod (result = mod + remainder), so subtract mod to
 * get remainder.
 */
static forceinline void btcp_mod(u64 *res, const u64 *product) noexcept
{
	u64 c = -secp256k1_p[0];
	u64 t[4 + 1];
	u64 r[4 + 2];
	vli_set<4>(r, product);
	r[4] = 0;
	vli_umult2<4>(t, product+4, c);
	r[5] = vli_add_to<5>(r, t);
	vli_clear<4>(t);
	vli_umult2<2>(t, r+4, c);
	r[4] = vli_add_to<4>(r, t);
	if (r[4]) {
		vli_sub<4>(res, r, secp256k1_p);
		return;
	}
	auto carry = vli_sub<4>(res, r, secp256k1_p);
	if (carry) vli_set<4>(res, r);
}


forceinline static void btc_reduction(u64 *result, const u64 *y) noexcept
{
	u64	s[4 + 1];
	u64	r[4];
	vli_set<4>(r, y);
	for (uint i=0; i < 4; i++) {
		u64	u = r[0] * secp256k1_p_k0;
		vli_btc_multP(s, u);
		s[4] += vli_add_to<4>(r, s);
		vli_rshift1w<4>(r, s[4]);
	}
	vli_mod<4>(result, r, secp256k1_p);
}

forceinline
static void btc_mont_mult(u64 *result, const u64 *x, const u64 *y) noexcept
{
	u64	s[4 + 1];
	u64	r[4 + 1];
	vli_clear<4 + 1>(r);
	for (uint i=0; i < 4;i++) {
		vli_umult2<4>(s, x, y[i]);
		vli_add_to<4 + 1>(r, s);
		u64	u = r[0] * secp256k1_p_k0;
		vli_btc_multP(s, u);
		u = vli_add_to<4 + 1>(r, s);
		vli_rshift1w<4 + 1>(r, u);
	}
	vli_mod<4>(result, r, secp256k1_p, r[4] != 0);
}


#ifdef	WITH_SM2_MULTSTEP
forceinline static void
sm2p_multStep(u64& r0, u64& r1, u64& r2, u64& r3, u64& r4, const u64& x0,
		const u64& x1, const u64& x2, const u64& x3, const u64 yi) noexcept
{
	uint128_t	pd;
	u64			cc, t0;	// t0 never overflow,  0xf * 0xf = 225 + 30 ... 255(ff)
	pd.mul_64_64(x0, yi);
	cc = 0;
	r0 = u64_addc(r0, pd.m_low(), cc);
	t0 = u64_addcz(pd.m_high(), cc);
	pd.mul_64_64(x1, yi);
	//cc = 0;
	r1 = u64_addc(r1, t0, cc);
	t0 = u64_addcz(pd.m_high(), cc);
	//cc = 0;
	r1 = u64_addc(r1, pd.m_low(), cc);
	t0 = u64_addcz(t0, cc);
	pd.mul_64_64(x2, yi);
	//cc = 0;
	r2 = u64_addc(r2, t0, cc);
	t0 = u64_addcz(pd.m_high(), cc);
	//cc = 0;
	r2 = u64_addc(r2, pd.m_low(), cc);
	t0 = u64_addcz(t0, cc);
	pd.mul_64_64(x3, yi);
	//cc = 0;
	r3 = u64_addc(r3, t0, cc);
	t0 = u64_addcz(pd.m_high(), cc);
	//cc = 0;
	r3 = u64_addc(r3, pd.m_low(), cc);
	r4 = u64_addc(r4, t0, cc);
}
#else
// r0..3 = r0..3 + s[0..3]
// return s[4] + carry (of above)
forceinline static u64
vli4_addc_to(u64& r0, u64& r1, u64& r2, u64& r3, const u64* s,
	const u64 carry=0) noexcept
{
	u64	cc = 0;
	r0 = u64_addc(r0, s[0], cc);
	r1 = u64_addc(r1, s[1], cc);
	r2 = u64_addc(r2, s[2], cc);
	r3 = u64_addc(r3, s[3], cc);
	return u64_addc(carry, s[4], cc);
}
#endif

forceinline
static void sm2p_multN(u64 *result, const u64 *x, const u64 *y) noexcept
{
#ifdef	WITH_SM2_MULTSTEP
	u64	r0=0, r1=0, r2=0, r3=0;
	u64	x0=x[0], x1=x[1], x2=x[2], x3=x[3];
	u64	carry=0;
	{
		sm2p_multStep(r0, r1, r2, r3, carry, x0, x1, x2, x3, y[0]);
		sm2p_reductionStep(r0, r1, r2, r3, carry);
		sm2p_multStep(r1, r2, r3, r0, carry, x0, x1, x2, x3, y[1]);
		sm2p_reductionStep(r1, r2, r3, r0, carry);
		sm2p_multStep(r2, r3, r0, r1, carry, x0, x1, x2, x3, y[2]);
		sm2p_reductionStep(r2, r3, r0, r1, carry);
		sm2p_multStep(r3, r0, r1, r2, carry, x0, x1, x2, x3, y[3]);
		sm2p_reductionStep(r3, r0, r1, r2, carry);
	}
#else
	u64	r0, r1, r2, r3;
	u64	carry;
	{
		u64	s[4+1];
		vli_umult2<4>(s, x, y[0]);
		vli4_load(s, r0, r1, r2, r3);
		carry = s[4];
		sm2p_reductionStep(r0, r1, r2, r3, carry);
		vli_umult2<4>(s, x, y[1]);
		carry = vli4_addc_to(r1, r2, r3, r0, s, carry);
		sm2p_reductionStep(r1, r2, r3, r0, carry);
		vli_umult2<4>(s, x, y[2]);
		carry = vli4_addc_to(r2, r3, r0, r1, s, carry);
		sm2p_reductionStep(r2, r3, r0, r1, carry);
		vli_umult2<4>(s, x, y[3]);
		carry = vli4_addc_to(r3, r0, r1, r2, s, carry);
		sm2p_reductionStep(r3, r0, r1, r2, carry);
	}
#endif
	// sm2p_mod
	{
		u64	cc=0;
		u64 s0 = u64_subc(r0, sm2_p[0], cc);
		u64 s1 = u64_subc(r1, sm2_p[1], cc);
		u64 s2 = u64_subc(r2, sm2_p[2], cc);
		u64 s3 = u64_subc(r3, sm2_p[3], cc);
		u64_subcz(carry, cc);
		if (cc != 0) vli4_save(result, r0, r1, r2, r3); else
			vli4_save(result, s0, s1, s2, s3);
	}
}


forceinline
static void sm2p_mult(u64 *result, const u64 *x, const u64 *y) noexcept
{
#if	defined(__x86_64__) || defined(__aarch64__)
	sm2p_mult_asm(result, x, y);
#else
	sm2p_multN(result, x, y);
#endif
}

template<const uint N, const u64 k0> forceinline
static void mont_sqr(u64 *result, const u64 *x, const u64 *prime) noexcept
{
	u64	s[N + 1];
	u64	r[N + 1];
	vli_clear<N + 1>(r);
	for (uint i=0; i < N;i++) {
		vli_umult2<N>(s, x, x[i]);
		vli_add_to<N + 1>(r, s);
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		u = vli_add_to<N + 1>(r, s);
		vli_rshift1w<N + 1>(r, u);
	}
	vli_mod<N>(result, r, prime, r[N] != 0);
}

template<const uint N, const u64 k0> forceinline
static void mont_sqrN(u64 *result, const u64 *x, const u64 *prime) noexcept
{
	u64	r[N * 2];
	u64	s[N + 1];
	vli_square<N>(r, x);
	// reduction r[0..N-1]
	for (uint i=0; i < N; i++) {
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		s[N] += vli_add_to<N>(r, s);
		vli_rshift1w<N>(r, s[N]);
	}
	// add r[N...2N-1], then mod P
	{
		s[N] = vli_add_to<N>(r, r + N);
		vli_mod<N>(result, r, prime, s[N] != 0);
	}
}

forceinline static void btc_sqrN(u64 *result, const u64 *x) noexcept
{
	u64	r[4 * 2];
	vli_square<4>(r, x);
	btcp_mod(result, r);
}

forceinline static void sm2p_sqrN(u64 *result, const u64 *x) noexcept
{
#if	defined(__x86_64__) || defined(__aarch64__)
	sm2p_sqr_asm(result, x);
#else
	u64	r[8];
	vli_square<4>(r, x);
	sm2p_reductionProd(result, r);
#endif
}

#endif	//	__MONT_HPP__
