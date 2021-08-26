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
#ifndef __MONT_ASM_H__
#define __MONT_ASM_H__

#include "cdefs.h"

#ifdef	__cplusplus
extern "C" {
#endif

// u * 2^256 mod sm2 prime
// p is 2^256 - 2^224 - 2^96 + 2^64 -1
// R(2^256) = 2^224 + 2^96 - 2^64 + 1 Mod p = 2^32 * 2^192 + (2^32-1) * 2^64 + 1
// u < 2^32
// used for carry_reduce
forceinline
static void vli_sm2_multR(u64 *result, const u64 uv)
{
	u64 u = uv & 0xffffffff;
	result[0] = u;
	result[1] = (u << 32) - u;
	result[2] = 0;
	result[3] = u << 32;
	return;
}


/* result > mod (result = mod + remainder), so subtract mod to
 * get remainder.
 */
static forceinline void
sm2p_mod_asm(u64 *res, const u64 *left, const bool carry)
{
#ifdef	__x86_64__
	asm volatile(
		"movq (%%rsi), %%r8\n"	// mov r8/9/10/11, right
		"movq 8(%%rsi), %%r9\n"
		"movq 16(%%rsi), %%r10\n"
		"movq 24(%%rsi), %%r11\n"
		"movq %%r8, %%r12\n"
		"movq %%r9, %%r13\n"
		"movq %%r10, %%r14\n"
		"movq %%r11, %%r15\n"
		"subq $-1, %%r12\n"
		"sbbq %[mod0], %%r13\n"
		"sbbq $-1, %%r14\n"
		"sbbq %[mod1], %%r15\n"
		"sbbq $0, %%rax\n"
		"cmovcq %%r8, %%r12\n"
		"cmovcq %%r9, %%r13\n"
		"cmovcq %%r10, %%r14\n"
		"cmovcq %%r11, %%r15\n"
		"movq %%r12, (%%rdi)\n"
		"movq %%r13, 8(%%rdi)\n"
		"movq %%r14, 16(%%rdi)\n"
		"movq %%r15, 24(%%rdi)\n"
		:
		: "S"(left), "D"(res), "a" ((u64)carry), [mod0] "m" (sm2_p[1]),
			[mod1] "m" (sm2_p[3])
		: "%r8", "%r9", "%r10", "%r11" , "%r12", "%r13", "%r14", "%r15", "cc", "memory");
#elif	defined(__aarch64__)
	asm volatile(
		//"mov x9, #-1\n"
		//"mov x10, %3\n"
		"mov x11, #-1\n"
		//"mov x12, %4\n"
		"ldp x4, x5, [%2]\n"
		"ldp x6, x7, [%2, 16]\n"
		"subs x9, x4, #-1\n"
		"sbcs x10, x5, %3\n"
		"sbcs x11, x6, x11\n"
		"sbcs x12, x7, %4\n"
		"sbcs %0, %0, xzr\n"
		"csel x4, x4, x9, cc\n"
		"csel x5, x5, x10, cc\n"
		"csel x6, x6 , x11, cc\n"
		"csel x7, x7, x12, cc\n"
		"stp x4, x5, [%1]\n"
		"stp x6, x7, [%1, 16]\n"
		//"adc %0, xzr, xzr\n"
		:
		: "r" ((u64)carry), "r" (res), "r" (left), "r" (sm2_p[1]), "r" (sm2_p[3])
		: "%x4", "%x5", "%x6", "%x7", "%x9", "%x10", "%x11", "%x12", "cc", "memory");
#endif
}


forceinline static void
sm2p_reductionS0(u64& r0, u64& r1, u64& r2, u64& r3)
{
#ifdef	__x86_64__
	register u64 res0 asm("r8")=r0;
	register u64 res1 asm("r9")=r1;
	register u64 res2 asm("r10")=r2;
	register u64 res3 asm("r11")=r3;
	asm volatile(
		//"MOVQ (8*0)(%%rsi), %%r8\n"
		//"MOVQ (8*1)(%%rsi), %%r9\n"
		//"MOVQ (8*2)(%%rsi), %%r10\n"
		//"MOVQ (8*3)(%%rsi), %%r11\n"

	// Only reduce, no multiplications are needed
	// First stage
		"MOVQ %%r8, %%rax\n"
		"MOVQ %%r8, %%rdx\n"
		"SHLQ $32, %%rax\n"
		"SHRQ $32, %%rdx\n"
		"ADDQ %%r8, %%r9\n"
		"ADCQ $0, %%r10\n"
		"ADCQ $0, %%r11\n"
		"ADCQ $0, %%r8\n"
		"subq %%rax, %%r9\n"
		"sbbq %%rdx, %%r10\n"
		"sbbq %%rax, %%r11\n"
		"sbbq %%rdx, %%r8\n"
	// Second stage
		"MOVQ %%r9, %%rax\n"
		"MOVQ %%r9, %%rdx\n"
		"SHLQ $32, %%rax\n"
		"SHRQ $32, %%rdx\n"
		"ADDQ %%r9, %%r10\n"
		"adcq $0, %%r11\n"
		"adcq $0, %%r8\n"
		"adcq $0, %%r9\n"
		"subq %%rax, %%r10\n"
		"sbbq %%rdx, %%r11\n"
		"sbbq %%rax, %%r8\n"
		"sbbq %%rdx, %%r9\n"
	// Third stage
		"MOVQ %%r10, %%rax\n"
		"MOVQ %%r10, %%rdx\n"
		"SHLQ $32, %%rax\n"
		"SHRQ $32, %%rdx\n"
		"ADDQ %%r10, %%r11\n"
		"adcq $0, %%r8\n"
		"adcq $0, %%r9\n"
		"adcq $0, %%r10\n"
		"subq %%rax, %%r11\n"
		"sbbq %%rdx, %%r8\n"
		"sbbq %%rax, %%r9\n"
		"sbbq %%rdx, %%r10\n"
	// Last stage
		"MOVQ %%r11, %%rax\n"
		"MOVQ %%r11, %%rdx\n"
		"SHLQ $32, %%rax\n"
		"SHRQ $32, %%rdx\n"
		"ADDQ %%r11, %%r8\n"
		"adcq $0, %%r9\n"
		"adcq $0, %%r10\n"
		"adcq $0, %%r11\n"
		"subq %%rax, %%r8\n"
		"sbbq %%rdx, %%r9\n"
		"sbbq %%rax, %%r10\n"
		"sbbq %%rdx, %%r11\n"
		: "+r" (res0), "+r" (res1), "+r" (res2), "+r" (res3)
		: //"S" (y)
		: "rax", "rdx", "cc");
	r0 = res0;
	r1 = res1;
	r2 = res2;
	r3 = res3;
#elif	defined(__aarch64__)
	register u64 res0 asm("x4") = r0;
	register u64 res1 asm("x5") = r1;
	register u64 res2 asm("x6") = r2;
	register u64 res3 asm("x7") = r3;
	asm volatile(
	// Only reduce, no multiplications are needed
	// First reduction step
		"LSR x13, x4, 32\n"
		"LSL x14, x4, 32\n"
		"ADDS x5, x5, x4\n"
		"ADCS x6, x6, XZR\n"
		"ADCS x7, x7, XZR\n"
		"ADCS x4, x4, xzr\n"
		"ADCS x0, XZR, XZR\n"
		"SUBS x5, x5, x14\n"
		"SBCS x6, x6, x13\n"
		"SBCS x7, x7, x14\n"
		"SBCS x4, x4, x13\n"
		"SBCS x0, x0, xzr\n"
	// Second reduction step
		"LSR x13, x5, 32\n"
		"LSL x14, x5, 32\n"
		"ADDS x6, x6, x5\n"
		"ADCS x7, x7, XZR\n"
		"ADCS x4, x4, XZR\n"
		"ADCS x5, x5, xzr\n"
		"ADCS x0, XZR, XZR\n"
		"SUBS x6, x6, x14\n"
		"SBCS x7, x7, x13\n"
		"SBCS x4, x4, x14\n"
		"SBCS x5, x5, x13\n"
		"SBCS x0, x0, xzr\n"
	// Third reduction step
		"LSR x13, x6, 32\n"
		"LSL x14, x6, 32\n"
		"ADDS x7, x7, x6\n"
		"ADCS x4, x4, XZR\n"
		"ADCS x5, x5, XZR\n"
		"ADCS x6, x6, xzr\n"
		"ADCS x0, XZR, XZR\n"
		"SUBS x7, x7, x14\n"
		"SBCS x4, x4, x13\n"
		"SBCS x5, x5, x14\n"
		"SBCS x6, x6, x13\n"
		"SBCS x0, x0, xzr\n"
	// Last reduction step
		"LSR x13, x7, 32\n"
		"LSL x14, x7, 32\n"
		"ADDS x4, x4, x7\n"
		"ADCS x5, x5, XZR\n"
		"ADCS x6, x6, XZR\n"
		"ADCS x7, x7, xzr\n"
		"ADCS x0, XZR, XZR\n"
		"SUBS x4, x4, x14\n"
		"SBCS x5, x5, x13\n"
		"SBCS x6, x6, x14\n"
		"SBCS x7, x7, x13\n"
		"SBCS x0, x0, xzr\n"
		: "+r" (res0), "+r" (res1), "+r" (res2), "+r" (res3)
		:
		: "x0", "x13", "x14", "cc");
	r0 = res0;
	r1 = res1;
	r2 = res2;
	r3 = res3;
#endif
}

forceinline static void sm2p_reductionS2(u64 *result, const u64 r0,
	const u64 r1, const u64 r2, const u64 r3, const u64 carry=0)
{
#ifdef	__x86_64__
	register u64 res0 asm("r8") = r0;
	register u64 res1 asm("r9") = r1;
	register u64 res2 asm("r10") = r2;
	register u64 res3 asm("r11") = r3;
	// mod prime
	asm volatile(
		"MOVQ %%r8, %%r12\n"
		"MOVQ %%r9, %%r13\n"
		"MOVQ %%r10, %%r14\n"
		"MOVQ %%r11, %%r15\n"

		"SUBQ $-1, %%r8\n"
		"SBBQ %[pr1], %%r9\n"
		"SBBQ $-1, %%r10\n"
		"SBBQ %[pr3], %%r11\n"
		"sbbq $0, %%rax\n"

		"CMOVCQ %%r12, %%r8\n"
		"CMOVCQ %%r13, %%r9\n"
		"CMOVCQ %%r14, %%r10\n"
		"CMOVCQ %%r15, %%r11\n"

		"MOVQ %%r8, (8*0)(%%rdi)\n"
		"MOVQ %%r9, (8*1)(%%rdi)\n"
		"MOVQ %%r10, (8*2)(%%rdi)\n"
		"MOVQ %%r11, (8*3)(%%rdi)\n"
		: 			// acc4/5/0/1
		: "D" (result), "a" (carry), [pr1] "m" (sm2_p[1]),
		[pr3] "m" (sm2_p[3]), "r" (res0), "r" (res1), "r" (res2), "r" (res3)
		: "r12", "r13", "r14", "r15", "cc", "memory");
#elif	defined(__aarch64__)
	register u64 res0 asm("x4") = r0;
	register u64 res1 asm("x5") = r1;
	register u64 res2 asm("x6") = r2;
	register u64 res3 asm("x7") = r3;
	register u64 _carry asm("x0") = carry;
	// mod prime
	asm volatile(
		//"MOV	x9, #-1\n"
		"MOV	x11, #-1\n"

		"SUBS	x9, x4, #-1\n"
		"SBCS	x10, x5, %1\n"
		"SBCS	x11, x6, x11\n"
		"SBCS	x12, x7, %2\n"
		"SBCS	x0, x0, XZR\n"

		"CSEL	x4, x4, x9, cc\n"
		"CSEL	x5, x5, x10, cc\n"
		"CSEL	x6, x6, x11, cc\n"
		"CSEL	x7, x7, x12, cc\n"

		"STP	x4, x5, [%0]\n"
		"STP	x6, x7, [%0, 16]\n"
		:
		: "r" (result), "r" (sm2_p[1]), "r" (sm2_p[3]), "r" (_carry),
			"r" (res0), "r" (res1), "r" (res2), "r" (res3)
		: "%x9", "%x10", "%x11", "%x12", "cc", "memory");
#endif
}


forceinline static void sm2p_mult_asm(u64 *result, const u64 *x, const u64 *y)
{
#ifdef	__x86_64__
	asm volatile(
	// x * y[0]
		"MOVQ (8*0)(%0), %%r14\n"

		"MOVQ (8*0)(%%rsi), %%RAX\n"
		"MULQ %%r14\n"
		"MOVQ %%RAX, %%r8\n"
		"MOVQ %%RDX, %%r9\n"

		"MOVQ (8*1)(%%rsi), %%RAX\n"
		"MULQ %%r14\n"
		"ADDQ %%RAX, %%r9\n"
		"ADCQ $0, %%RDX\n"
		"MOVQ %%RDX, %%r10\n"

		"MOVQ (8*2)(%%rsi), %%RAX\n"
		"MULQ %%r14\n"
		"ADDQ %%RAX, %%r10\n"
		"ADCQ $0, %%RDX\n"
		"MOVQ %%RDX, %%r11\n"

		"MOVQ (8*3)(%%rsi), %%RAX\n"
		"MULQ %%r14\n"
		"ADDQ %%RAX, %%r11\n"
		"ADCQ $0, %%RDX\n"
		"MOVQ %%RDX, %%r12\n"
		"XORQ %%r13, %%r13\n"
	// First reduction step
		"MOVQ %%r8, %%RAX\n"
		"MOVQ %%r8, %%r15\n"
		"SHLQ $32, %%r8\n"
		"SHRQ $32, %%r15\n"
		"ADDQ %%rax, %%r9\n"
		"ADCQ $0, %%r10\n"
		"ADCQ $0, %%r11\n"
		"ADCQ %%rax, %%r12\n"
		"ADCQ $0, %%r13\n"
		"SUBQ %%r8, %%r9\n"
		"SBBQ %%r15, %%r10\n"
		"SBBQ %%r8, %%r11\n"
		"SBBQ %%R15, %%r12\n"
		"SBBQ $0, %%r13\n"
		"XORQ %%r8, %%r8\n"
	// x * y[1]
		"MOVQ (8*1)(%0), %%r14\n"

		"MOVQ (8*0)(%%rsi), %%RAX\n"
		"MULQ %%r14\n"
		"ADDQ %%RAX, %%r9\n"
		"ADCQ $0, %%RDX\n"
		"MOVQ %%RDX, %%r15\n"

		"MOVQ (8*1)(%%rsi), %%RAX\n"
		"MULQ %%r14\n"
		"ADDQ %%r15, %%r10\n"
		"ADCQ $0, %%RDX\n"
		"ADDQ %%RAX, %%r10\n"
		"ADCQ $0, %%RDX\n"
		"MOVQ %%RDX, %%r15\n"

		"MOVQ (8*2)(%%rsi), %%RAX\n"
		"MULQ %%r14\n"
		"ADDQ %%r15, %%r11\n"
		"ADCQ $0, %%RDX\n"
		"ADDQ %%RAX, %%r11\n"
		"ADCQ $0, %%RDX\n"
		"MOVQ %%RDX, %%r15\n"

		"MOVQ (8*3)(%%rsi), %%RAX\n"
		"MULQ %%r14\n"
		"ADDQ %%r15, %%r12\n"
		"ADCQ $0, %%RDX\n"
		"ADDQ %%RAX, %%r12\n"
		"ADCQ %%RDX, %%r13\n"
		"ADCQ $0, %%r8\n"
	// Second reduction step
		"MOVQ %%r9, %%RAX\n"
		"MOVQ %%r9, %%r15\n"
		"SHLQ $32, %%r9\n"
		"SHRQ $32, %%r15\n"
		"ADDQ %%rax, %%r10\n"
		"ADCQ $0, %%r11\n"
		"ADCQ $0, %%r12\n"
		"ADCQ %%rax, %%r13\n"
		"ADCQ $0, %%r8\n"
		"SUBQ %%r9, %%r10\n"
		"SBBQ %%r15, %%r11\n"
		"SBBQ %%r9, %%r12\n"
		"SBBQ %%r15, %%r13\n"
		"SBBQ $0, %%r8\n"
		"XORQ %%r9, %%r9\n"
	// x * y[2]
		"MOVQ (8*2)(%0), %%r14\n"

		"MOVQ (8*0)(%%rsi), %%RAX\n"
		"MULQ %%r14\n"
		"ADDQ %%RAX, %%r10\n"
		"ADCQ $0, %%RDX\n"
		"MOVQ %%RDX, %%r15\n"

		"MOVQ (8*1)(%%rsi), %%RAX\n"
		"MULQ %%r14\n"
		"ADDQ %%r15, %%r11\n"
		"ADCQ $0, %%RDX\n"
		"ADDQ %%RAX, %%r11\n"
		"ADCQ $0, %%RDX\n"
		"MOVQ %%RDX, %%r15\n"

		"MOVQ (8*2)(%%rsi), %%RAX\n"
		"MULQ %%r14\n"
		"ADDQ %%r15, %%r12\n"
		"ADCQ $0, %%RDX\n"
		"ADDQ %%RAX, %%r12\n"
		"ADCQ $0, %%RDX\n"
		"MOVQ %%RDX, %%r15\n"

		"MOVQ (8*3)(%%rsi), %%RAX\n"
		"MULQ %%r14\n"
		"ADDQ %%r15, %%r13\n"
		"ADCQ $0, %%RDX\n"
		"ADDQ %%RAX, %%r13\n"
		"ADCQ %%RDX, %%r8\n"
		"ADCQ $0, %%r9\n"
	// Third reduction step
		"MOVQ %%r10, %%RAX\n"
		"MOVQ %%r10, %%r15\n"
		"SHLQ $32, %%r10\n"
		"SHRQ $32, %%r15\n"
		"ADDQ %%rax, %%r11\n"
		"ADCQ $0, %%r12\n"
		"ADCQ $0, %%r13\n"
		"ADCQ %%rax, %%r8\n"
		"ADCQ $0, %%r9\n"
		"SUBQ %%r10, %%r11\n"
		"SBBQ %%r15, %%r12\n"
		"SBBQ %%r10, %%r13\n"
		"SBBQ %%r15, %%r8\n"
		"SBBQ $0, %%r9\n"
		"XORQ %%r10, %%r10\n"
	// x * y[3]
		"MOVQ (8*3)(%0), %%r14\n"

		"MOVQ (8*0)(%%rsi), %%RAX\n"
		"MULQ %%r14\n"
		"ADDQ %%RAX, %%r11\n"
		"ADCQ $0, %%RDX\n"
		"MOVQ %%RDX, %%r15\n"

		"MOVQ (8*1)(%%rsi), %%RAX\n"
		"MULQ %%r14\n"
		"ADDQ %%r15, %%r12\n"
		"ADCQ $0, %%RDX\n"
		"ADDQ %%RAX, %%r12\n"
		"ADCQ $0, %%RDX\n"
		"MOVQ %%RDX, %%r15\n"

		"MOVQ (8*2)(%%rsi), %%RAX\n"
		"MULQ %%r14\n"
		"ADDQ %%r15, %%r13\n"
		"ADCQ $0, %%RDX\n"
		"ADDQ %%RAX, %%r13\n"
		"ADCQ $0, %%RDX\n"
		"MOVQ %%RDX, %%r15\n"

		"MOVQ (8*3)(%%rsi), %%RAX\n"
		"MULQ %%r14\n"
		"ADDQ %%r15, %%r8\n"
		"ADCQ $0, %%RDX\n"
		"ADDQ %%RAX, %%r8\n"
		"ADCQ %%RDX, %%r9\n"
		"ADCQ $0, %%r10\n"
	// Last reduction step
		"MOVQ %%r11, %%RAX\n"
		"MOVQ %%r11, %%r15\n"
		"SHLQ $32, %%r11\n"
		"SHRQ $32, %%r15\n"
		"ADDQ %%rax, %%r12\n"
		"ADCQ $0, %%r13\n"
		"ADCQ $0, %%r8\n"
		"ADCQ %%rax, %%r9\n"
		"ADCQ $0, %%r10\n"
		"SUBQ %%r11, %%r12\n"
		"SBBQ %%r15, %%r13\n"
		"SBBQ %%r11, %%r8\n"
		"SBBQ %%r15, %%r9\n"
		"SBB $0, %%r10\n"
	// Copy result [255:0]
		"MOVQ %%r12, %%rax\n"
		"MOVQ %%r13, %%r11\n"
		"MOVQ %%r8, %%r14\n"
		"MOVQ %%r9, %%r15\n"
	// Subtract sm2_p
		"SUBQ $-1, %%r12\n"
		"SBBQ %[pr1] ,%%r13\n"
		"SBBQ $-1, %%r8\n"
		"SBBQ %[pr3], %%r9\n"
		"SBBQ $0, %%r10\n"

		"CMOVCQ %%rax, %%r12\n"
		"CMOVCQ %%r11, %%r13\n"
		"CMOVCQ %%r14, %%r8\n"
		"CMOVCQ %%r15, %%r9\n"

		"MOVQ %%r12, (8*0)(%%rdi)\n"
		"MOVQ %%r13, (8*1)(%%rdi)\n"
		"MOVQ %%r8, (8*2)(%%rdi)\n"
		"MOVQ %%r9, (8*3)(%%rdi)\n"
		: 			// acc4/5/0/1
		: "r" (y), "D" (result), "S" (x), [pr1] "m" (sm2_p[1]),
		[pr3] "m" (sm2_p[3])
		: "rax", "rdx", "r8", "r9", "r10", "r11", "r12", "r13", "r14",
		"r15", "cc", "memory");
#elif	defined(__aarch64__)
	// x0 -- x3   register %%x4 -- %%x7
	register u64 x0 asm("x4") = x[0];
	register u64 x1 asm("x5") = x[1];
	register u64 x2 asm("x6") = x[2];
	register u64 x3 asm("x7") = x[3];
	asm volatile(
	// y[0] * x
	// load y0, y1
		"LDR x3, [%1]\n"
		"MUL	x9, x3, x4\n"
		"UMULH	x10, x3, x4\n"

		"MUL	x14, x3, x5\n"
		"ADDS	x10, x14, x10\n"
		"UMULH	x11, x3, x5\n"

		"MUL	x14, x3, x6\n"
		"ADCS	x11, x14, x11\n"
		"UMULH	x12, x3, x6\n"

		"MUL	x14, x3, x7\n"
		"ADCS	x12, x14, x12\n"
		"UMULH	x13, x3, x7\n"
		"ADC	x13, xzr, x13\n"
	// First reduction step
		"LSR	x14, x9, #32\n"
		"LSL	x15, x9, #32\n"
		"ADDS	x10, x10, x9\n"
		"ADCS	x11, x11, xzr\n"
		"ADCS	x12, x12, xzr\n"
		"ADCS	x13, x13, x9\n"
		"ADC	x9, xzr, xzr\n"
		"SUBS	x10, x10, x15\n"
		"SBCS	x11, x11, x14\n"
		"SBCS	x12, x12, x15\n"
		"SBCS	x13, x13, x14\n"
		"SBC	x9, x9, xzr\n"
	// y[1] * x
		"LDR x3, [%1, 8]\n"
		"MUL	x14, x3, x4\n"
		"ADDS	x10, x14, x10\n"
		"UMULH	x15, x3, x4\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x5\n"
		"ADDS	x11, x15, x11\n"
		"UMULH	x15, x3, x5\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x11, x14, x11\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x6\n"
		"ADDS	x12, x15, x12\n"
		"UMULH	x15, x3, x6\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x12, x14, x12\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x7\n"
		"ADDS	x13, x15, x13\n"
		"UMULH	x15, x3, x7\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x13, x14, x13\n"
		"ADC	x9, x15, x9\n"

	// Second reduction step
		"LSR	x14, x10, 32\n"
		"LSL	x15, x10, 32\n"
		"ADDS	x11, x11, x10\n"
		"ADCS	x12, x12, xzr\n"
		"ADCS	x13, x13, xzr\n"
		"ADCS	x9, x9, x10\n"
		"ADC	x10, xzr, xzr\n"
		"SUBS	x11, x11, x15\n"
		"SBCS	x12, x12, x14\n"
		"SBCS	x13, x13, x15\n"
		"SBCS	x9, x9, x14\n"
		"SBC	x10, x10, xzr\n"
	// y[2] * x
	// load y2, y3
		"LDR x3, [%1, 16]\n"
		"MUL	x14, x3, x4\n"
		"ADDS	x11, x14, x11\n"
		"UMULH	x15, x3, x4\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x5\n"
		"ADDS	x12, x15, x12\n"
		"UMULH	x15, x3, x5\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x12, x14, x12\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x6\n"
		"ADDS	x13, x15, x13\n"
		"UMULH	x15, x3, x6\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x13, x14, x13\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x7\n"
		"ADDS	x9, x15, x9\n"
		"UMULH	x15, x3, x7\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x9, x14, x9\n"
		"ADC	x10, x15, x10\n"

	// Third reduction step
		"LSR	x14, x11, 32\n"
		"LSL	x15, x11, 32\n"
		"ADDS	x12, x12, x11\n"
		"ADCS	x13, x13, xzr\n"
		"ADCS	x9, x9, xzr\n"
		"ADCS	x10, x10, x11\n"
		"ADC	x11, xzr, xzr\n"
		"SUBS	x12, x12, x15\n"
		"SBCS	x13, x13, x14\n"
		"SBCS	x9, x9, x15\n"
		"SBCS	x10, x10, x14\n"
		"SBC	x11, x11, xzr\n"
	// y[3] * x
		"LDR x3, [%1, 24]\n"
		"MUL	x14, x3, x4\n"
		"ADDS	x12, x14, x12\n"
		"UMULH	x15, x3, x4\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x5\n"
		"ADDS	x13, x15, x13\n"
		"UMULH	x15, x3, x5\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x13, x14, x13\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x6\n"
		"ADDS	x9, x15, x9\n"
		"UMULH	x15, x3, x6\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x9, x14, x9\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x7\n"
		"ADDS	x10, x15, x10\n"
		"UMULH	x15, x3, x7\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x10, x14, x10\n"
		"ADC	x11, x15, x11\n"

	// Last reduction step
		"LSR	x14, x12, 32\n"
		"LSL	x15, x12, 32\n"
		"ADDS	x13, x13, x12\n"
		"ADCS	x9, x9, xzr\n"
		"ADCS	x10, x10, xzr\n"
		"ADCS	x11, x11, x12\n"
		"ADC	x12, xzr, xzr\n"
		"SUBS	x13, x13, x15\n"
		"SBCS	x9, x9, x14\n"
		"SBCS	x10, x10, x15\n"
		"SBCS	x11, x11, x14\n"
		"SBC	x12, x12, xzr\n"

		"ldp	x4, x5, [%2]\n"
		"ldp	x6, x7, [%2, 16]\n"
		"SUBS	x4, x13, x4\n"
		"SBCS	x5, x9, x5\n"
		"SBCS	x6, x10, x6\n"
		"SBCS	x7, x11, x7\n"
		"SBCS	x12, x12, xzr\n"

		"CSEL	x4, x4, x13, cs\n"
		"CSEL	x5, x5, x9, cs\n"
		"CSEL	x6, x6, x10, cs\n"
		"CSEL	x7, x7, x11, cs\n"
		"stp	x4, x5, [%0]\n"
		"stp	x6, x7, [%0, 16]\n"
		:
		: "r" (result), "r" (y), "r" (sm2_p), "r" (x0), "r" (x1), "r" (x2),
		"r" (x3)
		: "%x3", "%x9", "%x10", "%x11", "%x12", "%x13", "%x14", "%x15",
		"cc", "memory");
#endif
}

forceinline static void sm2p_sqr_asm(u64 *result, const u64 *x)
{
#ifdef	__x86_64__
	asm volatile(
	// y[1:] * y[0]
	"MOVQ (8*0)(%%rsi), %%r14\n"

	"MOVQ (8*1)(%%rsi), %%rax\n"
	"MULQ %%r14\n"
	"MOVQ %%rax, %%r9\n"
	"MOVQ %%rdx, %%r10\n"

	"MOVQ (8*2)(%%rsi), %%rax\n"
	"MULQ %%r14\n"
	"ADDQ %%rax, %%r10\n"
	"ADCQ $0, %%rdx\n"
	"MOVQ %%rdx, %%r11\n"

	"MOVQ (8*3)(%%rsi), %%rax\n"
	"MULQ %%r14\n"
	"ADDQ %%rax, %%r11\n"
	"ADCQ $0, %%rdx\n"
	"MOVQ %%rdx, %%r12\n"
	// y[2:] * y[1]
	"MOVQ (8*1)(%%rsi), %%r14\n"

	"MOVQ (8*2)(%%rsi), %%rax\n"
	"MULQ %%r14\n"
	"ADDQ %%rax, %%r11\n"
	"ADCQ $0, %%rdx\n"
	"MOVQ %%rdx, %%r15\n"

	"MOVQ (8*3)(%%rsi), %%rax\n"
	"MULQ %%r14\n"
	"ADDQ %%r15, %%r12\n"
	"ADCQ $0, %%rdx\n"
	"ADDQ %%rax, %%r12\n"
	"ADCQ $0, %%rdx\n"
	"MOVQ %%rdx, %%r13\n"
	// y[3] * y[2]
	"MOVQ (8*2)(%%rsi), %%r14\n"

	"MOVQ (8*3)(%%rsi), %%rax\n"
	"MULQ %%r14\n"
	"ADDQ %%rax, %%r13\n"
	"ADCQ $0, %%rdx\n"
	"MOVQ %%rdx, %%rcx\n"
	"XORQ %%r15, %%r15\n"
	// *2
	"ADDQ %%r9, %%r9\n"
	"ADCQ %%r10, %%r10\n"
	"ADCQ %%r11, %%r11\n"
	"ADCQ %%r12, %%r12\n"
	"ADCQ %%r13, %%r13\n"
	"ADCQ %%rcx, %%rcx\n"
	"ADCQ $0, %%r15\n"
	// Missing products
	"MOVQ (8*0)(%%rsi), %%rax\n"
	"MULQ %%rax\n"
	"MOVQ %%rax, %%r8\n"
	"MOVQ %%rdx, %%r14\n"

	"MOVQ (8*1)(%%rsi), %%rax\n"
	"MULQ %%rax\n"
	"ADDQ %%r14, %%r9\n"
	"ADCQ %%rax, %%r10\n"
	"ADCQ $0, %%rdx\n"
	"MOVQ %%rdx, %%r14\n"

	"MOVQ (8*2)(%%rsi), %%rax\n"
	"MULQ %%rax\n"
	"ADDQ %%r14, %%r11\n"
	"ADCQ %%rax, %%r12\n"
	"ADCQ $0, %%rdx\n"
	"MOVQ %%rdx, %%r14\n"

	"MOVQ (8*3)(%%rsi), %%rax\n"
	"MULQ %%rax\n"
	"ADDQ %%r14, %%r13\n"
	"ADCQ %%rax, %%rcx\n"
	"ADCQ %%rdx, %%r15\n"
	"MOVQ %%r15, %%rbx\n"
	// First reduction step
	"MOVQ %%r8, %%rax\n"
	"MOVQ %%r8, %%rdx\n"
	"MOVQ %%r8, %%r15\n"
	"SHLQ $32, %%r8\n"
	"SHRQ $32, %%r15\n"
	"ADDQ %%rax, %%r9\n"
	"ADCQ $0, %%r10\n"
	"ADCQ $0, %%r11\n"
	"ADCQ $0, %%rdx\n"
	"SUBQ %%r8, %%r9\n"
	"SBBQ %%r15, %%r10\n"
	"SBBQ %%r8, %%r11\n"
	"SBBQ %%R15, %%rdx\n"
	"MOVQ %%rdx, %%r8\n"
	// Second reduction step
	"MOVQ %%r9, %%rax\n"
	"MOVQ %%r9, %%rdx\n"
	"MOVQ %%r9, %%r15\n"
	"SHLQ $32, %%r9\n"
	"SHRQ $32, %%r15\n"
	"ADDQ %%rax, %%r10\n"
	"ADCQ $0, %%r11\n"
	"ADCQ $0, %%r8\n"
	"ADCQ $0, %%rdx\n"
	"SUBQ %%r9, %%r10\n"
	"SBBQ %%r15, %%r11\n"
	"SBBQ %%r9, %%r8\n"
	"SBBQ %%r15, %%rdx\n"
	"MOVQ %%rdx, %%r9\n"
	// Third reduction step
	"MOVQ %%r10, %%rax\n"
	"MOVQ %%r10, %%rdx\n"
	"MOVQ %%r10, %%r15\n"
	"SHLQ $32, %%r10\n"
	"SHRQ $32, %%r15\n"
	"ADDQ %%rax, %%r11\n"
	"ADCQ $0, %%r8\n"
	"ADCQ $0, %%r9\n"
	"ADCQ $0, %%rdx\n"
	"SUBQ %%r10, %%r11\n"
	"SBBQ %%r15, %%r8\n"
	"SBBQ %%r10, %%r9\n"
	"SBBQ %%r15, %%rdx\n"
	"MOVQ %%rdx, %%r10\n"
	// Last reduction step
	"XORQ %%r14, %%r14\n"
	"MOVQ %%r11, %%rax\n"
	"MOVQ %%r11, %%rdx\n"
	"MOVQ %%r11, %%r15\n"
	"SHLQ $32, %%r11\n"
	"SHRQ $32, %%r15\n"
	"ADDQ %%rax, %%r8\n"
	"ADCQ $0, %%r9\n"
	"ADCQ $0, %%r10\n"
	"ADCQ $0, %%rdx\n"
	"SUBQ %%r11, %%r8\n"
	"SBBQ %%r15, %%r9\n"
	"SBBQ %%r11, %%r10\n"
	"SBBQ %%r15, %%rdx\n"
	"MOVQ %%rdx, %%r11\n"
	// Add bits [511:256] of the sqr result
	"ADCQ %%r12, %%r8\n"
	"ADCQ %%r13, %%r9\n"
	"ADCQ %%rcx, %%r10\n"
	"ADCQ %%rbx, %%r11\n"
	"ADCQ $0, %%r14\n"

	"MOVQ %%r8, %%r12\n"
	"MOVQ %%r9, %%r13\n"
	"MOVQ %%r10, %%rcx\n"
	"MOVQ %%r11, %%r15\n"
	// Subtract sm2_p
	"SUBQ $-1, %%r8\n"
	"SBBQ %[pr1] ,%%r9\n"
	"SBBQ $-1, %%r10\n"
	"SBBQ %[pr3], %%r11\n"
	"SBBQ $0, %%r14\n"

	"CMOVCQ %%r12, %%r8\n"
	"CMOVCQ %%r13, %%r9\n"
	"CMOVCQ %%rcx, %%r10\n"
	"CMOVCQ %%r15, %%r11\n"

	"MOVQ %%r8, (8*0)(%%rdi)\n"
	"MOVQ %%r9, (8*1)(%%rdi)\n"
	"MOVQ %%r10, (8*2)(%%rdi)\n"
	"MOVQ %%r11, (8*3)(%%rdi)\n"
		: 			// acc4/5/0/1
		: "D" (result), "S" (x), [pr1] "m" (sm2_p[1]), [pr3] "m" (sm2_p[3])
		: "rax", "rbx", "rcx", "rdx", "r8", "r9", "r12", "r13", "r10", "r11",
		"r14", "r15", "cc", "memory");
#elif	defined(__aarch64__)
	// x0 -- x3   registers for x[0] .. x[3]
	// x4 -- x7   registers for acc0 .. acc3
	// x9 -- x12  registers for acc4 .. acc7 
	register u64 x0 asm("x0") = x[0];
	register u64 x1 asm("x1") = x[1];
	register u64 x2 asm("x2") = x[2];
	register u64 x3 asm("x3") = x[3];
	register u64 acc0 asm("x4");
	register u64 acc1 asm("x5");
	register u64 acc2 asm("x6");
	register u64 acc3 asm("x7");
	register u64 acc4 asm("x9");
	register u64 acc5 asm("x10");
	register u64 acc6 asm("x11");
	register u64 acc7 asm("x12");
	asm volatile(
	// x[1:] * x[0]
		"MUL	x5, x0, x1\n"		// acc1
		"UMULH	x6, x0, x1\n"		// acc2

		"MUL	x13, x0, x2\n"
		"UMULH	x7, x0, x2\n"		//acc3

		"MUL	x14, x0, x3\n"
		"UMULH	x9, x0, x3\n"		//acc4
		"ADDS	x6, x13, x6\n"		//acc2
		"ADCS	x7, x14, x7\n"		//acc3
		"ADC	x9, xzr, x9\n"
	// x[2:] * x[1]
		"MUL	x13, x1, x2\n"
		"UMULH	x14, x1, x2\n"
		"ADDS	x7, x13, x7\n"		//acc3
		"ADCS	x9, x14, x9\n"		//acc4
		"ADC	x10, xzr, xzr\n"	//acc5

		"MUL	x13, x1, x3\n"
		"UMULH	x14, x1, x3\n"
		"ADDS	x9, x13, x9\n"		//acc4
		"ADC	x10, x14, x10\n"	//acc5
	// x[3] * x[2]
		"MUL	x13, x2, x3\n"
		"UMULH	x11, x2, x3\n"
		"ADDS	x10, x13, x10\n"	//acc5
		"ADC	x11, xzr, x11\n"	//acc6

	// *2
		"ADDS	x5, x5, x5\n"	// acc1
		"ADCS	x6, x6, x6\n"	//acc2
		"ADCS	x7, x7, x7\n"	//acc3
		"ADCS	x9, x9, x9\n"	//acc4
		"ADCS	x10, x10, x10\n"	//acc5
		"ADCS	x11, x11, x11\n"	//acc6
		"ADC	x12, xzr, xzr\n"	//acc7
	// Missing products
		"MUL	x4, x0, x0\n"	//acc0
		"UMULH	x13, x0, x0\n"
		"ADDS	x5, x13, x5\n"	//acc1

		"MUL	x13, x1, x1\n"
		"UMULH	x14, x1, x1\n"
		"ADCS	x6, x13, x6\n"	//acc2
		"ADCS	x7, x14, x7\n"	//acc3

		"MUL	x13, x2, x2\n"
		"UMULH	x14, x2, x2\n"
		"ADCS	x9, x13, x9\n"		//acc4
		"ADCS	x10, x14, x10\n"	//acc5

		"MUL	x13, x3, x3\n"
		"UMULH	x14, x3, x3\n"
		"ADCS	x11, x13, x11\n"	//acc6
		"ADC	x12, x14, x12\n"	//acc7
		: "=r" (acc0), "=r" (acc1), "=r" (acc2), "=r" (acc3), "=r" (acc4),
		"=r" (acc5), "=r" (acc6), "=r" (acc7)
		: "r" (x0), "r" (x1), "r" (x2), "r" (x3)
		: "x13", "x14", "cc");
	asm volatile(
	// First reduction step
		"LSL	x13, x4, 32\n"	// LSL $32, acc0, t0
		"LSR	x14, x4, 32\n"	// LSR $32, acc0, t1
		"ADDS	x5, x5, x4\n"	// ADDS acc0, acc1
		"ADCS	x6, x6, xzr\n"
		"ADCS	x7, x7, xzr\n"
		"ADC	x4, x4, xzr\n"
		"SUBS	x5, x5, x13\n"	// SUBS t0, acc1
		"SBCS	x6, x6, x14\n"
		"SBCS	x7, x7, x13\n"
		"SBC	x4, x4, x14\n"
	// Second reduction step
		"LSL	x13, x5, 32\n"	// LSL $32, acc1, t0
		"LSR	x14, x5, 32\n"	// LSR $32, acc1, t1
		"ADDS	x6, x6, x5\n"	// ADDS acc1, acc2
		"ADCS	x7, x7, xzr\n"
		"ADCS	x4, x4, xzr\n"
		"ADC	x5, x5, xzr\n"
		"SUBS	x6, x6, x13\n"	// SUBS t0, acc2
		"SBCS	x7, x7, x14\n"
		"SBCS	x4, x4, x13\n"	// SBCS t0, acc0
		"SBC	x5, x5, x14\n"
	// Third reduction step
		"LSL	x13, x6, 32\n"	// LSL $32, acc2, t0
		"LSR	x14, x6, 32\n"	// LSR $32, acc2, t1
		"ADDS	x7, x7, x6\n"	// ADDS acc2, acc3
		"ADCS	x4, x4, xzr\n"
		"ADCS	x5, x5, xzr\n"
		"ADC	x6, x6, xzr\n"
		"SUBS	x7, x7, x13\n"	// SUBS t0, acc3
		"SBCS	x4, x4, x14\n"
		"SBCS	x5, x5, x13\n"	// SBCS t0, acc1
		"SBC	x6, x6, x14\n"
	// Last reduction step
		"LSL	x13, x7, 32\n"	// LSL $32, acc3, t0
		"LSR	x14, x7, 32\n"	// LSR $32, acc3, t1
		"ADDS	x4, x4, x7\n"	// ADDS acc3, acc0
		"ADCS	x5, x5, xzr\n"
		"ADCS	x6, x6, xzr\n"
		"ADC	x7, x7, xzr\n"
		"SUBS	x4, x4, x13\n"	// SUBS t0, acc0
		"SBCS	x5, x5, x14\n"
		"SBCS	x6, x6, x13\n"	// SBCS t0, acc2
		"SBC	x7, x7, x14\n"
	// Add bits [511:256] of the sqr result
		"ADDS	x4, x9, x4\n"
		"ADCS	x5, x10, x5\n"
		"ADCS	x6, x11, x6\n"
		"ADCS	x7, x12, x7\n"
		"ADC	x0, xzr, xzr\n"	// t0 = carry
		: "+r" (x0), "+r" (acc0), "+r" (acc1), "+r" (acc2), "+r" (acc3)
		: "r" (acc4), "r" (acc5), "r" (acc6), "r" (acc7)
		: "x13", "x14", "cc", "memory");
	sm2p_reductionS2(result, acc0, acc1, acc2, acc3, x0);
#endif
}

#ifdef	__cplusplus
}
#endif

#endif	//	__MONT_ASM_H__
