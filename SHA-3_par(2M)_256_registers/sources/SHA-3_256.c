#include "keccak.h"

/*This file was implemented by Roberto Cabral; email:rbcabral@ic.unicamp.br 
 
In this file, we have the implementation of the SHA3-256 algorithm according to draft FIPS PUB 202.

This implementation is optimized to 128-bit registers and uses AVX-2 instructions. 
We have two different methods to read the data from the input. This is needed because, for this version,
we processes 17 word of 64 bits per time in the absorbing phase and , as 17 * 64 is not
multiple of 128, the data will not be aligned to do a efficient load in the next round. Thus, 
we use the other way to compute the load efficiently.

We have a define "round" that is a round of the permutations function. For this version,
the permutation function has 24 rounds.

The constants r0 to r11 are the constants used in to do the rotations in the rho step.

The constants RC are the constants used in to do the iota step.

To take advantage of the large registers we repesente the state as the follows

we have in the registers y0 to y12 the follows words
      y0 = 1  0    y7  = 16 15
      y1 = 3  2    y8  = 18 17
      y2 = 6  5    y9  = 19 14
      y3 = 8  7    y10 = 21 20
      y4 = 9  4    y11 = 23 22
      y5 = 11 10   y12 = 24 24
      y6 = 13 12

*/

#define and(X,Y)		_mm_and_si128(X,Y)

#define set_zero_256()		_mm256_setzero_si256()
#define set_zero_HI(X)		_mm_move_epi64(X)
#define alignr(X,Y,D)	 	_mm256_alignr_epi8(X,Y,D)
#define shiftR(Y,D)		_mm_srli_si128(Y,D)
#define shiftL(Y,D)		_mm_slli_si128(Y,D)


#define store(X,D,Y)		_mm256_store_si256((__m256i*)X+D,Y)
#define load(M,D) 		_mm256_load_si256((__m256i*)M+D)
#define unpackLO_64(X,Y)	_mm256_unpacklo_epi64(X,Y)
#define unpackHI_64(X,Y)	_mm256_unpackhi_epi64(X,Y)
#define xor(X,Y)		_mm256_xor_si256(X,Y)
#define and_not(X,Y)		_mm256_andnot_si256(X,Y)
#define load_256(M,d) 	  	_mm256_load_si256((__m256i*)M+d)
#define permut_2x128_avx(y,z,d) _mm256_permute2f128_si256(y,z,d)
#define xor_256(y,z) 	  	_mm256_xor_si256(y,z)
#define blend(y,z,d) 		_mm256_blend_epi32(y,z,d)
#define blend_32(y,z,d) 		_mm256_blend_epi32(y,z,d)
#define permut_64(y,d)  	_mm256_permute4x64_epi64(y,d)
#define shuffle_32(X,D)		_mm256_shuffle_epi32(X,D)
#define ROT(Y,D) 		_mm256_xor_si256(_mm256_add_epi64(Y,Y),_mm256_srli_epi64(Y,(63)))
#define ROTV(x,d) 		_mm256_xor_si256(_mm256_sllv_epi64(x,d),_mm256_srlv_epi64(x,_mm256_sub_epi64(_mm256_set1_epi64x(0x40),d)))

#ifdef __INTEL_COMPILER
#define BARRIER __memory_barrier()
#else
#define BARRIER asm volatile("" ::: "memory")
#endif

#define zero_state		y0 = set_zero_256();\
				y1 = set_zero_256();\
				y2 = set_zero_256();\
				y3 = set_zero_256();\
				y4 = set_zero_256();\
				y5 = set_zero_256();\
				y6 = set_zero_256();\
				y7 = set_zero_256();\
				y8 = set_zero_256();\
				y9 = set_zero_256();\
				y10= set_zero_256();\
				y11= set_zero_256();\
				y12= set_zero_256();
								
#define load_one(M1,M2)		\
				x0 = load_256(M1,0);\
				x1 = load_256(M2,0);\
				x4 = load_256(M1,2);\
				x5 = load_256(M2,2);\
				\
				x6 = permut_2x128_avx(x0,x1,0x20);/*b1 b0 a1 a0*/\
				x0 = permut_2x128_avx(x0,x1,0x31);/*b3 b2 a3 a2*/\
				x1 = permut_2x128_avx(x4,x5,0x31);/*b11 b10 a11 a10*/\
				x2 = load_256(M1,1);\
				x3 = load_256(M2,1);\
				\
				y0 = xor_256(y0,x6);\
				y1 = xor_256(y1,x0);\
				y5 = xor_256(y5,x1);\
				\
				x7 = blend(x2,x4,0x0C);\
				x9 = blend(x2,x4,0x03);\
				x8 = blend(x3,x5,0x0C);\
				x10= blend(x3,x5,0x03);\
				x0 = load_256(M1,3);\
				x1 = load_256(M2,3);\
				\
				x3 = permut_64(x9,0x39);\
				x4 = permut_64(x10,0x39);\
				x2 = permut_2x128_avx(x7,x8,0x20);/*b9 b4 a9 a4*/\
				x7 = permut_2x128_avx(x3,x4,0x20);/*b6 b5 a6 a5*/\
				x4 = permut_2x128_avx(x3,x4,0x31);/*b8 b7 a8 a7*/\
				x6 = permut_2x128_avx(x0,x1,0x20);/*b13 b12 a13 a12*/\
				y3 = xor_256(y3,x4);\
				\
				x3 = load_256(M1,4);\
				x4 = load_256(M2,4);\
				\
				x5 = set_zero_256();\
				y4 = xor_256(y4,x2);\
				y2 = xor_256(y2,x7);\
				\
				x3 = blend(x5,x3,0x03);\
				x4 = blend(x5,x4,0x03);\
				x1 = blend(x1,x4,0x0F);\
				x0 = blend(x0,x3,0x0F);\
				\
				y6 = xor_256(y6,x6);\
				\
				x3 = permut_64(x0,0x53);\
				x4 = permut_64(x1,0x35);\
				x5 = permut_64(x0,0x56);\
				x6 = permut_64(x1,0x65);\
				\
				x3 = xor_256(x3,x4);/*b16 b15 a16 a15*/\
				x6 = xor_256(x6,x5);/* 0  b14  0  a14*/\
				\
				y7 = xor_256(y7,x3);\
				y9 = xor_256(y9,x6);\

#define load_two(M1,M2)	\
				x0 = load_256(M1,0);\
				x1 = load_256(M2,0);\
				x2 = load_256(M1,1);\
				x3 = load_256(M2,1);\
				x4 = load_256(M1,2);\
				x5 = load_256(M2,2);\
				\
				x6 = blend_32(x0,x2,0x03);\
				x7 = blend_32(x1,x3,0x03);\
				x8 = blend_32(x2,x4,0xF0);\
				x9 = blend_32(x3,x5,0xF0);\
				\
				x6 = permut_64(x6,0x39);\
				x7 = permut_64(x7,0x39);\
				x8 = permut_64(x8,0x99);\
				x9 = permut_64(x9,0x99);\
				\
				x10= permut_2x128_avx(x6,x7,0x20);/*b1 b0 a1 a0*/\
				x6 = permut_2x128_avx(x6,x7,0x31);/*b3 b2 a3 a2*/\
				\
				x2 = permut_2x128_avx(x2,x3,0x31);/*b6 b5 a6 a5*/\
				x3 = permut_2x128_avx(x4,x5,0x20);/*b8 b7 a8 a7*/\
				/*x8 = permut_2x128_avx(x8,x9,0x30);b9 b4 a9 a4    cab be do with blend*/\
				\
				x8 = blend_32(x8,x9,0xF0);\
				y0 = xor_256(y0,x10);/*b1 b0 a1 a0*/\
				y1 = xor_256(y1,x6);/*b3 b2 a3 a2*/\
				y2 = xor_256(y2,x2);/*b6 b5 a6 a5*/\
				y3 = xor_256(y3,x3);/*b8 b7 a8 a7*/\
				y4 = xor_256(y4,x8);/*b9 b4 a9 a4*/\
				\
				x0 = load_256(M1,3);\
				x1 = load_256(M2,3);\
				x2 = load_256(M1,4);\
				x3 = load_256(M2,4);\
				\
				x8 = set_zero_256();\
				x6 = blend_32(x0,x4,0xC0);\
				x7 = blend_32(x1,x5,0xC0);\
				x9 = blend_32(x0,x8,0x3F);\
				x10= blend_32(x1,x8,0x3F);\
				\
				x6 = permut_64(x6,0x93);\
				x7 = permut_64(x7,0x93);\
				x9 = permut_64(x9,0x33);\
				x10= permut_64(x10,0x33);\
				\
				x4 = permut_2x128_avx(x6,x7,0x20);/*b11 b10 a11 a10*/\
				x7 = permut_2x128_avx(x6,x7,0x31);/*b13 b12 a13 a12*/\
				\
				x5 = permut_2x128_avx(x2,x3,0x20);/*b16 b15 a16 a15*/\
				/*x10= permut_2x128_avx(x9,x10,0x30);0 b14 0 a14*/\
				\
				x10= blend_32(x9,x10,0xF0);\
				y5 = xor_256(y5,x4);\
				y6 = xor_256(y6,x7);\
				y9 = xor_256(y9,x10);\
				y7 = xor_256(y7,x5);\
				\
				
#define load_three(M1,M2)	\
				x0 = load_256(M1,0);\
				x1 = load_256(M2,0);\
				x2 = load_256(M1,1);\
				x3 = load_256(M2,1);\
				x4 = load_256(M1,2);\
				x5 = load_256(M2,2);\
				\
				x0 = permut_2x128_avx(x0,x1,0x31);/*b1 b0 a1 a0*/\
				x1 = permut_2x128_avx(x2,x3,0x20);/*b3 b2 a3 a2*/\
				\
				x6 = blend(x2,x4,0xC0);\
				x7 = blend(x3,x5,0xC0);\
				x2 = blend(x2,x4,0x3F);\
				x3 = blend(x3,x5,0x3F);\
				\
				x2 = permut_64(x2,0x93);\
				x3 = permut_64(x3,0x93);\
				x6 = permut_2x128_avx(x6,x7,0x31);/*b9 b4 a9 a4*/\
				x7 = permut_2x128_avx(x2,x3,0x20);/*b6 b5 a6 a5*/\
				x4 = permut_2x128_avx(x2,x3,0x31);/*b8 b7 a8 a7*/\
				\
				y0 = xor_256(y0,x0);\
				y1 = xor_256(y1,x1);\
				y4 = xor_256(y4,x6);\
				y2 = xor_256(y2,x7);\
				y3 = xor_256(y3,x4);\
				\
				\
				x4 = set_zero_256();\
				x0 = load_256(M1,3);\
				x1 = load_256(M2,3);\
				x2 = load_256(M1,4);\
				x3 = load_256(M2,4);\
				\
				x5 = blend_32(x2,x4,0xC0);\
				x6 = blend_32(x3,x4,0xC0);\
				\
				x7 = permut_2x128_avx(x0,x1,0x20);\
				x8 = permut_2x128_avx(x0,x1,0x31);\
				x0 = permut_64(x5,0xCC);\
				x5 = permut_64(x5,0x99);\
				x2 = permut_64(x6,0xCC);\
				x6 = permut_64(x6,0x99);\
				x0 = permut_2x128_avx(x0,x2,0x30);\
				x2 = permut_2x128_avx(x5,x6,0x30);\
				\
				y5 = xor_256(y5,x7);\
				y6 = xor_256(y6,x8);\
				y7 = xor_256(y7,x2);\
				y9 = xor_256(y9,x0);\
				\

#define load_four(M1,M2)	\
				x0 = load_256(M1,0);\
				x1 = load_256(M2,0);\
				x2 = load_256(M1,1);\
				x3 = load_256(M2,1);\
				x4 = load_256(M1,3);\
				x5 = load_256(M2,3);\
				x8 = load_256(M1,4);\
				x9 = load_256(M2,4);\
				\
				x10= set_zero_256();\
				x6 = blend_32(x2,x4,0x3F);\
				x7 = blend_32(x3,x5,0x3F);\
				x4 = blend_32(x4,x10,0x3F);\
				x5 = blend_32(x5,x10,0x3F);\
				\
				x6 = permut_64(x6,0x39);\
				x7 = permut_64(x7,0x39);\
				\
				x2 = blend_32(x2,x0,0xC0);\
				x3 = blend_32(x3,x1,0xC0);\
				x4 = blend_32(x4,x8,0x0F);\
				x5 = blend_32(x5,x9,0x0F);\
				\
				x2 = permut_64(x2,0x93);\
				x3 = permut_64(x3,0x93);\
				x4 = permut_64(x4,0x93);\
				x5 = permut_64(x5,0x93);\
				\
				x0 = load_256(M1,2);\
				x10= load_256(M2,2);\
				\
				x8 = permut_2x128_avx(x8,x9,0x31);/*b16 b15 a16 a15*/\
				x1 = permut_2x128_avx(x2,x3,0x20);/*b1 b0 a1 a0*/\
				x2 = permut_2x128_avx(x2,x3,0x31);/*b3 b2 a3 a2*/\
				x3 = permut_2x128_avx(x6,x7,0x31);/*b9 b4 a9 a4*/\
				x6 = permut_2x128_avx(x6,x7,0x20);/*b11 b10 a11 a10*/\
				x7 = permut_2x128_avx(x4,x5,0x20);/*b13 b12 a13 a12*/\
				x4 = permut_2x128_avx(x4,x5,0x31);/*0 b14 0 a14*/\
				x5 = permut_2x128_avx(x0,x10,0x20);/*b6 b5 a6 a5*/\
				x10= permut_2x128_avx(x0,x10,0x31);/*b8 b7 a8 a7*/\
				\
				y0 = xor(y0,x1);\
				y1 = xor(y1,x2);\
				y2 = xor(y2,x5);\
				y3 = xor(y3,x10);\
				y4 = xor(y4,x3);\
				y5 = xor(y5,x6);\
				y6 = xor(y6,x7);\
				y7 = xor(y7,x8);\
				y9 = xor(y9,x4);\
			

 				
 
#define round			\
				/*After the first part of the step theta we'll have the follows values in the follows registers\
				  x0 = c0 ** \
				  x3 = c2 c1 \
				  x2 = c4 c3 \
				 */\
				 \
				x0 = xor(y5,y7);\
				x1 = xor(y1,y3);\
				x2 = xor(y4,y9);\
				\
				x0 = xor(x0,y0);\
				x1 = xor(x1,y6);\
				x3 = shuffle_32(x2,0x4E);\
				\
				x0 = xor(x0,y2);\
				x1 = xor(x1,y8);\
				x2 = xor(x2,x3);\
				\
				x0 = xor(x0,y10);\
				x1 = xor(x1,y11);\
				x2 = xor(x2,y12);\
				/*The secund part of the step theta\
				After this part we'll have the follows values in the follows registers\
				x0 = d0 d1\
				x1 = d2 d3\
				x2 = d4 d4\
				*/\
				\
				x4 = alignr(x0,x2,8);\
				x3 = xor(x0,ROT(x1,1));\
				x0 = xor(x2,ROT(x0,1));\
				x2 = xor(x1,ROT(x4,1));\
				\
				x0 = alignr(x3,x0,8);\
				x1 = alignr(x2,x3,8);\
				y2 = xor(y2,x0);\
				x2 = shuffle_32(x2,0xEE);/*Can be other way to do this. Here I'm doing a broadcast of the high word*/\
				\
				\
				y2 = ROTV(y2,r2);\
				y4 = xor(y4,x2);\
				y0 = xor(y0,x0);\
				y1 = xor(y1,x1);\
				\
				y4 = ROTV(y4,r4);\
				y12 = xor(y12,x2);\
				y6 = xor(y6,x1);\
				y5 = xor(y5,x0);\
				\
				y6 = ROTV(y6,r6);\
				y9 = xor(y9,x2);\
				y3 = xor(y3,x1);\
				y11 = xor(y11,x1);\
				\
				y5 = ROTV(y5,r5);\
				y7 = xor(y7,x0);\
				y8 = xor(y8,x1);\
				\
				x2 = alignr(y6,y5,8);\
				y8 = ROTV(y8,r8);\
				y3 = ROTV(y3,r3);\
				\
				x4 = blend_32(y0,y12,0x33);\
				y10 = xor(y10,x0);\
				y11 = ROTV(y11,r11);\
				\
				x5 = unpackLO_64(y4,y0);\
				y4 = unpackHI_64(y3,y4);\
				x0 = ROTV(x4,r0);\
				\
				x3 = and_not(y2,x2);\
				x3 = xor(x5,x3);\
				y10 = ROTV(y10,r10);\
				x1 = alignr(x0,y11,8);\
				\
				\
				y9 = ROTV(y9,r9);\
				x4 = and_not(x2,y8);\
				x4 = xor(y2,x4);\
				x2 = xor(x2,and_not(y8,x1));\
				\
				y5 = unpackLO_64(y9,y5);\
				y8 = xor(y8,and_not(x1,x5));\
				x1 = xor(x1,and_not(x5,y2));\
				\
				y1 = ROTV(y1,r1);\
				y11 = alignr(y11,y10,8);\
				\
				\
				y0 = xor(y1,and_not(y4,y5));\
				y7 = ROTV(y7,r7);\
				y3 = blend_32(y3,y6,0xCC);\
				y6 = unpackHI_64(y6,y9);\
				\
				y2 = xor(y4,and_not(y5,y7));\
				y4 = xor(y11,and_not(y1,y4));\
				y9 = xor(y5,and_not(y7,y11));\
				y7 = xor(y7,and_not(y11,y1));\
				\
				/*I'm not sure if this is the best way to compute these values*/\
				x5 = alignr(y3,x0,8);\
				x0 = blend_32(x0,y10,0x33);\
				y10 = alignr(y10,y6,8);\
				\
				y1 = unpackHI_64(x2,y8);\
				y11 = unpackLO_64(y9,y7);\
				y8 = unpackLO_64(x2,y8);\
				\
				x2 = xor(y10,and_not(x0,x5));\
				y5 = xor(x5,and_not(y3,y6));\
				y6 = xor(y6,and_not(x2,x0));\
				\
				y3 = unpackHI_64(y9,y7);\
				y9 = alignr(x1,x2,8);\
				\
				/*x2 = _mm256_slli_si256(y2,8);\
				x5 = _mm256_slli_si256(x4,8);\
				y10 = blend_32(y0,x2,0xCC);\
				y7  = blend_32(x3,x5,0xCC);\*/\
				\
				y10 = unpackLO_64(y0,y2);\
				y7 = unpackLO_64(x3,x4);\
				\
				y12 = shuffle_32(y4,0x44);\
				y4 = unpackHI_64(x1,y4);\
				y2 = unpackHI_64(y0,y2);\
				y0 = unpackHI_64(x3,x4);\

#define	other_rounds1  		\
				for(q=0;q<24;q++){\
				  round;\
				  x0 = load(RC,q);\
				  y0 = xor(y0,x0);\
				}\
				
#define other_rounds	        \
				round;\
				\
				x0 = load(RC,0);\
				y0 = xor(y0,x0);\
				\
				round;/*2 roudn*/\
				x0 = load(RC,1);\
				y0 = xor(y0,x0);\
				\
				round;/*3 round*/\
				x0 = load(RC,2);\
				y0 = xor(y0,x0);\
				\
				\
				round;/*4 round*/\
				x0 = load(RC,3);\
				y0 = xor(y0,x0);\
				\
				round;/*5 round*/\
				x0 = load(RC,4);\
				y0 = xor(y0,x0);\
				\
				round;/*6 round*/\
				x0 = load(RC,5);\
				y0 = xor(y0,x0);\
				\
				round;/*7 round*/\
				x0 = load(RC,6);\
				y0 = xor(y0,x0);\
				\
				round/*8 round*/\
				x0 = load(RC,7);\
				y0 = xor(y0,x0);\
				\
				round;/*9 round*/\
				x0 = load(RC,8);\
				y0 = xor(y0,x0);\
				\
				round;/*10 round*/\
				x0 = load(RC,9);\
				y0 = xor(y0,x0);\
				\
				round;/*11 round*/\
				x0 = load(RC,10);\
				y0 = xor(y0,x0);\
				\
				round;/*12 round*/\
				x0 = load(RC,11);\
				y0 = xor(y0,x0);\
				\
				round;/*13 round*/\
				x0 = load(RC,12);\
				y0 = xor(y0,x0);\
				\
				round;/*14 round*/\
				x0 = load(RC,13);\
				y0 = xor(y0,x0);\
				\
				round;/*15 round*/\
				x0 = load(RC,14);\
				y0 = xor(y0,x0);\
				\
				round;/*16 round*/\
				x0 = load(RC,15);\
				y0 = xor(y0,x0);\
				\
				round;/*17 round*/\
				x0 = load(RC,16);\
				y0 = xor(y0,x0);\
				\
				round;/*18 round*/\
				x0 = load(RC,17);\
				y0 = xor(y0,x0);\
				\
				round;/*19 round*/\
				x0 = load(RC,18);\
				y0 = xor(y0,x0);\
				\
				round;/*20 round*/\
				x0 = load(RC,19);\
				y0 = xor(y0,x0);\
				\
				round;/*21 round*/\
				x0 = load(RC,20);\
				y0 = xor(y0,x0);\
				\
				round;/*22 round*/\
				x0 = load(RC,21);\
				y0 = xor(y0,x0);\
				\
				round;/*23 round*/\
				x0 = load(RC,22);\
				y0 = xor(y0,x0);\
				\
				round;/*24 round*/\
				x0 = load(RC,23);\
				y0 = xor(y0,x0);				
				
			
int keccak_256(char **in_e, int inlen, uint8_t **md, int rsiz){

const ALIGN uint64_t RC[96] = {
    0x0000000000000001,0x0, 0x0000000000000001,0x0, 
    0x0000000000008082,0x0, 0x0000000000008082,0x0,
    0x800000000000808A,0x0, 0x800000000000808A,0x0,  
    0x8000000080008000,0x0, 0x8000000080008000,0x0,
    0x000000000000808B,0x0, 0x000000000000808B,0x0,  
    0x0000000080000001,0x0, 0x0000000080000001,0x0,
    0x8000000080008081,0x0, 0x8000000080008081,0x0,
    0x8000000000008009,0x0, 0x8000000000008009,0x0,
    0x000000000000008A,0x0, 0x000000000000008A,0x0,  
    0x0000000000000088,0x0, 0x0000000000000088,0x0,
    0x0000000080008009,0x0, 0x0000000080008009,0x0,  
    0x000000008000000A,0x0, 0x000000008000000A,0x0,
    0x000000008000808B,0x0, 0x000000008000808B,0x0,
    0x800000000000008B,0x0, 0x800000000000008B,0x0,
    0x8000000000008089,0x0, 0x8000000000008089,0x0,
    0x8000000000008003,0x0, 0x8000000000008003,0x0,
    0x8000000000008002,0x0, 0x8000000000008002,0x0,
    0x8000000000000080,0x0, 0x8000000000000080,0x0,
    0x000000000000800A,0x0, 0x000000000000800A,0x0,
    0x800000008000000A,0x0, 0x800000008000000A,0x0,
    0x8000000080008081,0x0, 0x8000000080008081,0x0,  
    0x8000000000008080,0x0, 0x8000000000008080,0x0,
    0x0000000080000001,0x0, 0x0000000080000001,0x0,  
    0x8000000080008008,0x0, 0x8000000080008008,0x0
  };
  
  
  const __m256i r2 = _mm256_set_epi64x(44,36,44,36);
  const __m256i r4 = _mm256_set_epi64x(20,27,20,27);
  const __m256i r6 = _mm256_set_epi64x(25,43,25,43);
  const __m256i r5 = _mm256_set_epi64x(10,3,10,3);
  const __m256i r8 = _mm256_set_epi64x(21,15,21,15);
  const __m256i r3 = _mm256_set_epi64x(55,6,55,6);
  const __m256i r11= _mm256_set_epi64x(56,61,56,61);
  const __m256i r0 = _mm256_set_epi64x(1,14,1,14);
  const __m256i r10 = _mm256_set_epi64x(2,18,2,18);
  const __m256i r9 = _mm256_set_epi64x(8,39,8,39);
  const __m256i r1 = _mm256_set_epi64x(28,62,28,62);
  const __m256i r7 = _mm256_set_epi64x(45,41,45,41);


  __m256i y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12;
  __m256i x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10;
   
  char *in_temp1 = in_e[0];
  char *in_temp2 = in_e[1];
  int l=0,q,k=1,it,i;
  ALIGN uint8_t temp1[144];
  ALIGN uint8_t temp2[144];
  
  zero_state;
  
  memset(temp1, 0, 144*sizeof(uint8_t));
  memset(temp2, 0, 144*sizeof(uint8_t));
  
  while(inlen-l >= rsiz){
    
    load_one(in_temp1,in_temp2);
    other_rounds(q);
     
    inlen -= 128;
    l=8;
    in_temp1 += 128;
    in_temp2 += 128;
    
    if(inlen-l < rsiz){
	  k=2;
	  break;
    }
    
    load_two(in_temp1,in_temp2);
    other_rounds(q);
    l=16;
    inlen -= 128;
    in_temp1 += 128;
    in_temp2 += 128;
    
    if(inlen-l < rsiz){
	  k=3;
	  break;
    }
    
    load_three(in_temp1,in_temp2);  
    other_rounds;
    l=24;
    inlen -= 128;
    in_temp1 += 128;
    in_temp2 += 128;
    

    
    if(inlen-l < rsiz){
	  k=4;
	  break;
    }    
    
    load_four(in_temp1,in_temp2);
    other_rounds;
    inlen -= 160;
    in_temp1 +=160;
    in_temp2 +=160;
    l=0; 

  }

   
    switch(k){
      case 2:
	inlen-=8;
	in_temp1 += 8;
	in_temp2 += 8;
	break;
      case 3:
	inlen-=16;
	in_temp1 += 16;
	in_temp2 += 16;
	break;
      case 4:
	inlen-=24;	
	in_temp1 += 24;
	in_temp2 += 24;
	break;
      case 1:
	break;
      default:
	printf("default\ninlen : %d\n k = %d\n",inlen,k);
	getchar();
	printf("ERRO!!!");
	return 0;
    }
   
   
   memcpy(temp1, in_temp1, inlen);
   memcpy(temp2, in_temp2, inlen);
   
   temp1[inlen] = 0x06;
   temp2[inlen++] = 0x06;
   temp1[rsiz - 1] |= 0x80;
   temp2[rsiz - 1] |= 0x80;
   
   load_one(temp1,temp2);
   other_rounds(q);
  
   x0 = permut_2x128_avx(y0,y1,0x20);
   x1 = permut_2x128_avx(y0,y1,0x31);
   
  store(md[0],0,x0);
  store(md[1],0,x1);
    
  return 0;
    
}
