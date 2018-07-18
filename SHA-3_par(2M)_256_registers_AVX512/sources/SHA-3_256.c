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

we have in the registers y[0] to y[12] the follows words
      y[0] = 1  0    y[7]  = 16 15
      y[1] = 3  2    y[8]  = 18 17
      y[2] = 6  5    y[9]  = 19 14
      y[3] = 8  7    y[10] = 21 20
      y[4] = 9  4    y[11] = 23 22
      y[5] = 11 10   y[12] = 24 24
      y[6] = 13 12

*/

#define set_zero_256()		_mm256_setzero_si256()
#define set_zero_HI(X)		_mm_move_epi64(X)
#define alignr(X,Y,D)	 	_mm256_alignr_epi8(X,Y,D)

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

#define zero_state		y[0] = set_zero_256();\
				y[1] = set_zero_256();\
				y[2] = set_zero_256();\
				y[3] = set_zero_256();\
				y[4] = set_zero_256();\
				y[5] = set_zero_256();\
				y[6] = set_zero_256();\
				y[7] = set_zero_256();\
				y[8] = set_zero_256();\
				y[9] = set_zero_256();\
				y[10]= set_zero_256();\
				y[11]= set_zero_256();\
				y[12]= set_zero_256();
								
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
				y[0] = xor_256(y[0],x6);\
				y[1] = xor_256(y[1],x0);\
				y[5] = xor_256(y[5],x1);\
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
				y[3] = xor_256(y[3],x4);\
				\
				x3 = load_256(M1,4);\
				x4 = load_256(M2,4);\
				\
				x5 = set_zero_256();\
				y[4] = xor_256(y[4],x2);\
				y[2] = xor_256(y[2],x7);\
				\
				x3 = blend(x5,x3,0x03);\
				x4 = blend(x5,x4,0x03);\
				x1 = blend(x1,x4,0x0F);\
				x0 = blend(x0,x3,0x0F);\
				\
				y[6] = xor_256(y[6],x6);\
				\
				x3 = permut_64(x0,0x53);\
				x4 = permut_64(x1,0x35);\
				x5 = permut_64(x0,0x56);\
				x6 = permut_64(x1,0x65);\
				\
				x3 = xor_256(x3,x4);/*b16 b15 a16 a15*/\
				x6 = xor_256(x6,x5);/* 0  b14  0  a14*/\
				\
				y[7] = xor_256(y[7],x3);\
				y[9] = xor_256(y[9],x6);\

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
				y[0] = xor_256(y[0],x10);/*b1 b0 a1 a0*/\
				y[1] = xor_256(y[1],x6);/*b3 b2 a3 a2*/\
				y[2] = xor_256(y[2],x2);/*b6 b5 a6 a5*/\
				y[3] = xor_256(y[3],x3);/*b8 b7 a8 a7*/\
				y[4] = xor_256(y[4],x8);/*b9 b4 a9 a4*/\
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
				y[5] = xor_256(y[5],x4);\
				y[6] = xor_256(y[6],x7);\
				y[9] = xor_256(y[9],x10);\
				y[7] = xor_256(y[7],x5);\
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
				y[0] = xor_256(y[0],x0);\
				y[1] = xor_256(y[1],x1);\
				y[4] = xor_256(y[4],x6);\
				y[2] = xor_256(y[2],x7);\
				y[3] = xor_256(y[3],x4);\
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
				y[5] = xor_256(y[5],x7);\
				y[6] = xor_256(y[6],x8);\
				y[7] = xor_256(y[7],x2);\
				y[9] = xor_256(y[9],x0);\
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
				y[0] = xor(y[0],x1);\
				y[1] = xor(y[1],x2);\
				y[2] = xor(y[2],x5);\
				y[3] = xor(y[3],x10);\
				y[4] = xor(y[4],x3);\
				y[5] = xor(y[5],x6);\
				y[6] = xor(y[6],x7);\
				y[7] = xor(y[7],x8);\
				y[9] = xor(y[9],x4);\
			

 				
 
#define keccakF(y,rnd)      \
                            /*After the first part of the step theta we'll have the follows values in the follows registers\\
                                    x0 = c0 -- \
                                    x3 = c2 c1 \
                                    x2 = c4 c3 \
                             */\
                            \
                            for(q=0;q<rnd;q+=2){\
                                    x2   = xor(y[4],y[9]);\
                                    x0 = _mm256_ternarylogic_epi64(y[0],y[7],y[5],0x96);\
                                    x0 = _mm256_ternarylogic_epi64(y[10],y[2],x0,0x96);\
                                    \
                                    x1 = _mm256_ternarylogic_epi64(y[1],y[3],y[6],0x96);\
                                    x1 = _mm256_ternarylogic_epi64(y[11],y[8],x1,0x96);\
                                    \
                                    x3 = shuffle_32(x2,0x4E);\
                                    x2 = _mm256_ternarylogic_epi64(x3,x2,y[12],0x96);\
                            \
                                    /*      The secund part of the step theta\
                                            After this part we'll have the follows values in the follows registers\
                                            x0 = d0 d1\\
                                            x1 = d2 d3\\
                                            x2 = d4 d4\\
                                    */\
	                            x4 = alignr(x0,x2,8);\
                                    x3 = xor(x0,ROT(x1,1));\
                                    x0 = xor(x2,ROT(x0,1));\
                                    x2 = xor(x1,ROT(x4,1));\
	                            \
                                    x0   = alignr(x3,x0,8);\
                                    x1   = alignr(x2,x3,8);\
                                    y[2] = xor(y[2],x0);\
                                    x2   = shuffle_32(x2,0xEE);\
                            \
                                    y[2] = ROTV(y[2],r2);/*1, 16*/\
                                    y[4] = xor(y[4],x2);\
                                    y[0] = xor(y[0],x0);\
                                    y[1] = xor(y[1],x1);\
                            \
                                    y[4]  = ROTV(y[4],r4);/*6, 15*/\
                                    y[12] = xor(y[12],x2);\
                                    y[6]  = xor(y[6],x1);\
                                    y[5]  = xor(y[5],x0);\
                                                            \
                                    y[6]  = ROTV(y[6],r6); /*12, 2*/\
                                    y[9]  = xor(y[9],x2);\
                                    y[3]  = xor(y[3],x1);\
                                    y[11] = xor(y[11],x1);\
                            \
                                    y[5] = ROTV(y[5],r5);/*17, 7*/\
                                    y[7] = xor(y[7],x0);\
                                    y[8] = xor(y[8],x1);\
                            \
                                    x2 = alignr(y[6],y[5],8);\
                                    y[8] = ROTV(y[8],r8);/*3, 18*/\
                                    y[3] = ROTV(y[3],r3);/*21, 11*/\
                            \
                                    x4    = blend_32(y[0],y[12],0x33); /*1, 24*/\
                                    y[10] = xor(y[10],x0);\
                                    y[11] = ROTV(y[11],r11);/*19, 9*/\
                            \
                                    x5   = unpackLO_64(y[4],y[0]);\
                                    y[4] = unpackHI_64(y[3],y[4]);\
                                    x0   = ROTV(x4,r0);/*10, 4*/\
                                                            \
                                    x3 = _mm256_ternarylogic_epi64(x5,y[2],x2,0xD2);\
                                    y[10] = ROTV(y[10],r10); /*24, 14*/\
                                    x1    = alignr(x0,y[11],8);\
                            \
                                    y[9] = ROTV(y[9],r9);/*13, 22*/\
                                    x4 = _mm256_ternarylogic_epi64(y[2],x2,y[8],0xD2);\
                                    x2 = _mm256_ternarylogic_epi64(x2,y[8],x1,0xD2);\
                            \
                                    y[5] = unpackLO_64(y[9],y[5]);\
                                    y[8] = _mm256_ternarylogic_epi64(y[8],x1,x5,0xD2);\
                                    x1 = _mm256_ternarylogic_epi64(x1,x5,y[2],0xD2);\
                            \
                                    y[1]  = ROTV(y[1],r1); /*5, 20*/\
                                    y[11] = alignr(y[11],y[10],8);\
                            \
                                    y[0] = _mm256_ternarylogic_epi64(y[1],y[4],y[5],0xD2);\
                                    y[7] = ROTV(y[7],r7); /*8, 23*/\
                                    y[3] = blend_32(y[3],y[6],0xCC);\
                                    y[6] = unpackHI_64(y[6],y[9]);\
                                                            \
                                    y[2] = _mm256_ternarylogic_epi64(y[4],y[5],y[7],0xD2);\
                                    y[4] = _mm256_ternarylogic_epi64(y[11],y[1],y[4],0xD2);\
                                    y[9] = _mm256_ternarylogic_epi64(y[5],y[7],y[11],0xD2);\
                                    y[7] = _mm256_ternarylogic_epi64(y[7],y[11],y[1],0xD2);\
                            \
                                    /*I'm not sure if this is the best way to compute these values*/\
                                    x5    = alignr(y[3],x0,8);\
                                    x0    = blend_32(x0,y[10],0x33);\
                                    y[10] = alignr(y[10],y[6],8);\
                            \
                                    y[1]  = unpackHI_64(x2,y[8]);\
                                    y[11] = unpackLO_64(y[9],y[7]);\
                                    y[8]  = unpackLO_64(x2,y[8]);\
                            \
                                    x2   = _mm256_ternarylogic_epi64(y[10],x0,x5,0xD2);\
                                    y[5] = _mm256_ternarylogic_epi64(x5,y[3],y[6],0xD2);\
                                    y[6] = _mm256_ternarylogic_epi64(y[6],x2,x0,0xD2);\
                            \
                                    y[3] = unpackHI_64(y[9],y[7]);\
                                    y[9] = alignr(x1,x2,8);\
                            \
                                   /*x2    = _mm_slli_si128(y[2],8);\
                                      x5    = _mm_slli_si128(x4,8);\
                                      y[10] = blend_32(y[0],x2,0xCC);\
                                      y[7]  = blend_32(x3,x5,0xCC);\
                                    */\
                                    \
                                    y[10] = unpackLO_64(y[0],y[2]);\
                                    y[7] = unpackLO_64(x3,x4);\
                            \
                                    y[12] = shuffle_32(y[4],0x44);\
                                    y[4] = unpackHI_64(x1,y[4]);\
                                    y[2] = unpackHI_64(y[0],y[2]);\
                                    y[0] = unpackHI_64(x3,x4);\
                                    x0 = load(RC,q);\
                                    y[0] = xor(y[0],x0);\
				                                     \
                                    x2   = xor(y[4],y[9]);\
                                    x0 = _mm256_ternarylogic_epi64(y[0],y[7],y[5],0x96);\
                                    x0 = _mm256_ternarylogic_epi64(y[10],y[2],x0,0x96);\
                                    \
                                    x1 = _mm256_ternarylogic_epi64(y[1],y[3],y[6],0x96);\
                                    x1 = _mm256_ternarylogic_epi64(y[11],y[8],x1,0x96);\
                                    \
                                    x3 = shuffle_32(x2,0x4E);\
                                    x2 = _mm256_ternarylogic_epi64(x3,x2,y[12],0x96);\
                            \
                                    /*      The secund part of the step theta\
                                            After this part we'll have the follows values in the follows registers\
                                            x0 = d0 d1\\
                                            x1 = d2 d3\\
                                            x2 = d4 d4\\
                                    */\
                                    x4 = alignr(x0,x2,8);\
                                    x3 = xor(x0,ROT(x1,1));\
                                    x0 = xor(x2,ROT(x0,1));\
                                    x2 = xor(x1,ROT(x4,1));\
                            \
                                   x0   = alignr(x3,x0,8);\
                                    x1   = alignr(x2,x3,8);\
                                    y[2] = xor(y[2],x0);\
                                    x2   = shuffle_32(x2,0xEE);\
                            \
                                    y[2] = ROTV(y[2],r2);/*1, 16*/\
                                    y[4] = xor(y[4],x2);\
                                    y[0] = xor(y[0],x0);\
                                    y[1] = xor(y[1],x1);\
                            \
                                    y[4]  = ROTV(y[4],r4);/*6, 15*/\
                                    y[12] = xor(y[12],x2);\
                                    y[6]  = xor(y[6],x1);\
                                    y[5]  = xor(y[5],x0);\
                                                            \
                                    y[6]  = ROTV(y[6],r6); /*12, 2*/\
                                    y[9]  = xor(y[9],x2);\
                                    y[3]  = xor(y[3],x1);\
                                    y[11] = xor(y[11],x1);\
                            \
                                    y[5] = ROTV(y[5],r5);/*17, 7*/\
                                    y[7] = xor(y[7],x0);\
                                    y[8] = xor(y[8],x1);\
                            \
                                    x2 = alignr(y[6],y[5],8);\
                                    y[8] = ROTV(y[8],r8);/*3, 18*/\
                                    y[3] = ROTV(y[3],r3);/*21, 11*/\
                            \
                                    x4    = blend_32(y[0],y[12],0x33); /*1, 24*/\
                                    y[10] = xor(y[10],x0);\
                                    y[11] = ROTV(y[11],r11);/*19, 9*/\
                            \
                                    x5   = unpackLO_64(y[4],y[0]);\
                                    y[4] = unpackHI_64(y[3],y[4]);\
                                    x0   = ROTV(x4,r0);/*10, 4*/\
                                                            \
                                    x3 = _mm256_ternarylogic_epi64(x5,y[2],x2,0xD2);\
                                    y[10] = ROTV(y[10],r10); /*24, 14*/\
                                    x1    = alignr(x0,y[11],8);\
                            \
                                    y[9] = ROTV(y[9],r9);/*13, 22*/\
                                    x4 = _mm256_ternarylogic_epi64(y[2],x2,y[8],0xD2);\
                                    x2 = _mm256_ternarylogic_epi64(x2,y[8],x1,0xD2);\
                            \
                                    y[5] = unpackLO_64(y[9],y[5]);\
                                    y[8] = _mm256_ternarylogic_epi64(y[8],x1,x5,0xD2);\
                                    x1 = _mm256_ternarylogic_epi64(x1,x5,y[2],0xD2);\
                            \
                                    y[1]  = ROTV(y[1],r1); /*5, 20*/\
                                    y[11] = alignr(y[11],y[10],8);\
                            \
                                    y[0] = _mm256_ternarylogic_epi64(y[1],y[4],y[5],0xD2);\
                                    y[7] = ROTV(y[7],r7); /*8, 23*/\
                                    y[3] = blend_32(y[3],y[6],0xCC);\
                                    y[6] = unpackHI_64(y[6],y[9]);\
                                                            \
                                    y[2] = _mm256_ternarylogic_epi64(y[4],y[5],y[7],0xD2);\
                                    y[4] = _mm256_ternarylogic_epi64(y[11],y[1],y[4],0xD2);\
                                    y[9] = _mm256_ternarylogic_epi64(y[5],y[7],y[11],0xD2);\
                                    y[7] = _mm256_ternarylogic_epi64(y[7],y[11],y[1],0xD2);\
                            \
                                    /*I'm not sure if this is the best way to compute these values*/\
                                    x5    = alignr(y[3],x0,8);\
                                    x0    = blend_32(x0,y[10],0x33);\
                                    y[10] = alignr(y[10],y[6],8);\
                            \
                                    y[1]  = unpackHI_64(x2,y[8]);\
                                    y[11] = unpackLO_64(y[9],y[7]);\
                                    y[8]  = unpackLO_64(x2,y[8]);\
                            \
                                    x2   = _mm256_ternarylogic_epi64(y[10],x0,x5,0xD2);\
                                    y[5] = _mm256_ternarylogic_epi64(x5,y[3],y[6],0xD2);\
                                    y[6] = _mm256_ternarylogic_epi64(y[6],x2,x0,0xD2);\
                                    y[3] = unpackHI_64(y[9],y[7]);\
                                    y[9] = alignr(x1,x2,8);\
                            \
                                    /*x2    = _mm_slli_si128(y[2],8);\
                                      x5    = _mm_slli_si128(x4,8);\
                                      y[10] = blend_32(y[0],x2,0xCC);\
                                      y[7]  = blend_32(x3,x5,0xCC);\
                                    */\
                                    \
                                    y[10] = unpackLO_64(y[0],y[2]);\
                                    y[7] = unpackLO_64(x3,x4);\
                            \
                                    y[12] = shuffle_32(y[4],0x44);\
                                    y[4] = unpackHI_64(x1,y[4]);\
                                    y[2] = unpackHI_64(y[0],y[2]);\
                                    y[0] = unpackHI_64(x3,x4);\
                                    x0 = load(RC,q+1);\
                                    y[0] = xor(y[0],x0);\
                            }\

void print_256(__m256i x){
     uint64_t c[4];
     store(c,0,x);
     printf("\n 2 = %8.16lx\t  2 = %8.16lx\t 1 = %8.16lx\t  0 = %8.16lx\t\n", c[3],c[2],c[1],c[0]);
}							
			
int keccak(char **in_e, int inlen, uint8_t **md, int rsiz){

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


  __m256i y[13];
  __m256i x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10;
   
  char *in_temp1 = in_e[0];
  char *in_temp2 = in_e[1];
  int l=0,q,k=1,it,i,o;
  ALIGN uint8_t temp1[144];
  ALIGN uint8_t temp2[144];
  
  zero_state;
  
  memset(temp1, 0, 144*sizeof(uint8_t));
  memset(temp2, 0, 144*sizeof(uint8_t));
  
  while(inlen-l >= rsiz){
    
    load_one(in_temp1,in_temp2);
     
    keccakF(y,24);
 
    inlen -= 128;
    l=8;
    in_temp1 += 128;
    in_temp2 += 128;
    
    if(inlen-l < rsiz){
	  k=2;
	  break;
    }
    
    load_two(in_temp1,in_temp2);
   
    keccakF(y,24);

    l=16;
    inlen -= 128;
    in_temp1 += 128;
    in_temp2 += 128;
    
    if(inlen-l < rsiz){
	  k=3;
	  break;
    }
    
    load_three(in_temp1,in_temp2);  
    keccakF(y,24);
    l=24;
    inlen -= 128;
    in_temp1 += 128;
    in_temp2 += 128;
    

    
    if(inlen-l < rsiz){
	  k=4;
	  break;
    }    
    
    load_four(in_temp1,in_temp2);
    keccakF(y,24);
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
   keccakF(y,24);
/*   for(o=0;o<13;o++)\
	print_256(y[o]);\
    getchar();\*/

   x0 = permut_2x128_avx(y[0],y[1],0x20);
   x1 = permut_2x128_avx(y[0],y[1],0x31);
   
  store(md[0],0,x0);
  store(md[1],0,x1);
    
  return 0;
    
}

