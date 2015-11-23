#include "keccak.h"

/*This file was implemented by Roberto Cabral; email:rbcabral@ic.unicamp.br 
 
In this file, we have the implementation of the SHA3-256 algorithm according to draft FIPS PUB 202.

This implementation is optimized to 256-bit registers and uses AVX-2 instructions. 
We have four different methods to read the data from the input. This is needed because, for this version,
we processes 17 word of 64 bits per time in the absorbing phase and , as 17 * 64 is not
multiple of 256, the data will not be aligned to do a efficient load in the next round. Thus, 
we use the other way to compute the load efficiently.

We have a define "round" that is a round of the permutations function. For this version,
the permutation function has 24 rounds.

The constants r0 to r5 are the constants used in to do the rotations in the rho step.

The constants RC are the constants used in to do the iota step.

To take advantage of the large registers we repesente the state as the follows

we have in the registers y0 to y6 the follows words
      y0 =  3  2  1  0
      y1 =  8  7  6  5
      y2 = 13 12 11 10
      y3 = 18 17 16 15
      y4 = 23 22 21 20
      y5 = 19 14  9  4
      y6 = 24 24 24 24

*/


#define load(x,m,d) 	  	x = _mm256_load_si256((__m256i*)m+d)
#define store(x,d,y)		_mm256_store_si256((__m256i*)x+d,y);
#define blend_32(x,y,z,d) 	x = _mm256_blend_epi32(y,z,d)
#define permut_64(x,y,d)  	x = _mm256_permute4x64_epi64(y,d)
#define xor(x,y,z) 	  	x = _mm256_xor_si256(y,z)
#define set_zero_256(x)		x = _mm256_setzero_si256()

#define broadcast_128L_256(x,y) x = _mm256_broadcastq_epi64(y)
#define xor_128(x,y,z)		x = _mm_xor_si128(y,z)
#define load_128(x,m,d) 	x = _mm_load_si128((__m128i*)m+d)
#define store_128(x,d,y)	_mm_store_si128((__m128i*)x+d,y);
#define unpackLO_128_64(x,y,z)	x = _mm256_unpacklo_epi64(y,z)
#define unpackHI_128_64(x,y,z)	x = _mm256_unpackhi_epi64(y,z)
#define cast_to_256(x,y)	x = _mm256_castsi128_si256(y)
#define cast_to_128(x,y)	x = _mm256_castsi256_si128(y)
#define extract_128(x,y,d)	x = _mm256_extracti128_si256(y,d)
#define insert_128(x,y,z,d)	x = _mm256_inserti128_si256(y,z,d)
#define blend_16_128(x,y,z,d)	x = _mm_blend_epi16(y,z,d)
#define xor_and_not(x,y,z)	x = _mm256_xor_si256(x,(_mm256_andnot_si256(y,z)));
#define permut_2x128_avx(x,y,z,d) 	x = _mm256_permute2f128_si256(y,z,d)
#define permut_2x128_avx2(x,y,z,d)	x = _mm256_permute2x128_si256(y,z,d)
#define alignr(y,z,d) 	_mm256_alignr_epi8(y,z,d)

#define RSIZ 136

#define ALIGN __attribute__ ((aligned (32)))
#define ROT(x,y)  _mm256_xor_si256(_mm256_add_epi64(x,x),_mm256_srli_epi64(x,(64-y)))
#define ROTV(x,d) _mm256_xor_si256(_mm256_sllv_epi64(x,d),_mm256_srlv_epi64(x,_mm256_sub_epi64(_mm256_set1_epi64x(0x40),d)))




#define zero_state		set_zero_256(y0);\
				set_zero_256(y1);\
				set_zero_256(y2);\
				set_zero_256(y3);\
				set_zero_256(y4);\
				set_zero_256(y5);\
				set_zero_256(y6);
				
#define load_last	  	\
				load(x1,temp,1);\
				load(x2,temp,2);\
				load(x3,temp,3);\
				load(x6,temp,4);\
				load(x0,temp,0);\
				\
				blend_32(a1,x1,x2,0x0C);\
				blend_32(a0,x1,x2,0x03);\
				blend_32(x6,x6,x3,0xF0);\
				blend_32(x4,a1,x3,0x30);\
				permut_2x128_avx(x2,x2,x3,0x21);\
				set_zero_256(x3);\
				blend_32(x4,x4,x3,0xC0);\
				permut_64(x6,x6,0x03);\
				xor(y6,y6,x4);\
				xor(y0,y0,x0);\
				\
				permut_64(x1,a0,0x39);\
				blend_32(x3,x6,x3,0xF0);\
				xor(y1,y1,x1);\
				extract_128(m0,y6,0);\
				extract_128(m1,y6,1);\
				xor_128(m0,m0,m1);\
				cast_to_256(x1,m0);\
				permut_64(x1,x1,0x55);\
				xor(c1,y0,y4);\
				xor(y2,y2,x2);\
				xor(c1,c1,y1);\
				\
				broadcast_128L_256(x0,m0);\
				xor(x0,x0,y5);\
				xor(x0,x0,x1);\
				xor(c1,c1,y2);\
				xor(y3,y3,x3);\
				xor(c1,c1,y3);\
				\
				xor(x1,x0,ROT(c1,0x01));\
				blend_32(x2,c1,x0,0x0C);\
				permut_64(x1,x1,0x55);\
				permut_64(x2,x2,0x1E);\
				xor(x2,c1,ROT(x2,0x01));\
				permut_64(x3,x2,0xFF);\
				\
				xor(y5,y5,x3);\
				xor(y6,y6,x3);\
				x5 = ROTV(y6,r5);\
				blend_32(x1,x2,x1,0xC0);\
				permut_64(x1,x1,0x93);\
				\
				xor(y0,y0,x1);\
				blend_32(x0,y0,y5,0x03);\
				x4 = ROTV(x0,r0);\
				xor(y1,y1,x1);\
				x6 = ROTV(y1,r1);\
				xor(y2,y2,x1);\
				x2 = ROTV(y2,r2);\
				xor(y3,y3,x1);\
				x3 = ROTV(y3,r3);\
				xor(y4,y4,x1);\
				x1 = ROTV(y4,r4);\
				\
				organize;\
				\
				load(x0,RC,0);\
				xor(y0,y0,x0);
				
#define load_one		\
				load(x1,in_temp,1);\
				load(x2,in_temp,2);\
				load(x3,in_temp,3);\
				load(x6,in_temp,4);\
				load(x0,in_temp,0);\
				\
				blend_32(a1,x1,x2,0x0C);\
				blend_32(a0,x1,x2,0x03);\
				blend_32(x6,x6,x3,0xF0);\
				permut_2x128_avx(x2,x2,x3,0x21);\
				blend_32(x4,a1,x3,0x30);\
				permut_64(x6,x6,0x03);\
				\
				permut_64(x1,a0,0x39);\
				set_zero_256(x3);\
				blend_32(x4,x4,x3,0xC0);\
				blend_32(x3,x6,x3,0xF0);\
				xor(y6,y6,x4);\
				xor(y0,y0,x0);\
				xor(y1,y1,x1);\
				extract_128(m0,y6,0);\
				\
				extract_128(m1,y6,1);\
				\
				xor_128(m0,m0,m1);\
				cast_to_256(x1,m0);\
				permut_64(x1,x1,0x55);\
				xor(c1,y0,y4);\
				xor(y2,y2,x2);\
				xor(c1,c1,y1);\
				broadcast_128L_256(x0,m0);\
				\
				xor(x0,x0,y5);\
				xor(x0,x0,x1);\
				xor(c1,c1,y2);\
				xor(y3,y3,x3);\
				xor(c1,c1,y3);\
				\
				xor(x1,x0,ROT(c1,0x01));\
				blend_32(x2,c1,x0,0x0C);\
				permut_64(x1,x1,0x55);\
				permut_64(x2,x2,0x1E);\
				xor(x2,c1,ROT(x2,0x01));\
				permut_64(x3,x2,0xFF);\
				\
				xor(y5,y5,x3);\
				xor(y6,y6,x3);\
				x5 = ROTV(y6,r5);\
				blend_32(x1,x2,x1,0xC0);\
				permut_64(x1,x1,0x93);\
				\
				xor(y0,y0,x1);\
				blend_32(x0,y0,y5,0x03);\
				x4 = ROTV(x0,r0);\
				xor(y1,y1,x1);\
				x6 = ROTV(y1,r1);\
				xor(y2,y2,x1);\
				x2 = ROTV(y2,r2);\
				xor(y3,y3,x1);\
				x3 = ROTV(y3,r3);\
				xor(y4,y4,x1);\
				x1 = ROTV(y4,r4);\
				\
				organize;\
				\
				load(x0,RC,0);\
				xor(y0,y0,x0);
				
#define load_two		\
				load(x0,in_temp,0);\
				load(x1,in_temp,1);\
				load(x2,in_temp,2);\
				load(x3,in_temp,3);\
				load(x4,in_temp,4);\
				\
				blend_32(x5,x1,x2,0x30);\
				blend_32(x0,x0,x1,0x03);\
				blend_32(x6,x2,x3,0x3F);\
				blend_32(x5,x5,x3,0xC0);\
				permut_2x128_avx(x1,x1,x2,0x21);\
				permut_64(x0,x0,0x39);\
				permut_64(x6,x6,0x93);\
				permut_64(x5,x5,0x39);\
				xor(y0,y0,x0);\
				xor(y1,y1,x1);\
				xor(c1,y4,y0);\
				set_zero_256(a1);\
				blend_32(x4,x4,a1,0xF0);\
				blend_32(x5,x5,a1,0xC0);\
				\
				xor(c1,c1,y1);\
				xor(y2,y2,x6);\
				xor(y6,y6,x5);\
				extract_128(m0,y6,0);\
				extract_128(m1,y6,1);\
				xor_128(m0,m0,m1);\
				cast_to_256(x1,m0);\
				broadcast_128L_256(x0,m0);\
				permut_64(x1,x1,0x55);\
				xor(y3,y3,x4);\
				xor(c1,c1,y2);\
				xor(c1,c1,y3);\
				xor(x0,x0,y5);\
				xor(x0,x0,x1);\
				\
				xor(x1,x0,ROT(c1,0x01));\
				blend_32(x2,c1,x0,0x0C);\
				permut_64(x1,x1,0x55);\
				permut_64(x2,x2,0x1E);\
				xor(x2,c1,ROT(x2,0x01));\
				permut_64(x3,x2,0xFF);\
				\
				xor(y5,y5,x3);\
				xor(y6,y6,x3);\
				x5 = ROTV(y6,r5);\
				blend_32(x1,x2,x1,0xC0);\
				permut_64(x1,x1,0x93);\
				\
				xor(y0,y0,x1);\
				blend_32(x0,y0,y5,0x03);\
				x4 = ROTV(x0,r0);\
				xor(y1,y1,x1);\
				xor(y2,y2,x1);\
				xor(y3,y3,x1);\
				xor(y4,y4,x1);\
				x6 = ROTV(y1,r1);\
				x2 = ROTV(y2,r2);\
				x3 = ROTV(y3,r3);\
				x1 = ROTV(y4,r4);\
				\
				organize;\
				\
				load(x0,RC,0);\
				xor(y0,y0,x0);
				
#define load_three		\
				load(x0,in_temp,0);\
				load(x1,in_temp,1);\
				load(x4,in_temp,3);\
				load(x5,in_temp,4);\
				\
				load(x2,in_temp,2);\
				blend_32(x3,x2,x1,0xC0);\
				blend_32(a1,x1,x2,0xC0);\
				blend_32(a1,a1,x5,0x03);\
				permut_2x128_avx(x0,x0,x1,0x21);\
				permut_64(x3,x3,0x93);\
				permut_64(a1,a1,0x0E);\
				permut_64(x5,x5,0x99);\
				xor(y0,y0,x0);\
				xor(y2,y2,x4);\
				xor(x4,y0,y2);\
				xor(y1,y1,x3);\
				xor(x4,x4,y1);\
				\
				xor(x4,x4,y4);\
				set_zero_256(x6);\
				blend_32(x1,a1,x6,0xC0);\
				xor(y6,y6,x1);\
				\
				extract_128(m1,y6,0);\
				extract_128(m2,y6,1);\
				xor_128(m1,m1,m2);\
				m2 = _mm_shuffle_epi32(m1,0x4E);\
				xor_128(m1,m1,m2);\
				cast_to_256(x2,m1);\
				xor(x2,x2,y5);\
				\
				extract_128(m0,x5,0);\
				insert_128(x6,x6,m0,0);\
				xor(y3,y3,x6);\
				xor(x1,x4,y3);\
				\
				xor(x0,x2,ROT(x1,0x01));\
				blend_32(x3,x1,x2,0x0C);\
				permut_64(x3,x3,0x1E);\
				xor(x3,x1,ROT(x3,0x01));\
				permut_64(x0,x0,0x55);\
				permut_64(x1,x3,0xFF);\
				xor(y5,y5,x1);\
				xor(y6,y6,x1);\
				x5 = ROTV(y6,r5);\
				permut_64(x2,x3,0x90);\
				blend_32(x0,x0,x2,0xFC);\
				\
				xor(y0,y0,x0);\
				blend_32(x1,y0,y5,0x03);\
				x4 = ROTV(x1,r0);\
				xor(y1,y1,x0);\
				xor(y4,y4,x0);\
				xor(y2,y2,x0);\
				xor(y3,y3,x0);\
				x6 = ROTV(y1,r1);\
				x1 = ROTV(y4,r4);\
				x2 = ROTV(y2,r2);\
				x3 = ROTV(y3,r3);\
				\
				organize;\
				load(x0,RC,0);\
				xor(y0,y0,x0);
				
				
#define load_four		\
				load(x0,in_temp,0);\
				load(x1,in_temp,1);\
				load(x2,in_temp,2);\
				load(x3,in_temp,3);\
				load(x4,in_temp,4);\
				\
				blend_32(x5,x1,x3,0x03);\
				blend_32(x0,x1,x0,0xC0);\
				xor(y1,y1,x2);\
				extract_128(m0,x4,1);\
				cast_to_256(x1,m0);\
				blend_32(x6,x3,x4,0x03);\
				blend_32(x5,x5,x4,0x0C);\
				permut_64(x0,x0,0x93);\
				permut_64(x3,x6,0x39);\
				permut_64(x5,x5,0x93);\
				blend_32(x5,x5,x1,0xC0);\
				xor(y0,y0,x0);\
				xor(y3,y3,x1);\
				xor(c1,y1,y4);\
				xor(y6,y6,x5);\
				\
				xor(x6,y3,y0);\
				extract_128(m1,y6,0);\
				extract_128(m2,y6,1);\
				xor(y2,y2,x3);\
				xor(c1,c1,x6);\
				xor_128(m1,m1,m2);\
				m2 = _mm_shuffle_epi32(m1,0x4E);\
				xor_128(m1,m1,m2);\
				cast_to_256(x2,m1);\
				xor(x0,x2,y5);\
				xor(c1,c1,y2);\
				\
				xor(x1,x0,ROT(c1,0x01));\
				blend_32(x2,c1,x0,0x0C);\
				permut_64(x1,x1,0x55);\
				permut_64(x2,x2,0x1E);\
				xor(x2,c1,ROT(x2,0x01));\
				permut_64(x3,x2,0xFF);\
				\
				xor(y5,y5,x3);\
				xor(y6,y6,x3);\
				x5 = ROTV(y6,r5);\
				blend_32(x1,x2,x1,0xC0);\
				permut_64(x1,x1,0x93);\
				\
				xor(y0,y0,x1);\
				blend_32(x0,y0,y5,0x03);\
				x4 = ROTV(x0,r0);\
				xor(y1,y1,x1);\
				x6 = ROTV(y1,r1);\
				xor(y2,y2,x1);\
				x2 = ROTV(y2,r2);\
				xor(y3,y3,x1);\
				x3 = ROTV(y3,r3);\
				xor(y4,y4,x1);\
				x1 = ROTV(y4,r4);\
				\
				organize;\
				\
				load(x0,RC,0);\
				xor(y0,y0,x0);
				

#define organize 		\
				blend_32(y1,x4,x6,0xCC);\
				blend_32(x0,x2,x3,0xCC);\
				blend_32(y6,x4,x1,0x33);\
				permut_64(y4,x5,0x4C);\
				extract_128(m0,y1,1);\
				\
				blend_32(y2,x0,y1,0x0F);\
				blend_32(x0,x0,y6,0xF0);\
				blend_32(a1,y4,x6,0x3F);\
				permut_64(y1,y2,0x39);\
				\
				blend_32(y5,y0,y1,0xFC);\
				blend_32(y0,y0,y2,0xFC);\
				blend_32(c1,x3,x1,0xCC);\
				permut_64(y5,y5,0x39);\
				\
				a0 = y0;\
				xor_and_not(y0,y1,y5);\
				\
				blend_32(y5,x6,x2,0xCC);\
				blend_32(y2,x0,a1,0xC0);\
				blend_32(x2,c1,y5,0x0F);\
				blend_32(x3,x2,x5,0x03);\
				permut_64(y2,y2,0x93);\
				permut_64(x3,x3,0x39);\
				\
				blend_32(y1,y2,x0,0xC0);\
				blend_32(y3,x2,x3,0xC0);\
				permut_64(y1,y1,0x93);\
				\
				xor_and_not(y1,y2,x0);\
				\
				permut_64(y3,y3,0x93);\
				\
				xor_and_not(y3,x2,x3);\
				blend_32(x2,y5,y6,0x0F);\
				permut_64(a1,a1,0x2D);\
				permut_64(y6,x1,0xC9);\
				permut_64(x0,x2,0x4B);\
				\
				blend_32(x0,x0,y4,0x0C);\
				blend_32(y6,y6,x4,03);\
				blend_32(x2,x0,x2,0xF3);\
				permut_64(x4,x4,0x9E);\
				permut_64(x2,x2,0x1E);\
				\
				blend_32(y2,x2,x0,0xC0);\
				blend_32(x1,x4,x3,0xC0);\
				blend_32(x3,x4,x5,0x3F);\
				permut_64(y2,y2,0x93);\
				blend_32(x1,x1,a0,0x03);\
				\
				xor_and_not(y2,x2,x0);\
				blend_32(x5,x3,c1,0x0F);\
				permut_64(x5,x5,0xD2);\
				xor_and_not(y6,x1,a1);\
				blend_32(x4,x5,x6,0xC0);\
				\
				extract_128(m1,x5,0);\
				permut_64(x4,x4,0x93);\
				cast_to_256(y4,m0);\
				insert_128(y4,y4,m1,1);\
				\
				xor_and_not(y4,x4,x5);\
				xor_and_not(x4,x5,x6);\
				permut_64(y5,x4,0xFF);\
				
				
				
#define new_round  		\
				extract_128(m0,y6,0);\
				extract_128(m1,y6,1);\
				xor(x0,y0,y1);\
				xor(x1,y2,y3);\
				xor_128(m0,m0,m1);\
				\
				m1 = _mm_shuffle_epi32(m0,0x4E);\
				xor(x0,x0,y4);\
				xor_128(m1,m1,m0);\
				cast_to_256(x2,m1);\
				xor(x1,x0,x1);\
				xor(x2,x2,y5);\
				blend_32(x3,x1,x2,0x0C);\
				\
				permut_64(x3,x3,0x1E);\
				xor(x3,x1,ROT(x3,0x01));\
				xor(x0,x2,ROT(x1,0x01));\
				\
				permut_64(x5,x3,0xFF);\
				permut_64(x0,x0,0x55);\
				xor(y6,y6,x5);\
				xor(y5,y5,x5);\
				permut_64(x2,x3,0x90);\
				blend_32(a1,x0,x2,0xFC);\
				\
				xor(y4,y4,a1);\
				xor(y1,y1,a1);\
				xor(y2,y2,a1);\
				xor(y0,y0,a1);\
				xor(y3,y3,a1);\
				y4 = ROTV(y4,r15);\
				\
				blend_32(x0,y0,y6,0x03);\
				blend_32(x5,y1,y2,0x33);\
				blend_32(x2,y2,y6,0x30);\
				blend_32(x1,y1,y6,0x0C);\
				\
				permut_64(x4,x0,0x39);\
				permut_64(y1,x1,0x1E);\
				permut_64(y2,x2,0x4B);\
				\
				x4 = ROTV(x4,r11);\
				blend_32(x6,y5,y3,0xC0);\
				blend_32(x3,y3,y6,0xC0);\
				blend_32(x5,x5,x6,0xC3);\
				\
				y1 = ROTV(y1,r12);\
				\
				\
				y5 = ROTV(x5,r10);\
				\
				\
				blend_32(a0,y5,y0,0x03);\
				permut_64(x6,y5,0x39);\
				y2 = ROTV(y2,r13);\
				x1 = (_mm256_andnot_si256(y1,y2));\
				blend_32(a1,a0,x6,0xFC);\
				\
				permut_64(y3,x3,0x93);\
				\
				permut_64(a1,a1,0x39);\
				y3 = ROTV(y3,r14);\
				\
				x1 = _mm256_xor_si256(x4,x1);\
				x2 = (_mm256_andnot_si256(y2,y3));\
				x2 = _mm256_xor_si256(y1,x2);\
				\
				unpackLO_128_64(c1,x1,x2);\
				\
				unpackHI_128_64(x1,x1,x2);\
				x3 = _mm256_xor_si256(y2,_mm256_andnot_si256(y3,y4));\
				x5 = _mm256_xor_si256(y4,_mm256_andnot_si256(x4,y1));\
				x4 = _mm256_xor_si256(y3,_mm256_andnot_si256(y4,x4));\
				\
				unpackLO_128_64(x2,x3,x4);\
				extract_128(m0,c1,1);\
				cast_to_256(y1,m0);\
				blend_32(y1,y1,x2,0xF0);\
				\
				unpackHI_128_64(x3,x3,x4);\
				y0 = (_mm256_andnot_si256(x6,a1));\
				y0 = _mm256_xor_si256(a0,y0);\
				\
				permut_2x128_avx(y2,c1,x2,0x20);\
				x6 = _mm256_xor_si256(y5,_mm256_andnot_si256(a0,x6));\
				\
				extract_128(m0,x1,1);\
				cast_to_256(y3,m0);\
				blend_32(y3,y3,x3,0xF0);\
				\
				permut_2x128_avx(y4,x1,x3,0x20);\
				y5 = _mm256_shuffle_epi32(x5,0xEE);\
				\
				permut_64(y6,x5,0xC9);\
				blend_32(y6,y6,x6,0x03);\
				
				
								

#ifdef __INTEL_COMPILER
#define	other_rounds  		for(q=1;q<24;q++){\
				  round;\
				  load(x0,RC,q);\
				  xor(y0,y0,x0);}
				  
#else
#define other_rounds	        new_round;/*2 roudn*/\
				load(x0,RC,1);\
				xor(y0,y0,x0);\
				\
				new_round;/*3 round*/\
				load(x0,RC,2);\
				xor(y0,y0,x0);\
				\
				\
				new_round;/*4 round*/\
				load(x0,RC,3);\
				xor(y0,y0,x0);\
				\
				new_round;/*5 round*/\
				load(x0,RC,4);\
				xor(y0,y0,x0);\
				\
				new_round;/*6 round*/\
				load(x0,RC,5);\
				xor(y0,y0,x0);\
				\
				new_round;/*7 round*/\
				load(x0,RC,6);\
				xor(y0,y0,x0);\
				\
				new_round/*8 round*/\
				load(x0,RC,7);\
				xor(y0,y0,x0);\
				\
				new_round;/*9 round*/\
				load(x0,RC,8);\
				xor(y0,y0,x0);\
				\
				new_round;/*10 round*/\
				load(x0,RC,9);\
				xor(y0,y0,x0);\
				\
				new_round;/*11 round*/\
				load(x0,RC,10);\
				xor(y0,y0,x0);\
				\
				new_round;/*12 round*/\
				load(x0,RC,11);\
				xor(y0,y0,x0);\
				\
				new_round;/*13 round*/\
				load(x0,RC,12);\
				xor(y0,y0,x0);\
				\
				new_round;/*14 round*/\
				load(x0,RC,13);\
				xor(y0,y0,x0);\
				\
				new_round;/*15 round*/\
				load(x0,RC,14);\
				xor(y0,y0,x0);\
				\
				new_round;/*16 round*/\
				load(x0,RC,15);\
				xor(y0,y0,x0);\
				\
				new_round;/*17 round*/\
				load(x0,RC,16);\
				xor(y0,y0,x0);\
				\
				new_round;/*18 round*/\
				load(x0,RC,17);\
				xor(y0,y0,x0);\
				\
				new_round;/*19 round*/\
				load(x0,RC,18);\
				xor(y0,y0,x0);\
				\
				new_round;/*20 round*/\
				load(x0,RC,19);\
				xor(y0,y0,x0);\
				\
				new_round;/*21 round*/\
				load(x0,RC,20);\
				xor(y0,y0,x0);\
				\
				new_round;/*22 round*/\
				load(x0,RC,21);\
				xor(y0,y0,x0);\
				\
				new_round;/*23 round*/\
				load(x0,RC,22);\
				xor(y0,y0,x0);\
				\
				new_round;/*24 round*/\
				load(x0,RC,23);\
				xor(y0,y0,x0);
#endif				
			
#define xor_and_not(x,y,z)	x = _mm256_xor_si256(x,(_mm256_andnot_si256(y,z)));


int keccak(const uint8_t *in, int inlen, uint8_t *md){

const ALIGN uint64_t RC[96] = {
    0x0000000000000001,0x0,0x0,0x0,
    0x0000000000008082,0x0,0x0,0x0,
    0x800000000000808A,0x0,0x0,0x0,
    0x8000000080008000,0x0,0x0,0x0,
    0x000000000000808B,0x0,0x0,0x0,
    0x0000000080000001,0x0,0x0,0x0,
    0x8000000080008081,0x0,0x0,0x0,
    0x8000000000008009,0x0,0x0,0x0,
    0x000000000000008A,0x0,0x0,0x0,
    0x0000000000000088,0x0,0x0,0x0,
    0x0000000080008009,0x0,0x0,0x0,
    0x000000008000000A,0x0,0x0,0x0,
    0x000000008000808B,0x0,0x0,0x0,
    0x800000000000008B,0x0,0x0,0x0,
    0x8000000000008089,0x0,0x0,0x0,
    0x8000000000008003,0x0,0x0,0x0,
    0x8000000000008002,0x0,0x0,0x0,
    0x8000000000000080,0x0,0x0,0x0,
    0x000000000000800A,0x0,0x0,0x0,
    0x800000008000000A,0x0,0x0,0x0,
    0x8000000080008081,0x0,0x0,0x0,
    0x8000000000008080,0x0,0x0,0x0,
    0x0000000080000001,0x0,0x0,0x0,
    0x8000000080008008,0x0,0x0,0x0
  };
  __m256i x0,x1,x2,x3,x4,x5,x6;
  __m256i y0,y1,y2,y3,y4,y5,y6;
  __m256i a0,a1,c1;
  __m128i m0,m1,m2,m3;
  char *in_temp = (char*)in;
  int k=1,l=0,q,i;
  ALIGN uint8_t temp[144];
  
    
   const __m256i r0 = _mm256_set_epi64x(0x1C,0x3E,0x1,0xE);
   const __m256i r1 = _mm256_set_epi64x(0x37,0x6,0x2C,0x24);
   const __m256i r2 = _mm256_set_epi64x(0x19,0x2B,0xA,0x3);
   const __m256i r3 = _mm256_set_epi64x(0x15,0xF,0x2D,0x29);
   const __m256i r4 = _mm256_set_epi64x(0x38,0x3D,0x2,0x12);
   const __m256i r5 = _mm256_set_epi64x(0x8,0x27,0x14,0x1B);
   
   const __m256i r10 = _mm256_set_epi64x(0x15,0x2B,0x2C,0xE);
   const __m256i r11 = _mm256_set_epi64x(0x1B,0x1C,0x3E,0x1);
   const __m256i r12 = _mm256_set_epi64x(0x24,0x14,0x37,0x6);
   const __m256i r13 = _mm256_set_epi64x(0xA,0x3,0x27,0x19);
   const __m256i r14 = _mm256_set_epi64x(0xF,0x2D,0x29,0x8);
   const __m256i r15 = _mm256_set_epi64x(0x38,0x3D,0x2,0x12);
  
  memset(temp, 0, 144*sizeof(uint8_t));
  zero_state;

  while(inlen-l >= RSIZ){
    
    load_one;

    other_rounds;
    inlen -= 128;
    l=8;
    in_temp += 128;
    
    if(inlen-l < RSIZ){
	  k=2;
	  break;
    }
    
    load_two;
    other_rounds;
    l=16;
    inlen -= 128;
    in_temp += 128;
    
    if(inlen-l < RSIZ){
	  k=3;
	  break;
    }
    
    load_three;
    other_rounds;
    l=24;
    inlen -= 128;
    in_temp += 128;
    
    if(inlen-l < RSIZ){
	  k=4;
	  break;
    }    
    
    load_four;
    other_rounds;
    inlen -= 160;
    in_temp +=160;
    l=0;
  }
 
  switch(k){
      case 2:
	inlen-=8;
	in_temp += 8;
	break;
      case 3:
	inlen-=16;
	in_temp += 16;
	break;
      case 4:
	inlen-=24;
	in_temp += 24;
	break;
      case 1:
	break;
      default:
	printf("ERRO!!!");
	return 0;
    }

    memcpy(temp, in_temp, inlen);
    temp[inlen++] = 0x06;
    temp[RSIZ - 1] |= 0x80;
    
    load_last;
    other_rounds;
    
    store(md,0,y0);
    
    return 0;

}
