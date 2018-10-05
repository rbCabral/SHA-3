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


#define xor(X,Y)		_mm_xor_si128(X,Y)
#define and(X,Y)		_mm_and_si128(X,Y)
#define load(M,D) 		_mm_load_si128((__m128i*)M+D)
#define store(X,D,Y)		_mm_store_si128((__m128i*)X+D,Y)
#define unpackLO_64(X,Y)	_mm_unpacklo_epi64(X,Y)
#define unpackHI_64(X,Y)	_mm_unpackhi_epi64(X,Y)
#define set_zero_128()		_mm_setzero_si128()
#define set_zero_HI(X)		_mm_move_epi64(X)
#define alignr(X,Y,D)	 	_mm_alignr_epi8 (X,Y,D)
//#define shiftR(Y,D)		_mm_srli_si128(Y,D)
//#define shiftL(Y,D)		_mm_slli_si128(Y,D)
#define shuffle_32(X,D)		_mm_shuffle_epi32(X,D)
#define blend_32(X,Y,D)		_mm_blend_epi32(X,Y,D)
//#define blend_LH(X,Y)		xor(and(X,bb),and(Y,aa));
//#define blend_HL(X,Y)		xor(and(X,aa),and(Y,bb));
#define and_not(X,Y)		_mm_andnot_si128(X,Y)


#define ROT(Y,D) 		_mm_xor_si128(_mm_add_epi64(Y,Y),_mm_srli_epi64(Y,(63)))
#define ROTV(x,d) 		_mm_xor_si128(_mm_sllv_epi64(x,d),_mm_srlv_epi64(x,_mm_sub_epi64(_mm_set1_epi64x(0x40),d)))

#ifdef __INTEL_COMPILER
#define BARRIER __memory_barrier()
#else
#define BARRIER asm volatile("" ::: "memory")
#endif



#define zero_state		y0 = set_zero_128();\
				y1 = set_zero_128();\
				y2 = set_zero_128();\
				y3 = set_zero_128();\
				y4 = set_zero_128();\
				y5 = set_zero_128();\
				y6 = set_zero_128();\
				y7 = set_zero_128();\
				y8 = set_zero_128();\
				y9 = set_zero_128();\
				y10 = set_zero_128();\
				y11 = set_zero_128();\
				y12 = set_zero_128();\

#define load_one_224(M)		\
				x0 = load(M,0);\
				x1 = load(M,1);\
				y0 = xor(x0,y0);\
				y1 = xor(x1,y1);\
				x2 = load(M,2);\
				x3 = load(M,3);\
				x4 = load(M,4);\
				\
				x0 = alignr(x3,x2,8);\
				x1 = alignr(x4,x3,8);\
				x3 = blend_32(x2,x4,0xC);\
				\
				x4 = load(M,7);\
				x2 = load(M,8);\
				y2 = xor(x0,y2);\
				y3 = xor(x1,y3);\
				\
				x1 = load(M,9);\
				x0 = alignr(x1,x2,8);\
				x0 = set_zero_HI(x0);\
				y4 = xor(x3,y4);\
				\
				y8 = xor(x0,y8);\
				x2 = alignr(x2,x4,8);\
				x4 = set_zero_HI(x4);\
				\
				x0 = load(M,5);\
				x1 = load(M,6);\
				y7 = xor(x2,y7);\
				y9 = xor(x4,y9);\
				y5 = xor(x0,y5);\
				y6 = xor(x1,y6);\

#define load_one_256(M)		\
				x0 = load(M,0);\
				x1 = load(M,1);\
				y0 = xor(x0,y0);\
				y1 = xor(x1,y1);\
				x2 = load(M,2);\
				x3 = load(M,3);\
				x4 = load(M,4);\
				\
				x0 = alignr(x3,x2,8);\
				x1 = alignr(x4,x3,8);\
				x3 = blend_32(x2,x4,0xC);\
				\
				x4 = load(M,7);\
				x2 = load(M,8);\
				y2 = xor(x0,y2);\
				y3 = xor(x1,y3);\
				y4 = xor(x3,y4);\
				\
				x2 = alignr(x2,x4,8);\
				x4 = set_zero_HI(x4);\
				\
				x0 = load(M,5);\
				x1 = load(M,6);\
				y7 = xor(x2,y7);\
				y9 = xor(x4,y9);\
				y5 = xor(x0,y5);\
				y6 = xor(x1,y6);\

#define load_one_384(M)		\
				x0 = load(M,0);\
				x1 = load(M,1);\
				y0 = xor(x0,y0);\
				y1 = xor(x1,y1);\
				x2 = load(M,2);\
				x3 = load(M,3);\
				x4 = load(M,4);\
				\
				x0 = alignr(x3,x2,8);\
				x1 = alignr(x4,x3,8);\
				x3 = blend_32(x2,x4,0xC);\
				\
				x4 = load(M,5);\
				x2 = load(M,6);\
				y2 = xor(x0,y2);\
				y3 = xor(x1,y3);\
				y4 = xor(x3,y4);\
				\
				y5 = xor(x4,y5);\
				x2 = set_zero_HI(x2);\
				y6 = xor(x2,y6);\
				
#define load_one_512(M)		\
				x0 = load(M,0);\
				x1 = load(M,1);\
				y0 = xor(x0,y0);\
				y1 = xor(x1,y1);\
				x2 = load(M,2);\
				x3 = load(M,3);\
				x4 = load(M,4);\
				\
				x0 = alignr(x3,x2,8);\
				x1 = alignr(x4,x3,8);\
				x3 = set_zero_HI(x2);\
				\
				y2 = xor(x0,y2);\
				y3 = xor(x1,y3);\
				y4 = xor(x3,y4);\
								
#define load_two_256(M)		\
				\
				x0 = load(M,0);\
				x1 = load(M,1);\
				x2 = load(M,2);\
				x3 = load(M,3);\
				x4 = load(M,4);\
				\
				\
				y2 = xor(y2,x3);\
				y3 = xor(y3,x4);\
				x0 = alignr(x1,x0,8);\
				x1 = alignr(x2,x1,8);\
				y0 = xor(y0,x0);\
				y1 = xor(y1,x1);\
				\
				x0 = load(M,5);\
				x1 = load(M,6);\
				x3 = load(M,7);\
				x4 = load(M,8);\
				\
				\
				x2 = alignr(x0,x2,8);\
				x0 = alignr(x1,x0,8);\
				x1 = alignr(x3,x1,8);\
				x3 = _mm_srli_si128(x3,8);\
				y4 = xor(y4,x2);\
				y5 = xor(y5,x0);\
				y6 = xor(y6,x1);\
				y9 = xor(y9,x3);\
				y7 = xor(y7,x4);\

#define load_two_384(M)		\
				\
				x0 = load(M,0);\
				x1 = load(M,1);\
				x2 = load(M,2);\
				x3 = load(M,3);\
				x4 = load(M,4);\
				\
				y2 = xor(y2,x3);\
				y3 = xor(y3,x4);\
				x0 = alignr(x1,x0,8);\
				x1 = alignr(x2,x1,8);\
				y0 = xor(y0,x0);\
				y1 = xor(y1,x1);\
				\
				x0 = load(M,5);\
				x1 = load(M,6);\
				x2 = alignr(x0,x2,8);\
				x0 = alignr(x1,x0,8);\
				\
				x1 = _mm_srli_si128(x1,8);\
				y4 = xor(y4,x2);\
				y5 = xor(y5,x0);\
				y6 = xor(y6,x1);\

#define load_two_512(M)		\
				x0 = load(M,0);\
				x1 = load(M,1);\
				x2 = load(M,2);\
				x3 = load(M,3);\
				x4 = load(M,4);\
				\
				y2 = xor(y2,x3);\
				y3 = xor(y3,x4);\
				x0 = alignr(x1,x0,8);\
				x1 = alignr(x2,x1,8);\
				y0 = xor(y0,x0);\
				y1 = xor(y1,x1);\
				\
				x2 = _mm_srli_si128(x2,8);\
				y4 = xor(y4,x2);\
		

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
				x4 = blend_32(y0,y12,0x3);\
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
				y0 = xor(y1,and_not(y4,y5));\
				y7 = ROTV(y7,r7);\
				y3 = blend_32(y3,y6,0xC);\
				y6 = unpackHI_64(y6,y9);\
				\
				y2 = xor(y4,and_not(y5,y7));\
				y4 = xor(y11,and_not(y1,y4));\
				y9 = xor(y5,and_not(y7,y11));\
				y7 = xor(y7,and_not(y11,y1));\
				\
				/*I'm not sure if this is the best way to compute these values*/\
				x5 = alignr(y3,x0,8);\
				x0 = blend_32(x0,y10,0x3);\
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
				/*x2 = _mm_slli_si128(y2,8);\
				x5 = _mm_slli_si128(x4,8);\
				y10 = blend_32(y0,x2,0xC);\
				y7  = blend_32(x3,x5,0xC);\*/\
				\
				y10 = unpackLO_64(y0,y2);\
				y7 = unpackLO_64(x3,x4);\
				\
				y12 = shuffle_32(y4,0x44);\
				y4 = unpackHI_64(x1,y4);\
				y2 = unpackHI_64(y0,y2);\
				y0 = unpackHI_64(x3,x4);\
				\


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
			

void print_128(__m128i x){
     uint64_t c[2];
     store(c,0,x);
     printf("\n 1 = %8.16lx\t  0 = %8.16lx\t\n", c[1],c[0]);
}			
			
int keccak(const uint8_t *in, int inlen, uint8_t *md, int r){

const ALIGN uint64_t RC[48] = {
    0x0000000000000001,0x0, 0x0000000000008082,0x0,
    0x800000000000808A,0x0, 0x8000000080008000,0x0,
    0x000000000000808B,0x0, 0x0000000080000001,0x0,
    0x8000000080008081,0x0, 0x8000000000008009,0x0,
    0x000000000000008A,0x0, 0x0000000000000088,0x0,
    0x0000000080008009,0x0, 0x000000008000000A,0x0,
    0x000000008000808B,0x0, 0x800000000000008B,0x0,
    0x8000000000008089,0x0, 0x8000000000008003,0x0,
    0x8000000000008002,0x0, 0x8000000000000080,0x0,
    0x000000000000800A,0x0, 0x800000008000000A,0x0,
    0x8000000080008081,0x0, 0x8000000000008080,0x0,
    0x0000000080000001,0x0, 0x8000000080008008,0x0
  };
  
  
  const __m128i r2 = _mm_set_epi64x(44,36);
  const __m128i r4 = _mm_set_epi64x(20,27);
  const __m128i r6 = _mm_set_epi64x(25,43);
  const __m128i r5 = _mm_set_epi64x(10,3);
  const __m128i r8 = _mm_set_epi64x(21,15);
  const __m128i r3 = _mm_set_epi64x(55,6);
  const __m128i r11= _mm_set_epi64x(56,61);
  const __m128i r0 = _mm_set_epi64x(1,14);
  const __m128i r10 = _mm_set_epi64x(2,18);
  const __m128i r9 = _mm_set_epi64x(8,39);
  const __m128i r1 = _mm_set_epi64x(28,62);
  const __m128i r7 = _mm_set_epi64x(45,41);

  __m128i y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12;
  __m128i x0,x1,x2,x3,x4,x5,x6;
  
  const uint8_t *in_temp = in;
  int q,it=0,i;
  ALIGN uint8_t temp[144];
       
  memset(temp, 0, 144*sizeof(uint8_t));
  
  zero_state;
  
  if(inlen >= r ){
    it = inlen / r;
  }

  switch(r){
    
    case 144:
      for(i=0;i<it;i++){
	load_one_224(in_temp);
	other_rounds(q);
	in_temp +=144;
      }

      inlen = inlen - it * r;
      memcpy(temp, in_temp, inlen);
      temp[inlen++] = 0x06;
      temp[r - 1] |= 0x80;
  
      load_one_224(temp);
      other_rounds(q);
      store(md,0,y0);
      store(md,1,y1);
      break;
    
    case 136:
      for(i=0;i<it/2;i++){
	load_one_256(in_temp);
	other_rounds(q);
      
	in_temp += 128;
	
	load_two_256(in_temp);
	other_rounds(q);
	in_temp +=144;
      }
    
      if(it%2 != 0 ){
	load_one_256(in_temp);
	other_rounds(q);
	in_temp += 136;
      }
      inlen = inlen - it * r;
      memcpy(temp, in_temp, inlen);
      temp[inlen++] = 0x06;
      temp[r - 1] |= 0x80;
      load_one_256(temp);
      other_rounds(q);
      store(md,0,y0);
      store(md,1,y1);
      break;
	
    case 104:
      for(i=0;i<it/2;i++){
	load_one_384(in_temp);
	other_rounds(q);
	in_temp += 96;
    
	load_two_384(in_temp);
	other_rounds(q);
	in_temp +=112;
      }
  
      if(it%2 != 0 ){
	load_one_384(in_temp);
	other_rounds(q);
	in_temp += 104;
      }

      inlen = inlen - it * r;
      memcpy(temp, in_temp, inlen);
      temp[inlen++] = 0x06;
      temp[r - 1] |= 0x80;
      load_one_384(temp);
      other_rounds(q);
      store(md,0,y0);
      store(md,1,y1);
      store(md,2,unpackLO_64(y4,y2));
      break;
    case 72:
      for(i=0;i<it/2;i++){
	load_one_512(in_temp);
	other_rounds(q);
	in_temp += 64;
	
	load_two_512(in_temp);
	other_rounds(q);
	in_temp +=80;
      }
  
      if(it%2 != 0 ){
	load_one_512(in_temp);
	other_rounds(q);
	in_temp += 72;
      }

      inlen = inlen - it * r;
      memcpy(temp, in_temp, inlen);
      temp[inlen++] = 0x06;
      temp[r - 1] |= 0x80;
      load_one_512(temp);
      other_rounds(q);
      store(md,0,y0);
      store(md,1,y1);
      store(md,2,unpackLO_64(y4,y2));
      store(md,3,alignr(y3,y2,8));
      break;
    default:
      printf("The value of bit rate is not defined.\n");
  }
  
 
  return 0;
}
