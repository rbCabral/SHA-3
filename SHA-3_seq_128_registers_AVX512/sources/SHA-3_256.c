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
#define load(M,D) 		_mm_loadu_si128((__m128i*)M+D)
#define store(X,D,Y)		_mm_storeu_si128((__m128i*)X+D,Y)
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


//#define ROT(Y,D) 		_mm_xor_si128(_mm_add_epi64(Y,Y),_mm_srli_epi64(Y,(63)))
#define ROT(Y,D) 		_mm_rol_epi64(Y,D)
//#define ROTV(x,d) 		_mm_xor_si128(_mm_sllv_epi64(x,d),_mm_srlv_epi64(x,_mm_sub_epi64(_mm_set1_epi64x(0x40),d)))
#define ROTV(x,d) 		_mm_rolv_epi64(x,d)

#ifdef __INTEL_COMPILER
#define BARRIER __memory_barrier()
#else
#define BARRIER asm volatile("" ::: "memory")
#endif


#define zero_state		y[0] = set_zero_128();\
				y[1] = set_zero_128();\
				y[2] = set_zero_128();\
				y[3] = set_zero_128();\
				y[4] = set_zero_128();\
				y[5] = set_zero_128();\
				y[6] = set_zero_128();\
				y[7] = set_zero_128();\
				y[8] = set_zero_128();\
				y[9] = set_zero_128();\
				y[10] = set_zero_128();\
				y[11] = set_zero_128();\
				y[12] = set_zero_128();\

#define loadM(M)		\
				x0 = load(M,0);\
				x1 = load(M,1);\
				y[0] = xor(x0,y[0]);\
				y[1] = xor(x1,y[1]);\
				x2 = load(M,2);\
				x3 = load(M,3);\
				x4 = load(M,4);\
				\
				x0 = alignr(x3,x2,8);\
				x1 = alignr(x4,x3,8);\
				x3 = blend_32(x2,x4,0xC);\
				\
				y[2] = xor(x0,y[2]);\
				y[3] = xor(x1,y[3]);\
				y[4] = xor(x3,y[4]);\
				if(rsiz==104){\
					x4 = load(M,5);\
					x2 = load(M,6);\
					y[5] = xor(x4,y[5]);\
					x2 = set_zero_HI(x2);\
					y[6] = xor(x2,y[6]);\
				}else if (rsiz==136){\
					x4 = load(M,7);\
                	x2 = load(M,8);\
                	x2 = alignr(x2,x4,8);\
                	x4 = set_zero_HI(x4);\
                	\
                	x0 = load(M,5);\
                	x1 = load(M,6);\
                	y[7] = xor(x2,y[7]);\
                	y[9] = xor(x4,y[9]);\
                	y[5] = xor(x0,y[5]);\
                	y[6] = xor(x1,y[6]);\
                }else if(rsiz == 144){\
                	x4 = load(M,7);\
                	x2 = load(M,8);\
                	x1 = load(M,9);\
                	x0 = alignr(x1,x2,8);\
                	x0 = set_zero_HI(x0);\
                	\
	                y[8] = xor(x0,y[8]);\
	                x2 = alignr(x2,x4,8);\
	                x4 = set_zero_HI(x4);\
	                \
	                x0 = load(M,5);\
	                x1 = load(M,6);\
	                y[7] = xor(x2,y[7]);\
	                y[9] = xor(x4,y[9]);\
	                y[5] = xor(x0,y[5]);\
	                y[6] = xor(x1,y[6]);\
                }\


#define keccakF(y,rnd)\
	/*After the first part of the step theta we'll have the follows values in the follows registers\\
		x0 = c0 -- \
		x3 = c2 c1 \
		x2 = c4 c3 \
	 */\
\
	for(q=0;q<rnd;q+=2){\
		x2   = xor(y[4],y[9]);\
		x0 = _mm_ternarylogic_epi64(y[0],y[7],y[5],0x96);\
		x0 = _mm_ternarylogic_epi64(y[10],y[2],x0,0x96);\
		\
		x1 = _mm_ternarylogic_epi64(y[1],y[3],y[6],0x96);\
		x1 = _mm_ternarylogic_epi64(y[11],y[8],x1,0x96);\
		\
		x3 = shuffle_32(x2,0x4E);\
		x2 = _mm_ternarylogic_epi64(x3,x2,y[12],0x96);\
\
		/*	The secund part of the step theta\
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
		x4    = blend_32(y[0],y[12],0x3); /*1, 24*/\
		y[10] = xor(y[10],x0);\
		y[11] = ROTV(y[11],r11);/*19, 9*/\
\
		x5   = unpackLO_64(y[4],y[0]);\
		y[4] = unpackHI_64(y[3],y[4]);\
		x0   = ROTV(x4,r0);/*10, 4*/\
					\
		x3 = _mm_ternarylogic_epi64(x5,y[2],x2,0xD2);\
		y[10] = ROTV(y[10],r10); /*24, 14*/\
		x1    = alignr(x0,y[11],8);\
\
		y[9] = ROTV(y[9],r9);/*13, 22*/\
		x4 = _mm_ternarylogic_epi64(y[2],x2,y[8],0xD2);\
		x2 = _mm_ternarylogic_epi64(x2,y[8],x1,0xD2);\
\
		y[5] = unpackLO_64(y[9],y[5]);\
		y[8] = _mm_ternarylogic_epi64(y[8],x1,x5,0xD2);\
		x1 = _mm_ternarylogic_epi64(x1,x5,y[2],0xD2);\
\
		y[1]  = ROTV(y[1],r1); /*5, 20*/\
		y[11] = alignr(y[11],y[10],8);\
\
		y[0] = _mm_ternarylogic_epi64(y[1],y[4],y[5],0xD2);\
		y[7] = ROTV(y[7],r7); /*8, 23*/\
		y[3] = blend_32(y[3],y[6],0xC);\
		y[6] = unpackHI_64(y[6],y[9]);\
					\
		y[2] = _mm_ternarylogic_epi64(y[4],y[5],y[7],0xD2);\
		y[4] = _mm_ternarylogic_epi64(y[11],y[1],y[4],0xD2);\
		y[9] = _mm_ternarylogic_epi64(y[5],y[7],y[11],0xD2);\
		y[7] = _mm_ternarylogic_epi64(y[7],y[11],y[1],0xD2);\
\
		/*I'm not sure if this is the best way to compute these values*/\
		x5    = alignr(y[3],x0,8);\
		x0    = blend_32(x0,y[10],0x3);\
		y[10] = alignr(y[10],y[6],8);\
\
		y[1]  = unpackHI_64(x2,y[8]);\
		y[11] = unpackLO_64(y[9],y[7]);\
		y[8]  = unpackLO_64(x2,y[8]);\
\
		x2   = _mm_ternarylogic_epi64(y[10],x0,x5,0xD2);\
		y[5] = _mm_ternarylogic_epi64(x5,y[3],y[6],0xD2);\
		y[6] = _mm_ternarylogic_epi64(y[6],x2,x0,0xD2);\
\
		y[3] = unpackHI_64(y[9],y[7]);\
		y[9] = alignr(x1,x2,8);\
\
		/*x2    = _mm_slli_si128(y[2],8);\
		  x5    = _mm_slli_si128(x4,8);\
		  y[10] = blend_32(y[0],x2,0xC);\
		  y[7]  = blend_32(x3,x5,0xC);\
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
		x0 = _mm_ternarylogic_epi64(y[0],y[7],y[5],0x96);\
		x0 = _mm_ternarylogic_epi64(y[10],y[2],x0,0x96);\
		\
		x1 = _mm_ternarylogic_epi64(y[1],y[3],y[6],0x96);\
		x1 = _mm_ternarylogic_epi64(y[11],y[8],x1,0x96);\
		\
		x3 = shuffle_32(x2,0x4E);\
		x2 = _mm_ternarylogic_epi64(x3,x2,y[12],0x96);\
\
		/*	The secund part of the step theta\
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
		x4    = blend_32(y[0],y[12],0x3); /*1, 24*/\
		y[10] = xor(y[10],x0);\
		y[11] = ROTV(y[11],r11);/*19, 9*/\
\
		x5   = unpackLO_64(y[4],y[0]);\
		y[4] = unpackHI_64(y[3],y[4]);\
		x0   = ROTV(x4,r0);/*10, 4*/\
					\
		x3 = _mm_ternarylogic_epi64(x5,y[2],x2,0xD2);\
		y[10] = ROTV(y[10],r10); /*24, 14*/\
		x1    = alignr(x0,y[11],8);\
\
		y[9] = ROTV(y[9],r9);/*13, 22*/\
		x4 = _mm_ternarylogic_epi64(y[2],x2,y[8],0xD2);\
		x2 = _mm_ternarylogic_epi64(x2,y[8],x1,0xD2);\
\
		y[5] = unpackLO_64(y[9],y[5]);\
		y[8] = _mm_ternarylogic_epi64(y[8],x1,x5,0xD2);\
		x1 = _mm_ternarylogic_epi64(x1,x5,y[2],0xD2);\
\
		y[1]  = ROTV(y[1],r1); /*5, 20*/\
		y[11] = alignr(y[11],y[10],8);\
\
		y[0] = _mm_ternarylogic_epi64(y[1],y[4],y[5],0xD2);\
		y[7] = ROTV(y[7],r7); /*8, 23*/\
		y[3] = blend_32(y[3],y[6],0xC);\
		y[6] = unpackHI_64(y[6],y[9]);\
					\
		y[2] = _mm_ternarylogic_epi64(y[4],y[5],y[7],0xD2);\
		y[4] = _mm_ternarylogic_epi64(y[11],y[1],y[4],0xD2);\
		y[9] = _mm_ternarylogic_epi64(y[5],y[7],y[11],0xD2);\
		y[7] = _mm_ternarylogic_epi64(y[7],y[11],y[1],0xD2);\
\
		/*I'm not sure if this is the best way to compute these values*/\
		x5    = alignr(y[3],x0,8);\
		x0    = blend_32(x0,y[10],0x3);\
		y[10] = alignr(y[10],y[6],8);\
\
		y[1]  = unpackHI_64(x2,y[8]);\
		y[11] = unpackLO_64(y[9],y[7]);\
		y[8]  = unpackLO_64(x2,y[8]);\
\
		x2   = _mm_ternarylogic_epi64(y[10],x0,x5,0xD2);\
		y[5] = _mm_ternarylogic_epi64(x5,y[3],y[6],0xD2);\
		y[6] = _mm_ternarylogic_epi64(y[6],x2,x0,0xD2);\
\
		y[3] = unpackHI_64(y[9],y[7]);\
		y[9] = alignr(x1,x2,8);\
\
		/*x2    = _mm_slli_si128(y[2],8);\
		  x5    = _mm_slli_si128(x4,8);\
		  y[10] = blend_32(y[0],x2,0xC);\
		  y[7]  = blend_32(x3,x5,0xC);\
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



void print_128(__m128i x){
     uint64_t c[2];
     store(c,0,x);
     printf("\n 1 = %8.16lx\t  0 = %8.16lx\t\n", c[1],c[0]);
}			
			
int keccak(const uint8_t *in, int inlen, uint8_t *md, int r){

  __m128i y[13];
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
	__m128i x0,x1,x2,x3,x4,x5,x6;

  uint8_t *in_temp = (uint8_t*)in;
  int i,q,it=0;
  ALIGN uint8_t temp[144];

  int rsiz = 200 - 2*r; 
  zero_state;
  
  for ( ; inlen >= rsiz; inlen -= rsiz) { 
            loadM(in_temp); 
            keccakF(y, 24); 
            in_temp +=rsiz;
   } 

  // last block and padding 
   memcpy(temp, in_temp, inlen); 
   temp[inlen++] = 0x06;           // XXX Padding Changed from Keccak 3.0 
   memset(temp + inlen, 0, rsiz - inlen); 
   temp[rsiz - 1] |= 0x80; 

   uint8_t *t = temp;
   loadM(t);
   keccakF(y,24);

   if(rsiz==144){
   		store(md,0,y[0]);
    	store(md,1,y[1]);
   }else if(rsiz==136){
   		store(md,0,y[0]);
   		store(md,1,y[1]);	
   }else if(rsiz==104){
   		store(md,0,y[0]);
      	store(md,1,y[1]);
      	store(md,2,unpackLO_64(y[4],y[2]));
   }else if(rsiz == 72){
   		store(md,0,y[0]);
      	store(md,1,y[1]);
      	store(md,2,unpackLO_64(y[4],y[2]));
      	store(md,3,alignr(y[3],y[2],8));
   }
   
 
  return 0;
}

