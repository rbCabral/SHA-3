//#include "immintrin.h"
//#include "emmintrin.h"
//#include <string.h>
//#include <stdio.h>
//#include <stdint.h>
#include "keccak.h"

/*#ifdef __INTEL_COMPILER
#define ALIGN __declspec(align(32))
#else
#define ALIGN __attribute__ ((aligned (64)))
#endif
*/
#define load(M)	\
	x[0] = _mm512_xor_si512(x[0],_mm512_maskz_loadu_epi64(0x1F,(__m512i*)M)+0);\
	M=M + 40;\
	x[1] = _mm512_xor_si512(x[1],_mm512_maskz_loadu_epi64(0x1F,(__m512i*)M)+0);	\
	M+=40;\
	x[2] = _mm512_xor_si512(x[2],_mm512_maskz_loadu_epi64(0x1F,(__m512i*)M)+0);	\
	M+=40;\
	x[3] = _mm512_xor_si512(x[3],_mm512_maskz_loadu_epi64(0x03,(__m512i*)M)+0);	\
	M+=16;\
	

void print_512(__m512i x){
    uint64_t c[8];
    _mm512_store_si512((__m512*)c,x);
    printf("\n 0 = %8.16lx\t  1 = %8.16lx\t  2 = %8.16lx\t  3 = %8.16lx\t\n 4 = %8.16lx\t  5 = %8.16lx\t  6 = %8.16lx\t  7 = %8.16lx\t\n", c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7]);
}


void keccakF(__m512i* x, int rnd){
	const uint64_t RC[24] = {
	    0x0000000000000001, 0x0000000000008082, 0x800000000000808a,
	    0x8000000080008000, 0x000000000000808b, 0x0000000080000001,
	    0x8000000080008081, 0x8000000000008009, 0x000000000000008a,
	    0x0000000000000088, 0x0000000080008009, 0x000000008000000a,
	    0x000000008000808b, 0x800000000000008b, 0x8000000000008089,
	    0x8000000000008003, 0x8000000000008002, 0x8000000000000080,
	    0x000000000000800a, 0x800000008000000a, 0x8000000080008081,
	    0x8000000000008080, 0x0000000080000001, 0x8000000080008008
	};
	const __m512i r0 = _mm512_set_epi64(0, 0, 0, 27, 28, 62, 1, 0);
	const __m512i r1 = _mm512_set_epi64(0, 0, 0, 20, 55, 6, 44, 36);
	const __m512i r2 = _mm512_set_epi64(0, 0, 0, 39, 25, 43, 10, 3);
	const __m512i r3 = _mm512_set_epi64(0, 0, 0, 8, 21, 15, 45, 41);
	const __m512i r4 = _mm512_set_epi64(0, 0, 0, 14, 56, 61, 2, 18);
	const __m512i p1 = _mm512_set_epi64(0, 0, 0, 2, 1, 0, 4, 3);
	const __m512i p2 = _mm512_set_epi64(0, 0, 0, 3, 2, 1, 0, 4);
	const __m512i p9 = _mm512_set_epi64(0, 0, 0, 1, 0, 4, 3, 2);
	const __m512i p10 = _mm512_set_epi64(0, 0, 0, 0, 4, 3, 2, 1);
	const __m512i p11 = _mm512_set_epi64(0, 0, 0, 3, 2, 1, 0, 4);
	const __m512i p12 = _mm512_set_epi64(0, 0, 0, 0, 4, 3, 2, 1);

	
	uint64_t* rc = (uint64_t*)RC;
	int i;
	__m512i c1,c2,c3,c4;
/*	print_512(x[0]);
	print_512(x[1]);
	print_512(x[2]);
	print_512(x[3]);
	print_512(x[4]);*/

	for(i=0;i<24;i++){

	/*theta step*/
	//c1 = _mm512_xor_si512(x[0],x[1]);
	//c2 = _mm512_xor_si512(x[2],x[3]);
	//c1 = _mm512_xor_si512(c1,x[4]);
	c1 = _mm512_ternarylogic_epi64(x[0],x[1],x[2],0x96);
	c1 = _mm512_ternarylogic_epi64(c1,x[3],x[4],0x96);

	//c1 = _mm512_xor_si512(c1,c2);
	c2 = _mm512_permutexvar_epi64(p11, c1);
        c1 = _mm512_permutexvar_epi64(p12, c1);
        c1 = _mm512_xor_si512(c2,_mm512_rol_epi64(c1,0x01));

	x[0] = _mm512_xor_si512(c1,x[0]);
	x[1] = _mm512_xor_si512(c1,x[1]);
	x[2] = _mm512_xor_si512(c1,x[2]);
	x[3] = _mm512_xor_si512(c1,x[3]);
	x[4] = _mm512_xor_si512(c1,x[4]);
	/*theta step*/


	/*rho and pi step*/
	x[0] = _mm512_rolv_epi64(x[0],r0);/*15-5-20-10-0*/
	x[1] = _mm512_rolv_epi64(x[1],r1);/*6-21-11-1-16*/
	x[2] = _mm512_rolv_epi64(x[2],r2);/*22-12-2-17-7*/
	x[3] = _mm512_rolv_epi64(x[3],r3);/*13-3-18-8-23*/
	x[4] = _mm512_rolv_epi64(x[4],r4);/*4-19-9-24-14*/


	c1   = _mm512_mask_blend_epi64(0x12,x[1],x[0]);/*15-x-11-10-16*/
	c2   = _mm512_mask_blend_epi64(0x08,x[0],x[1]);/*x-21-22-x-x*/
	x[0] = _mm512_mask_blend_epi64(0x12,x[0],x[1]);/*6-5-x-1-0*/
	x[1] = _mm512_mask_blend_epi64(0x0A,x[2],x[3]);/*x-3-2-8-7*/
	c3   = _mm512_mask_blend_epi64(0x1A,x[3],x[2]);/*22-x-18-17-23*/
	c4   = _mm512_mask_blend_epi64(0x10,x[2],x[3]);/*13-12-x-x-x*/
	
	c2   = _mm512_mask_blend_epi64(0x11,c2,c3);/*22-21-20-x-23*/
	c3   = _mm512_mask_blend_epi64(0x11,c3,c1);/*15-x-18-17-16*/
	c1   = _mm512_mask_blend_epi64(0x18,c1,c4);/*13-12-11-10-x*/
	c4   = _mm512_mask_blend_epi64(0x03,x[0],x[1]);/*6-5-x-8-7*/
	x[0] = _mm512_mask_blend_epi64(0x0C,x[0],x[1]);/*x-3-2-1-0*/

	x[0] = _mm512_mask_blend_epi64(0x10,x[0],x[4]);/*4-3-2-1-0*/
	x[1] = _mm512_mask_blend_epi64(0x04,c4,x[4]);/*6-5-9-8-7*/
	x[2] = _mm512_mask_blend_epi64(0x01,c1,x[4]);/*13-12-11-10-14*/
	x[3] = _mm512_mask_blend_epi64(0x08,c3,x[4]);/*15-19-18-17-16*/
	x[4] = _mm512_mask_blend_epi64(0x02,c2,x[4]);/*22-21-20-24-23*/
	/*rho and pi step*/
	

	/*chi step*/
	c1 = _mm512_permutexvar_epi64(p10,x[0]);
	c2 = _mm512_permutexvar_epi64(p10,c1);
	//x[0] = _mm512_xor_si512(x[0],_mm512_andnot_si512(c1,c2));
	x[0] = _mm512_ternarylogic_epi64(x[0],c1,c2,0xD2);

	c1 = _mm512_permutexvar_epi64(p1,x[1]);
	c2 = _mm512_permutexvar_epi64(p2,x[1]);
	//x[1] = _mm512_xor_si512(c1,_mm512_andnot_si512(c2,x[1]));
	x[1] = _mm512_ternarylogic_epi64(c1,c2,x[1],0xD2);

	x[2] = _mm512_permutexvar_epi64(p10,x[2]);
	c1 = _mm512_permutexvar_epi64(p10,x[2]);
	c2 = _mm512_permutexvar_epi64(p10,c1);
	//x[2] = _mm512_xor_si512(x[2],_mm512_andnot_si512(c1,c2));
	x[2] = _mm512_ternarylogic_epi64(x[2],c1,c2,0xD2);

	c1 = _mm512_permutexvar_epi64(p2,x[3]);
	c2 = _mm512_permutexvar_epi64(p10,x[3]);
	//x[3] = _mm512_xor_si512(c1,_mm512_andnot_si512(x[3],c2));
	x[3] = _mm512_ternarylogic_epi64(c1,x[3],c2,0xD2);

	x[4] = _mm512_permutexvar_epi64(p9,x[4]);
	c1 = _mm512_permutexvar_epi64(p10,x[4]);
	c2 = _mm512_permutexvar_epi64(p10,c1);
	//x[4] = _mm512_xor_si512(x[4],_mm512_andnot_si512(c1,c2));
	x[4] = _mm512_ternarylogic_epi64(x[4],c1,c2,0xD2);
	/*chi step*/
/*	print_512(x[0]);
	print_512(x[1]);
	print_512(x[2]);
	print_512(x[3]);
	print_512(x[4]);
	getchar();*/

	/*iota step*/
	c1 = _mm512_maskz_loadu_epi64(0x01,(__m512i*)rc);
	rc+=1;
	x[0] = _mm512_xor_si512(x[0],c1);
	/*iota step*/
	}
/*	print_512(x[0]);
	print_512(x[1]);
	print_512(x[2]);
	print_512(x[3]);
	print_512(x[4]);*/
}

int keccak(const uint8_t *in, int inlen, uint8_t *md, int r){
	
    uint8_t * in_temp = (uint8_t*)in;
    uint8_t *t;
    __m512i x[5];

    ALIGN uint8_t temp[144];
    memset(temp, 0, 144*sizeof(uint8_t));
    x[0] = _mm512_setzero_si512();
    x[1] = _mm512_setzero_si512();
    x[2] = _mm512_setzero_si512();
    x[3] = _mm512_setzero_si512();
    x[4] = _mm512_setzero_si512();
   int rsiz = 200 - 2*32;
   for ( ; inlen >= rsiz; inlen -= rsiz) {
            load(in_temp);
        	keccakF(x, 24);
//	exit(1);
   }
/*   	print_512(x[0]);
	print_512(x[1]);
	print_512(x[2]);
	print_512(x[3]);
	print_512(x[4]);*/


   // last block and padding
    memcpy(temp, in_temp, inlen);
    temp[inlen++] = 0x06;           // XXX Padding Changed from Keccak 3.0
    memset(temp + inlen, 0, rsiz - inlen);
    temp[rsiz - 1] |= 0x80;

    t = temp;
    load(t);
    keccakF(x, 24);
   
    _mm512_storeu_si512((__m512*)md,x[0]);
    
    return 0;
}


