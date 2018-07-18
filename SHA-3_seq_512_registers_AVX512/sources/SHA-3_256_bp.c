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
	const __m512i p1 = _mm512_set_epi64(8, 4, 10, 1, 12, 3, 9, 0);
	const __m512i p2 = _mm512_set_epi64(10, 1, 12, 3, 11, 2, 9, 0);
	const __m512i p3 = _mm512_set_epi64(9, 8, 3, 2, 11, 10, 1, 0);
	const __m512i p4 = _mm512_set_epi64(15, 14, 7, 6, 13, 12, 5, 4);
	const __m512i p5 = _mm512_set_epi64(2, 1, 0, 12, 3, 2, 1, 0);
	const __m512i p6 = _mm512_set_epi64(6, 5, 4, 10, 7, 6, 5, 4);
	const __m512i p7 = _mm512_set_epi64(2, 1, 0, 8, 3, 2, 1, 0);
	const __m512i p8 = _mm512_set_epi64(6, 5, 4, 11, 7, 6, 5, 4);
	const __m512i p9 = _mm512_set_epi64(4, 3, 2, 9, 0, 4, 3, 2);
	const __m512i p10 = _mm512_set_epi64(1, 1, 1, 0, 4, 3, 2, 1);
	const __m512i p11 = _mm512_set_epi64(4, 0, 1, 2, 3, 0, 0, 0);
	const __m512i p12 = _mm512_set_epi64(1, 2, 3, 4, 0, 0, 0, 0);

	
	uint64_t* rc = RC;
	int i;
	__m512i c1,c2,c3;
/*	print_512(x[0]);
	print_512(x[1]);
	print_512(x[2]);
	print_512(x[3]);
	print_512(x[4]);*/

	for(i=0;i<24;i++){
	/*theta step*/
	c1 = _mm512_xor_si512(x[0],x[1]);
	c2 = _mm512_xor_si512(x[2],x[3]);
	c1 = _mm512_xor_si512(c1,x[4]);
	c1 = _mm512_xor_si512(c1,c2);
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

	/*print_512(x[0]);
	print_512(x[1]);
	print_512(x[2]);
	print_512(x[3]);
	print_512(x[4]);*/

	c1 = _mm512_permutex2var_epi64(x[0],p1,x[1]);/*16-15-11-10-6-5-1-0*/
	c2 = _mm512_permutex2var_epi64(x[2],p2,x[3]);/*18-17-13-12-3-2-8-7*/
	x[0] = _mm512_mask_blend_epi64(0x08,x[0],x[1]);/*x-21-20-x-x*/
	x[1] = _mm512_mask_blend_epi64(0x01,x[2],x[3]);/*22-x-x-x-23*/
	c3 = _mm512_mask_blend_epi64(0x11,x[0],x[1]);/*22-21-20-x-23*/

	x[2] = _mm512_permutex2var_epi64(c1,p3,c2);/*8-7-6-5-3-2-1-0*/
	x[3] = _mm512_permutex2var_epi64(c1,p4,c2);/*18-17-16-15-13-12-11-10*/


	x[0] = _mm512_permutex2var_epi64(x[2],p5,x[4]);/*2-1-0-4-3-2-1-0*/
	x[1] = _mm512_permutex2var_epi64(x[2],p6,x[4]);/*7-6-5-9-8-7-6-5*/
	x[2] = _mm512_permutex2var_epi64(x[3],p7,x[4]);/*12-11-10-14-13-12-11-10*/
	x[3] = _mm512_permutex2var_epi64(x[3],p8,x[4]);/*17-16-15-19-18-17-16-15*/
	x[4] = _mm512_permutex2var_epi64(c3,p9,x[4]);/*22-21-20-24-23-22-21-20*/
	/*rho and pi step*/

	/*chi step*/
	c1 = _mm512_permutexvar_epi64(p10,x[0]);
	c2 = _mm512_permutexvar_epi64(p10,c1);
	x[0] = _mm512_xor_si512(x[0],_mm512_andnot_si512(c1,c2));

	c1 = _mm512_permutexvar_epi64(p10,x[1]);
	c2 = _mm512_permutexvar_epi64(p10,c1);
	x[1] = _mm512_xor_si512(x[1],_mm512_andnot_si512(c1,c2));

	c1 = _mm512_permutexvar_epi64(p10,x[2]);
	c2 = _mm512_permutexvar_epi64(p10,c1);
	x[2] = _mm512_xor_si512(x[2],_mm512_andnot_si512(c1,c2));

	c1 = _mm512_permutexvar_epi64(p10,x[3]);
	c2 = _mm512_permutexvar_epi64(p10,c1);
	x[3] = _mm512_xor_si512(x[3],_mm512_andnot_si512(c1,c2));

	c1 = _mm512_permutexvar_epi64(p10,x[4]);
	c2 = _mm512_permutexvar_epi64(p10,c1);
	x[4] = _mm512_xor_si512(x[4],_mm512_andnot_si512(c1,c2));
	/*chi step*/

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
	
    uint8_t * in_temp = in;
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


*   // last block and padding
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

/*
int main(int argc, char **argv)
{
   ALIGN char *msgstr;
   int inlen,i;
   ALIGN uint8_t md[64];
   ALIGN uint8_t md1[64];
      

   if(argc != 2){
     printf("ERROR!, You need to pass the message size\n\n");
     exit(1);
   }

   inlen = atoi(argv[1]);

   msgstr = (char*)_mm_malloc(inlen*sizeof(char),32);
   
   for(i=0;i<inlen;i++){
      msgstr[i] = (char)(rand()%256);
    }
   

    printf("\nSequential implementation using 128-bit registers.\n\n");    
    printf("<------------------------------------------------------>\n");    
      
   printf("Keccak_256 \t\t\t");
    
   memset(md, 0, 64*sizeof(uint8_t));
   memset(md1, 0, 64*sizeof(uint8_t));
   
   keccak((uint8_t*)msgstr, inlen,md1,32);
   keccak_std(msgstr, inlen,md,32);
   
   if(memcmp(md,md1,32) == 0){
     printf("ok!\n");
    }else{
      printf("Error!!\n");
    }
    
    printf("<------------------------------------------------------>\n\n\n");    
    

      
  _mm_free(msgstr);

  return 0;

}
*/



/*
int main(){
 ALIGN uint64_t RC[60] = {
    0x0,0x1, 0x2,0x3,
    0x4,0x5, 0x6,0x7,
    0x8,0x9, 0x10,0x11,
    0x12,0x13, 0x14,0x15,
    0x16,0x17, 0x18,0x19,
    0x20,0x21, 0x22,0x23,
    0x24,0x25, 0x26,0x27,
    0x28,0x29, 0x30,0x31,
    0x32,0x33, 0x34,0x35,
    0x36,0x37, 0x38,0x39,
    0x40,0x41, 0x42,0x43,
    0x44,0x45, 0x46,0x47
  };

	uint8_t * temp = (uint8_t*)RC;
__m512i x[5];
x[0] = _mm512_setzero_si512();
x[1] = _mm512_setzero_si512();
x[2] = _mm512_setzero_si512();
x[3] = _mm512_setzero_si512();
x[4] = _mm512_setzero_si512();



	__m512i x[5];  
	const uint8_t *in_temp = in;

	ALIGN uint8_t temp[144];
       

  	memset(temp, 0, 144*sizeof(uint8_t));

  	x[4] = _mm512_setzero_si512();
	for ( ; inlen >= r; inlen -= r) {
	    load(in_temp);
	    getchar();
        keccakF(x, 24);
    }

    // last block and padding
    memcpy(temp, in, inlen);
    temp[inlen++] = 0x06;           // XXX Padding Changed from Keccak 3.0
    memset(temp + inlen, 0, r - inlen);
    temp[r - 1] |= 0x80;

    load(in_temp);
	keccakF(x, 24);

    //memcpy(md, st, mdlen);
    return 0;
}
*/

