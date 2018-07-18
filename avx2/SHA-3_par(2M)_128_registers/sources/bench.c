#include <stdio.h>
#include "bench.h"
#include "keccak.h"

/*====================================================================*/
/* Private definitions                                                */
/*====================================================================*/

/**
 * Timer type.
 */

typedef unsigned long long bench_t;

static inline bench_t cycles(void) {
	unsigned int hi, lo;
	__asm volatile ("rdtsc\n\t":"=a" (lo), "=d"(hi));
	return ((bench_t) lo) | (((bench_t) hi) << 32);
}

/**
 * Stores the time measured before the execution of the benchmark.
 */
static bench_t before;

/**
 * Stores the time measured after the execution of the benchmark.
 */
static bench_t after;

/**
 * Stores the sum of timings for the current benchmark.
 */
long long total;

/*====================================================================*/
/* Public definitions                                                 */
/*====================================================================*/

void bench_reset() {
	total = 0;
}

void bench_before() {
	before = cycles();
}

void bench_after() {
	long long result;
	after = cycles();
	result = (after - before);
	total += result;
}

void bench_compute(int benches) {
	total = total / benches;
}

void bench_print() {
	printf("%lld cycles\n", total);
	printf("\n");
}

unsigned long long bench_get_total() {
	return total;
}



int main(int argc, char **argv)
{
   int inlen,i,j;
   ALIGN char *msgstr[2];
   ALIGN uint8_t *md[2];    

  
    if(argc != 2){
      printf("ERROR!, You need to pass the message size\n\n");
      exit(1);
    }

    inlen = atoi(argv[1]);
   
   for(i=0;i<2;i++){
      msgstr[i] = (char*)_mm_malloc(inlen*sizeof(char),32);
      md[i]  = (uint8_t*)_mm_malloc(64*sizeof(uint8_t),32);
   }    
    
    printf("\n 4-way implementation using 256-bit registers.\n\n");    
    printf("<------------------------------------------------------>\n");    
    
    printf("\nSHA3-256.\n");    
    BENCH_FUNCTION(keccak,msgstr, inlen,md,136);
    printf("cycles per bytes :%f\n",(double)total/inlen);
    
    printf("\nSHA3-384.\n");    
    BENCH_FUNCTION(keccak,msgstr, inlen,md,104);
    printf("cycles per bytes :%f\n",(double)total/inlen);
    
    printf("\nSHA3-512.\n");    
    BENCH_FUNCTION(keccak,msgstr, inlen,md,72);
    printf("cycles per bytes :%f\n",(double)total/inlen);
    
    printf("<------------------------------------------------------>\n\n\n");    

   for(i=0;i<2;i++){
	  _mm_free(msgstr[i]);
	  _mm_free(md[i]);
    }     
   return 0;
   
}

