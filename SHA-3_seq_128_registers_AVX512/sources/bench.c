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
   ALIGN uint8_t *msgstr;
   int inlen,i,j;
   ALIGN uint8_t md[64];

  
    if(argc != 2){
      printf("ERROR!, You need to pass the message size\n\n");
      exit(1);
    }

    inlen = atoi(argv[1]);
   
    msgstr = (uint8_t*)_mm_malloc(inlen*sizeof(uint8_t),32);
    
    for(i=0;i<inlen;i++){
      msgstr[i] = (char)(rand()%256);
    }
    
    printf("\nSequential implementation using 128-bit registers.\n\n");    
    printf("<------------------------------------------------------>\n");    
    
    printf("\nSHA3-224.\n");    
    BENCH_FUNCTION(keccak,msgstr, inlen,md,28);
    printf("cycles per bytes :%f\n",(double)total/inlen);

    printf("\nSHA3-256.\n");    
    BENCH_FUNCTION(keccak,msgstr, inlen,md,32);
    printf("cycles per bytes :%f\n",(double)total/inlen);
    
    printf("\nSHA3-384.\n");    
    BENCH_FUNCTION(keccak,msgstr, inlen,md,48);
    printf("cycles per bytes :%f\n",(double)total/inlen);
    
  /*  printf("\nSHA3-512.\n");    
    BENCH_FUNCTION(keccak,msgstr, inlen,md,64);
    printf("cycles per bytes :%f\n",(double)total/inlen);*/
    
    printf("<------------------------------------------------------>\n\n\n");    


    
    _mm_free(msgstr);
    
   return 0;
   
}

