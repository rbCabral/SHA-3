#include "immintrin.h"
#include "emmintrin.h"
#include <string.h>
#include <stdio.h>
#include <stdint.h>

#define BENCH 100
#define ALIGN __attribute__ ((aligned (32)))

#define KECCAK256

typedef unsigned long long uInt256[4];
typedef unsigned long long uInt512[8];

// compute a keccak hash (md) of given byte length from "in"
int keccak_std(const uint8_t *in, int inlen, uint8_t *md, int mdlen);

int keccak_256(char **in_e, int inlen, uint8_t *md,uint8_t *md1,uint8_t *md2,uint8_t *md3);
int keccak_384(char **in_e, int inlen, uint8_t *md,uint8_t *md1,uint8_t *md2,uint8_t *md3);
int keccak_512	(char **in_e, int inlen, uint8_t *md,uint8_t *md1,uint8_t *md2,uint8_t *md3);





