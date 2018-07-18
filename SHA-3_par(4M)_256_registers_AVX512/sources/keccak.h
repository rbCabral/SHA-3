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

int keccak(char **in_e, int inlen, uint8_t **md, int rsiz);






