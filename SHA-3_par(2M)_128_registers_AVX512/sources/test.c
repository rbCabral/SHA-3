#include "keccak.h"
#include <stdint.h>

int main(int argc, char **argv)
{
  int inlen,i,j; 
  ALIGN char *msgstr[2];
  ALIGN uint8_t *md[2];    
  ALIGN uint8_t *md1[64];


   if(argc != 2){
     printf("ERROR!, You need to pass the message size\n\n");
     exit(1);
   }

   inlen = atoi(argv[1]);

   for(i=0;i<2;i++){
      msgstr[i] = (char*)_mm_malloc(inlen*sizeof(char),32);
      md[i]  = (uint8_t*)_mm_malloc(64*sizeof(uint8_t),32);
      md1[i]  = (uint8_t*)_mm_malloc(64*sizeof(uint8_t),32);
    }
    
    for(i=0;i<inlen;i++){
      msgstr[0][i] = (char)(rand()%256);
      msgstr[1][i] = (char)(rand()%256);
    }
    
    printf("\n 2-way implementation using 128-bit registers.\n\n");    
    printf("<------------------------------------------------------>\n");    
   
  printf("SHA3-224 2 2-way \t\t\t");

   keccak(msgstr, inlen,md,28,28);

   keccak_std((uint8_t*)msgstr[0], inlen,md1[0],28);
   keccak_std((uint8_t*)msgstr[1], inlen,md1[1],28);

    if((memcmp(md[0],md1[0],28) !=0 || memcmp(md[1],md1[1],28)!= 0))
       printf("Error!!\n");
    else
       printf("Ok!!\n");

   printf("SHA3-256 2 way \t\t\t");

   keccak(msgstr, inlen,md,32,32);

   keccak_std((uint8_t*)msgstr[0], inlen,md1[0],32);
   keccak_std((uint8_t*)msgstr[1], inlen,md1[1],32);

    if((memcmp(md[0],md1[0],32) !=0 || memcmp(md[1],md1[1],32)!= 0))
       printf("Error!!\n");
    else
       printf("Ok!!\n");
    
    printf("SHA3-384 2 way \t\t\t");
    
    keccak(msgstr, inlen,md,48,48);
   
    keccak_std((uint8_t*)msgstr[0], inlen,md1[0],48);
    keccak_std((uint8_t*)msgstr[1], inlen,md1[1],48);
    
    if((memcmp(md[0],md1[0],48) !=0 || memcmp(md[1],md1[1],48)!= 0))
      printf("Error!!\n");
    else
      printf("Ok!!\n");
    
    printf("SHA3-512 2 way \t\t\t");

    keccak(msgstr, inlen,md,64,64);
   
    keccak_std((uint8_t*)msgstr[0], inlen,md1[0],64);
    keccak_std((uint8_t*)msgstr[1], inlen,md1[1],64);
    
    if((memcmp(md[0],md1[0],64) != 0 || memcmp(md[1],md1[1],64)!= 0))
      printf("Error!!\n");
    else    
      printf("Ok!!\n");
    
    printf("SHAKE128 2 way \t\t\t");

    keccak(msgstr, inlen,md,16,16);
   
    keccak_std((uint8_t*)msgstr[0], inlen,md1[0],16);
    keccak_std((uint8_t*)msgstr[1], inlen,md1[1],16);
    
    if((memcmp(md[0],md1[0],16) != 0 || memcmp(md[1],md1[1],16)!= 0))
      printf("Error!!\n");
    else    
      printf("Ok!!\n");

    printf("SHAKE256 2 way \t\t\t");

    keccak(msgstr, inlen,md,32,32);
   
    keccak_std((uint8_t*)msgstr[0], inlen,md1[0],32);
    keccak_std((uint8_t*)msgstr[1], inlen,md1[1],32);
    
    if((memcmp(md[0],md1[0],32) != 0 || memcmp(md[1],md1[1],32)!= 0))
      printf("Error!!\n");
    else    
      printf("Ok!!\n");

    printf("<------------------------------------------------------>\n");    

   for(i=0;i<2;i++){
	  _mm_free(msgstr[i]);
	  _mm_free(md[i]);
	  _mm_free(md1[i]);
    } 
 
 return 0;
}

