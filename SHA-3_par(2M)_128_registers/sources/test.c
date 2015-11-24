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
    }
    
    printf("\n 4-way implementation using 256-bit registers.\n\n");    
    printf("<------------------------------------------------------>\n");    
    
   printf("Keccak_256 4 way \t\t\t");

   keccak(msgstr, inlen,md,136);

   keccak_std(msgstr[0], inlen,md1[0],32);
   keccak_std(msgstr[1], inlen,md1[1],32);
   if((memcmp(md[0],md1[0],32) && memcmp(md[1],md1[1],32)!= 0)){
     printf("Error!!\n");
     exit(1);
   }
    
    printf("Ok!!\n");
    
    printf("Keccak_384 4 way \t\t\t");
    
    keccak(msgstr, inlen,md,104);
   
    keccak_std(msgstr[0], inlen,md1[0],48);
    keccak_std(msgstr[1], inlen,md1[1],48);
    
   if((memcmp(md[0],md1[0],48) && memcmp(md[1],md1[1],48)!= 0)){
    printf("Error!!\n");
    exit(1);
   }
    printf("Ok!!\n");
    
    printf("Keccak_512 4 way \t\t\t");

    keccak(msgstr, inlen,md,72);
   
    keccak_std(msgstr[0], inlen,md1[0],64);
    keccak_std(msgstr[1], inlen,md1[1],64);
    
   if((memcmp(md[0],md1[0],64) && memcmp(md[1],md1[1],64)!= 0)){
    printf("Error!!\n");
    exit(1);
   }
    
    printf("Ok!!\n");
    
    printf("<------------------------------------------------------>\n");    

   for(i=0;i<2;i++){
	  _mm_free(msgstr[i]);
	  _mm_free(md[i]);
	  _mm_free(md1[i]);
    } 
 
 return 0;
}

