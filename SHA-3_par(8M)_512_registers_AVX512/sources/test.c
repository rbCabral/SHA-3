#include "keccak.h"
#include <stdint.h>

int main(int argc, char **argv)
{
  int inlen,i,j; 
  ALIGN char *msgstr[8];
  ALIGN uint8_t *md[8];    
  ALIGN uint8_t *md1[8];


   if(argc != 2){
     printf("ERROR!, You need to pass the message size\n\n");
     exit(1);
   }

   inlen = atoi(argv[1]);

   for(i=0;i<8;i++){
      msgstr[i] = (char*)_mm_malloc(inlen*sizeof(char),32);
      md[i]  = (uint8_t*)_mm_malloc(65*sizeof(uint8_t),32);
      md1[i]  = (uint8_t*)_mm_malloc(65*sizeof(uint8_t),32);
    }
    
    for(j=0;j<8;j++){
      for(i=0;i<inlen;i++){
	       msgstr[j][i] = (char)(rand()%256);
/*	      char c =  (char)(rand()%256);

	      msgstr[0][i] = c;
	      msgstr[1][i] = c;
	      msgstr[2][i] = c;
	      msgstr[3][i] = c;
	      msgstr[4][i] = c;
	      msgstr[5][i] = c;
	      msgstr[6][i] = c;
	      msgstr[7][i] = c;*/
      }
    }
    
    printf("\n 4-way implementation using 256-bit registers.\n\n");    
    printf("<------------------------------------------------------>\n");    
    
   printf("Keccak_256 4 way \t\t\t");

   keccak(msgstr, inlen,md,32);

    keccak_std((uint8_t*)msgstr[0], inlen,md1[0],32);
    keccak_std((uint8_t*)msgstr[1], inlen,md1[1],32);
    keccak_std((uint8_t*)msgstr[2], inlen,md1[2],32);
    keccak_std((uint8_t*)msgstr[3], inlen,md1[3],32);
    keccak_std((uint8_t*)msgstr[4], inlen,md1[4],32);
    keccak_std((uint8_t*)msgstr[5], inlen,md1[5],32);
    keccak_std((uint8_t*)msgstr[6], inlen,md1[6],32);
    keccak_std((uint8_t*)msgstr[7], inlen,md1[7],32);
    
    if(memcmp(md[0],md1[0],32)!=0 || memcmp(md[1],md1[1],32) !=0 || memcmp(md[2],md1[2],32)!=0 || memcmp(md[3],md1[3],32) != 0 \
		    || memcmp(md[4],md1[4],32)!=0 || memcmp(md[5],md1[5],32) !=0 || memcmp(md[6],md1[6],32)!=0 || memcmp(md[7],md1[7],32) != 0){
      printf("Error!!\n");
      exit(1);
    }
    
    printf("Ok!!\n");
   /*
    printf("Keccak_384 4 way \t\t\t");
    
    keccak(msgstr, inlen,md,104);
   
    keccak_std((uint8_t*)msgstr[0], inlen,md1[0],48);
    keccak_std((uint8_t*)msgstr[1], inlen,md1[1],48);
    keccak_std((uint8_t*)msgstr[2], inlen,md1[2],48);
    keccak_std((uint8_t*)msgstr[3], inlen,md1[3],48);

    if((memcmp(md[0],md1[0],48) || memcmp(md[1],md1[1],48) || memcmp(md[2],md1[2],48) || memcmp(md[3],md1[3],48))  != 0){
      printf("Error!!\n");
      exit(1);
    }
    printf("Ok!!\n");
    
    printf("Keccak_512 4 way \t\t\t");

    keccak(msgstr, inlen,md,72);
   
    keccak_std((uint8_t*)msgstr[0], inlen,md1[0],64);
    keccak_std((uint8_t*)msgstr[1], inlen,md1[1],64);
    keccak_std((uint8_t*)msgstr[2], inlen,md1[2],64);
    keccak_std((uint8_t*)msgstr[3], inlen,md1[3],64);

    if((memcmp(md[0],md1[0],64) || memcmp(md[1],md1[1],64) || memcmp(md[2],md1[2],64) || memcmp(md[3],md1[3],64))  != 0){
      printf("Error!!\n");
      exit(1);
    }
    
    printf("Ok!!\n");
    */
    printf("<------------------------------------------------------>\n");    

   for(i=0;i<8;i++){
	  _mm_free(msgstr[i]);
	  _mm_free(md[i]);
	  _mm_free(md1[i]);
    } 
 
 return 0;
}

