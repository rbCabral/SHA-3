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
      md[i]  = (uint8_t*)_mm_malloc(200*sizeof(uint8_t),32);
      md1[i]  = (uint8_t*)_mm_malloc(200*sizeof(uint8_t),32);
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
    
   printf("SHA3-224 4 way \t\t\t");

   keccak(msgstr, inlen,md,28,28);

    keccak_std((uint8_t*)msgstr[0], inlen,md1[0],28);
    keccak_std((uint8_t*)msgstr[1], inlen,md1[1],28);
    keccak_std((uint8_t*)msgstr[2], inlen,md1[2],28);
    keccak_std((uint8_t*)msgstr[3], inlen,md1[3],28);
    keccak_std((uint8_t*)msgstr[4], inlen,md1[4],28);
    keccak_std((uint8_t*)msgstr[5], inlen,md1[5],28);
    keccak_std((uint8_t*)msgstr[6], inlen,md1[6],28);
    keccak_std((uint8_t*)msgstr[7], inlen,md1[7],28);
    
    if(memcmp(md[0],md1[0],28)!=0 || memcmp(md[1],md1[1],28) !=0 || memcmp(md[2],md1[2],28)!=0 || memcmp(md[3],md1[3],28) != 0 \
		    || memcmp(md[4],md1[4],28)!=0 || memcmp(md[5],md1[5],28) !=0 || memcmp(md[6],md1[6],28)!=0 || memcmp(md[7],md1[7],28) != 0)
	    printf("Error!!\n");
    else
	    printf("Ok!!\n");
    
    printf("SHA3-256 4 way \t\t\t");

   keccak(msgstr, inlen,md,32,32);

    keccak_std((uint8_t*)msgstr[0], inlen,md1[0],32);
    keccak_std((uint8_t*)msgstr[1], inlen,md1[1],32);
    keccak_std((uint8_t*)msgstr[2], inlen,md1[2],32);
    keccak_std((uint8_t*)msgstr[3], inlen,md1[3],32);
    keccak_std((uint8_t*)msgstr[4], inlen,md1[4],32);
    keccak_std((uint8_t*)msgstr[5], inlen,md1[5],32);
    keccak_std((uint8_t*)msgstr[6], inlen,md1[6],32);
    keccak_std((uint8_t*)msgstr[7], inlen,md1[7],32);
    
    if(memcmp(md[0],md1[0],32)!=0 || memcmp(md[1],md1[1],32) !=0 || memcmp(md[2],md1[2],32)!=0 || memcmp(md[3],md1[3],32) != 0 \
		    || memcmp(md[4],md1[4],32)!=0 || memcmp(md[5],md1[5],32) !=0 || memcmp(md[6],md1[6],32)!=0 || memcmp(md[7],md1[7],32) != 0)
	    printf("Error!!\n");
    else
	    printf("Ok!!\n");

    printf("SHA3-384 4 way \t\t\t");

   keccak(msgstr, inlen,md,48,48);

    keccak_std((uint8_t*)msgstr[0], inlen,md1[0],48);
    keccak_std((uint8_t*)msgstr[1], inlen,md1[1],48);
    keccak_std((uint8_t*)msgstr[2], inlen,md1[2],48);
    keccak_std((uint8_t*)msgstr[3], inlen,md1[3],48);
    keccak_std((uint8_t*)msgstr[4], inlen,md1[4],48);
    keccak_std((uint8_t*)msgstr[5], inlen,md1[5],48);
    keccak_std((uint8_t*)msgstr[6], inlen,md1[6],48);
    keccak_std((uint8_t*)msgstr[7], inlen,md1[7],48);
    
    if(memcmp(md[0],md1[0],48)!=0 || memcmp(md[1],md1[1],48) !=0 || memcmp(md[2],md1[2],48)!=0 || memcmp(md[3],md1[3],48) != 0 \
		    || memcmp(md[4],md1[4],48)!=0 || memcmp(md[5],md1[5],48) !=0 || memcmp(md[6],md1[6],48)!=0 || memcmp(md[7],md1[7],48) != 0)
	    printf("Error!!\n");
    else
	    printf("Ok!!\n");

   printf("Keccak_512 4 way \t\t\t");

   keccak(msgstr, inlen,md,64,64);

    keccak_std((uint8_t*)msgstr[0], inlen,md1[0],64);
    keccak_std((uint8_t*)msgstr[1], inlen,md1[1],64);
    keccak_std((uint8_t*)msgstr[2], inlen,md1[2],64);
    keccak_std((uint8_t*)msgstr[3], inlen,md1[3],64);
    keccak_std((uint8_t*)msgstr[4], inlen,md1[4],64);
    keccak_std((uint8_t*)msgstr[5], inlen,md1[5],64);
    keccak_std((uint8_t*)msgstr[6], inlen,md1[6],64);
    keccak_std((uint8_t*)msgstr[7], inlen,md1[7],64);
    
    if(memcmp(md[0],md1[0],64)!=0 || memcmp(md[1],md1[1],64) !=0 || memcmp(md[2],md1[2],64)!=0 || memcmp(md[3],md1[3],64) != 0 \
		    || memcmp(md[4],md1[4],64)!=0 || memcmp(md[5],md1[5],64) !=0 || memcmp(md[6],md1[6],64)!=0 || memcmp(md[7],md1[7],64) != 0)
	    printf("Error!!\n");
    else
	    printf("Ok!!\n");

    printf("SHAKE128 4 way \t\t\t");

   keccak(msgstr, inlen,md,16,16);

    keccak_std((uint8_t*)msgstr[0], inlen,md1[0],16);
    keccak_std((uint8_t*)msgstr[1], inlen,md1[1],16);
    keccak_std((uint8_t*)msgstr[2], inlen,md1[2],16);
    keccak_std((uint8_t*)msgstr[3], inlen,md1[3],16);
    keccak_std((uint8_t*)msgstr[4], inlen,md1[4],16);
    keccak_std((uint8_t*)msgstr[5], inlen,md1[5],16);
    keccak_std((uint8_t*)msgstr[6], inlen,md1[6],16);
    keccak_std((uint8_t*)msgstr[7], inlen,md1[7],16);
    
    if(memcmp(md[0],md1[0],16)!=0 || memcmp(md[1],md1[1],16) !=0 || memcmp(md[2],md1[2],16)!=0 || memcmp(md[3],md1[3],16) != 0 \
		    || memcmp(md[4],md1[4],16)!=0 || memcmp(md[5],md1[5],16) !=0 || memcmp(md[6],md1[6],16)!=0 || memcmp(md[7],md1[7],16) != 0)
	    printf("Error!!\n");
    else
	    printf("Ok!!\n");

    printf("SHAKE256 4 way \t\t\t");

   keccak(msgstr, inlen,md,32,32);

    keccak_std((uint8_t*)msgstr[0], inlen,md1[0],32);
    keccak_std((uint8_t*)msgstr[1], inlen,md1[1],32);
    keccak_std((uint8_t*)msgstr[2], inlen,md1[2],32);
    keccak_std((uint8_t*)msgstr[3], inlen,md1[3],32);
    keccak_std((uint8_t*)msgstr[4], inlen,md1[4],32);
    keccak_std((uint8_t*)msgstr[5], inlen,md1[5],32);
    keccak_std((uint8_t*)msgstr[6], inlen,md1[6],32);
    keccak_std((uint8_t*)msgstr[7], inlen,md1[7],32);
    
    if(memcmp(md[0],md1[0],32)!=0 || memcmp(md[1],md1[1],32) !=0 || memcmp(md[2],md1[2],32)!=0 || memcmp(md[3],md1[3],32) != 0 \
		    || memcmp(md[4],md1[4],32)!=0 || memcmp(md[5],md1[5],32) !=0 || memcmp(md[6],md1[6],32)!=0 || memcmp(md[7],md1[7],32) != 0)
	    printf("Error!!\n");
    else
	    printf("Ok!!\n");

       
    printf("<------------------------------------------------------>\n");    

   for(i=0;i<8;i++){
	  _mm_free(msgstr[i]);
	  _mm_free(md[i]);
	  _mm_free(md1[i]);
    } 
 
 return 0;
}

