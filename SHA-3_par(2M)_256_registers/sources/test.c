#include "keccak.h"

int main(int argc, char **argv)
{
   ALIGN char *msgstr[2];
   int inlen,i;
   ALIGN uint8_t *md[2];
   ALIGN uint8_t *md1[2];
      

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
   

    printf("\n 2-way implementation using 256-bit registers.\n\n");    
    printf("<------------------------------------------------------>\n");    

   
   printf("Keccak_256 \t\t\t");
   
   for(inlen=0;inlen<10000;inlen++){   
   keccak_256(msgstr,inlen,md,136);
   
   keccak_std(msgstr[0], inlen,md1[0],32);
   keccak_std(msgstr[1], inlen,md1[1],32);
   
   if((memcmp(md[0],md1[0],32) && memcmp(md[1],md1[1],32)!= 0)){
     printf("Error!!\n");
     exit(1);
   }
   }
    printf("Ok!!\n");

    
    printf("<------------------------------------------------------>\n\n\n");    
    
   for(i=0;i<2;i++){
	  _mm_free(msgstr[i]);
	  _mm_free(md[i]);
	  _mm_free(md1[i]);
    } 

  return 0;

}
