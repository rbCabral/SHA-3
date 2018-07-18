#include "keccak.h"

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
   keccak_std((uint8_t*)msgstr, inlen,md,32);
   
   if(memcmp(md,md1,32) == 0){
     printf("ok!\n");
    }else{
      printf("Error!!\n");
    }
    
    printf("<------------------------------------------------------>\n\n\n");    
    

      
  _mm_free(msgstr);

  return 0;

}
