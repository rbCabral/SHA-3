#include "keccak.h"

#define ROT(x,y)  _mm_rol_epi64(x,y)
#define ROT1(x)  _mm_rol_epi64(x,0x1)


#define load(x,m,d) 	  	x = _mm_load_si128((__m128i*)m+d)
#define store(x,d,y)		_mm_store_si128((__m128i*)x+d,y);
#define unpackLO_64(x,y,z)	x = _mm_unpacklo_epi64(y,z)
#define unpackHI_64(x,y,z)	x = _mm_unpackhi_epi64(y,z)
#define set_zero_128(x)		x = _mm_setzero_si128()
#define xor(x,y,z) 	  	x = _mm_xor_si128(y,z)
#define andnot(y,z)		_mm_andnot_si128(y,z)



#define init_S_zeros(y)		for(y=0;y<25;y++){\
				  set_zero_128(s[y]);}


#define join_last(y,p)		for(y=0;y<p;y++){\
				  load(ain[4],temp[0],y);\
				  load(ain[5],temp[1],y);\
				  \
				  unpackLO_64(ain[0],ain[4],ain[5]);\
				  unpackHI_64(ain[1],ain[4],ain[5]);\
				  \
				  xor(s[0+(y*2)],s[0+(y*2)],ain[0]);\
				  xor(s[1+(y*2)],s[1+(y*2)],ain[1]);\
				}\
				\
				load(ain[4],temp[0],y);\
				load(ain[5],temp[1],y);\
				\
				unpackLO_64(ain[0],ain[4],ain[5]);\
				\
				xor(s[0+(y*2)],s[0+(y*2)],ain[0]);
				  

#define join1(x,y,p)		for(y=0;y<p;y++){\
				  load(ain[4],x[0],y);\
				  load(ain[5],x[1],y);\
				  \
				  unpackLO_64(ain[0],ain[4],ain[5]);\
				  unpackHI_64(ain[1],ain[4],ain[5]);\
				  \
				  xor(s[0+(y*2)],s[0+(y*2)],ain[0]);\
				  xor(s[1+(y*2)],s[1+(y*2)],ain[1]);\
				}\
				\
				load(ain[4],x[0],y);\
				load(ain[5],x[1],y);\
				\
				unpackLO_64(ain[0],ain[4],ain[5]);\
				\
				xor(s[0+(y*2)],s[0+(y*2)],ain[0]);

#define join2(x,y,p)		load(ain[4],x[0],0);\
				load(ain[5],x[1],0);\
				unpackHI_64(ain[0],ain[4],ain[5]);\
				xor(s[0],s[0],ain[0]);\
				\
				for(y=0;y<p;y++){\
				  load(ain[4],x[0],y+1);\
				  load(ain[5],x[1],y+1);\
				  \
				  unpackLO_64(ain[0],ain[4],ain[5]);\
				  unpackHI_64(ain[1],ain[4],ain[5]);\
				  \
				  xor(s[1+(y*2)],s[1+(y*2)],ain[0]);\
				  xor(s[2+(y*2)],s[2+(y*2)],ain[1]);\
				}\
				
			
#define back(y,p)		for(y=0;y<p;y++){\
				  unpackLO_64(ain[0],s[0+(y*4)],s[1+(y*4)]);\
				  unpackHI_64(ain[1],s[0+(y*4)],s[1+(y*4)]);\
				  unpackLO_64(ain[2],s[2+(y*4)],s[3+(y*4)]);\
				  unpackHI_64(ain[3],s[2+(y*4)],s[3+(y*4)]);\
				  \
				  s[0+(y*4)] = ain[0];\
				  s[1+(y*4)] = ain[1];\
				  s[2+(y*4)] = ain[2];\
				  s[3+(y*4)] = ain[3];\
				  }


#define keccakf(s)      \
                        for(i=0;i<24;i=i+2){\
                        \
                        aux[0] = _mm_ternarylogic_epi64(s[0],s[5],s[10],0x96);\
                        aux[0] = _mm_ternarylogic_epi64(aux[0],s[15],s[20],0x96);\
                        \
                        aux[2] = _mm_ternarylogic_epi64(s[2],s[7],s[12],0x96);\
                        aux[2] = _mm_ternarylogic_epi64(aux[2],s[17],s[22],0x96);\
                        \
                        xor(aux[6],aux[0],ROT1(aux[2]));/*aux[6] = d1*/\
                        \
                        aux[1] = _mm_ternarylogic_epi64(s[1],s[6],s[11],0x96);\
                        aux[1] = _mm_ternarylogic_epi64(aux[1],s[16],s[21],0x96);\
                        \
                        aux[3] = _mm_ternarylogic_epi64(s[3],s[8],s[13],0x96);\
                        aux[3] = _mm_ternarylogic_epi64(aux[3],s[18],s[23],0x96);\
                        \
                        xor(aux[7],aux[1],ROT1(aux[3]));/*aux[7] = d2*/\
                        xor(aux[9],aux[3],ROT1(aux[0]));/*aux[9] = d4*/\
                        \
                        aux[4] = _mm_ternarylogic_epi64(s[4],s[9],s[14],0x96);\
                        aux[4] = _mm_ternarylogic_epi64(aux[4],s[19],s[24],0x96);\
                        \
                        xor(aux[5],aux[4],ROT1(aux[1]));/*aux[5] = d0*/\
                        xor(aux[8],aux[2],ROT1(aux[4]));/*aux[8] = d3*/\
                        \
                        xor(s[0],s[0],aux[5]);\
                        xor(s[5],s[5],aux[5]);\
                        xor(s[10],s[10],aux[5]);\
                        xor(s[15],s[15],aux[5]);\
                        xor(s[20],s[20],aux[5]);\
                        \
                        xor(s[1],s[1],aux[6]);\
                        xor(s[6],s[6],aux[6]);\
                        xor(s[11],s[11],aux[6]);\
                        xor(s[16],s[16],aux[6]);\
                        xor(s[21],s[21],aux[6]);\
                        \
                        xor(s[2],s[2],aux[7]);\
                        xor(s[7],s[7],aux[7]);\
                        xor(s[12],s[12],aux[7]);\
                        xor(s[17],s[17],aux[7]);\
                        xor(s[22],s[22],aux[7]);\
                        \
                        xor(s[3],s[3],aux[8]);\
                        xor(s[8],s[8],aux[8]);\
                        xor(s[13],s[13],aux[8]);\
                        xor(s[18],s[18],aux[8]);\
                        xor(s[23],s[23],aux[8]);\
                        \
                        xor(s[4],s[4],aux[9]);\
                        xor(s[9],s[9],aux[9]);\
                        xor(s[14],s[14],aux[9]);\
                        xor(s[19],s[19],aux[9]);\
                        xor(s[24],s[24],aux[9]);\
                        \
                        \
                        aux[1] = ROT(s[6],0x2C);\
                        aux[2] = ROT(s[12],0x2B);\
                        aux[3] = ROT(s[18],0x15);\
                        aux[4] = ROT(s[24],0x0E);\
                        \
                        \
                        s[18] = _mm_ternarylogic_epi64(aux[4],s[0],aux[1],0xD2);/*s[18] = s4*/\
                        s[6] = _mm_ternarylogic_epi64(aux[3],aux[4],s[0],0xD2);/*s[6] = s3*/\
                        s[0] = _mm_ternarylogic_epi64(s[0],aux[1],aux[2],0xD2);/*s[0] = s0*/\
                        s[12] = _mm_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[12] = s1*/\
                        load(aux[6],RC,i);\
                        s[24] = _mm_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[24] = s2*/\
                        \
                        xor(s[0],s[0],aux[6]);\
                        \
                        \
                        aux[0] = ROT(s[3],0x1C);\
                        aux[1] = ROT(s[9],0x14);\
                        aux[2] = ROT(s[10],0x3);\
                        aux[3] = ROT(s[16],0x2D);\
                        aux[4] = ROT(s[22],0x3D);\
                        \
                        s[16] = _mm_ternarylogic_epi64(aux[0],aux[1],aux[2],0xD2);/*s[16] = s5*/\
                        s[3]  = _mm_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[3] = s6*/\
                        s[10] = _mm_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[10] = s7*/\
                        s[22] = _mm_ternarylogic_epi64(aux[3],aux[4],aux[0],0xD2);/*s[22] = s8*/\
                        s[9]  = _mm_ternarylogic_epi64(aux[4],aux[0],aux[1],0xD2);/*s[9] = s9*/\
                        \
                        aux[0] = ROT(s[1],0x01);\
                        aux[1] = ROT(s[7],0x06);\
                        aux[2] = ROT(s[13],0x19);\
                        aux[3] = ROT(s[19],0x8);\
                        aux[4] = ROT(s[20],0x12);\
                        \
                        s[7]   = _mm_ternarylogic_epi64(aux[0],aux[1],aux[2],0xD2);/*s[7] = s10*/\
                        s[19]  = _mm_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[19] = s11*/\
                        s[1]   = _mm_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[1] = s12*/\
                        s[13]  = _mm_ternarylogic_epi64(aux[3],aux[4],aux[0],0xD2);/*s[13] = s13*/\
                        s[20]  = _mm_ternarylogic_epi64(aux[4],aux[0],aux[1],0xD2);/*s[20] = s14*/\
                        \
                        aux[0] = ROT(s[4],0x1B);\
                        aux[1] = ROT(s[5],0x24);\
                        aux[2] = ROT(s[11],0xA);\
                        aux[3] = ROT(s[17],0xF);\
                        aux[4] = ROT(s[23],0x38);\
                        \
                        s[23] = _mm_ternarylogic_epi64(aux[0],aux[1],aux[2],0xD2);/*s[23] = s15*/\
                        s[5]  = _mm_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[5] = s16*/\
                        s[17] = _mm_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[17] = s17*/\
                        s[4]  = _mm_ternarylogic_epi64(aux[3],aux[4],aux[0],0xD2);/*s[4] = s18*/\
                        s[11] = _mm_ternarylogic_epi64(aux[4],aux[0],aux[1],0xD2);/*s[11] = s19*/\
                        \
                        aux[0] = ROT(s[2],0x3E);\
                        aux[1] = ROT(s[8],0x37);\
                        aux[2] = ROT(s[14],0x27);\
                        aux[3] = ROT(s[15],0x29);\
                        aux[4] = ROT(s[21],0x2);\
                        \
                        s[14] = _mm_ternarylogic_epi64(aux[0],aux[1],aux[2],0xD2);/*s[14] = s20*/\
                        s[21] = _mm_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[21] = s21*/\
                        s[8]  = _mm_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[8] = s22*/\
                        s[15] = _mm_ternarylogic_epi64(aux[3],aux[4],aux[0],0xD2);/*s[15] = s23*/\
                        s[2]  = _mm_ternarylogic_epi64(aux[4],aux[0],aux[1],0xD2);/*s[2] = s24*/\
                        \
                        aux[0] = _mm_ternarylogic_epi64(s[0],s[16],s[7],0x96);\
                        aux[0] = _mm_ternarylogic_epi64(aux[0],s[23],s[14],0x96);\
                        \
                        aux[2] = _mm_ternarylogic_epi64(s[24],s[10],s[1],0x96);\
                        aux[2] = _mm_ternarylogic_epi64(aux[2],s[17],s[8],0x96);\
                        \
                        xor(aux[6],aux[0],ROT(aux[2],0x01));/*aux[6] = d1*/\
                        \
                        aux[1] = _mm_ternarylogic_epi64(s[12],s[3],s[19],0x96);\
                        aux[1] = _mm_ternarylogic_epi64(aux[1],s[5],s[21],0x96);\
                        \
                        aux[3] = _mm_ternarylogic_epi64(s[6],s[22],s[13],0x96);\
                        aux[3] = _mm_ternarylogic_epi64(aux[3],s[4],s[15],0x96);\
                        \
                        xor(aux[7],aux[1],ROT(aux[3],0x01));/*aux[7] = d2*/\
                        xor(aux[9],aux[3],ROT(aux[0],0x01));/*aux[9] = d4*/\
                        \
                        aux[4] = _mm_ternarylogic_epi64(s[18],s[9],s[20],0x96);\
                        aux[4] = _mm_ternarylogic_epi64(aux[4],s[11],s[2],0x96);\
                        \
                        xor(aux[5],aux[4],ROT(aux[1],0x01));/*aux[5] = d0*/\
                        xor(aux[8],aux[2],ROT(aux[4],0x01));/*aux[8] = d3*/\
                        \
                        \
                        xor(s[0],s[0],aux[5]);\
                        xor(s[16],s[16],aux[5]);\
                        xor(s[7],s[7],aux[5]);\
                        xor(s[23],s[23],aux[5]);\
                        xor(s[14],s[14],aux[5]);\
                        \
                        xor(s[12],s[12],aux[6]);\
                        xor(s[3],s[3],aux[6]);\
                        xor(s[19],s[19],aux[6]);\
                        xor(s[5],s[5],aux[6]);\
                        xor(s[21],s[21],aux[6]);\
                        \
                        xor(s[24],s[24],aux[7]);\
                        xor(s[10],s[10],aux[7]);\
                        xor(s[1],s[1],aux[7]);\
                        xor(s[17],s[17],aux[7]);\
                        xor(s[8],s[8],aux[7]);\
                        \
                        xor(s[6],s[6],aux[8]);\
                        xor(s[22],s[22],aux[8]);\
                        xor(s[13],s[13],aux[8]);\
                        xor(s[4],s[4],aux[8]);\
                        xor(s[15],s[15],aux[8]);\
                        \
                        xor(s[18],s[18],aux[9]);\
                        xor(s[9],s[9],aux[9]);\
                        xor(s[20],s[20],aux[9]);\
                        xor(s[11],s[11],aux[9]);\
                        xor(s[2],s[2],aux[9]);\
                        \
                        \
                        aux[1] = ROT(s[3],0x2C);\
                        aux[2] = ROT(s[1],0x2B);\
                        aux[3] = ROT(s[4],0x15);\
                        aux[4] = ROT(s[2],0x0E);\
                        \
                        s[3] = _mm_ternarylogic_epi64(aux[3],aux[4],s[0],0xD2);/*s[3] = s3*/\
                        s[4] = _mm_ternarylogic_epi64(aux[4],s[0],aux[1],0xD2);/*s[4] = s4*/\
                        s[0] = _mm_ternarylogic_epi64(s[0],aux[1],aux[2],0xD2);/*s[0] = s0*/\
                        s[1] = _mm_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[1] = s1*/\
                        load(aux[6],RC,(i+1));\
                        s[2] = _mm_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[2] = s2*/\
                        \
                        xor(s[0],s[0],aux[6]);\
                        \
                        aux[0] = ROT(s[6],0x1C);\
                        aux[1] = ROT(s[9],0x14);\
                        aux[2] = ROT(s[7],0x3);\
                        aux[3] = ROT(s[5],0x2D);\
                        aux[4] = ROT(s[8],0x3D);\
                        \
                        s[5] = _mm_ternarylogic_epi64(aux[0],aux[1],aux[2],0xD2);/*s[5] = s5*/\
                        s[6] = _mm_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[6] = s6*/\
                        s[7] = _mm_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[7] = s7*/\
                        s[8] = _mm_ternarylogic_epi64(aux[3],aux[4],aux[0],0xD2);/*s[8] = s8*/\
                        s[9] = _mm_ternarylogic_epi64(aux[4],aux[0],aux[1],0xD2);/*s[9] = s9*/\
                        \
                        aux[0] = ROT(s[12],0x1);\
                        aux[1] = ROT(s[10],0x6);\
                        aux[2] = ROT(s[13],0x19);\
                        aux[3] = ROT(s[11],0x8);\
                        aux[4] = ROT(s[14],0x12);\
                        \
                        s[10]  = _mm_ternarylogic_epi64(aux[0],aux[1],aux[2],0xD2);/*s[10] = s10*/\
                        s[11]  = _mm_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[11] = s11*/\
                        s[12]  = _mm_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[12] = s12*/\
                        s[13]  = _mm_ternarylogic_epi64(aux[3],aux[4],aux[0],0xD2);/*s[13] = s13*/\
                        s[14]  = _mm_ternarylogic_epi64(aux[4],aux[0],aux[1],0xD2);/*s[14] = s14*/\
                        \
                        aux[0] = ROT(s[18],0x1B);\
                        aux[1] = ROT(s[16],0x24);\
                        aux[2] = ROT(s[19],0xA);\
                        aux[3] = ROT(s[17],0xF);\
                        aux[4] = ROT(s[15],0x38);\
                        \
                        s[15] = _mm_ternarylogic_epi64(aux[0],aux[1],aux[2],0xD2);/*s[15] = s15*/\
                        s[16]  = _mm_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[16] = s16*/\
                        s[17] = _mm_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[17] = s17*/\
                        s[18]  = _mm_ternarylogic_epi64(aux[3],aux[4],aux[0],0xD2);/*s[18] = s18*/\
                        s[19] = _mm_ternarylogic_epi64(aux[4],aux[0],aux[1],0xD2);/*s[19] = s19*/\
                        \
                        aux[0] = ROT(s[24],0x3E);\
                        aux[1] = ROT(s[22],0x37);\
                        aux[2] = ROT(s[20],0x27);\
                        aux[3] = ROT(s[23],0x29);\
                        aux[4] = ROT(s[21],0x2);\
                        \
                        s[20] = _mm_ternarylogic_epi64(aux[0],aux[1],aux[2],0xD2);/*s[20] = s20*/\
                        s[21] = _mm_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[21] = s21*/\
                        s[22] = _mm_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[22] = s22*/\
                        s[23] = _mm_ternarylogic_epi64(aux[3],aux[4],aux[0],0xD2);/*s[23] = s23*/\
                        s[24] = _mm_ternarylogic_epi64(aux[4],aux[0],aux[1],0xD2);/*s[24] = s24*/\
                        }\
                        \


void print_128(__m128i x){
     uint64_t c[2];
     store(c,0,x);
     printf("\n 1 = %8.16lx\t  0 = %8.16lx\t\n", c[1],c[0]);
}			
	


int keccak(char **in_e, int inlen, uint8_t **md, int rsiz)
{
   __m128i aux[10];
   int i;
    const ALIGN uint64_t RC[100] = {
    0x0000000000000001,0x0000000000000001,
    0x0000000000008082,0x0000000000008082,
    0x800000000000808A,0x800000000000808A,
    0x8000000080008000,0x8000000080008000,
    0x000000000000808B,0x000000000000808B,
    0x0000000080000001,0x0000000080000001,
    0x8000000080008081,0x8000000080008081,
    0x8000000000008009,0x8000000000008009,
    0x000000000000008A,0x000000000000008A,
    0x0000000000000088,0x0000000000000088,
    0x0000000080008009,0x0000000080008009,
    0x000000008000000A,0x000000008000000A,
    0x000000008000808B,0x000000008000808B,
    0x800000000000008B,0x800000000000008B,
    0x8000000000008089,0x8000000000008089,
    0x8000000000008003,0x8000000000008003,
    0x8000000000008002,0x8000000000008002,
    0x8000000000000080,0x8000000000000080,
    0x000000000000800A,0x000000000000800A,
    0x800000008000000A,0x800000008000000A,
    0x8000000080008081,0x8000000080008081,
    0x8000000000008080,0x8000000000008080,
    0x0000000080000001,0x0000000080000001,
    0x8000000080008008,0x8000000080008008
  };


  ALIGN uint8_t temp[2][144];
    
  int j,l=0,it=0,k,o=0;
  __m128i s[28],ain[8];
  __m128i auxl[5];
  char* con[32];
  char* in[2];
  int t,g,gf;
  
  in[0]=in_e[0];
  in[1]=in_e[1];

  init_S_zeros(j);
  
  if(inlen >= rsiz ){
	it = inlen / rsiz;
  }
      
  switch (rsiz){
    case 136:
      t  =   8;
      g  = 128;
      gf = 144;
      break;
    case 104: 
      t  =   6;
      g  =  96;
      gf = 112;
      break;
    case 72: 
      t  =   4;
      g  =  64;
      gf =  80;
      break;
  }    
      
  for(k=0;k<it/2;k++){
    join1(in,j,t);


    keccakf(s);
    
    
    in[0] += g;
    in[1] += g;
    join2(in,j,t);
      
    keccakf(s);
    in[0] += gf;
    in[1] += gf;  
  }
      
  if(it%2 != 0 ){
    join1(in,j,t);
    
    keccakf(s);
    in[0] += rsiz;
    in[1] += rsiz;
  }   
   inlen = inlen - it * rsiz;
  for(i=0;i<2;i++){
    memset(temp[i], 0, 144*sizeof(uint8_t));
    memcpy(temp[i], in[i], inlen);
    temp[i][inlen] = 0x06;
    temp[i][rsiz - 1] |= 0x80;
  }

  join_last(j,t); 
  keccakf(s);
  


  switch (rsiz){
    case 136: 
      back(j,1);
      store(md[0],0,s[0]);
      store(md[0],1,s[2]);
      store(md[1],0,s[1]);
      store(md[1],1,s[3]);
      break;
    case 104:
      back(j,2);
      store(md[0],0,s[0]);
      store(md[0],1,s[2]);
      store(md[1],0,s[1]);
      store(md[1],1,s[3]);
      store(md[0],2,s[4]);
      store(md[0],3,s[6]);
      store(md[1],2,s[5]);
      store(md[1],3,s[7]);
      break;
    case 72:
      back(j,2);
      store(md[0],0,s[0]);
      store(md[0],1,s[2]);
      store(md[1],0,s[1]);
      store(md[1],1,s[3]);
      store(md[0],2,s[4]);
      store(md[0],3,s[6]);
      store(md[1],2,s[5]);
      store(md[1],3,s[7]);
      break;
  } 
  

  return 0;
}
