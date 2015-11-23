#include "keccak.h"

#define ROT(x,y)  _mm256_or_si256(_mm256_slli_epi64(x,y),_mm256_srli_epi64(x,(64-y)))
#define ROT1(x)  _mm256_or_si256(_mm256_add_epi64(x,x),_mm256_srli_epi64(x,(63)))

#define load(x,m,d) 	  	x = _mm256_load_si256((__m256i*)m+d)
#define store(x,d,y)		_mm256_store_si256((__m256i*)x+d,y);
#define blend_32(x,y,z,d) 	x = _mm256_blend_epi32(y,z,d)
#define permut_64(x,y,d)  	x = _mm256_permute4x64_epi64(y,d)
#define xor_pd(x,y,z) 	  	x = _mm256_xor_pd(y,z)
#define xor(x,y,z) 	  	x = _mm256_xor_si256(y,z)
#define andnot(y,z)		_mm256_andnot_si256(y,z)
#define unpackLO_64(x,y,z)	x = _mm256_unpacklo_epi64(y,z)
#define unpackHI_64(x,y,z)	x = _mm256_unpackhi_epi64(y,z)
#define extract_128(x,y,d)	x = _mm256_extracti128_si256(y,d)
#define insert_128(x,y,z,d)	x = _mm256_inserti128_si256(y,z,d)
#define set_zero_256(x)		x = _mm256_setzero_si256()

#define load_128(x,m,d) 	x = _mm_load_si128((__m128i*)m+d)
#define store_128(x,d,y)	_mm_store_si128((__m128i*)x+d,y);
#define unpackLO_128_64(x,y,z)	x = _mm_unpacklo_epi64(y,z)
#define unpackHI_128_64(x,y,z)	x = _mm_unpackhi_epi64(y,z)
#define cast_to_256(x,y)	x = _mm256_castsi128_si256(y)

#define init_S_zeros(y)		for(y=0;y<25;y++){\
				  set_zero_256(s[y]);}

#define join_last(y,p)		\
				for(y=0;y<p;y++){\
				  load(ain[4],temp0,y);\
				  load(ain[5],temp1,y);\
				  load(ain[6],temp2,y);\
				  load(ain[7],temp3,y);\
				  \
				  unpackLO_64(ain[0],ain[4],ain[5]);\
				  unpackHI_64(ain[1],ain[4],ain[5]);\
				  unpackLO_64(ain[2],ain[6],ain[7]);\
				  unpackHI_64(ain[3],ain[6],ain[7]);\
				  \
				  extract_128(auxl[0],ain[0],1);\
				  extract_128(auxl[1],ain[1],1);\
				  extract_128(auxl[2],ain[2],0);\
				  extract_128(auxl[3],ain[3],0);\
				  cast_to_256(ain[4],auxl[1]);\
				  cast_to_256(ain[5],auxl[0]);\
				  \
				  blend_32(ain[3],ain[3],ain[4],0x0F);\
				  blend_32(ain[2],ain[2],ain[5],0x0F);\
				  \
				  /*insert_128(ain[2],ain[2],auxl[0],0);\
				  insert_128(ain[3],ain[3],auxl[1],0);\*/\
				  insert_128(ain[0],ain[0],auxl[2],1);\
				  insert_128(ain[1],ain[1],auxl[3],1);\
				  \
				  xor(s[2+(y*4)],s[2+(y*4)],ain[2]);\
				  xor(s[3+(y*4)],s[3+(y*4)],ain[3]);\
				  xor(s[0+(y*4)],s[0+(y*4)],ain[0]);\
				  xor(s[1+(y*4)],s[1+(y*4)],ain[1]);\
				  }\
				load_128(auxl[0],temp0,y*2);\
				load_128(auxl[1],temp1,y*2);\
				load_128(auxl[2],temp2,y*2);\
				load_128(auxl[3],temp3,y*2);\
				\
				unpackLO_128_64(auxl[4],auxl[0],auxl[1]);\
				unpackLO_128_64(auxl[0],auxl[2],auxl[3]);\
				cast_to_256(ain[0],auxl[4]);\
				insert_128(ain[1],ain[0],auxl[0],1);\
				xor(s[0+(y*4)],s[0+(y*4)],ain[1]);
				  

#define join1(x,y,p)		   for(y=0;y<p;y++){\
				  load(ain[4],x[0],y);\
				  load(ain[5],x[1],y);\
				  load(ain[6],x[2],y);\
				  load(ain[7],x[3],y);\
				  \
				  unpackLO_64(ain[0],ain[4],ain[5]);\
				  unpackHI_64(ain[1],ain[4],ain[5]);\
				  unpackLO_64(ain[2],ain[6],ain[7]);\
				  unpackHI_64(ain[3],ain[6],ain[7]);\
				  \
				  extract_128(auxl[0],ain[0],1);\
				  extract_128(auxl[1],ain[1],1);\
				  extract_128(auxl[2],ain[2],0);\
				  extract_128(auxl[3],ain[3],0);\
				  cast_to_256(ain[4],auxl[1]);\
				  cast_to_256(ain[5],auxl[0]);\
				  \
				  blend_32(ain[3],ain[3],ain[4],0x0F);\
				  blend_32(ain[2],ain[2],ain[5],0x0F);\
				  \
				  /*insert_128(ain[2],ain[2],auxl[0],0);\
				  insert_128(ain[3],ain[3],auxl[1],0);\*/\
				  insert_128(ain[0],ain[0],auxl[2],1);\
				  insert_128(ain[1],ain[1],auxl[3],1);\
				  \
				  xor(s[2+(y*4)],s[2+(y*4)],ain[2]);\
				  xor(s[3+(y*4)],s[3+(y*4)],ain[3]);\
				  xor(s[0+(y*4)],s[0+(y*4)],ain[0]);\
				  xor(s[1+(y*4)],s[1+(y*4)],ain[1]);\
				  }\
				load_128(auxl[0],x[0],y*2);\
				load_128(auxl[1],x[1],y*2);\
				load_128(auxl[2],x[2],y*2);\
				load_128(auxl[3],x[3],y*2);\
				\
				unpackLO_128_64(auxl[4],auxl[0],auxl[1]);\
				unpackLO_128_64(auxl[0],auxl[2],auxl[3]);\
				cast_to_256(ain[0],auxl[4]);\
				insert_128(ain[1],ain[0],auxl[0],1);\
				xor(s[0+(y*4)],s[0+(y*4)],ain[1]);

#define join2(x,y,p)		load(ain[4],x[0],0);\
				load(ain[5],x[1],0);\
				load(ain[6],x[2],0);\
				load(ain[7],x[3],0);\
				\
				unpackLO_64(ain[0],ain[4],ain[5]);\
				unpackHI_64(ain[1],ain[4],ain[5]);\
				unpackLO_64(ain[2],ain[6],ain[7]);\
				unpackHI_64(ain[3],ain[6],ain[7]);\
				\
				extract_128(auxl[0],ain[0],1);\
				extract_128(auxl[1],ain[1],1);\
				extract_128(auxl[3],ain[3],0);\
				\
				insert_128(ain[2],ain[2],auxl[0],0);\
				insert_128(ain[3],ain[3],auxl[1],0);\
				insert_128(ain[1],ain[1],auxl[3],1);\
				\
				xor(s[1],s[1],ain[2]);\
				xor(s[2],s[2],ain[3]);\
				xor(s[0],s[0],ain[1]);\
				\
				for(y=0;y<p;y++){\
				  load(ain[4],x[0],y+1);\
				  load(ain[5],x[1],y+1);\
				  load(ain[6],x[2],y+1);\
				  load(ain[7],x[3],y+1);\
				  \
				  unpackLO_64(ain[0],ain[4],ain[5]);\
				  unpackHI_64(ain[1],ain[4],ain[5]);\
				  unpackLO_64(ain[2],ain[6],ain[7]);\
				  unpackHI_64(ain[3],ain[6],ain[7]);\
				  \
				  extract_128(auxl[0],ain[0],1);\
				  extract_128(auxl[1],ain[1],1);\
				  extract_128(auxl[2],ain[2],0);\
				  extract_128(auxl[3],ain[3],0);\
				  cast_to_256(ain[4],auxl[1]);\
				  cast_to_256(ain[5],auxl[0]);\
				  \
				  blend_32(ain[3],ain[3],ain[4],0x0F);\
				  blend_32(ain[2],ain[2],ain[5],0x0F);\
				  \
				  /*insert_128(ain[2],ain[2],auxl[0],0);\
				  insert_128(ain[3],ain[3],auxl[1],0);\*/\
				  insert_128(ain[0],ain[0],auxl[2],1);\
				  insert_128(ain[1],ain[1],auxl[3],1);\
				  \
				  xor(s[5+(y*4)],s[5+(y*4)],ain[2]);\
				  xor(s[6+(y*4)],s[6+(y*4)],ain[3]);\
				  xor(s[3+(y*4)],s[3+(y*4)],ain[0]);\
				  xor(s[4+(y*4)],s[4+(y*4)],ain[1]);\
				}\
				load_128(auxl[0],x[0],(y+1)*2);\
				load_128(auxl[1],x[1],(y+1)*2);\
				load_128(auxl[2],x[2],(y+1)*2);\
				load_128(auxl[3],x[3],(y+1)*2);\
				\
				unpackLO_128_64(auxl[4],auxl[0],auxl[1]);\
				unpackHI_128_64(auxl[1],auxl[0],auxl[1]);\
				unpackLO_128_64(auxl[0],auxl[2],auxl[3]);\
				unpackHI_128_64(auxl[2],auxl[2],auxl[3]);\
				\
				cast_to_256(ain[0],auxl[4]);\
				cast_to_256(ain[1],auxl[1]);\
				insert_128(ain[0],ain[0],auxl[0],1);\
				insert_128(ain[1],ain[1],auxl[2],1);\
				\
				xor(s[3+(y*4)],s[3+(y*4)],ain[0]);\
				xor(s[4+(y*4)],s[4+(y*4)],ain[1]);
				
				
				
#define join3(x,y,p)		load_128(auxl[0],x[0],1);\
				load_128(auxl[1],x[1],1);\
				load_128(auxl[2],x[2],1);\
				load_128(auxl[3],x[3],1);\
				\
				unpackLO_128_64(auxl[4],auxl[0],auxl[1]);\
				unpackHI_128_64(auxl[1],auxl[0],auxl[1]);\
				unpackLO_128_64(auxl[0],auxl[2],auxl[3]);\
				unpackHI_128_64(auxl[2],auxl[2],auxl[3]);\
				\
				cast_to_256(ain[0],auxl[4]);\
				cast_to_256(ain[1],auxl[1]);\
				insert_128(ain[0],ain[0],auxl[0],1);\
				insert_128(ain[1],ain[1],auxl[2],1);\
				\
				xor(s[0],s[0],ain[0]);\
				xor(s[1],s[1],ain[1]);\
				\
				for(y=0;y<p;y++){\
				  load(ain[4],x[0],y+1);\
				  load(ain[5],x[1],y+1);\
				  load(ain[6],x[2],y+1);\
				  load(ain[7],x[3],y+1);\
				  \
				  unpackLO_64(ain[0],ain[4],ain[5]);\
				  unpackHI_64(ain[1],ain[4],ain[5]);\
				  unpackLO_64(ain[2],ain[6],ain[7]);\
				  unpackHI_64(ain[3],ain[6],ain[7]);\
				  \
				  extract_128(auxl[0],ain[0],1);\
				  extract_128(auxl[1],ain[1],1);\
				  extract_128(auxl[2],ain[2],0);\
				  extract_128(auxl[3],ain[3],0);\
				  cast_to_256(ain[4],auxl[1]);\
				  cast_to_256(ain[5],auxl[0]);\
				  \
				  blend_32(ain[3],ain[3],ain[4],0x0F);\
				  blend_32(ain[2],ain[2],ain[5],0x0F);\
				  \
				  /*insert_128(ain[2],ain[2],auxl[0],0);\
				  insert_128(ain[3],ain[3],auxl[1],0);\*/\
				  insert_128(ain[0],ain[0],auxl[2],1);\
				  insert_128(ain[1],ain[1],auxl[3],1);\
				  \
				  xor(s[4+(y*4)],s[4+(y*4)],ain[2]);\
				  xor(s[5+(y*4)],s[5+(y*4)],ain[3]);\
				  xor(s[2+(y*4)],s[2+(y*4)],ain[0]);\
				  xor(s[3+(y*4)],s[3+(y*4)],ain[1]);\
				}\
				load(ain[4],x[0],y+1);\
				load(ain[5],x[1],y+1);\
				load(ain[6],x[2],y+1);\
				load(ain[7],x[3],y+1);\
				\
				unpackLO_64(ain[0],ain[4],ain[5]);\
				unpackHI_64(ain[1],ain[4],ain[5]);\
				unpackLO_64(ain[2],ain[6],ain[7]);\
				unpackHI_64(ain[3],ain[6],ain[7]);\
				\
				extract_128(auxl[0],ain[0],1);\
				extract_128(auxl[2],ain[2],0);\
				extract_128(auxl[3],ain[3],0);\
				\
				insert_128(ain[0],ain[0],auxl[2],1);\
				insert_128(ain[1],ain[1],auxl[3],1);\
				insert_128(ain[2],ain[2],auxl[0],0);\
				\
				xor(s[2+(y*4)],s[2+(y*4)],ain[0]);\
				xor(s[3+(y*4)],s[3+(y*4)],ain[1]);\
				xor(s[4+(y*4)],s[4+(y*4)],ain[2]);

#define join4(x,y,p)		load_128(auxl[0],x[0],1);\
				load_128(auxl[1],x[1],1);\
				load_128(auxl[2],x[2],1);\
				load_128(auxl[3],x[3],1);\
				\
				unpackHI_128_64(auxl[1],auxl[0],auxl[1]);\
				unpackHI_128_64(auxl[2],auxl[2],auxl[3]);\
				\
				cast_to_256(ain[0],auxl[1]);\
				insert_128(ain[0],ain[0],auxl[2],1);\
				xor(s[0],s[0],ain[0]);\
				\
				for(y=0;y<p;y++){\
				  load(ain[4],x[0],y+1);\
				  load(ain[5],x[1],y+1);\
				  load(ain[6],x[2],y+1);\
				  load(ain[7],x[3],y+1);\
				  \
				  unpackLO_64(ain[0],ain[4],ain[5]);\
				  unpackHI_64(ain[1],ain[4],ain[5]);\
				  unpackLO_64(ain[2],ain[6],ain[7]);\
				  unpackHI_64(ain[3],ain[6],ain[7]);\
				  \
				  extract_128(auxl[0],ain[0],1);\
				  extract_128(auxl[1],ain[1],1);\
				  extract_128(auxl[2],ain[2],0);\
				  extract_128(auxl[3],ain[3],0);\
				  cast_to_256(ain[4],auxl[1]);\
				  cast_to_256(ain[5],auxl[0]);\
				  \
				  blend_32(ain[3],ain[3],ain[4],0x0F);\
				  blend_32(ain[2],ain[2],ain[5],0x0F);\
				  \
				  /*insert_128(ain[2],ain[2],auxl[0],0);\
				  insert_128(ain[3],ain[3],auxl[1],0);\*/\
				  insert_128(ain[0],ain[0],auxl[2],1);\
				  insert_128(ain[1],ain[1],auxl[3],1);\
				  \
				  xor(s[3+(y*4)],s[3+(y*4)],ain[2]);\
				  xor(s[4+(y*4)],s[4+(y*4)],ain[3]);\
				  xor(s[1+(y*4)],s[1+(y*4)],ain[0]);\
				  xor(s[2+(y*4)],s[2+(y*4)],ain[1]);\
				}
				
#define back(y,p)		for(y=0;y<p;y++){\
				  unpackLO_64(ain[0],s[0+(y*4)],s[1+(y*4)]);\
				  unpackHI_64(ain[1],s[0+(y*4)],s[1+(y*4)]);\
				  unpackLO_64(ain[2],s[2+(y*4)],s[3+(y*4)]);\
				  unpackHI_64(ain[3],s[2+(y*4)],s[3+(y*4)]);\
				  \
				  extract_128(auxl[0],ain[0],1);\
				  extract_128(auxl[1],ain[1],1);\
				  extract_128(auxl[2],ain[2],0);\
				  extract_128(auxl[3],ain[3],0);\
				  cast_to_256(ain[4],auxl[1]);\
				  cast_to_256(ain[5],auxl[0]);\
				  \
				  blend_32(s[3+(y*4)],ain[3],ain[4],0x0F);\
				  blend_32(s[2+(y*4)],ain[2],ain[5],0x0F);\
				  \
				  /*insert_128(s[2+(y*4)],ain[2],auxl[0],0);\
				  insert_128(s[3+(y*4)],ain[3],auxl[1],0);\*/\
				  insert_128(s[0+(y*4)],ain[0],auxl[2],1);\
				  insert_128(s[1+(y*4)],ain[1],auxl[3],1);\
				  }

#define keccakf(s)	\
			for(i=0;i<24;i=i+2){\
      \
	xor(aux[0],s[0],s[5]);\
	xor(aux[0],aux[0],s[10]);\
	xor(aux[0],aux[0],s[15]);\
	xor(aux[0],aux[0],s[20]);\
	\
	xor(aux[2],s[2],s[7]);\
	xor(aux[2],aux[2],s[12]);\
	xor(aux[2],aux[2],s[17]);\
	xor(aux[2],aux[2],s[22]);\
	\
	xor(aux[6],aux[0],ROT1(aux[2]));/*aux[6] = d1*/\
	\
	xor(aux[1],s[1],s[6]);\
	xor(aux[1],aux[1],s[11]);\
	xor(aux[1],aux[1],s[16]);\
	xor(aux[1],aux[1],s[21]);\
	\
	xor(aux[3],s[3],s[8]);\
	xor(aux[3],aux[3],s[13]);\
	xor(aux[3],aux[3],s[18]);\
	xor(aux[3],aux[3],s[23]);\
	\
	xor(aux[7],aux[1],ROT1(aux[3]));/*aux[7] = d2*/\
	xor(aux[9],aux[3],ROT1(aux[0]));/*aux[9] = d4*/\
	\
	xor(aux[4],s[4],s[9]);\
	xor(aux[4],aux[4],s[14]);\
	xor(aux[4],aux[4],s[19]);\
	xor(aux[4],aux[4],s[24]);\
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
	aux[1] = ROT(s[6],44);\
	aux[2] = ROT(s[12],43);\
	aux[3] = ROT(s[18],21);\
	aux[4] = ROT(s[24],14);\
	\
	\
	xor(s[18],aux[4],andnot(s[0],aux[1]));/*s[18] = s4*/\
	xor(s[6],aux[3],andnot(aux[4],s[0]));/*s[6] = s3*/\
	xor(s[0] ,s[0],andnot(aux[1],aux[2]));/*s[0] = s0*/\
	xor(s[12] ,aux[1],andnot(aux[2],aux[3]));/*s[12] = s1*/\
	load(aux[6],RC,i);\
	xor(s[24],aux[2],andnot(aux[3],aux[4]));/*s[24] = s2*/\
	\
	\
	\
	xor(s[0],s[0],aux[6]);\
	\
	\
	aux[0] = ROT(s[3],28);\
	aux[1] = ROT(s[9],20);\
	aux[2] = ROT(s[10],3);\
	aux[3] = ROT(s[16],45);\
	aux[4] = ROT(s[22],61);\
	\
	xor(s[16] ,aux[0],andnot(aux[1],aux[2]));/*s[16] = s5*/\
	xor(s[3] ,aux[1],andnot(aux[2],aux[3]));/*s[3] = s6*/\
	xor(s[10],aux[2],andnot(aux[3],aux[4]));/*s[10]= s7*/\
	xor(s[22],aux[3],andnot(aux[4],aux[0]));/*s[22]= s8*/\
	xor(s[9],aux[4],andnot(aux[0],aux[1]));/*s[9]= s9*/\
	\
	aux[0] = ROT1(s[1]);\
	aux[1] = ROT(s[7],6);\
	aux[2] = ROT(s[13],25);\
	aux[3] = ROT(s[19],8);\
	aux[4] = ROT(s[20],18);\
	\
	xor(s[7] ,aux[0],andnot(aux[1],aux[2]));/*s[7] = s10*/\
	xor(s[19] ,aux[1],andnot(aux[2],aux[3]));/*s[19] = s11*/\
	xor(s[1],aux[2],andnot(aux[3],aux[4]));/*s[1]= s12*/\
	xor(s[13],aux[3],andnot(aux[4],aux[0]));/*s[13]= s13*/\
	xor(s[20],aux[4],andnot(aux[0],aux[1]));/*s[20]= s14*/\
	\
	aux[0] = ROT(s[4],27);\
	aux[1] = ROT(s[5],36);\
	aux[2] = ROT(s[11],10);\
	aux[3] = ROT(s[17],15);\
	aux[4] = ROT(s[23],56);\
	\
	xor(s[23] ,aux[0],andnot(aux[1],aux[2]));/*s[23] = s15*/\
	xor(s[5] ,aux[1],andnot(aux[2],aux[3]));/*s[5] = s16*/\
	xor(s[17],aux[2],andnot(aux[3],aux[4]));/*s[17]= s17*/\
	xor(s[4],aux[3],andnot(aux[4],aux[0]));/*s[4]= s18*/\
	xor(s[11],aux[4],andnot(aux[0],aux[1]));/*s[11]= s19*/\
	\
	aux[0] = ROT(s[2],62);\
	aux[1] = ROT(s[8],55);\
	aux[2] = ROT(s[14],39);\
	aux[3] = ROT(s[15],41);\
	aux[4] = ROT(s[21],2);\
	\
	xor(s[14] ,aux[0],andnot(aux[1],aux[2]));/*s[14] = s20*/\
	xor(s[21] ,aux[1],andnot(aux[2],aux[3]));/*s[21] = s21*/\
	xor(s[8],aux[2],andnot(aux[3],aux[4]));/*s[8]= s22*/\
	xor(s[15],aux[3],andnot(aux[4],aux[0]));/*s[15]= s23*/\
	xor(s[2],aux[4],andnot(aux[0],aux[1]));/*s[2]= s24*/\
	\
	xor(aux[0],s[0],s[16]);\
	xor(aux[0],aux[0],s[7]);\
	xor(aux[0],aux[0],s[23]);\
	xor(aux[0],aux[0],s[14]);\
	\
	xor(aux[2],s[24],s[10]);\
	xor(aux[2],aux[2],s[1]);\
	xor(aux[2],aux[2],s[17]);\
	xor(aux[2],aux[2],s[8]);\
	\
	xor(aux[6],aux[0],ROT(aux[2],0x01));/*aux[6] = d1*/\
	\
	xor(aux[1],s[12],s[3]);\
	xor(aux[1],aux[1],s[19]);\
	xor(aux[1],aux[1],s[5]);\
	xor(aux[1],aux[1],s[21]);\
	\
	xor(aux[3],s[6],s[22]);\
	xor(aux[3],aux[3],s[13]);\
	xor(aux[3],aux[3],s[4]);\
	xor(aux[3],aux[3],s[15]);\
	\
	xor(aux[7],aux[1],ROT(aux[3],0x01));/*aux[7] = d2*/\
	xor(aux[9],aux[3],ROT(aux[0],0x01));/*aux[9] = d4*/\
	\
	xor(aux[4],s[18],s[9]);\
	xor(aux[4],aux[4],s[20]);\
	xor(aux[4],aux[4],s[11]);\
	xor(aux[4],aux[4],s[2]);\
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
	aux[1] = ROT(s[3],44);\
	aux[2] = ROT(s[1],43);\
	aux[3] = ROT(s[4],21);\
	aux[4] = ROT(s[2],14);\
	\
	\
	xor(s[3],aux[3],andnot(aux[4],s[0]));/*s[3] = s3*/\
	xor(s[4],aux[4],andnot(s[0],aux[1]));/*s[4] = s4*/\
	xor(s[0] ,s[0],andnot(aux[1],aux[2]));/*s[0] = s0*/\
	xor(s[1] ,aux[1],andnot(aux[2],aux[3]));/*s[1] = s1*/\
	load(aux[6],RC,(i+1));\
	xor(s[2],aux[2],andnot(aux[3],aux[4]));/*s[2] = s2*/\
	\
	\
	xor(s[0],s[0],aux[6]);\
	\
	aux[0] = ROT(s[6],28);\
	aux[1] = ROT(s[9],20);\
	aux[2] = ROT(s[7],3);\
	aux[3] = ROT(s[5],45);\
	aux[4] = ROT(s[8],61);\
	\
	xor(s[5] ,aux[0],andnot(aux[1],aux[2]));/*s[5] = s5*/\
	xor(s[6] ,aux[1],andnot(aux[2],aux[3]));/*s[6] = s6*/\
	xor(s[7],aux[2],andnot(aux[3],aux[4]));/*s[7]= s7*/\
	xor(s[8],aux[3],andnot(aux[4],aux[0]));/*s[8]= s8*/\
	xor(s[9],aux[4],andnot(aux[0],aux[1]));/*s[9]= s9*/\
	\
	aux[0] = ROT(s[12],1);\
	aux[1] = ROT(s[10],6);\
	aux[2] = ROT(s[13],25);\
	aux[3] = ROT(s[11],8);\
	aux[4] = ROT(s[14],18);\
	\
	xor(s[10] ,aux[0],andnot(aux[1],aux[2]));/*s[10] = s10*/\
	xor(s[11] ,aux[1],andnot(aux[2],aux[3]));/*s[11] = s11*/\
	xor(s[12],aux[2],andnot(aux[3],aux[4]));/*s[12]= s12*/\
	xor(s[13],aux[3],andnot(aux[4],aux[0]));/*s[13]= s13*/\
	xor(s[14],aux[4],andnot(aux[0],aux[1]));/*s[14]= s14*/\
	\
	aux[0] = ROT(s[18],27);\
	aux[1] = ROT(s[16],36);\
	aux[2] = ROT(s[19],10);\
	aux[3] = ROT(s[17],15);\
	aux[4] = ROT(s[15],56);\
	\
	xor(s[15] ,aux[0],andnot(aux[1],aux[2]));/*s[15] = s15*/\
	xor(s[16] ,aux[1],andnot(aux[2],aux[3]));/*s[16] = s16*/\
	xor(s[17],aux[2],andnot(aux[3],aux[4]));/*s[17]= s17*/\
	xor(s[18],aux[3],andnot(aux[4],aux[0]));/*s[18]= s18*/\
	xor(s[19],aux[4],andnot(aux[0],aux[1]));/*s[19]= s19*/\
	\
	aux[0] = ROT(s[24],62);\
	aux[1] = ROT(s[22],55);\
	aux[2] = ROT(s[20],39);\
	aux[3] = ROT(s[23],41);\
	aux[4] = ROT(s[21],2);\
	\
	xor(s[20] ,aux[0],andnot(aux[1],aux[2]));/*s[20] = s20*/\
	xor(s[21] ,aux[1],andnot(aux[2],aux[3]));/*s[21] = s21*/\
	xor(s[22],aux[2],andnot(aux[3],aux[4]));/*s[22]= s22*/\
	xor(s[23],aux[3],andnot(aux[4],aux[0]));/*s[23]= s23*/\
	xor(s[24],aux[4],andnot(aux[0],aux[1]));/*s[24]= s24*/\
	}\
	\




int keccak_256(char **in_e, int inlen, uint8_t *md,uint8_t *md1,uint8_t *md2,uint8_t *md3)
{
   __m256i aux[10];
   int i;
  const ALIGN uint64_t RC[100] = {
    0x0000000000000001,0x0000000000000001,0x0000000000000001,0x0000000000000001,
    0x0000000000008082,0x0000000000008082,0x0000000000008082,0x0000000000008082,
    0x800000000000808A,0x800000000000808A,0x800000000000808A,0x800000000000808A,  
    0x8000000080008000,0x8000000080008000,0x8000000080008000,0x8000000080008000, 
    0x000000000000808B,0x000000000000808B,0x000000000000808B,0x000000000000808B,   
    0x0000000080000001,0x0000000080000001,0x0000000080000001,0x0000000080000001, 
    0x8000000080008081,0x8000000080008081,0x8000000080008081,0x8000000080008081, 
    0x8000000000008009,0x8000000000008009,0x8000000000008009,0x8000000000008009, 
    0x000000000000008A,0x000000000000008A,0x000000000000008A,0x000000000000008A, 
    0x0000000000000088,0x0000000000000088,0x0000000000000088,0x0000000000000088, 
    0x0000000080008009,0x0000000080008009,0x0000000080008009,0x0000000080008009, 
    0x000000008000000A,0x000000008000000A,0x000000008000000A,0x000000008000000A, 
    0x000000008000808B,0x000000008000808B,0x000000008000808B,0x000000008000808B, 
    0x800000000000008B,0x800000000000008B,0x800000000000008B,0x800000000000008B, 
    0x8000000000008089,0x8000000000008089,0x8000000000008089,0x8000000000008089, 
    0x8000000000008003,0x8000000000008003,0x8000000000008003,0x8000000000008003, 
    0x8000000000008002,0x8000000000008002,0x8000000000008002,0x8000000000008002, 
    0x8000000000000080,0x8000000000000080,0x8000000000000080,0x8000000000000080, 
    0x000000000000800A,0x000000000000800A,0x000000000000800A,0x000000000000800A, 
    0x800000008000000A,0x800000008000000A,0x800000008000000A,0x800000008000000A, 
    0x8000000080008081,0x8000000080008081,0x8000000080008081,0x8000000080008081, 
    0x8000000000008080,0x8000000000008080,0x8000000000008080,0x8000000000008080, 
    0x0000000080000001,0x0000000080000001,0x0000000080000001,0x0000000080000001, 
    0x8000000080008008,0x8000000080008008,0x8000000080008008,0x8000000080008008 
  };
  
  const int rsiz = 136;
  const int g = 128;
  const int gf = 160;

  ALIGN uint8_t temp0[144];
  ALIGN uint8_t temp1[144];
  ALIGN uint8_t temp2[144];
  ALIGN uint8_t temp3[144];
  int k=1,j,l=0;
  __m256i s[28],ain[8];
  __m128i auxl[5];
  char* con[32];
  char* in[4];
  in[0]=in_e[0];
  in[1]=in_e[1];
  in[2]=in_e[2];
  in[3]=in_e[3];

  init_S_zeros(j);
      
    
  while(inlen -l >= rsiz){

    join1(in,j,4);
    keccakf(s);
  
    inlen -= g;
    l=8;
    in[0] += g;
    in[1] += g;
    in[2] += g;
    in[3] += g;
	    
    if(inlen -l < rsiz){
      k=2;
      break;
    }

    
    join2(in,j,3);
    
    keccakf(s);
    inlen -= g;
    l=16;
    in[0] += g;
    in[1] += g;
    in[2] += g;
    in[3] += g;
    if(inlen -l < rsiz){
      k=3;
      break;
    }

    join3(in,j,3);

    keccakf(s);
    inlen -= g;
    l=24;
    in[0] += g;
    in[1] += g;
    in[2] += g;
    in[3] += g;
    
    
    
    if(inlen -l < rsiz){
      k=4;
      break;
    }

    join4(in,j,4);

    keccakf(s);
    inlen -= gf;
    l=0;
    in[0] += gf;
    in[1] += gf;
    in[2] += gf;
    in[3] += gf;
    
    
    if(inlen -l < rsiz){
      k=1;
      break;
    }
  }
    

  switch(k){
    case 2:
      inlen-=8;
      in[0] += 8;
      in[1] += 8;
      in[2] += 8;
      in[3] += 8;
      break;
    case 3:
      inlen-=16;
      in[0] += 16;
      in[1] += 16;
      in[2] += 16;
      in[3] += 16;
      break;
    case 4:
      inlen-=24;
      in[0] += 24;
      in[1] += 24;
      in[2] += 24;
      in[3] += 24;
      break;
    case 1:
      break;
    default:
      printf("ERRO!!!");
      return 0;
  }
    
 
  memset(temp0, 0, 144*sizeof(uint8_t));
  memset(temp1, 0, 144*sizeof(uint8_t));
  memset(temp2, 0, 144*sizeof(uint8_t));
  memset(temp3, 0, 144*sizeof(uint8_t));
  
  // last block and padding
  memcpy(temp0, in[0], inlen);
  memcpy(temp1, in[1], inlen);
  memcpy(temp2, in[2], inlen);
  memcpy(temp3, in[3], inlen);
  temp0[inlen] = 0x06;
  temp1[inlen] = 0x06;
  temp2[inlen] = 0x06;
  temp3[inlen++] = 0x06;
  temp0[rsiz - 1] |= 0x80;
  temp1[rsiz - 1] |= 0x80;
  temp2[rsiz - 1] |= 0x80;
  temp3[rsiz - 1] |= 0x80;

  join_last(j,4); 
  
  keccakf(s);
  
//     for(i=0;i<25;i++){
//       print_256(s[i]);
//     }
//      getchar();
//   
  
//   print_256(s[0]);
//   print_256(s[1]);
//   print_256(s[2]);
//   print_256(s[3]);
//   getchar();
  
  
  back(j,1);
  
//   print_256(s[0]);
//   print_256(s[1]);
//   print_256(s[2]);
//   print_256(s[3]);
//   getchar();
  
  store(md,0,s[0]);
  store(md1,0,s[1]);
  store(md2,0,s[2]);
  store(md3,0,s[3]);

  return 0;
}

int keccak_384(char **in_e, int inlen, uint8_t *md,uint8_t *md1,uint8_t *md2,uint8_t *md3)
{
   __m256i aux[10];
   int i;
  const ALIGN uint64_t RC[100] = {
    0x0000000000000001,0x0000000000000001,0x0000000000000001,0x0000000000000001,
    0x0000000000008082,0x0000000000008082,0x0000000000008082,0x0000000000008082,
    0x800000000000808A,0x800000000000808A,0x800000000000808A,0x800000000000808A,  
    0x8000000080008000,0x8000000080008000,0x8000000080008000,0x8000000080008000, 
    0x000000000000808B,0x000000000000808B,0x000000000000808B,0x000000000000808B,   
    0x0000000080000001,0x0000000080000001,0x0000000080000001,0x0000000080000001, 
    0x8000000080008081,0x8000000080008081,0x8000000080008081,0x8000000080008081, 
    0x8000000000008009,0x8000000000008009,0x8000000000008009,0x8000000000008009, 
    0x000000000000008A,0x000000000000008A,0x000000000000008A,0x000000000000008A, 
    0x0000000000000088,0x0000000000000088,0x0000000000000088,0x0000000000000088, 
    0x0000000080008009,0x0000000080008009,0x0000000080008009,0x0000000080008009, 
    0x000000008000000A,0x000000008000000A,0x000000008000000A,0x000000008000000A, 
    0x000000008000808B,0x000000008000808B,0x000000008000808B,0x000000008000808B, 
    0x800000000000008B,0x800000000000008B,0x800000000000008B,0x800000000000008B, 
    0x8000000000008089,0x8000000000008089,0x8000000000008089,0x8000000000008089, 
    0x8000000000008003,0x8000000000008003,0x8000000000008003,0x8000000000008003, 
    0x8000000000008002,0x8000000000008002,0x8000000000008002,0x8000000000008002, 
    0x8000000000000080,0x8000000000000080,0x8000000000000080,0x8000000000000080, 
    0x000000000000800A,0x000000000000800A,0x000000000000800A,0x000000000000800A, 
    0x800000008000000A,0x800000008000000A,0x800000008000000A,0x800000008000000A, 
    0x8000000080008081,0x8000000080008081,0x8000000080008081,0x8000000080008081, 
    0x8000000000008080,0x8000000000008080,0x8000000000008080,0x8000000000008080, 
    0x0000000080000001,0x0000000080000001,0x0000000080000001,0x0000000080000001, 
    0x8000000080008008,0x8000000080008008,0x8000000080008008,0x8000000080008008 
  };
  

  const int rsiz = 104;
  const int g = 96;
  const int gf = 128;
  ALIGN uint8_t temp0[144];
  ALIGN uint8_t temp1[144];
  ALIGN uint8_t temp2[144];
  ALIGN uint8_t temp3[144];
  int k=1,j,l=0;
  __m256i s[28],ain[8];
  __m128i auxl[5];
  char* con[32];
  char* in[4];
  in[0]=in_e[0];
  in[1]=in_e[1];
  in[2]=in_e[2];
  in[3]=in_e[3];

  init_S_zeros(j);
      
    
  while(inlen -l >= rsiz){


    join1(in,j,3);
	
    keccakf(s);
    inlen -= g;
    l=8;
    in[0] += g;
    in[1] += g;
    in[2] += g;
    in[3] += g;
	    
    if(inlen -l < rsiz){
      k=2;
      break;
    }
    join2(in,j,2);

    keccakf(s);
    inlen -= g;
    l=16;
    in[0] += g;
    in[1] += g;
    in[2] += g;
    in[3] += g;
    if(inlen -l < rsiz){
      k=3;
      break;
    }

    join3(in,j,2);
    
    keccakf(s);
    inlen -= g;
    l=24;
    in[0] += g;
    in[1] += g;
    in[2] += g;
    in[3] += g;
    
    if(inlen -l < rsiz){
      k=4;
      break;
    }

    join4(in,j,3);
    
    keccakf(s);
    inlen -= gf;
    l=0;
    in[0] += gf;
    in[1] += gf;
    in[2] += gf;
    in[3] += gf;
    
    
    if(inlen -l < rsiz){
      k=1;
      break;
    }
  }
    

    switch(k){
      case 2:
	inlen-=8;
	in[0] += 8;
	in[1] += 8;
	in[2] += 8;
	in[3] += 8;
	break;
      case 3:
	inlen-=16;
	in[0] += 16;
	in[1] += 16;
	in[2] += 16;
	in[3] += 16;
	break;
      case 4:
	inlen-=24;
	in[0] += 24;
	in[1] += 24;
	in[2] += 24;
	in[3] += 24;
	break;
      case 1:
	break;
      default:
	printf("ERRO!!!");
	return 0;
    }
    
 
    memset(temp0, 0, 144*sizeof(uint8_t));
    memset(temp1, 0, 144*sizeof(uint8_t));
    memset(temp2, 0, 144*sizeof(uint8_t));
    memset(temp3, 0, 144*sizeof(uint8_t));
    
    // last block and padding
    memcpy(temp0, in[0], inlen);
    memcpy(temp1, in[1], inlen);
    memcpy(temp2, in[2], inlen);
    memcpy(temp3, in[3], inlen);
    temp0[inlen] = 0x06;
    temp1[inlen] = 0x06;
    temp2[inlen] = 0x06;
    temp3[inlen++] = 0x06;
    temp0[rsiz - 1] |= 0x80;
    temp1[rsiz - 1] |= 0x80;
    temp2[rsiz - 1] |= 0x80;
    temp3[rsiz - 1] |= 0x80;
    

    join_last(j,3); 
    
    keccakf(s);
    

    back(j,2);
    store(md,0,s[0]);
    store(md,1,s[4]);
    store(md1,0,s[1]);
    store(md1,1,s[5]);
    store(md2,0,s[2]);
    store(md2,1,s[6]);
    store(md3,0,s[3]);
    store(md3,1,s[7]);

    return 0;
}

int keccak_512	(char **in_e, int inlen, uint8_t *md,uint8_t *md1,uint8_t *md2,uint8_t *md3)
{
   __m256i aux[10];
   int i;
  const ALIGN uint64_t RC[100] = {
    0x0000000000000001,0x0000000000000001,0x0000000000000001,0x0000000000000001,
    0x0000000000008082,0x0000000000008082,0x0000000000008082,0x0000000000008082,
    0x800000000000808A,0x800000000000808A,0x800000000000808A,0x800000000000808A,  
    0x8000000080008000,0x8000000080008000,0x8000000080008000,0x8000000080008000, 
    0x000000000000808B,0x000000000000808B,0x000000000000808B,0x000000000000808B,   
    0x0000000080000001,0x0000000080000001,0x0000000080000001,0x0000000080000001, 
    0x8000000080008081,0x8000000080008081,0x8000000080008081,0x8000000080008081, 
    0x8000000000008009,0x8000000000008009,0x8000000000008009,0x8000000000008009, 
    0x000000000000008A,0x000000000000008A,0x000000000000008A,0x000000000000008A, 
    0x0000000000000088,0x0000000000000088,0x0000000000000088,0x0000000000000088, 
    0x0000000080008009,0x0000000080008009,0x0000000080008009,0x0000000080008009, 
    0x000000008000000A,0x000000008000000A,0x000000008000000A,0x000000008000000A, 
    0x000000008000808B,0x000000008000808B,0x000000008000808B,0x000000008000808B, 
    0x800000000000008B,0x800000000000008B,0x800000000000008B,0x800000000000008B, 
    0x8000000000008089,0x8000000000008089,0x8000000000008089,0x8000000000008089, 
    0x8000000000008003,0x8000000000008003,0x8000000000008003,0x8000000000008003, 
    0x8000000000008002,0x8000000000008002,0x8000000000008002,0x8000000000008002, 
    0x8000000000000080,0x8000000000000080,0x8000000000000080,0x8000000000000080, 
    0x000000000000800A,0x000000000000800A,0x000000000000800A,0x000000000000800A, 
    0x800000008000000A,0x800000008000000A,0x800000008000000A,0x800000008000000A, 
    0x8000000080008081,0x8000000080008081,0x8000000080008081,0x8000000080008081, 
    0x8000000000008080,0x8000000000008080,0x8000000000008080,0x8000000000008080, 
    0x0000000080000001,0x0000000080000001,0x0000000080000001,0x0000000080000001, 
    0x8000000080008008,0x8000000080008008,0x8000000080008008,0x8000000080008008 
  };
  
  const int rsiz = 72;
  const int g = 64;
  const int gf = 96;
  ALIGN uint8_t temp0[144];
  ALIGN uint8_t temp1[144];
  ALIGN uint8_t temp2[144];
  ALIGN uint8_t temp3[144];
  int k=1,j,l=0;
  __m256i s[28],ain[8];
  __m128i auxl[5];
  char* con[32];
  char* in[4];
  in[0]=in_e[0];
  in[1]=in_e[1];
  in[2]=in_e[2];
  in[3]=in_e[3];

  init_S_zeros(j);
    
    
  while(inlen -l >= rsiz){
    
    join1(in,j,2);
          
    keccakf(s);
    inlen -= g;
    l=8;
    in[0] += g;
    in[1] += g;
    in[2] += g;
    in[3] += g;
    
    
	    
    if(inlen -l < rsiz){
      k=2;
      break;
    }

    join2(in,j,1);
    
    keccakf(s);
    inlen -= g;
    l=16;
    in[0] += g;
    in[1] += g;
    in[2] += g;
    in[3] += g;
    if(inlen -l < rsiz){
      k=3;
      break;
    }

    join3(in,j,1);
    
    keccakf(s);
    inlen -= g;
    l=24;
    in[0] += g;
    in[1] += g;
    in[2] += g;
    in[3] += g;
    
    if(inlen -l < rsiz){
      k=4;
      break;
    }

    join4(in,j,2);
    
    keccakf(s);
    inlen -= gf;
    l=0;
    in[0] += gf;
    in[1] += gf;
    in[2] += gf;
    in[3] += gf;
    
    
    if(inlen -l < rsiz){
      k=1;
      break;
    }
  }
  

  switch(k){
    case 2:
      inlen-=8;
      in[0] += 8;
      in[1] += 8;
      in[2] += 8;
      in[3] += 8;
      break;
    case 3:
      inlen-=16;
      in[0] += 16;
      in[1] += 16;
      in[2] += 16;
      in[3] += 16;
      break;
    case 4:
      inlen-=24;
      in[0] += 24;
      in[1] += 24;
      in[2] += 24;
      in[3] += 24;
      break;
    case 1:
      break;
    default:
      printf("ERRO!!!");
      return 0;
  }
  

  memset(temp0, 0, 144*sizeof(uint8_t));
  memset(temp1, 0, 144*sizeof(uint8_t));
  memset(temp2, 0, 144*sizeof(uint8_t));
  memset(temp3, 0, 144*sizeof(uint8_t));
  
  // last block and padding
  memcpy(temp0, in[0], inlen);
  memcpy(temp1, in[1], inlen);
  memcpy(temp2, in[2], inlen);
  memcpy(temp3, in[3], inlen);
  temp0[inlen] = 0x06;
  temp1[inlen] = 0x06;
  temp2[inlen] = 0x06;
  temp3[inlen++] = 0x06;
  temp0[rsiz - 1] |= 0x80;
  temp1[rsiz - 1] |= 0x80;
  temp2[rsiz - 1] |= 0x80;
  temp3[rsiz - 1] |= 0x80;
  
  
  join_last(j,2); 
    
  keccakf(s);
  
  back(j,2);
  store(md,0,s[0]);
  store(md,1,s[4]);
  store(md1,0,s[1]);
  store(md1,1,s[5]);
  store(md2,0,s[2]);
  store(md2,1,s[6]);
  store(md3,0,s[3]);
  store(md3,1,s[7]);

  return 0;
}


int openmp_256(char ** in, int inlen, uint8_t **md,uint8_t **md1,uint8_t **md2,uint8_t **md3, int t)
{
  int i;
  #pragma omp parallel for num_threads(t)
     for(i=0;i<t;i++){
        keccak_256(in, inlen,md[i],md1[i],md2[i],md3[i]);
     }
}

int openmp_384(char ** in, int inlen, uint8_t **md,uint8_t **md1,uint8_t **md2,uint8_t **md3)
{
  int i;
  #pragma omp parallel for num_threads(4)
     for(i=0;i<4;i++){
        keccak_384(in, inlen,md[i],md1[i],md2[i],md3[i]);
     }
}


int openmp_512(char ** in, int inlen, uint8_t **md,uint8_t **md1,uint8_t **md2,uint8_t **md3)
{
  int i;
  #pragma omp parallel for num_threads(4)
     for(i=0;i<4;i++){
        keccak_512(in, inlen,md[i],md1[i],md2[i],md3[i]);
     }
}