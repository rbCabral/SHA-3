#include "keccak.h"

#define ROT(x,y)  _mm_or_si128(_mm_slli_epi64(x,y),_mm_srli_epi64(x,(64-y)))
#define ROT1(x)   _mm_or_si128(_mm_add_epi64(x,x),_mm_srli_epi64(x,(63)))

#define load(x,m,d) 	  	x = _mm_loadu_si128((__m128i*)m+d)
#define store(x,d,y)		_mm_storeu_si128((__m128i*)x+d,y);
#define unpackLO_64(x,y,z)	x = _mm_unpacklo_epi64(y,z)
#define unpackHI_64(x,y,z)	x = _mm_unpackhi_epi64(y,z)
#define set_zero_128(x)		x = _mm_setzero_si128()
#define xor(x,y,z) 	  	x = _mm_xor_si128(y,z)
#define andnot(y,z)		_mm_andnot_si128(y,z)


#define blend_32(x,y,z,d) 	x = _mm256_blend_epi32(y,z,d)
#define permut_64(x,y,d)  	x = _mm256_permute4x64_epi64(y,d)
#define xor_pd(x,y,z) 	  	x = _mm256_xor_pd(y,z)


#define extract_128(x,y,d)	x = _mm256_extracti128_si256(y,d)
#define insert_128(x,y,z,d)	x = _mm256_inserti128_si256(y,z,d)


#define load_128(x,m,d) 	x = _mm_load_si128((__m128i*)m+d)
#define store_128(x,d,y)	_mm_store_si128((__m128i*)x+d,y);
#define unpackLO_128_64(x,y,z)	x = _mm_unpacklo_epi64(y,z)
#define unpackHI_128_64(x,y,z)	x = _mm_unpackhi_epi64(y,z)
#define cast_to_256(x,y)	x = _mm256_castsi128_si256(y)

#define init_S_zeros(y)		\
							for(y=0;y<25;y++)\
				  				set_zero_128(s[y]);

#define loadM(j,p,M)	\
				for(j=0;j<p;j++){\
					load(t[1],M[0],j);\
				  	load(t[2],M[1],j);\
				  	\
				  	unpackLO_64(t[0],t[1],t[2]);\
				  	unpackHI_64(t[1],t[1],t[2]);\
				  	\
				  	xor(s[0+(j*2)],s[0+(j*2)],t[0]);\
				  	xor(s[1+(j*2)],s[1+(j*2)],t[1]);\
				}\
				if(rsiz != 144){\
					load(t[0],M[0],j);\
					load(t[1],M[1],j);\
					\
					unpackLO_64(t[0],t[0],t[1]);\
					\
					xor(s[0+(j*2)],s[0+(j*2)],t[0]);\
				}\

#define back(j,p,M)	\
				for(j=0;j<p;j++){\
				  	unpackLO_64(t[0],s[0+(j*2)],s[1+(j*2)]);\
				  	unpackHI_64(t[1],s[0+(j*2)],s[1+(j*2)]);\
				  	\
				  	store(M[0],j,t[0]);\
				  	store(M[1],j,t[1]);\
				}\
				if(rsiz != 144){\
					unpackLO_64(t[0],s[0+(j*2)],s[1+(j*2)]);\
					unpackHI_64(t[1],s[0+(j*2)],s[1+(j*2)]);\
					store(M[0],j,t[0]);\
				  	store(M[1],j,t[1]);\
				}\

#if defined (AVX2)

inline void keccakF(__m128i *s){\
	int i;
	__m128i aux[10];
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
	for(i=0;i<24;i=i+2){\
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
		xor(s[5],s[5],aux[5]);\
		xor(s[10],s[10],aux[5]);\
		xor(s[15],s[15],aux[5]);\
		xor(s[20],s[20],aux[5]);\
		\
		xor(s[11],s[11],aux[6]);\
		xor(s[16],s[16],aux[6]);\
		xor(s[21],s[21],aux[6]);\
		\
		xor(s[7],s[7],aux[7]);\
		xor(s[17],s[17],aux[7]);\
		xor(s[22],s[22],aux[7]);\
		\
		xor(s[8],s[8],aux[8]);\
		xor(s[13],s[13],aux[8]);\
		xor(s[23],s[23],aux[8]);\
		\
		xor(s[9],s[9],aux[9]);\
		xor(s[14],s[14],aux[9]);\
		xor(s[19],s[19],aux[9]);\
		\
		\
		xor(s[6],s[6],aux[6]);\
		xor(s[12],s[12],aux[7]);\
		xor(s[18],s[18],aux[8]);\
		xor(s[24],s[24],aux[9]);\
		xor(s[0],s[0],aux[5]);\
		xor(s[1],s[1],aux[6]);\
		xor(s[2],s[2],aux[7]);\
		xor(s[3],s[3],aux[8]);\
		xor(s[4],s[4],aux[9]);\
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
	}
#endif

#if defined (AVX512)
		void keccakF(__m128i *s){   
		int i;   
			__m128i aux[10];
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
                    }
#endif

void print_128(__m128i x){
     uint64_t c[2];
     store(c,0,x);
     printf("\n 1 = %8.16lx\t  0 = %8.16lx\t\n", c[1],c[0]);
}

int keccak(char **in_e, int inlen, uint8_t **md, int ra, int out)
{
	

  ALIGN uint8_t temP[2][200];
    
  int j,i;
  __m128i s[25],t[3];
  char* in[2];
  int r;
  
  in[0]=in_e[0];
  in[1]=in_e[1];

  init_S_zeros(j);
  int rsiz = 200-2*ra;
  r = (rsiz/8)/2;

  for(;inlen >=rsiz; inlen -=rsiz){
          loadM(j,r,in);
          keccakF(s);
          in[0] += rsiz;
          in[1] += rsiz;
  }

  for(i=0;i<2;i++){
    memset(temP[i], 0, 168*sizeof(uint8_t));
    memcpy(temP[i], in[i], inlen);
    temP[i][inlen] = 0x06;
    temP[i][rsiz - 1] |= 0x80;
  }
  loadM(j,r,temP);
  keccakF(s);

  back(j,r,temP);
  if(rsiz < out){
  	memcpy(md[0], temP[0], rsiz);
  	memcpy(md[1], temP[1], rsiz);
  	out -= rsiz;
  	back(j,r,temP);
  	keccakF(s);
  }
  memcpy(md[0], temP[0], out);
  memcpy(md[1], temP[1], out);
  
  // printf("\n");
  // for(i=0;i<32;i++){
  // 	printf("%x", md[0][i]);
  // }
  // printf("\n");
  // for(i=0;i<32;i++){
  // 	printf("%x", md[1][i]);
  // }
  // printf("\n");
  

  return 0;
}
