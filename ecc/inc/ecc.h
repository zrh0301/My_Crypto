#include "../bn/inc/wlc_bn.h"

#define ECC_LEN 8
#define WORD_BITS 32
#define ECC_BITS  256

typedef struct{
	dig_t x[ECC_LEN];
	dig_t y[ECC_LEN];
} ecc_point;

typedef struct{
	dig_t p[ECC_LEN];
	dig_t n[ECC_LEN];
	dig_t a[ECC_LEN];
	dig_t b[ECC_LEN];
} ecc_param;

void ecc_copy(ecc_point*& r,const ecc_point* p);
void ecc_make_point(ecc_point*& r,const dig_t* x,const dig_t* y,const ecc_param* param);
void ecc_add(ecc_point*& r,const ecc_point* p,const ecc_point* q,const ecc_param* param);
void ecc_double(ecc_point*& r,const ecc_point* p,const ecc_param* param);
void ecc_mul(ecc_point*& pk, const dig_t* sk, const ecc_param* param);