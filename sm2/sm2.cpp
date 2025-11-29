#include "../bn/inc/wlc_bn.h"
#include "../ecc/inc/ecc.h"

//dig_t 32bits
//p[0]->lowest 32 bits
const dig_t p[8] = {
	0xffffffff,0xffffffff,0x00000000,0xffffffff,
	0xffffffff,0xffffffff,0xffffffff,0xfffffffe
};
const dig_t a[8] = {
	0xfffffffc,0xffffffff,0x00000000,0xffffffff,
	0xffffffff,0xffffffff,0xffffffff,0xfffffffe
};
const dig_t b[8] = {
	0x4d940e93,0xddbcbd41,0x15ab8f92,0xf39789f5,
	0xcf6509a7,0x4d5a9e4b,0x9d9f5e34,0x28e9fa9e
};
const dig_t n[8] = {
	0x39d54123,0x53bbf409,0x21c6052b,0x7203df6b,
	0xffffffff,0xffffffff,0xffffffff,0xfffffffe
};
const dig_t g_x[8] = {
	0x334c747c,0x715a4859,0xf2660be1,0x8fe30bbf,
	0x6a39c994,0x5f990446,0x1f198119,0x32c4ae2c
};
const dig_t g_y[8] = {
	0x2139f0a0,0x02df32e5,0xc62a4740,0xd0a9877c,
	0x6b692153,0x59bdcee3,0xf4f6779c,0xbc3736a2
};

void init(ecc_param*& param,ecc_point*& G){
	bn_copy(param->a,a,ECC_LEN);
	bn_copy(param->b,b,ECC_LEN);
	bn_copy(param->n,n,ECC_LEN);
	bn_copy(param->p,p,ECC_LEN);
	bn_copy(G->x,g_x,ECC_LEN);
	bn_copy(G->y,g_y,ECC_LEN);
}

void key_gen(dig_t*& sk,ecc_point*& pk,const ecc_point* G,const ecc_param* param){
	do{
		bn_rand(sk,8);
	}while(bn_cmp(sk,param->n,ECC_LEN) == 1 or bn_is_zero(sk,ECC_LEN));
	ecc_copy(pk,G);
	ecc_mul(pk,sk,param);
}

int main(){
	ecc_param* param = (ecc_param*)calloc(1,sizeof(ecc_param));
	ecc_point* G = (ecc_point*)calloc(1,sizeof(ecc_point));
	ecc_point* pk = (ecc_point*)calloc(1,sizeof(ecc_point));
	dig_t* sk = (dig_t*)calloc(8,sizeof(dig_t));
	init(param,G);
	key_gen(sk,pk,G,param);
	bn_print(sk,ECC_LEN);
	bn_print(pk->x,ECC_LEN);
	bn_print(pk->y,ECC_LEN);
	return 0;
}