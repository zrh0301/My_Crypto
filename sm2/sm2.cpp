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

void key_gen(dig_t*& sk,ecc_point*& pk,
	    	 const ecc_point* G,const ecc_param* param){
	do{
		bn_rand(sk,8);
	}while(bn_cmp(sk,param->n,ECC_LEN)==1 or 
		   bn_is_zero(sk,ECC_LEN));
	// ecc_copy(pk,G);
	ecc_mul(pk,G,sk,param);
}
/*
void ecc_sign(dig_t*& sig_r,dig_t*& sig_s,
			  dig_t* sk,const ecc_point* pk,
			  const dig_t* e,
			  const ecc_point* G,const ecc_param* param){
	dig_t* k = (dig_t*)calloc(8,sizeof(dig_t));
	do{
		bn_rand(k,8);
	}while(bn_cmp(k,param->n,ECC_LEN)==1 or 
		   bn_is_zero(k,ECC_LEN));

	ecc_point* g1 = (ecc_point*)calloc(1,sizeof(ecc_point));
	ecc_mul(g1,G,k,param);
	dig_t* x1 = (dig_t*)calloc(8,sizeof(dig_t));
	dig_t* y1 = (dig_t*)calloc(8,sizeof(dig_t));
	bn_copy(x1,g1->x,ECC_LEN);
	bn_copy(y1,g1->y,ECC_LEN);

	bn_mod_add(sig_r,e,x1,param->n,ECC_LEN);
	printf("sig_r:");
	bn_print(sig_r,ECC_LEN);	

	dig_t* chk = (dig_t*)calloc(8,sizeof(dig_t));
	bn_mod_add(chk,sig_r,k,param->n,ECC_LEN);
	if(bn_is_zero(chk,ECC_LEN) or bn_is_zero(sig_r,ECC_LEN)) 
		ecc_sign(sig_r,sig_s,sk,pk,e,G,param);

	dig_t* tmp1 = (dig_t*)calloc(8,sizeof(dig_t));

	bn_mod_mul(tmp1,sig_r,sk,param->n,ECC_LEN);
	bn_mod_sub(tmp1,k,tmp1,param->n,ECC_LEN);
	dig_t* tmp2 = (dig_t*)calloc(8,sizeof(dig_t));
	tmp2[0] = 1;
	for(int i=1;i<=7;++i) tmp2[i] = 0;
	bn_mod_add(tmp2,tmp2,sk,param->n,ECC_LEN);

	dig_t* sk_inv = (dig_t*)calloc(8,sizeof(dig_t));
	bn_mod_inv(sk_inv,tmp2,param->n,ECC_LEN);
	bn_mod_mul(sig_s,tmp1,sk_inv,param->n,ECC_LEN);
	printf("sig_s:");
	bn_print(sig_s,ECC_LEN);
}

void ecc_verify(const dig_t* r,const dig_t* s,
				const ecc_point* pk,
				const dig_t* e,
				const ecc_point* G,const ecc_param* param){
	dig_t* tmp = (dig_t*)calloc(8,sizeof(dig_t));
	bn_mod_add(tmp,r,s,param->n,ECC_LEN);

	ecc_point* tmp1 = (ecc_point*)calloc(1,sizeof(ecc_point));
	ecc_mul(tmp1,G,s,param);
	ecc_point* tmp2 = (ecc_point*)calloc(1,sizeof(ecc_point));
	ecc_mul(tmp2,pk,tmp,param);
	ecc_add(tmp2,tmp1,tmp2,param);
	bn_mod_add(tmp,e,tmp2->x,param->n,ECC_LEN);
	bn_print(tmp,ECC_LEN);
}
*/
// 辅助宏，方便释放
#define FREE_ALL(...) do { \
    void* ptrs[] = {__VA_ARGS__}; \
    for(int i=0; i<sizeof(ptrs)/sizeof(void*); i++) free(ptrs[i]); \
} while(0)

void ecc_sign(dig_t* sig_r, dig_t* sig_s, // 去掉引用 &
              const dig_t* sk, const ecc_point* pk,
              const dig_t* e,
              const ecc_point* G, const ecc_param* param) {
    
    dig_t* k = (dig_t*)calloc(ECC_LEN, sizeof(dig_t));
    ecc_point* g1 = (ecc_point*)calloc(1, sizeof(ecc_point));
    dig_t* r_val = (dig_t*)calloc(ECC_LEN, sizeof(dig_t)); // 内部用的 r
    dig_t* chk = (dig_t*)calloc(ECC_LEN, sizeof(dig_t));
    dig_t* tmp1 = (dig_t*)calloc(ECC_LEN, sizeof(dig_t));
    dig_t* tmp2 = (dig_t*)calloc(ECC_LEN, sizeof(dig_t));
    dig_t* sk_inv = (dig_t*)calloc(ECC_LEN, sizeof(dig_t));
    dig_t* one = (dig_t*)calloc(ECC_LEN, sizeof(dig_t));
    one[0] = 1;

    while (1) { // 使用循环代替递归
        // 1. 产生随机数 k in [1, n-1]
        do {
            bn_rand(k, ECC_LEN);
        } while (bn_cmp(k, param->n, ECC_LEN) >= 0 || bn_is_zero(k, ECC_LEN));

        // 2. 计算点 (x1, y1) = k * G
        ecc_mul(g1, G, k, param);

        // 3. r = (e + x1) mod n
        // 注意：虽然 x1 是 mod p，但这里要 mod n
        bn_mod_add(r_val, e, g1->x, param->n, ECC_LEN);

        // 4. 检查 r != 0 且 r + k != n
        if (bn_is_zero(r_val, ECC_LEN)) continue;

        bn_mod_add(chk, r_val, k, param->n, ECC_LEN);
        if (bn_is_zero(chk, ECC_LEN)) continue;

        // 5. s = (1 + d)^-1 * (k - r*d) mod n
        
        // tmp2 = 1 + d
        bn_mod_add(tmp2, one, sk, param->n, ECC_LEN);
        
        // sk_inv = (1 + d)^-1
        // 注意：如果 (1+d) 是 0 (即 d = n-1)，则求逆失败。
        // 虽然私钥通常不会是 n-1，但健壮的代码应该处理
        if (bn_is_zero(tmp2, ECC_LEN)) continue; 
        bn_mod_inv(sk_inv, tmp2, param->n, ECC_LEN);

        // tmp1 = r * d
        bn_mod_mul(tmp1, r_val, sk, param->n, ECC_LEN);
        // tmp1 = k - r*d
        bn_mod_sub(tmp1, k, tmp1, param->n, ECC_LEN);

        // s = sk_inv * tmp1
        bn_mod_mul(sig_s, sk_inv, tmp1, param->n, ECC_LEN);

        // 6. 检查 s != 0
        if (bn_is_zero(sig_s, ECC_LEN)) continue;

        // 成功，复制 r 并跳出循环
        bn_copy(sig_r, r_val, ECC_LEN);
        break;
    }

    // 打印调试
    printf("g->x:"); bn_print(g1->x,ECC_LEN);
    printf("g->y:"); bn_print(g1->y,ECC_LEN);
    printf("sig_r:"); bn_print(sig_r, ECC_LEN);
    printf("sig_s:"); bn_print(sig_s, ECC_LEN);

    // 释放内存
    FREE_ALL(k, g1, r_val, chk, tmp1, tmp2, sk_inv, one);
}

int ecc_verify(const dig_t* r, const dig_t* s,
               const ecc_point* pk,
               const dig_t* e,
               const ecc_point* G, const ecc_param* param) {
    
    // 检查 r, s 范围 [1, n-1]
    if (bn_is_zero(r, ECC_LEN) || bn_cmp(r, param->n, ECC_LEN) >= 0) return 0;
    if (bn_is_zero(s, ECC_LEN) || bn_cmp(s, param->n, ECC_LEN) >= 0) return 0;

    dig_t* t = (dig_t*)calloc(ECC_LEN, sizeof(dig_t));
    ecc_point* sg = (ecc_point*)calloc(1, sizeof(ecc_point));
    ecc_point* tpk = (ecc_point*)calloc(1, sizeof(ecc_point));
    ecc_point* p = (ecc_point*)calloc(1, sizeof(ecc_point));
    dig_t* R = (dig_t*)calloc(ECC_LEN, sizeof(dig_t));
    int result = 0;

    // t = (r + s) mod n
    bn_mod_add(t, r, s, param->n, ECC_LEN);

    if (bn_is_zero(t, ECC_LEN)) {
        goto cleanup; // t=0 验证失败
    }

    // P = s*G + t*PA
    // 1. sg = s*G
    ecc_mul(sg, G, s, param);
    
    // 2. tpk = t*PA
    ecc_mul(tpk, pk, t, param);

    // 3. p = sg + tpk
    ecc_add(p, sg, tpk, param);

    // R = (e + x1) mod n
    // x1 取模 n 比较安全，虽然直接加也行(依赖 mod_add 实现)
    bn_mod_add(R, e, p->x, param->n, ECC_LEN);

    printf("p->x:"); bn_print(p->x,ECC_LEN);
    printf("p->y:");bn_print(p->y,ECC_LEN);
    printf("calc_R:"); bn_print(R, ECC_LEN);

    // 比较 R 和 r
    if (bn_cmp(R, r, ECC_LEN) == 0) {
        result = 1; // 成功
    }

cleanup:
    FREE_ALL(t, sg, tpk, p, R);
    return result;
}

int main(){
	ecc_param* param = (ecc_param*)calloc(1,sizeof(ecc_param));
	ecc_point* G = (ecc_point*)calloc(1,sizeof(ecc_point));
	ecc_point* pk = (ecc_point*)calloc(1,sizeof(ecc_point));
	dig_t* sk = (dig_t*)calloc(8,sizeof(dig_t));
	
	init(param,G);
	printf("a:"); bn_print(param->a,ECC_LEN);
	printf("b:"); bn_print(param->b,ECC_LEN);
	printf("p:"); bn_print(param->p,ECC_LEN);
	
	key_gen(sk,pk,G,param);
	
	printf("sk:"); bn_print(sk,ECC_LEN);
	printf("pk->x"); bn_print(pk->x,ECC_LEN);
	printf("pk->y"); bn_print(pk->y,ECC_LEN);
	
	/*
	ecc_double(pk,G,param);
	bn_print(pk->x,ECC_LEN);
	bn_print(pk->y,ECC_LEN);
	*/
	
	dig_t* e = (dig_t*)calloc(8,sizeof(dig_t));//e=sm3(Z,m)
	bn_rand(e,ECC_LEN);
	printf("m:"); bn_print(e,ECC_LEN);
	dig_t* sig_r = (dig_t*)calloc(8,sizeof(dig_t));
	dig_t* sig_s = (dig_t*)calloc(8,sizeof(dig_t));
	ecc_sign(sig_r,sig_s,sk,pk,e,G,param);

	ecc_verify(sig_r,sig_s,pk,e,G,param);
	
	return 0;
}