#include "ecc.h"

void ecc_copy(ecc_point*& r,const ecc_point* p){
	memcpy(r,p,sizeof(ecc_point));
}

void ecc_make_point(ecc_point*& r,const dig_t* x,const dig_t* y,const ecc_param* param){
	bn_copy(r->x,x,ECC_LEN);
	bn_copy(r->y,y,ECC_LEN);
}

void ecc_add(ecc_point*& r,const ecc_point* p,const ecc_point* q,const ecc_param* param){
	
	if(bn_is_zero(p->x,ECC_LEN) and bn_is_zero(p->y,ECC_LEN)){
		ecc_copy(r,q);
		return ;
	}
	if(bn_is_zero(q->x,ECC_LEN) and bn_is_zero(q->y,ECC_LEN)){
		ecc_copy(r,p);
		return ;
	}

	dig_t* t1 = (dig_t*)calloc(8,sizeof(dig_t));//t1=x1-x2
	dig_t* t2 = (dig_t*)calloc(8,sizeof(dig_t));//t2=y1-y2

	bn_mod_sub(t1,p->x,q->x,param->p,ECC_LEN);
	bn_mod_sub(t2,p->y,q->y,param->p,ECC_LEN);

	if(bn_is_zero(t1,ECC_LEN)){
		if(bn_is_zero(t2,ECC_LEN)){
			ecc_double(r,p,param);
		}
		else{
			memset(r,0,sizeof(ecc_point));
		}
		return ;
	}

	dig_t* t1_inv = (dig_t*)calloc(8,sizeof(dig_t));//t1^{-1}
	bn_mod_inv(t1_inv,t1,param->p,ECC_LEN);
	dig_t* lambda = (dig_t*)calloc(8,sizeof(dig_t));
	bn_mod_mul(lambda,t1_inv,t2,param->p,ECC_LEN);

	dig_t* x3 = (dig_t*)calloc(8,sizeof(dig_t));
	dig_t* y3 = (dig_t*)calloc(8,sizeof(dig_t));

	bn_mod_mul(x3,lambda,lambda,param->p,ECC_LEN);
	bn_mod_sub(x3,x3,p->x,param->p,ECC_LEN);
	bn_mod_sub(x3,x3,q->x,param->p,ECC_LEN);
	
	dig_t* t3 = (dig_t*)calloc(8,sizeof(dig_t));//t3=x1-x3
	bn_mod_sub(t3,p->x,x3,param->p,ECC_LEN);
	bn_mod_mul(y3,lambda,t3,param->p,ECC_LEN);
	bn_mod_sub(y3,y3,p->y,param->p,ECC_LEN);

	ecc_make_point(r,x3,y3,param);

	free(t1);free(t2);free(t1_inv);free(lambda);
	free(x3);free(y3);free(t3);
}

void ecc_double(ecc_point*& r, const ecc_point* p, const ecc_param* param) {
    
    if (bn_is_zero(p->y, ECC_LEN) or 
    	bn_is_zero(p->x, ECC_LEN)){
        memset(r, 0, sizeof(ecc_point));
        return;
    }

    dig_t* t1 = (dig_t*)calloc(ECC_LEN, sizeof(dig_t)); // x^2, 3x^2
    dig_t* t2 = (dig_t*)calloc(ECC_LEN, sizeof(dig_t)); // 2y, (2y)^-1
    dig_t* lambda = (dig_t*)calloc(ECC_LEN, sizeof(dig_t)); // 斜率
    dig_t* x3 = (dig_t*)calloc(ECC_LEN, sizeof(dig_t));
    dig_t* y3 = (dig_t*)calloc(ECC_LEN, sizeof(dig_t));
    
    
    // t1 = x^2
    bn_mod_mul(t1, p->x, p->x, param->p, ECC_LEN);

    bn_mod_add(lambda, t1, t1, param->p, ECC_LEN); // 2x^2
    bn_mod_add(lambda, lambda, t1, param->p, ECC_LEN); // 3x^2
    
    // lambda_fenzi = 3x^2 + a
    bn_mod_add(lambda, lambda, param->a, param->p, ECC_LEN);

    // t2 = y + y
    bn_mod_add(t2, p->y, p->y, param->p, ECC_LEN);

    // t1 = (2y)^-1
    bn_mod_inv(t1, t2, param->p, ECC_LEN);

    // lambda_final
    bn_mod_mul(lambda, lambda, t1, param->p, ECC_LEN);

    // x3 = lambda^2 - 2x1 
    bn_mod_mul(x3, lambda, lambda, param->p, ECC_LEN); // x3 = lambda^2
    bn_mod_sub(x3, x3, p->x, param->p, ECC_LEN);
    bn_mod_sub(x3, x3, p->x, param->p, ECC_LEN);       // - x1 (即 -2x1)
    

    // y3 = lambda(x1 - x3) - y1 
    // t1 = (x1 - x3)
    bn_mod_sub(t1, p->x, x3, param->p, ECC_LEN);       
    bn_mod_mul(y3, t1, lambda, param->p, ECC_LEN);     // y3 = lambda(x1 - x3)
    bn_mod_sub(y3, y3, p->y, param->p, ECC_LEN);       // - y1

    bn_copy(r->x, x3, ECC_LEN);
    bn_copy(r->y, y3, ECC_LEN);
    
    free(t1); free(t2); free(lambda); free(x3); free(y3);
}

int get_bit(const dig_t* sk,int index) {
    int word_idx = index/WORD_BITS;
    int bit_idx  = index%WORD_BITS;
    
    return (sk[word_idx]>>bit_idx)&1;
}

void ecc_mul(ecc_point*& r,const ecc_point* g, const dig_t* sk, const ecc_param* param){
    
    ecc_point* R = (ecc_point*)calloc(1,sizeof(ecc_point)); // 累加器 (Result)
    ecc_point* P = (ecc_point*)calloc(1,sizeof(ecc_point)); // 基点副本 (Base Point)
    memcpy(P,g,sizeof(ecc_point));

    
    int i = ECC_BITS-1;
    for(;i>=0;--i){
        if(get_bit(sk,i)) break;
    }

    if(i>=0){
        // memcpy(R, P, sizeof(ecc_point)); // R = P
        --i; 
        ecc_copy(R,P);
    } 
    else{
        memset(r, 0, sizeof(ecc_point));
        free(R); free(P);
        return;
    }

    for(;i>=0;--i) {
        // Double: R = 2 * R
        ecc_double(R,R,param);

        // Add: if get_bit == 1，then R = R + P
        if (get_bit(sk,i)){
            ecc_add(R,R,P,param);
        }
    }

    ecc_copy(r,R);

    free(R);free(P);
}