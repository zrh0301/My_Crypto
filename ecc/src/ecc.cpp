#include "ecc.h"

void ecc_copy(ecc_point*& r,const ecc_point* p){
	memcpy(r,p,sizeof(ecc_point));
}

void ecc_make_point(ecc_point*& r,const dig_t* x,const dig_t* y,const ecc_param* param){
	bn_copy(r->x,x,ECC_LEN);
	bn_copy(r->y,y,ECC_LEN);
}

void ecc_add(ecc_point*& r,const ecc_point* p,const ecc_point* q,const ecc_param* param){
	//calc r = p+q
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
}

void ecc_double(ecc_point*& r,const ecc_point* p,const ecc_param* param){
	if(bn_is_zero(p->y,ECC_LEN)){
		memset(r,0,sizeof(ecc_point));
		return ;
	}

	dig_t* t1 = (dig_t*)calloc(8,sizeof(dig_t));//3x^2+a
	dig_t* t2 = (dig_t*)calloc(8,sizeof(dig_t));//2y
	dig_t* lambda = (dig_t*)calloc(8,sizeof(dig_t));

	bn_mod_mul(t1,p->x,p->x,param->p,ECC_LEN);
	for(int i=1;i<=3;++i) bn_mod_add(lambda,t1,lambda,param->p,ECC_LEN);
	bn_mod_add(lambda,lambda,param->a,param->p,ECC_LEN);
	
	dig_t* t2_inv = (dig_t*)calloc(8,sizeof(dig_t));
	bn_mod_inv(t2_inv,t2,param->p,ECC_LEN);
	bn_mod_add(t2_inv,t2_inv,t2_inv,param->p,ECC_LEN);

	bn_mod_mul(lambda,lambda,t2_inv,param->p,ECC_LEN);

	dig_t* x3 = (dig_t*)calloc(8,sizeof(dig_t));//x3=lambda^2-2x1
	dig_t* y3 = (dig_t*)calloc(8,sizeof(dig_t));//y3=lambda(x1-x3)-y1

	bn_mod_mul(x3,lambda,lambda,param->p,ECC_LEN);
	bn_mod_sub(x3,x3,p->x,param->p,ECC_LEN);
	bn_mod_sub(x3,x3,p->x,param->p,ECC_LEN);
	dig_t* t3 = (dig_t*)calloc(9,sizeof(dig_t));//t3=x1-x3
	bn_mod_sub(t3,p->x,x3,param->p,ECC_LEN);
	bn_mod_mul(y3,t3,lambda,param->p,ECC_LEN);
	bn_mod_sub(y3,y3,p->y,param->p,ECC_LEN);
	ecc_make_point(r,x3,y3,param);
}

// 检查点是否为无穷远点 (零点)
int ecc_is_infinity(const ecc_point* P) {
    // 简单判断：如果 x 和 y 都是 0，视为无穷远点（具体视初始化而定）
    // 需要配合大数库的 bn_is_zero 使用
    // return bn_is_zero(P->x) && bn_is_zero(P->y);
    return 0; // 占位
}

// ==========================================
// 获取大数 sk 第 index 位的比特值 (0 或 1)
// ==========================================
int get_bit(const dig_t* sk, int index) {
    int word_idx = index / WORD_BITS; // 第几个 dig_t
    int bit_idx  = index % WORD_BITS; // dig_t 内部的第几位
    
    // 假设 dig_t 是 uint32_t，且数组是小端序 (sk[0] 是最低位字)
    // 如果是大端序 (sk[0] 是最高位字)，计算方式需调整
    return (sk[word_idx] >> bit_idx) & 1;
}

// ==========================================
// 标量乘法实现: pk = sk * pk
// 算法：Left-to-Right Binary Method (MSB first)
// ==========================================
void ecc_mul(ecc_point*& pk, const dig_t* sk, const ecc_param* param) {
    // 1. 分配临时变量，防止修改 pk 过程中破坏基点数据
    ecc_point* R = (ecc_point*)calloc(1, sizeof(ecc_point)); // 累加器 (Result)
    ecc_point* P = (ecc_point*)calloc(1, sizeof(ecc_point)); // 基点副本 (Base Point)
    
    // 2. 保存输入的基点 P (因为计算过程中 pk 会变，或者我们需要 P 用于加法)
    if (pk == NULL) {
        // 错误处理
        free(R); free(P);
        return;
    }
    memcpy(P, pk, sizeof(ecc_point));

    // 3. 寻找 sk 中最高位的 '1' (MSB)
    // 通常从 255 位开始向下找
    int i = ECC_BITS - 1;
    for (; i >= 0; i--) {
        if (get_bit(sk, i)) {
            break;
        }
    }

    // 4. 初始化 R
    // 为了避免处理“无穷远点”的加法逻辑，我们将 R 初始化为基点 P，
    // 然后循环从 MSB 的下一位开始。
    if (i >= 0) {
        memcpy(R, P, sizeof(ecc_point)); // R = P
        i--; // 处理下一位
    } else {
        // sk 为 0 的情况，结果为无穷远点 (通常全0)
        memset(pk, 0, sizeof(ecc_point));
        free(R); free(P);
        return;
    }

    // 5. Double-and-Add 循环
    for (; i >= 0; i--) {
        // Double: R = 2 * R
        ecc_double(R, R, param);

        // Add: 如果当前位是 1，则 R = R + P
        if (get_bit(sk, i)) {
            ecc_add(R, R, P, param);
        }
    }

    // 6. 将结果回写到 pk
    memcpy(pk, R, sizeof(ecc_point));

    // 7. 清理内存
    free(R);
    free(P);
}