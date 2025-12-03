#include "../inc/wlc_bn.h"

void bn_set_zero(dig_t *x, int digs){
	for (int i = 0; i < digs; ++i){
		x[i] = 0;
	}
}

void bn_rand(dig_t *x, int digs){
    ssize_t n = getrandom(x,digs*sizeof(dig_t),0);
    if(n<0){
        if(errno == EINTR) bn_rand(x,digs);
    }
}

void bn_copy(dig_t *r, const dig_t *x, int digs){
	for (int i = 0; i < digs; ++i) r[i] = x[i];
}

void bn_print(const dig_t* r,int digs){
    for(int i=digs-1;i>=0;--i) printf("%08x",r[i]);
    printf("\n");
}

int bn_get_bits(const dig_t *x, int digs){
	int bits = digs << 5;
	for (int i = digs - 1; i >= 0; --i){
		if (x[i] == 0) bits -= 32;
		else{
			register dig_t cur = x[i];
			for (int j = 0; j < 32; ++j){
				cur >>= 1;
				if (cur == 0){
					bits -= (31 - j);
					return bits;
				}
			}
		}
		return 0;
	}
}



dig_t bn_add(dig_t *r, const dig_t *x, const dig_t *y, int digs){
	int i = 0;
	register udi_t carry = 0;
	for (i = 0; i < digs; ++i){
		carry += (udi_t)x[i] + (udi_t)y[i];
		r[i] = (dig_t)carry;
		carry >>= 32;
	}
	return (dig_t)carry;
}

dig_t bn_sub(dig_t *r, const dig_t *x, const dig_t *y, int digs){
	register dig_t carry = 0;
	for (int i = 0; i < digs; ++i){
		sdi_t tpx = (sdi_t)x[i];
		sdi_t tpy = (sdi_t)y[i];
		r[i] = x[i] - y[i] - carry;
		carry = ((tpx - carry) < tpy);
	}
	return carry;
}

void bn_mul(dig_t* r,const dig_t* x,const dig_t* y,int digs){
    bn_set_zero(r,digs<<1);
    for(int i=0;i<digs;++i){
        for(int j=0;j<digs;++j){
            register udi_t tmp=(udi_t)x[i]*(udi_t)y[j];
            dig_t tmp_high=tmp>>32;
            dig_t tmp_low=(dig_t)tmp;
            r[i+j+1]+=bn_add(r+i+j,r+i+j,&tmp_low,1);
            r[i+j+2]+=bn_add(r+i+j+1,r+i+j+1,&tmp_high,1);
        }
    }
}




int bn_cmp(const dig_t *x, const dig_t *y, int digs){
	for (int i = digs - 1; i >= 0; --i){
		if (x[i] < y[i]) return -1;
		else if (x[i] > y[i]) return 1;
	}
	return 0;
}



dig_t bn_lsh_low(dig_t *r, const dig_t *x, int bits, int digs){
	if (bits >= 32){
		while (bits >= 31){
			bn_lsh_low(x, x, 31, digs);
			bits -= 31;
		}
		return bn_lsh_low(x, x, bits, digs);
	}

	dig_t carry, carry_new = 0;

	for (int i = 0; i < digs; ++i){
		carry = carry_new;
		carry_new = x[i] >> (32 - bits);
		r[i] = (x[i] << bits) + carry;
	}
	return carry_new;
}

dig_t bn_rsh_low(dig_t *r, const dig_t *x, int bits, int digs){
	if (bits >= 32){
		while (bits >= 31){
			bn_rsh_low(x, x, 31, digs);
			bits -= 31;
		}
		return bn_rsh_low(x, x, bits, digs);
	}

	dig_t carry, carry_new = 0;

	for (int i = digs - 1; i >= 0; --i){
		carry = carry_new;
		carry_new = (*(x + i)) << 32 - bits;
		r[i] = (x[i] >> bits) + carry;
	}
	return carry_new >> (32 - bits);
}



void bn_mod_add(dig_t *r, const dig_t *x, const dig_t *y, const dig_t *m, int digs){
	register dig_t *temp = calloc((digs + 1), sizeof(dig_t));
	register dig_t *m_new = calloc((digs + 1), sizeof(dig_t));
	bn_copy(m_new, m, digs);

	temp[digs] = bn_add(temp, x, y, digs);
	while (bn_cmp(temp, m_new, digs + 1) >= 0) bn_sub(temp, temp, m_new, digs + 1);

	bn_copy(r, temp, digs);
	free(temp);
	free(m_new);
}

void bn_mod_sub(dig_t *r, const dig_t *x, const dig_t *y, const dig_t *m, int digs) {
    // 比较 x 和 y
    int cmp = 0;
    for (int i = digs - 1; i >= 0; --i) {
        if (x[i] > y[i]) {
            cmp = 1;
            break;
        } else if (x[i] < y[i]) {
            cmp = -1;
            break;
        }
    }
    
    if (cmp >= 0) {
        // x >= y，直接减法
        bn_sub(r, x, y, digs);
    } else {
        // x < y，计算 r = x - y + m
        dig_t carry = 0;
        for (int i = 0; i < digs; ++i) {
            // 先计算 m[i] - y[i]
            dig_t diff = m[i] - y[i];
            
            // 加上 x[i] 和进位
            dig_t sum = diff + x[i] + carry;
            
            // 更新进位
            if (diff < y[i]) {
                carry = 1;  // 借位
            } else if (diff == y[i] && carry) {
                carry = 1;
            } else {
                // 检查加法是否溢出
                carry = (sum < diff) || (carry && sum == diff) ? 1 : 0;
            }
            
            r[i] = sum;
        }
        
        // 如果还有进位，需要再次减去 m
        if (carry) {
            bn_sub(r, r, m, digs);
        }
    }
}

void bn_mod_hlv(dig_t *r, const dig_t *x, const dig_t *m, int digs){
	if(x[0] & 1 == 1) bn_add(r, x, m, digs);
	bn_rsh_low(r, x, 1, digs);
	while(bn_cmp(r, m, digs) >= 0) bn_sub(r, r, m, digs);
}

void bn_mod_rdc(dig_t *r, const dig_t *x, const dig_t *m, int xdigs, int mdigs)
{
	dig_t *carry = calloc(mdigs + 1, sizeof(dig_t));
	dig_t *m_new = calloc(mdigs + 1, sizeof(dig_t));
	dig_t *carry_now = calloc(mdigs + 1, sizeof(dig_t));
	dig_t *r_new = calloc(mdigs + 1, sizeof(dig_t));
	dig_t *temp = calloc(mdigs + 1, sizeof(dig_t));
	dig_t *x_new = calloc(xdigs / mdigs * mdigs + ((xdigs % mdigs == 0) ? 0 : 1) * mdigs, sizeof(dig_t));
	bn_copy(m_new, m, mdigs);
	bn_copy(x_new, x, xdigs);

	carry[mdigs] = 1;
	carry_now[0] = 1;

	for (int i = 0; i < xdigs; i += mdigs){
		bn_mod_mul(temp, x_new + i, carry_now, m, mdigs);
		bn_add(r_new, temp, r_new, mdigs + 1);

		while (bn_cmp(r_new, m_new, mdigs + 1) >= 0) bn_sub(r_new, r_new, m_new, mdigs + 1);
		bn_mod_mul(carry_now, carry_now, carry, m_new, mdigs + 1);
	}

	bn_copy(r, r_new, mdigs);
	free(carry);
	free(m_new);
	free(carry_now);
	free(r_new);
	free(temp);
	free(x_new);
}

void bn_mod_mul(dig_t *r, const dig_t *x, const dig_t *y, const dig_t *m, int digs){
	int digs_new = digs + 1;
	dig_t *zero_bn = calloc(digs, sizeof(dig_t));
	dig_t *x_new = calloc(digs_new, sizeof(dig_t));
	dig_t *r_new = calloc(digs_new, sizeof(dig_t));
	dig_t *y_new = calloc(digs, sizeof(dig_t));
	dig_t *m_new = calloc(digs_new, sizeof(dig_t));
	bn_copy(x_new, x, digs);
	bn_copy(y_new, y, digs);
	bn_copy(m_new, m, digs);

	while (bn_cmp(x_new, m_new, digs) >= 0) bn_sub(x_new, x_new, m_new, digs);
	while (bn_cmp(y_new, m_new, digs) >= 0) bn_sub(y_new, y_new, m_new, digs);

	while (bn_cmp(y_new, zero_bn, digs)){
		if (y_new[0] & 1 == 1){
			bn_add(r_new, r_new, x_new, digs + 1);
			while (bn_cmp(r_new, m_new, digs + 1) >= 0) bn_sub(r_new, r_new, m_new, digs + 1);
		}

		bn_lsh_low(x_new, x_new, 1, digs_new);
		while (bn_cmp(x_new, m_new, digs_new) > 0) bn_sub(x_new, x_new, m_new, digs_new);
		bn_rsh_low(y_new, y_new, 1, digs);
	}

	bn_copy(r, r_new, digs);
	free(x_new);
	free(y_new);
	free(m_new);
	free(r_new);
	free(zero_bn);
}

/*
// Barrett 预计算结构
typedef struct {
    dig_t *mu;      // μ = floor(2^(2k) / m)
    dig_t *m;       // 模数
    int k;          // 模数的位数(以dig_t为单位)
    int mu_len;     // μ的长度
} barrett_ctx_t;

// 初始化 Barrett 上下文（只需计算一次）
void barrett_init(barrett_ctx_t *ctx, const dig_t *m, int digs) {
    ctx->k = digs;
    ctx->m = calloc(digs, sizeof(dig_t));
    ctx->mu = calloc(digs + 2, sizeof(dig_t));  // μ 最多 digs+1 位
    ctx->mu_len = digs + 1;
    
    bn_copy(ctx->m, m, digs);
    
    // 计算 μ = floor(2^(2k*DIGIT_BITS) / m)
    // 即计算 2^(2k*DIGIT_BITS) 然后除以 m
    int double_len = 2 * digs + 1;
    dig_t *dividend = calloc(double_len, sizeof(dig_t));
    dividend[2 * digs] = 1;  // 2^(2k*DIGIT_BITS)
    
    // 大数除法: μ = dividend / m
    bn_div(ctx->mu, NULL, dividend, double_len, m, digs);
    
    free(dividend);
}

// 释放 Barrett 上下文
void barrett_free(barrett_ctx_t *ctx) {
    free(ctx->mu);
    free(ctx->m);
}

// Barrett 模约减: r = x mod m
// 输入 x 的长度最多为 2*k
void barrett_reduce(dig_t *r, const dig_t *x, int x_len, const barrett_ctx_t *ctx) {
    int k = ctx->k;
    
    // 如果 x 已经小于 m，直接返回
    if (x_len <= k && bn_cmp(x, ctx->m, k) < 0) {
        bn_copy(r, x, k);
        return;
    }
    
    // 步骤1: q1 = floor(x / 2^((k-1)*DIGIT_BITS))
    // 即取 x 的高位部分（右移 k-1 个 digit）
    int q1_len = x_len - (k - 1);
    if (q1_len <= 0) {
        bn_copy(r, x, k);
        while (bn_cmp(r, ctx->m, k) >= 0) {
            bn_sub(r, r, ctx->m, k);
        }
        return;
    }
    const dig_t *q1 = x + (k - 1);  // 指向 x 的高位
    
    // 步骤2: q2 = q1 * μ
    int q2_len = q1_len + ctx->mu_len;
    dig_t *q2 = calloc(q2_len, sizeof(dig_t));
    bn_mul(q2, q1, ctx->mu, ctx->mu_len);
    
    // 步骤3: q3 = floor(q2 / 2^((k+1)*DIGIT_BITS))
    // 即取 q2 的高位（右移 k+1 个 digit）
    int shift = k + 1;
    dig_t *q3;
    int q3_len;
    if (q2_len > shift) {
        q3 = q2 + shift;
        q3_len = q2_len - shift;
    } else {
        // q3 = 0，直接用减法处理
        bn_copy(r, x, k + 1 < x_len ? k + 1 : x_len);
        while (bn_cmp(r, ctx->m, k) >= 0) {
            bn_sub(r, r, ctx->m, k);
        }
        free(q2);
        return;
    }
    
    // 步骤4: r1 = x mod 2^((k+1)*DIGIT_BITS)
    // 即取 x 的低 k+1 个 digit
    int r1_len = (x_len < k + 1) ? x_len : k + 1;
    dig_t *r1 = calloc(k + 2, sizeof(dig_t));
    bn_copy(r1, x, r1_len);
    
    // 步骤5: r2 = (q3 * m) mod 2^((k+1)*DIGIT_BITS)
    dig_t *q3m = calloc(q3_len + k + 1, sizeof(dig_t));
    bn_mul(q3m, q3, ctx->m, k);
    
    dig_t *r2 = calloc(k + 2, sizeof(dig_t));
    int r2_copy_len = (q3_len + k < k + 1) ? q3_len + k : k + 1;
    bn_copy(r2, q3m, r2_copy_len);
    
    // 步骤6: r = r1 - r2
    dig_t *result = calloc(k + 2, sizeof(dig_t));
    if (bn_cmp(r1, r2, k + 1) >= 0) {
        bn_sub(result, r1, r2, k + 1);
    } else {
        // r1 < r2，需要借位：r1 + 2^((k+1)*DIGIT_BITS) - r2
        dig_t *borrow = calloc(k + 2, sizeof(dig_t));
        borrow[k + 1] = 1;  // 2^((k+1)*DIGIT_BITS)
        bn_add(r1, r1, borrow, k + 2);
        bn_sub(result, r1, r2, k + 2);
        free(borrow);
    }
    
    // 步骤7: 最多需要减去 2 次 m
    while (bn_cmp(result, ctx->m, k + 1) >= 0) {
        bn_sub(result, result, ctx->m, k + 1);
    }
    
    bn_copy(r, result, k);
    
    // 清理
    free(q2);
    free(r1);
    free(r2);
    free(q3m);
    free(result);
}

// 优化后的模乘法
void bn_mod_mul_barrett(dig_t *r, const dig_t *x, const dig_t *y, 
                        const barrett_ctx_t *ctx) {
    int k = ctx->k;
    int prod_len = 2 * k;
    
    // 步骤1: 计算完整乘积 z = x * y
    dig_t *z = calloc(prod_len, sizeof(dig_t));
    bn_mul(z, x, y, k);
    
    // 步骤2: Barrett 约减
    barrett_reduce(r, z, prod_len, ctx);
    
    free(z);
}

// 简化接口（兼容原函数签名）
void bn_mod_mul(dig_t *r, const dig_t *x, const dig_t *y, const dig_t *m, int digs) {
    barrett_ctx_t ctx;
    barrett_init(&ctx, m, digs);
    bn_mod_mul_barrett(r, x, y, &ctx);
    barrett_free(&ctx);
}
*/
/*
void bn_mod_exp(dig_t *r, const dig_t *x, const dig_t *e, const dig_t *m, int digs){
	dig_t *x_new = calloc(digs, sizeof(dig_t));
	dig_t *e_new = calloc(digs, sizeof(dig_t));
	dig_t *zero_bn = calloc(digs, sizeof(dig_t));
	bn_copy(x_new, x, digs);
	bn_copy(e_new, e, digs);

	bn_set_zero(r, digs);
	r[0] = 1;

	while (bn_cmp(x_new, m, digs) >= 0) bn_sub(x_new, x_new, m, digs);

	register dig_t carry = 0;
	for (int count = 0; count < digs << 5; ++count){
		carry = bn_lsh_low(e_new, e_new, 1, digs);
		bn_mod_mul(r, r, r, m, digs);
		if (carry) bn_mod_mul(r, r, x_new, m, digs);
	}

	free(x_new);
	free(e_new);
	free(zero_bn);
}
*/

// ==========================================
// 1. 辅助工具函数
// ==========================================

// 计算大数的实际比特长度 (用于跳过前导零)
int bn_bit_len(const dig_t *a, int digs) {
    for (int i = digs - 1; i >= 0; i--) {
        if (a[i] != 0) {
            int bits = 0;
            dig_t val = a[i];
            while (val) { bits++; val >>= 1; }
            return i * 32 + bits;
        }
    }
    return 0;
}

// 检查大数第 bit_index 位是否为 1
int bn_check_bit(const dig_t *a, int bit_index) {
    int word_idx = bit_index / 32;
    int bit_pos = bit_index % 32;
    return (a[word_idx] >> bit_pos) & 1;
}

// 计算蒙哥马利常数 n0 = -m^(-1) mod 2^32
// 用于约减步骤
dig_t bn_calc_monty_n0(dig_t m0) {
    dig_t y = 1;
    // 牛顿迭代法求逆 (针对 2^32)
    y = y * (2 - m0 * y);
    y = y * (2 - m0 * y);
    y = y * (2 - m0 * y);
    y = y * (2 - m0 * y);
    y = y * (2 - m0 * y);
    return (dig_t)(0 - y);
}

// 简单的比较函数: 1 if a>b, -1 if a<b, 0 if a==b
int bn_cmp_safe(const dig_t *a, const dig_t *b, int digs) {
    for (int i = digs - 1; i >= 0; i--) {
        if (a[i] > b[i]) return 1;
        if (a[i] < b[i]) return -1;
    }
    return 0;
}

// 简单减法: r = a - b (假设 a >= b)
void bn_sub_safe(dig_t *r, const dig_t *a, const dig_t *b, int digs) {
    dig_t borrow = 0;
    for (int i = 0; i < digs; i++) {
        udi_t temp = (udi_t)a[i] - b[i] - borrow;
        r[i] = (dig_t)temp;
        borrow = (temp >> 32) & 1;
    }
}

// ==========================================
// 2. 核心: 蒙哥马利乘法 (A * B * R^-1 mod m)
// ==========================================
void bn_mon_mul(dig_t *r, const dig_t *a, const dig_t *b, const dig_t *m, dig_t n0, int digs) {
    // T 需要 2 * digs + 1 的空间来存储中间乘积结果
    int t_size = 2 * digs + 1;
    dig_t *T = (dig_t *)calloc(t_size, sizeof(dig_t));

    // --- 步骤 1: 计算 T = a * b (普通大数乘法) ---
    for (int i = 0; i < digs; i++) {
        dig_t carry = 0;
        for (int j = 0; j < digs; j++) {
            udi_t val = (udi_t)T[i + j] + (udi_t)a[i] * b[j] + carry;
            T[i + j] = (dig_t)val;
            carry = (dig_t)(val >> 32);
        }
        T[i + digs] = carry;
    }

    // --- 步骤 2: 蒙哥马利约减 (Reduction) ---
    // 这一步让 T 能够被 R (2^(32*digs)) 整除，同时 mod m 不变
    for (int i = 0; i < digs; i++) {
        // u = T[i] * n0 mod 2^32
        dig_t u = T[i] * n0;
        
        // T += m * u * (2^(32*i))
        // 我们只需要从 T[i] 开始加，因为低位已经处理完了
        dig_t carry = 0;
        for (int j = 0; j < digs; j++) {
            udi_t val = (udi_t)T[i + j] + (udi_t)u * m[j] + carry;
            T[i + j] = (dig_t)val; // 这里 T[i] 会变成 0 (理论上，实际可忽略)
            carry = (dig_t)(val >> 32);
        }
        
        // 处理进位传播
        for (int j = digs; j < digs + 2 && i + j < t_size; j++) {
            udi_t val = (udi_t)T[i + j] + carry;
            T[i + j] = (dig_t)val;
            carry = (dig_t)(val >> 32);
            if (carry == 0) break;
        }
    }

    // --- 步骤 3: 结果是 T / R ---
    // 相当于直接取 T 的高 digs 位
    // 也就是 T[digs] 到 T[2*digs-1]
    for (int i = 0; i < digs; i++) {
        r[i] = T[i + digs];
    }
    
    // --- 步骤 4: 最后的修正 ---
    // 如果结果 >= m，需要减去 m
    // 检查第 digs 位是否有溢出(T[2*digs]) 或者 r >= m
    if (T[2 * digs] || bn_cmp_safe(r, m, digs) >= 0) {
        bn_sub_safe(r, r, m, digs);
    }

    free(T);
}

// ==========================================
// 3. 预计算 R^2 mod m
// ==========================================
// 这里的 R = 2^(32 * digs)
// 我们通过移位和取模来计算 R^2 mod m，避免大数除法
void bn_calc_R2_mod_m(dig_t *R2, const dig_t *m, int digs) {
    // 初始化 R2 为 1 (或者直接从 R 开始，这里从 1 开始最通用)
    memset(R2, 0, digs * sizeof(dig_t));
    R2[0] = 1;

    // R^2 相当于 1 左移 (2 * digs * 32) 位
    // 我们每次左移 1 位，如果溢出(>=m)就减去 m
    int total_bits = 2 * digs * 32;
    
    for (int i = 0; i < total_bits; i++) {
        // R2 = R2 * 2
        dig_t carry = 0;
        for (int j = 0; j < digs; j++) {
            dig_t next_carry = (R2[j] >> 31) & 1;
            R2[j] = (R2[j] << 1) | carry;
            carry = next_carry;
        }
        
        // 如果 R2 >= m，则 R2 = R2 - m
        if (carry || bn_cmp_safe(R2, m, digs) >= 0) {
            bn_sub_safe(R2, R2, m, digs);
        }
    }
}

// ==========================================
// 4. 主函数: 蒙哥马利模幂
// ==========================================
void bn_mod_exp(dig_t *r, const dig_t *x, const dig_t *e, const dig_t *m, int digs){
    if (!x || !e || !m || !r || digs <= 0) return;

    // 1. 准备工作
    dig_t *R2 = (dig_t*)calloc(digs, sizeof(dig_t));
    dig_t *x_mon = (dig_t*)calloc(digs, sizeof(dig_t));
    dig_t *t_res = (dig_t*)calloc(digs, sizeof(dig_t)); // 临时结果
    dig_t *one = (dig_t*)calloc(digs, sizeof(dig_t));
    one[0] = 1;

    // 计算魔法常数 n0
    dig_t n0 = bn_calc_monty_n0(m[0]);

    // 2. 预计算 R^2 mod m
    bn_calc_R2_mod_m(R2, m, digs);

    // 3. 将输入转换到蒙哥马利域
    // x_mon = x * R mod m = MonMul(x, R^2)
    bn_mon_mul(x_mon, x, R2, m, n0, digs);
    
    // r 的初始值 (对应 1)
    // r_mon = 1 * R mod m = MonMul(1, R^2)
    bn_mon_mul(t_res, one, R2, m, n0, digs);

    // 4. 核心循环: 从左向右二进制幂 (Square and Multiply)
    int bits = bn_bit_len(e, digs);
    
    for (int i = bits - 1; i >= 0; i--) {
        // Square: res = res * res
        bn_mon_mul(t_res, t_res, t_res, m, n0, digs);

        if (bn_check_bit(e, i)) {
            // Multiply: res = res * x
            bn_mon_mul(t_res, t_res, x_mon, m, n0, digs);
        }
    }

    // 5. 将结果转回普通域
    // r = t_res * R^-1 mod m = MonMul(t_res, 1)
    bn_mon_mul(r, t_res, one, m, n0, digs);

    // 6. 清理内存
    free(R2);
    free(x_mon);
    free(t_res);
    free(one);
}

// 检查大数是否为 0
// 返回值: 1 是 0, 0 不是 0
int bn_is_zero(const dig_t *a, int digs) {
    // 遍历每一个 limb (字)
    for (int i = 0; i < digs; i++) {
        // 只要发现任意一位不是 0，那就肯定不是 0
        if (a[i] != 0) {
            return 0; 
        }
    }
    // 全部检查通过，确实是 0
    return 1;
}

// 检查大数是否为 1
// 返回值: 1 是 1, 0 不是 1
int bn_is_one(const dig_t *a, int digs) {
    if (digs <= 0) return 0; // 防御性编程

    // 条件1: 最低位必须是 1
    if (a[0] != 1) {
        return 0;
    }

    // 条件2: 除了最低位，其他高位必须全是 0
    for (int i = 1; i < digs; i++) {
        if (a[i] != 0) {
            return 0;
        }
    }

    return 1;
}

// 辅助函数：带符号的大数减法或者模减法
// x = (x - y) mod m
// 注意：在这个算法里，x 和 y 的差可能是负数，如果是负数需要 +m
void bn_mod_sub_optimize(dig_t *x, const dig_t *y, const dig_t *m, int digs) {
    // 如果 x >= y，直接减
    if (bn_cmp(x, y, digs) >= 0) {
        bn_sub(x, x, y, digs);
    } else {
        // 如果 x < y，结果是负数，需要 x - y + m
        // 等价于 m - (y - x)
        dig_t *tmp = calloc(digs, sizeof(dig_t));
        bn_sub(tmp, y, x, digs); // tmp = y - x
        bn_sub(x, m, tmp, digs); // x = m - tmp
        free(tmp);
    }
}

// 辅助函数：模 m 除以 2
// r = a / 2 mod m (前提：m 是奇数)
void bn_mod_div2(dig_t *r, dig_t *a, const dig_t *m, int digs) {
    // 检查 a 是否为奇数 (看最低位)
    if (a[0] & 1) {
        // 是奇数：r = (a + m) >> 1
        // 注意：a+m 可能会进位溢出最高位，需要小心处理
        dig_t carry = bn_add(r, a, m, digs); 
        bn_rsh_low(r, r, 1, digs);
        if (carry) {
            // 如果有进位，进位的那 1 也是要右移进来的，变成了最高位的 0x8000...
            r[digs-1] |= ((dig_t)1 << 31); 
        }
    } else {
        // 是偶数：r = a >> 1
        bn_rsh_low(r, a, 1, digs);
    }
}

// ==========================================
// 二进制求模逆 (Binary Inversion)
// 计算 r = a^-1 mod m
// 限制：m 必须是奇数
// ==========================================
void bn_mod_inv(dig_t *r, const dig_t *a, const dig_t *m, int digs) {
    dig_t *u = calloc(digs, sizeof(dig_t));
    dig_t *v = calloc(digs, sizeof(dig_t));
    dig_t *x1 = calloc(digs, sizeof(dig_t));
    dig_t *x2 = calloc(digs, sizeof(dig_t));

    // 初始化
    bn_copy(u, a, digs);
    bn_copy(v, m, digs);
    x1[0] = 1; // x1 = 1
    bn_set_zero(x2, digs); // x2 = 0

    // 循环直到 u == 0 或 v == 0 (实际上 u=1 或 v=1 时就差不多了)
    // 这里为了通用，循环直到其中一个为 0，非 0 的那个就是 GCD
    while (!bn_is_zero(u, digs) && !bn_is_zero(v, digs)) {
        
        // 1. 消除 u 的偶因子
        while ((u[0] & 1) == 0) { // while u is even
            bn_rsh_low(u, u, 1, digs); // u = u / 2
            bn_mod_div2(x1, x1, m, digs); // x1 = x1 / 2 mod m
        }

        // 2. 消除 v 的偶因子
        while ((v[0] & 1) == 0) { // while v is even
            bn_rsh_low(v, v, 1, digs); // v = v / 2
            bn_mod_div2(x2, x2, m, digs); // x2 = x2 / 2 mod m
        }

        // 3. 减法步骤
        if (bn_cmp(u, v, digs) >= 0) {
            bn_sub(u, u, v, digs);        // u = u - v
            bn_mod_sub_optimize(x1, x2, m, digs); // x1 = x1 - x2 mod m
        } else {
            bn_sub(v, v, u, digs);        // v = v - u
            bn_mod_sub_optimize(x2, x1, m, digs); // x2 = x2 - x1 mod m
        }
    }

    // 结果处理
    // 如果 GCD(a, m) == 1，那么 u 和 v 中最后剩下的那个非零值一定是 1
    if (bn_is_one(u, digs)) {
        bn_copy(r, x1, digs);
    } else if (bn_is_one(v, digs)) {
        bn_copy(r, x2, digs);
    } else {
        // 错误：a 和 m 不互质，不存在逆元
        bn_set_zero(r, digs); 
    }

    free(u); free(v); free(x1); free(x2);
}

