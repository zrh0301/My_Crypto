#include <stdlib.h>
#include <stdint.h>
#include <string.h>

// 假设 dig_t 是 32 位无符号整数
typedef uint32_t dig_t;
typedef uint64_t ddig_t; // 双精度，用于乘法中间结果

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
        ddig_t temp = (ddig_t)a[i] - b[i] - borrow;
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
            ddig_t val = (ddig_t)T[i + j] + (ddig_t)a[i] * b[j] + carry;
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
            ddig_t val = (ddig_t)T[i + j] + (ddig_t)u * m[j] + carry;
            T[i + j] = (dig_t)val; // 这里 T[i] 会变成 0 (理论上，实际可忽略)
            carry = (dig_t)(val >> 32);
        }
        
        // 处理进位传播
        for (int j = digs; j < digs + 2 && i + j < t_size; j++) {
            ddig_t val = (ddig_t)T[i + j] + carry;
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