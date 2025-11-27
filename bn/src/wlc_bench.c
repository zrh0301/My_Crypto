#include "../inc/wlc_bench.h"

#define ANSI 5
#define WORDLEN 8
#define TIMER ANSI

typedef clock_t bench_t;

typedef struct _bench_st
{
    /** Stores the time measured before the execution of the benchmark. */
    bench_t before;
    /** Stores the time measured after the execution of the benchmark. */
    bench_t after;
    /** Stores the sum of timings for the current benchmark. */
    long long total;

	double result;
} bench_ctx;

static bench_ctx g_bench;

void bench_reset(void)
{
    g_bench.total = 0;
}

void bench_before(void)
{
    g_bench.before = clock();
}

void bench_after(void)
{
    long long result;

    g_bench.after = clock();
    result = (g_bench.after - g_bench.before) * 1000000 / CLOCKS_PER_SEC;
    g_bench.total += result;
}

void bench_compute(int benches)
{
    g_bench.result = (double)g_bench.total / benches;
}

void bench_print(void)
{
    util_print("%9.2f us", g_bench.result);

    if (g_bench.result < 0)
    {
        util_print(" (overflow or bad overhead estimation)\n");
    }
    else
    {
        util_print("\n");
    }
}

unsigned long long bench_total(void)
{
    return g_bench.total;
}

void bench_bn()
{
	//n=8542D69E4C044F18E8B92435BF6FF7DD297720630485628D5AE74EE7C32E79B7
	unsigned char byte_m[] = { 0x85, 0x42, 0xD6, 0x9E, 0x4C, 0x04, 0x4F, 0x18, 0xE8, 0xB9, 0x24, 0x35, 0xBF, 0x6F, 0xF7, 0xDD,
		0x29, 0x77, 0x20, 0x63, 0x04, 0x85, 0x62, 0x8D, 0x5A, 0xE7, 0x4E, 0xE7, 0xC3, 0x2E, 0x79, 0xB7 };

	dig_t a[MAX_BN_DIGS] = { 0 };
	dig_t b[MAX_BN_DIGS] = { 0 };
	dig_t c[MAX_BN_DIGS] = { 0 };
	dig_t d[MAX_BN_DIGS] = { 0 };
	
	dig_t m[WORDLEN];
	ByteToBN(byte_m, 32, m, WORDLEN);

	int benchs = 10;
	int rounds = 100;

	printf("+------------------------------------------------------------------------------+\n");
	printf("+---    Bench module  --->  bn\n");

	rand_seed(0, 0);

	BENCH_START("bn_set_zero", benchs);
	BENCH_ITEM(bn_set_zero(a, WORDLEN), rounds);
	BENCH_FINAL(benchs, rounds);

	BENCH_START("bn_rand", benchs);
	BENCH_ITEM(bn_rand(a, WORDLEN), rounds);
	BENCH_FINAL(benchs, rounds);

	bn_rand(a, WORDLEN);
	BENCH_START("bn_copy", benchs);
	BENCH_ITEM(bn_copy(c, a, WORDLEN), rounds);
	BENCH_FINAL(benchs, rounds);

	BENCH_START("bn_get_bits", benchs);
	BENCH_ITEM(bn_get_bits(a, WORDLEN), rounds);
	BENCH_FINAL(benchs, rounds);

	bn_rand(a, WORDLEN);
	bn_rand(b, WORDLEN);
	BENCH_START("bn_add", benchs);
	BENCH_ITEM(bn_add(c, a, b, WORDLEN), rounds);
	BENCH_FINAL(benchs, rounds);

	bn_rand(a, WORDLEN);
	bn_rand(b, WORDLEN);
	BENCH_START("bn_sub", benchs);
	BENCH_ITEM(bn_sub(c, a, b, WORDLEN), rounds);
	BENCH_FINAL(benchs, rounds);

	bn_rand(a, WORDLEN);
	bn_rand(b, WORDLEN);
	BENCH_START("bn_mul", benchs);
	BENCH_ITEM(bn_mul(c, a, b, WORDLEN), rounds);
	BENCH_FINAL(benchs, rounds);

	bn_rand(a, WORDLEN);
	bn_rand(b, WORDLEN);
	BENCH_START("bn_cmp", benchs);
	BENCH_ITEM(bn_cmp(a, b, WORDLEN), rounds);
	BENCH_FINAL(benchs, rounds);

	bn_rand(a, WORDLEN);
	BENCH_START("bn_lsh_low", benchs);
	BENCH_ITEM(bn_lsh_low(a, a, 85, WORDLEN), rounds);
	BENCH_FINAL(benchs, rounds);

	bn_rand(a, WORDLEN);
	BENCH_START("bn_rsh_low", benchs);
	BENCH_ITEM(bn_rsh_low(a, a, 85, WORDLEN), rounds);
	BENCH_FINAL(benchs, rounds);

	bn_rand(a, WORDLEN);
	bn_rand(b, WORDLEN);
	if (bn_cmp(a, m, WORDLEN) >=0)
	{
		bn_sub(a, a, m, WORDLEN);
	}
	if (bn_cmp(b, m, WORDLEN) >=0)
	{
		bn_sub(b, b, m, WORDLEN);
	}
	BENCH_START("bn_mod_add", benchs);
	BENCH_ITEM(bn_mod_add(c, a, b, m, WORDLEN), rounds);
	BENCH_FINAL(benchs, rounds);

	bn_rand(a, WORDLEN);
	bn_rand(b, WORDLEN);
	if (bn_cmp(a, m, WORDLEN) >=0)
	{
		bn_sub(a, a, m, WORDLEN);
	}
	if (bn_cmp(b, m, WORDLEN) >=0)
	{
		bn_sub(b, b, m, WORDLEN);
	}
	BENCH_START("bn_mod_sub", benchs);
	BENCH_ITEM(bn_mod_sub(c, a, b, m, WORDLEN), rounds);
	BENCH_FINAL(benchs, rounds);

	bn_rand(a, WORDLEN);
	BENCH_START("bn_mod_hlv", benchs);
	BENCH_ITEM(bn_mod_hlv(c, a, m, WORDLEN), rounds);
	BENCH_FINAL(benchs, rounds);

	 bn_rand(a, WORDLEN << 1);
	 BENCH_START("bn_mod_rdc", benchs);
	 BENCH_ITEM(bn_mod_rdc(c, a, m, WORDLEN << 1, WORDLEN), 1);
	 BENCH_FINAL(benchs, 1);

	 bn_rand(a, WORDLEN);
	 bn_rand(b, WORDLEN);
	 bn_mod_rdc(a, a, m, WORDLEN, WORDLEN);
	 bn_mod_rdc(b, b, m, WORDLEN, WORDLEN);
	 BENCH_START("bn_mod_mul", benchs);
	 BENCH_ITEM(bn_mod_mul(c, a, b, m, WORDLEN), rounds);
	 BENCH_FINAL(benchs, rounds);

	 bn_rand(a, WORDLEN);
	 bn_rand(b, WORDLEN);
	 bn_mod_rdc(a, a, m, WORDLEN, WORDLEN);
	 BENCH_START("bn_mul_exp", benchs);
	 BENCH_ITEM(bn_mod_exp(c, a, b, m, WORDLEN), rounds);
	 BENCH_FINAL(benchs, rounds);

	printf("+------------------------------------------------------------------------------+\n");
}




