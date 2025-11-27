#include "../inc/wlc_test.h"

#define WORDLEN 8

void test_bn(void)
{

	//m=8542D69E4C044F18E8B92435BF6FF7DD297720630485628D5AE74EE7C32E79B7
	unsigned char byte_m[] = { 0x85, 0x42, 0xD6, 0x9E, 0x4C, 0x04, 0x4F, 0x18, 0xE8, 0xB9, 0x24, 0x35, 0xBF, 0x6F, 0xF7, 0xDD,
		0x29, 0x77, 0x20, 0x63, 0x04, 0x85, 0x62, 0x8D, 0x5A, 0xE7, 0x4E, 0xE7, 0xC3, 0x2E, 0x79, 0xB7 };

	dig_t a[WORDLEN] = { 0 };
	dig_t b[WORDLEN] = { 0 };
	dig_t c[WORDLEN << 1] = { 0 };
	dig_t d[WORDLEN << 1] = { 0 };
	dig_t e[WORDLEN] = { 0 };
	dig_t zero[WORDLEN];

	//Zero
	for (int i = 0; i < WORDLEN; i++)
		zero[i] = 0;

	dig_t m[WORDLEN];
	ByteToBN(byte_m, 32, m, WORDLEN);

	int tests = 1000;

	printf("+------------------------------------------------------------------------------+\n");
	printf("+---    Test module  --->  bn\n");

	rand_seed(0, 0);


	TEST_START("copy and comparison is correct", tests)
	{
		bn_rand(a, WORDLEN);
		bn_rand(b, WORDLEN);
		bn_rand(c, WORDLEN);

		bn_copy(c, a, WORDLEN);
		TEST_ITEM(bn_cmp(a, c, WORDLEN) == 0, END);
		bn_copy(b, a, WORDLEN);
		TEST_ITEM(bn_cmp(a, b, WORDLEN) == 0, END);
	}
	TEST_FINAL;

	TEST_START("bn_lsh_low is correct", 1)
	{
		char* charbuf_a = "E7933E845253E808990A81E00D2C0443F6F505392A29E6916A53128A5B6E5E9C";
		char* charbuf_b = "D08A4A7D011321503C01A580887EDEA0A725453CD22D4A62514B6DCBD3800000";

		unsigned char bytebuf[200];
		int charlen = 0;
		int bytelen = 0;
		int result = 0;

		charlen = 64;
		result = CharToByte(charbuf_a, charlen, bytebuf, &bytelen);
		result = ByteToBN(bytebuf, bytelen, a, WORDLEN);

		charlen = 64;
		result = CharToByte(charbuf_b, charlen, bytebuf, &bytelen);
		result = ByteToBN(bytebuf, bytelen, b, WORDLEN);

		bn_lsh_low(a, a, 21, WORDLEN);
		TEST_ITEM(bn_cmp(a, b, WORDLEN) == 0, END);
	}
	TEST_FINAL;

	TEST_START("bn_rsh_low is correct", 1)
	{
		char* charbuf_a = "E7933E845253E808990A81E00D2C0443F6F505392A29E6916A53128A5B6E5E9C";
		char* charbuf_b = "0000073C99F422929F4044C8540F006960221FB7A829C9514F348B52989452DB";

		unsigned char bytebuf[200];
		int charlen = 0;
		int bytelen = 0;
		int result = 0;

		charlen = 64;
		result = CharToByte(charbuf_a, charlen, bytebuf, &bytelen);
		result = ByteToBN(bytebuf, bytelen, a, WORDLEN);

		charlen = 64;
		result = CharToByte(charbuf_b, charlen, bytebuf, &bytelen);
		result = ByteToBN(bytebuf, bytelen, b, WORDLEN);

		bn_rsh_low(a, a, 21, WORDLEN);
		TEST_ITEM(bn_cmp(a, b, WORDLEN) == 0, END);
	}
	TEST_FINAL;

	TEST_START("mod_addition & bn_sub is correct", 1)
	{
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

		bn_mod_add(c, a, b, m, WORDLEN);
		bn_mod_add(d, b, a, m, WORDLEN);
		TEST_ITEM(bn_cmp(c, d, WORDLEN) == 0, END);
	}
	TEST_FINAL;

	TEST_START("mod_subtraction is correct", 1)
	{
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

		bn_mod_sub(c, a, b, m, WORDLEN);
		bn_mod_sub(d, b, a, m, WORDLEN);
		bn_mod_sub(e, zero, c, m, WORDLEN);
		TEST_ITEM(bn_cmp(e, d, WORDLEN) == 0, END);
	}
	TEST_FINAL;

	TEST_START("mod_hlv is correct", 1)
	{
		bn_rand(a, WORDLEN);
		if (bn_cmp(a, m, WORDLEN) >=0)
		{
			bn_sub(a, a, m, WORDLEN);
		}
		
		bn_mod_hlv(b, a, m, WORDLEN);
		bn_mod_add(c, b, b, m, WORDLEN);
		if (bn_cmp(c, m, WORDLEN) >=0)
		{
			bn_sub(c, c, m, WORDLEN);
		}
		TEST_ITEM(bn_cmp(c, a, WORDLEN) == 0, END);
	}
	TEST_FINAL;

	TEST_START("bn_multiply is correct", 1)
	{
		char* charbuf_a = "E7933E845253E808990A81E00D2C0443F6F505392A29E6916A53128A5B6E5E9C";
		char* charbuf_b = "5B3472D55740B0938241A245B0F95667DEA6A679C2A41C10E6BFE4F511780CE2";
		char* charbuf_s = "5280C9024A6C0C1150B4D1F1376002EA123AEB696D32C9BDAD2EDCF7E96E7C8C0141CF47ED8E3B5E9659D30569C6547EF56E59AD6A980A782E6985BE98FED5B8";

		unsigned char bytebuf[200];
		int charlen = 0;
		int bytelen = 0;
		int result = 0;

		charlen = 64;
		result = CharToByte(charbuf_a, charlen, bytebuf, &bytelen);
		result = ByteToBN(bytebuf, bytelen, a, WORDLEN);

		charlen = 64;
		result = CharToByte(charbuf_b, charlen, bytebuf, &bytelen);
		result = ByteToBN(bytebuf, bytelen, b, WORDLEN);

		charlen = 128;
		result = CharToByte(charbuf_s, charlen, bytebuf, &bytelen);
		result = ByteToBN(bytebuf, bytelen, d, (WORDLEN << 1));

		//c = a*b mod M
		bn_mul(c, a, b, WORDLEN);

		TEST_ITEM(bn_cmp(c, d, (WORDLEN << 1)) == 0, END);
	}
	TEST_FINAL;

	 TEST_START("mod_multiply is correct", 1)
	 {
	 	char* charbuf_a = "E7933E845253E808990A81E00D2C0443F6F505392A29E6916A53128A5B6E5E9C";
	 	char* charbuf_b = "5B3472D55740B0938241A245B0F95667DEA6A679C2A41C10E6BFE4F511780CE2";
	 	char* charbuf_s = "7150FC3BC4E380963B59A304957119F4D6D573EE4893FD44CB5C38D32D943606";

	 	unsigned char bytebuf[200];
	 	int charlen = 0;
	 	int bytelen = 0;
	 	int result = 0;

	 	charlen = 64;
	 	result = CharToByte(charbuf_a, charlen, bytebuf, &bytelen);
	 	result = ByteToBN(bytebuf, bytelen, a, WORDLEN);

	 	charlen = 64;
 		result = CharToByte(charbuf_b, charlen, bytebuf, &bytelen);
	 	result = ByteToBN(bytebuf, bytelen, b, WORDLEN);

	 	charlen = 64;
	 	result = CharToByte(charbuf_s, charlen, bytebuf, &bytelen);
	 	result = ByteToBN(bytebuf, bytelen, d, WORDLEN);

	 	//c = a*b mod M
	 	bn_mod_mul(c, a, b, m, WORDLEN);

	 	TEST_ITEM(bn_cmp(c, d, WORDLEN) == 0, END);
	 }
	 TEST_FINAL;

	 TEST_START("mod_exponent is correct", 1)
	 {
	 	char* charbuf_a = "E7933E845253E808990A81E00D2C0443F6F505392A29E6916A53128A5B6E5E9C";

	 	char* charbuf_b = "20268CCD1BA7CD9B4AC3CD107DD0417ECD964936A3B0B50A65EBA12A0BD66A85";
	 	char* charbuf_s = "0C76CD59CAA8AE970B279C2C3C4F317168C5126AE4A23F0E51A6D6285370471D";

	 	// char* charbuf_b = "0000000000000000000000000000000000000000000000000000000000000005";
	 	// char* charbuf_s = "6EB9EF75E672ACD9AB5F9A993D07A8BB03C596C14FAF887BD1454BC1269B2010";

	 	unsigned char bytebuf[200];
	 	int charlen = 0;
	 	int bytelen = 0;
	 	int result = 0;

	 	charlen = 64;
	 	result = CharToByte(charbuf_a, charlen, bytebuf, &bytelen);
	 	result = ByteToBN(bytebuf, bytelen, a, WORDLEN);

	 	charlen = 64;
	 	result = CharToByte(charbuf_b, charlen, bytebuf, &bytelen);
	 	result = ByteToBN(bytebuf, bytelen, b, WORDLEN);

	 	charlen = 64;
	 	result = CharToByte(charbuf_s, charlen, bytebuf, &bytelen);
	 	result = ByteToBN(bytebuf, bytelen, d, WORDLEN);

	 	//c = a^b mod M
	 	bn_mod_exp(c, a, b, m, WORDLEN);

	 	TEST_ITEM(bn_cmp(c, d, WORDLEN) == 0, END);
	 }
	 TEST_FINAL;

	END:

		printf("+------------------------------------------------------------------------------+\n");
	}

