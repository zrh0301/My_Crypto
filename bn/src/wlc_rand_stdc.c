
#include "../inc/wlc_rand.h"

void rand_seed(uint8_t* buf, int size)
{
	unsigned int seed;
	int i;

	for (seed = (unsigned int)time(NULL), i = 0; i < size; i++)
	{
		seed += buf[i] << ((i * 3) % 0x1F);
	}
	srand(seed);
}