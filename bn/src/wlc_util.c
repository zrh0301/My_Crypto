#include "../inc/wlc_util.h"

int util_hex2bin(unsigned char *dst, int dst_size, unsigned char *src, int src_size)
{
    int i, j;
    unsigned char c, t;

    for (i = 0, j = 0; i < src_size; i++)
    {
        c = src[i];
        if (c >= '0' && c <= '9')
            t = c - '0';
        else if (c >= 'a' && c <= 'f')
            t = c - 'a' + 10;
        else if (c >= 'A' && c <= 'F')
            t = c - 'A' + 10;
        else
            continue;
        if (j & 1)
            dst[j >> 1] |= t;
        else
            dst[j >> 1] = t << 4;
        j++;
        if (j >= (dst_size << 1))
            break;
    }

    return ((j + 1) >> 1);
}

int util_memcmp(unsigned char *a, unsigned char *b, int size)
{
    int i;

    for (i = 0; i < size; i++)
    {
        if (a[i] > b[i])
            return 1;
        if (a[i] < b[i])
            return -1;
    }
    return 0;
}