#ifndef _WLC_UTIL_H_
#define _WLC_UTIL_H_

#include <stdio.h>
#include <stdlib.h>

/*============================================================================*/
/* Function definitions                                                       */
/*============================================================================*/

#ifdef __cplusplus
extern "C"{	/* start of __cplusplus */
#endif

int util_hex2bin(unsigned char *dst, int dst_size, unsigned char *src, int src_size);

int util_memcmp(unsigned char *a, unsigned char *b, int size);

#define util_print  printf

#ifdef __cplusplus
}			/* end of __cplusplus */
#endif

#endif  /* !_WLC_UTIL_H_ */