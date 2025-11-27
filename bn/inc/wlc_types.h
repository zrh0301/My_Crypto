/**
 * @file
 *
 * Elementary types.
 *
 * @ingroup libwhucrypt
 */

#ifndef _WLC_TYPES_H_
#define _WLC_TYPES_H_

//#include "wlc_config.h"

#include <stdint.h>
#include <memory.h>

/*============================================================================*/
/* Constant definitions                                                       */
/*============================================================================*/

/* Bits of a word (or digit) */
#define WBITS   32

/* Bytes of a word (or digit) */
#define WBYTS  (WBITS/8)

/*============================================================================*/
/* Type definitions                                                           */
/*============================================================================*/

/**
 * Represent the multiple precision integer.
 */
#if WBITS == 8
typedef uint8_t dig_t;
#elif WBITS == 16
typedef uint16_t dig_t;
#elif WBITS == 32
typedef uint32_t    dig_t; /* The digit of multiple precision integer */
typedef uint32_t    usi_t; /* Unsigned Single-precision Integer */
typedef int32_t     ssi_t; /* Signed Single-precision Integer */
typedef uint64_t    udi_t; /* Unsigned Double-precision Integer */
typedef int64_t     sdi_t; /* Signed Double-precision Integer */

#define WMASK   0xFFFFFFFF
#define WLGBITS 5

#elif WBITS == 64
typedef uint64_t dig_t;
#else
#error "Unsupported WBITS"
#endif

#define WTODIGS(bits)    (((bits) + WBITS - 1) >> WLGBITS)


 /*
  * Parameter identifiers.
  */
typedef enum _pid_en
{
	SM2_P256 = 1,
	SM9_P256,
	BLS12_P381
} pid_st;


#define	WLC_OK		0

#define	WLC_TRUE	1
#define WLC_FLASE	0

#define WLC_ERR_VRF  -1

#endif