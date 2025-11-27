/**
 * @defgroup rand Pseudo-random number generators.
 */

/**
 * @file
 *
 * Interface of the module for pseudo-random number generation.
 *
 * @ingroup rand
 */
#ifndef _WLC_RAND_H_
#define _WLC_RAND_H_

#include "wlc_types.h"

#include <stdlib.h>
#include <time.h>

 /*============================================================================*/
 /* Function definitions                                                        */
 /*============================================================================*/

#ifdef __cplusplus
extern "C"
{ /* start of __cplusplus */
#endif

/**
 * Sets the initial state of the pseudo-random number generator.
 *
 * @param[in] buf			- the buffer that represents the initial state.
 * @param[in] size			- the number of bytes.
 */
void rand_seed(uint8_t* buf, int size);

#ifdef __cplusplus
} /* end of __cplusplus */
#endif

#endif /* !_WLC_RAND_H_ */
