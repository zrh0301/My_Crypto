#ifndef _WLC_TEST_H_
#define _WLC_TEST_H_

#include "wlc_bn.h"
#include "wlc_rand.h"
#include "common.h"

#include <string.h>
#include <stdio.h>

/*============================================================================*/
/* Macro definitions                                                          */
/*============================================================================*/

/**
 * Tests a sequence of commands to see if they respect some property.
 *
 * @param[in] _DESCRIPTION		- the property description.
 */
#define TEST_START(_DESCRIPTION, _TESTS)                                                    \
    printf("Testing if " _DESCRIPTION "...%*c", (int)(64 - strlen(_DESCRIPTION)), ' '); \
    for (int _i = 0; _i < _TESTS; _i++)                                                     \
    {

/**
 * Asserts a condition.
 *
 * If the condition is not satisfied, a unconditional jump is made to the passed
 * label.
 *
 * @param[in] _CONDITION		- the condition to assert.
 * @param[in] _LABEL			- the label to jump if the condition is no satisfied.
 */
#define TEST_ITEM(_CONDITION, _LABEL) \
    if (!(_CONDITION))                \
    {                                 \
        printf("[FAIL]!!!!!!!!\n");       \
        goto _LABEL;                  \
    }

/**
 * Finalizes a test printing the test result.
 */
#define TEST_FINAL \
    }              \
    printf("[PASS]\n");

/*============================================================================*/
/* Function definitions                                                       */
/*============================================================================*/

void test_bn(void);

#endif /* !_WLC_TEST_H_ */