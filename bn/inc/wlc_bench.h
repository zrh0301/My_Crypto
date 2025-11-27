#ifndef _WLC_BENCH_H_
#define _WLC_BENCH_H_

#include "wlc_types.h"
#include "wlc_util.h"
#include "wlc_bn.h"
#include "wlc_rand.h"
#include "common.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*============================================================================*/
/* Macro definitions                                                          */
/*============================================================================*/

/** Number of times each benchmark is run. */

/**
 * Runs a new benchmark.
 *
 * @param[in] LABEL			- the label for this benchmark.
 */
#define BENCH_START(_LABEL, _BENCHS)                                       \
   bench_reset();                                                          \
   util_print("BENCH: " _LABEL "%*c = ", (int)(32 - strlen(_LABEL)), ' '); \
   for (int _b = 0; _b < _BENCHS; _b++)                                    \
   {

/**
 * Prints the average timing of each execution in the chosen metric.
 */
#define BENCH_FINAL(_BENCHS, _ROUNDS) \
   }                                  \
   bench_compute(_BENCHS * _ROUNDS);   \
   bench_print();

/**
 * Measures the time of one execution and adds it to the benchmark total.
 *
 * @param[in] FUNCTION		- the function executed.
 */
#define BENCH_ITEM(_FUNCTION, _ROUNDS)  \
   _FUNCTION;                           \
   bench_before();                      \
   for (int _r = 0; _r < _ROUNDS; _r++) \
   {                                    \
      _FUNCTION;                        \
   }                                    \
   bench_after();

/*============================================================================*/
/* Function definitions                                                        */
/*============================================================================*/

#ifdef __cplusplus
extern "C"
{ /* start of __cplusplus */
#endif

   /**
    * Resets the benchmark data.
    *
    * @param[in] label			- the benchmark label.
    */
   void bench_reset(void);

   /**
    * Measures the time before a benchmark is executed.
    */
   void bench_before(void);

   /**
    * Measures the time after a benchmark was started and adds it to the total.
    */
   void bench_after(void);

   /**
    * Computes the mean elapsed time between the start and the end of a benchmark.
    *
    * @param benches			- the number of executed benchmarks.
    */
   void bench_compute(int benches);

   /**
    * Prints the last benchmark.
    */
   void bench_print(void);

   /**
    * Returns the result of the last benchmark.
    *
    * @return the last benchmark.
    */
   unsigned long long bench_total(void);

   void bench_bn(void);

#ifdef __cplusplus
} /* end of __cplusplus */
#endif

#endif /* !_WLC_BENCH_H_ */