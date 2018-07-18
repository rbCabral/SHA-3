#ifndef _BENCH_H_
#define _BENCH_H_
/**
 * Bench a function FUNC, 
 * that receives many 
 * __VA_ARGS__ parameters.
 */
#define BENCH_FUNCTION(FUNC,...) \
do{ \
	int _bench_i_ = 0,_bench_j_ = 0;\
	BENCH_BEGIN( #FUNC, BENCH) { \
		BENCH_ADD(FUNC(__VA_ARGS__), BENCH); \
	} \
	BENCH_END(BENCH); \
}while(0)

/**
 * Runs a new benchmark once.
 *
 * @param[in] LABEL			- the label for this benchmark.
 * @param[in] FUNCTION		- the function to benchmark.
 */
#define BENCH_ONCE(LABEL, FUNCTION)										\
	bench_reset();														\
	printf("Initializes the %s ... \n", LABEL);							\
	bench_before();														\
	FUNCTION;															\

/**
 * Runs a new benchmark.
 *
 * @param[in] LABEL			- the label for this benchmark.
 */
#define BENCH_BEGIN(LABEL, BENCH)										\
	bench_reset();														\
	printf("Runs the benchmark of the %s ... \n", LABEL);				\
	for (_bench_i_ = 0; _bench_i_ < BENCH; _bench_i_++)	{				\

/**
 * Prints the mean timing of each execution in nanoseconds.
 */
#define BENCH_END(BENCH)												\
	}																	\
	bench_compute(BENCH * BENCH);										\
	bench_print()														\

/**
 * Measures the time of one execution and adds it to the benchmark total.
 *
 * @param[in] FUNCTION		- the function executed.
 */
#define BENCH_ADD(FUNCTION, BENCH)										\
	FUNCTION;															\
	bench_before();														\
	for (_bench_j_ = 0; _bench_j_ < BENCH; _bench_j_++) {				\
		FUNCTION;														\
	}																	\
	bench_after();														\


/*====================================================================*/
/* Function prototypes                                                */
/*====================================================================*/

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
unsigned long long bench_get_total(void);

#endif
