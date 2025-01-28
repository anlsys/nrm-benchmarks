/*******************************************************************************
 * Copyright 2021 UChicago Argonne, LLC.
 * (c.f. AUTHORS, LICENSE)
 *
 * This file is part of the NRM Benchmarks project.
 * For more info, see https://github.com/anlsys/nrm-benchmarks
 *
 * SPDX-License-Identifier: BSD-3-Clause
 ******************************************************************************/

#include "config.h"

#include "nrm-benchmarks.h"

#include <nrm.h>
#include <math.h>

static int *key_array, *key_buff1, *key_buff2;

#pragma omp threadprivate(bucket_ptrs)
static int **bucket_size, *bucket_ptrs;

static int  TOTAL_KEYS        ;
static int  TOTAL_KS1         ;
static int  TOTAL_KS2         ;
static int  MAX_KEY           ;
static int
static int  TOTAL_KEYS        ;
static int  TOTAL_KS1         ;
static int  TOTAL_KS2         ;
static int  MAX_KEY           ;
static int
static int  NUM_BUCKETS       ;
static int  NUM_KEYS          ;
static int  SIZE_OF_BUFFERS   ;
                                    
/*****************************************************************/
/************   F  I  N  D  _  M  Y  _  S  E  E  D    ************/
/************                                         ************/
/************ returns parallel random number seq seed ************/
/*****************************************************************/

/*
 * Create a random number sequence of total length nn residing
 * on np number of processors.  Each processor will therefore have a
 * subsequence of length nn/np.  This routine returns that random
 * number which is the first random number for the subsequence belonging
 * to processor rank kn, and which is used as seed for proc kn ran # gen.
 */

double   find_my_seed( int kn,        /* my processor rank, 0<=kn<=num procs */
                       int np,        /* np = num procs                      */
                       long nn,       /* total num of ran numbers, all procs */
                       double s,      /* Ran num seed, for ex.: 314159265.00 */
                       double a )     /* Ran num gen mult, try 1220703125.00 */
{

      double t1,t2;
      long   mq,nq,kk,ik;

      if ( kn == 0 ) return s;

      mq = (nn/4 + np - 1) / np;
      nq = mq * 4 * kn;               /* number of rans to be skipped */

      t1 = s;
      t2 = a;
      kk = nq;
      while ( kk > 1 ) {
      	 ik = kk / 2;
         if( 2 * ik ==  kk ) {
            (void)randlc( &t2, &t2 );
	    kk = ik;
	 }
	 else {
            (void)randlc( &t1, &t2 );
	    kk = kk - 1;
	 }
      }
      (void)randlc( &t1, &t2 );

      return( t1 );

}



/*****************************************************************/
/*************      C  R  E  A  T  E  _  S  E  Q      ************/
/*****************************************************************/

void	create_seq( double seed, double a )
{
	double x, s;
	INT_TYPE i, k;

#pragma omp parallel private(x,s,i,k)
    {
	INT_TYPE k1, k2;
	double an = a;
	int myid = 0, num_threads = 1;
        INT_TYPE mq;

#ifdef _OPENMP
	myid = omp_get_thread_num();
	num_threads = omp_get_num_threads();
#endif

	mq = (NUM_KEYS + num_threads - 1) / num_threads;
	k1 = mq * myid;
	k2 = k1 + mq;
	if ( k2 > NUM_KEYS ) k2 = NUM_KEYS;

	KS = 0;
	s = find_my_seed( myid, num_threads,
			  (long)4*NUM_KEYS, seed, an );

        k = MAX_KEY/4;

	for (i=k1; i<k2; i++)
	{
	    x = randlc(&s, &an);
	    x += randlc(&s, &an);
    	    x += randlc(&s, &an);
	    x += randlc(&s, &an);  

            key_array[i] = k*x;
	}
    } /*omp parallel*/
}


/*****************************************************************/
/*************             R  A  N  K             ****************/
/*****************************************************************/


void is_kernel( int iteration )
{

    INT_TYPE    i, k;
    INT_TYPE    *key_buff_ptr, *key_buff_ptr2;

    int shift = MAX_KEY_LOG_2 - NUM_BUCKETS_LOG_2;
    INT_TYPE num_bucket_keys = (1L << shift);


    key_array[iteration] = iteration;
    key_array[iteration+MAX_ITERATIONS] = MAX_KEY - iteration;


/*  Determine where the partial verify test keys are, load into  */
/*  top of array bucket_size                                     */
    for( i=0; i<TEST_ARRAY_SIZE; i++ )
        partial_verify_vals[i] = key_array[test_index_array[i]];


/*  Setup pointers to key buffers  */
    key_buff_ptr2 = key_buff2;
    key_buff_ptr = key_buff1;


#pragma omp parallel private(i, k)
  {
    INT_TYPE *work_buff, m, k1, k2;
    int myid = 0, num_threads = 1;

    myid = omp_get_thread_num();
    num_threads = omp_get_num_threads();

/*  Bucket sort is known to improve cache performance on some   */
/*  cache based systems.  But the actual performance may depend */
/*  on cache size, problem size. */

    work_buff = bucket_size[myid];

/*  Initialize */
    for( i=0; i<NUM_BUCKETS; i++ )  
        work_buff[i] = 0;

/*  Determine the number of keys in each bucket */
    #pragma omp for schedule(static)
    for( i=0; i<NUM_KEYS; i++ )
        work_buff[key_array[i] >> shift]++;

/*  Accumulative bucket sizes are the bucket pointers.
    These are global sizes accumulated upon to each bucket */
    bucket_ptrs[0] = 0;
    for( k=0; k< myid; k++ )  
        bucket_ptrs[0] += bucket_size[k][0];

    for( i=1; i< NUM_BUCKETS; i++ ) { 
        bucket_ptrs[i] = bucket_ptrs[i-1];
        for( k=0; k< myid; k++ )
            bucket_ptrs[i] += bucket_size[k][i];
        for( k=myid; k< num_threads; k++ )
            bucket_ptrs[i] += bucket_size[k][i-1];
    }


/*  Sort into appropriate bucket */
    #pragma omp for schedule(static)
    for( i=0; i<NUM_KEYS; i++ )  
    {
        k = key_array[i];
        key_buff2[bucket_ptrs[k >> shift]++] = k;
    }

/*  The bucket pointers now point to the final accumulated sizes */
    if (myid < num_threads-1) {
        for( i=0; i< NUM_BUCKETS; i++ )
            for( k=myid+1; k< num_threads; k++ )
                bucket_ptrs[i] += bucket_size[k][i];
    }


/*  Now, buckets are sorted.  We only need to sort keys inside
    each bucket, which can be done in parallel.  Because the distribution
    of the number of keys in the buckets is Gaussian, the use of
    a dynamic schedule should improve load balance, thus, performance     */

    #pragma omp for schedule(dynamic)
    for( i=0; i< NUM_BUCKETS; i++ ) {

/*  Clear the work array section associated with each bucket */
        k1 = i * num_bucket_keys;
        k2 = k1 + num_bucket_keys;
        for ( k = k1; k < k2; k++ )
            key_buff_ptr[k] = 0;

/*  Ranking of all keys occurs in this section:                 */

/*  In this section, the keys themselves are used as their 
    own indexes to determine how many of each there are: their
    individual population                                       */
        m = (i > 0)? bucket_ptrs[i-1] : 0;
        for ( k = m; k < bucket_ptrs[i]; k++ )
            key_buff_ptr[key_buff_ptr2[k]]++;  /* Now they have individual key   */
                                       /* population                     */

/*  To obtain ranks of each key, successively add the individual key
    population, not forgetting to add m, the total of lesser keys,
    to the first key population                                          */
        key_buff_ptr[k1] += m;
        for ( k = k1+1; k < k2; k++ )
            key_buff_ptr[k] += key_buff_ptr[k-1];

    }

  } /*omp parallel*/
}      

int main(int argc, char **argv)
{
	/* configuration parameters:
	 * - T is the only user-input parameter, it determines the size of the
	 *   problem. All other parameters are derived from it.
	 *   From the NAS NPB documentation:
	 *   - Class: S  W  A  B  C  D  E 
	 *   - T val: 16 20 23 25 27 31 35 
	 */

	size_t T;
	size_t total_keys_log_2, max_key_log_2, num_buckets_log_2 = 10;
	long int times;

	/* needed for performance measurement */
	int64_t sumtime = 0, mintime = INT64_MAX, maxtime = 0;
	nrm_time_t start, end;
	int num_threads;

	/* retrieve the size of the problem and initialize the rest of the
	 * configuration from it
	 */
	assert(argc == 3);
	errno = 0;
	T = strtoull(argv[1], NULL, 0);
	assert(!errno);
	/* we don't support the type switching required to make class D & E work
	 */
	assert(T < 27);
	errno = 0;
	times = strtol(argv[2], NULL, 0);
	assert(!errno);

	total_keys_log_2 = T;
	max_key_log_2 = T - 4;

	/* ensure that OpenMP is giving us the right number of threads */
#pragma omp parallel
#pragma omp master
	num_threads = omp_get_num_threads();
	int err = 0;
#pragma omp parallel
#pragma omp atomic
	err++;
	assert(num_threads == err);
	err = 0;

	/* initialization: random number generator and private array
	 * dum arrays are there to avoid dead code elimination
	 */
	/*  Generate random number sequence and subsequent keys on all procs */
	create_seq( 314159265.00,                    /* Random number gen seed */
		    1220703125.00 );                 /* Random number gen mult */

	bucket_size = (int **)calloc(sizeof(int *), num_threads);
	for (i = 0; i < num_threads; i++) {
		bucket_size[i] = (int *)calloc(sizeof(int), num_buckets);
	}

	TOTAL_KEYS = (1L << total_keys_log_2);
	TOTAL_KS1 = TOTAL_KEYS;
	TOTAL_KS2 = 1;
	MAX_KEY = (1 << max_key_log_2);
	NUM_BUCKETS = (1 << num_buckets_log_2);
	NUM_KEYS = TOTAL_KEYS;
	SIZE_OF_BUFFERS = NUM_KEYS;
                                           
	key_array = (int *)calloc(sizeof(int), SIZE_OF_BUFFERS);
	key_buff1 = (int *)calloc(sizeof(int), MAX_KEY);
	key_buff2 = (int *)calloc(sizeof(int), SIZE_OF_BUFFERS);
        **key_buff1_aptr = NULL;


#pragma omp parallel for
	for( i=0; i<num_keys; i++ )
		key_buff2[i] = 0;

	/*  Do one interation for free (i.e., untimed) to guarantee initialization of  
	    all data and code pages and respective tables */
	is_kernel( 1 );  

	/* NRM Context init */
	nrmb_init(argv[0]);

	/* this version of the benchmarks reports one progress each time it goes
	 * through the entire array.
	 */
	for(long int iter = 0; iter < times; iter++)
	{
		int64_t time;
		nrm_time_gettime(&start);

		/* the actual benchmark is quite involved,
		 * so we put it in a separate function
		 */
		is_kernel(iter);
		nrm_time_gettime(&end);

		nrmb_send_progress(1.0);

		time = nrm_time_diff(&start, &end);
		sumtime += time;
		mintime = NRMB_MIN(time, mintime);
		maxtime = NRMB_MAX(time, maxtime);
	}

	nrmb_finalize();

	/* report the configuration and timings */
	fprintf(stdout, "NRM Benchmarks:      %s\n", argv[0]);
	fprintf(stdout, "Version:             %s\n", PACKAGE_VERSION);
	fprintf(stdout, "Description: one progress per iteration, NPB IS benchmark\n");
	fprintf(stdout, "Problem size:        %zu.\n", T);
	fprintf(stdout, "Kernel was executed: %ld times.\n", times);
	fprintf(stdout, "Number of threads:   %d\n", num_threads);
	fprintf(stdout, "Time (s): avg:       %11.6f min:  %11.6f max: %11.6f\n",
		1.0E-09 * sumtime/times, 1.0E-09 * mintime, 1.0E-09 * maxtime);

	fprintf(stdout, "VALIDATION disabled\n");
	return 0;
}
