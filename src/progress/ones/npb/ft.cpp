/*
MIT License

Copyright (c) 2021 Parallel Applications Modelling Group - GMAP 
	GMAP website: https://gmap.pucrs.br
	
	Pontifical Catholic University of Rio Grande do Sul (PUCRS)
	Av. Ipiranga, 6681, Porto Alegre - Brazil, 90619-900

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

------------------------------------------------------------------------------

The original NPB 3.4.1 version was written in Fortran and belongs to: 
	http://www.nas.nasa.gov/Software/NPB/

Authors of the Fortran code:
	D. Bailey
	W. Saphir
	H. Jin

------------------------------------------------------------------------------

The serial C++ version is a translation of the original NPB 3.4.1
Serial C++ version: https://github.com/GMAP/NPB-CPP/tree/master/NPB-SER

Authors of the C++ code: 
	Dalvan Griebler <dalvangriebler@gmail.com>
	Gabriell Araujo <hexenoften@gmail.com>
 	Júnior Löff <loffjh@gmail.com>

------------------------------------------------------------------------------

The OpenMP version is a parallel implementation of the serial C++ version
OpenMP version: https://github.com/GMAP/NPB-CPP/tree/master/NPB-OMP

Authors of the OpenMP code:
	Júnior Löff <loffjh@gmail.com>

*/

#include "config.h"
#include "omp.h"
#include "npb-CPP.hpp"

extern "C" {
#include "nrm-benchmarks.h"
}

#define NX               512
#define NY               256
#define NZ               256
#define MAXDIM           512
#define NITER_DEFAULT    20
#define NXP              513
#define NYP              256
#define NTOTAL           33554432
#define NTOTALP          33619968
#define DEFAULT_BEHAVIOR 1
/*
 * ---------------------------------------------------------------------
 * u0, u1, u2 are the main arrays in the problem. 
 * depending on the decomposition, these arrays will have different 
 * dimensions. to accomodate all possibilities, we allocate them as 
 * one-dimensional arrays and pass them to subroutines for different 
 * views
 * - u0 contains the initial (transformed) initial condition
 * - u1 and u2 are working arrays
 * - twiddle contains exponents for the time evolution operator. 
 * ---------------------------------------------------------------------
 * large arrays are in common so that they are allocated on the
 * heap rather than the stack. this common block is not
 * referenced directly anywhere else. padding is to avoid accidental 
 * cache problems, since all array sizes are powers of two.
 * ---------------------------------------------------------------------
 * we need a bunch of logic to keep track of how
 * arrays are laid out. 
 *
 * note: this serial version is the derived from the parallel 0D case
 * of the ft NPB.
 * the computation proceeds logically as
 *
 * set up initial conditions
 * fftx(1)
 * transpose (1->2)
 * ffty(2)
 * transpose (2->3)
 * fftz(3)
 * time evolution
 * fftz(3)
 * transpose (3->2)
 * ffty(2)
 * transpose (2->1)
 * fftx(1)
 * compute residual(1)
 * 
 * for the 0D, 1D, 2D strategies, the layouts look like xxx
 *
 *            0D        1D        2D
 * 1:        xyz       xyz       xyz
 * 2:        xyz       xyz       yxz
 * 3:        xyz       zyx       zxy
 * the array dimensions are stored in dims(coord, phase)
 * ---------------------------------------------------------------------
 * if processor array is 1x1 -> 0D grid decomposition
 * 
 * cache blocking params. these values are good for most
 * RISC processors.  
 * FFT parameters:
 * fftblock controls how many ffts are done at a time. 
 * the default is appropriate for most cache-based machines
 * on vector machines, the FFT can be vectorized with vector
 * length equal to the block size, so the block size should
 * be as large as possible. this is the size of the smallest
 * dimension of the problem: 128 for class A, 256 for class B
 * and 512 for class C.
 * ---------------------------------------------------------------------
 */
#define	FFTBLOCK_DEFAULT DEFAULT_BEHAVIOR
#define	FFTBLOCKPAD_DEFAULT DEFAULT_BEHAVIOR
#define FFTBLOCK FFTBLOCK_DEFAULT
#define FFTBLOCKPAD FFTBLOCKPAD_DEFAULT
#define	SEED 314159265.0
#define	A 1220703125.0
#define	PI 3.141592653589793238
#define	ALPHA 1.0e-6
#define	T_TOTAL 1
#define	T_SETUP 2
#define	T_FFT 3
#define	T_EVOLVE 4
#define	T_CHECKSUM 5
#define	T_FFTX 6
#define	T_FFTY 7
#define	T_FFTZ 8
#define T_MAX 8

/* global variables */
#if defined(DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION)
static dcomplex sums[NITER_DEFAULT+1];
static double twiddle[NTOTAL];
static dcomplex u[MAXDIM];
static dcomplex u0[NTOTAL];
static dcomplex u1[NTOTAL];
static int dims[3];
#else
static dcomplex (*sums)=(dcomplex*)malloc(sizeof(dcomplex)*(NITER_DEFAULT+1));
static double (*twiddle)=(double*)malloc(sizeof(double)*(NTOTAL));
static dcomplex (*u)=(dcomplex*)malloc(sizeof(dcomplex)*(MAXDIM));
static dcomplex (*u0)=(dcomplex*)malloc(sizeof(dcomplex)*(NTOTAL));
static dcomplex (*u1)=(dcomplex*)malloc(sizeof(dcomplex)*(NTOTAL));
static int (*dims)=(int*)malloc(sizeof(int)*(3));
#endif
static int niter;

/* function prototypes */
static void cffts1(int is,
		int d1,
		int d2,
		int d3,
		void* pointer_x,
		void* pointer_xout,
		dcomplex y1[][FFTBLOCKPAD],
		dcomplex y2[][FFTBLOCKPAD]);
static void cffts2(int is,
		int d1,
		int d2,
		int d3,
		void* pointer_x,
		void* pointer_xout,
		dcomplex y1[][FFTBLOCKPAD],
		dcomplex y2[][FFTBLOCKPAD]);
static void cffts3(int is,
		int d1,
		int d2,
		int d3,
		void* pointer_x,
		void* pointer_xout,
		dcomplex y1[][FFTBLOCKPAD],
		dcomplex y2[][FFTBLOCKPAD]);
static void cfftz(int is,
		int m,
		int n,
		dcomplex x[][FFTBLOCKPAD],
		dcomplex y[][FFTBLOCKPAD]);
static void checksum(int i,
		void* pointer_u1,
		int d1,
		int d2,
		int d3);
static void compute_indexmap(void* pointer_twiddle,
		int d1,
		int d2,
		int d3);
static void compute_initial_conditions(void* pointer_u0,
		int d1,
		int d2,
		int d3);
static void evolve(void* pointer_u0,
		void* pointer_u1,
		void* pointer_twiddle,
		int d1,
		int d2,
		int d3);
static void fft(int dir,
		void* pointer_x1,
		void* pointer_x2);
static void fft_init(int n);
static void fftz2(int is,
		int l,
		int m,
		int n,
		int ny,
		int ny1,
		dcomplex u[],
		dcomplex x[][FFTBLOCKPAD],
		dcomplex y[][FFTBLOCKPAD]);
static int ilog2(int n);
static void init_ui(void* pointer_u0,
		void* pointer_u1,
		void* pointer_twiddle,
		int d1,
		int d2,
		int d3);
static void ipow46(double a,
		int exponent,
		double* result);

/* ft */
int main(int argc, char **argv){
#if defined(DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION)
	printf(" DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION mode on\n");
#endif
	int i;	
	int iter;
	double total_time, mflops;
	boolean verified;

	int64_t sumtime = 0, mintime = INT64_MAX, maxtime = 0;
	nrm_time_t start, end;
	int num_threads;

	#pragma omp parallel
	#pragma omp master
	num_threads = omp_get_num_threads();

	assert(argc == 2);
	errno = 0;
	niter = atoi(argv[1]);
	assert(!errno);

	/*
	 * ---------------------------------------------------------------------
	 * run the entire problem once to make sure all data is touched. 
	 * this reduces variable startup costs, which is important for such a 
	 * short benchmark. the other NPB 2 implementations are similar. 
	 * ---------------------------------------------------------------------
	 */
	dims[0] = NX;
	dims[1] = NY;
	dims[2] = NZ;

	init_ui(u0, u1, twiddle, dims[0], dims[1], dims[2]);
	compute_indexmap(twiddle, dims[0], dims[1], dims[2]);
	compute_initial_conditions(u1, dims[0], dims[1], dims[2]);
	fft_init(MAXDIM);
	#pragma omp parallel
		fft(1, u1, u0);
	/*
	 * ---------------------------------------------------------------------
	 * start over from the beginning. note that all operations must
	 * be timed, in contrast to other benchmarks. 
	 * ---------------------------------------------------------------------
	 */
	nrmb_init(argv[0]);

	compute_indexmap(twiddle, dims[0], dims[1], dims[2]);

	compute_initial_conditions(u1, dims[0], dims[1], dims[2]);

	fft_init(MAXDIM);

	#pragma omp parallel private(iter) firstprivate(niter)
    {
		fft(1, u1, u0);
		for(iter=1; iter<=niter; iter++){
			int64_t time;
			nrm_time_gettime(&start);

			evolve(u0, u1, twiddle, dims[0], dims[1], dims[2]);
			fft(-1, u1, u1);
			#pragma omp barrier

			nrm_time_gettime(&end);
			
			#pragma omp master
			{
				nrmb_send_progress(1.0);

				time = nrm_time_diff(&start, &end);
				sumtime += time;
				mintime = NRMB_MIN(time, mintime);
				maxtime = NRMB_MAX(time, maxtime);
			}
		}
	} /* end parallel */

	nrmb_finalize();

	/* report the configuration and timings */
	fprintf(stdout, "NRM Benchmarks:      %s\n", argv[0]);
	fprintf(stdout, "Version:             %s\n", PACKAGE_VERSION);
	fprintf(stdout,
		"Description: one progress per iteration, NPB FT benchmark\n");
	fprintf(stdout, "Problem size:        %zu.\n", (size_t)NX);
	fprintf(stdout, "Kernel was executed: %d times.\n", niter);
	fprintf(stdout, "Number of threads:   %d\n", num_threads);
	fprintf(stdout, "Time (s): avg:       %11.6f min:  %11.6f max: %11.6f\n",
		1.0E-09 * sumtime / niter, 1.0E-09 * mintime, 1.0E-09 * maxtime);

	fprintf(stdout, "VALIDATION disabled\n");

	return 0;
}

static void cffts1(int is,
		int d1,
		int d2,
		int d3,
		void* pointer_x,
		void* pointer_xout,
		dcomplex y1[][FFTBLOCKPAD],
		dcomplex y2[][FFTBLOCKPAD]){
	dcomplex (*x)[NY][NX] = (dcomplex(*)[NY][NX])pointer_x;
	dcomplex (*xout)[NY][NX] = (dcomplex(*)[NY][NX])pointer_xout;

	int logd1;
	int i, j, k, jj;

	logd1 = ilog2(d1);

	
	#pragma omp for	
	for(k=0; k<d3; k++){
		for(jj=0; jj<=d2-FFTBLOCK; jj+=FFTBLOCK){
			for(j=0; j<FFTBLOCK; j++){
				for(i=0; i<d1; i++){
					y1[i][j] = x[k][j+jj][i];
				}
			}
			cfftz(is, logd1, d1, y1, y2);
			for(j=0; j<FFTBLOCK; j++){
				for(i=0; i<d1; i++){
					xout[k][j+jj][i] = y1[i][j];
				}
			}
		}
	}
}

static void cffts2(int is,
		int d1,
		int d2,
		int d3,
		void* pointer_x,
		void* pointer_xout,
		dcomplex y1[][FFTBLOCKPAD],
		dcomplex y2[][FFTBLOCKPAD]){
	dcomplex (*x)[NY][NX] = (dcomplex(*)[NY][NX])pointer_x;
	dcomplex (*xout)[NY][NX] = (dcomplex(*)[NY][NX])pointer_xout;

	int logd2;
	int i, j, k, ii;

	logd2 = ilog2(d2);

	#pragma omp for	
	for(k=0; k<d3; k++){
		for(ii=0; ii<=d1-FFTBLOCK; ii+=FFTBLOCK){
			for(j=0; j<d2; j++){
				for(i=0; i<FFTBLOCK; i++){
					y1[j][i] = x[k][j][i+ii];
				}
			}
			cfftz(is, logd2, d2, y1, y2);
			for(j=0; j<d2; j++){
				for(i=0; i<FFTBLOCK; i++){
					xout[k][j][i+ii] = y1[j][i];
				}
			}
		}
	}

}

static void cffts3(int is,
		int d1,
		int d2,
		int d3,
		void* pointer_x,
		void* pointer_xout,
		dcomplex y1[][FFTBLOCKPAD],
		dcomplex y2[][FFTBLOCKPAD]){
	dcomplex (*x)[NY][NX] = (dcomplex(*)[NY][NX])pointer_x;
	dcomplex (*xout)[NY][NX] = (dcomplex(*)[NY][NX])pointer_xout;

	int logd3;
	int i, j, k, ii;

	logd3 = ilog2(d3);

	#pragma omp for
	for(j=0; j<d2; j++){
		for(ii=0; ii<=d1-FFTBLOCK; ii+=FFTBLOCK){
			for(k=0; k<d3; k++){
				for(i=0; i<FFTBLOCK; i++){
					y1[k][i] = x[k][j][i+ii];
				}
			}
			cfftz(is, logd3, d3, y1, y2);
			for(k=0; k<d3; k++){
				for(i=0; i<FFTBLOCK; i++){
					xout[k][j][i+ii] = y1[k][i];
				}
			}
		}
	}

}

/*
 * ---------------------------------------------------------------------
 * computes NY N-point complex-to-complex FFTs of X using an algorithm due
 * to swarztrauber. X is both the input and the output array, while Y is a 
 * scratch array. it is assumed that N = 2^M. before calling CFFTZ to 
 * perform FFTs, the array U must be initialized by calling CFFTZ with is 
 * set to 0 and M set to MX, where MX is the maximum value of M for any 
 * subsequent call.
 * ---------------------------------------------------------------------
 */
static void cfftz(int is,
		int m,
		int n,
		dcomplex x[][FFTBLOCKPAD],
		dcomplex y[][FFTBLOCKPAD]){
	int i,j,l,mx;

	/*
	 * ---------------------------------------------------------------------
	 * check if input parameters are invalid.
	 * ---------------------------------------------------------------------
	 */
	mx = (int)(u[0].real);
	if((is != 1 && is != -1) || m < 1 || m > mx){
		printf("CFFTZ: Either U has not been initialized, or else\n"    
				"one of the input parameters is invalid%5d%5d%5d\n", is, m, mx);
		exit(EXIT_FAILURE); 
	}

	/*
	 * ---------------------------------------------------------------------
	 * perform one variant of the Stockham FFT.
	 * ---------------------------------------------------------------------
	 */
	for(l=1; l<=m; l+=2){
		fftz2(is, l, m, n, FFTBLOCK, FFTBLOCKPAD, u, x, y);
		if(l==m){
			/*
			 * ---------------------------------------------------------------------
			 * copy Y to X.
			 * ---------------------------------------------------------------------
			 */
			for(j=0; j<n; j++){
				for(i=0; i<FFTBLOCK; i++){
					x[j][i] = y[j][i];
				}
			}
			break;
		}
		fftz2(is, l + 1, m, n, FFTBLOCK, FFTBLOCKPAD, u, y, x);
	}	
}

/*
 * compute function from local (i,j,k) to ibar^2+jbar^2+kbar^2 
 * for time evolution exponent. 
 */
static void compute_indexmap(void* pointer_twiddle,
		int d1,
		int d2,
		int d3){
	double (*twiddle)[NY][NX] = (double(*)[NY][NX])pointer_twiddle;

	int i, j, k, kk, kk2, jj, kj2, ii;
	double ap;

	/*
	 * ---------------------------------------------------------------------
	 * basically we want to convert the fortran indices 
	 * 1 2 3 4 5 6 7 8 
	 * to 
	 * 0 1 2 3 -4 -3 -2 -1
	 * the following magic formula does the trick:
	 * mod(i-1+n/2, n) - n/2
	 * ---------------------------------------------------------------------
	 */
	ap = - 4.0 * ALPHA * PI * PI;
	#pragma omp parallel for private(i,j,kk,kk2,jj,kj2,ii)
	for(k=0; k<d3; k++){
		kk = ((k+NZ/2) % NZ) - NZ/2;
		kk2 = kk*kk;
		for(j=0; j<d2; j++){
			jj = ((j+NY/2) % NY) - NY/2;
			kj2 = jj*jj+kk2;
			for(i=0; i<d1; i++){
				ii = ((i+NX/2) % NX) - NX/2;
				twiddle[k][j][i] = exp(ap*(double)(ii*ii+kj2));
			}
		}
	}
}

/*
 * ---------------------------------------------------------------------
 * fill in array u0 with initial conditions from 
 * random number generator 
 * ---------------------------------------------------------------------
 */
static void compute_initial_conditions(void* pointer_u0,
		int d1,
		int d2,
		int d3){
	dcomplex (*u0)[NY][NX] = (dcomplex(*)[NY][NX])pointer_u0;

	int k, j;
	double x0, start, an, starts[NZ];
	start = SEED;

	/*
	 * ---------------------------------------------------------------------
	 * jump to the starting element for our first plane.
	 * ---------------------------------------------------------------------
	 */
	ipow46(A, 0, &an);
	randlc(&start, an);
	ipow46(A, 2*NX*NY, &an);

	starts[0] = start;
	for(int k=1; k<dims[2]; k++){
		randlc(&start, an);
		starts[k] = start;
	}

	/*
	 * ---------------------------------------------------------------------
	 * go through by z planes filling in one square at a time.
	 * ---------------------------------------------------------------------
	 */
	#pragma omp parallel for private(k,j,x0)
	for(k=0; k<dims[2]; k++){
		x0 = starts[k];
		for(j=0; j<dims[1]; j++){			
			vranlc(2*NX, &x0, A, (double*)&u0[k][j][0]);
		}
	}
}

/*
 * ---------------------------------------------------------------------
 * evolve u0 -> u1 (t time steps) in fourier space
 * ---------------------------------------------------------------------
 */
static void evolve(void* pointer_u0,
		void* pointer_u1,
		void* pointer_twiddle,
		int d1,
		int d2,
		int d3){
	dcomplex (*u0)[NY][NX] = (dcomplex(*)[NY][NX])pointer_u0;
	dcomplex (*u1)[NY][NX] = (dcomplex(*)[NY][NX])pointer_u1;
	double (*twiddle)[NY][NX] = (double(*)[NY][NX])pointer_twiddle;

	int i, j, k;
	#pragma omp for
	for(k=0; k<d3; k++){
		for(j=0; j<d2; j++){
			for(i=0; i<d1; i++){
				u0[k][j][i] = dcomplex_mul2(u0[k][j][i], twiddle[k][j][i]);
				u1[k][j][i] = u0[k][j][i];
			}
		}
	}
}

static void fft(int dir,
		void* pointer_x1,
		void* pointer_x2){
	dcomplex y1[MAXDIM*FFTBLOCKPAD];
	dcomplex y2[MAXDIM*FFTBLOCKPAD];

	/*
	 * ---------------------------------------------------------------------
	 * note: args x1, x2 must be different arrays
	 * note: args for cfftsx are (direction, layout, xin, xout, scratch)
	 * xin/xout may be the same and it can be somewhat faster
	 * if they are
	 * ---------------------------------------------------------------------
	 */
	if(dir==1){
		cffts1(1, dims[0], dims[1], dims[2], pointer_x1, pointer_x1,
				(dcomplex(*)[FFTBLOCKPAD])(void*)y1,
				(dcomplex(*)[FFTBLOCKPAD])(void*)y2);
		cffts2(1, dims[0], dims[1], dims[2], pointer_x1, pointer_x1,
				(dcomplex(*)[FFTBLOCKPAD])(void*)y1,
				(dcomplex(*)[FFTBLOCKPAD])(void*)y2);
		cffts3(1, dims[0], dims[1], dims[2], pointer_x1, pointer_x2,
				(dcomplex(*)[FFTBLOCKPAD])(void*)y1,
				(dcomplex(*)[FFTBLOCKPAD])(void*)y2);
	}else{
		cffts3(-1, dims[0], dims[1], dims[2], pointer_x1, pointer_x1,
				(dcomplex(*)[FFTBLOCKPAD])(void*)y1,
				(dcomplex(*)[FFTBLOCKPAD])(void*)y2);
		cffts2(-1, dims[0], dims[1], dims[2], pointer_x1, pointer_x1,
				(dcomplex(*)[FFTBLOCKPAD])(void*)y1,
				(dcomplex(*)[FFTBLOCKPAD])(void*)y2);
		cffts1(-1, dims[0], dims[1], dims[2], pointer_x1, pointer_x2,
				(dcomplex(*)[FFTBLOCKPAD])(void*)y1,
				(dcomplex(*)[FFTBLOCKPAD])(void*)y2);
	}
}

/*
 * ---------------------------------------------------------------------
 * compute the roots-of-unity array that will be used for subsequent FFTs. 
 * ---------------------------------------------------------------------
 */
static void fft_init(int n){
	int m,nu,ku,i,j,ln;
	double t, ti;
	/*
	 * ---------------------------------------------------------------------
	 * initialize the U array with sines and cosines in a manner that permits
	 * stride one access at each FFT iteration.
	 * ---------------------------------------------------------------------	
	 */	
	nu = n;
	m = ilog2(n);
	u[0] = dcomplex_create((double)m, 0.0);
	ku = 2;
	ln = 1;

	for(j=1; j<=m; j++){
		t = PI / ln;

		for(i=0; i<=ln-1; i++){
			ti = i * t;
			u[i+ku-1] = dcomplex_create(cos(ti), sin(ti));
		}

		ku = ku + ln;
		ln = 2 * ln;
	}
}

/*
 * ---------------------------------------------------------------------
 * performs the l-th iteration of the second variant of the stockham FFT
 * ---------------------------------------------------------------------
 */
static void fftz2(int is,
		int l,
		int m,
		int n,
		int ny,
		int ny1,
		dcomplex u[],
		dcomplex x[][FFTBLOCKPAD],
		dcomplex y[][FFTBLOCKPAD]){
	int k,n1,li,lj,lk,ku,i,j,i11,i12,i21,i22;
	dcomplex u1,x11,x21;

	/*
	 * ---------------------------------------------------------------------
	 * set initial parameters.
	 * ---------------------------------------------------------------------
	 */
	n1 = n / 2;
	lk = 1 << (l - 1);
	li = 1 << (m - l);
	lj = 2 * lk;
	ku = li;

	for(i=0; i<=li-1; i++){
		i11 = i * lk;
		i12 = i11 + n1;
		i21 = i * lj;
		i22 = i21 + lk;
		if(is>=1){
			u1 = u[ku+i];
		}else{
			u1 = dconjg(u[ku+i]);
		}

		/*
		 * ---------------------------------------------------------------------
		 * this loop is vectorizable.
		 * ---------------------------------------------------------------------
		 */
		for(k=0; k<=lk-1; k++){
			for(j=0; j<ny; j++){
				x11 = x[i11+k][j];
				x21 = x[i12+k][j];
				y[i21+k][j] = dcomplex_add(x11, x21);
				y[i22+k][j] = dcomplex_mul(u1, dcomplex_sub(x11, x21));
			}
		}
	}
}

static int ilog2(int n){
	int nn, lg;
	if(n==1){return 0;}
	lg = 1;
	nn = 2;
	while(nn<n){
		nn*=2;
		lg+=1;
	}
	return lg;
}

/*
 * ---------------------------------------------------------------------
 * touch all the big data
 * ---------------------------------------------------------------------  
 */
static void init_ui(void* pointer_u0,
		void* pointer_u1,
		void* pointer_twiddle,
		int d1,
		int d2,
		int d3){
	dcomplex (*u0)[NY][NX] = (dcomplex(*)[NY][NX])pointer_u0;
	dcomplex (*u1)[NY][NX] = (dcomplex(*)[NY][NX])pointer_u1;
	double (*twiddle)[NY][NX] = (double(*)[NY][NX])pointer_twiddle;

	int i, j, k;
	#pragma omp parallel for private(i,j,k)
	for(k=0; k<d3; k++){
		for(j=0; j<d2; j++){
			for(i=0; i<d1; i++){
				u0[k][j][i] = dcomplex_create(0.0, 0.0);
				u1[k][j][i] = dcomplex_create(0.0, 0.0);
				twiddle[k][j][i] = 0.0;
			}
		}
	}
}

/*
 * ---------------------------------------------------------------------
 * compute a^exponent mod 2^46
 * ---------------------------------------------------------------------
 */
static void ipow46(double a,
		int exponent,
		double* result){
	double q, r;
	int n, n2;

	/*
	 * --------------------------------------------------------------------
	 * use
	 * a^n = a^(n/2)*a^(n/2) if n even else
	 * a^n = a*a^(n-1)       if n odd
	 * -------------------------------------------------------------------
	 */
	*result = 1;
	if(exponent==0){return;}
	q = a;
	r = 1;
	n = exponent;

	while(n > 1){
		n2 = n/2;
		if(n2*2==n){
			randlc(&q, q);
			n = n2;
		}else{
			randlc(&r, q);
			n = n-1;
		}
	}
	randlc(&r, q);
	*result = r;
}
