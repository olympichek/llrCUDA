/**************************************************************
 *
 *	gwypnum.c
 *
 *  Modulo k*2^n+/-1 DWFFT multiplications and squarings prototype
 *  source code, GPU version ; must be used linked with cuda and cufft
 *  libraries (versions 8.0.44 work fine!).
 *  This code is built from the one included in Shoichiro Yamada's
 *  llrcuda.0.931 that was released the 07/11/2015 on
 *  www.mersenneforum.org, and, indeed, also from llrp version 3.8.1
 *  (the last released portable LLR version).
 * 
 *  Thanks to the nice work of Shoichiro, it was not difficult for me to
 *  extend his code to rational bases DWT, and also, to generic modular
 *  reduction ; that is what is done here!
 * 
 *  This code is fully C and C++ written, no Assembler code.
 *  Large numbers (at least 1 mega digits) benefit more from the GPU
 *  parallelism, but this program may also be used on smaller positive
 *  results for verification...
 * 
 *  Below is a bit of history of the llrp program :
 *  Nothing original here ; my goal was to have a code portable on any
 *  system having a C / C++ compiler,
 *  and, indeed a processor with a sufficiently powerful floating point
 *  unit!
 *  First : 14/09/2005 : uses George Woltman's 1/k IBDWT method and
 *  cyclic (c = -1), negacyclic (c = +1) real convolutions.
 *  Drawback in the negacyclic case : k must be small due to the 1/cos
 *  factor in the inverse DFFT.
 *  (this factor becomes large near the middle of the FFT array!)
 *  Updates:
 *  14/04/2008 : Full complex, half length convolution used when
 *  computing modulo k*2^n+1 .
 *  Nov. 2010 : zero-padded FFT implemented for k's up to 45 bits large.
 *  Dec. 2010 : generic modular reduction implemented for k's larger
 *  than 45 bits or for general form moduli.
 *  May. 2011 : This code must be linked with  Matteo Frigo and Steven
 *  G. Johnson's FFTW library.
 *  Thanks to this FFTW usage, a power of two FFT length is no more
 *  required.
 *  January 2018 : In this GPU version , CUFFT 8.0.44 is used for all
 *  Fourier transforms.
 *  Jean Penne  02/01/2018, E-mail : jpenne@free.fr
 *
 **************************************************************/

/* Include Files */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "fftw3.h"
#if defined (__linux__) || defined (__FreeBSD__) || defined (__APPLE__)
#include <sys/time.h>
#define _timeb		timeb
#define _ftime		ftime
#else
#include <time.h>
#endif
#include "giants.h"
#include "gwdbldbl.h"
#include "gwypnum.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include "cuda_safecalls.h"


/* definitions */

#define BITS 16
#ifndef _MIPSBUILD_
#define MAXBITSPERDOUBLE (double)26
#define MAXKBITS (double)19
#else
#define MAXBITSPERDOUBLE (double)26
#define MAXKBITS (double)19
#endif

/* The maximum value k * mulbyconst that can be in a zero pad FFT.  Larger */
/* values must use generic modular reduction. */

#if 0 //defined (__linux__) || defined (__FreeBSD__) || defined (__APPLE__)
#define MAX_ZEROPAD_K	35184372088831.0	/* 45-bit k's seem to be OK. */
#else
#define MAX_ZEROPAD_K	68719476735.0		/* 36-bit k's seem to be OK. */
#endif

#define EB 10	// Extra bits of precision for generic reduction

#ifndef WIN32
#define LINE_FEED "\r"
#elif defined (_CONSOLE)
#define LINE_FEED "\r"
#else
#define LINE_FEED "\n"
#endif


/* global variables */

double gwyptimers[2*NBTIMERS] = {0.0};		/* Up to NBTIMERS separate timers */

int     cufftonly = FALSE;
int	E_CHK = 0;   // JP 20/06/17
int	CUMULATIVE_TIMING = 0;
int     error_log = 0;   // JP 20/06/17
double	MAXERR;
giant	gmodulus = NULL;
giant	grecip = NULL;
giant	gtmp = NULL;
gwypnum	GWP_RANDOM = NULL;
gwypnum	modulus = NULL;
gwypnum	recip = NULL;
gwypnum	gwyptmp = NULL;
void	(*printfunction)(char*) = NULL;
void	(*screen_output)(char*) = NULL;
void	(*both_output)(char*) = NULL;

//int	plus = 0, compl = 0, zp = 0, generic = 0, zcomplex = 0, verbose = 0;
int	plus = 0;
int	compl2 = 0;
int	zp = 0;
int	generic = 0;
int	zcomplex = TRUE;
int	verbose = 0;

double	*cn=NULL, *sn=NULL, *cnp=NULL, *snp=NULL, *two_to_phi=NULL, *two_to_minusphi=NULL,
		*w=NULL, *permutedx=NULL, *permutedy=NULL, *invlimit=NULL, *flimit=NULL,
		*hlimit=NULL, *limitbv=NULL;
double 	r, high, low, highinv, lowinv, last, lastinv, addinvalue, wrapfactor;
double	BIGVAL, kg, SMALLMULCONST = 1.0, MAXMULCONST = 1.0;
double	ttmp, avg_num_b_per_word;
double		*xin=NULL, *yin=NULL;
fftw_complex	*cxin=NULL, *cyin=NULL, *cxout=NULL, *cyout=NULL;
double	*cuda_xin=NULL;
double	*cuda_cxin=NULL;
double	*cuda_cxout=NULL;
double	*cuda_yin=NULL;
double	*cuda_cyin=NULL;
double	*cuda_cyout=NULL;
double *cuda_x=NULL;
double *cuda_y=NULL;
double *cuda_m=NULL;
double *cuda_r=NULL;
double *cuda_cm=NULL;
double *cuda_cr=NULL;
double *cuda_tmp=NULL;
double *cuda_tmp_g=NULL;
double *cuda_two_to_phi=NULL;
double *cuda_two_to_minusphi=NULL;
double *cuda_cnp=NULL;
double *cuda_snp=NULL;
double *g_limitbv=NULL,*g_invlimit=NULL,*g_carry=NULL;
float  *g_err=NULL;
float  *l_err=NULL;
double	BIGVAL2 = 6755399441055744.0; 
int     g_fftlen = 0;

fftw_plan		fwpx, fwpy, bwpx, bwpy;
cufftHandle 		cuda_fwpx;
cufftHandle 		cuda_bwpx;
unsigned long 	ng, bit_length, zerowordslow, zerowordshigh; 
int 	 *ip=NULL, *permute=NULL, *fftbase=NULL, addinindex, wrapindex;
int		FFTLEN = 0, debug = 0, MULBYCONST = 0, FFTINC = 0;
int s_FFTLEN = 0; //cuda
char	gwypbuf[256];

// Variables used for modular reduction in zero-padded mode

unsigned long temp = 0, rem, hwcount, lwcount, hwoffset, bits;
int		inc;
double	mult, invmult, shift, limit_high, limit_inverse_high, limit_high_bigval;
double	*scr=NULL, *scral=NULL;


/**************************************************************
 *
 *	Functions
 *
 **************************************************************/

/* rint is not ANSI compatible, so we need a definition for 
 * WIN32 and other platforms with rint.
 */

double
RINT(double x)
{
	return floor(x + 0.5);
}

#define RINTP(x) ((x)-BIGVAL)+BIGVAL

// Allocation routine

gwypnum gwypalloc()
{
	return((gwypnum)malloc(FFTLEN*sizeof(double)));
}

// Utilities

void gwyptrace (int n) {
	printfunction = (verbose)? both_output : screen_output;
	sprintf (gwypbuf, "OK until number %d\n", n);
	if (printfunction != NULL)
		(*printfunction)(gwypbuf);
}

void gwypclearline (int size) {
	char buf[256];
	int i;

	for (i=0; i<256; i++)
		buf[i] = '\0';
	for (i=0; i<size; i++)
		buf[i] = ' ';
	buf[size-1] = '\r';
#if !defined(WIN32) || defined(_CONSOLE)
	printf("%s", buf);
#endif
}


/* Routines used to time code chunks */

void gwypclear_timers () {
	int	i;
	for (i = 0; i < 2*NBTIMERS; i++) gwyptimers[i] = 0.0;
}

void gwypclear_timer (
	int	i)
{
	gwyptimers[i] = 0.0;
}

void gwypstart_timer ( 
	int	i) 
{ 
	struct _timeb timeval; 
	if (i >= NBTIMERS)
		return;
	if (gwyptimers[i+NBTIMERS] != 0.0)			// to avoid double start...
		return;
/*	if (HIGH_RES_TIMER) { 
		gwyptimers[i] -= getHighResTimer (); 
	} else { */
		_ftime (&timeval); 
		gwyptimers[i] -= (double) timeval.time * 1000.0 + timeval.millitm; 
//	} 
	gwyptimers[i+NBTIMERS] = 1.0;				// to show that gwyptimers[i] is already started
} 
 
void gwypend_timer ( 
	int	i) 
{ 
	struct _timeb timeval; 
	if (i >= NBTIMERS)
		return;
	if (gwyptimers[i+NBTIMERS] == 0.0)			// to avoid double end...
		return;
/*	if (HIGH_RES_TIMER) { 
		gwyptimers[i] += getHighResTimer (); 
	} else { */
		_ftime (&timeval); 
		gwyptimers[i] += (double) timeval.time * 1000.0 + timeval.millitm; 
//	} 
	gwyptimers[i+NBTIMERS] = 0.0;				// to show that gwyptimers[i] is ended
} 
 
void gwypdivide_timer (
	int	i,
	int	j)
{
	gwyptimers[i] = gwyptimers[i] / j;
}

double gwyptimer_value ( 
	int	i) 
{ 
/*	if (HIGH_RES_TIMER) 
		return (gwyptimers[i] / getHighResTimerFrequency ()); 
	else */
		return (gwyptimers[i] / 1000.0); 
} 
 
void gwypprint_timer (
	int	i,
	int	flags)
{ 
	char	buf[40]; 
	double	t; 
 
	t = gwyptimer_value (i); 
	if (flags & TIMER_NL)
		if (t >= 1.0)  
			sprintf (buf, "%.3f sec."LINE_FEED"", t); 
		else 
			sprintf (buf, "%.3f ms."LINE_FEED"", t * 1000.0);
	else
		if (t >= 1.0)  
			sprintf (buf, "%.3f sec.", t); 
		else 
			sprintf (buf, "%.3f ms.", t * 1000.0);


	printfunction = screen_output;
	if (printfunction != NULL)
		(*printfunction)(buf);


	if (flags & TIMER_CLR) gwyptimers[i] = 0.0; 
	if ((flags & TIMER_OPT_CLR) && !CUMULATIVE_TIMING) gwyptimers[i] = 0.0; 
} 

void gwypwrite_timer (		// JP 23/11/07
	char* buf,
	int	i, 
	int	flags) 
{ 
	double	t; 
 
	t = gwyptimer_value (i); 
	if (flags & TIMER_NL)
		if (t >= 1.0)  
			sprintf (buf, "%.3f sec.\n", t); 
		else 
			sprintf (buf, "%.3f ms.\n", t * 1000.0);
	else
		if (t >= 1.0)  
			sprintf (buf, "%.3f sec.", t); 
		else 
			sprintf (buf, "%.3f ms.", t * 1000.0);
 
 
	if (flags & TIMER_CLR) gwyptimers[i] = 0.0; 
	if ((flags & TIMER_OPT_CLR) && !CUMULATIVE_TIMING) gwyptimers[i] = 0.0; 
} 

void
print(
	double 	*x,
	int  	N
)
{
	int  	printed = 0;

	printfunction = (verbose)? both_output : screen_output;

	while (N >= 0)
	{
		if ((x[N]==0) && (!printed))
			goto decr;
		sprintf(gwypbuf, "%g  ",x[N]);
		if (printfunction != NULL)
			(*printfunction)(gwypbuf);
		printed=1;
decr:		N--;
	}
	sprintf(gwypbuf, "||\n");
	if (printfunction != NULL)
		(*printfunction)(gwypbuf);
}

/* Routine that copy a gwypnum from */
/* source to dest while zeroing some lower FFT words */

__global__ void
cuda_gwypcopyzero_kernel (
    double *s,
    double *d,
    unsigned long n,
    unsigned long len)
{
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	register double zero = 0.0;
	register unsigned long 	i;
        
        i = threadID;
        if (i<n)
            d[i] = zero;
        else if (i<len)
            d[i] = s[i];
}

void gwypcopyzero (
	gwypnum	s,
	gwypnum	d,
	unsigned long n)
{
	register double zero = 0.0;
	register double *sptr = s + n;
	register double *dptr = d;
	register double *maxptr;

	if (debug)
		gwypstart_timer (4);
	maxptr = d + n;
	while (dptr < maxptr)
		*dptr++ = zero;

	maxptr = d + FFTLEN;
	while (dptr < maxptr)
		*dptr++ = *sptr++;
	if (debug)
		gwypend_timer (4);
}

/* Set a gwypnum to zero */

void gwypzero (gwypnum s) {
	long j;

	for(j=0; j<FFTLEN; ++j)
		s[j] = 0;
	return;
}

/* Routine that zero some high words in a gwypnum */

__global__ void
cuda_gwypsetzero_kernel (
    double *s,
    unsigned long n,
    unsigned long len)
{
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	register double zero = 0.0;
	register unsigned long 	i;
        
        i = threadID+len-n;
        if (i<len)
            s[i] = zero;
}

void gwypsetzero (
	gwypnum s,
	unsigned long n)
{
	register double zero = 0.0;
	register double *sptr = s + FFTLEN - n;
	register double *maxptr = s + FFTLEN;

	if (debug)
		gwypstart_timer (4);
while (sptr < maxptr)
	*sptr++ = zero;
	if (debug)
		gwypend_timer (4);
}

// User side large integers arithmetic operations

__global__ void
cuda_gwypcopy_kernel (
    double *s,
    double *d,
    int n)
{
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	register int 	i;
        
        i = threadID;
        if (i<n)
            d[i] = s[i];
}

void gwypcopy (
	gwypnum s,
	gwypnum d)
{
	int i;

	for (i=0; i<FFTLEN; i++)
		d[i] = s[i];
}

__global__ void
cuda_gwypaddquick_kernel (
    double *s,
    double *d,
    int n)
{
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	register int 	i;
        
        i = threadID;
        if (i<n)
            d[i] += s[i];
}

void gwypaddquick (
	gwypnum s,
	gwypnum d)
{
	int i;

	for (i=0; i<FFTLEN; i++)
		d[i] += s[i];
}

__global__ void
cuda_gwypsubquick_kernel (
    double *s,
    double *d,
    int n)
{
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	register int 	i;
        
        i = threadID;
        if (i<n)
            d[i] -= s[i];
}

void gwypsubquick (
	gwypnum s,
	gwypnum d)
{
	int i;

	for (i=0; i<FFTLEN; i++)
		d[i] -= s[i];
}

// These functions do the relevant dyadic multiplications or squarings on Fourier transformed data

__global__ void
cuda_mul_complex_kernel(
	double 	*a,
	double 	*b,
	int 			n
)
{
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	register int 	k;
	register double Reb;

	k=threadID;
	if (k<n) {
            k=k*2;
            Reb = a[k]*b[k]-a[k+1]*b[k+1];
            b[k+1] = a[k+1]*b[k]+a[k]*b[k+1];
            b[k] = Reb;
	}
}

void
_mul_complex(
	fftw_complex 	*a,
	fftw_complex 	*b,
	int 			n
)
{
	register int 	k;
	register double Reb;

	for (k=0; k<n; k++) {
            Reb = a[k][0]*b[k][0]-a[k][1]*b[k][1];
            b[k][1] = a[k][1]*b[k][0]+a[k][0]*b[k][1];
            b[k][0] = Reb;
	}
}

__global__ void
cuda_square_complex_kernel(
	double *b,
	int n
)
{
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	register int k;
	register double Reb;
	k=threadID;
	if(k<n)
	{
            k=k*2;
            Reb = b[k]*b[k]-b[k+1]*b[k+1];
            b[k+1] = 2*b[k+1]*b[k];
            b[k] = Reb;
	}

}

void
cuda_square_complex(
	fftw_complex *b,
	int n
)
{
	register int k;
	register double Reb;

	for (k=0; k<n; k++) {
            Reb = b[k][0]*b[k][0]-b[k][1]*b[k][1];
            b[k][1] = 2*b[k][1]*b[k][0];
            b[k][0] = Reb;
	}

}



void
_square_complex(
	fftw_complex *b,
	int n
)
{
	register int k;
	register double Reb;

	for (k=0; k<n; k++) {
            Reb = b[k][0]*b[k][0]-b[k][1]*b[k][1];
            b[k][1] = 2*b[k][1]*b[k][0];
            b[k][0] = Reb;
	}

}

// These function do the general multiplication or squaring of large integers, using DFFT

__global__ void
cuda_cnp_m_snp_kernel(
	double *cxin,
	double *cnp,
	double *snp,
	int n
)
{
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	register int k,j;
	register double ReX;
	j=threadID;
	if(j<n)
	{
            k=j*2;
            ReX = cnp[j]*cxin[k]-snp[j]*cxin[k+1];
            cxin[k+1] = cnp[j]*cxin[k+1]+snp[j]*cxin[k];
            cxin[k] = ReX;
	}
}

__global__ void
cuda_cnp_p_snp_kernel(
	double *cxin,
	double *cnp,
	double *snp,
	int n
)
{
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	register int k,j;
	register double ReX;
	j=threadID;
	if(j<n)
	{
            k=j*2;
            ReX = cnp[j]*cxin[k]+snp[j]*cxin[k+1];
            cxin[k+1] = cnp[j]*cxin[k+1]-snp[j]*cxin[k];
            cxin[k] = ReX;
	}
}

// These function do the general multiplication or squaring of large integers, using DFFT

void cuda_fftwsquare_g (
	int size)
{
//	register int j;
//	register double ReX;

        if (compl2) 	
	{    // Full complex, half size DWFFT
		cuda_cnp_m_snp_kernel<<<(size/2+127)/128,128>>>(cuda_cxin,cuda_cnp,cuda_snp, size/2);
		cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cxin,(cufftDoubleComplex *)cuda_cxout,CUFFT_FORWARD));
		cuda_square_complex_kernel<<<(size/2+127)/128,128>>>(cuda_cxout, size/2);
		cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cxout,(cufftDoubleComplex *)cuda_cxin,CUFFT_INVERSE));
		cuda_cnp_p_snp_kernel<<<(size/2+127)/128,128>>>(cuda_cxin,cuda_cnp,cuda_snp, size/2);
	}
	else  // Real to complex, full size DWFFT
	{
		cufftSafeCall(cufftExecD2Z(cuda_fwpx,(cufftDoubleReal *)cuda_xin,(cufftDoubleComplex *)cuda_cxout));
		cuda_square_complex_kernel<<<(size/2+1+127)/128,128>>>(cuda_cxout, size/2+1);
		cufftSafeCall(cufftExecZ2D(cuda_bwpx,(cufftDoubleComplex *)cuda_cxout,(cufftDoubleReal *)cuda_xin));
	}	
}

// These function do the general multiplication or squaring of large integers, using DFFT

void fftwsquare_g (
	int size)
{
	register int j;
	register double ReX;

        if (compl2) {    // Full complex, half size DWFFT
                                        // Multiply cxin by exp(i*j*pi/size) to prepare a right-angle convolution
                for (j=0; j<size/2; j++) {
                        ReX = cnp[j]*cxin[j][0]-snp[j]*cxin[j][1];
                        cxin[j][1] = cnp[j]*cxin[j][1]+snp[j]*cxin[j][0];
                        cxin[j][0] = ReX;
                }

        }

	//fftw_execute (fwpx);	// Execute the relevant forward DFFT
        if (compl2)     
	{
		cutilSafeCall(cudaMemcpy(cuda_cxin,cxin,sizeof(double)*size,cudaMemcpyHostToDevice));
		cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cxin,(cufftDoubleComplex *)cuda_cxout,CUFFT_FORWARD));
		cutilSafeCall(cudaMemcpy(cxout,cuda_cxout,sizeof(double)*size,cudaMemcpyDeviceToHost));
	}
	else
	{
		cutilSafeCall(cudaMemcpy(cuda_xin,xin,sizeof(double)*size,cudaMemcpyHostToDevice));
		cufftSafeCall(cufftExecD2Z(cuda_fwpx,(cufftDoubleReal *)cuda_xin,(cufftDoubleComplex *)cuda_cxout));
		cutilSafeCall(cudaMemcpy(cxout,cuda_cxout,sizeof(double)*2*(size/2+1),cudaMemcpyDeviceToHost));
	}

        if (compl2)                              // Compute the relevant Dyadic squaring
                _square_complex (cxout, size/2);
        else
                _square_complex(cxout, size/2+1);

	//fftw_execute (bwpx);	// Execute the relevant backward DFFT

        if (compl2)     
	{
		cutilSafeCall(cudaMemcpy(cuda_cxout,cxout,sizeof(double)*size,cudaMemcpyHostToDevice));
		cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cxout,(cufftDoubleComplex *)cuda_cxin,CUFFT_INVERSE));
		cutilSafeCall(cudaMemcpy(cxin,cuda_cxin,sizeof(double)*size,cudaMemcpyDeviceToHost));
	}
	else
	{
		cutilSafeCall(cudaMemcpy(cuda_cxout,cxout,sizeof(double)*2*(size/2+1),cudaMemcpyHostToDevice));
		cufftSafeCall(cufftExecZ2D(cuda_bwpx,(cufftDoubleComplex *)cuda_cxout,(cufftDoubleReal *)cuda_xin));
		cutilSafeCall(cudaMemcpy(xin,cuda_xin,sizeof(double)*size,cudaMemcpyDeviceToHost));
	}

        if (compl2) {                    // Full complex, half size DWFFT
                                // Multiply cxin by exp(-i*j*pi/size) to complete the right-angle convolution
                for (j=0; j<size/2; ++j)
                {
                        ReX = cnp[j]*cxin[j][0]+snp[j]*cxin[j][1];
                        cxin[j][1] = cnp[j]*cxin[j][1]-snp[j]*cxin[j][0];
                        cxin[j][0] = ReX;
                }
        }

}

void cuda_fftwmul_g (
	int size)
{
//	register int j;
//	register double ReX;

        if (compl2) 	
	{    // Full complex, half size DWFFT
		cuda_cnp_m_snp_kernel<<<(size/2+127)/128,128>>>(cuda_cxin,cuda_cnp,cuda_snp, size/2);
		cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cxin,(cufftDoubleComplex *)cuda_cxout,CUFFT_FORWARD));
		cuda_cnp_m_snp_kernel<<<(size/2+127)/128,128>>>(cuda_cyin,cuda_cnp,cuda_snp, size/2);
		cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cyin,(cufftDoubleComplex *)cuda_cyout,CUFFT_FORWARD));
		cuda_mul_complex_kernel<<<(size/2+127)/128,128>>>(cuda_cxout, cuda_cyout, size/2);
//		cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cxout,(cufftDoubleComplex *)cuda_cxin,CUFFT_INVERSE));
//		cuda_cnp_p_snp_kernel<<<(size/2+127)/128,128>>>(cuda_cxin,cuda_cnp,cuda_snp, size/2);
		cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cyout,(cufftDoubleComplex *)cuda_cyin,CUFFT_INVERSE));
		cuda_cnp_p_snp_kernel<<<(size/2+127)/128,128>>>(cuda_cyin,cuda_cnp,cuda_snp, size/2);
	}
	else
	{
		cufftSafeCall(cufftExecD2Z(cuda_fwpx,(cufftDoubleReal *)cuda_xin,(cufftDoubleComplex *)cuda_cxout));
		cufftSafeCall(cufftExecD2Z(cuda_fwpx,(cufftDoubleReal *)cuda_yin,(cufftDoubleComplex *)cuda_cyout));
		cuda_mul_complex_kernel<<<(size/2+1+127)/128,128>>>(cuda_cxout, cuda_cyout, size/2+1);
//		cufftSafeCall(cufftExecZ2D(cuda_bwpx,(cufftDoubleComplex *)cuda_cxout,(cufftDoubleReal *)cuda_xin));
		cufftSafeCall(cufftExecZ2D(cuda_bwpx,(cufftDoubleComplex *)cuda_cyout,(cufftDoubleReal *)cuda_yin));
	}	
}

void cuda_fftwmulbym_g (
	int size)
{
//	register int j;
//	register double ReX;

        if (compl2) 	
	{    // Full complex, half size DWFFT
//		cuda_cnp_m_snp_kernel<<<(size/2+127)/128,128>>>(cuda_cxin,cuda_cnp,cuda_snp, size/2);
//		cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cxin,(cufftDoubleComplex *)cuda_cxout,CUFFT_FORWARD));
		cuda_cnp_m_snp_kernel<<<(size/2+127)/128,128>>>(cuda_cyin,cuda_cnp,cuda_snp, size/2);
		cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cyin,(cufftDoubleComplex *)cuda_cyout,CUFFT_FORWARD));
		cuda_mul_complex_kernel<<<(size/2+127)/128,128>>>(cuda_cm, cuda_cyout, size/2);
//		cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cxout,(cufftDoubleComplex *)cuda_cxin,CUFFT_INVERSE));
//		cuda_cnp_p_snp_kernel<<<(size/2+127)/128,128>>>(cuda_cxin,cuda_cnp,cuda_snp, size/2);
		cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cyout,(cufftDoubleComplex *)cuda_cyin,CUFFT_INVERSE));
		cuda_cnp_p_snp_kernel<<<(size/2+127)/128,128>>>(cuda_cyin,cuda_cnp,cuda_snp, size/2);
	}
	else
	{
//		cufftSafeCall(cufftExecD2Z(cuda_fwpx,(cufftDoubleReal *)cuda_xin,(cufftDoubleComplex *)cuda_cxout));
		cufftSafeCall(cufftExecD2Z(cuda_fwpx,(cufftDoubleReal *)cuda_yin,(cufftDoubleComplex *)cuda_cyout));
		cuda_mul_complex_kernel<<<(size/2+1+127)/128,128>>>(cuda_cm, cuda_cyout, size/2+1);
//		cufftSafeCall(cufftExecZ2D(cuda_bwpx,(cufftDoubleComplex *)cuda_cxout,(cufftDoubleReal *)cuda_xin));
		cufftSafeCall(cufftExecZ2D(cuda_bwpx,(cufftDoubleComplex *)cuda_cyout,(cufftDoubleReal *)cuda_yin));
	}	
}

void cuda_fftwmulbyr_g (
	int size)
{
//	register int j;
//	register double ReX;

        if (compl2) 	
	{    // Full complex, half size DWFFT
//		cuda_cnp_m_snp_kernel<<<(size/2+127)/128,128>>>(cuda_cxin,cuda_cnp,cuda_snp, size/2);
//		cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cxin,(cufftDoubleComplex *)cuda_cxout,CUFFT_FORWARD));
		cuda_cnp_m_snp_kernel<<<(size/2+127)/128,128>>>(cuda_cyin,cuda_cnp,cuda_snp, size/2);
		cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cyin,(cufftDoubleComplex *)cuda_cyout,CUFFT_FORWARD));
		cuda_mul_complex_kernel<<<(size/2+127)/128,128>>>(cuda_cr, cuda_cyout, size/2);
//		cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cxout,(cufftDoubleComplex *)cuda_cxin,CUFFT_INVERSE));
//		cuda_cnp_p_snp_kernel<<<(size/2+127)/128,128>>>(cuda_cxin,cuda_cnp,cuda_snp, size/2);
		cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cyout,(cufftDoubleComplex *)cuda_cyin,CUFFT_INVERSE));
		cuda_cnp_p_snp_kernel<<<(size/2+127)/128,128>>>(cuda_cyin,cuda_cnp,cuda_snp, size/2);
	}
	else
	{
//		cufftSafeCall(cufftExecD2Z(cuda_fwpx,(cufftDoubleReal *)cuda_xin,(cufftDoubleComplex *)cuda_cxout));
		cufftSafeCall(cufftExecD2Z(cuda_fwpx,(cufftDoubleReal *)cuda_yin,(cufftDoubleComplex *)cuda_cyout));
		cuda_mul_complex_kernel<<<(size/2+1+127)/128,128>>>(cuda_cr, cuda_cyout, size/2+1);
//		cufftSafeCall(cufftExecZ2D(cuda_bwpx,(cufftDoubleComplex *)cuda_cxout,(cufftDoubleReal *)cuda_xin));
		cufftSafeCall(cufftExecZ2D(cuda_bwpx,(cufftDoubleComplex *)cuda_cyout,(cufftDoubleReal *)cuda_yin));
	}	
}

void fftwmul_g (int size) {
	register int j;
	register double ReX, ReY;

        if (compl2) {    // Full complex, half size DWFFT
                        // Multiply cxin and cyin by exp(i*j*pi/size) to prepare a right-angle convolution
                for (j=0; j<size/2; j++) {
                        ReX = cnp[j]*cxin[j][0]-snp[j]*cxin[j][1];
                        cxin[j][1] = cnp[j]*cxin[j][1]+snp[j]*cxin[j][0];
                        cxin[j][0] = ReX;
                        ReY = cnp[j]*cyin[j][0]-snp[j]*cyin[j][1];
                        cyin[j][1] = cnp[j]*cyin[j][1]+snp[j]*cyin[j][0];
                        cyin[j][0] = ReY;
                }
        }
    if (compl2) {
	//fftw_execute (fwpx);	// Execute the two relevant forward DFFTs
	cutilSafeCall(cudaMemcpy(cuda_cxin,cxin,sizeof(double)*size,cudaMemcpyHostToDevice));
	cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cxin,(cufftDoubleComplex *)cuda_cxout,CUFFT_FORWARD));
	cutilSafeCall(cudaMemcpy(cxout,cuda_cxout,sizeof(double)*2*(size/2+1),cudaMemcpyDeviceToHost));
	//fftw_execute (fwpy);
	cutilSafeCall(cudaMemcpy(cuda_cyin,cyin,sizeof(double)*size,cudaMemcpyHostToDevice));
	cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cyin,(cufftDoubleComplex *)cuda_cyout,CUFFT_FORWARD));
	cutilSafeCall(cudaMemcpy(cyout,cuda_cyout,sizeof(double)*2*(size/2+1),cudaMemcpyDeviceToHost));
    }
    else {
	//fftw_execute (fwpx);	// Execute the two relevant forward DFFTs
	cutilSafeCall(cudaMemcpy(cuda_xin,xin,sizeof(double)*size,cudaMemcpyHostToDevice));
	cufftSafeCall(cufftExecD2Z(cuda_fwpx,(cufftDoubleReal *)cuda_xin,(cufftDoubleComplex *)cuda_cxout));
	cutilSafeCall(cudaMemcpy(cxout,cuda_cxout,sizeof(double)*2*(size/2+1),cudaMemcpyDeviceToHost));
	//fftw_execute (fwpy);
	cutilSafeCall(cudaMemcpy(cuda_yin,yin,sizeof(double)*size,cudaMemcpyHostToDevice));
	cufftSafeCall(cufftExecD2Z(cuda_fwpx,(cufftDoubleReal *)cuda_yin,(cufftDoubleComplex *)cuda_cyout));
	cutilSafeCall(cudaMemcpy(cyout,cuda_cyout,sizeof(double)*2*(size/2+1),cudaMemcpyDeviceToHost));
    }

	if (compl2)				// Compute the relevant Dyadic product
		_mul_complex (cxout, cyout, size/2);
	else
		_mul_complex(cxout, cyout, size/2+1);
    if (compl2) {
            cutilSafeCall(cudaMemcpy(cuda_cyout,cyout,sizeof(double)*size,cudaMemcpyHostToDevice));
            cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cyout,(cufftDoubleComplex *)cuda_cyin,CUFFT_INVERSE));
            cutilSafeCall(cudaMemcpy(cyin,cuda_cyin,sizeof(double)*size,cudaMemcpyDeviceToHost));
    }
    else {
	//fftw_execute (bwpy);	// Execute the relevant backward DFFT
	cutilSafeCall(cudaMemcpy(cuda_cyout,cyout,sizeof(double)*2*(size/2+1),cudaMemcpyHostToDevice));
	cufftSafeCall(cufftExecZ2D(cuda_bwpx,(cufftDoubleComplex *)cuda_cyout,(cufftDoubleReal *)cuda_yin));
	cutilSafeCall(cudaMemcpy(yin,cuda_yin,sizeof(double)*size,cudaMemcpyDeviceToHost));
    }

        if (compl2) {                    // Full complex, half size DWFFT
                                // Multiply cyin by exp(-i*j*pi/size) to complete the right-angle convolution
                for (j=0; j<size/2; ++j)
                {
                        ReY = cnp[j]*cyin[j][0]+snp[j]*cyin[j][1];
                        cyin[j][1] = cnp[j]*cyin[j][1]-snp[j]*cyin[j][0];
                        cyin[j][0] = ReY;
                }
        }



}

#define STRIDE_DIM 256
/**************************************************************
 *
 *      Functions
 *
 **************************************************************/
#define BLOCK_DIM 16

// This kernel is optimized to ensure all global reads and writes are coalesced,
// and to avoid bank conflicts in shared memory.  This kernel is up to 11x faster
// than the naive kernel below.  Note that the shared memory array is sized to 
// (BLOCK_DIM+1)*BLOCK_DIM.  This pads each row of the 2D block in shared memory 
// so that bank conflicts do not occur when threads address the array column-wise.
__global__ void transpose(double *odata, double *idata, int width, int height)
{
        __shared__ double block[BLOCK_DIM][BLOCK_DIM+1];

        // read the matrix tile into shared memory
        unsigned int xIndex = blockIdx.x * BLOCK_DIM + threadIdx.x;
        unsigned int yIndex = blockIdx.y * BLOCK_DIM + threadIdx.y;
        if((xIndex < width) && (yIndex < height))
        {
                unsigned int index_in = yIndex * width + xIndex;
                block[threadIdx.y][threadIdx.x] = idata[index_in];
        }
        __syncthreads();

        // write the transposed matrix tile to global memory
        xIndex = blockIdx.y * BLOCK_DIM + threadIdx.x;
        yIndex = blockIdx.x * BLOCK_DIM + threadIdx.y;
        if((xIndex < height) && (yIndex < width))
        {
                unsigned int index_out = yIndex * height + xIndex;
                odata[index_out] = block[threadIdx.x][threadIdx.y];
        }
}

__global__ void mul_const_transpose(double *odata, double *idata, double c, int width, int height)
{
        __shared__ double block[BLOCK_DIM][BLOCK_DIM+1];

        // read the matrix tile into shared memory
        unsigned int xIndex = blockIdx.x * BLOCK_DIM + threadIdx.x;
        unsigned int yIndex = blockIdx.y * BLOCK_DIM + threadIdx.y;
        if((xIndex < width) && (yIndex < height))
        {
                unsigned int index_in = yIndex * width + xIndex;
                block[threadIdx.y][threadIdx.x] = idata[index_in]*c;
        }
        __syncthreads();

        // write the transposed matrix tile to global memory
        xIndex = blockIdx.y * BLOCK_DIM + threadIdx.x;
        yIndex = blockIdx.x * BLOCK_DIM + threadIdx.y;
        if((xIndex < height) && (yIndex < width))
        {
                unsigned int index_out = yIndex * height + xIndex;
                odata[index_out] = block[threadIdx.x][threadIdx.y];
        }
}

__global__ void mul_0_transpose(double *odata, double *idata, double *mul, int width, int height)
{
        __shared__ double block[BLOCK_DIM][BLOCK_DIM+1];

        // read the matrix tile into shared memory
        unsigned int xIndex = blockIdx.x * BLOCK_DIM + threadIdx.x;
        unsigned int yIndex = blockIdx.y * BLOCK_DIM + threadIdx.y;
        if((xIndex < width) && (yIndex < height))
        {
                unsigned int index_in = yIndex * width + xIndex;
                block[threadIdx.y][threadIdx.x] = idata[index_in]*mul[index_in];
        }
        __syncthreads();

        // write the transposed matrix tile to global memory
        xIndex = blockIdx.y * BLOCK_DIM + threadIdx.x;
        yIndex = blockIdx.x * BLOCK_DIM + threadIdx.y;
        if((xIndex < height) && (yIndex < width))
        {
                unsigned int index_out = yIndex * height + xIndex;
                odata[index_out] = block[threadIdx.x][threadIdx.y];
        }
}

__global__ void mul_1_transpose(double *odata, double *idata, double *mul, int width, int height)
{
        __shared__ double block[BLOCK_DIM][BLOCK_DIM+1];

        // read the matrix tile into shared memory
        unsigned int xIndex = blockIdx.x * BLOCK_DIM + threadIdx.x;
        unsigned int yIndex = blockIdx.y * BLOCK_DIM + threadIdx.y;
        if((xIndex < width) && (yIndex < height))
        {
                unsigned int index_in = yIndex * width + xIndex;
                block[threadIdx.y][threadIdx.x] = idata[index_in];
        }
        __syncthreads();

        // write the transposed matrix tile to global memory
        xIndex = blockIdx.y * BLOCK_DIM + threadIdx.x;
        yIndex = blockIdx.x * BLOCK_DIM + threadIdx.y;
        if((xIndex < height) && (yIndex < width))
        {
                unsigned int index_out = yIndex * height + xIndex;
                odata[index_out] = block[threadIdx.x][threadIdx.y]*mul[index_out];
        }
}

//#define IDX(i) ((((i) >> 14) << 14) + (((i) & (128*128-1)) >> 7)  + (((i) & 127) << 7))
#define IDX(i) ((((i) >> 16) << 16) + (((i) & (256*256-1)) >> 8)  + (((i) & 255) << 8))
//#define IDX(i) ((((i) >> 18) << 18) + (((i) & (512*512-1)) >> 9)  + (((i) & 511) << 9))

//      Functions that do the normalization of the result of a multiplication or squaring

__global__ void cuda_inormalize_kernel(     // Used for irrational bases DWFFT
        double          *x,
        int             N,
        int             error_log,
        int             noadd,
        int             nomul,
        double          *g_limitbv,
        double          *g_limitbv2,
        double          *g_invlimit,
        int             STRIDE,
        double          bv,
        double          bv2,
        double          *g_carry,
        float           *g_err,
        int             wrapindex,
        double          wrapfactor,
        int             MULBYCONST,
        double          SMALLMULCONST,
        int             addinindex,
        double          addinvalue,
        int             plus
)
{
	register int    	j;
        register double         xx,zz,tx;
        register double         carry,err,maxerr;
        const int 		threadID = blockIdx.x * blockDim.x + threadIdx.x;
        if(addinindex==0)
	{
                if( threadID==0)
                        if (!noadd) 
				x[IDX(addinindex)] += addinvalue;   // Add the optional small constant
        }
        else
	{
        	if( threadID*STRIDE < N && threadID*STRIDE+STRIDE >= N)
                        if (!noadd) 
				x[IDX(addinindex)] += addinvalue;   // Add the optional small constant
        }
        carry = 0;
        maxerr = 0;
        for (j = threadID*STRIDE; j < threadID*STRIDE+STRIDE && j < N;++j)
        {
                tx = x[IDX(j)];
                if(MULBYCONST && !nomul)
                        tx *= SMALLMULCONST; // Optionaly multiply by a small constant
                xx = (tx + bv) - bv2; // Round to the nearest integer
                if(error_log) 
		{                    // Compute the rounding error if required
                        if(xx<0)
                                err = fabs(-xx + tx);
                        else
                                err = fabs(xx  - tx);
                        if(err > maxerr)
                                maxerr = err;
                }
                xx += carry;
                zz = (xx+g_limitbv[IDX(j)])-g_limitbv2[IDX(j)];
                carry = zz*g_invlimit[IDX(j)]; // Compute the carry on next word
                x[IDX(j)] = xx - zz;           // And the balanced remainder in current word
        }
        g_carry[threadID]=carry;
        if(g_err[threadID] < maxerr)
		g_err[threadID] = maxerr;
}

//      Functions that do the normalization of the result of a multiplication or squaring
__global__ void cuda_inormalize2_kernel(        // Used for irrational bases DWFFT
        double          *x,
        int             N,
        int             error_log,
        int             noadd,
        int             nomul,
        double          *g_limitbv,
        double          *g_limitbv2,
        double          *g_invlimit,
        int             STRIDE,
        double          bv,
        double          bv2,
        double          *g_carry,
        float           *g_err,
        int             wrapindex,
        double          wrapfactor,
        int             MULBYCONST,
        double          SMALLMULCONST,
        int             addinindex,
        double          addinvalue,
        int             plus
)
{
        register int    	j;
        register double         tx,xx,zz;
        register double         carry,carry2;
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

        if(threadID*STRIDE < N && threadID*STRIDE+STRIDE >= N)
        {
                carry = g_carry[threadID];
		carry2 = 0;
                if (carry)
                {
                        j = 0;
                        if (wrapindex)
                                carry2 = carry*wrapfactor;
                        if (plus) 
				carry = -carry;
                        while (carry||carry2)
                        {       
				if (wrapindex && !carry) // Skip already normalized words
                                 	j = wrapindex;
                                tx = x[IDX(j)];
                                xx = tx + carry;
                                if (wrapindex && j==wrapindex) 
				{
                                        xx += carry2;
                                        carry2 = 0.0;
                                }
                                zz = (xx+g_limitbv[IDX(j)])-g_limitbv2[IDX(j)];
                                carry = zz*g_invlimit[IDX(j)]; // Compute the carry on next word
                                x[IDX(j)] = xx - zz; // And the balanced remainder in current word
                                if (++j == N)
                                {
                                        j = 0;
                                        if (wrapindex)
                                                carry2 = carry*wrapfactor;
                                        if (plus) 
						carry = -carry;
                                }
                        }
                }
        }
        else
        {
                carry = g_carry[threadID];
		carry2 = 0;
                j = threadID*STRIDE+STRIDE;
                while (carry)
                {
                        tx = x[IDX(j)];
                        xx = tx + carry;
                        zz = (xx+g_limitbv[IDX(j)])-g_limitbv2[IDX(j)];
                        carry = zz*g_invlimit[IDX(j)];          // Compute the carry on next word
                        x[IDX(j)] = xx - zz;                      // And the balanced remainder in current word
                        j++;
                }
        }
}

//	Functions that do the normalization of the result of a multiplication or squaring

double
inormalize(                     // Used for irrational bases DWFFT
	double 			*x,
	int 			N,
	int				error_log,
	int				noadd,
	int				nomul
)
{
	register int 	j;
	register double 	*px = x, xx, zz, bv = BIGVAL2;    // JP 08/07/17
	register double  	carry = 0.0, carry2 = 0.0, err, maxerr = 0.0;

	if (!noadd)
		x[addinindex] += addinvalue;            // Add the optional small constant
	for (j=0; j<N; ++j)
	{
		if (MULBYCONST && !nomul)
			*px *= SMALLMULCONST;		// Optionaly multiply by a small constant
		xx = (*px + bv) - bv;			// Round to the nearest integer
		if (error_log ) {			// Compute the rounding error if required
			if (xx<0)
          			err = fabs(-xx + *px);  
     			else 
     				err = fabs(xx  - *px);
	 		if (err > maxerr) 
	 			maxerr = err;
     		}

		xx += carry;
		zz = (xx+limitbv[j])-limitbv[j];
		carry = zz*invlimit[j];			// Compute the carry on next word
		*(px++) = xx - zz;			// And the balanced remainder in current word

	}

	if (carry)
	{
		j = 0;
		px = x;
		if (wrapindex)
			carry2 = carry*wrapfactor;
		if (plus)
			carry = -carry;
		while (carry||carry2)
		{	if (wrapindex && !carry) {	// Skip already normalized words
				j = wrapindex;
				px = x + wrapindex;
			}
			xx = *px + carry;
			if (wrapindex && j==wrapindex) {
				xx += carry2;
				carry2 = 0.0;
			}

			zz = (xx+limitbv[j])-limitbv[j];
			carry = zz*invlimit[j];		// Compute the carry on next word
			*(px++) = xx - zz;			// And the balanced remainder in current word

			if (++j == N)
			{
				j = 0;
				px = x;
				if (wrapindex)
					carry2 = carry*wrapfactor;
				if (plus)
					carry = -carry;
			}
		}
	}
	return(maxerr);
}

//      Functions that do the normalization of the result of a multiplication or squaring

__global__ void cuda_rnormalize_kernel(     // Used for rational bases DWFFT
        double          *x,
        int             N,
        int             error_log,
        int             noadd,
        int             nomul,
        double          limitbv,
        double          limitbv2,
        double          invlimit,
        int             STRIDE,
        double          bv,
        double          bv2,
        double          *g_carry,
        float           *g_err,
        int             MULBYCONST,
        double          SMALLMULCONST,
        int             addinindex,
        double          addinvalue,
        int             plus
)
{
	register int    	j;
        register double         xx, zz, tx;
        register double         carry = 0.0, err, maxerr = 0.0;
        const int 		threadID = blockIdx.x * blockDim.x + threadIdx.x;
        if(addinindex==0)
	{
                if( threadID==0)
                        if (!noadd) 
				x[IDX(addinindex)] += addinvalue;   // Add the optional small constant
        }
        else
	{
        	if( threadID*STRIDE < N && threadID*STRIDE+STRIDE >= N)
                        if (!noadd) 
				x[IDX(addinindex)] += addinvalue;   // Add the optional small constant
        }
        carry = 0;
        maxerr = 0;
        for (j = threadID*STRIDE; j < threadID*STRIDE+STRIDE && j < N;++j)
        {
                tx = x[IDX(j)];
                if(MULBYCONST && !nomul)
                        tx *= SMALLMULCONST; // Optionaly multiply by a small constant
                xx = (tx + bv) - bv2; // Round to the nearest integer
                if(error_log) 
		{                     // Compute the rounding error if required
                        if(xx<0)
                                err = fabs(-xx + tx);
                        else
                                err = fabs(xx  - tx);
                        if(err > maxerr)
                                maxerr = err;
                }
                xx += carry;
                zz = (xx+limitbv)-limitbv2;
                carry = zz*invlimit;            // Compute the carry on next word
                x[IDX(j)] = xx - zz;            // And the balanced remainder in current word
        }
        g_carry[threadID]=carry;
        if(g_err[threadID] < maxerr)
		g_err[threadID] = maxerr;
}

//      Functions that do the normalization of the result of a multiplication or squaring
__global__ void cuda_rnormalize2_kernel(        // Used for rational bases DWFFT
        double          *x,
        int             N,
        int             error_log,
        int             noadd,
        int             nomul,
        double          limitbv,
        double          limitbv2,
        double          invlimit,
        int             STRIDE,
        double          bv,
        double          bv2,
        double          *g_carry,
        float           *g_err,
        int             MULBYCONST,
        double          SMALLMULCONST,
        int             addinindex,
        double          addinvalue,
        int             plus
)
{
        register int    	j;
        register double         xx, zz, tx;
        register double         carry = 0.0;
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

        if(threadID*STRIDE < N && threadID*STRIDE+STRIDE >= N)
        {
                carry = g_carry[threadID];
                if (carry)
                {
                        j = 0;
                        if (plus) 
				carry = -carry;
                        while (carry)
                        {       
                                tx = x[IDX(j)];
                                xx = tx + carry;
                                zz = (xx+limitbv)-limitbv2;
                                carry = zz*invlimit; // Compute the carry on next word
                                x[IDX(j)] = xx - zz; // And the balanced remainder in current word
                                if (++j == N)
                                {
                                        j = 0;
                                        if (plus) 
						carry = -carry;
                                }
                        }
                }
        }
        else
        {
                carry = g_carry[threadID];
                j = threadID*STRIDE+STRIDE;
                while (carry)
                {
                        tx = x[IDX(j)];
                        xx = tx + carry;
                        zz = (xx+limitbv)-limitbv2;
                        carry = zz*invlimit;    // Compute the carry on next word
                        x[IDX(j)] = xx - zz;    // And the balanced remainder in current word
                        j++;
                }
        }
}

double
rnormalize(			// Used for rational bases DWFFT
	double 			*x,
	int 				N,
	int				error_log,
	int				noadd,
	int				nomul
)
{
	register int 	j;
	register double 	*px = x, xx, zz, bv = BIGVAL2, invlimit = limit_inverse_high; // JP 08/07/17
	register double  	carry = 0.0, limitbv = (limit_high*BIGVAL2)-BIGVAL2, err, maxerr = 0.0; // JP 08/07/17

	if (!noadd)
		x[addinindex] += addinvalue;	// Add the optional small constant
	for (j=0; j<N; ++j)
	{
		if (MULBYCONST && !nomul)
			*px *= SMALLMULCONST;	// Optionaly multiply by a small constant
		xx = (*px + bv) - bv;		// Round to the nearest integer
		if (error_log ) {		// Compute the rounding error if required
			if (xx<0)
          			err = fabs(-xx + *px);  
     			else 
     				err = fabs(xx  - *px);
	 		if (err > maxerr) 
	 			maxerr = err;
     		}

		xx += carry;
		zz = (xx+limitbv)-limitbv;
		carry = zz*invlimit;			// Compute the carry on next word
		*(px++) = xx - zz;			// And the balanced remainder in current word

	}

	if (carry)
	{
		j = 0;
		px = x;
		if (plus)
			carry = -carry;
		while (carry)
		{
			xx = *px + carry;
			zz = (xx+limitbv)-limitbv;
			carry = zz*invlimit;		// Compute the carry on next word
			*(px++) = xx - zz;		// And the balanced remainder in current word

			if (++j == N)
			{
				j = 0;
				px = x;
				if (plus)
					carry = -carry;
			}
		}
	}
	return(maxerr);
}
					        // Modular reduction of a zero padded integer
/*
__global__ void cuda_modred_kernel (
    double *s,
    double *scral,
    int len,
    int hwcount,
    int lwcount,
    int hindex,
    double bv,
    double limit_high,
    double mult,
    double invmult,
    double g_limitbv
    
)
{
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	register int 	i, j, liminf;
	register double carry = 0.0, xx, yy, zz, q;
        
        i = threadID;
        j = threadID + len -1;
        liminf = threadID + len -hindex;
        if (i<12)
            scral[i] = carry;   // Zero 12 double words in scratch area
            
}*/

void
modred (
	double *x
)
{
	register double *pscr = scral;		// Init. ; Point to scratch area
	register double *px = x+FFTLEN-1; 	// Point to last FFT word
	register long hindex = hwcount; 	// Count of upper FFT words
	register double carry = 0.0, bv = BIGVAL2, limitbv = (limit_high*bv)-bv, xx, yy, zz, q;

	if (debug)
		gwypstart_timer (3);

	*pscr++ = carry;					// Zero 12 double words in scratch area
	*pscr++ = carry;
	*pscr++ = carry;
	*pscr++ = carry;
	*pscr++ = carry;
	*pscr++ = carry;
	*pscr++ = carry;
	*pscr++ = carry;
	*pscr++ = carry;
	*pscr++ = carry;
	*pscr++ = carry;
	*pscr++ = carry;

	while (hindex-- > 0) {	// divide by mult the upper FFT words and save them in scratch area
		xx = *px+carry*limit_high;
		q = (xx*invmult+bv)-bv;		// q = xx/mult integer
		carry = xx-mult*q;			// carry = xx%h
		*px-- = 0.0;				// zero the high word
		*pscr++ = shift*q;			// save -inc*q*2^(ng-q)%bits
	}
	px++;
	pscr--;
	*px = carry;
	px = x;
	hindex = lwcount;
	carry = bv;
	while (hindex-- > 0) {	// Add or subtract saved words into the lower part of FFT
		xx = *px + *pscr-- + carry;			// xx = x + *pscr + carry
		yy = xx+limitbv;					// y = xx/limit_high*limit_high + bv
		zz = yy-limitbv;					// z = xx/limit_high*limit_high
		carry = limit_inverse_high*yy;		// carry = y / limit_high + bv
		*px++ = xx-zz;						// new x = xx%limit_high
	}
	px = x + hwoffset + 4;					// Point to 5th upper FFT word;
	carry = *px;
	*px-- = 0.0;

	xx = *px+carry*limit_high;
	q = (xx*invmult+bv)-bv;					// q = xx/mult integer
	carry = xx-mult*q;
	*px-- = 0.0;
	*pscr = shift*q;

	xx = *px+carry*limit_high;
	q = (xx*invmult+bv)-bv;					// q = xx/mult integer
	carry = xx-mult*q;
	*px-- = 0.0;
	*(pscr+1) = shift*q;

	xx = *px + carry*limit_high;
	q = (xx*invmult+bv)-bv;					// q = xx/mult integer
	carry = xx-mult*q;
	*px-- = 0.0;	// added
	*(pscr+2) = shift*q;

	xx = *px + carry*limit_high;
	q = (xx*invmult+bv)-bv;					// q = xx/mult integer
	carry = xx-mult*q;
	*px = 0.0;	// added
	*(pscr+3) = shift*q;
	
	while ((carry>(limit_high/2)) || (carry<(-limit_high/2))) {
		xx = carry+bv;						// Split the remainder
		yy = xx + limitbv;					// y = xx/limit_high*limit_high + bv
		zz = yy - limitbv;					// z = xx/limit_high*limit_high
		carry = limit_inverse_high*yy-bv;	// carry = y / limit_high
		*px++ = xx-zz;						// Save lower bits
	}
	*(px) = carry;
	px = x;								// Reload source pointer

	xx = *px + *(pscr+3) + bv;
	yy = xx + limitbv;					// y = xx/limit_high*limit_high + bv
	zz = yy - limitbv;					// z = xx/limit_high*limit_high
	carry = limit_inverse_high*yy;		// carry = y / limit_high
	*px++ = xx-zz;						// Save new value

	xx = *px + *(pscr+2) + carry;
	yy = xx + limitbv;					// y = xx/limit_high*limit_high + bv
	zz = yy - limitbv;					// z = xx/limit_high*limit_high
	carry = limit_inverse_high*yy;		// carry = y / limit_high
	*px++ = xx-zz;						// Save new value

	xx = *px + *(pscr+1) + carry;
	yy = xx + limitbv;					// y = xx/limit_high*limit_high + bv
	zz = yy - limitbv;					// z = xx/limit_high*limit_high
	carry = limit_inverse_high*yy;		// carry = y / limit_high
	*px++ = xx-zz;						// Save new value

	xx = *px + *(pscr+0) + carry;
	yy = xx + limitbv;					// y = xx/limit_high*limit_high + bv
	zz = yy - limitbv;					// z = xx/limit_high*limit_high
	carry = limit_inverse_high*yy;		// carry = y / limit_high
	*px++ = xx-zz;						// Save new value

	xx = *px + carry;
	yy = xx + limitbv;					// y = xx/limit_high*limit_high + bv
	zz = yy - limitbv;					// z = xx/limit_high*limit_high
	carry = limit_inverse_high*yy;		// carry = y / limit_high
	*px++ = xx-zz;						// Save new value

	*px	+= carry - bv;					// Adjust final carry

	if (debug)
		gwypend_timer (3);

}

void
check_balanced(						// Check if the balanced form of a result is correct
	double 	*x,
	int 	N
)
{
	int 	j;
	double 	lim, *ptrx = x;

	for (j=0; j<N; ++j)
	{
		lim = hlimit[j];
		assert ((*ptrx<=lim) && (*ptrx>=-lim));
		++ptrx;
	}
}

__global__ void
cuda_mulbyconst_kernel(
	double *out,
	double *in,
	double c,
	int n
)
{
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	if(threadID<n)
		out[threadID]=in[threadID]*c;
}

__global__ void
cuda_mul_kernel(
	double *out,
	double *in,
	double *in2,
	int n
)
{
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	if(threadID<n)
		out[threadID]=in[threadID]*in2[threadID];
}

__global__ void
cuda_mul_two_to_phi_kernel(
	double *x,
	double *cxin,
	double *two_to_phi,
	int hn
)
{
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	if(threadID<hn)
		cxin[threadID*2]=x[threadID]*two_to_phi[threadID];
		cxin[threadID*2+1]=x[threadID+hn]*two_to_phi[threadID+hn];
}

__global__ void
cuda_fold_kernel(
	double *x,
	double *cxin,
	int hn
)
{
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	if(threadID<hn)
		cxin[threadID*2]=x[threadID];
		cxin[threadID*2+1]=x[threadID+hn];
}

__global__ void
cuda_mul_two_to_minusphi_kernel(
	double *x,
	double *cxin,
	double *two_to_minusphi,
	int hn
)
{
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	if(threadID<hn)
		x[threadID]=cxin[threadID*2]*two_to_minusphi[threadID];
		x[threadID+hn]=cxin[threadID*2+1]*two_to_minusphi[threadID+hn];
}

__global__ void
cuda_unfold_kernel(
	double *x,
	double *cxin,
	double c,
	int hn
)
{
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	if(threadID<hn)
		x[threadID]=cxin[threadID*2]*c;
		x[threadID+hn]=cxin[threadID*2+1]*c;
}

//	Funtions that do squaring and multiplications of large integers ; inputs and output are normalized

double
cuda_lucas_square(	// Squaring of a large integer ; input and output normalized
	double			*x,
	int 			N,
	int 			error_log,
	int				noadd,
	int				nomul,
	int 			flag
)
{
	register int 		j, hn = N/2;
	register double 	err;
	int STRIDE;
	STRIDE=STRIDE_DIM;
	err=0;
       	dim3 grid(STRIDE_DIM/BLOCK_DIM,STRIDE_DIM/BLOCK_DIM, 1);
       	dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);

	if (zp || generic) {
            if(compl2)
            {
		if(flag&1){
			cutilSafeCall(cudaMemcpy(cuda_x,x,sizeof(double)*N,cudaMemcpyHostToDevice));
		}
		cuda_fold_kernel<<<(hn+127)/128,128>>>(cuda_x,cuda_cxin,hn);
		cuda_fftwsquare_g(N);	// DWT squaring
		cuda_unfold_kernel<<<(hn+127)/128,128>>>(cuda_x,cuda_cxin,(double)ttmp,hn);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_x[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
		cuda_rnormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_rnormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_x[j],(double *)&cuda_tmp[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
#ifdef linux
		cutilSafeCall(cudaMemcpy(x,cuda_x,sizeof(double),cudaMemcpyDeviceToHost));
#endif
        	if (error_log)
		{
			cutilSafeCall(cudaMemcpy(x,cuda_x,sizeof(double)*N,cudaMemcpyDeviceToHost));
                	check_balanced(x, N);
			cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
			for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
		}
		if(!error_log && flag&2){
			cutilSafeCall(cudaMemcpy(x,cuda_x,sizeof(double)*N,cudaMemcpyDeviceToHost));
		} 
		if (!generic)
                    modred(x);
        	return(err);
            }
            else
            {
		if(flag&1){
			cutilSafeCall(cudaMemcpy(cuda_xin,x,sizeof(double)*N,cudaMemcpyHostToDevice));
		}
		cuda_fftwsquare_g(N);	// DWT squaring

        	//dim3 grid(STRIDE_DIM/BLOCK_DIM,STRIDE_DIM/BLOCK_DIM, 1);
        	//dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		mul_const_transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_xin[j],(double)ttmp,(int)  STRIDE_DIM,(int) STRIDE_DIM);

		cuda_rnormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_rnormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);

		if (error_log){
        		for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           			transpose<<<grid, threads>>>((double *)&cuda_x[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
			cutilSafeCall(cudaMemcpy(x,cuda_x,sizeof(double)*N,cudaMemcpyDeviceToHost));
			cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
			for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
			check_balanced(x, N);
		}

		if(!error_log && flag&2) {
        		for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           			transpose<<<grid, threads>>>((double *)&cuda_x[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
			cutilSafeCall(cudaMemcpy(x,cuda_x,sizeof(double)*N,cudaMemcpyDeviceToHost));
		} 
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
        		transpose<<<grid, threads>>>((double *)&cuda_xin[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
#ifdef linux
		cutilSafeCall(cudaMemcpy(xin,cuda_xin,sizeof(double),cudaMemcpyDeviceToHost));
#endif
                if (!generic)
                    modred(x);
		return(err);
            }
        }
        else
            if(compl2)
            {
		if(flag&1){
			cutilSafeCall(cudaMemcpy(cuda_x,x,sizeof(double)*N,cudaMemcpyHostToDevice));
		}
		cuda_mul_two_to_phi_kernel<<<(hn+127)/128,128>>>(cuda_x,cuda_cxin,cuda_two_to_phi,hn);
		cuda_fftwsquare_g(N);	// DWT squaring
		cuda_mul_two_to_minusphi_kernel<<<(hn+127)/128,128>>>(cuda_x,cuda_cxin,cuda_two_to_minusphi,hn);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_x[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
		cuda_inormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,g_limitbv,g_limitbv,g_invlimit,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,wrapindex,wrapfactor,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_inormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,g_limitbv,g_limitbv,g_invlimit,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,wrapindex,wrapfactor,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_x[j],(double *)&cuda_tmp[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
#ifdef linux
		cutilSafeCall(cudaMemcpy(x,cuda_x,sizeof(double),cudaMemcpyDeviceToHost));
#endif
        	if (error_log)
		{
			cutilSafeCall(cudaMemcpy(x,cuda_x,sizeof(double)*N,cudaMemcpyDeviceToHost));
                	check_balanced(x, N);
			cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
			for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
		}
		if(!error_log && flag&2){
			cutilSafeCall(cudaMemcpy(x,cuda_x,sizeof(double)*N,cudaMemcpyDeviceToHost));
		} 
        	return(err);
            }
            else
            {
		if(flag&1){
			cutilSafeCall(cudaMemcpy(cuda_x,x,sizeof(double)*N,cudaMemcpyHostToDevice));
			cuda_mul_kernel<<<(N+127)/128,128>>>(cuda_xin,cuda_x,cuda_two_to_phi, N);
		}
		cuda_fftwsquare_g(N);	// DWT squaring

        	//dim3 grid(STRIDE_DIM/BLOCK_DIM,STRIDE_DIM/BLOCK_DIM, 1);
        	//dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		mul_0_transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_xin[j],(double *)&cuda_two_to_minusphi[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);

		cuda_inormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,g_limitbv,g_limitbv,g_invlimit,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,wrapindex,wrapfactor,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_inormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,g_limitbv,g_limitbv,g_invlimit,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,wrapindex,wrapfactor,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);

		if (error_log){
        		for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           			transpose<<<grid, threads>>>((double *)&cuda_x[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
			cutilSafeCall(cudaMemcpy(x,cuda_x,sizeof(double)*N,cudaMemcpyDeviceToHost));
			cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
			for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
			check_balanced(x, N);
		}

		if(!error_log && flag&2) {
        		for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           			transpose<<<grid, threads>>>((double *)&cuda_x[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
			cutilSafeCall(cudaMemcpy(x,cuda_x,sizeof(double)*N,cudaMemcpyDeviceToHost));
		} 
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
        		mul_1_transpose<<<grid, threads>>>((double *)&cuda_xin[j],(double *)&cuda_tmp[j],(double *)&cuda_two_to_phi[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
#ifdef linux
		cutilSafeCall(cudaMemcpy(xin,cuda_xin,sizeof(double),cudaMemcpyDeviceToHost));
#endif

		return(err);
            }
}

double
cuda_lucas_square_generic(	// Squaring of a large integer ; input and output normalized, generic reduction included.
	double			*x,
	int 			N,
	int 			error_log,
	int				noadd,
	int				nomul,
	int 			flag
)
{
	register int 		j, hn = N/2;
	register double 	err;
	int STRIDE;
	STRIDE=STRIDE_DIM;
	err=0;
       	dim3 grid(STRIDE_DIM/BLOCK_DIM,STRIDE_DIM/BLOCK_DIM, 1);
       	dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);

            if(compl2)
            {
                if(flag&1){
                    cutilSafeCall(cudaMemcpy(cuda_x,x,sizeof(double)*N,cudaMemcpyHostToDevice));
		}
		cuda_fold_kernel<<<(hn+127)/128,128>>>(cuda_x,cuda_cxin,hn);
		cuda_fftwsquare_g(N);	// DWT squaring
		cuda_unfold_kernel<<<(hn+127)/128,128>>>(cuda_x,cuda_cxin,(double)ttmp,hn);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_x[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
		cuda_rnormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_rnormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_x[j],(double *)&cuda_tmp[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
        	if (error_log)
		{
			cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
			for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
		}
		cuda_gwypcopyzero_kernel<<<(N+127)/128,128>>>(cuda_x, cuda_tmp_g, zerowordslow, FFTLEN);
		cuda_fold_kernel<<<(hn+127)/128,128>>>(cuda_tmp_g,cuda_cyin,hn);
                cuda_fftwmulbyr_g(N);
		cuda_unfold_kernel<<<(hn+127)/128,128>>>(cuda_tmp_g,cuda_cyin,(double)ttmp,hn);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_tmp_g[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
		cuda_rnormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, 1, 1,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_rnormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, 1, 1,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_tmp_g[j],(double *)&cuda_tmp[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
        	if (error_log)
		{
			cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
			for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
		}
		cuda_gwypsetzero_kernel<<<(N+127)/128,128>>>(cuda_tmp_g, zerowordshigh, FFTLEN);
		cuda_fold_kernel<<<(hn+127)/128,128>>>(cuda_tmp_g,cuda_cyin,hn);
                cuda_fftwmulbym_g(N);
		cuda_unfold_kernel<<<(hn+127)/128,128>>>(cuda_tmp_g,cuda_cyin,(double)ttmp,hn);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_tmp_g[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
		cuda_rnormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, 1, 1,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_rnormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, 1, 1,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_tmp_g[j],(double *)&cuda_tmp[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
        	if (error_log)
		{
			cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
			for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
		}
		cuda_gwypaddquick_kernel<<<(N+127)/128,128>>>(cuda_tmp_g, cuda_x, FFTLEN);
                cuda_gwypsubquick_kernel<<<(N+127)/128,128>>>(cuda_m, cuda_x, FFTLEN);
                cutilSafeCall(cudaMemcpy(x,cuda_x,sizeof(double)*N,cudaMemcpyDeviceToHost));
        	return(err);
            }
            else
            {
		if(flag&1){
			cutilSafeCall(cudaMemcpy(cuda_xin,x,sizeof(double)*N,cudaMemcpyHostToDevice));
		}
		cuda_fftwsquare_g(N);	// DWT squaring

        	//dim3 grid(STRIDE_DIM/BLOCK_DIM,STRIDE_DIM/BLOCK_DIM, 1);
        	//dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		mul_const_transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_xin[j],(double)ttmp,(int)  STRIDE_DIM,(int) STRIDE_DIM);

		cuda_rnormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_rnormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
                for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
                    transpose<<<grid, threads>>>((double *)&cuda_x[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
		if (error_log){
                    cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
                    for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
                }
		cuda_gwypcopyzero_kernel<<<(N+127)/128,128>>>(cuda_x, cuda_yin, zerowordslow, FFTLEN);
                cuda_fftwmulbyr_g(N);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		mul_const_transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_yin[j],(double)ttmp,(int)  STRIDE_DIM,(int) STRIDE_DIM);

		cuda_rnormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, 1, 1,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_rnormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, 1, 1,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
                for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
                    transpose<<<grid, threads>>>((double *)&cuda_yin[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
		if (error_log){
                    cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
                    for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
                }
		cuda_gwypsetzero_kernel<<<(N+127)/128,128>>>(cuda_yin, zerowordshigh, FFTLEN);
                cuda_fftwmulbym_g(N);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		mul_const_transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_yin[j],(double)ttmp,(int)  STRIDE_DIM,(int) STRIDE_DIM);

		cuda_rnormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, 1, 1,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_rnormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, 1, 1,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
                for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
                    transpose<<<grid, threads>>>((double *)&cuda_tmp_g[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
		if (error_log){
                    cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
                    for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
                }
                cuda_gwypsubquick_kernel<<<(N+127)/128,128>>>(cuda_tmp_g, cuda_x, FFTLEN);
                cutilSafeCall(cudaMemcpy(x,cuda_x,sizeof(double)*N,cudaMemcpyDeviceToHost));
		return(err);
            }
}

//	Funtions that do squaring and multiplications of large integers ; inputs and output are normalized

double
lucas_square(	// Squaring of a large integer ; input and output normalized
	double			*x,
	int 			N,
	int 			error_log,
	int				noadd,
	int				nomul
)
{
	register int 		j, hn = N/2;
	register double 	err;

	if (debug)
		gwypstart_timer (5);

	if (zp || generic)
		if (compl2) {	// Transform an N sized real FFT array into an N/2 sized complex FFT one
			for (j=0; j<hn; ++j)
			{
				cxin[j][0] = x[j];
				cxin[j][1] = x[j+hn];
			}
		}
		else {
			for (j=0; j<N; ++j)
			{
				xin[j] = x[j];
			}
		}
	else
		if (compl2) {	// Transform an N sized real FFT array into an N/2 sized complex FFT one
			for (j=0; j<hn; ++j)
			{
				cxin[j][0] = x[j] * two_to_phi[j];
				cxin[j][1] = x[j+hn] * two_to_phi[j+hn];
			}
		}
		else {
			for (j=0; j<N; ++j)
			{
				xin[j] = x[j] * two_to_phi[j];
			}
		}

	fftwsquare_g(N);	// DWT squaring

	if (zp || generic)
		if (compl2) {	// Unfold the N/2 sized complex ouput array into a N sized real one
			for (j=0; j<hn; ++j)
			{
				x[j] = cxin[j][0] * ttmp;
				x[j+hn] = cxin[j][1] * ttmp;
			}
		}
		else {
			for (j=0; j<N; ++j)
			{
				x[j] = xin[j] * ttmp;
			}
		}
	else
		if (compl2) {	// Unfold the N/2 sized complex ouput array into a N sized real one
			for (j=0; j<hn; ++j)
			{
				x[j] = cxin[j][0] *  two_to_minusphi[j];
				x[j+hn] = cxin[j][1] *  two_to_minusphi[j+hn];
			}
		}
		else {
			for (j=0; j<N; ++j)
			{
				x[j] = xin[j] *  two_to_minusphi[j];
			}
		}

	if (debug)
		gwypend_timer (5);

	if (debug)
		gwypstart_timer (2);
	err = (zp || generic)? rnormalize(x, N, error_log, noadd, nomul) : inormalize(x, N, error_log, noadd, nomul);
//	if (error_log)
//		check_balanced(x, N);
	if (debug)
		gwypend_timer (2);

	if (zp)
		modred (x);
	if (error_log)
		check_balanced(x, N);

	return(err);
}

double
cuda_lucas_mul(	// Multiplication of large integers ; inputs and output normalized
	double			*x,
        double                  *y,
	int 			N,
	int 			error_log,
	int				noadd,
	int				nomul,
	int 			flag
)
{
	register int 		j, hn = N/2;
	register double 	err;
	int STRIDE;
	STRIDE=STRIDE_DIM;
	err=0;
       	dim3 grid(STRIDE_DIM/BLOCK_DIM,STRIDE_DIM/BLOCK_DIM, 1);
       	dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);

	if (zp || generic)
            if(compl2)
            {
		if(flag&1){
			cutilSafeCall(cudaMemcpy(cuda_x,x,sizeof(double)*N,cudaMemcpyHostToDevice));
		}
                cutilSafeCall(cudaMemcpy(cuda_y,y,sizeof(double)*N,cudaMemcpyHostToDevice));
		cuda_fold_kernel<<<(hn+127)/128,128>>>(cuda_x,cuda_cxin,hn);
		cuda_fold_kernel<<<(hn+127)/128,128>>>(cuda_y,cuda_cyin,hn);
		cuda_fftwmul_g(N);	// DWT multipication
		cuda_unfold_kernel<<<(hn+127)/128,128>>>(cuda_y,cuda_cyin,(double)ttmp,hn);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_y[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
		cuda_rnormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_rnormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_y[j],(double *)&cuda_tmp[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
#ifdef linux
		cutilSafeCall(cudaMemcpy(y,cuda_y,sizeof(double),cudaMemcpyDeviceToHost));
#endif
        	if (error_log)
		{
			cutilSafeCall(cudaMemcpy(y,cuda_y,sizeof(double)*N,cudaMemcpyDeviceToHost));
                	check_balanced(y, N);
			cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
			for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
		}
		if(!error_log && flag&2){
			cutilSafeCall(cudaMemcpy(y,cuda_y,sizeof(double)*N,cudaMemcpyDeviceToHost));
		} 
		if (!generic)
                    modred (y);
        	return(err);
            }
            else
            {
		if(flag&1){
			cutilSafeCall(cudaMemcpy(cuda_xin,x,sizeof(double)*N,cudaMemcpyHostToDevice));
		}
                cutilSafeCall(cudaMemcpy(cuda_yin,y,sizeof(double)*N,cudaMemcpyHostToDevice));
		cuda_fftwmul_g(N);	// DWT multiplication

        	//dim3 grid(STRIDE_DIM/BLOCK_DIM,STRIDE_DIM/BLOCK_DIM, 1);
        	//dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		mul_const_transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_yin[j],(double)ttmp,(int)  STRIDE_DIM,(int) STRIDE_DIM);

		cuda_rnormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_rnormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);

		if (error_log){
        		for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           			transpose<<<grid, threads>>>((double *)&cuda_y[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
			cutilSafeCall(cudaMemcpy(y,cuda_y,sizeof(double)*N,cudaMemcpyDeviceToHost));
			cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
			for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
			check_balanced(y, N);
		}

		if(!error_log && flag&2) {
        		for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           			transpose<<<grid, threads>>>((double *)&cuda_y[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
			cutilSafeCall(cudaMemcpy(y,cuda_y,sizeof(double)*N,cudaMemcpyDeviceToHost));
		} 
//        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
//        		mul_1_transpose<<<grid, threads>>>((double *)&cuda_yin[j],(double *)&cuda_tmp[j],(double *)&cuda_two_to_phi[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
#ifdef linux
		cutilSafeCall(cudaMemcpy(yin,cuda_yin,sizeof(double),cudaMemcpyDeviceToHost));
#endif
                if (!generic)
                    modred (y);
		return(err);
            }
        else
            if(compl2)
            {
		if(flag&1){
			cutilSafeCall(cudaMemcpy(cuda_x,x,sizeof(double)*N,cudaMemcpyHostToDevice));
		}
                cutilSafeCall(cudaMemcpy(cuda_y,y,sizeof(double)*N,cudaMemcpyHostToDevice));
		cuda_mul_two_to_phi_kernel<<<(hn+127)/128,128>>>(cuda_x,cuda_cxin,cuda_two_to_phi,hn);
		cuda_mul_two_to_phi_kernel<<<(hn+127)/128,128>>>(cuda_y,cuda_cyin,cuda_two_to_phi,hn);
		cuda_fftwmul_g(N);	// DWT multipication
		cuda_mul_two_to_minusphi_kernel<<<(hn+127)/128,128>>>(cuda_y,cuda_cyin,cuda_two_to_minusphi,hn);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_y[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
		cuda_inormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,g_limitbv,g_limitbv,g_invlimit,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,wrapindex,wrapfactor,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_inormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,g_limitbv,g_limitbv,g_invlimit,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,wrapindex,wrapfactor,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_y[j],(double *)&cuda_tmp[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
#ifdef linux
		cutilSafeCall(cudaMemcpy(y,cuda_y,sizeof(double),cudaMemcpyDeviceToHost));
#endif
        	if (error_log)
		{
			cutilSafeCall(cudaMemcpy(y,cuda_y,sizeof(double)*N,cudaMemcpyDeviceToHost));
                	check_balanced(y, N);
			cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
			for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
		}
		if(!error_log && flag&2){
			cutilSafeCall(cudaMemcpy(y,cuda_y,sizeof(double)*N,cudaMemcpyDeviceToHost));
		} 
        	return(err);
            }
            else
            {
		if(flag&1){
			cutilSafeCall(cudaMemcpy(cuda_x,x,sizeof(double)*N,cudaMemcpyHostToDevice));
		}
                cutilSafeCall(cudaMemcpy(cuda_y,y,sizeof(double)*N,cudaMemcpyHostToDevice));
                cuda_mul_kernel<<<(N+127)/128,128>>>(cuda_xin,cuda_x,cuda_two_to_phi, N);
                cuda_mul_kernel<<<(N+127)/128,128>>>(cuda_yin,cuda_y,cuda_two_to_phi, N);
		cuda_fftwmul_g(N);	// DWT multiplication

        	//dim3 grid(STRIDE_DIM/BLOCK_DIM,STRIDE_DIM/BLOCK_DIM, 1);
        	//dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		mul_0_transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_yin[j],(double *)&cuda_two_to_minusphi[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);

		cuda_inormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,g_limitbv,g_limitbv,g_invlimit,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,wrapindex,wrapfactor,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_inormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,g_limitbv,g_limitbv,g_invlimit,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,wrapindex,wrapfactor,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);

		if (error_log){
        		for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           			transpose<<<grid, threads>>>((double *)&cuda_y[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
			cutilSafeCall(cudaMemcpy(y,cuda_y,sizeof(double)*N,cudaMemcpyDeviceToHost));
			cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
			for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
			check_balanced(y, N);
		}

		if(!error_log && flag&2) {
        		for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           			transpose<<<grid, threads>>>((double *)&cuda_y[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
			cutilSafeCall(cudaMemcpy(y,cuda_y,sizeof(double)*N,cudaMemcpyDeviceToHost));
		} 
//        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
//        		mul_1_transpose<<<grid, threads>>>((double *)&cuda_yin[j],(double *)&cuda_tmp[j],(double *)&cuda_two_to_phi[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
#ifdef linux
		cutilSafeCall(cudaMemcpy(yin,cuda_yin,sizeof(double),cudaMemcpyDeviceToHost));
#endif

		return(err);
	}
}

double
cuda_lucas_mul_generic(	// Multiplication of large integers ; inputs and output normalized
	double			*x,
        double                  *y,
	int 			N,
	int 			error_log,
	int				noadd,
	int				nomul,
	int 			flag
)
{
	register int 		j, hn = N/2;
	register double 	err;
	int STRIDE;
	STRIDE=STRIDE_DIM;
	err=0;
       	dim3 grid(STRIDE_DIM/BLOCK_DIM,STRIDE_DIM/BLOCK_DIM, 1);
       	dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);

            if(compl2)
            {
		if(flag&1){
			cutilSafeCall(cudaMemcpy(cuda_x,x,sizeof(double)*N,cudaMemcpyHostToDevice));
		}
                cutilSafeCall(cudaMemcpy(cuda_y,y,sizeof(double)*N,cudaMemcpyHostToDevice));
		cuda_fold_kernel<<<(hn+127)/128,128>>>(cuda_x,cuda_cxin,hn);
		cuda_fold_kernel<<<(hn+127)/128,128>>>(cuda_y,cuda_cyin,hn);
		cuda_fftwmul_g(N);	// DWT multipication
		cuda_unfold_kernel<<<(hn+127)/128,128>>>(cuda_y,cuda_cyin,(double)ttmp,hn);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_y[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
		cuda_rnormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_rnormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_y[j],(double *)&cuda_tmp[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
		cuda_gwypcopyzero_kernel<<<(N+127)/128,128>>>(cuda_y, cuda_tmp_g, zerowordslow, FFTLEN);
		cuda_fold_kernel<<<(hn+127)/128,128>>>(cuda_tmp_g,cuda_cyin,hn);
                cuda_fftwmulbyr_g(N);
		cuda_unfold_kernel<<<(hn+127)/128,128>>>(cuda_tmp_g,cuda_cyin,(double)ttmp,hn);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_tmp_g[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
		cuda_rnormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, 1, 1,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_rnormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, 1, 1,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_tmp_g[j],(double *)&cuda_tmp[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
        	if (error_log)
		{
			cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
			for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
		}
		cuda_gwypsetzero_kernel<<<(N+127)/128,128>>>(cuda_tmp_g, zerowordshigh, FFTLEN);
		cuda_fold_kernel<<<(hn+127)/128,128>>>(cuda_tmp_g,cuda_cyin,hn);
                cuda_fftwmulbym_g(N);
		cuda_unfold_kernel<<<(hn+127)/128,128>>>(cuda_tmp_g,cuda_cyin,(double)ttmp,hn);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_tmp_g[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
		cuda_rnormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, 1, 1,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_rnormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, 1, 1,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		transpose<<<grid, threads>>>((double *)&cuda_tmp_g[j],(double *)&cuda_tmp[j],(int)  STRIDE_DIM,(int) STRIDE_DIM);
        	if (error_log)
		{
			cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
			for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
		}
		cuda_gwypaddquick_kernel<<<(N+127)/128,128>>>(cuda_tmp_g, cuda_y, FFTLEN);
                cuda_gwypsubquick_kernel<<<(N+127)/128,128>>>(cuda_m, cuda_y, FFTLEN);
                cutilSafeCall(cudaMemcpy(y,cuda_y,sizeof(double)*N,cudaMemcpyDeviceToHost));
        	return(err);
            }
            else
            {
		if(flag&1){
			cutilSafeCall(cudaMemcpy(cuda_xin,x,sizeof(double)*N,cudaMemcpyHostToDevice));
		}
                cutilSafeCall(cudaMemcpy(cuda_yin,y,sizeof(double)*N,cudaMemcpyHostToDevice));
		cuda_fftwmul_g(N);	// DWT multiplication

        	//dim3 grid(STRIDE_DIM/BLOCK_DIM,STRIDE_DIM/BLOCK_DIM, 1);
        	//dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		mul_const_transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_yin[j],(double)ttmp,(int)  STRIDE_DIM,(int) STRIDE_DIM);

		cuda_rnormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_rnormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, noadd, nomul,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);

                for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
                    transpose<<<grid, threads>>>((double *)&cuda_y[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
		if (error_log){
                    cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
                    for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
		}
		cuda_gwypcopyzero_kernel<<<(N+127)/128,128>>>(cuda_y, cuda_yin, zerowordslow, FFTLEN);
                cuda_fftwmulbyr_g(N);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		mul_const_transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_yin[j],(double)ttmp,(int)  STRIDE_DIM,(int) STRIDE_DIM);

		cuda_rnormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, 1, 1,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_rnormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, 1, 1,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
                for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
                    transpose<<<grid, threads>>>((double *)&cuda_yin[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
		if (error_log){
                    cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
                    for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
                }
		cuda_gwypsetzero_kernel<<<(N+127)/128,128>>>(cuda_yin, zerowordshigh, FFTLEN);
                cuda_fftwmulbym_g(N);
        	for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
           		mul_const_transpose<<<grid, threads>>>((double *)&cuda_tmp[j],(double *)&cuda_yin[j],(double)ttmp,(int)  STRIDE_DIM,(int) STRIDE_DIM);

		cuda_rnormalize_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, 1, 1,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
		cuda_rnormalize2_kernel<<<(N+STRIDE*128-1)/STRIDE/128,128>>>(cuda_tmp, N, error_log, 1, 1,limit_high_bigval,limit_high_bigval,limit_inverse_high,STRIDE,BIGVAL2,BIGVAL2,g_carry,g_err,MULBYCONST,SMALLMULCONST,addinindex,addinvalue,plus);
                for(j=0;j<N;j+=(STRIDE_DIM*STRIDE_DIM))
                    transpose<<<grid, threads>>>((double *)&cuda_tmp_g[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
		if (error_log){
                    cutilSafeCall(cudaMemcpy(l_err,g_err,sizeof(float)*(N+STRIDE-1)/STRIDE,cudaMemcpyDeviceToHost));
                    for(j=0;j<(N+STRIDE-1)/STRIDE;j++)if(err<l_err[j])err=l_err[j];
                }
                cuda_gwypsubquick_kernel<<<(N+127)/128,128>>>(cuda_tmp_g, cuda_y, FFTLEN);
                cutilSafeCall(cudaMemcpy(y,cuda_y,sizeof(double)*N,cudaMemcpyDeviceToHost));
		return(err);
            }
}

double
lucas_mul(	// Multiplication of large integers ; inputs and output normalized
	double	*x,
	double	*y,
	int 	N,
	int 	error_log,
	int		noadd,
	int		nomul
)
{
	register int 		j, hn = N/2;
	register double 	err;

	if (debug)
		gwypstart_timer (5);

	if (zp || generic)
		if (compl2) {	// Transform N sized real FFT arrays into N/2 sized complex FFT ones
			for (j=0; j<hn; ++j)
			{
				cxin[j][0] = x[j];
				cxin[j][1] = x[j+hn];
				cyin[j][0] = y[j];
				cyin[j][1] = y[j+hn];
			}
		}
		else {
			for (j=0; j<N; ++j)
			{
				xin[j] = x[j];
				yin[j] = y[j];
			}
		}
	else
		if (compl2) {	// Transform N sized real FFT arrays into N/2 sized complex FFT ones
			for (j=0; j<hn; ++j)
			{
				cxin[j][0] = x[j] * two_to_phi[j];
				cxin[j][1] = x[j+hn] * two_to_phi[j+hn];
				cyin[j][0] = y[j] * two_to_phi[j];
				cyin[j][1] = y[j+hn] * two_to_phi[j+hn];
			}
		}
		else {
			for (j=0; j<N; ++j)
			{
				xin[j] = x[j] * two_to_phi[j];
				yin[j] = y[j] * two_to_phi[j];
			}
		}

	fftwmul_g(N);	// DWT mutiplication

	if (zp || generic)
		if (compl2) {	// Unfold the N/2 sized complex ouput array into a N sized real one
			for (j=0; j<hn; ++j)
			{
				y[j] = cyin[j][0] * ttmp;
				y[j+hn] = cyin[j][1] * ttmp;
			}
		}
		else {
			for (j=0; j<N; ++j)
			{
				y[j] = yin[j] * ttmp;
			}
		}
	else
		if (compl2) {	// Unfold the N/2 sized complex ouput array into a N sized real one
			for (j=0; j<hn; ++j)
			{
				y[j] = cyin[j][0] *  two_to_minusphi[j];
				y[j+hn] = cyin[j][1] *  two_to_minusphi[j+hn];
			}
		}
		else {
			for (j=0; j<N; ++j)
			{
				y[j] = yin[j] *  two_to_minusphi[j];
			}
		}

	if (debug)
		gwypend_timer (5);

	if (debug)
		gwypstart_timer (2);
	err = (zp || generic)? rnormalize(y, N, error_log, noadd, nomul) : inormalize(y, N, error_log, noadd, nomul);
	if (debug)
		gwypend_timer (2);

	if (zp)
		modred (y);
	if (error_log)
		check_balanced(y, N);

	return(err);
}

/* ------------ Initialization routines ------------------- */


int set_fftlen (double k, unsigned long b, unsigned long n) {
    double bpw, bpwg, tbpw, tbpwg, kbits = 0, rdwt, rzpad, rgeneric; 
    double log2k, log2b, bitsmc = log(MAXMULCONST)/log(2.0);
    int incfft, fftwincr, fftwincrg, fftdwt = 2, fftzpad = 2, fftgeneric = 2 , zpad = 0;
    unsigned long len = bitlen (gmodulus), ki = (unsigned long)k;

    printfunction = (verbose)? both_output : screen_output;
    log2b = log(b)/log(2.0);
        
    while (ki) {    // Compute bit lengh of k
        kbits += 1;
        ki >>= 1;
    }
    
    if (k > 0.0) {
        log2k = log(k)/log(2.0);
        rdwt = n*log2b + log2k;	// Exponent for DWT mode
    }
    else {
        rdwt = n*log2b;	// Exponent for DWT mode
    }

    rzpad = len + len + 32; // Exponent for zero padded or generic modes
    rgeneric = rzpad + 2*EB;

    while (1) {	// Compute zero padded FFT length (raw)
        bpw = RINT((rzpad+fftzpad-1)/fftzpad);
        tbpw = 2.0*bpw + bitsmc;
        if (tbpw <= MAXBITSPERDOUBLE)
            break;
        else
            fftzpad *= 2;
    }

    while (1) {	// Compute generic FFT length (raw)
        bpwg = RINT((rgeneric+fftgeneric-1)/fftgeneric);
        tbpwg = 2.0*bpwg + bitsmc + 2.0;
        if (tbpwg <= MAXBITSPERDOUBLE)
            break;
        else
            fftgeneric *= 2;
    }

    incfft = FFTINC; // Compute zero padded FFT length (fine)
    while ((fftzpad <= 64) && (incfft-- > 0))
        fftzpad *= 2;

    if (fftzpad > 64) {
        fftzpad /= 2;
        fftwincr = fftzpad / 16;
        while (1) {
            bpw = RINT((rzpad+fftzpad-1)/fftzpad);
            tbpw = 2.0*bpw + bitsmc;
            if (tbpw <= MAXBITSPERDOUBLE)
                break;
            else
                fftzpad += fftwincr;
        }
        while (incfft-- > 0)
            fftzpad += fftwincr;
    }

    incfft = FFTINC;	// Compute generic FFT length (fine)
    while ((fftgeneric <= 64) && (incfft-- > 0))
        fftgeneric *= 2;

    if (fftgeneric > 64) {
        fftgeneric /= 2;
        fftwincrg = fftgeneric / 16;
        while (1) {
            bpwg = RINT((rgeneric+fftgeneric-1)/fftgeneric);
            tbpwg = 2.0*bpwg + bitsmc + 2.0;
            if (tbpwg <= MAXBITSPERDOUBLE)
                break;
            else
                fftgeneric += fftwincrg;
        }
        while (incfft-- > 0)
            fftgeneric += fftwincrg;
    }

//    if (k > 0.0) {
        while (1) {     // Compute DWT FFT length (raw)
            bpw = RINT((rdwt+fftdwt-1)/fftdwt);
            tbpw = 2.0*bpw + kbits + bitsmc + 0.6*log((double)fftdwt)/log(2.0);
            if (tbpw <= MAXBITSPERDOUBLE || bpw < 4.0)
                break;
            else
                fftdwt *= 2;
        }
        incfft = FFTINC; // Compute DWT FFT length (fine)
        while ((fftdwt <= 64) && (incfft-- > 0))
            fftdwt *= 2;
        if (fftdwt > 64) {
            fftdwt /= 2;
            fftwincr = fftdwt / 16;
            while (1) {
                bpw = RINT((rdwt+fftdwt-1)/fftdwt);
                tbpw = 2.0*bpw + kbits + bitsmc + 0.6*log((double)fftdwt)/log(2.0);
                if (tbpw <= MAXBITSPERDOUBLE || bpw < 4.0)
                    break;
                else
                    fftdwt += fftwincr;
            }
            while (incfft-- > 0)
                fftdwt += fftwincr;
        }
//    }
    
    if(fftgeneric < STRIDE_DIM*2) {
        cufftonly = 1;
//      generic = 1;
        if (!generic)
            zpad = 1;
    }


    if (generic) {
        zpad = 0;
        FFTLEN = fftgeneric;
        bpwg = RINT((rgeneric+fftgeneric-1)/fftgeneric);
        ng = (unsigned long)bpwg*FFTLEN;
        kg = 1.0;
        grecip = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);
        itog (1, grecip); 
        gshiftleft (len + len + EB, grecip); 
        divg (gmodulus, grecip);// compute len+EB+1 bits of reciprocal 
        gshiftleft (ng - len - len - EB, grecip);
        // shift so gwmul routines wrap quotient to lower end of fft 
        zerowordslow = (len - EB) / (unsigned long)bpwg; 
        zerowordshigh = FFTLEN - (len + (unsigned long)bpwg - 1) / (unsigned long)bpwg - 1; 
    }
    else if (zpad || (/*b != 2 || */kbits > MAXKBITS /*|| kbits  > MAXBITSPERDOUBLE || fftdwt > fftzpad ||  bpw < 5.0*/)) {
        zpad = 1;
        FFTLEN = fftzpad;
        bpw = RINT((rzpad+fftzpad-1)/fftzpad);
        ng = (unsigned long)bpw*FFTLEN;
        kg = 1.0;
    }
    else {
        zpad = 0;
        bpw = RINT((rdwt+fftdwt-1)/fftdwt);
        FFTLEN = fftdwt;
//	}       Use the code below only for IBDWT, JP 08/07/17
	
// 	{       And not for rational base !
        int ii,jj;
        if(g_fftlen != 0)
            FFTLEN = g_fftlen * 2 + 2 ;
        for(jj = ii = 64;ii < (FFTLEN >> 1) ;ii*=2)
            for(jj=ii;jj < (FFTLEN >> 1) && jj < ii*2 ;jj+=ii/4);
        FFTLEN = jj;
        g_fftlen = jj;
        {
            if(g_fftlen != 0)
                FFTLEN = g_fftlen * 2 + 2 ;
            for(jj = ii = 64;ii < (FFTLEN >> 1) ;ii*=2)
                for(jj=ii;jj < (FFTLEN >> 1) && jj < ii*2 ;jj+=ii/4);
            FFTLEN = jj;
            g_fftlen = jj;
        }

        if(s_FFTLEN > 0)
            {
                FFTLEN = s_FFTLEN;
                g_fftlen = s_FFTLEN;
                s_FFTLEN = 0;
            }
    }      // End special IBDWT code.
	
    tbpw = (zpad || generic) ? 2.0*bpw + bitsmc : 2.0*bpw + kbits + bitsmc + 0.6*log((double)FFTLEN)/log(2.0);
        
    if (debug) {
        sprintf (gwypbuf, "FFTLEN = %d, bpw = %f, Bits per double = %f, Maxbpd = %f\n",
        FFTLEN, bpw, tbpw, MAXBITSPERDOUBLE);
        if (printfunction != NULL)
            (*printfunction)(gwypbuf);
    }
    return (zpad);
}

void gwypset_larger_fftlen_count(
	int count
)
{
	FFTINC = count;
}

void	init_fftw (int n) {
	int j, hn = n/2;


        if (compl2) {
                if (cnp==NULL)
                    cnp = (double *)malloc(FFTLEN*sizeof(double));
                if (snp==NULL)
                    snp = (double *)malloc(FFTLEN*sizeof(double));
                for (j=0;j<FFTLEN;j++)
                {
                        cnp[j] = fftcosp(j);
                        snp[j] = fftsinp(j);
                }
                if (!cufftonly) {
                    if (cuda_cnp == NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_cnp, sizeof(double)*(n+STRIDE_DIM*STRIDE_DIM)));
                    if (cuda_snp == NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_snp, sizeof(double)*(n+STRIDE_DIM*STRIDE_DIM)));
                    cutilSafeCall(cudaMemcpy(cuda_cnp,cnp,n*sizeof(double),cudaMemcpyHostToDevice));
                    cutilSafeCall(cudaMemcpy(cuda_snp,snp,n*sizeof(double),cudaMemcpyHostToDevice));
                }
                if (cxin==NULL)
                    cxin = (fftw_complex *) malloc (hn * sizeof(fftw_complex));
                if (cyin==NULL)
                    cyin = (fftw_complex *) malloc (hn * sizeof(fftw_complex));
                if (cxout==NULL)
                    cxout = (fftw_complex *) malloc (hn * sizeof(fftw_complex));
                if (cyout==NULL)
                    cyout = (fftw_complex *) malloc (hn * sizeof(fftw_complex));
                if (cuda_cxin == NULL)
                    cutilSafeCall(cudaMalloc((void**)&cuda_cxin, sizeof(fftw_complex)*(hn+STRIDE_DIM*STRIDE_DIM)));
                if (cuda_cxout == NULL)
                    cutilSafeCall(cudaMalloc((void**)&cuda_cxout, sizeof(fftw_complex)*(hn+STRIDE_DIM*STRIDE_DIM)));
                if (cuda_cyin == NULL)
                    cutilSafeCall(cudaMalloc((void**)&cuda_cyin, sizeof(fftw_complex)*(hn+STRIDE_DIM*STRIDE_DIM)));
                if (cuda_cyout == NULL)
                    cutilSafeCall(cudaMalloc((void**)&cuda_cyout, sizeof(fftw_complex)*(hn+STRIDE_DIM*STRIDE_DIM)));
                cufftSafeCall(cufftPlan1d(&cuda_fwpx,FFTLEN/2, CUFFT_Z2Z, 1));
                if (generic && !cufftonly) {
                    if (cuda_tmp_g==NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_tmp_g, sizeof(fftw_complex)*(hn+STRIDE_DIM*STRIDE_DIM)));
                    if (cuda_m==NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_m, sizeof(fftw_complex)*(hn+STRIDE_DIM*STRIDE_DIM)));
                   cutilSafeCall(cudaMemcpy(cuda_m,modulus,sizeof(double)*n,cudaMemcpyHostToDevice));
                    if (cuda_r==NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_r, sizeof(fftw_complex)*(hn+STRIDE_DIM*STRIDE_DIM)));
                    cutilSafeCall(cudaMemcpy(cuda_r,recip,sizeof(double)*n,cudaMemcpyHostToDevice));
                    if (cuda_cm==NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_cm, sizeof(fftw_complex)*(hn+STRIDE_DIM*STRIDE_DIM)));
                    cuda_fold_kernel<<<(hn+127)/128,128>>>(cuda_m,cuda_cxin,hn);
                    cuda_cnp_m_snp_kernel<<<(hn+127)/128,128>>>(cuda_cxin,cuda_cnp,cuda_snp, hn);
                    cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cxin,(cufftDoubleComplex *)cuda_cm,CUFFT_FORWARD));
                    if (cuda_cr==NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_cr, sizeof(fftw_complex)*(hn+STRIDE_DIM*STRIDE_DIM)));
                    cuda_fold_kernel<<<(hn+127)/128,128>>>(cuda_r,cuda_cxin,hn);
                    cuda_cnp_m_snp_kernel<<<(hn+127)/128,128>>>(cuda_cxin,cuda_cnp,cuda_snp, hn);
                    cufftSafeCall(cufftExecZ2Z(cuda_fwpx,(cufftDoubleComplex *)cuda_cxin,(cufftDoubleComplex *)cuda_cr,CUFFT_FORWARD));
                }
                if (!cufftonly) {
                    if (cuda_x==NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_x, sizeof(fftw_complex)*(hn+STRIDE_DIM*STRIDE_DIM)));
                    if (cuda_y==NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_y, sizeof(fftw_complex)*(hn+STRIDE_DIM*STRIDE_DIM)));
                }
	}
        else 
	{
                if (xin==NULL)
                    xin = (double *) malloc (n * sizeof(double));
                if (yin==NULL)
                    yin = (double *) malloc (n * sizeof(double));
                if (cxout==NULL)
                    cxout = (fftw_complex *) malloc ((hn+1) * sizeof(fftw_complex));
                if (cyout==NULL)
                    cyout = (fftw_complex *) malloc ((hn+1) * sizeof(fftw_complex));
		cufftSafeCall(cufftPlan1d(&cuda_fwpx,FFTLEN, CUFFT_D2Z, 1));
		cufftSafeCall(cufftPlan1d(&cuda_bwpx,FFTLEN, CUFFT_Z2D, 1));
                if (generic && !cufftonly) {
                    if (cuda_tmp_g==NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_tmp_g, sizeof(double)*(n+STRIDE_DIM*STRIDE_DIM)));
                    if (cuda_m==NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_m, sizeof(double)*(n+STRIDE_DIM*STRIDE_DIM)));
                    cutilSafeCall(cudaMemcpy(cuda_m,modulus,sizeof(double)*n,cudaMemcpyHostToDevice));
                    if (cuda_r==NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_r, sizeof(double)*(n+STRIDE_DIM*STRIDE_DIM)));
                    cutilSafeCall(cudaMemcpy(cuda_r,recip,sizeof(double)*n,cudaMemcpyHostToDevice));
                    if (cuda_cm==NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_cm, sizeof(double)*2*(hn+1)));
                    cufftSafeCall(cufftExecD2Z(cuda_fwpx,(cufftDoubleReal *)cuda_m,(cufftDoubleComplex *)cuda_cm));
                    if (cuda_cr==NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_cr, sizeof(double)*2*(hn+1)));
                    cufftSafeCall(cufftExecD2Z(cuda_fwpx,(cufftDoubleReal *)cuda_r,(cufftDoubleComplex *)cuda_cr));
                }
                if (!cufftonly) {
                    if (cuda_x==NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_x, sizeof(double)*(n+STRIDE_DIM*STRIDE_DIM)));
                    if (cuda_y==NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_y, sizeof(double)*(n+STRIDE_DIM*STRIDE_DIM)));
                    if (cuda_cxin==NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_cxin, sizeof(fftw_complex)*(hn+STRIDE_DIM*STRIDE_DIM)));
                    if (cuda_cyin==NULL)
                        cutilSafeCall(cudaMalloc((void**)&cuda_cyin, sizeof(fftw_complex)*(hn+STRIDE_DIM*STRIDE_DIM)));
                }
                if (cuda_xin==NULL)
                    cutilSafeCall(cudaMalloc((void**)&cuda_xin, sizeof(double)*(n+STRIDE_DIM*STRIDE_DIM)));
                if (cuda_yin==NULL)
                    cutilSafeCall(cudaMalloc((void**)&cuda_yin, sizeof(double)*(n+STRIDE_DIM*STRIDE_DIM)));
                if (cuda_cxout==NULL)
                    cutilSafeCall(cudaMalloc((void**)&cuda_cxout, sizeof(double)*2*(hn+1)));
                if (cuda_cyout==NULL)
                    cutilSafeCall(cudaMalloc((void**)&cuda_cyout, sizeof(double)*2*(hn+1)));
	}
}

// Get the pointers to the user output functions

void gwypsetoutputs (
	void(*screenf)(char *),
	void(*bothf)(char *)
)
{
	screen_output = screenf;
	both_output = bothf;
}

// Initialize the gwypnum system

int
gwypsetup(
	double	k,		// The multiplier
	unsigned long	b,	// The base (force generic reduction if not two)
	unsigned long 	n,	// The exponent
	signed long 	c,	// c, in k*b^n+c (force generic reduction if not +1 or -1)
	giant modulus_arg	// The modulus of the modular reduction
)
{
	long j, len, g, limit;
	double 	log2;
	double	log2k, log2b;
	double	tc1 = 12345.6789, tc2 = 6789.12345;

	printfunction = (verbose)? both_output : screen_output;
	gmodulus = modulus_arg;

	if (abs(c) != 1)
		generic = 1;
	else
		plus = compl2 = (c==1)? 1 : 0;

        log2 = log(2.0);
	log2b = log(b)/log2;

#define logb(x) log(x)/log(b)

	if (k > 0.0) {
		log2k = log(k)/log2;
		r = (double)n*log2b + log2k;
	}
	else
		r = (double)n*log2b;

	if (b != 2 || k == 0.0 || k*MAXMULCONST > MAX_ZEROPAD_K) {
		generic = 1;
	}

	MAXERR = 0.0;

	BIGVAL = 3.0;
	while (RINTP(tc1) != RINT(tc1) || RINTP(tc2) != RINT(tc2)) {
		BIGVAL *= 2.0;
	}
	while (RINTP(tc1) == RINT(tc1) && RINTP(tc2) == RINT(tc2)) {
		BIGVAL *= 2.0;
	}
	BIGVAL /= 2.0;

/* Calculate the number of bits in k*2^n.  This will be helpful in */
/* determining how much memory to allocate for giants. */

	bit_length = bitlen (gmodulus);

	if (debug) {
		sprintf (gwypbuf, "k = %14.1f, b = %lu, n = %lu, c = %ld, bit_length = %lu\n", k, b, n, c, bit_length);
		if (printfunction != NULL)
			(*printfunction)(gwypbuf);
	}


	ng = n;
	kg = k;
	zp = set_fftlen(k, b, n);

	if(FFTLEN < 128 && ((k!=1) || (b!=4) || (c!=1)))
	{
            return(1);  // Make only an APRCL test...
 	}

	avg_num_b_per_word = ((zp || generic) ? n * 2.0 : (logb(k) + n)) / FFTLEN;



	if (zp) {
		inc = c;		// copy for zp
		bits = ng/FFTLEN;
		limit_high = (double)(1<<bits);
		limit_inverse_high = 1.0/limit_high;
		limit_high_bigval = (limit_high*BIGVAL2)-BIGVAL2;
		mult = k;		// copy for zp
		invmult = 1.0/mult;
		rem = (ng - n)%bits;
		shift = (double) (1L << rem) * -inc;
		hwcount = (ng - n)/bits;
		lwcount = (((ng - n)/bits + 9)/8)*8;
		hwoffset = (n + bits - 1) / bits;
		if (debug) {
			sprintf (gwypbuf, "bits = %lu, ng = %lu, hwcount = %lu, lwcount = %lu, hwoffset = %lu\n", bits, ng, hwcount, lwcount, hwoffset);
			if (printfunction != NULL)
				(*printfunction)(gwypbuf);
			sprintf (gwypbuf, "limit_high = %g, limit_inverse_high = %g, limit_high_bigval = %g, mult = %g, invmult = %g, rem = %lu, shift = %g\n",
				limit_high, limit_inverse_high, limit_high_bigval, mult, invmult, rem, shift);
			if (printfunction != NULL)
				(*printfunction)(gwypbuf);
		}
		temp = (long) malloc ((hwcount + 24) * sizeof (double) + 256);
		scr = (double *) temp;					// address of scratch area
		scral = (double *) ((temp + 7) & 0xFFFFFFFFFFFFFFF8);	// id. double word aligned
		k = kg;
		n = ng;
		c = (zcomplex)? 1 : -1;     // Optional, defaulted to TRUE.
		plus = compl2 = (c == 1)? 1 : 0;
	}
	else if (generic) {
		k = kg;
		n = ng;
		c = (zcomplex)? 1 : -1;     // Optional, defaulted to TRUE.
		plus = compl2 = (c == 1)? 1 : 0;
		bits = ng/FFTLEN;
		limit_high = (double)(1<<bits);
		limit_inverse_high = 1.0/limit_high;
		limit_high_bigval = (limit_high*BIGVAL2)-BIGVAL2;
	}
	
	gwfft_weight_setup (0, k, n, c, FFTLEN);
	len = (FFTLEN+1)*sizeof(double);
	if (fftbase == NULL)
		fftbase = (int *)malloc((FFTLEN+1)*sizeof(int));
	if (two_to_phi == NULL)
		two_to_phi = (double *)malloc(len);
	if (two_to_minusphi == NULL)
		two_to_minusphi = (double *)malloc(len);
	if (invlimit == NULL)
		invlimit = (double *)malloc(len);
	if (flimit == NULL)
		flimit = (double *)malloc(len);
	if (hlimit == NULL)
		hlimit = (double *)malloc(len);
	if (limitbv == NULL)
		limitbv = (double *)malloc(len);
        if (!cufftonly) {
            if (cuda_two_to_phi==NULL)
                cutilSafeCall(cudaMalloc((void**)&cuda_two_to_phi, len+STRIDE_DIM*STRIDE_DIM*sizeof(double)));
            if (cuda_two_to_minusphi==NULL)
                cutilSafeCall(cudaMalloc((void**)&cuda_two_to_minusphi, len+STRIDE_DIM*STRIDE_DIM*sizeof(double)));
            if (g_invlimit==NULL)
                cutilSafeCall(cudaMalloc((void**)&g_invlimit, len+STRIDE_DIM*STRIDE_DIM*sizeof(double)));
            if (g_limitbv==NULL)
                cutilSafeCall(cudaMalloc((void**)&g_limitbv, len+STRIDE_DIM*STRIDE_DIM*sizeof(double)));
            if (cuda_tmp==NULL)
                cutilSafeCall(cudaMalloc((void**)&cuda_tmp, len+STRIDE_DIM*STRIDE_DIM*sizeof(double)));

            if (g_carry==NULL)
                cutilSafeCall(cudaMalloc((void**)&g_carry,len));
            if (g_err==NULL)
                cutilSafeCall(cudaMalloc((void**)&g_err,sizeof(float)*FFTLEN));
        }

        l_err = (float *)malloc(sizeof(float)*FFTLEN);

	high = (double)(1<<gwfft_base(1));
	low = 0.5*high;
	last = 1<<(gwfft_base(FFTLEN)-gwfft_base(FFTLEN-1));
	highinv = 1.0/high;
	lowinv = 1.0/low;
	lastinv = 1.0/last;
        
	if (debug) {
		sprintf (gwypbuf, "clog2k = %d, log2k = %7.4f\n", (int)ceil(log2k), log2k);
		if (printfunction != NULL)
			(*printfunction)(gwypbuf);
	}
	
	wrapindex = 0;
	wrapfactor = 1.0;
	g = 1;
        
	if (k != 1.0) {
		g = (1<<(int)ceil(log2k))-(int)k;
		while (n > gwfft_base(wrapindex))
			wrapindex++;
		wrapindex--;
		wrapfactor = g;
		for (j=0;j+gwfft_base(wrapindex)<n;j++)
			wrapfactor *= 2.0;
	}
	
	if (debug) {
		sprintf (gwypbuf, "wrapindex = %d, wrapfactor = %f\n", wrapindex, wrapfactor);
		if (printfunction != NULL)
			(*printfunction)(gwypbuf);
		sprintf(gwypbuf, "INIT : r = %7.4f, low = %7.4f, high = %7.4f, last = %7.4f, g = %ld, B = %g\n",
		r, low, high, last, g, BIGVAL2); // JP 08/07/17
		if (printfunction != NULL)
			(*printfunction)(gwypbuf);
	}
	
	for(j=0; j<FFTLEN; ++j) {
		two_to_phi[j] = gwfft_weight (j);
		two_to_minusphi[j] = (compl2) ? 2.0*gwfft_weight_inverse_over_fftlen (j) :
						gwfft_weight_inverse_over_fftlen (j);
		fftbase[j] = gwfft_base(j+1) - gwfft_base(j);
	}
	
	if (!cufftonly) {
            cutilSafeCall(cudaMemcpy(cuda_two_to_phi,two_to_phi,len,cudaMemcpyHostToDevice));
            cutilSafeCall(cudaMemcpy(cuda_two_to_minusphi,two_to_minusphi,len,cudaMemcpyHostToDevice));
            for(j=0; j<FFTLEN; ++j) l_err[j]=0;
            cutilSafeCall(cudaMemcpy(g_err,l_err,FFTLEN*sizeof(float),cudaMemcpyHostToDevice));
        }

        ttmp = two_to_minusphi[0];
	flimit[0] = high;
	invlimit[0] = highinv;
	hlimit[0] = low;
	limitbv[0] = high*BIGVAL2-BIGVAL2;
	flimit[FFTLEN-1] = last;
	invlimit[FFTLEN-1] = lastinv;
	hlimit[FFTLEN-1] = 0.5*last;
	limitbv[FFTLEN-1] = last*BIGVAL2-BIGVAL2;
        
	for(j=1; j<FFTLEN-1; ++j) {
		limit = 1<<(gwfft_base(j+1) - gwfft_base(j));
		flimit[j] = limit;
		invlimit[j] = 1.0/limit;
		hlimit[j] = 0.5*limit;
		limitbv[j] = limit*BIGVAL2-BIGVAL2;
	}

	if (generic) {
		modulus = gwypalloc ();
		recip = gwypalloc ();
		gwyptmp = gwypalloc ();
		gianttogwyp (grecip, recip);
		gianttogwyp (gmodulus, modulus);
	}
	
        dim3 grid(STRIDE_DIM/BLOCK_DIM,STRIDE_DIM/BLOCK_DIM, 1);
        dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);
        
	if (!cufftonly) {
            cutilSafeCall(cudaMemcpy(cuda_tmp,invlimit,len,cudaMemcpyHostToDevice));
//            dim3 grid(STRIDE_DIM/BLOCK_DIM,STRIDE_DIM/BLOCK_DIM, 1);
//            dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);
            for(j=0;j<FFTLEN;j+=(STRIDE_DIM*STRIDE_DIM))
           	transpose<<<grid, threads>>>((double *)&g_invlimit[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
            cutilSafeCall(cudaMemcpy(cuda_tmp,limitbv,len,cudaMemcpyHostToDevice));
            for(j=0;j<FFTLEN;j+=(STRIDE_DIM*STRIDE_DIM))
           	transpose<<<grid, threads>>>((double *)&g_limitbv[j],(double *)&cuda_tmp[j],(int) STRIDE_DIM,(int) STRIDE_DIM);
        }
        
	limitbv[0] = high*BIGVAL2-BIGVAL2; // remplacé BIGVAL par BIGVAL2 JP 04/07/17
	limitbv[FFTLEN-1] = last*BIGVAL2-BIGVAL2;
        
	for(j=1; j<FFTLEN-1; ++j) {
		limit = 1<<(gwfft_base(j+1) - gwfft_base(j));
		limitbv[j] = limit*BIGVAL2-BIGVAL2;
	}

        init_fftw((FFTLEN < STRIDE_DIM*2)?STRIDE_DIM*2:FFTLEN);
        
	addinindex = 0;
	addinvalue = 0.0;
	SMALLMULCONST = 1.0;
	MULBYCONST = 0;
        
	return 0;
}

int gwypsetup_general_mod_giant (
	giant modulus_arg	// The modulus of the modular reduction
)
{
	return gwypsetup (0.0, 1, 1, -1, modulus_arg);
}

//	Miscellanous utility routines

void bitaddr (
	unsigned long bit,
	unsigned long *word,
	unsigned long *bit_in_word,
	unsigned long FFTLEN)
{

/* What word is the bit in? */

	*word = (unsigned long) ((double) bit / r * FFTLEN);
	if (*word >= FFTLEN) *word = FFTLEN - 1;

/* Compute the bit within the word. */

	*bit_in_word = bit - (int)ceil(*word*r/FFTLEN);
}

void setaddin (
	long	value, int N)
{
	unsigned long word, bit_in_word;

	if (zp || generic) {
		addinvalue = (double) value;
		addinindex = 0;
		return;
	}

/* If value is even, shift it right and increment bit number.  This */
/* will ensure that we modify the proper FFT word. */

	for (bit_in_word = 0; value && (value & 1) == 0; value >>= 1)
		bit_in_word++;

/* Convert the input value to 1/k format.  Case 1 (2^n+/-1: Inverse of k */
/* is 1.  Case 2 (k*2^n-1): Inverse of k is 2^n.  Case 3 (k*2^n+1): Inverse */
/* of k is -2^n.  No other cases can be handled. */

	if (kg == 1.0) {
		bitaddr (bit_in_word, &word, &bit_in_word, N);
	}
	else {
		bitaddr (ng + bit_in_word, &word, &bit_in_word, N);
	}
	addinvalue = (plus && (kg != 1.0)) ? -(double)(value<<bit_in_word) : (double)(value<<bit_in_word);
	addinindex = word;
}

/* Add a small constant at the specified power of b after the */
/* next multiplication.  That is, value*b^power_of_b is added to */
/* the next multiplication result.  This only works if k=1. */

void gwypsetaddinatpowerofb (
	long	value,
	unsigned long power_of_b_arg)
{
	unsigned long word, b_in_word, power_of_b;

/* If value is even, shift it right and increment bit number.  This */
/* will ensure that we modify the proper FFT word. */

	for (power_of_b = power_of_b_arg; value && (value & 1) == 0; value >>= 1)
		power_of_b++;

	bitaddr (power_of_b, &word, &b_in_word, FFTLEN);

	addinvalue = (double)(value<<b_in_word);
	addinindex = word;
}

/* Test if a gwypnum is zero */

#define	MAX_NZ_COUNT 10

int
gwypiszero(
	gwypnum 	gg
)
{
	register long 		j, count = 0;
	register double 	*gp = gg;
	long	result;
	giant	gtest;

	if (!generic)
		for(j=0; j<FFTLEN; j++) {

			if (count > MAX_NZ_COUNT)
				return 0;		// Too much non zero words, the gwypnum is not zero.
			else if (*gp++)
				count++;		// Update the count of non-zero words.
		}
	if (count || generic) {		// The gwypnum is likely to be zero but needs a more accurate test...
		gtest = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);// Allocate memory for the giant
		gwyptogiant (gg, gtest);
		result = isZero (gtest);
		gwypfree (gtest);			// Free memory
		return (result);
	}
	else
		return 1;			// The gwypnum is zero
}

/* Test two gwypnums for equality */

int	gwypequal (
	gwypnum gw1, 
	gwypnum gw2
) 
{
	gwypnum gwdiff;
	int result;

	gwdiff = gwypalloc ();			// Reserve memory for the difference
	gwypsub3 (gw1, gw2, gwdiff);		// Normalized subtract...
	result = gwypiszero (gwdiff);	// Test for zero difference
	gwypfree (gwdiff);				// Free memory
	return (result);
}


//	User side functions...

double
gwypnormalize(
	gwypnum s
)
{
	double err;

	if (debug)
		gwypstart_timer (2);
	err = (zp || generic)? rnormalize(s, FFTLEN, E_CHK, 0, 0) : inormalize(s, FFTLEN, E_CHK, 0, 0);
	if (debug)
		gwypend_timer (2);
	return (err);
}


double
gwyprawnormalize(
	gwypnum s
)
{
	double err;

	if (debug)
		gwypstart_timer (2);
	err = (zp || generic)? rnormalize(s, FFTLEN, E_CHK, 1, 1) : inormalize(s, FFTLEN, E_CHK, 1, 1);
	if (debug)
		gwypend_timer (2);
	return (err);
}

void cuda_generic_modred(
	gwypnum s,
        int flag
)
{
	double err;

	gwypcopyzero (s, gwyptmp, zerowordslow);
	err = cuda_lucas_mul (recip, gwyptmp, FFTLEN, E_CHK, 1, 1, flag);
	if (err > MAXERR)
		MAXERR = err;
	gwypsetzero (gwyptmp, zerowordshigh);
	err = cuda_lucas_mul (modulus, gwyptmp, FFTLEN, E_CHK, 1, 1, flag);
	if (err > MAXERR)
		MAXERR = err;
	if (debug)
		gwypstart_timer (4);
	if (compl2) {
		gwypaddquick (gwyptmp, s);
		gwypsubquick (modulus, s);
	}
	else
		gwypsubquick (gwyptmp, s);
	if (debug)
		gwypend_timer (4);
}

void generic_modred(
	gwypnum s
)
{
	double err;

	gwypcopyzero (s, gwyptmp, zerowordslow);
	err = lucas_mul (recip, gwyptmp, FFTLEN, E_CHK, 1, 1);
	if (err > MAXERR)
		MAXERR = err;
	gwypsetzero (gwyptmp, zerowordshigh);
	err = lucas_mul (modulus, gwyptmp, FFTLEN, E_CHK, 1, 1);
	if (err > MAXERR)
		MAXERR = err;
	if (debug)
		gwypstart_timer (4);
	if (compl2) {
		gwypaddquick (gwyptmp, s);
		gwypsubquick (modulus, s);
	}
	else
		gwypsubquick (gwyptmp, s);
	if (debug)
		gwypend_timer (4);
}

void gwypfft_description (
	char *buf)
{
	if (zp)
            if (compl2)
		sprintf (buf, "Using complex zero-padded rational base DWT, FFT length = %d", FFTLEN);
            else
		sprintf (buf, "Using real zero-padded rational base DWT, FFT length = %d", FFTLEN);
	else if (generic)
            if (compl2)
		sprintf (buf, "Using complex rational base DWT and generic reduction, FFT length = %d", FFTLEN);
            else
		sprintf (buf, "Using real rational base DWT and generic reduction, FFT length = %d", FFTLEN);
	else if (compl2)
            sprintf (buf, "Using complex irrational base DWT, FFT length = %d", FFTLEN);
	else
            sprintf (buf, "Using real irrational base DWT, FFT length = %d", FFTLEN);
}



// User side large integers arithmetic operations


void gwypadd (
	gwypnum s,
	gwypnum d)
{
	int i;

	for (i=0; i<FFTLEN; i++)
		d[i] += s[i];
	gwyprawnormalize (d);
//	if (zp)
//		modred (d);
//	else if (generic)
//		generic_modred (d);
}


void gwypsub (
	gwypnum s,
	gwypnum d)
{
	int i;

	for (i=0; i<FFTLEN; i++) {
		d[i] -= s[i];
	}
	gwyprawnormalize (d);
//	if (zp)
//		modred (d);
//	else if (generic)
//		generic_modred (d);
}

void gwypadd3 (
	gwypnum s1,
	gwypnum s2,
	gwypnum d)
{
	int i;

	for (i=0; i<FFTLEN; i++)
		d[i] = s1[i] + s2[i];
	gwyprawnormalize (d);
//	if (zp)
//		modred (d);
//	else if (generic)
//		generic_modred (d);
}


void gwypsub3 (
	gwypnum s1,
	gwypnum s2,
	gwypnum d)
{
	int i;

	for (i=0; i<FFTLEN; i++)
		d[i] = s1[i] - s2[i];
	gwyprawnormalize (d);
//	if (zp)
//		modred (d);
//	else if (generic)
//		generic_modred (d);
}

void
cuda_gwypsquare (gwypnum s, int flag) {
	double err;
        if (generic)
            err = cuda_lucas_square_generic (s, FFTLEN, E_CHK, 0, 0, flag);
        else
            err = cuda_lucas_square (s, FFTLEN, E_CHK, 0, 0, flag);
	if (err > MAXERR)
		MAXERR = err;
//	if (generic)
//		cuda_generic_modred (s, flag);
//                generic_modred (s);
}

void
gwypsquare2 (gwypnum s) {
	double err;

	err = lucas_square (s, FFTLEN, E_CHK, 0, 0);
	if (err > MAXERR)
		MAXERR = err;
	if (generic)
		generic_modred (s);
}

void
gwypsquare (gwypnum s) {
	double err;

	err = lucas_square (s, FFTLEN, E_CHK, 0, 0);
	if (err > MAXERR)
		MAXERR = err;
	if (generic)
		generic_modred (s);
}

void
cuda_gwypmul (gwypnum s, gwypnum d, int flag) {
	double err;
        if (generic)
            err = cuda_lucas_mul_generic (s, d, FFTLEN, E_CHK, 0, 0, flag);
        else
            err = cuda_lucas_mul (s, d, FFTLEN, E_CHK, 0, 0, flag);
	if (err > MAXERR)
		MAXERR = err;
//	if (generic)
//		cuda_generic_modred (d, flag);
//                generic_modred (d);
}

void
gwypmul (gwypnum s, gwypnum d) {
	double err;

	err = lucas_mul (s, d, FFTLEN, E_CHK, 0, 0);
	if (err > MAXERR)
		MAXERR = err;
	if (generic)
		generic_modred (d);
}

/* Generate random FFT data */

void gwyp_random_number (
	gwypnum	x)
{
	giant	g;
	unsigned long i, len;

/* Generate the random number */

	srand ((unsigned int) time (NULL));
	len = (unsigned long) (bit_length / 16) + 1;
	g = newgiant (len);
	for (i = 0; i < len; i++) {
		g->n[i] = ((unsigned short) rand() << 10) +
			  ((unsigned short) rand() << 5) +
			  (unsigned short) rand();
	}
	g->sign = len;
	modg (gmodulus, g);
	gianttogwyp (g, x);
	gwyprawnormalize (x);
	gwypfree(g);
}


/* Square a number using a slower method that will have reduced */
/* round-off error on non-random input data.*/

void gwypsquare_carefully (
	gwypnum	s)		/* Source and destination */
{
	gwypnum	tmp1, tmp2;
	double	err;

/* Generate a random number, if we have't already done so */

	if (GWP_RANDOM == NULL) {
		GWP_RANDOM = gwypalloc ();
		gwyp_random_number (GWP_RANDOM);
	}

/* Now do the squaring using three multiplies and adds */

	tmp1 = gwypalloc ();
	tmp2 = gwypalloc ();
	gwypadd3 (s, GWP_RANDOM, tmp1);		/* Compute s+random */
	gwypcopy (GWP_RANDOM, tmp2);
	err = lucas_mul (tmp2, s, FFTLEN, E_CHK, 1, 0);	/* Compute s*random without adding*/
	if (err > MAXERR)
		MAXERR = err;
	if (generic)
		generic_modred (s);
	err = lucas_square (tmp2, FFTLEN, E_CHK, 1, 0);	/* Compute random^2 without addin*/
	if (err > MAXERR)
		MAXERR = err;
	if (generic)
		generic_modred (tmp2);
	gwypsquare (tmp1);					/* Compute (s+random)^2  + addinvalue */
	gwypsubquick (tmp2, tmp1);			/* Calc s^2 from 3 results */
	gwypaddquick (s, s);
	gwypsub3 (tmp1, s, s);

/* Free memory and return */

	gwypfree (tmp1);
	gwypfree (tmp2);
}


/* Multiply numbers using a slower method that will have reduced */
/* round-off error on non-random input data.*/

void gwypmul_carefully (
	gwypnum s, gwypnum t)					/* Source and destination */
{
	gwypnum	tmp1, tmp2, tmp3, tmp4;
	double	err;

/* Generate a random number, if we have't already done so */

	if (GWP_RANDOM == NULL) {
		GWP_RANDOM = gwypalloc ();
		gwyp_random_number (GWP_RANDOM);
	}

/* Now do the multiply using four multiplies and adds */

	tmp1 = gwypalloc ();
	tmp2 = gwypalloc ();
	tmp3 = gwypalloc ();
	tmp4 = gwypalloc ();
	gwypcopy (s, tmp4);
	gwypadd3 (s, GWP_RANDOM, tmp1);		/* Compute s+random */
	gwypadd3 (t, GWP_RANDOM, tmp3);		/* Compute t+random */
	gwypcopy (GWP_RANDOM, tmp2);
	err = lucas_mul (tmp2, tmp4, FFTLEN, E_CHK, 1, 0);	/* Compute s*random without adding*/
	if (err > MAXERR)
		MAXERR = err;
	if (generic)
		generic_modred (tmp4);
	err = lucas_mul (tmp2, t, FFTLEN, E_CHK, 1, 0);		/* Compute t*random without adding*/
	if (err > MAXERR)
		MAXERR = err;
	if (generic)
		generic_modred (t);
	err = lucas_square (tmp2, FFTLEN, E_CHK, 1, 0);		/* Compute random^2 without addin*/
	if (err > MAXERR)
		MAXERR = err;
	if (generic)
		generic_modred (tmp2);
	err = lucas_mul (tmp1, tmp3, FFTLEN, E_CHK, 0, 0);	/* Compute (s+random)*(t+random) + addinvalue */
	if (err > MAXERR)
		MAXERR = err;
	if (generic)
		generic_modred (tmp3);
	gwypsubquick (tmp2, tmp3);				/* Subtract random^2 */
	gwypsubquick (t, tmp3);
	gwypsub3 (tmp3, tmp4, t);

/* Free memory and return */

	gwypfree (tmp1);
	gwypfree (tmp2);
	gwypfree (tmp3);
	gwypfree (tmp4);
}


/* Do a squaring very carefully.  This is done after a normal */
/* iteration gets a roundoff error above 0.40.  This careful iteration */
/* will not generate a roundoff error. */

void gwypcareful_squaring (
	gwypnum s
)		
{
	gwypnum	hi, lo;
	long i;
	double	err;

/* Copy the data to hi and lo.  Zero out half the FFT data in each. */

	hi = gwypalloc ();
	lo = gwypalloc ();
	gwypcopy (s, hi);
	gwypcopy (s, lo);
	for (i = 0; i < FFTLEN/2; i++)
		hi[i] = 0;
	for ( ; i < FFTLEN; i++)
		lo[i] = 0;

/* Now do the squaring using three multiplies and adds */

	gwypcopy (hi, s);
	err = lucas_mul (lo, s, FFTLEN, E_CHK, 1, 0);	// Compute hi*lo without adding
	if (err > MAXERR)
		MAXERR = err;
	if (generic)
		generic_modred (s);
	gwypsquare (hi);									// Compute hi^2 + addin value
	err = lucas_square (lo, FFTLEN, E_CHK, 1, 0);	// Compute lo^2 without addin
	if (err > MAXERR)
		MAXERR = err;
	if (generic)
		generic_modred (lo);
	gwypaddquick (s, s);				// 2*hi*lo
	gwypaddquick (hi, s);			// + hi^2
	gwypadd (lo, s);					// + lo^2

/* Free memory and return */

	gwypfree (hi);
	gwypfree (lo);
}


void gwypcareful_multiply (
	gwypnum s,
	gwypnum t
)		
{
	gwypnum	hi, lo, hi2, lo2, u;
	long i;
	double	err;

/* Copy the data to hi, lo, hi2 and lo2.  Zero out half the FFT data in each. */

	hi = gwypalloc ();
	lo = gwypalloc ();
	hi2 = gwypalloc ();
	lo2 = gwypalloc ();
	u = gwypalloc ();

	gwypcopy (s, hi);
	gwypcopy (s, lo);
	gwypcopy (t, hi2);
	gwypcopy (t, lo2);

	for (i = 0; i < FFTLEN/2; i++) {
		hi[i] = 0;
		hi2[i] = 0;
	}

	for ( ; i < FFTLEN; i++) {
		lo[i] = 0;
		lo2[i] = 0;
	}

/* Now do the multiply using four multiplies and adds */

	gwypcopy (hi2, t);
	gwypmul (hi, t);									// Compute hi*hi2 + addin value
	gwypcopy (lo2, u);
	gwypmul (lo, u);
	err = lucas_mul (lo, u, FFTLEN, E_CHK, 1, 0);	// Compute lo*lo2 without adding
	if (err > MAXERR)
		MAXERR = err;
	if (generic)
		generic_modred (u);
	gwypmul (lo2, hi);
	err = lucas_mul (lo2, hi, FFTLEN, E_CHK, 1, 0);	// Compute lo2*hi without adding
	if (err > MAXERR)
		MAXERR = err;
	if (generic)
		generic_modred (hi);
	gwypmul (hi2, lo);
	err = lucas_mul (hi2, lo, FFTLEN, E_CHK, 1, 0);	// Compute hi2*lo without adding
	if (err > MAXERR)
		MAXERR = err;
	if (generic)
		generic_modred (lo);
	gwypaddquick (u, t);
	gwypaddquick (hi, t);
	gwypadd (lo, t);

/* Free memory and return */

	gwypfree (hi);
	gwypfree (lo);
	gwypfree (hi2);
	gwypfree (lo2);
	gwypfree (u);
}

// Set small add-in constant

void gwypsetaddin(
	long s
)
{
	setaddin (s, FFTLEN);
}

// Set small multiplicative constant

void gwypsetmulbyconst(
	long s
)
{
	SMALLMULCONST = (double)s;
}

// Set the maximum of the multiplicative constant

void gwypsetmaxmulbyconst(
	long s
)
{
	MAXMULCONST = (double)s;
}

// Conversion routines

void itogwyp(					// Set a gwypnum to a small value
	int s,
	gwypnum d
)
{
	int j, saveindex;
	double savevalue;

	for (j=0; j<FFTLEN; j++)	// Init large integer to zero
		d[j] = 0.0;
	saveindex = addinindex;		// Save setaddin internal status
	savevalue = addinvalue;
	setaddin(s, FFTLEN);
	d[addinindex] = addinvalue;	// Set the large integer value in 1/k format
	addinindex = saveindex;		// Restore setaddin
	addinvalue = savevalue;
}

void gwypaddsmall(				// Add a small value to a gwypnum
	gwypnum d,
	int s
)
{
	int saveindex;
	double savevalue;

	saveindex = addinindex;		// Save setaddin internal status
	savevalue = addinvalue;
	setaddin (s, FFTLEN);
	gwypnormalize(d);			// Normalize the number while adding the value
	addinindex = saveindex;		// Restore setaddin
	addinvalue = savevalue;
}

/* Convert a giant to gwypnum FFT format.  Giant must be a positive number. */

void gianttogwyp (
	giant	a,
	gwypnum	g)
{
	giant	newg;
	unsigned e1len;
	int	i, bits_in_next_binval;
	unsigned long binval, carry;
	unsigned short *e1;

/* To make the mod k*b^n+c step faster, gwypnum's are pre-multiplied by 1/k */
/* If k is greater than 1, then we calculate the inverse of k, multiply */
/* the giant by the inverse of k, and do a mod k*b^n+c. */

	if (kg > 1.0) {

		newg = newgiant (((unsigned long)( bit_length / 16) + 1) * 2);

		/* Easy case 1 (k*2^n-1): Inverse of k is 2^n */

		if (!plus) {
			gtog (a, newg);
			gshiftleft (ng, newg);
		}

		/* Easy case 2 (k*2^n+1): Inverse of k is -2^n */

		else {
			gtog (gmodulus, newg);	// replace -a by gmodulus-a to have positive numbers!!
			subg (a, newg);
			gshiftleft (ng, newg);
		}

		modg (gmodulus, newg);

		a = newg;
	}

/* Now convert the giant to FFT format */

	ASSERTG (a->sign >= 0);
	e1len = a->sign;
	e1 = a->n;

	if (e1len) {binval = *e1++; e1len--; bits_in_next_binval = 16;}
	else binval = 0;
	carry = 0;
	for (i = 0; i < FFTLEN; i++) {
		int	bits;
		long	value, mask;
		bits = (zp || generic)? ng/FFTLEN : fftbase [i];
		mask = (1L << bits) - 1;
		if (i == FFTLEN - 1) value = binval;
		else value = binval & mask;
		value = value + carry;
		if (value > (mask >> 1) && bits > 1 && i != (FFTLEN - 1)) {
			value = value - (mask + 1);
			carry = 1;
		} else {
			carry = 0;
		}
		g[i] = (double)value;

		binval >>= bits;
		if (e1len == 0) continue;
		if (bits_in_next_binval < bits) {
			if (bits_in_next_binval)
				binval |= (*e1 >> (16 - bits_in_next_binval)) << (16 - bits);
			bits -= bits_in_next_binval;
			e1++; e1len--; bits_in_next_binval = 16;
			if (e1len == 0) continue;
		}
		if (bits) {
			binval |= (*e1 >> (16 - bits_in_next_binval)) << (16 - bits);
			bits_in_next_binval -= bits;
		}
	}

/* Free allocated memory */

	if (kg > 1.0)
		gwypfree(newg);
}


int gwyptogiant (
	gwypnum	gg,
	giant	v)
{
	long	val;
	int	i, j, limit, bits, bitsout, carry;
	unsigned short *outptr;

	limit = FFTLEN;

/* Collect bits until we have all of them */

	carry = 0;
	bitsout = 0;
	outptr = v->n;
	*outptr = 0;
	if (zp || generic)
		bits = ng / FFTLEN;

	for (i = 0; i < limit; i++) {
		val = (long) gg[i];
		if (!zp && !generic)
			bits = fftbase[i];
		val += carry;
		for (j = 0; j < bits; j++) {
			*outptr >>= 1;
			if (val & 1) *outptr += 0x8000;
			val >>= 1;
			bitsout++;
			if (bitsout == 16) {
				outptr++;
				bitsout = 0;
			}
		}
		carry = val;
	}

/* Finish outputting the last word and any carry data */

	while (bitsout || (carry != -1 && carry != 0)) {
		*outptr >>= 1;
		if (carry & 1) *outptr += 0x8000;
		carry >>= 1;
		bitsout++;
		if (bitsout == 16) {
			outptr++;
			bitsout = 0;
		}
	}

/* Set the length */

	v->sign = (long) (outptr - v->n);
	while (v->sign && (v->n[v->sign-1] == 0)) v->sign--;

/* If carry is -1, the gwnum is negative.  Ugh.  Flip the bits and sign. */
	
	if (carry == -1) {
		for (j = 0; j < v->sign; j++) v->n[j] = ~v->n[j];
		while (v->sign && (v->n[v->sign-1] == 0)) v->sign--;
		iaddg (1, v);
		v->sign = -v->sign;
	}

/* The gwnum is not guaranteed to be smaller than k*b^n+c.  Handle this */
/* possibility.  This also converts negative values to positive. */

	modg (gmodulus, v);

/* Since all gwnums are premultiplied by the inverse of k, we must now */
/* multiply by k to get the true result. */

	if (kg > 1.0) {
		giant	newg;
		newg = newgiant ((unsigned long) (bit_length / 16) + 3);
		itog ((int)kg, newg);
		mulg (v, newg);
		modg (gmodulus, newg);
		gtog (newg, v);
		gwypfree (newg);
	}

/* Return success */

	return (0);
}

char	timebuf[40];

void	gwypdone (void) {				// Free all the memory used by this code.
	printfunction = (verbose)? both_output : screen_output;
	MAXMULCONST = 1.0;
	gwypfree (fftbase);
	gwypfree (two_to_phi);
	gwypfree (two_to_minusphi);
	gwypfree (invlimit);
	gwypfree (flimit);
	gwypfree (hlimit);
	gwypfree (limitbv);
	gwypfree (cn);
	gwypfree (sn);
	gwypfree (permute);
	gwypfree (permutedx);
	gwypfree (permutedy);	
	gwypfree (ip);
	gwypfree (w);
	gwypfree (cnp);
	gwypfree (snp);
	if (zp)
		gwypfree (scr);
	gwypfree (GWP_RANDOM);
	gwypfree (modulus);
	gwypfree (recip);
	gwypfree (gwyptmp);
	gwypfree (grecip);
        cufftSafeCall(cufftDestroy(cuda_fwpx));
	if(!compl2) {
        	cufftSafeCall(cufftDestroy(cuda_bwpx));
        }
        if (!cufftonly) {
            cutilSafeCall(cudaFree((char *)cuda_x));
            cuda_x = NULL;
            cutilSafeCall(cudaFree((char *)cuda_y));
            cuda_y = NULL;
            cutilSafeCall(cudaFree((char *)cuda_tmp));
            cuda_tmp = NULL;
        }
        if (!compl2) {
            cutilSafeCall(cudaFree((char *)cuda_xin));
            cuda_xin=NULL;
            cutilSafeCall(cudaFree((char *)cuda_yin));
            cuda_yin=NULL;
        }
        else {
            cutilSafeCall(cudaFree((char *)cuda_cxin));
            cuda_cxin=NULL;
            cutilSafeCall(cudaFree((char *)cuda_cyin));
            cuda_cyin=NULL;
        }
        cutilSafeCall(cudaFree((char *)cuda_cxout));
        cuda_cxout=NULL;
        cutilSafeCall(cudaFree((char *)cuda_cyout));
        cuda_cyout=NULL;
        if (!cufftonly) {
            cutilSafeCall(cudaFree((char *)cuda_two_to_phi));            cuda_two_to_phi = NULL;
            cutilSafeCall(cudaFree((char *)cuda_two_to_minusphi));
            cuda_two_to_minusphi = NULL;
            cutilSafeCall(cudaFree((char *)g_carry));
            g_carry = NULL;
            cutilSafeCall(cudaFree((char *)g_err));
            g_err = NULL;
            cutilSafeCall(cudaFree((char *)cuda_cnp));
            cuda_cnp = NULL;
            cutilSafeCall(cudaFree((char *)cuda_snp));
            cuda_snp = NULL;
            free (l_err);
            if (generic) {
                cutilSafeCall(cudaFree((char *)cuda_tmp_g));
                cuda_tmp_g = NULL;
                cutilSafeCall(cudaFree((char *)cuda_m));
                cuda_m = NULL;
                cutilSafeCall(cudaFree((char *)cuda_r));
                cuda_r = NULL;
                cutilSafeCall(cudaFree((char *)cuda_cm));
                cuda_cm = NULL;
                cutilSafeCall(cudaFree((char *)cuda_cr));
                cuda_cr = NULL;
            }
        }

	if (debug) {
		if (MAXERR != 0.0) {
			sprintf (gwypbuf, "Maximum Round off : %10.10f\n", MAXERR);
			if (printfunction != NULL)
				(*printfunction)(gwypbuf);
		}
		gwypwrite_timer (timebuf, 5, TIMER_CLR);
		sprintf (gwypbuf, "Squaring and/or Mult. time : %s\n", timebuf); 
		if (printfunction != NULL)
			(*printfunction)(gwypbuf);
		gwypwrite_timer (timebuf, 2, TIMER_CLR);
		sprintf (gwypbuf, "Normalization time : %s\n", timebuf); 
		if (printfunction != NULL)
			(*printfunction)(gwypbuf);
		if (zp) {
			gwypwrite_timer (timebuf, 3, TIMER_CLR); 
			sprintf (gwypbuf, "Modular reduction time : %s\n", timebuf); 
			if (printfunction != NULL)
				(*printfunction)(gwypbuf);
		}
		if (generic) {
			gwypwrite_timer (timebuf, 4, TIMER_CLR); 
			sprintf (gwypbuf, "Generic copyzero + setzero time : %s\n", timebuf); 
			if (printfunction != NULL)
				(*printfunction)(gwypbuf);
		}
	}
	zp = 0;
	generic = 0;
}
