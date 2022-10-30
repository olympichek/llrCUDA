#ifndef _GWPNUMI_
#define _GWPNUMI_

/* Include Files */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#if defined (__linux__) || defined (__FreeBSD__) || defined (__APPLE__)
#include <sys/time.h>
#else
#include <time.h>
#endif
#include <sys/timeb.h>


/* compiler options */

#ifdef _WIN32
#pragma warning( disable : 4127 4706 ) /* disable conditional is constant warning */
#endif

#define TIMER_NL	0x1
#define TIMER_CLR	0x2
#define TIMER_OPT_CLR	0x4
#define NBTIMERS 10

// Global flags and data defined in gwypnum.cu file

extern	int FFTLEN;
extern	int s_FFTLEN;//cuda
extern  int g_fftlen;//cuda
extern  unsigned long maxbitsinfftlen, maxbitsinfftword;
                                    // JP 20/09/20
extern	int MULBYCONST;
extern	int zp;
extern	int generic;
extern	int inc;
extern	int debug;
extern	int tdebug;
extern	int verbose;
extern	int E_CHK;
extern  int balerr;
extern	double MAXERR;
extern	double gwyptimers[];

typedef double* gwypnum;

/* function prototypes */

#define gwypsetnormroutine(z,e,c) {E_CHK=(e);MULBYCONST=(c);}
#define gwypfree(ptr) if ((void*)ptr != NULL) {free (ptr); ptr = NULL;}
#define gwypcudaFree(ptr) if ((void*)ptr != NULL) {cutilSafeCall(cudaFree (ptr)); ptr = NULL;}
#define gwypswap(s,d)	{gwypnum t; t = s; s = d; d = t;}

void gwypclear_timers (void);

void gwypclear_timer (
	int
);

void gwypstart_timer ( 
	int
); 

void gwypend_timer ( 
	int
);
 
void gwypdivide_timer (
	int,
	int
);

double gwyptimer_value (
	int
);

void gwypprint_timer (
	int,
	int
);

void gwypwrite_timer (		// JP 23/11/07
	char*,
	int, 
	int
	);

void gwypclearline (int
);

/* ------------ gwypnum - specific routines ------------------- */

void gwypsetoutputs (	// Get the pointers to user output functions
	void(*screenf)(char *), 
	void(*bothf)(char *)
);

int
gwypsetup(		// Initialize the gwypnum system
	double,		// The multiplier
	unsigned long,	// The base (generic mode forced if not two)
	unsigned long,	// The exponent
	signed long,	// c, in k*b^n+c (force generic reduction if neither +1 nor -1)
	giant		// modulus
);

int gwypsetup_general_mod_giant (
	giant		// The modulus of the modular reduction
);

void gwypset_larger_fftlen_count(
	int
);

void gwypfft_description (
	char *
);

gwypnum gwypalloc(	// Memory allocation
);

void itogwyp(		// Integer to gwypnum
	int,
	gwypnum
);

void gwypaddsmall(	// Add a small value to a gwypnum
	gwypnum,
	int
);

void gwypsetaddin(	// Set the addin constant before normalizing
	long
);

void gwypsetaddinatpowerofb (
	long,
	unsigned long
);

void gwypsetmaxmulbyconst(// Set the maximum of the multiplicative constant
	unsigned long
);

void gwypsetmulbyconst(	// Set the multiplicative constant before normalizing
	long
);

int gwyptogiant (
	gwypnum,
	giant
);

void gianttogwyp (
	giant,
	gwypnum
);

// User side large integers arithmetic operations

void gwypcopy (
	gwypnum,
	gwypnum
);

void			// Square a large integer
cuda_gwypsquare(
	gwypnum,
	int
);

void			// Square a large integer
gwypsquare(
	gwypnum
);

void			// Multiply two large integers
cuda_gwypmul(
	gwypnum,
	gwypnum,
        int
);

void			// Multiply two large integers
gwypmul(
	gwypnum,
	gwypnum
);

void gwypaddquick (
	gwypnum,
	gwypnum
);

void gwypsubquick (
	gwypnum,
	gwypnum
);

void gwypadd (
	gwypnum s,
	gwypnum d
);

void gwypsub (
	gwypnum s,
	gwypnum d
);

void gwypadd3 (
	gwypnum s1,
	gwypnum s2,
	gwypnum d
);

void gwypsub3 (
	gwypnum s1,
	gwypnum s2,
	gwypnum d
);

/* Square a number using a slower method that will have reduced */
/* round-off error on non-random input data.*/

void gwypsquare_carefully (
	gwypnum
);

/* Multiply numbers using a slower method that will have reduced */
/* round-off error on non-random input data.*/

void gwypmul_carefully (
	gwypnum,
	gwypnum
);

int			// Test is the large integer is zero
gwypiszero(
	gwypnum
);

int	gwypequal (	// Test two gwypnums for equality
	gwypnum, 
	gwypnum
);

/********************* Internal functions ***************************/

double gwyp_get_maxerr ();
/*{
	return (MAXERR);
}*/

void gwyp_clear_maxerr ();
/*{
	MAXERR = 0.0;
}*/

double			// Normalize a large integer
gwypnormalize(
	gwypnum
);

void gwypcopyzero (
	gwypnum,
	gwypnum,
	unsigned long
);

void gwypsetzero (
	gwypnum,
	unsigned long
);

void	gwypdone (
void
);			// Free all the memory used by this code.

#endif
