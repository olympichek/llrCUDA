/*----------------------------------------------------------------------
| gwdbldbl.h
|
| This file contains the headers for the gwnum helper routines that use
| extended-precision floats.
| 
|  Copyright 2005 Just For Fun Software, Inc.
|  All Rights Reserved.
+---------------------------------------------------------------------*/

#ifndef _GWDBLDBL_H
#define _GWDBLDBL_H

/* This is a C library.  If used in a C++ program, don't let the C++ */
/* compiler mangle names. */

#ifdef __cplusplus
extern "C" {
#endif

/* Include common definitions */

#include "common.h"

/* Extended precision helper routines */

double fftcos (unsigned long);
double fftsin (unsigned long);
double fftcosp (unsigned long);
double fftsinp (unsigned long);
double fftinv2cosp (unsigned long j, unsigned long fftlen);
double fftmul (double a, double b);
double fftadd (double a, double b);
double fftsub (double a, double b);
void gwfft_weight_setup (int, double, unsigned long, signed long, unsigned long);
double gwfft_weight (unsigned long);
double gwfft_weight_sloppy (unsigned long);
double gwfft_weight_inverse (unsigned long);
double gwfft_weight_inverse_sloppy (unsigned long);
double gwfft_weight_inverse_over_fftlen (unsigned long/*, unsigned long */);
void gwfft_weights3 (unsigned long, double *, double *, double *);
double gwfft_weight_exponent (unsigned long);
unsigned long gwfft_base (unsigned long);

#ifdef __cplusplus
}
#endif

#endif
