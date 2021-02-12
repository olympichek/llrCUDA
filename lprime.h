/* Handy definitions */

#define FALSE	0
#define TRUE	1

/* This controls whether we want to pause computation if the load average */
/* becomes too great.  This does not apply to OS/2. */

#if defined (__linux__) || defined (__FreeBSD__) || defined (__APPLE__)
#define MPRIME_LOADAVG
#define max(x, y) (x > y)? x : y
void Sleep (long);

/* Handle differences between Windows and Linux runtime libraries */

#define stricmp(x,y)	strcasecmp(x,y)
#define _commit(f)	fsync(f)
#define _open		open
#define _close		close
#define _read		read
#define _write		write
#define _lseek		lseek
#define _unlink		unlink
#define _creat		creat
#define _chdir		chdir
#define _ftime		ftime
#define _timeb		timeb
#define IsCharAlphaNumeric(c) isalnum(c)
#define _O_APPEND	O_APPEND
#define _O_RDONLY	O_RDONLY
#define _O_WRONLY	O_WRONLY
#define _O_RDWR		O_RDWR
#define _O_CREAT	O_CREAT
#define _O_TRUNC	O_TRUNC
#define _O_BINARY 	0
#define _O_TEXT		0

#endif

#define EXTERNC

/* The common include files */

#include <time.h>
#include <assert.h>
extern int NO_GUI;
#include "./giants.h"
#include "./giantext.h"
#include "./gwypnum.h"
#include "./common.h"
#include "Llr.h"

/* Global variables */

extern int volatile THREAD_STOP;	/* TRUE if thread should stop */
extern int volatile THREAD_KILL;	/* TRUE if program should terminate */
extern int MENUING;			/* TRUE when main menu active */

/* Internal routines */

void main_menu ();
void linuxContinue (char *);

