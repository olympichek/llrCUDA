#ifndef LLRDOTH
#define LLRDOTH

/* Constants */

#define VERSION		"3.8.9"

short default_work_type (void);
#define WORK_FACTOR		0
#define WORK_TEST		1
#define WORK_ADVANCEDTEST	2
#define WORK_DBLCHK		3
#define WORK_ECM		4
#define WORK_PMINUS1		5
#define WORK_PFACTOR		6
#define WORK_ADVANCEDFACTOR	7

struct work_unit {		/* One line from the worktodo file */
	int	work_type;	/* Type of work to do */
	char	assignment_uid[33]; /* Primenet assignment ID */
	char	extension[9];	/* Optional save file extension */
	double	k;		/* K in k*b^n+c */
	unsigned long b;	/* B in k*b^n+c */
	unsigned long n;	/* N in k*b^n+c */
	signed long c;		/* C in k*b^n+c */
	unsigned long minimum_fftlen;/* Minimum FFT length to use.  Zero means default fftlen */
	double	sieve_depth;	/* How far it has been trial factored */
	double	factor_to;	/* How far we should trial factor to */
	int	pminus1ed;	/* TRUE if has been P-1 factored */
	double	B1;		/* ECM and P-1 - Stage 1 bound */
	double	B2_start;	/* ECM and P-1 - Stage #2 start */
	double	B2;		/* ECM and P-1 - Stage #2 end */
	int	nth_run;	/* P+1 - 1 for start 2/7, 2 for start 6/5, 3+ for random start  JP */
	unsigned int curves_to_do; /* ECM - curves to try */
	double	curve;		/* ECM - Specific curve to test (debug tool) */
	double	tests_saved;	/* Pfactor - primality tests saved if a factor is found */
	unsigned int prp_base;	/* PRP base to use */	
	int	prp_residue_type; /* PRP residue to output -- see primenet.h */
	int	prp_dblchk;	/* True if this is a doublecheck of a previous PRP */
	int	cert_squarings; /* Number of squarings required for PRP proof certification JP */
	char	*known_factors;	/* ECM, P-1, PRP - list of known factors */
	char	*comment;	/* Comment line in worktodo.ini */
		/* Runtime variables */
	struct work_unit *next; /* Next in doubly-linked list */
	struct work_unit *prev; /* Previous in doubly-linked list */
	int	in_use_count;	/* Count of threads accessing this work unit */
	int	high_memory_usage;/* Set if we are using a lot of memory */
				/* If user changes the available memory */
				/* settings, then we should stop and */
				/* restart our computations */
	char	stage[11];	/* Test stage (e.g. TF,P-1,LL) */
	double	pct_complete;	/* Percent complete (misnomer as value is */
				/* between 0.0 and 1.0) */
	unsigned long fftlen;	/* FFT length in use */
	int	ra_failed;	/* Set when register assignment fails, tells */
				/* us not to try registering it again. */
};

struct work_unit_array {	/* All the lines for one worker thread */
	struct work_unit *first; /* First work unit */
	struct work_unit *last;	/* Last work unit */
};


/* Global variables */

extern char INI_FILE[80];		/* Name of the prime INI file */

extern int ERRCHK;					/* 1 to turn on error checking */
extern unsigned int PRIORITY;		/* Desired priority level */
extern unsigned int CPU_AFFINITY;	/* NT Processor affinity */


extern unsigned long volatile ITER_OUTPUT;/* Iterations between outputs */
extern unsigned long volatile ITER_OUTPUT_RES;/* Iterations between results */
					/* file outputs */
extern unsigned long volatile DISK_WRITE_TIME;
					/* Number of minutes between writing */
					/* intermediate results to disk */
extern int TWO_BACKUP_FILES;		/* TRUE for 2 backup files(qXXXXXXX) */
extern int RUN_ON_BATTERY;		/* Run program even on battery power */
extern int TRAY_ICON;			/* Display tiny tray icon */
extern int HIDE_ICON;			/* Display no icon */
extern unsigned int PRECISION;		/* Number of decimal places to output*/
					/* in percent complete lines */
extern int CUMULATIVE_TIMING;		/* True if outputting cumulative time*/

/* Common routines */


// void getCpuInfo (); 
 
int isPrime (unsigned long p);

void nameIniFiles (int named_ini_files);
void readIniFiles ();

void IniGetString (char *, char *, char *, unsigned int, char *);
long IniGetInt (char *, char *, long);
void IniWriteString (char *, char *, char *);
void IniWriteInt (char *, char *, long);

void IniFileOpen (char *, int);
void processTimedIniFile (char *);
void IniFileClose (char *);
unsigned int IniGetNumLines (char *);
void IniGetLineAsString (char *, unsigned int, char *, unsigned int,
			 char *, unsigned int);
void IniGetLineAsInt (char *, unsigned int, char *, unsigned int, long *);
void IniReplaceLineAsString (char *, unsigned int, char *, char *);
void IniReplaceLineAsInt (char *, unsigned int, char *, long);
void IniInsertLineAsString (char *, unsigned int, char *, char *);
void IniInsertLineAsInt (char *, unsigned int, char *, long);
void IniAppendLineAsString (char *, char *, char *);
void IniAppendLineAsInt (char *, char *, long);
void IniDeleteLine (char *, unsigned int);
void IniDeleteAllLines (char *);

void OutputBoth (char *);
void OutputSomewhere (char *);
void LogMsg (char *);
void ReplaceableLine (int);

unsigned long pick_fft_length (unsigned long);

void tempFileName (char	*, char, giant);
int fileExists (char *);
int readFileHeader (char *, int *, short *, unsigned long *);
int writeResults (char	*);


/* Routines called by common routines */

void OutputStr (char *);
// int isHighResTimerAvailable (void); 
// double getHighResTimer (void); 
// double getHighResTimerFrequency (void); 
unsigned long num_cpus ();
#define stopCheck escapeCheck
int escapeCheck ();
#define	WORKING_ICON	0
#define	IDLE_ICON	1
void ChangeIcon (int);
void BlinkIcon (int);


/* Common routines */

int primeContinue ();

/* Routines used to time code chunks */

extern double __timers[10];		/* 10 timers are available */

/* Routines called by common routines */

void title (char *);
void flashWindowAndBeep ();
void SetPriority ();


#endif
