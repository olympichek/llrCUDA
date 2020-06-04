#ifndef LLRDOTH
#define LLRDOTH

/* Constants */

#define VERSION		"3.8.0"

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
