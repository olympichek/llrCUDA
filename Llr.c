/*----------------------------------------------------------------------
| This file contains routines and global variables that are common for
| all operating systems the program has been ported to.  It is included
| in one of the source code files of each port.  See Llr.h for the
| common #defines and common routine definitions.
+---------------------------------------------------------------------*/
 
#define CHECK_IF_ANY_ERROR(X,J,N,K) \
		checknumber = K;\
\
/* Check for excessive roundoff error  */\
\
		if (MAXERR > maxroundoff) {\
			lasterr_point = J;\
			if (J == last_bit[K] &&\
			    MAXERR == last_maxerr[K] && !abonroundoff && !will_try_larger_fft) {\
				clearline(100);\
				OutputBoth (ERROK);\
				GWERROR = 0;\
				clearline(100);\
				OutputBoth (ERRMSG6);\
				maxerr_recovery_mode[K] = TRUE;\
				sleep5 = FALSE;\
				goto error;\
			} else {\
				char	msg[80];\
				sprintf (msg, ERRMSG1C, MAXERR, maxroundoff);\
				sprintf (buf, ERRMSG0L, J, N, msg);\
				clearline(100);\
				OutputBoth (buf);\
				/*g_roundoff=1;*//*cuda*/\
				will_try_larger_fft = TRUE;/*cuda*/\
				gwypdone();/*cuda*/\
				if (J == last_bit[K])\
					will_try_larger_fft = TRUE;\
				if (will_try_larger_fft) {\
					last_bit[K]  = 0;\
					last_maxerr[K]  = 0.0;\
				}\
				else {\
					last_bit[K] = J;\
					last_maxerr[K] = MAXERR;\
				}\
				sleep5 = FALSE;\
				goto error;\
			}\
		}\
\
		if (ERRCHK) {\
			if (MAXERR < reallyminerr && J > 30)\
				reallyminerr = MAXERR;\
			if (MAXERR > reallymaxerr)\
				reallymaxerr = MAXERR;\
		} 

int it = 0; //cuda ; init to zero, JP 22/06/17

// Some ABC format strings

char ckstring[] = "(2^$a$b)^2-2";	// Carol/Kynea
char cwstring[] = "$a*$b^$a$c";		// Cullen/Woodall
char ffstring[] = "$a*2^$b+1";		// FermFact output
char ffmstring[] = "$a*2^$b-1";		// Lei output
char gmstring[] = "4^$a+1";		// Gaussian Mersenne norms
char gfstring[] = "$a^$b+1";		// Special primality test for generalized Fermat numbers
char spstring[] = "(2^$a+1)/3";		// Special SPRP test for Wagstaff numbers
char repustring[] = "(10^$a-1)/9";	// PRP test for repunits numbers
char grepustring[] = "($a^$b-1)/($a-1)";// PRP test for generalized repunits numbers
char diffnumpstring[] = "$a^$b-$a^$c+%d";// If $b>$c, it is [$a^($b-$c)-1]*$a^$c+%d so, form K*B^N+C
char diffnummstring[] = "$a^$b-$a^$c-%d";// If $b>$c, it is [$a^($b-$c)-1]*$a^$c-%d so, form K*B^N+C
char diffnumstring[] = "$a^$b-$a^$c$d";	// General diffnum format

// Fixed k and c forms for k*b^n+c

char fkpstring[] = "%d*$a^$b+%d";
char fkmstring[] = "%d*$a^$b-%d";

// Fixed b and c forms for k*b^n+c

char fbpstring[]  = "$a*"$LLF"^$b+%d";
char fbmstring[]  = "$a*"$LLF"^$b-%d";

// Fixed n and c forms for k*b^n+c

char fnpstring[] = "$a*$b^%lu+%d";
char fnmstring[] = "$a*$b^%lu-%d";

// Fixed k forms for k*b^n+c

char fkastring[]  = ""$LLF"*$a^$b$c";

// Fixed b forms for k*b^n+c

char fbastring[] = "$a*"$LLF"^$b$c";

// Fixed n forms for k*b^n+c

char fnastring[]  = "$a*$b^%lu$c";

// Fixed c forms for k*b^n+c

char abcpstring[]  = "$a*$b^$c+%d";
char abcmstring[]  = "$a*$b^$c-%d";

// General k*b^n+c format 

char abcastring[] = "$a*$b^$c$d";

// General (k*b^n+c)/d format 

char abcadstring[] = "($a*$b^$c$d)/$e";

// Test the primality of a number given as a string

char numberstring[] = "$a";

// Test if $a is a base $b Wieferich prime

char wftstring[] = "$a$b";

// Search for base $c Wieferich primes in the range $a to $b

char wfsstring[] = "$a$b$c";

/* Process a number from newpgen output file */
/* NEWPGEN output files use the mask as defined below: */

#define MODE_PLUS    0x01	/* k.b^n+1 */
#define MODE_MINUS   0x02	/* k.b^n-1 */
#define MODE_2PLUS   0x04	/* k.b^(n+1)+1 (*) */
#define MODE_2MINUS  0x08	/* k.b^(n+1)-1 (*) */
#define MODE_4PLUS   0x10	/* k.b^(n+2)+1 (*) */
#define MODE_4MINUS  0x20	/* k.b^(n+2)-1 (*) */
#define MODE_PRIMORIAL 0x40	/* PRIMORIAL - can't handle this */
#define MODE_PLUS5  0x80	/* k.b^n+5 */
#define MODE_2MINUS3 0x100	/* 2k.b^n-3 JP 23/08/17 */
#define MODE_AP	    0x200	/* 2^n+2k-1 */
#define MODE_PLUS7  0x800	/* k.b^n+7 */
#define MODE_2PLUS3 0x1000	/* 2k.b^n+3 */
#define MODE_DUAL 0x8000
#define MODE_PLUS_DUAL 0x8001	/* b^n+k */
#define MODE_MINUS_DUAL 0x8002	/* b^n-k */
#define MODE_NOTGENERALISED 0x400


/* Define the world/group/owner read/write attributes for creating files */
/* I've always used 0666 in Unix (everyone gets R/W access), but MSVC 8 */
/* now refuses to work with that setting -- insisting on 0600 instead. */

#ifdef _WIN32
#define	CREATE_FILE_ACCESS	0600
#else
#define	CREATE_FILE_ACCESS	0666
#endif

#define IBSIZE 300

char	greatbuf[10001] = {0};
char	INI_FILE[80] = {0};
char	SVINI_FILE[80] = {0};
char	RESFILE[80] = {0};
char	LOGFILE[80] = {0};
char	EXTENSION[8] = {0};

int fftlen = 0;
int ERRCHK = 0;

unsigned int PRIORITY = 1;
unsigned int CPU_AFFINITY = 99;
unsigned int CPU_TYPE = 0;
unsigned long CPU_SPEED = 25;
int      nbdg, nbdg1, nbdg2;
int      showdigits = FALSE;
int      maxaprcl = 400;
int ZERO_PADDED_FFT;
int GWERROR = 0;
unsigned long volatile ITER_OUTPUT = 0;
unsigned long volatile ITER_OUTPUT_RES = 99999999;
unsigned long volatile DISK_WRITE_TIME = 30;
int	TWO_BACKUP_FILES = 1;
int	RUN_ON_BATTERY = 1;
int	TRAY_ICON = TRUE;
int	HIDE_ICON = FALSE;
unsigned int PRECISION = 2;
// int	CUMULATIVE_TIMING = 0;
int	HIGH_RES_TIMER = 0; 

// double MAXERR = 0.0;
double maxdiffmult = 1.0;

/* PRP and LLR global variables */

#define	sgkbufsize 20000

giant	N = NULL;		/* Number being tested */
giant	NP = NULL;		/* Number being tested */
giant	M = NULL;		/* Gaussian Mersenne modulus = N*NP */
giant	gk = NULL;		/* k multiplier */
giant	gb = NULL;		/* Generalized Fermat base may be a large integer... */

unsigned long Nlen = 0;	/* Bit length of number being LLRed or PRPed */
unsigned long klen = 0;	/* Number of bits of k multiplier */
long OLDFFTLEN = 0; /* previous value of FFTLEN, used by setuponly option */
unsigned long ndiff = 0;/* used for b^n-b^m+c number processing */
unsigned long gformat;	/* used for b^n-b^m+c number processing */
unsigned long globalb = 2;	// base of the candidate in a global
double	globalk = 1.0;   // k value of the candidate in a global
							

/* Other gwypnum globals */

giant testn, testnp;
unsigned long facn = 0, facnp = 0;
int resn = 0, resnp = 0;
char facnstr[80], facnpstr[80];
char sgd[sgkbufsize];
// char sgq[sgkbufsize];
char sgb[sgkbufsize];
char bpfstring[sgkbufsize];

static unsigned long last_bit[10] = {0};
static double last_suminp[10] = {0.0};
static double last_sumout[10] = {0.0};
static double last_maxerr[10] = {0.0};
static double maxroundoff = 0.40;

static unsigned long mask;

extern int zcomplex;
extern int cufftonly;

unsigned int Fermat_only = FALSE;
unsigned int strong = TRUE;
unsigned int quotient = FALSE;
unsigned int vrbareix = FALSE;
unsigned int dualtest = FALSE;
unsigned int setuponly = FALSE;
unsigned int abonillsum = FALSE;
unsigned int abonmismatch = FALSE;
unsigned int aborted = FALSE;
unsigned int testgm = FALSE;
unsigned int testgq = FALSE;
unsigned int testfac = FALSE;
unsigned int nofac = FALSE;
unsigned int abonroundoff = FALSE;
unsigned int will_try_larger_fft = FALSE;
unsigned int checknumber = 0;
unsigned int sleep5 = FALSE;
unsigned int maxerr_recovery_mode [10] = {0};
unsigned int lasterr_point = 0;
unsigned long interimFiles, interimResidues, throttle, facfrom, facto;
unsigned long factored = 0, eliminated = 0;
unsigned long pdivisor = 1000000, pquotient = 1;
unsigned long bpf[30], bpc[30], vpf[30];
// Base prime factors, cofactors, power of p.f.
giantstruct*	gbpf[30] = {NULL}; // Large integer prime factors
giantstruct*	gbpc[30] = {NULL}; // Large integer cofactors
unsigned long nrestarts = 0;
// Nb. of restarts for an N+1 or N-1 prime test

double smargin = 0.0;
double	pcfftlim = 0.5;	

int genProthBase(giant, uint32_t);
long generalLucasBase(giant , uint32_t *, uint32_t *);
unsigned long gcd (
	unsigned long x,
	unsigned long y);

/* Utility output routines */

void LineFeed ()
{
#ifndef WIN32
	OutputStr ((char*)"\r");
#elif defined (_CONSOLE)
	OutputStr ((char*)"\r");
#else
	OutputStr ((char*)"\n");
#endif
}

void OutputNum (
	unsigned long num)
{
	char	buf[20];
	sprintf (buf, "%lu", num);
	OutputStr (buf);
}

/* Sleep five minutes before restarting */

char ERROK[] = "Disregard last error.  Result is reproducible and thus not a hardware problem.\n";
char ERRMSG0L[] = "Iter: %ld/%ld, %s";
char ERRMSG0[] = "Bit: %ld/%ld, %s"; 
char ERRMSG1A[] = "ERROR: ILLEGAL SUMOUT\n";
char ERRMSG1B[] = "ERROR: SUM(INPUTS) != SUM(OUTPUTS), %.16g != %.16g\n";
char ERRMSG1C[] = "ERROR: ROUND OFF (%.10g) > %.10g\n";
char ERRMSG2[] = "Possible hardware failure, consult the readme file.\n";
char ERRMSG3[] = "Continuing from last save file.\n";
char ERRMSG4[] = "Waiting five minutes before restarting.\n";
char ERRMSG5[] = "Fatal Error, Check Number = %d, test of %s aborted\n";
char ERRMSG6[] = "For added safety, redoing iteration using a slower, more reliable method.\n";
char ERRMSG7[] = "Fatal Error, divisibility test of %d^(N-1)-1 aborted\n";
char ERRMSG8[] = "Unrecoverable error, Restarting with next larger FFT length...\n";
char WRITEFILEERR[] = "Error writing intermediate file: %s\n";


void	trace(int n) {			// Debugging tool...
	char buf[100];
	sprintf(buf, "OK until number %d\n", n);
	OutputBoth (buf); 	
}

void	strace(int n) {			// Debugging tool...
	char buf[100];
	sprintf(buf, "OK until number %d\n", n);
	OutputStr (buf); 	
}

void clearline (int size) {
	char buf[256];
	int i;

	for (i=0; i<256; i++)
		buf[i] = '\0';
	for (i=0; i<size; i++)
		buf[i] = ' ';
	buf[size-1] = '\r';
#if !defined(WIN32) || defined(_CONSOLE)
	OutputStr(buf);
#endif
}

int	digitstrcmp (const char *s1, const char *s2) {
	if (strlen (s1) < strlen (s2))
		return (-1);
	else if (strlen (s1) > strlen (s2))
		return (1);
	else 
		return (strcmp (s1, s2));
}

int SleepFive ()
{
	int	i;

	OutputStr (ERRMSG4);
	BlinkIcon (10);			/* Blink icon for 10 seconds */
	Sleep (10000);
	ChangeIcon (IDLE_ICON);		/* Idle icon for rest of 5 minutes */
	for (i = 0; i < 290; i++) {
		Sleep (1000);
		if (escapeCheck ()) return (FALSE);
	}
	ChangeIcon (WORKING_ICON);	/* And back to the working icon */
	return (TRUE);
}

/* Truncate a percentage to the requested number of digits. */
/* Truncating prevents 99.5% from showing up as 100% complete. */

double trunc_percent (
	double	percent)
{
	if (percent > 100.0) percent = 100.0;
	percent -= 0.5 * pow (10.0, - (double) PRECISION);
	if (percent < 0.0) return (0.0);
	return (percent);
}

 
//  Test if a string contains only valid digits. 
 
int isDigitString(char *s) { 
    while (*s) { 
	if (!isdigit(*s)) return (FALSE); 
	s++; 
    } 
    return (TRUE); 
} 
 
void OutputTimeStamp ()
{
	time_t	this_time;
	char	tbuf[40], buf[40];

//	if (TIMESTAMPING) {
		time (&this_time);
		strcpy (tbuf, ctime (&this_time)+4);
		tbuf[12] = 0;
		sprintf (buf, "[%s] ", tbuf);
		OutputStr (buf);
//	}
}
 
/* Determine if a small number is prime */

int isPrime (
	unsigned long p)
{
	unsigned long i;
	if (p < 2)	// 1 is not prime.
		return (FALSE);
	for (i = 2; (i < 65536) && (i*i <= p); i = (i + 1) | 1)
            // Avoid an overflow when computing i*i ...
            if (p % i == 0) return (FALSE);
	return (TRUE);
}

/* Determine the names of the INI files */

void nameIniFiles (
	int	named_ini_files)
{
	char	buf[120];

	if (named_ini_files < 0) {
		strcpy (INI_FILE, "llr.ini");
		strcpy (RESFILE, "lresults.txt");
		strcpy (LOGFILE, "lprime.log");
		strcpy (EXTENSION, "");
	} else {
		sprintf (INI_FILE, "llr%04d.ini", named_ini_files);
		sprintf (RESFILE, "lresu%04d.txt", named_ini_files);
		sprintf (LOGFILE, "lprim%04d.log", named_ini_files);
		sprintf (EXTENSION, ".%03d", named_ini_files);
	}

/* Let the user rename these files and pick a different working directory */

	IniGetString (INI_FILE, (char*)"WorkingDir", buf, sizeof(buf), NULL);
	IniGetString (INI_FILE, (char*)"results.txt", RESFILE, 80, RESFILE);
	IniGetString (INI_FILE, (char*)"prime.log", LOGFILE, 80, LOGFILE);
	IniGetString (INI_FILE, (char*)"prime.ini", INI_FILE, 80, INI_FILE);
	if (buf[0]) {
		_chdir (buf);
		IniFileOpen (INI_FILE, 0);
	}
}

/* Read the INI files */

void readIniFiles () 
{ 
	int	temp; 
 
//	getCpuInfo (); 
 
	PRECISION = (unsigned int) IniGetInt (INI_FILE, (char*)"PercentPrecision", 2); 
	if (PRECISION > 6) PRECISION = 6; 
 
	ITER_OUTPUT = IniGetInt (INI_FILE, (char*)"OutputIterations", 10000); 
	if (ITER_OUTPUT <= 0) ITER_OUTPUT = 1; 
	ITER_OUTPUT_RES = IniGetInt (INI_FILE, (char*)"ResultsFileIterations", 
				     99999999); 
	if (ITER_OUTPUT_RES < 1000) ITER_OUTPUT_RES = 1000; 
	DISK_WRITE_TIME = IniGetInt (INI_FILE, (char*)"DiskWriteTime", 30); 
	TWO_BACKUP_FILES = (int) IniGetInt (INI_FILE, (char*)"TwoBackupFiles", 1); 
	RUN_ON_BATTERY = (int) IniGetInt (INI_FILE, (char*)"RunOnBattery", 1); 
 
	temp = (int) IniGetInt (INI_FILE, (char*)"ErrorCheck", 0); 
	ERRCHK = (temp != 0); 
	PRIORITY = (unsigned int) IniGetInt (INI_FILE, (char*)"Priority", 1); 
	CPU_AFFINITY = (unsigned int) IniGetInt (INI_FILE, (char*)"Affinity", 99); 
	HIDE_ICON = (int) IniGetInt (INI_FILE, (char*)"HideIcon", 0); 
	TRAY_ICON = (int) IniGetInt (INI_FILE, (char*)"TrayIcon", 1); 
 
/* Guess the CPU type if it isn't known.  Otherwise, validate it. */ 
 
//	getCpuInfo (); 
 
/* Other oddball options */ 
 
	CUMULATIVE_TIMING = IniGetInt (INI_FILE, (char*)"CumulativeTiming", 0); 
//	HIGH_RES_TIMER = isHighResTimerAvailable (); 
} 
 
/*----------------------------------------------------------------------
| Portable routines to read and write ini files!  NOTE:  These only
| work if you open no more than 5 ini files.  Also you must not
| change the working directory at any time during program execution.
+---------------------------------------------------------------------*/

struct IniLine {
	char	*keyword;
	char	*value;
	int	active;
};
struct IniCache {
	char	*filename;
	int	immediate_writes;
	int	dirty;
	unsigned int num_lines;
	unsigned int array_size;
	struct IniLine **lines;
};

void growIniLineArray (
	struct IniCache *p)
{
	struct IniLine **newlines;

	if (p->num_lines != p->array_size) return;

	newlines = (struct IniLine **)
		malloc ((p->num_lines + 100) * sizeof (struct IniLine **));
	if (p->num_lines) {
		memcpy (newlines, p->lines, p->num_lines * sizeof (struct IniLine *));
		gwypfree (p->lines);
	}
	p->lines = newlines;
	p->array_size = p->num_lines + 100;
}

struct IniCache *openIniFile (
	char	*filename,
	int	forced_read)
{
static	struct IniCache *cache[10] = {0};
	struct IniCache *p;
	FILE	*fd;
	unsigned int i;
	char	line[80];
	char	*val;

/* See if file is cached */

	for (i = 0; i < 10; i++) {
		p = cache[i];
		if (p == NULL) {
			p = (struct IniCache *) malloc (sizeof (struct IniCache));
			p->filename = (char *) malloc (strlen (filename) + 1);
			strcpy (p->filename, filename);
			p->immediate_writes = 1;
			p->dirty = 0;
			p->num_lines = 0;
			p->array_size = 0;
			p->lines = NULL;
			forced_read = 1;
			cache[i] = p;
			break;
		}
		if (strcmp (filename, p->filename) == 0)
			break;
	}

/* Skip reading the ini file if appropriate */

	if (!forced_read) return (p);
	if (p->dirty) return (p);

/* Free the data if we've already read some in */

	for (i = 0; i < p->num_lines; i++) {
		gwypfree (p->lines[i]->keyword);
		gwypfree (p->lines[i]->value);
		gwypfree (p->lines[i]);
	}
	p->num_lines = 0;

/* Read the IniFile */
	
	fd = fopen (filename, "r");
	if (fd == NULL) return (p);

	while (fgets (line, 80, fd)) {
		if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = 0;
		if (line[0] == 0) continue;
		if (line[strlen(line)-1] == '\r') line[strlen(line)-1] = 0;
		if (line[0] == 0) continue;

		val = strchr (line, '=');
		if (val == NULL) {
			char	buf[130];
			sprintf (buf, "Illegal line in INI file: %s\n", line);
			OutputSomewhere (buf);
			continue;
		}
		*val++ = 0;

		growIniLineArray (p);
		
/* Allocate and fill in a new line structure */

		i = p->num_lines++;
		p->lines[i] = (struct IniLine *) malloc (sizeof (struct IniLine));
		p->lines[i]->keyword = (char *) malloc (strlen (line) + 1);
		p->lines[i]->value = (char *) malloc (strlen (val) + 1);
		p->lines[i]->active = TRUE;
		strcpy (p->lines[i]->keyword, line);
		strcpy (p->lines[i]->value, val);
	}
	fclose (fd);

	return (p);
}

void writeIniFile (
	struct IniCache *p)
{
	int	fd;
	unsigned int j;
	char	buf[100];

/* Delay writing the file unless this INI file is written */
/* to immediately */

	if (!p->immediate_writes) {
		p->dirty = 1;
		return;
	}

/* Create and write out the INI file */

	fd = _open (p->filename, _O_CREAT | _O_TRUNC | _O_WRONLY | _O_TEXT, 0666);
	if (fd < 0) return;
	for (j = 0; j < p->num_lines; j++) {
		strcpy (buf, p->lines[j]->keyword);
		strcat (buf, "=");
		strcat (buf, p->lines[j]->value);
		strcat (buf, "\n");
		_write (fd, buf, strlen (buf));
	}
	p->dirty = 0;
	_close (fd);
}

void save_IniFile (char *filename, char *savedfilename) {
	struct IniCache *p;	
	p = openIniFile (filename, 1);		// open and read the source IniFile.
	p->filename = savedfilename;		// put the target filename in the structure.
	writeIniFile (p);					// Write the target.
	p->filename = filename;				// Restore the structure in memory.
}

void truncated_strcpy (
	char	*buf,
	unsigned int bufsize,
	char	*val)
{
	if (strlen (val) >= bufsize) {
		memcpy (buf, val, bufsize-1);
		buf[bufsize-1] = 0;
	} else {
		strcpy (buf, val);
	}
}

void IniGetString (
	char	*filename,
	char	*keyword,
	char	*val,
	unsigned int val_bufsize,
	char	*default_val)
{
	struct IniCache *p;
	unsigned int i;

/* Open ini file */

	p = openIniFile (filename, 1);

/* Look for the keyword */

	for (i = 0; ; i++) {
		if (i == p->num_lines) {
			if (default_val == NULL) {
				val[0] = 0;
			} else {
				truncated_strcpy (val, val_bufsize, default_val);
			}
			return;
		}
		if (p->lines[i]->active &&
		    stricmp (keyword, p->lines[i]->keyword) == 0) break;
	}

/* Copy info from the line structure to the user buffers */

	truncated_strcpy (val, val_bufsize, p->lines[i]->value);
}

long IniGetInt (
	char	*filename,
	char	*keyword,
	long	default_val)
{
	char	buf[20], defval[20];
	sprintf (defval, "%ld", default_val);
	IniGetString (filename, keyword, buf, 20, defval);
	return (atol (buf));
}

void IniWriteString (
	char	*filename,
	char	*keyword,
	char	*val)
{
	struct IniCache *p;
	unsigned int i, j;

/* Open ini file */
	p = openIniFile (filename, 1);
/* Look for the keyword */

	for (i = 0; ; i++) {
		if (i == p->num_lines ||
		    stricmp (p->lines[i]->keyword, "Time") == 0) {

/* Ignore request if we are deleting line */

			if (val == NULL) return;

/* Make sure the line array has room for the new line */

			growIniLineArray (p);

/* Shuffle entries down to make room for this entry */

			for (j = p->num_lines; j > i; j--)
				p->lines[j] = p->lines[j-1];

/* Allocate and fill in a new line structure */

			p->lines[i] = (struct IniLine *) malloc (sizeof (struct IniLine));
			p->lines[i]->keyword = (char *) malloc (strlen (keyword) + 1);
			strcpy (p->lines[i]->keyword, keyword);
			p->lines[i]->value = NULL;
			p->num_lines++;
			break;
		}
		if (p->lines[i]->active &&
		    stricmp (keyword, p->lines[i]->keyword) == 0) {
			if (val != NULL && strcmp (val, p->lines[i]->value) == 0) return;
			break;
		}
	}
/* Delete the line if requested */

	if (val == NULL) {
		IniDeleteLine (filename, i+1);
		return;
	}
/* Replace the value associated with the keyword */

	gwypfree (p->lines[i]->value);
	p->lines[i]->value = (char *) malloc (strlen (val) + 1);
	strcpy (p->lines[i]->value, val);

/* Write the INI file back to disk */
	writeIniFile (p);
}

void IniWriteInt (
	char	*filename,
	char	*keyword,
	long	val)
{
	char	buf[40];
	sprintf (buf, "%ld", val);
	IniWriteString (filename, keyword, buf);
}

void IniFileOpen (
	char	*filename,
	int	immediate_writes)
{
	struct IniCache *p;
	p = openIniFile (filename, 1);
	p->immediate_writes = immediate_writes;
}

void IniFileClose (
	char	*filename)
{
	struct IniCache *p;
	p = openIniFile (filename, 0);
	if (p->dirty) {
		p->immediate_writes = 1;
		writeIniFile (p);
		p->immediate_writes = 0;
	}
}

int IniFileWritable (
	char	*filename)
{
	struct IniCache *p;
	int	fd;
	unsigned int j;
	char	buf[100];

/* Create and write out the INI file */

	p = openIniFile (filename, 0);
	fd = _open (p->filename, _O_CREAT | _O_TRUNC | _O_WRONLY | _O_TEXT, 0666);
	if (fd < 0) return (FALSE);
	for (j = 0; j < p->num_lines; j++) {
		strcpy (buf, p->lines[j]->keyword);
		strcat (buf, "=");
		strcat (buf, p->lines[j]->value);
		strcat (buf, "\n");
		if (_write (fd, buf, strlen (buf)) != (int) strlen (buf)) {
			_close (fd);
			return (FALSE);
		}
	}
	if (p->num_lines == 0) {
		if (_write (fd, "DummyLine=XXX\n", 14) != 14) {
			_close (fd);
			return (FALSE);
		}
		p->dirty = 1;
	}
	_close (fd);
	return (TRUE);
}

unsigned int IniGetNumLines (
	char	*filename)
{
	struct IniCache *p;
	p = openIniFile (filename, 0);
	return (p->num_lines);
}

void IniGetLineAsString (
	char	*filename,
	unsigned int line,
	char	*keyword,
	unsigned int keyword_bufsize,
	char	*val,
	unsigned int val_bufsize)
{
	struct IniCache *p;

/* Open ini file */

	p = openIniFile (filename, 0);

/* Copy info from the line structure to the user buffers */

	truncated_strcpy (keyword, keyword_bufsize, p->lines[line-1]->keyword);
	truncated_strcpy (val, val_bufsize, p->lines[line-1]->value);
}

void IniGetLineAsInt (
	char	*filename,
	unsigned int line,
	char	*keyword,
	unsigned int keyword_bufsize,
	long	*val)
{
	char	buf[20];
	IniGetLineAsString (filename, line, keyword, keyword_bufsize, buf, 20);
	*val = atol (buf);
}

void IniReplaceLineAsString (
	char	*filename,
	unsigned int line,
	char	*keyword,
	char	*val)
{
	IniDeleteLine (filename, line);
	IniInsertLineAsString (filename, line, keyword, val);
}

void IniReplaceLineAsInt (
	char	*filename,
	unsigned int line,
	char	*keyword,
	long	val)
{
	char	buf[20];
	sprintf (buf, "%ld", val);
	IniReplaceLineAsString (filename, line, keyword, buf);
}

void IniInsertLineAsString (
	char	*filename,
	unsigned int line,
	char	*keyword,
	char	*val)
{
	struct IniCache *p;
	unsigned int i;

/* Open ini file, do not reread it as that could change the line numbers! */

	p = openIniFile (filename, 0);

/* Adjust line number if it doesn't make sense */

	if (line == 0) line = 1;
	if (line > p->num_lines+1) line = p->num_lines+1;

/* Make sure the line array has room for the new line */

	growIniLineArray (p);

/* Shuffle lines down in the array to make room for the new line */

	for (i = p->num_lines; i >= line; i--) p->lines[i] = p->lines[i-1];
	p->num_lines++;

/* Allocate and fill in a new line structure */

	p->lines[line-1] = (struct IniLine *) malloc (sizeof (struct IniLine));
	p->lines[line-1]->keyword = (char *) malloc (strlen (keyword) + 1);
	p->lines[line-1]->value = (char *) malloc (strlen (val) + 1);
	p->lines[line-1]->active = TRUE;
	strcpy (p->lines[line-1]->keyword, keyword);
	strcpy (p->lines[line-1]->value, val);

/* Write the INI file back to disk */

	writeIniFile (p);
}

void IniInsertLineAsInt (
	char	*filename,
	unsigned int line,
	char	*keyword,
	long	val)
{
	char	buf[20];
	sprintf (buf, "%ld", val);
	IniInsertLineAsString (filename, line, keyword, buf);
}

void IniAppendLineAsString (
	char	*filename,
	char	*keyword,
	char	*val)
{
	struct IniCache *p;
	p = openIniFile (filename, 0);
	IniInsertLineAsString (filename, p->num_lines+1, keyword, val);
}

void IniAppendLineAsInt (
	char	*filename,
	char	*keyword,
	long	val)
{
	char	buf[20];
	sprintf (buf, "%ld", val);
	IniAppendLineAsString (filename, keyword, buf);
}

void IniDeleteLine (
	char	*filename,
	unsigned int line)
{
	struct IniCache *p;
	unsigned int i;

/* Open ini file, do not reread it as that could change the line numbers! */

	p = openIniFile (filename, 0);
	if (line == 0 || line > p->num_lines) return;

/* Free the data associated with the given line */

	gwypfree (p->lines[line-1]->keyword);
	gwypfree (p->lines[line-1]->value);
	gwypfree (p->lines[line-1]);

/* Delete the line from the lines array */

	for (i = line; i < p->num_lines; i++) p->lines[i-1] = p->lines[i];
	p->num_lines--;

/* Write the INI file back to disk */

	writeIniFile (p);
}

void IniDeleteAllLines (
	char	*filename)
{
	struct IniCache *p;
	unsigned int i;

/* Open ini file! */

	p = openIniFile (filename, 0);

/* Free the data associated with the given line */

	for (i = 0; i < p->num_lines; i++) {
		gwypfree (p->lines[i]->keyword);
		gwypfree (p->lines[i]->value);
		gwypfree (p->lines[i]);
	}
	p->num_lines = 0;

/* Write the INI file back to disk */

	writeIniFile (p);
}

/* Output string to screen or results file */

void OutputSomewhere (
	char	*buf)
{
	if (NO_GUI) writeResults (buf);
	else OutputStr (buf);
}

/* Output string to both the screen and results file */

void OutputBoth (
	char	*buf)
{
	OutputStr (buf);
	writeResults (buf);
}

void OutputAuxTimes (void)
{
	char Auxtimebuf [256];
	if (debug) {
		sprintf (Auxtimebuf, "Squaring and/or Mult. time : "); 
		gwypwrite_timer (Auxtimebuf+strlen(Auxtimebuf), 5, TIMER_CLR | TIMER_NL);
		OutputBoth (Auxtimebuf); 
		sprintf (Auxtimebuf, "Normalization time : "); 
		gwypwrite_timer (Auxtimebuf+strlen(Auxtimebuf), 2, TIMER_CLR | TIMER_NL);
		OutputBoth (Auxtimebuf); 
		if (zp) {
			sprintf (Auxtimebuf, "Modular reduction time : "); 
			gwypwrite_timer (Auxtimebuf+strlen(Auxtimebuf), 3, TIMER_CLR | TIMER_NL); 
			OutputBoth (Auxtimebuf); 
		}
		if (generic) {
			sprintf (Auxtimebuf, "Generic copyzero + setzero time : "); 
			gwypwrite_timer (Auxtimebuf+strlen(Auxtimebuf), 4, TIMER_CLR | TIMER_NL); 
			OutputBoth (Auxtimebuf); 
		}
	}
}

/* Output message to screen and prime.log file */

void LogMsg (
	char	*str)
{
	int	fd;
	unsigned long filelen;
static	time_t	last_time = 0;
	time_t	this_time;

/* Output it to the screen */

	OutputStr (str);

/* Open the log file and position to the end */

	fd = _open (LOGFILE, _O_TEXT | _O_RDWR | _O_CREAT, 0666);
	if (fd < 0) {
		OutputStr ((char*)"Unable to open log file.\n");
		return;
	}
	filelen = _lseek (fd, 0L, SEEK_END);

/* Output to the log file only if it hasn't grown too big */

	if (filelen < 250000) {

/* If it has been at least 5 minutes since the last time stamp */
/* was output, then output a new timestamp */

		time (&this_time);
		if (this_time - last_time > 300) {
			char	buf[48];
			last_time = this_time;
			buf[0] = '[';
			strcpy (buf+1, ctime (&this_time));
			sprintf (buf+25, " - ver %s]\n", VERSION);
			_write (fd, buf, strlen (buf));
		}

/* Output the message */

		_write (fd, str, strlen (str));
	}

/* Display message about full log file */
	
	else {
		char	*fullmsg = (char*)"Prime.log file full.  Please delete it.\n";
		OutputStr (fullmsg);
		if (filelen < 251000)
			_write (fd, fullmsg, strlen (fullmsg));
	}
	_close (fd);
}

int gmodi (uint32_t, giant);


/* Generate temporary file name */
 
void tempFileName ( 
	char	*buf, char c, giant NN) 
{ 
	int remainder;
 
	remainder = gmodi(19999981, NN);
	sprintf (buf, "%1c%07i", c, remainder % 10000000); 
} 

/* See if the given file exists */

int fileExists (
	char	*filename)
{
	int	fd;
	fd = _open (filename, _O_RDONLY | _O_BINARY);
	if (fd < 0) return (0);
	_close (fd);
	return (1);
}

/* Open the results file and write a line to the end of it. */

int writeResults (
	char	*msg)
{
static	time_t	last_time = 0;
	time_t	this_time;
	int	fd;

/* Open file, position to end */

	fd = _open (RESFILE, _O_TEXT | _O_RDWR | _O_CREAT | _O_APPEND, 0666);
	if (fd < 0) {
		LogMsg ((char*)"Error opening the results file ; see result below :\n");
		LogMsg (msg);			// Log the unwrited message (Darren Bedwell'request)
		return (FALSE);
	}

/* If it has been at least 5 minutes since the last time stamp */
/* was output, then output a new timestamp */

	time (&this_time);
	if (this_time - last_time > 300) {
		char	buf[32];
		last_time = this_time;
		buf[0] = '[';
		strcpy (buf+1, ctime (&this_time));
		buf[25] = ']';
		buf[26] = '\n';
		if (verbose)
			_write (fd, buf, 27);
	}

/* Output the message */

	if (_write (fd, msg, strlen (msg)) < 0) goto fail;
	_close (fd);
	return (TRUE);

/* On a write error, close file and return error flag */

fail:	_close (fd);
	return (FALSE);
}


/* Read and write intermediate results to a file */
/*
int read_gwypnum (
	int	fd,
	gwypnum	g,
	long	*sum)
{
	giant	tmp;
	long	i, len, bytes;

	tmp = newgiant(FFTLEN*sizeof(double)/sizeof(short) + 16);
	if (_read (fd, &len, sizeof (long)) != sizeof (long)) return (FALSE);
	bytes = len * sizeof (short);
	if (_read (fd, tmp->n, bytes) != bytes) return (FALSE);
	tmp->sign = len;
	*sum += len;
	for (i = 0; i < len; i++) *sum += tmp->n[i];
	gianttogwyp (tmp, g);
	gwypfree (tmp);
	return (TRUE);
}
*/
int read_gwypnum (
        int     fd,
        gwypnum  g,
        long    *sum)
{
 //       giant   tmp;
        long    len, bytes;

        if (_read (fd, &len, sizeof (long)) != sizeof (long)) return (FALSE);
        bytes = len * sizeof (double);
        if (_read (fd, g, bytes) != bytes) return (FALSE);
 	*sum = 0; 
        return (TRUE);
}

/*
int write_gwypnum (
        int     fd,
        gwypnum  g,
        long    *sum)
{
        giant   tmp;
        long    i, len, bytes;

        tmp = newgiant(FFTLEN*sizeof(double)/sizeof(short) + 16);
        gwyptogiant (g, tmp);
        len = tmp->sign;
        if (_write (fd, &len, sizeof (long)) != sizeof (long)) return (FALSE);
        bytes = len * sizeof (short);
        if (_write (fd, tmp->n, bytes) != bytes) return (FALSE);
        *sum += len;
        for (i = 0; i < len; i++) *sum += tmp->n[i];
        gwypfree (tmp);
        return (TRUE);
}
*/
int write_gwypnum (
	int	fd,
	gwypnum	g,
	long	*sum)
{
//	giant	tmp;
	long	len, bytes;

	len = FFTLEN;
	if (_write (fd, &len, sizeof (long)) != sizeof (long)) return (FALSE);
	bytes = len * sizeof (double);
	if (_write (fd, g, bytes) != bytes) return (FALSE);
	*sum = 0;
	return (TRUE);
}

int read_long (
	int	fd,
	unsigned long *val,
	long	*sum)
{
	if (_read (fd, val, sizeof (long)) != sizeof (long)) return (FALSE);
	*sum += *val;
	return (TRUE);
}

int write_long (
	int	fd,
	unsigned long val,
	long	*sum)
{
	if (_write (fd, &val, sizeof (long)) != sizeof (long)) return (FALSE);
	*sum += val;
	return (TRUE);
}

int read_double (
	int	fd,
	double *val,
	long	*sum)
{
	if (_read (fd, val, sizeof (double)) != sizeof (double)) return (FALSE);
	*sum += (long)floor(*val);
	return (TRUE);
}

int write_double (
	int	fd,
	double val,
	long	*sum)
{
	if (_write (fd, &val, sizeof (double)) != sizeof (double)) return (FALSE);
	*sum += (long)floor(val);
	return (TRUE);
}

int writeToFile (
	char	*filename,
	unsigned long j,
	gwypnum	x,
	gwypnum	y)
{
	char	newfilename[16], errmsg[100];
	int	fd;
	unsigned long magicnum, version;
	long	sum = 0, i;


	if (debug) {
		sprintf (errmsg, "Writing %s intermediate file, j = %lu\n",filename, j);
		OutputBoth (errmsg);
	}


/* If we are allowed to create multiple intermediate files, then */
/* write to a file called yNNNNNNN. */

	strcpy (newfilename, filename);
	if (TWO_BACKUP_FILES) newfilename[0] = 'y';

/* Create the intermediate file */

	fd = _open (newfilename, _O_BINARY|_O_WRONLY|_O_TRUNC|_O_CREAT, 0666);
	if (fd < 0) return (FALSE);
/* Write the file header. */

	magicnum = 0x9f2b3cd4;
	if (_write (fd, &magicnum, sizeof (long)) != sizeof (long))
		goto writeerr;

	version = 1;
	if (_write (fd, &version, sizeof (long)) != sizeof (long))
		goto writeerr;

	if (_write (fd, &FFTLEN, sizeof (int)) != sizeof (int)) //cuda
		goto writeerr; //cuda
/* Write the file data */

	if (!write_long (fd, j, &sum)) goto writeerr;

/* Write the data values */

	if (!write_gwypnum (fd, x, &sum)) goto writeerr;
	if (y != NULL && !write_gwypnum (fd, y, &sum)) goto writeerr; 

/* Write the checksum */

	if (_write (fd, &sum, sizeof (long)) != sizeof (long)) goto writeerr;

/* Save the timers */

	for (i=0; i<NBTIMERS; i++) {
		if (gwyptimers[i+NBTIMERS] != 0.0) {// if the timer was running
			gwypend_timer (i);			// update and save it
			if (! write_double (fd, gwyptimers[i], &sum)) {
				gwypstart_timer (i);	// and then, restart it, even if write is in error!
				goto writeerr;
			}
			if (! write_double (fd, gwyptimers[i+NBTIMERS], &sum)) { // save the timer status
				gwypstart_timer (i);	// and then, restart it, even if write is in error!
				goto writeerr;
			}
			gwypstart_timer (i);	// and then, restart it!
		}
		else {
			if (! write_double (fd, gwyptimers[i], &sum)) goto writeerr;	// save the timer
			if (! write_double (fd, gwyptimers[i+NBTIMERS], &sum)) goto writeerr;	// save its status
		}
	}

	_commit (fd);
	_close (fd);

/* Now rename the intermediate files */

	if (TWO_BACKUP_FILES) {
		_unlink (filename);
		rename (newfilename, filename);
	}
	return (TRUE);

/* An error occured.  Close and delete the current file. */

writeerr:
	_close (fd);
	_unlink (newfilename);
	return (FALSE);
}
int readFFTLENFromFile (
        char    *filename,
        unsigned long *j,
        gwypnum  x,
        gwypnum  y)
{
        int     fd;
        unsigned long magicnum, version;
 //       long    sum = 0, i;
        char    errmsg[100];

/* Open the intermediate file */

        fd = _open (filename, _O_BINARY | _O_RDONLY);
        if (fd < 0) goto error;

/* Read the file header */

        if (_read (fd, &magicnum, sizeof (long)) != sizeof (long))
                goto readerr;
        if (magicnum != 0x9f2b3cd4) goto readerr;

        if (_read (fd, &version, sizeof (long)) != sizeof (long)) goto readerr;
        if (version != 1 && version != 2) goto readerr;
        if (_read (fd, &s_FFTLEN, sizeof (int)) != sizeof (int)) goto readerr; //cuda
	_close (fd);
	return (TRUE);
readerr:
	sprintf (errmsg,"Error reading %s intermediate file, j = %lu\n", filename, *j);
	OutputStr (errmsg);
	_close (fd);
error:
	_unlink (filename);
	return (FALSE);
}

int readFromFile (
	char	*filename,
	unsigned long *j,
	gwypnum	x,
	gwypnum	y)
{
	int	fd;
	unsigned long magicnum, version;
	long	sum = 0, i;
	char	errmsg[100];

/* Open the intermediate file */

	fd = _open (filename, _O_BINARY | _O_RDONLY);
	if (fd < 0) goto error;

/* Read the file header */

	if (_read (fd, &magicnum, sizeof (long)) != sizeof (long))
		goto readerr;
	if (magicnum != 0x9f2b3cd4) goto readerr;

	if (_read (fd, &version, sizeof (long)) != sizeof (long)) goto readerr;
	if (version != 1 && version != 2) goto readerr;
	if (_read (fd, &s_FFTLEN, sizeof (int)) != sizeof (int)) goto readerr; //cuda

/* Read the file data */

	if (!read_long (fd, j, &sum)) goto readerr;

/* Read the values */

	if (!read_gwypnum (fd, x, &sum)) goto readerr;
	if (y != NULL && !read_gwypnum (fd, y, &sum)) goto readerr; 

/* Read and compare the checksum */

	if (_read (fd, &i, sizeof (long)) != sizeof (long)) goto readerr;
	if (i != sum) goto readerr;

/* Read the timers and their status */

	for (i=0; i<NBTIMERS; i++) {
		if (!read_double (fd, &gwyptimers[i], &sum)) goto readerr;
		if (!read_double (fd, &gwyptimers[i+NBTIMERS], &sum)) goto readerr;
	}

	_close (fd);

	if (debug) {
		sprintf (errmsg, "Reading %s intermediate file, j = %lu\n", filename, *j);
		OutputBoth (errmsg);
	}

	return (TRUE);

/* An error occured.  Delete the current intermediate file. */
/* Set stage to -1 to indicate an error. */

readerr:
	sprintf (errmsg,"Error reading %s intermediate file, j = %lu\n", filename, *j);
	OutputStr (errmsg);
	_close (fd);
error:
	_unlink (filename);
	return (FALSE);
}

int writeToFileB (
	char	*filename,
	unsigned long j,
	unsigned long B,
	unsigned long nr,
	unsigned long *bpf,
	gwypnum	x,
	gwypnum	y)
{
	char	newfilename[16],errmsg[100];
	int	fd;
	unsigned long magicnum, version;
	long	sum = 0, i;


	if (debug) {
		sprintf (errmsg, "Writing %s intermediate file, j = %lu\n",filename, j);
		OutputBoth (errmsg);
	}

/* If we are allowed to create multiple intermediate files, then */
/* write to a file called yNNNNNNN. */

	strcpy (newfilename, filename);
	if (TWO_BACKUP_FILES) newfilename[0] = 'y';

/* Create the intermediate file */

	fd = _open (newfilename, _O_BINARY|_O_WRONLY|_O_TRUNC|_O_CREAT, 0666);
	if (fd < 0) return (FALSE);

/* Write the file header. */

	magicnum = 0x9f2b3cd4;
	if (_write (fd, &magicnum, sizeof (long)) != sizeof (long))
		goto writeerr;
	version = 1;
	if (_write (fd, &version, sizeof (long)) != sizeof (long))
		goto writeerr;

/* Write the file data */

	if (! write_long (fd, j, &sum)) goto writeerr;
	if (! write_long (fd, B, &sum)) goto writeerr;
	if (! write_long (fd, nr, &sum)) goto writeerr;
	for (i=0; i<10; i++) {
		if (! write_long (fd, bpf[i], &sum)) goto writeerr;
	}

/* Write the data values */

	if (! write_gwypnum (fd, x, &sum)) goto writeerr;
	if (y != NULL && ! write_gwypnum (fd, y, &sum)) goto writeerr; 

/* Write the checksum */

	if (_write (fd, &sum, sizeof (long)) != sizeof (long)) goto writeerr;

/* Save the timers */

	for (i=0; i<NBTIMERS; i++) {
		if (gwyptimers[i+NBTIMERS] != 0.0) {// if the timer was running
			gwypend_timer (i);			// update and save it
			if (! write_double (fd, gwyptimers[i], &sum)) {
				gwypstart_timer (i);	// and then, restart it, even if write is in error!
				goto writeerr;
			}
			if (! write_double (fd, gwyptimers[i+NBTIMERS], &sum)) { // save the timer status
				gwypstart_timer (i);	// and then, restart it, even if write is in error!
				goto writeerr;
			}
			gwypstart_timer (i);	// and then, restart it!
		}
		else {
			if (! write_double (fd, gwyptimers[i], &sum)) goto writeerr;	// save the timer
			if (! write_double (fd, gwyptimers[i+NBTIMERS], &sum)) goto writeerr;	// save its status
		}
	}

	_commit (fd);
	_close (fd);

/* Now rename the intermediate files */

	if (TWO_BACKUP_FILES) {
		_unlink (filename);
		rename (newfilename, filename);
	}
	return (TRUE);

/* An error occured.  Close and delete the current file. */

writeerr:
	_close (fd);
	_unlink (newfilename);
	return (FALSE);
}

int readFromFileB (
	char	*filename,
	unsigned long *j,
	unsigned long *B,
	unsigned long *nr,
	unsigned long *bpf,
	gwypnum	x,
	gwypnum	y)
{
	int	fd;
	unsigned long magicnum, version;
	long	sum = 0, i;
	char	errmsg[100];

/* Open the intermediate file */

	fd = _open (filename, _O_BINARY | _O_RDONLY);
	if (fd < 0) goto error;

/* Read the file header */

	if (_read (fd, &magicnum, sizeof (long)) != sizeof (long))
		goto readerr;
	if (magicnum != 0x9f2b3cd4) goto readerr;

	if (_read (fd, &version, sizeof (long)) != sizeof (long)) goto readerr;
	if (version != 1 && version != 2) goto readerr;

/* Read the file data */

	if (! read_long (fd, j, &sum)) goto readerr;
	if (! read_long (fd, B, &sum)) goto readerr;
	if (! read_long (fd, nr, &sum)) goto readerr;
	for (i=0; i<10; i++) {
		if (! read_long (fd, &bpf[i], &sum)) goto readerr;
	}

/* Read the values */

	if (! read_gwypnum (fd, x, &sum)) goto readerr;
	if (y != NULL && ! read_gwypnum (fd, y, &sum)) goto readerr; 

/* Read and compare the checksum */

	if (_read (fd, &i, sizeof (long)) != sizeof (long)) goto readerr;
	if (i != sum) goto readerr;

/* Read the timers and their status */

	for (i=0; i<NBTIMERS; i++) {
		if (! read_double (fd, &gwyptimers[i], &sum)) goto readerr;
		if (! read_double (fd, &gwyptimers[i+NBTIMERS], &sum)) goto readerr;
	}

	_close (fd);

	if (debug) {
		sprintf (errmsg, "Reading %s intermediate file, j = %lu\n", filename, *j);
		OutputBoth (errmsg);
	}

	return (TRUE);

/* An error occured.  Delete the current intermediate file. */
/* Set stage to -1 to indicate an error. */

readerr:
	sprintf (errmsg,"Error reading %s intermediate file j = %lu\n", filename, *j);
	OutputStr (errmsg);
	_close (fd);
error:
	_unlink (filename);
	return (FALSE);
}

int gmwriteToFile (
	char	*filename,
	unsigned long j,
	unsigned long ubx,
	unsigned long uby,
	gwypnum	x,
	gwypnum	y)
{
	char	newfilename[16], errmsg[100];
	int	fd;
	unsigned long magicnum, version;
	long	sum = 0, i;


	if (debug) {
		sprintf (errmsg, "Writing %s intermediate file, j = %lu\n",filename, j);
		OutputBoth (errmsg);
	}

/* If we are allowed to create multiple intermediate files, then */
/* write to a file called yNNNNNNN. */

	strcpy (newfilename, filename);
	if (TWO_BACKUP_FILES) newfilename[0] = 'y';

/* Create the intermediate file */

	fd = _open (newfilename, _O_BINARY|_O_WRONLY|_O_TRUNC|_O_CREAT, 0666);
	if (fd < 0) return (FALSE);

/* Write the file header. */

	magicnum = 0x9f2b3cd4;
	if (_write (fd, &magicnum, sizeof (long)) != sizeof (long))
		goto writeerr;
	version = 1;
	if (_write (fd, &version, sizeof (long)) != sizeof (long))
		goto writeerr;

/* Write the file data */

	if (! write_long (fd, j, &sum)) goto writeerr;
	if (! write_long (fd, ubx, &sum)) goto writeerr;
	if (! write_long (fd, uby, &sum)) goto writeerr;

/* Write the data values */

	if (! write_gwypnum (fd, x, &sum)) goto writeerr;
	if (y != NULL && ! write_gwypnum (fd, y, &sum)) goto writeerr; 

/* Write the checksum */

	if (_write (fd, &sum, sizeof (long)) != sizeof (long)) goto writeerr;

/* Save the timers */

	for (i=0; i<NBTIMERS; i++) {
		if (gwyptimers[i+NBTIMERS] != 0.0) {// if the timer was running
			gwypend_timer (i);			// update and save it
			if (! write_double (fd, gwyptimers[i], &sum)) {
				gwypstart_timer (i);	// and then, restart it, even if write is in error!
				goto writeerr;
			}
			if (! write_double (fd, gwyptimers[i+NBTIMERS], &sum)) { // save the timer status
				gwypstart_timer (i);	// and then, restart it, even if write is in error!
				goto writeerr;
			}
			gwypstart_timer (i);	// and then, restart it!
		}
		else {
			if (! write_double (fd, gwyptimers[i], &sum)) goto writeerr;	// save the timer
			if (! write_double (fd, gwyptimers[i+NBTIMERS], &sum)) goto writeerr;	// save its status
		}
	}

	_commit (fd);
	_close (fd);

/* Now rename the intermediate files */

	if (TWO_BACKUP_FILES) {
		_unlink (filename);
		rename (newfilename, filename);
	}
	return (TRUE);

/* An error occured.  Close and delete the current file. */

writeerr:
	_close (fd);
	_unlink (newfilename);
	return (FALSE);
}

int gmreadFromFile (
	char	*filename,
	unsigned long *j,
	unsigned long *ubx,
	unsigned long *uby,
	gwypnum	x,
	gwypnum	y)
{
	int	fd;
	unsigned long magicnum, version;
	long	sum = 0, i;
	char	errmsg[100];

/* Open the intermediate file */

	fd = _open (filename, _O_BINARY | _O_RDONLY);
	if (fd < 0) goto error;

/* Read the file header */

	if (_read (fd, &magicnum, sizeof (long)) != sizeof (long))
		goto readerr;
	if (magicnum != 0x9f2b3cd4) goto readerr;

	if (_read (fd, &version, sizeof (long)) != sizeof (long)) goto readerr;
	if (version != 1 && version != 2) goto readerr;

/* Read the file data */

	if (! read_long (fd, j, &sum)) goto readerr;
	if (! read_long (fd, ubx, &sum)) goto readerr;
	if (! read_long (fd, uby, &sum)) goto readerr;

/* Read the values */

	if (! read_gwypnum (fd, x, &sum)) goto readerr;
	if (y != NULL && ! read_gwypnum (fd, y, &sum)) goto readerr; 

/* Read and compare the checksum */

	if (_read (fd, &i, sizeof (long)) != sizeof (long)) goto readerr;
	if (i != sum) goto readerr;

/* Read the timers and their status */

	for (i=0; i<NBTIMERS; i++) {
		if (! read_double (fd, &gwyptimers[i], &sum)) goto readerr;
		if (! read_double (fd, &gwyptimers[i+NBTIMERS], &sum)) goto readerr;
	}

	_close (fd);

	if (debug) {
		sprintf (errmsg, "Reading %s intermediate file, j = %lu\n", filename, *j);
		OutputBoth (errmsg);
	}

	return (TRUE);

/* An error occured.  Delete the current intermediate file. */
/* Set stage to -1 to indicate an error. */

readerr:
	sprintf (errmsg,"Error reading %s intermediate file, j = %lu\n", filename, *j);
	OutputStr (errmsg);
	_close (fd);
error:
	_unlink (filename);
	return (FALSE);
}

int LwriteToFile (					// To save a Lucas sequence matrix and its Discriminant
	char	*filename,
	unsigned long j,
	unsigned long D,
	unsigned long nr,
	unsigned long *bpf,
	gwypnum	x,
	gwypnum	y,
	gwypnum	z,
	gwypnum	t)
{
	char	newfilename[16], errmsg[100];
	int	fd;
	unsigned long magicnum, version;
	long	sum = 0, i;


	if (debug) {
		sprintf (errmsg, "Writing %s intermediate file, j = %lu\n",filename, j);
		OutputBoth (errmsg);
	}

/* If we are allowed to create multiple intermediate files, then */
/* write to a file called yNNNNNNN. */

	strcpy (newfilename, filename);
	if (TWO_BACKUP_FILES) newfilename[0] = 'y';

/* Create the intermediate file */

	fd = _open (newfilename, _O_BINARY|_O_WRONLY|_O_TRUNC|_O_CREAT, 0666);
	if (fd < 0) return (FALSE);

/* Write the file header. */

	magicnum = 0x9f2b3cd4;
	if (_write (fd, &magicnum, sizeof (long)) != sizeof (long))
		goto writeerr;
	version = 1;
	if (_write (fd, &version, sizeof (long)) != sizeof (long))
		goto writeerr;

/* Write the file data */

	if (! write_long (fd, j, &sum)) goto writeerr;
	if (! write_long (fd, D, &sum)) goto writeerr;
	if (! write_long (fd, nr, &sum)) goto writeerr;
	for (i=0; i<10; i++) {
		if (! write_long (fd, bpf[i], &sum)) goto writeerr;
	}

/* Write the data values */

	if (! write_gwypnum (fd, x, &sum)) goto writeerr;
	if (y != NULL && ! write_gwypnum (fd, y, &sum)) goto writeerr; 
	if (z != NULL && ! write_gwypnum (fd, z, &sum)) goto writeerr; 
	if (t != NULL && ! write_gwypnum (fd, t, &sum)) goto writeerr; 

/* Write the checksum */

	if (_write (fd, &sum, sizeof (long)) != sizeof (long)) goto writeerr;

/* Save the five gwyptimers */

	for (i=0; i<NBTIMERS; i++) {
		if (gwyptimers[i+NBTIMERS] != 0.0) {// if the timer was running
			gwypend_timer (i);			// update and save it
			if (! write_double (fd, gwyptimers[i], &sum)) {
				gwypstart_timer (i);	// and then, restart it, even if write is in error!
				goto writeerr;
			}
			if (! write_double (fd, gwyptimers[i+NBTIMERS], &sum)) { // save the timer status
				gwypstart_timer (i);	// and then, restart it, even if write is in error!
				goto writeerr;
			}
			gwypstart_timer (i);	// and then, restart it!
		}
		else {
			if (! write_double (fd, gwyptimers[i], &sum)) goto writeerr;	// save the timer
			if (! write_double (fd, gwyptimers[i+NBTIMERS], &sum)) goto writeerr;	// save its status
		}
	}

	_commit (fd);
	_close (fd);

/* Now rename the intermediate files */

	if (TWO_BACKUP_FILES) {
		_unlink (filename);
		rename (newfilename, filename);
	}

	return (TRUE);

/* An error occured.  Close and delete the current file. */

writeerr:
	_close (fd);
	_unlink (newfilename);
	return (FALSE);
}

int LreadFromFile (						// To restore a Lucas sequence matrix
	char	*filename,
	unsigned long *j,
	unsigned long *D,
	unsigned long *nr,
	unsigned long *bpf,
	gwypnum	x,
	gwypnum	y,
	gwypnum	z,
	gwypnum	t)
{
	int	fd;
	unsigned long magicnum, version;
	long	sum = 0, i;
	char	errmsg[100];

/* Open the intermediate file */

	fd = _open (filename, _O_BINARY | _O_RDONLY);
	if (fd < 0) goto error;

/* Read the file header */

	if (_read (fd, &magicnum, sizeof (long)) != sizeof (long))
		goto readerr;
	if (magicnum != 0x9f2b3cd4) goto readerr;

	if (_read (fd, &version, sizeof (long)) != sizeof (long)) goto readerr;
	if (version != 1 && version != 2) goto readerr;

/* Read the file data */

	if (! read_long (fd, j, &sum)) goto readerr;
	if (! read_long (fd, D, &sum)) goto readerr;
	if (! read_long (fd, nr, &sum)) goto readerr;
	for (i=0; i<10; i++) {
		if (! read_long (fd, &bpf[i], &sum)) goto readerr;
	}

/* Read the values */

	if (! read_gwypnum (fd, x, &sum)) goto readerr;
	if (y != NULL && ! read_gwypnum (fd, y, &sum)) goto readerr; 
	if (z != NULL && ! read_gwypnum (fd, z, &sum)) goto readerr; 
	if (t != NULL && ! read_gwypnum (fd, t, &sum)) goto readerr; 

/* Read and compare the checksum */

	if (_read (fd, &i, sizeof (long)) != sizeof (long)) goto readerr;
	if (i != sum) goto readerr;

/* Read the timers and their status */

	for (i=0; i<NBTIMERS; i++) {
		if (! read_double (fd, &gwyptimers[i], &sum)) goto readerr;
		if (! read_double (fd, &gwyptimers[i+NBTIMERS], &sum)) goto readerr;
	}

	_close (fd);

	if (debug) {
		sprintf (errmsg, "Reading %s intermediate file, j = %lu\n", filename, *j);
		OutputBoth (errmsg);
	}

	return (TRUE);

/* An error occured.  Delete the current intermediate file. */
/* Set stage to -1 to indicate an error. */

readerr:
	sprintf (errmsg,"Error reading %s intermediate file, j = %lu\n", filename, *j);
	OutputStr (errmsg);
	_close (fd);
error:
	_unlink (filename);
	return (FALSE);
}


/* Print some words of a gwypnum */

int
gwyprint(
	gwypnum 	gg,
	int 	N
)
{
	int 	j;
	long	val;
	char buf[20];

	OutputStr ((char*)"\n");
	for(j=0; j<N; ++j)
	{
		val = (long)gg[j];
		if (val) {
			sprintf (buf, "%ld ", val);
			OutputBoth (buf);
		}
	}
	OutputBoth ((char*)"\n");
	return 0;
}

void writeresidue (
	gwypnum s,			// The gwypnum data
	giant m,			// The current external modulus
	giant t,			// A temporary giant file	
	char *b,			// The output buffer
	const char *str,	// The tested number as a string
	const int bit		// The iteration bit number
)
{
	char restr[20];

	gwyptogiant (s, t);		// The modulo reduction is done here
	modg (m, t);			// External modulus and gwypnums one may be different...
	if (abs(t->sign) < 1)	// make a 64 bit residue correct !!
		sprintf (restr, "%04X%04X%04X%04X", 0, 0, 0, 0);
	else if (abs(t->sign) < 2)
		sprintf (restr, "%04X%04X%04X%04X", 0, 0, 0, t->n[0]);
	else if (abs(t->sign) < 3)
		sprintf (restr, "%04X%04X%04X%04X", 0, 0, t->n[1], t->n[0]);
	else if (abs(t->sign) < 4)
		sprintf (restr, "%04X%04X%04X%04X", 0, t->n[2], t->n[1], t->n[0]);
	else
		sprintf (restr, "%04X%04X%04X%04X", t->n[3], t->n[2], t->n[1], t->n[0]);
	sprintf (b, "%s interim residue %s at bit %d\n", str, restr, bit);
	OutputBoth (b);
}

char res64[17]; /* VP : This variable has been made global */

// Compute the number of digits of a large integer.

int gnbdg (giant nb, int digitbase) {
	giant gnbmax;
	int templ;
	if (digitbase==2)
		return (bitlen (nb));
	else {
		templ = (int)(floor((double)bitlen (nb) * log (2.0) /log ((double)digitbase))); // default value
		gnbmax = newgiant (abs (nb->sign) + 2);
                itog (digitbase, gnbmax);
                power (gnbmax, templ);
		while (gcompg (gnbmax, nb) < 0) {
			templ++;                                      // adjust the value if it is too small
			ulmulg ((unsigned long) digitbase, gnbmax);   // and adjust the comparand
		}
		free (gnbmax);
		return (templ);
	}
}

#if defined (__linux__) || defined (__FreeBSD__) || defined (__APPLE__)

int aprcltest (int prptest, int verbose)		// Primality test using compiled APRCL code
{
	int retcode;
	mpz_t n1;
	mpz_init (n1);
	if (mpz_set_str (n1, greatbuf, 10) != 0)
		return (6);								// Invalid numeric string...
	if (prptest)
		retcode = mpz_strongbpsw_prp (n1);		// Strong BPSW PRP test
	else
		retcode = mpz_aprtcle (n1, verbose);	// Prime test possibly with print out.
	return ((retcode == -1)? 7 : retcode);		// Result code returned (-1 if in error...)
}


#else

	STARTUPINFO aprstinfo;
	PROCESS_INFORMATION aprprinfo;

int aprcltest (int prptest, int verbose)		// Primality test using external APRCL program as a child process
{
	int ofd;
	unsigned long exitcode, errcode;
	char errbuf[100];
	char titre[] = "APRT-CLE Primality Test";
	char line[100] = {0};
	ofd = _open ("__pc__", _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
	if (ofd) {
		writelg = _write (ofd, greatbuf, strlen (greatbuf));
		_close (ofd);
	}
	else
		return (8);				// Could not create the file...

//	Initialisation de la structure wfstinfo avec les seules données requises

	aprstinfo.cb = sizeof (STARTUPINFO);
	aprstinfo.lpDesktop = NULL;
	aprstinfo.lpTitle = titre;
	aprstinfo.dwFlags = STARTF_USESHOWWINDOW;
	aprstinfo.lpReserved2 = NULL;

	exitcode = STILL_ACTIVE;

// Construire la ligne de commande

	if (prptest)
		sprintf (line, "aprcl __pc__ prp");
	else if (verbose == 1)
		sprintf (line, "aprcl __pc__ prime progress");
	else if (verbose == 2)
		sprintf (line, "aprcl __pc__ prime details");
	else
		sprintf (line, "aprcl __pc__");

// Creer le processus fils

	if (!CreateProcess (				// Tentative de creation du processus fils
		"aprcl.exe",
		line,
		NULL,
		NULL,
		FALSE,
		0L,
		NULL,
		NULL,
		&aprstinfo,
		&aprprinfo))

	{
		errcode = GetLastError ();		// Echec...
		sprintf (errbuf, "Error %lu while trying to create new process\n", errcode);
		OutputStr (errbuf);
		_unlink ("__pc__");
		return (9);
	}

	while (exitcode == STILL_ACTIVE) {
		GetExitCodeProcess(aprprinfo.hProcess, &exitcode);
	}
	_unlink ("__pc__");
	return (exitcode);
}

#endif

int saprcltest (char *str, int prptest, int verbose) {
	int i;
	for (i=0;i<10000;i++) {
		greatbuf[i] = str[i];
		if (!str[i])
			break;
	}
	return (aprcltest(prptest, verbose));
}

int gaprcltest (giant N, int prptest, int verbose) {
	if (gnbdg (N,10) > 10000)
		return (-1);				// Number too large...
	if (N->sign == 1) {				// Trial divisions test is sufficient for this small number...
		return (isPrime (N->n[0]) ? 12 : 10);
	}
	gtoc (N, greatbuf, 10000);
	return (aprcltest(prptest, verbose));
}

int setupok (int errcode, giant bignumber, char *string, int *resultat)	                                            // Test if the call to gwypsetup is successful
{
	char buff[256];
        int resaprcl;
	if (!errcode)
		return TRUE;
        else if (errcode == 1) {    // Number too small...
            nbdg = gnbdg (bignumber, 10);
            gwypstart_timer(1);
            if (nbdg > maxaprcl)    // Make only a Strong BPSW PRP test
                resaprcl = gaprcltest (bignumber, 1, 0);
            else if (debug)    // Primality test while showing progress 
                resaprcl = gaprcltest (bignumber, 0, 2);
            else                    // Primality test silently done
                resaprcl = gaprcltest (bignumber, 0, 0);
            gwypend_timer (1);
            if (resaprcl == 10) {
                sprintf (buff,"%s is not prime. (Trial divisions)", string);
            }
            else if (resaprcl == 12)
                sprintf (buff,"%s is prime! (%d decimal digits, Trial divisions)", string, nbdg);
            else if (resaprcl == 0)
			sprintf (buff,"%s is not prime. (APRCL test) ", string);
		else if (resaprcl == 1)
			sprintf (buff,"%s is a probable BPSW prime! (%d decimal digits, APRCL test) ", string, nbdg);
		else if (resaprcl == 2)
			sprintf (buff,"%s is prime! (%d decimal digits, APRCL test)", string, nbdg);
		else if (resaprcl == 6)
			sprintf (buff,"Invalid numerical string in %s\n", string);
		else if (resaprcl == 7)
			sprintf (buff,"APRCL error while testing %s...\n", string);
		else {
			if (resaprcl == 9)
				sprintf (buff, "APRCL primality test not available for %s\n", string);
			else
				sprintf (buff,"Unexpected return value : %d, while APRCL testing %s...\n", resaprcl, string);
		}
        *resultat = ((resaprcl==1)||(resaprcl==2)||(resaprcl==12));

#if defined(WIN32) && !defined(_CONSOLE)

        sprintf (buff+strlen(buff), "  Time : "); 

#else

        clearline(100);

        if (*resultat) {
#ifdef _CONSOLE
		hConsole = GetStdHandle(STD_OUTPUT_HANDLE);	// Access to Console attributes
		SetConsoleTextAttribute(hConsole, BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_RED);
		OutputBoth(buff);
		SetConsoleTextAttribute(hConsole, FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED);
#else
		OutputStr((char*)"\033[7m");
		OutputBoth(buff);
		OutputStr((char*)"\033[0m");
#endif
        }
        else
                OutputBoth(buff);

            sprintf (buff, "  Time : "); 

#endif

/* Output the final timings for APRCL test */

            gwypwrite_timer (buff+strlen(buff), 1, TIMER_CLR | TIMER_NL); 
            OutputBoth (buff);
//	retval = TRUE; 
            return FALSE;
        }
	else {
            sprintf (buff, "%s : Fatal error at setup.", string);
//		gwerror_text (errcode, buff+strlen(buff), 255-strlen(buff));
            strcpy(buff+strlen(buff), "\n");
            OutputBoth (buff);
            return FALSE;
	}
}

/* Test if M divides a^(N-1) - 1 -- gwypsetup has already been called. */

int isexpdiv (
	long a,
	giant N,
	giant M,
	int	*res)
{
//	unsigned long bit, bitpos, firstone = 0, iters;
	unsigned long bit, iters;
	gwypnum	x;
	giant	tmp;
	char	filename[20], buf[sgkbufsize+256], fft_desc[256], oldres64[17];
	long	write_time = DISK_WRITE_TIME * 60;
	int	echk, saving, stopping;
	time_t	start_time, current_time;
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;


/* Init, subtract 1 from N to compute a^(N-1) mod M */

	iaddg (-1, N);

	Nlen = bitlen (N);
	nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.

	*res = TRUE;		/* Assume the residue is one */

/* Init filename */

	tempFileName (filename, 'z', N);

/* Allocate memory */

	x = gwypalloc ();

	tmp = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);

/* Optionally resume from save file and output a message */
/* indicating we are resuming a test */

	if (fileExists (filename) && readFromFile (filename, &bit, x, NULL)) {
		char	fmt_mask[80];
		double	pct;
		pct = trunc_percent (bit * 100.0 / Nlen);
		sprintf (fmt_mask,
			 "Resuming divisibility test of %%d^(N-1)-1 at bit %%ld [%%.%df%%%%]\n",
			 PRECISION);
		sprintf (buf, fmt_mask, a, bit, pct);
		OutputStr (buf);
		if (verbose)
			writeResults (buf);
	}

/* Otherwise, output a message indicating we are starting test */

	else {
		gwypclear_timers ();	// Init. the timers
		sprintf (buf, "Starting divisibility test of %ld^(N-1)-1\n", a);
		OutputStr (buf);
		if (verbose)
			writeResults (buf);
		bit = 1;
//		dbltogw ((double) a, x);
		itogwyp (a, x);
	}

/* Get the current time */

	gwypstart_timer (0);
	gwypstart_timer (1);
	time (&start_time);

/* Output a message about the FFT length */

	gwypfft_description (fft_desc);
#ifdef WIN32
	sprintf (buf, "%s\n", fft_desc);
#else
	sprintf (buf, "%s", fft_desc);
#endif
	OutputStr (buf);
	LineFeed ();
	if (verbose) {
#if !defined(WIN32) 
		strcat (buf, "\n");
#endif
		writeResults (buf);
	}
	ReplaceableLine (1);	/* Remember where replaceable line is */

/* Init the title */

	title ((char*)"Divisibility test in progress...");

/* Do the PRP test */

	gwypsetmulbyconst (a);
	iters = 0;
	while (bit < Nlen) {

/* Error check the last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		stopping = stopCheck ();
		echk = stopping || ERRCHK || (bit <= 50) || (bit >= Nlen-50);
		if (((bit & 127) == 0) || (bit == 1) || (bit == (lasterr_point-1))) {
			echk = 1;
			time (&current_time);
			saving = ((current_time - start_time > write_time) || (bit == 1) || (bit == (lasterr_point-1)));
		} else
			saving = 0;

/* Process this bit */

		if (bitval (N, Nlen-bit-1)) {
			gwypsetnormroutine (0, echk, 1);
		} else {
			gwypsetnormroutine (0, echk, 0);
		}

		if (/*(bit+25 < Nlen) && (bit > 25) && */((bit != lasterr_point) || !maxerr_recovery_mode[6])) {
                    if (cufftonly)
                        gwypsquare (x);
                    else if (zp || generic)
                        cuda_gwypsquare (x, 3);
                    else if(bit==1 || it==0)  
                        {cuda_gwypsquare (x,1);it=1;}
                    else  if(bit != (lasterr_point-1)) 
                        cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                    else
                        cuda_gwypsquare (x,2);
                }
		else {
			gwypsquare_carefully (x);
			will_try_larger_fft = TRUE;
			if (bit == lasterr_point)
				maxerr_recovery_mode[6] = FALSE;
		}

		CHECK_IF_ANY_ERROR (x, (bit), Nlen, 6);

/* That iteration succeeded, bump counters */

		if (will_try_larger_fft && (bit == lasterr_point))
			saving = 1;					// Be sure to restart after this recovery iteration!
		will_try_larger_fft = FALSE;
		bit++;
		iters++;

/* Print a message every so often */

		if (bit % ITER_OUTPUT == 0) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (bit * 100.0 / Nlen);
			sprintf (fmt_mask, "%%.%df%%%% of %%ld", PRECISION);
			sprintf (buf, fmt_mask, pct, Nlen);
			title (buf);
			ReplaceableLine (2);	/* Replace line */
			sprintf (fmt_mask,
				 "%%d^(N-1)-1, bit: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, a, bit, Nlen, pct);
			OutputStr (buf);
			if (ERRCHK && bit > 30) {
				OutputStr ((char*)".  Round off: ");
				sprintf (buf, "%10.10f", reallyminerr);
				OutputStr (buf);
				sprintf (buf, " to %10.10f", reallymaxerr);
				OutputStr (buf);
			}
			gwypend_timer (0);
			if (CUMULATIVE_TIMING) {
				OutputStr ((char*)".  Time thusfar: ");
			} else {
				OutputStr ((char*)".  Time per bit: ");
				gwypdivide_timer (0, iters);
				iters = 0;
			}
			gwypprint_timer (0, TIMER_NL | TIMER_OPT_CLR);
			gwypstart_timer (0);
		}

/* Print a results file message every so often */

		if (bit % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
			sprintf (buf, "Bit %ld / %ld\n", bit, Nlen);
			writeResults (buf);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving || stopping) {
			write_time = DISK_WRITE_TIME * 60;
			if (! writeToFile (filename, bit, x, NULL)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (buf);
				if (write_time > 600) write_time = 600;
			}	
			time (&start_time);

/* If an escape key was hit, write out the results and return */

			if (stopping) {
				iaddg (1, N);	// Restore the modulus
				gwypfree (tmp);
//				gwypfree (x);
//				gwypdone ();
				*res = FALSE;
				return (FALSE);
			}
		}

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next iteration so that we can compare our */
/* residues to programs that start counter at zero or one. */

		if (interimResidues && bit % interimResidues < 2) {
				gwyptogiant (x, tmp);		// The modulo reduction is done here
			if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
			else if (abs(tmp->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
			else if (abs(tmp->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
			sprintf (buf, "%ld^(N-1) interim residue %s at bit %ld\n", a, res64, bit);
			OutputBoth (buf);
		}

/* Write a save file every "interimFiles" iterations. */

		if (interimFiles && bit % interimFiles == 0) {
			char	interimfile[20];
			sprintf (interimfile, "%.8s.%03lu",
				 filename, bit / interimFiles);
			if (! writeToFile (interimfile, bit, x, NULL)) {
				sprintf (buf, WRITEFILEERR, interimfile);
				OutputBoth (buf);
			}
		}
	}

/* See if the residue is one.  If not, format a 64-bit residue. */
/* Old versions of PRP used a non-standard 64-bit residue, computing */
/* a^N-a mod N rather than the more standard a^(N-1) mod N.  Since */
/* some projects recorded these non-standard residues, output that */
/* residue too.  Note that some really old versions printed out the */
/* 32-bit chunks of the non-standard residue in reverse order. */

	clearline (100);

	gwyptogiant (x, tmp);
	if (!isone (tmp)) {
		*res = FALSE;	/* Residue not one */
			if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
			else if (abs(tmp->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
			else if (abs(tmp->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
		smulg ((unsigned short)a, tmp); modg (N, tmp); iaddg (-a, tmp);
			if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
				sprintf (oldres64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
			else if (abs(tmp->sign) < 3)
				sprintf (oldres64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
			else if (abs(tmp->sign) < 4)
				sprintf (oldres64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
			else
				sprintf (oldres64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
	}

/* Print results.  Do not change the format of this line as Jim Fougeron of */
/* PFGW fame automates his QA scripts by parsing this line. */

	sprintf (buf, "End of divisibility test of %ld^(N-1)-1\n", a);

	gwypfree (tmp);
	gwypfree (x);

/* Output the final timings */

	gwypend_timer (1);
	sprintf (buf+strlen(buf)-1, "  Time: ");
	ReplaceableLine (2);	/* Replace line */
	gwypwrite_timer (buf+strlen(buf), 1, TIMER_CLR | TIMER_NL); 
	if (verbose)
		OutputBoth (buf);
	else
		OutputStr (buf);

/* Cleanup and return */

	iaddg (1, N);					// Restore the modulus
	Nlen = bitlen (N);
//	gwypdone ();
	_unlink (filename);
	lasterr_point = 0;
	return (TRUE);

/* An error occured, sleep, then try restarting at last save point. */

error:
	iaddg (1, N);					// Restore the value of N
	Nlen = bitlen (N);
	gwypfree (tmp);
	gwypfree (x);
//	gwypdone ();
	*res = FALSE;					// To avoid credit mesage...

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
		sprintf (buf,ERRMSG7,a);
		OutputBoth (buf);
		gwypdone ();
		will_try_larger_fft = FALSE;
		_unlink (filename);
		return (TRUE);
	}

//	gwypdone ();

/* Output a message saying we are restarting */
	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
		will_try_larger_fft = FALSE;
		return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
		OutputBoth (ERRMSG8);
		IniWriteInt(INI_FILE, (char*)"FFT_Increment", IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
		_unlink (filename);
		will_try_larger_fft = FALSE;
	}
	return (-1);
}

/* Test for an N+1 Frobenius probable prime -- gwypsetup has already been called. */

int commonFrobeniusPRP (
	unsigned long P,
	unsigned long Q,
	int	*res, char *str)
{
//	unsigned long bit, firstone = 0, iters, D, bitv;
	unsigned long bit, iters, D;
	gwypnum x, y, gwA, gw2;
	giant	tmp, tmp2, tmp3, A;
	char	filename[20], buf[sgkbufsize+256], fft_desc[256];
	long	write_time = DISK_WRITE_TIME * 60;
	int	echk, saving, stopping;
	time_t	start_time, current_time;
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;


/* Allocate memory */

	x = gwypalloc ();
	y = gwypalloc ();
	gwA = gwypalloc ();
	gw2 = gwypalloc ();
//	dbltogw (2.0, gw2);
	itogwyp (2, gw2);

	tmp = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);
	tmp2 = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);
	tmp3 = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);
	A = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);

	D = P*P - 4*Q;
	*res = TRUE;		/* Assume it is a probable prime */
        nbdg = gnbdg (N, 10);
/* Init file name */

	tempFileName (filename, 'F', N);

/* Optionally resume from save file and output a message */
/* indicating we are resuming a Frobenius test */

	if (fileExists (filename) && readFromFile (filename, &bit, x, y)) {
		char	fmt_mask[80];
		double	pct;
		pct = trunc_percent (bit * 100.0 / Nlen);
		sprintf (fmt_mask,
			 "Resuming Frobenius sequence at bit %%ld [%%.%df%%%%]\n",
			 PRECISION);
		sprintf (buf, fmt_mask, bit, pct);
		gwypstart_timer (0);
		gwypstart_timer (1);
		goto Frobeniusresume;
	}
	else {

/* Init, compute (N+1)/2 to compute x and y mod N */

		gtog (N, tmp3);
		iaddg (1, tmp3);
		gshiftright (1, tmp3);
		Nlen = bitlen (tmp3);

		itog (Q, A);		// Compute A = P*P*Q^-1 - 2 mod N
		invg (N, A);
		ulmulg (P*P, A);
		iaddg (-2, A);
		modg (N, A);
		gianttogwyp (A, gwA);

		tempFileName (filename, 'L', N);	// Reinit file name

/* Optionally resume from save file and output a message */
/* indicating we are resuming a Lucas test */

		if (fileExists (filename) && readFromFile (filename, &bit, x, y)) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (bit * 100.0 / Nlen);
			sprintf (fmt_mask,
				"Resuming Lucas sequence at bit %%ld [%%.%df%%%%]\n",
				PRECISION);
			sprintf (buf, fmt_mask, bit, pct);
			OutputStr (buf);
			if (verbose)
				writeResults (buf);
		}

/* Otherwise, output a message indicating we are starting a Lucas test */

		else {
//			gwypclear_timers ();	// Init. timers
                        if (showdigits)
                            sprintf (buf, "Starting Lucas sequence (%d decimal digits)\n", nbdg);
                        else
                            sprintf (buf, "Starting Lucas sequence\n");
			OutputStr (buf);
			if (verbose)
				writeResults (buf);

			bit = 1;
//			dbltogw (2.0, x);			// Initial values
			itogwyp (2, x);		// Initial values
			gwypcopy (gwA, y);
//			if (! writeToFile (filename, bit, x, y)) {
//				sprintf (buf, WRITEFILEERR, filename);
//				OutputBoth (buf);
//			}	
		}
	}


/* Get the current time */

	gwypstart_timer (0);
	gwypstart_timer (1);
	time (&start_time);

/* Output a message about the FFT length */

	gwypfft_description (fft_desc);
#ifdef WIN32
	sprintf (buf, "%s, P = %lu, Q = %lu\n", fft_desc, P, Q);
#else
	sprintf (buf, "%s, P = %lu, Q = %lu", fft_desc, P, Q);
#endif
	OutputStr (buf);
	LineFeed ();
	if (verbose) {
#if !defined(WIN32) 
		strcat (buf, "\n");
#endif
		writeResults (buf);
	}
	ReplaceableLine (1);	/* Remember where replaceable line is */

/* Init the title */

	title ((char*)"Lucas PRP test in progress...");

/* Do the PRP test */

	will_try_larger_fft = FALSE;
	iters = 0;
	while (bit <= Nlen) {

/* Error check the first and last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		stopping = stopCheck ();
		echk = stopping || ERRCHK || (bit <= 50) || (bit >= Nlen-50);
		if (((bit & 127) == 0) || (bit == 1) || (bit == (lasterr_point-1))) {
			echk = 1;
			time (&current_time);
			saving = ((current_time - start_time > write_time) || (bit == 1) || (bit == (lasterr_point-1)));
		} else
			saving = 0;

/* Process this bit */

		gwypsetnormroutine (0, echk, 0);

		if (bitval (tmp3, Nlen-bit)) {
			gwypsetaddin (0);
			if (/*(bit+26 < Nlen) && (bit > 26) && */
				((bit != lasterr_point) || (!maxerr_recovery_mode[1] && !maxerr_recovery_mode[2]))) {
                            if (cufftonly)
                                gwypmul (y, x);
                            else
                                cuda_gwypmul (y, x, 3);
                            will_try_larger_fft = FALSE;
			}
			else {
                            gwypmul_carefully (y, x);
                            if (bit == lasterr_point)
                                    maxerr_recovery_mode[1] = FALSE;
                            will_try_larger_fft = TRUE;
			}
			CHECK_IF_ANY_ERROR(x, (bit), Nlen, 1)
			gwypsub3 (x, gwA, x);
			if (abs(inc)==1)
				gwypsetaddin (-2);
			if (/*(bit+26 < Nlen) && (bit > 26) && */
				((bit != lasterr_point) || !maxerr_recovery_mode[2])) {
                            if (cufftonly)
                                gwypsquare (y);
                            else if (zp || generic)
                                cuda_gwypsquare (y, 3);
                            else if(bit==1 || it==0)  
                                {cuda_gwypsquare (y,1);it=1;}
                            else  if(bit != (lasterr_point-1)) 
                                cuda_gwypsquare (y,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (y,2);
                            will_try_larger_fft = FALSE;
			}
			else {
                            gwypsquare_carefully (y);
                            if (bit == lasterr_point)
                                    maxerr_recovery_mode[2] = FALSE;
                            will_try_larger_fft = TRUE;
			}
			CHECK_IF_ANY_ERROR(y, (bit), Nlen, 2)
			if (abs(inc)!=1)
				gwypsubquick (gw2, y);
		}
		else {
			gwypsetaddin (0);
			if (/*(bit+26 < Nlen) && (bit > 26) && */
				((bit != lasterr_point) || (!maxerr_recovery_mode[3] && !maxerr_recovery_mode[4]))) {
                            if (cufftonly)
                                gwypmul (x, y);
                            else
                                cuda_gwypmul (x, y, 3);
                            will_try_larger_fft = FALSE;
			}
			else {
				gwypmul_carefully (x, y);
				if (bit == lasterr_point)
					maxerr_recovery_mode[3] = FALSE;
				will_try_larger_fft = TRUE;
			}
			CHECK_IF_ANY_ERROR(y, (bit), Nlen, 3)
			gwypsub3 (y, gwA, y);
			if (abs(inc)==1)
				gwypsetaddin (-2);
			if (/*(bit+26 < Nlen) && (bit > 26) && */
				((bit != lasterr_point) || !maxerr_recovery_mode[4])) {
                            if (cufftonly)
                                gwypsquare (x);
                            else if (zp || generic)
                                cuda_gwypsquare (x, 3);
                            else if(bit==1 || it==0)  
                                {cuda_gwypsquare (x,1);it=1;}
                            else  if(bit != (lasterr_point-1)) 
                                cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (x,2);
                            will_try_larger_fft = FALSE;
			}
			else {
				gwypsquare_carefully (x);
				if (bit == lasterr_point)
					maxerr_recovery_mode[4] = FALSE;
				will_try_larger_fft = TRUE;
			}
			CHECK_IF_ANY_ERROR(x, (bit), Nlen, 4)
			if (abs(inc)!=1)
				gwypsubquick (gw2, x);
		}

 /* That iteration succeeded, bump counters */

		if (will_try_larger_fft && (bit == lasterr_point))
			saving = 1;					// Be sure to restart after this recovery iteration!
		will_try_larger_fft = FALSE;
		bit++;
		iters++;

/* Print a message every so often */

		if (bit % ITER_OUTPUT == 0) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (bit * 100.0 / Nlen);
			sprintf (fmt_mask, "%%.%df%%%% of %%ld", PRECISION);
			sprintf (buf, fmt_mask, pct, Nlen);
			title (buf);
			ReplaceableLine (2);	/* Replace line */
			sprintf (fmt_mask,
				 "%%s, bit: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, str, bit, Nlen, pct);
			OutputStr (buf);
			if (ERRCHK && bit > 30) {
				OutputStr ((char*)".  Round off: ");
				sprintf (buf, "%10.10f", reallyminerr);
				OutputStr (buf);
				sprintf (buf, " to %10.10f", reallymaxerr);
				OutputStr (buf);
			}
			gwypend_timer (0);
			if (CUMULATIVE_TIMING) {
				OutputStr ((char*)".  Time thusfar: ");
			} else {
				OutputStr ((char*)".  Time per bit: ");
				gwypdivide_timer (0, iters);
				iters = 0;
			}
			gwypprint_timer (0, TIMER_NL | TIMER_OPT_CLR);
			gwypstart_timer (0);
		}

/* Print a results file message every so often */

		if (bit % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
			sprintf (buf, "Bit %ld / %ld\n", bit, Nlen);
			writeResults (buf);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving || stopping) {
			write_time = DISK_WRITE_TIME * 60;
			saving = FALSE;
			if (! writeToFile (filename, bit, x, y)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (buf);
				if (write_time > 600) write_time = 600;
			}	
			time (&start_time);

/* If an escape key was hit, write out the results and return */

			if (stopping) {
				gwypfree (tmp);
				gwypfree (tmp2);
				gwypfree (tmp3);
				gwypfree (A);
				gwypfree (x);				// Clean up
				gwypfree (y);
				gwypfree (gwA);
				gwypfree (gw2);
				*res = FALSE;
				return (FALSE);
			}
		}

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next iteration so that we can compare our */
/* residues to programs that start counter at zero or one. */

		if (interimResidues && bit % interimResidues < 2) {
				gwyptogiant (x, tmp);		// The modulo reduction is done here
				modg (N, tmp);				// External modulus and gwypnums one may be different...
			if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
			else if (abs(tmp->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
			else if (abs(tmp->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
			sprintf (buf, "%s interim residue %s at bit %ld\n", str, res64, bit);
			OutputBoth (buf);
		}

/* Write a save file every "interimFiles" iterations. */

		if (interimFiles && bit % interimFiles == 0) {
			char	interimfile[20];
			sprintf (interimfile, "%.8s.%03lu",
				 filename, bit / interimFiles);
			if (! writeToFile (interimfile, bit, x, y)) {
				sprintf (buf, WRITEFILEERR, interimfile);
				OutputBoth (buf);
			}
		}
	}

	will_try_larger_fft = FALSE;

/* See if we've found a Lucas probable prime.  If not, format a 64-bit residue. */

	clearline (100);

	gwyptogiant (x, tmp);			// V(m)
	gwyptogiant (y, tmp2);			// V(m+1)
	mulg (A, tmp);					// A*V(m)
	gshiftleft (1, tmp2);			// 2*V(m+1)
	subg (tmp2, tmp);				// A*V(m)-2*V(m+1)
	modg (N, tmp);

	if (!isZero (tmp)) {
		*res = FALSE;				// N is composite.
			if (abs(tmp->sign) < 2)	// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
			else if (abs(tmp->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
			else if (abs(tmp->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);

		if (IniGetInt(INI_FILE, (char*)"LucasPRPtest", 0))
			sprintf (buf, "%s is not prime(P = %lu, Q = %lu), Lucas RES64: %s", str, P, Q, res64);
		else
			sprintf (buf, "%s is strong-Fermat PSP, but composite!! (P = %lu, Q = %lu), Lucas RES64: %s", str, P, Q, res64);
	}
	if (*res) {						// N may be prime ; do now the Frobenius PRP test
		_unlink (filename);			// Remove Lucas save file

		tempFileName (filename, 'F', N);	// Reinit file name

		bit = 1;
//		dbltogw ((double)Q, y);
		itogwyp (Q, y);
//		lasterr_point = 0;			// Reset a possible Lucas roundoff error point
		sprintf (buf, "%s is Lucas PRP, Starting Frobenius test sequence\n", str);
		if (! writeToFile (filename, bit, x, y)) {
			sprintf (buf, WRITEFILEERR, filename);
			OutputBoth (buf);
		}			// Make a save file to avoid a false restart after a roundoff error...

Frobeniusresume:

	clearline (100);

		OutputStr (buf);
		if (verbose)
			writeResults (buf);
		gtog (N, tmp3);				// compute (N-1)/2 to compute y mod N
		iaddg (-1, tmp3);
		gshiftright (1, tmp3);
		Nlen = bitlen (tmp3);

/* Get the current time */

		time (&start_time);

/* Output a message about the FFT length */

		gwypfft_description (fft_desc);
#ifdef WIN32
		sprintf (buf, "%s, Q = %lu\n", fft_desc, Q);
#else
		sprintf (buf, "%s, Q = %lu", fft_desc, Q);
#endif
		OutputStr (buf);
		LineFeed();
		if (verbose) {
#if !defined(WIN32) 
			strcat (buf, "\n");
#endif
			writeResults (buf);
		}
		ReplaceableLine (1);	/* Remember where replaceable line is */

/* Init the title */

		title ((char*)"Frobenius PRP test in progress...");

/* Do the PRP test */
		
//		asm_data = (struct gwasm_data *) ->asm_data;	// Get the struct pointer
//		asm_data->SPREAD_CARRY_OVER_4_WORDS = TRUE;
		gwypsetaddin (0);
		gwypsetmulbyconst (Q);
		iters = 0;
		while (bit < Nlen) {

/* Error check the first and last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

			stopping = stopCheck ();
			echk = stopping || ERRCHK || (bit <= 50) || (bit >= Nlen-50);
			if (((bit & 127) == 0) || (bit == 1) || (bit == (lasterr_point-1))) {
				echk = 1;
				time (&current_time);
				saving = ((current_time - start_time > write_time) || (bit == 1) || (bit == (lasterr_point-1)));
			} else
				saving = 0;

/* Process this bit */


			if (bitval (tmp3, Nlen-bit-1)) {
				gwypsetnormroutine (0, echk, 1);
			} else {
				gwypsetnormroutine (0, echk, 0);
			}

			if (/*(bit+25 < Nlen) && (bit > 25) && */((bit != lasterr_point) || !maxerr_recovery_mode[6])) {
                            if (cufftonly)
                                gwypsquare (y);
                            else if (zp || generic)
                                cuda_gwypsquare (y, 3);
                            else if(bit==1 || it==0)  
                                {cuda_gwypsquare (y,1);it=1;}
                            else  if(bit != (lasterr_point-1)) 
                                cuda_gwypsquare (y,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (y,2);
                        }
			else {
				gwypsquare_carefully (y);
				will_try_larger_fft = TRUE;
				if (bit == lasterr_point)
					maxerr_recovery_mode[6] = FALSE;
			}

			CHECK_IF_ANY_ERROR (y, (bit), Nlen, 6);

/* That iteration succeeded, bump counters */

			if (will_try_larger_fft && (bit == lasterr_point))
				saving = 1;					// Be sure to restart after this recovery iteration!
			will_try_larger_fft = FALSE;
			bit++;
			iters++;

/* Print a message every so often */

			if (bit % ITER_OUTPUT == 0) {
				char	fmt_mask[80];
				double	pct;
				pct = trunc_percent (bit * 100.0 / Nlen);
				sprintf (fmt_mask, "%%.%df%%%% of %%ld", PRECISION);
				sprintf (buf, fmt_mask, pct, Nlen);
				title (buf);
				ReplaceableLine (2);	/* Replace line */
				sprintf (fmt_mask,
					 "%%s, bit: %%ld / %%ld [%%.%df%%%%]",
					 PRECISION);
				sprintf (buf, fmt_mask, str, bit, Nlen, pct);
				OutputStr (buf);
				if (ERRCHK && bit > 30) {
					OutputStr ((char*)".  Round off: ");
					sprintf (buf, "%10.10f", reallyminerr);
					OutputStr (buf);
					sprintf (buf, " to %10.10f", reallymaxerr);
					OutputStr (buf);
				}
				gwypend_timer (0);
				if (CUMULATIVE_TIMING) {
					OutputStr ((char*)".  Time thusfar: ");
				} else {
					OutputStr ((char*)".  Time per bit: ");
					gwypdivide_timer (0, iters);
					iters = 0;
				}
				gwypprint_timer (0, TIMER_NL | TIMER_OPT_CLR);
				gwypstart_timer (0);
			}

/* Print a results file message every so often */

			if (bit % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
				sprintf (buf, "Bit %ld / %ld\n", bit, Nlen);
				writeResults (buf);
			}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

			if (saving || stopping) {
				write_time = DISK_WRITE_TIME * 60;
				saving = 0;
				if (! writeToFile (filename, bit, x, y)) {
					sprintf (buf, WRITEFILEERR, filename);
					OutputBoth (buf);
					if (write_time > 600) write_time = 600;
				}	
				time (&start_time);

/* If an escape key was hit, write out the results and return */

				if (stopping) {
					gwypfree (tmp);
					gwypfree (tmp2);
					gwypfree (tmp3);
					gwypfree (A);
					gwypfree (x);
					gwypfree (y);
					gwypfree (gwA);
					gwypfree (gw2);
//					gwypdone ();
					*res = FALSE;
					return (FALSE);
				}
			}

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next iteration so that we can compare our */
/* residues to programs that start counter at zero or one. */

			if (interimResidues && bit % interimResidues < 2) {
				gwyptogiant (y, tmp);	// The modulo reduction is done here
			if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
			else if (abs(tmp->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
			else if (abs(tmp->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
			sprintf (buf, "%s interim residue %s at bit %ld\n", str, res64, bit);
				OutputBoth (buf);
			}

/* Write a save file every "interimFiles" iterations. */

			if (interimFiles && bit % interimFiles == 0) {
				char	interimfile[20];
				sprintf (interimfile, "%.8s.%03lu",
				 filename, bit / interimFiles);
				if (! writeToFile (interimfile, bit, x, y)) {
					sprintf (buf, WRITEFILEERR, interimfile);
					OutputBoth (buf);
				}
			}
		}
		gwypsetaddin (0);
		gwypsetnormroutine (0, 1, 0);
		will_try_larger_fft = TRUE;
//		gwypmul_carefully (x, y);	// y = B*V(m)-2
                if (cufftonly)
                    gwypmul (x, y);
                else
                    cuda_gwypmul (x, y, 3);
		CHECK_IF_ANY_ERROR (y, (Nlen), Nlen, 6);
		will_try_larger_fft = FALSE;
		gwypsubquick (gw2, y);
		gwyptogiant (y, tmp);
		modg (N, tmp);					// modulo N
		if (!isZero (tmp)) {
			*res = FALSE;				// N is Lucas PSP, but composite!!
			if (abs(tmp->sign) < 2)	// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
			else if (abs(tmp->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
			else if (abs(tmp->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);

			if (IniGetInt(INI_FILE, (char*)"LucasPRPtest", 0))
				sprintf (buf, "%s is Lucas PSP (P = %lu, Q = %lu), but composite!!. Frobenius RES64: %s", str, P, Q, res64);
			else
				sprintf (buf, "%s is strong-Fermat and Lucas PSP (P = %lu, Q = %lu), but composite!!. Frobenius RES64: %s", str, P, Q, res64);
		}
	}

/* Print results.  */

	clearline (100);

	if (*res)
		sprintf (buf, "%s is Frobenius PRP! (P = %lu, Q = %lu, D = %lu)", str, P, Q, D);

	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (tmp3);
	gwypfree (A);
	gwypfree (x);				// Clean up
	gwypfree (y);
	gwypfree (gwA);
	gwypfree (gw2);

/* Update the output file */

	if ((*res && IniGetInt (INI_FILE, (char*)"OutputPrimes", 0)) ||
	    (!*res && IniGetInt (INI_FILE, (char*)"OutputComposites", 0)))
		writeResults (buf);


#if defined(WIN32) && !defined(_CONSOLE)

	sprintf (buf+strlen(buf), "  Time : "); 
	ReplaceableLine (2);	/* Replace line */ 

#else

	clearline(100);

#ifdef _CONSOLE
	OutputBoth(buf);
#else
	if (*res) {
		OutputStr((char*)"\033[7m");
		OutputBoth(buf);
		OutputStr((char*)"\033[0m");
	}
	else
		OutputBoth(buf);
#endif

	sprintf (buf, "  Time : "); 

#endif

/* Output the final timings */

	gwypend_timer (1);
	gwypwrite_timer (buf+strlen(buf), 1, TIMER_CLR | TIMER_NL); 
	if ((*res && IniGetInt (INI_FILE, (char*)"OutputPrimes", 0)) ||
	    (!*res && IniGetInt (INI_FILE, (char*)"OutputComposites", 0)))
		OutputStr (buf);
	else
		OutputBoth (buf);

/* Cleanup and return */

	Nlen = bitlen (N);
//	gwypdone ();
	_unlink (filename);
	lasterr_point = 0;
	return (TRUE);

/* An error occured, sleep, then try restarting at last save point. */

error:
	Nlen = bitlen (N);
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (tmp3);
	gwypfree (A);
	gwypfree (x);				// Clean up
	gwypfree (y);
	gwypfree (gwA);
	gwypfree (gw2);
//	gwypdone ();
	*res = FALSE;

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
		sprintf (buf, ERRMSG5, checknumber, str);
		OutputBoth (buf);
//		gwypdone ();
		will_try_larger_fft = FALSE;
		_unlink (filename);
		return (TRUE);
	}

//	gwypdone ();

/* Output a message saying we are restarting */
	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
		will_try_larger_fft = FALSE;
		return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
		OutputBoth (ERRMSG8);
		IniWriteInt(INI_FILE, (char*)"FFT_Increment", IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
		_unlink (filename);
		will_try_larger_fft = FALSE;
	}
	return (-1);
}

/* Test for a strong Fermat probable prime -- gwypsetup has already been called. */

int commonPRP (
	long a,
	int	*res, char *str)
{
	unsigned long bit, bitpos, firstone = 0, iters;
	gwypnum	x, y, gwypminusone, gwypone;
	giant	tmp;
	char	filename[20], buf[sgkbufsize+256], fft_desc[256], oldres64[17];
	long	write_time = DISK_WRITE_TIME * 60;
	int	echk, saving, stopping;
	time_t	start_time, current_time;
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;


/* Init, subtract 1 from N to compute a^(N-1) mod N */

	iaddg (-1, N);
	while (bitval (N, firstone) == 0)	// position of first one bit in N-1
		firstone++;
	Nlen = bitlen (N);
	nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.
	*res = TRUE;		/* Assume it is a probable prime */

/* Init filename */

	tempFileName (filename, 'z', N);

/* Allocate memory */

	x = gwypalloc ();
	y = gwypalloc ();
	gwypminusone = gwypalloc ();
	gwypone = gwypalloc ();

//	dbltogw (1.0, gwypone);
//	gianttogwyp (N, gwypminusone);
	itogwyp (1, gwypone);
	if (abs(inc) == 1)
		itogwyp (-1, gwypminusone);
	else
		gianttogwyp (N, gwypminusone);
	tmp = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);

/* Optionally resume from save file and output a message */
/* indicating we are resuming a test */

	if (fileExists (filename) && readFromFile (filename, &bit, x, NULL)) {
		char	fmt_mask[80];
		double	pct;
		pct = trunc_percent (bit * 100.0 / Nlen);
		sprintf (fmt_mask,
			 "Resuming probable prime test of %%s at bit %%ld [%%.%df%%%%]\n",
			 PRECISION);
		sprintf (buf, fmt_mask, str, bit, pct);
		OutputStr (buf);
		if (verbose)
			writeResults (buf);
	}

/* Otherwise, output a message indicating we are starting test */

	else {
		gwypclear_timers ();	// Init. timers
		if (showdigits)
			sprintf (buf, "Starting probable prime test of %s (%d decimal digits)\n", str, nbdg);
		else
			sprintf (buf, "Starting probable prime test of %s\n", str);
		OutputStr (buf);
		if (verbose)
			writeResults (buf);
		bit = 1;
		itogwyp (a, x);
	}

/* Get the current time */

	gwypstart_timer (0);
	gwypstart_timer (1);
	time (&start_time);

/* Output a message about the FFT length */

	gwypfft_description (fft_desc);
#ifdef WIN32
	sprintf (buf, "%s, a = %ld\n", fft_desc, a);
#else
	sprintf (buf, "%s, a = %ld", fft_desc, a);
#endif
	OutputStr (buf);
	LineFeed();
	if (verbose) {
#if !defined(WIN32) 
		strcat (buf, "\n");
#endif
		writeResults (buf);
	}
	ReplaceableLine (1);	/* Remember where replaceable line is */

/* Init the title */

	title ((char*)"Strong Fermat PRP test in progress...");

/* Do the PRP test */

	gwypsetmulbyconst (a);
	iters = 0;
	while (bit < Nlen) {

/* Error check the first and last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		stopping = stopCheck ();
		echk = stopping || ERRCHK || (bit <= 50) || (bit >= Nlen-50);
		if (((bit & 127) == 0) || (bit == 1) || (bit == (lasterr_point-1))) {
			echk = 1;
			time (&current_time);
			saving = ((current_time - start_time > write_time) || (bit == 1) || (bit == (lasterr_point-1)));
		} else
			saving = 0;

/* Process this bit */


		if (bitval (N, bitpos = Nlen-bit-1)) {
			gwypsetnormroutine (0, echk, 1);
		} else {
			gwypsetnormroutine (0, echk, 0);
		}

		if (/*(bit+25 < Nlen) && (bit > 25) && */((bit != lasterr_point) || !maxerr_recovery_mode[6])) {
                    if (cufftonly)
                        gwypsquare (x);
                    else if (zp || generic)
                        cuda_gwypsquare (x, 3);
                    else if(bit==1 || it==0)  
                        {cuda_gwypsquare (x,1);it=1;}
                    else  if(bit != (lasterr_point-1)) 
                        cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                    else
                        cuda_gwypsquare (x,2);
                }
		else {
			gwypsquare_carefully (x);
//			gwypcareful_squaring (x);
			will_try_larger_fft = TRUE;
			if (bit == lasterr_point)
				maxerr_recovery_mode[6] = FALSE;
		}

		CHECK_IF_ANY_ERROR (x, (bit), Nlen, 6);


		if (bitpos == firstone)
			gwypcopy (x, y);			// Keep this value for the strong PRP test

/* That iteration succeeded, bump counters */

		if (will_try_larger_fft && (bit == lasterr_point))
			saving = 1;					// Be sure to restart after this recovery iteration!
		will_try_larger_fft = FALSE;
		bit++;
		iters++;

/* Print a message every so often */

		if (bit % ITER_OUTPUT == 0) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (bit * 100.0 / Nlen);
			sprintf (fmt_mask, "%%.%df%%%% of %%ld", PRECISION);
			sprintf (buf, fmt_mask, pct, Nlen);
			title (buf);
			ReplaceableLine (2);	/* Replace line */
			sprintf (fmt_mask,
				 "%%s, bit: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, str, bit, Nlen, pct);
			OutputStr (buf);
			if (ERRCHK && bit > 30) {
				OutputStr ((char*)".  Round off: ");
				sprintf (buf, "%10.10f", reallyminerr);
				OutputStr (buf);
				sprintf (buf, " to %10.10f", reallymaxerr);
				OutputStr (buf);
			}
			gwypend_timer (0);
			if (CUMULATIVE_TIMING) {
				OutputStr ((char*)".  Time thusfar: ");
			} else {
				OutputStr ((char*)".  Time per bit: ");
				gwypdivide_timer (0, iters);
				iters = 0;
			}
			gwypprint_timer (0, TIMER_NL | TIMER_OPT_CLR);
			gwypstart_timer (0);
		}

/* Print a results file message every so often */

		if (bit % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
			sprintf (buf, "Bit %ld / %ld\n", bit, Nlen);
			writeResults (buf);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving || stopping) {
			write_time = DISK_WRITE_TIME * 60;
			if (! writeToFile (filename, bit, x, NULL)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (buf);
				if (write_time > 600) write_time = 600;
			}	
			time (&start_time);

/* If an escape key was hit, write out the results and return */

			if (stopping) {
				gwypfree (tmp);
				gwypfree (x);
				gwypfree (y);
				gwypfree (gwypminusone);
				gwypfree (gwypone);
//				gwypdone ();
				*res = FALSE;
				return (FALSE);
			}
		}

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next iteration so that we can compare our */
/* residues to programs that start counter at zero or one. */

		if (interimResidues && bit % interimResidues < 2) {
				gwyptogiant (x, tmp);	// The modulo reduction is done here
			if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
			else if (abs(tmp->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
			else if (abs(tmp->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
			sprintf (buf, "%s interim residue %s at bit %ld\n", str, res64, bit);
			OutputBoth (buf);
		}

/* Write a save file every "interimFiles" iterations. */

		if (interimFiles && bit % interimFiles == 0) {
			char	interimfile[20];
			sprintf (interimfile, "%.8s.%03lu",
				 filename, bit / interimFiles);
			if (! writeToFile (interimfile, bit, x, NULL)) {
				sprintf (buf, WRITEFILEERR, interimfile);
				OutputBoth (buf);
			}
		}
	}

/* See if we've found a probable prime.  If not, format a 64-bit residue. */
/* Old versions of PRP used a non-standard 64-bit residue, computing */
/* 3^N-3 mod N rather than the more standard 3^(N-1) mod N.  Since */
/* some projects recorded these non-standard residues, output that */
/* residue too.  Note that some really old versions printed out the */
/* 32-bit chunks of the non-standard residue in reverse order. */

	clearline (100);

	iaddg (1, N);	// Restore the modulus

	gwyptogiant (x, tmp);
	modg (N, tmp);

	if (!isone (tmp)) {
		*res = FALSE;	/* Not a prime */
		if (abs(tmp->sign) < 2)	// make a 64 bit residue correct !!
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
		else if (abs(tmp->sign) < 3)
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
		else if (abs(tmp->sign) < 4)
			sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
		else
			sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
		smulg ((unsigned short)a, tmp); modg (N, tmp); iaddg (-a, tmp);
		if (abs(tmp->sign) < 2)	// make a 64 bit residue correct !!
			sprintf (oldres64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
		else if (abs(tmp->sign) < 3)
			sprintf (oldres64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
		else if (abs(tmp->sign) < 4)
			sprintf (oldres64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
		else
			sprintf (oldres64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
		if (IniGetInt (INI_FILE, (char*)"OldRes64", 1))
			sprintf (buf, "%s is not prime.  RES64: %s.  OLD64: %s", str, res64, oldres64);
		else
			sprintf (buf, "%s is not prime.  RES64: %s", str, res64);
	}
	else {
		if (strong) {
			iters = 1;
			gwyptogiant (y, tmp);
			if (!isone (tmp)) {
				for (bitpos = 0; bitpos < firstone; bitpos++) {
					iaddg (1, tmp);
					iters++;
					if (gcompg (N, tmp) == 0)	// success!
						break;
					if (firstone - bitpos == 1) {
						*res = FALSE;			// N is PSP, but composite!!
						sprintf (buf, "%s is not prime, although %ld-PSP!! (%lu more test(s))", str, a, iters);
						break;
					}
                                        if (cufftonly)
                                            gwypsquare (y);
                                        else if (zp || generic)
                                            cuda_gwypsquare (y, 3);
                                        else if(bit==1 || it==0)  
                                            {cuda_gwypsquare (y,1);it=1;}
                                        else  if(bit != (lasterr_point-1)) 
                                            cuda_gwypsquare (y,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                                        else
                                            cuda_gwypsquare (y,2);
					gwyptogiant (y, tmp);
				}
			}
			if (*res)
			sprintf (buf, "%s is base %ld-Strong Fermat PRP! (%d decimal digits, %lu more test(s))", str, a, nbdg, iters);
		}
		else
			sprintf (buf, "%s is base %ld-Fermat PRP! (%d decimal digits)", str, a, nbdg);
	}

/* Print results.  Do not change the format of this line as Jim Fougeron of */
/* PFGW fame automates his QA scripts by parsing this line. */

	gwypfree (tmp);
	gwypfree (x);
	gwypfree (y);
	gwypfree (gwypminusone);
	gwypfree (gwypone);

/* Update the output file */

	if ((*res && IniGetInt (INI_FILE, (char*)"OutputPrimes", 0)) ||
	    (!*res && IniGetInt (INI_FILE, (char*)"OutputComposites", 0)))
		writeResults (buf);


#if defined(WIN32) && !defined(_CONSOLE)

	sprintf (buf+strlen(buf), "  Time : "); 
	ReplaceableLine (2);	/* Replace line */ 

#else

	clearline(100);

#ifdef _CONSOLE
	OutputBoth(buf);
#else
	if (*res) {
		OutputStr((char*)"\033[7m");
		OutputBoth(buf);
		OutputStr((char*)"\033[0m");
	}
	else
		OutputBoth(buf);
#endif

	sprintf (buf, "  Time : "); 

#endif

/* Output the final timings */

	gwypend_timer (1);
	gwypwrite_timer (buf+strlen(buf), 1, TIMER_CLR | TIMER_NL); 
	if ((*res && IniGetInt (INI_FILE, (char*)"OutputPrimes", 0)) ||
	    (!*res && IniGetInt (INI_FILE, (char*)"OutputComposites", 0)))
		OutputStr (buf);
	else
		OutputBoth (buf);

/* Cleanup and return */

//	iaddg (1, N);					// Restore the value of N
	Nlen = bitlen (N);
//	gwypdone ();
	_unlink (filename);
	lasterr_point = 0;
	return (TRUE);

/* An error occured, sleep, then try restarting at last save point. */

error:
	iaddg (1, N);					// Restore the value of N
	Nlen = bitlen (N);
	gwypfree (tmp);
	gwypfree (x);
	gwypfree (y);
	gwypfree (gwypone);
	gwypfree (gwypminusone);
//	gwypdone ();
	*res = FALSE;					// To avoid credit mesage...

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
		sprintf (buf, ERRMSG5, checknumber, str);
		OutputBoth (buf);
//		gwypdone ();
		will_try_larger_fft = FALSE;
		_unlink (filename);
		return (TRUE);
	}

//	gwypdone ();

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
		will_try_larger_fft = FALSE;
		return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
		OutputBoth (ERRMSG8);
		IniWriteInt(INI_FILE, (char*)"FFT_Increment", IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
		_unlink (filename);
		will_try_larger_fft = FALSE;
	}
	return (-1);
}

/* Test for 2*N+1 prime, knowing that N is prime -- gwypsetup has already been called. */

int commonCC1P (
	unsigned long a,
	int	*res, char *str)
{
	unsigned long bit, iters, nreduced, gcdn;
	gwypnum	x;
	giant	exponent, tmp;
	char	filename[20], buf[sgkbufsize+256], fft_desc[256], oldres64[17];
	long	write_time = DISK_WRITE_TIME * 60;
	int	echk, saving, stopping;
	time_t	start_time, current_time;
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;

/* First, test if gcd ((a^2-1), N) is one... */

	nreduced = gmodi (a*a-1, N);
	if (!nreduced)
		gcdn = a*a-1;
	else
		gcdn = gcd (a*a-1, nreduced);
	if (gcdn != 1) {
		sprintf (buf, "%s has a small factor : %lu!!\n", str, gcdn);
		OutputBoth (buf);
//		gwypdone ();
		*res = FALSE;
		return (TRUE);
	}


	exponent = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);
	tmp = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);
        
/* Init, subtract 1 from N to compute a^(N-1) mod N */

        gtog (N, exponent);
	iaddg (-1, exponent);
	Nlen = bitlen (N);
	nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.
	*res = TRUE;		/* Assume it is a prime */

/* Init filename */

	tempFileName (filename, 'z', N);

/* Allocate memory */

	x = gwypalloc ();

/* Optionally resume from save file and output a message */
/* indicating we are resuming a test */

	if (fileExists (filename) && readFromFile (filename, &bit, x, NULL)) {
		char	fmt_mask[80];
		double	pct;
		pct = trunc_percent (bit * 100.0 / Nlen);
		sprintf (fmt_mask,
			 "Resuming Pocklington prime test of %%s at bit %%ld [%%.%df%%%%]\n",
			 PRECISION);
		sprintf (buf, fmt_mask, str, bit, pct);
		OutputStr (buf);
		if (verbose)
			writeResults (buf);
	}

/* Otherwise, output a message indicating we are starting test */

	else {
		gwypclear_timers ();	// Init. timers
		if (setuponly) {
			if (FFTLEN != OLDFFTLEN) {
				OutputBoth (str); 
				OutputBoth ((char*)" : "); 
			}
		}
		else {
			if (showdigits)
				sprintf (buf, "Starting Pocklington prime test of %s (%d decimal digits)\n", str, nbdg);
			else
				sprintf (buf, "Starting Pocklington prime test of %s\n", str);
			OutputStr (buf);
		if (verbose)
			writeResults (buf);
		}
		bit = 1;
		itogwyp (a, x);
	}

/* Get the current time */

	gwypstart_timer (0);
	gwypstart_timer (1);
	time (&start_time);

/* Output a message about the FFT length */

	gwypfft_description (fft_desc);
#ifdef WIN32
	sprintf (buf, "%s, a = %lu\n", fft_desc, a);
#else
	sprintf (buf, "%s, a = %lu", fft_desc, a);
#endif
	OutputStr (buf);
	LineFeed();
	if (verbose) {
#if !defined(WIN32) 
		strcat (buf, "\n");
#endif
		writeResults (buf);
	}
	ReplaceableLine (1);	/* Remember where replaceable line is */

/* Init the title */

	title ((char*)"Pocklington prime test in progress...");

/* Do the PRP test */

	gwypsetmulbyconst (a);

	iters = 0;
	while (bit < Nlen) {

/* Error check the first and last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		stopping = stopCheck ();
		echk = stopping || ERRCHK || (bit <= 50) || (bit >= Nlen-50);
		if (((bit & 127) == 0) || (bit == 1) || (bit == (lasterr_point-1))) {
			echk = 1;
			time (&current_time);
			saving = ((current_time - start_time > write_time) || (bit == 1) || (bit == (lasterr_point-1)));
		} else
			saving = 0;

/* Process this bit */


		if (bitval (exponent, Nlen-bit-1)) {
			gwypsetnormroutine (0, echk, 1);
		} else {
			gwypsetnormroutine (0, echk, 0);
		}
		if (/*(bit+25 < Nlen) && (bit > 25) && */((bit != lasterr_point) || !maxerr_recovery_mode[6])) {
                    if (cufftonly)
                        gwypsquare (x);
                    else if (zp || generic)
                        cuda_gwypsquare (x, 3);
                    else if(bit==1 || it==0)  
                        {cuda_gwypsquare (x,1);it=1;}
                    else  if(bit != (lasterr_point-1)) 
                        cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                    else
                        cuda_gwypsquare (x,2);
                }
		else {
			gwypsquare_carefully (x);
			will_try_larger_fft = TRUE;
			if (bit == lasterr_point)
				maxerr_recovery_mode[6] = FALSE;
		}

		CHECK_IF_ANY_ERROR (x, (bit), Nlen, 6);

/* That iteration succeeded, bump counters */

		if (will_try_larger_fft && (bit == lasterr_point))
			saving = 1;					// Be sure to restart after this recovery iteration!
		will_try_larger_fft = FALSE;
		bit++;
		iters++;

/* Print a message every so often */

		if (bit % ITER_OUTPUT == 0) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (bit * 100.0 / Nlen);
			if (strlen (str) < 40) {
				sprintf (fmt_mask, "%%.%df%%%% of %%s", PRECISION);
				sprintf (buf, fmt_mask, pct, str);
			}
			else {
				sprintf (fmt_mask, "%%.%df%%%% of %%ld", PRECISION);
				sprintf (buf, fmt_mask, pct, Nlen);
			}
			title (buf);
			ReplaceableLine (2);	/* Replace line */
			sprintf (fmt_mask,
				 "%%s, bit: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, str, bit, Nlen, pct);
			OutputStr (buf);
			if (ERRCHK && bit > 30) {
				OutputStr ((char*)".  Round off: ");
				sprintf (buf, "%10.10f", reallyminerr);
				OutputStr (buf);
				sprintf (buf, " to %10.10f", reallymaxerr);
				OutputStr (buf);
			}
			gwypend_timer (0);
			if (CUMULATIVE_TIMING) {
				OutputStr ((char*)".  Time thusfar: ");
			} else {
				OutputStr ((char*)".  Time per bit: ");
				gwypdivide_timer (0, iters);
				iters = 0;
			}
			gwypprint_timer (0, TIMER_NL | TIMER_OPT_CLR);
			gwypstart_timer (0);
		}

/* Print a results file message every so often */

		if (bit % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
			sprintf (buf, "Bit %ld / %ld\n", bit, Nlen);
			writeResults (buf);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving || stopping) {
			write_time = DISK_WRITE_TIME * 60;
			if (! writeToFile (filename, bit, x, NULL)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (buf);
				if (write_time > 600) write_time = 600;
			}
			time (&start_time);

/* If an escape key was hit, write out the results and return */

			if (stopping) {
				gwypfree (tmp);
				gwypfree (x);
//				gwypdone ();
				*res = FALSE;
				return (FALSE);
			}
		}

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next iteration so that we can compare our */
/* residues to programs that start counter at zero or one. */

		if (interimResidues && bit % interimResidues < 2) {
				gwyptogiant (x, tmp);			// The modulo reduction is done here
			if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
			else if (abs(tmp->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
			else if (abs(tmp->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
			sprintf (buf, "%s interim residue %s at bit %ld\n", str, res64, bit);
			OutputBoth (buf);
		}

/* Write a save file every "interimFiles" iterations. */

		if (interimFiles && bit % interimFiles == 0) {
			char	interimfile[20];
			sprintf (interimfile, "%.8s.%03lu",
				 filename, bit / interimFiles);
			if (! writeToFile (interimfile, bit, x, NULL)) {
				sprintf (buf, WRITEFILEERR, interimfile);
				OutputBoth (buf);
			}
		}
	}

/* See if we've found a prime.  If not, format a 64-bit residue. */
/* Old versions of PRP used a non-standard 64-bit residue, computing */
/* 3^N-3 mod N rather than the more standard 3^(N-1) mod N.  Since */
/* some projects recorded these non-standard residues, output that */
/* residue too.  Note that some really old versions printed out the */
/* 32-bit chunks of the non-standard residue in reverse order. */

	clearline (100);

	gwyptogiant (x, tmp);
	if (!isone (tmp)) {
		*res = FALSE;	/* Not a prime */
		if (abs(tmp->sign) < 2)	// make a 64 bit residue correct !!
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
		else if (abs(tmp->sign) < 3)
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
		else if (abs(tmp->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
		else
			sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
		smulg (3, tmp); modg (N, tmp); iaddg (-3, tmp);
		if (abs(tmp->sign) < 2)	// make a 64 bit residue correct !!
			sprintf (oldres64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
		else if (abs(tmp->sign) < 3)
			sprintf (oldres64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
		else if (abs(tmp->sign) < 4)
			sprintf (oldres64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
		else
			sprintf (oldres64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
	}

/* Print results.  Do not change the format of this line as Jim Fougeron of */
/* PFGW fame automates his QA scripts by parsing this line. */

	if (*res)
		sprintf (buf, "%s is prime! (%d decimal digits", str, nbdg);
	else if (IniGetInt (INI_FILE, (char*)"OldRes64", 0))
		sprintf (buf, "%s is not prime.  RES64: %s.  OLD64: %s", str, res64, oldres64);
	else
		sprintf (buf, "%s is not prime.  RES64: %s", str, res64);

#if defined(WIN32) && !defined(_CONSOLE)

	sprintf (buf+strlen(buf), "  Time : "); 
	ReplaceableLine (2);	/* Replace line */ 

#else

	clearline(100);

#ifdef _CONSOLE
	OutputBoth(buf);
#else
	if (*res) {
		OutputStr((char*)"\033[7m");
		OutputBoth(buf);
		OutputStr((char*)"\033[0m");
	}
	else
		OutputBoth(buf);
#endif

	sprintf (buf, "  Time : "); 

#endif


/* Output the final timings */

	gwypend_timer (1);
	gwypwrite_timer (buf+strlen(buf), 1, TIMER_CLR | TIMER_NL); 
	OutputBoth (buf);

/* Cleanup and return */

        gwypfree (exponent);
        gwypfree (tmp);
	gwypfree (x);
//	gwypdone ();
	_unlink (filename);
	lasterr_point = 0;
	return (TRUE);

/* An error occured, sleep, then try restarting at last save point. */

error:
        gwypfree (exponent);
	gwypfree (tmp);
	gwypfree (x);
//	gwypdone ();
	*res = FALSE;					// To avoid credit mesage...

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
		sprintf (buf, ERRMSG5, checknumber, str);
		OutputBoth (buf);
//		gwypdone ();
		will_try_larger_fft = FALSE;
		_unlink (filename);
		return (TRUE);
	}

//	gwypdone ();

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
		will_try_larger_fft = FALSE;
		return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
		OutputBoth (ERRMSG8);
		IniWriteInt(INI_FILE, (char*)"FFT_Increment", IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
		_unlink (filename);
		will_try_larger_fft = FALSE;
	}
	return (-1);
}

/* Test for 2*N-1 prime, knowing that N is prime -- gwypsetup has already been called. */

int commonCC2P (
	unsigned long P,
	int *res,
	char *str)
{
	char	filename[20], buf[sgkbufsize+256], fft_desc[256]; 
//	unsigned long bits, explen, iters, bit, bitv, D, mask=0x80000000, frestart=FALSE;
	unsigned long explen, iters, bit, D;
	unsigned long nreduced, gcdn;
	long	write_time = DISK_WRITE_TIME * 60;
	int	echk, saving, stopping;
	time_t	start_time, current_time;
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;
	giant exponent, tmp, tmp2, tmp3;
	gwypnum x, y, gwypinvD;

/* First, test if gcd (U(2), N) is one... */

	nreduced = gmodi (P, N);
	if (!nreduced)
		gcdn = P;
	else
		gcdn = gcd (P, nreduced);
	if (gcdn != 1) {
		sprintf (buf, "%s has a small factor : %lu!!\n", str, gcdn);
		OutputBoth (buf);
//		gwypdone ();
		*res = FALSE;
		return (TRUE);
	}


	x = gwypalloc ();					// allocate memory for the gwypnums
	y = gwypalloc ();
	gwypinvD = gwypalloc ();

//	bits = bitlen (N);

	exponent = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);	// Allocate memory for exponent
	tmp = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);	// Allocate memory for tmp
	tmp2 = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);		// Allocate memory for tmp2
	tmp3 = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);		// Allocate memory for tmp3

	gtog (N, exponent);
	iaddg (1, exponent);					// exponent = modulus + 1

	explen = bitlen (exponent);
        nbdg = gnbdg (N, 10);


/* Init filename */

	tempFileName (filename, 'L', N);

/* Optionally resume from save file and output a message */
/* indicating we are resuming a test */

	if (fileExists (filename) && readFromFile (filename, &bit, x, y)) {
		char	fmt_mask[80];
		double	pct;
		pct = trunc_percent (bit * 100.0 / explen);
		sprintf (fmt_mask,
			 "Resuming Lucas sequence at bit %%ld [%%.%df%%%%]\n",
			 PRECISION);
		sprintf (buf, fmt_mask, bit, pct);
		OutputStr (buf);
		if (verbose)
			writeResults (buf);
		D = P*P-4;

	}

/* Otherwise, output a message indicating we are starting test */

	else {
		gwypclear_timers ();	// Init. timers
		D = P*P-4;
		sprintf (buf, "Starting Morrison prime test of %s\n", str);
		OutputStr (buf);
		if (verbose)
			writeResults (buf);
		bit = 1;
		itogwyp (2, x);
		itogwyp (P, y);
//		dbltogw (2.0, x);
//		dbltogw ((double)P, y);

	}

	itog (D, tmp3);    // Compute the inverse of D modulo N
	invg (N, tmp3);
	gianttogwyp (tmp3, gwypinvD);  // Convert it to gwypnum

/* Get the current time */

	gwypstart_timer (0);
	gwypstart_timer (1);
	time (&start_time);

/* Output a message about the FFT length */

	gwypfft_description (fft_desc);
#ifdef WIN32
	sprintf (buf, "%s, P = %lu\n", fft_desc, P);
#else
	sprintf (buf, "%s, P = %lu", fft_desc, P);
#endif
	OutputStr (buf);
	LineFeed();
	if (verbose) {
#if !defined(WIN32) 
		strcat (buf, "\n");
#endif
		writeResults (buf);
	}

	title ((char*)"Morrison prime test in progress...");

	ReplaceableLine (1);	/* Remember where replaceable line is */

	iters = 0;
	gwypsetnormroutine (0, 1, 0);
	will_try_larger_fft = FALSE;

	while (bit <= explen) {

/* Error check the first and last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		stopping = stopCheck ();
		echk = stopping || ERRCHK || (bit <= 50) || (bit >= Nlen-50);
		if (((bit & 127) == 0) || (bit == 1) || (bit == (lasterr_point-1))) {
			echk = 1;
			time (&current_time);
			saving = ((current_time - start_time > write_time) || (bit == 1) || (bit == (lasterr_point-1)));
		} else
			saving = 0;

/* Process this bit */

		gwypsetnormroutine (0, echk, 0);


		if (bitval (exponent, explen-bit)) {
			gwypsetaddin (-(int)P);
			if (/*(bit+26 < explen) && (bit > 26) && */
				((bit != lasterr_point) || (!maxerr_recovery_mode[1] && !maxerr_recovery_mode[2]))) {
                            if (cufftonly)
                                gwypmul (y, x);
                            else
                                cuda_gwypmul (y, x, 3);
                            will_try_larger_fft = FALSE;
			}
			else {
                            gwypmul_carefully (y, x);
                            if (bit == lasterr_point)
                                    maxerr_recovery_mode[1] = FALSE;
                            will_try_larger_fft = TRUE;
			}
			CHECK_IF_ANY_ERROR(x, (bit), explen, 1)
			gwypsetaddin (-2);
			if (/*(bit+26 < explen) && (bit > 26) && */
				((bit != lasterr_point) || !maxerr_recovery_mode[2])) {
                            if (cufftonly)
                                gwypsquare (y);
                            else if (zp || generic)
                                cuda_gwypsquare (y, 3);
                            else if(bit==1 || it==0)  
                                {cuda_gwypsquare (y,1);it=1;}
                            else  if(bit != (lasterr_point-1)) 
                                cuda_gwypsquare (y,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (y,2);
                            will_try_larger_fft = FALSE;
			}
			else {
				gwypsquare_carefully (y);
				if (bit == lasterr_point)
					maxerr_recovery_mode[2] = FALSE;
				will_try_larger_fft = TRUE;
			}
			CHECK_IF_ANY_ERROR(y, (bit), explen, 2)
		}
		else {
			gwypsetaddin (-(int)P);
			if (/*(bit+26 < explen) && (bit > 26) && */
				((bit != lasterr_point) || (!maxerr_recovery_mode[3] && !maxerr_recovery_mode[4]))) {
                            if (cufftonly)
                                gwypmul (x, y);
                            else
                                cuda_gwypmul (x, y, 3);
                            will_try_larger_fft = FALSE;
			}
			else {
                            gwypmul_carefully (x, y);
                            if (bit == lasterr_point)
                                    maxerr_recovery_mode[3] = FALSE;
                            will_try_larger_fft = TRUE;
			}
			CHECK_IF_ANY_ERROR(y, (bit), explen, 3)
			gwypsetaddin (-2);
			if (/*(bit+26 < explen) && (bit > 26) && */
				((bit != lasterr_point) || !maxerr_recovery_mode[4])) {
                            if (cufftonly)
                                gwypsquare (x);
                            else if (zp || generic)
                                cuda_gwypsquare (x, 3);
                            else if(bit==1 || it==0)  
                                {cuda_gwypsquare (x,1);it=1;}
                            else  if(bit != (lasterr_point-1)) 
                                cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (x,2);
                            will_try_larger_fft = FALSE;
			}
			else {
				gwypsquare_carefully (x);
				if (bit == lasterr_point)
					maxerr_recovery_mode[4] = FALSE;
				will_try_larger_fft = TRUE;
			}
			CHECK_IF_ANY_ERROR(x, (bit), explen, 4)
		}

 /* That iteration succeeded, bump counters */

		if (will_try_larger_fft && (bit == lasterr_point))
			saving = 1;					// Be sure to restart after this recovery iteration!
		will_try_larger_fft = FALSE;
		bit++;
		iters++;

/* Print a message every so often */

		if (bit % ITER_OUTPUT == 0) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (bit * 100.0 / explen);
			sprintf (fmt_mask, "%%.%df%%%% of %%ld", PRECISION);
			sprintf (buf, fmt_mask, pct, Nlen);
			title (buf);
			ReplaceableLine (2);	/* Replace line */
			sprintf (fmt_mask,
				 "%%s, bit: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, str, bit, explen, pct);
			OutputStr (buf);
			if (ERRCHK && bit > 30) {
				OutputStr ((char*)".  Round off: ");
				sprintf (buf, "%10.10f", reallyminerr);
				OutputStr (buf);
				sprintf (buf, " to %10.10f", reallymaxerr);
				OutputStr (buf);
			}
			gwypend_timer (0);
			if (CUMULATIVE_TIMING) {
				OutputStr ((char*)".  Time thusfar: ");
			} else {
				OutputStr ((char*)".  Time per bit: ");
				gwypdivide_timer (0, iters);
				iters = 0;
			}
			gwypprint_timer (0, TIMER_NL | TIMER_OPT_CLR);
			gwypstart_timer (0);
		}

/* Print a results file message every so often */

		if (bit % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
			sprintf (buf, "Bit %ld / %ld\n", bit, explen);
			writeResults (buf);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving || stopping) {
			write_time = DISK_WRITE_TIME * 60;
			if (! writeToFile (filename, bit, x, y)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (buf);
				if (write_time > 600) write_time = 600;
			}	
			time (&start_time);

/* If an escape key was hit, write out the results and return */

			if (stopping) {
				gwypfree (tmp);
				gwypfree (tmp2);
				gwypfree (tmp3);
				gwypfree (x);				// Clean up
				gwypfree (y);
				gwypfree (gwypinvD);
//				gwypdone ();
				*res = FALSE;
				return (FALSE);
			}
		}

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next iteration so that we can compare our */
/* residues to programs that start counter at zero or one. */

		if (interimResidues && bit % interimResidues < 2) {
			gwyptogiant (x, tmp);		// The modulo reduction is done here
			modg (N, tmp);				// External modulus and gwypnums one may be different...
			if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
			else if (abs(tmp->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
			else if (abs(tmp->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
			sprintf (buf, "%s interim residue %s at bit %ld\n", str, res64, bit);
			OutputBoth (buf);
		}

/* Write a save file every "interimFiles" iterations. */

		if (interimFiles && bit % interimFiles == 0) {
			char	interimfile[20];
			sprintf (interimfile, "%.8s.%03lu",
				 filename, bit / interimFiles);
			if (! writeToFile (interimfile, bit, x, y)) {
				sprintf (buf, WRITEFILEERR, interimfile);
				OutputBoth (buf);
			}
		}
	}

	clearline (100);

	will_try_larger_fft = TRUE;

	gwypsetaddin (0);			// Reset addin constant.
	gwypsetnormroutine (0, 1, 1);	// set mul. by const.
	gwypsetmulbyconst (2);
//	gwypmul (gwypinvD, y);			// y = D^-1*2*V(N+2) modulo N
        if (cufftonly)
            gwypmul (gwypinvD, y);
        else
            cuda_gwypmul (gwypinvD, y, 3);
	CHECK_IF_ANY_ERROR(y, (explen), explen, 4)
	gwypsetmulbyconst (P);
//	gwypmul (gwypinvD, x);			// x = D^-1*P*V(N+1) modulo N
        if (cufftonly)
            gwypmul (gwypinvD, x);
        else
            cuda_gwypmul (gwypinvD, x, 3);
	CHECK_IF_ANY_ERROR(x, (explen), explen, 4)
	gwypsub (x, y);				// y = D^-1*(2*V(N+2)-P*V(N+1)) = U(N+1) modulo N
	gwypsetnormroutine (0, 1, 0);	// reset mul by const

	gwyptogiant (y, tmp);		// tmp = U(N+1) modulo N

	will_try_larger_fft = FALSE;

	if (!isZero (tmp)) {
		*res = FALSE;				// Not a prime.
		if (abs(tmp->sign) < 2)	// make a 64 bit residue correct !!
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
		else if (abs(tmp->sign) < 3)
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
		else if (abs(tmp->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
		else
			sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
		sprintf (buf, "%s is not prime. Lucas RES64: %s", str, res64);
	}
	else
		sprintf (buf, "%s is prime! (%d decimal digits)", str, nbdg);


#if defined(WIN32) && !defined(_CONSOLE)

	sprintf (buf+strlen(buf), "  Time : "); 
	ReplaceableLine (2);	/* Replace line */ 

#else

	clearline(100);

#ifdef _CONSOLE
	OutputBoth(buf);
#else
	if (*res) {
		OutputStr((char*)"\033[7m");
		OutputBoth(buf);
		OutputStr((char*)"\033[0m");
	}
	else
		OutputBoth(buf);
#endif

	sprintf (buf, "  Time : "); 

#endif

/* Output the final timings */

	gwypend_timer (1);
	gwypwrite_timer (buf+strlen(buf), 1, TIMER_CLR | TIMER_NL); 
	OutputBoth (buf);

	gwypfree (x);				// Clean up
	gwypfree (y);
	gwypfree (gwypinvD);
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (tmp3);
	gwypfree (exponent);
//	gwypdone ();
	_unlink (filename);
	lasterr_point = 0;
	return TRUE;

/* An error occured, sleep, then try restarting at last save point. */

error:

	gwypfree (x);				// Clean up
	gwypfree (y);
	gwypfree (gwypinvD);
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (tmp3);
	gwypfree (exponent);
//	gwypdone ();
	*res = FALSE;

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
		sprintf (buf, ERRMSG5, checknumber, str);
		OutputBoth (buf);
//		gwypdone ();
		will_try_larger_fft = FALSE;
		_unlink (filename);
		return (TRUE);
	}

//	gwypdone ();

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
		will_try_larger_fft = FALSE;
		return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
		OutputBoth (ERRMSG8);
		IniWriteInt(INI_FILE, (char*)"FFT_Increment", IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
		_unlink (filename);
		will_try_larger_fft = FALSE;
	}
	return (-1);
}

/* Test if k*2^n+c is a probable prime. */

int fastIsPRP (
	double	k,				/* k in k*b^n+c */
	unsigned long b,		/* b in k*b^n+c */
	unsigned long n,		/* n in k*b^n+c */
	signed long c,			/* c in k*b^n+c */
	char *str,
	int	*res)
{
	int	retval, a;


	a = IniGetInt (INI_FILE, (char*)"FBase", 3);

	do {
		gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
		gwypsetmaxmulbyconst (a);
		if (!setupok (gwypsetup (k, b, n, c, N), N, str, res)) {
			return TRUE;
		}

/* Do the PRP test */

		retval = commonPRP (a, res, str);
	} while (retval == -1);

/* Clean up and return */

	if (retval == TRUE)
		IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
	gwypdone ();
	return (retval);
}


/* Test if k*b^n+c is the next prime in a Cunningham chain of the first kind. */

int fastIsCC1P (
	double	k,				/* k in k*b^n+c */
	unsigned long b,		/* b in k*b^n+c */
	unsigned long n,		/* n in k*b^n+c */
	signed long c,			/* c in k*b^n+c */
	char *str,
	int	*res)
{
	int	retval, a;

	a = IniGetInt (INI_FILE, (char*)"FBase", 3);

	do {
		gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
		gwypsetmaxmulbyconst (a);
		if (!setupok (gwypsetup (k, b, n, c, N), N, str, res)) {
			return TRUE;
		}

/* Do the Pocklington test */

		retval = commonCC1P (a, res, str);
	} while (retval == -1);

/* Clean up and return */

	if (retval == TRUE)
		IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
	gwypdone ();
	return (retval);
}

/* Test if k*2^n+c is the next prime in a Cunningham chain of the second kind. */

int fastIsCC2P (
	double	k,				/* k in k*b^n+c */
	unsigned long b,		/* b in k*b^n+c */
	unsigned long n,		/* n in k*b^n+c */
	signed long c,			/* c in k*b^n+c */
	char *str,
	int	*res)
{
	char	buf[sgkbufsize+256]; 
	int	retval, P;

/* Setup the assembly code. */


//	Compute P for the Morrison test (Q = 1)

	P = genLucasBaseP (N, IniGetInt (INI_FILE, (char*)"PBase", 3));
	if (P < 0) {
		if (P == -1)
			sprintf (buf, "Cannot compute P to test %s...\nThis is surprising, please, let me know that!!\nMy E-mail is jpenne@free.fr\n", str);
		else
			sprintf (buf, "%s has a small factor : %d !!\n", str, abs(P));
		OutputBoth (buf);
		return (TRUE); 
	}

	do {
		gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
		gwypsetmaxmulbyconst(P);
		if (!setupok (gwypsetup (k, b, n, c, N), N, str, res)) {
			return TRUE;
		}

/* Do the Morrison test */

		retval = commonCC2P (P, res, str);
	} while (retval == -1);

/* Clean up and return */

	if (retval == TRUE)
		IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
	gwypdone ();
	return (retval);
}

/* Test if k*b^n+c (or a factor of it) is a Frobenius probable prime. */

int fastIsFrobeniusPRP (
	double	k,				/* k in k*b^n+c */
	unsigned long b,		/* b in k*b^n+c */
	unsigned long n,		/* n in k*b^n+c */
	signed long c,			/* c in k*b^n+c */
	char *str,
	int	*res)
{
	char	buf[sgkbufsize+256]; 
	int	retval;
	uint32_t P = 3, Q = 0;
	long D;

/* Setup the assembly code. */

	P = IniGetInt (INI_FILE, (char*)"PBase", 3);

	D = generalLucasBase (N, &P, &Q);
	if (D < 0) {
		if (D == -1)
			sprintf (buf, "Cannot compute D to test %s...\nThis is surprising, please, let me know that!!\nMy E-mail is jpenne@free.fr\n", str);
		else
			sprintf (buf, "%s has a small factor : %ld !!\n", str, abs(D));
		OutputBoth (buf);
		return (TRUE); 
	}

	do {
		gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
		gwypsetmaxmulbyconst (max (3, Q));
		if (!setupok (gwypsetup (k, b, n, c, N), N, str, res)) {
			return TRUE;
		}

/* Do the Frobenius PRP test */

		retval = commonFrobeniusPRP (P, Q, res, str);

	} while (retval == -1);

/* Clean up and return */
	
	if (retval == TRUE)
		IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
	gwypdone ();
	return (retval);
}

/* Test if N is a Frobenius probable prime.  The number N can be of ANY form. */

int slowIsFrobeniusPRP (
	char	*str,		/* string representation of N */
	int	*res)
{
	char	buf[sgkbufsize+256]; 
	int	retval;

	uint32_t P = 3, Q = 0;
	long D;

/* Setup the assembly code. */

	P = IniGetInt (INI_FILE, (char*)"PBase", 3);

	D = generalLucasBase (N, &P, &Q);
	if (D < 0) {
		if (D == -1)
			sprintf (buf, "Cannot compute D to test %s...\nThis is surprising, please, let me know that!!\nMy E-mail is jpenne@free.fr\n", str);
		else
			sprintf (buf, "%s has a small factor : %ld !!\n", str, abs(D));
		OutputBoth (buf);
		return (TRUE); 
	}

	do {
		gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
		gwypsetmaxmulbyconst (max (3, Q));
		if (!setupok (gwypsetup_general_mod_giant (N), N, str, res)) {
			return TRUE;
		}

/* Do the Frobenius PRP test */

		retval = commonFrobeniusPRP (P, Q, res, str);

	} while (retval == -1);

/* Clean up and return */

	if (retval == TRUE)
		IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
	gwypdone ();
	return (retval);
}
/* Test if N is a Wieferich prime.  The number N can be of ANY form. */

int slowIsWieferich (
	char	*str,		/* string representation of N */
	int	*res)
{
	int	a,retval;
	char	buf[sgkbufsize+256]; 

	M = newgiant ((bitlen (N) >> 3) + 8);

	gtog (N, M);
	squareg (M);

/* Setup the gwypnum code */

	a = IniGetInt (INI_FILE, (char*)"FBase", 2);

	do {
		gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
		gwypsetmaxmulbyconst (a);
		if (!setupok (gwypsetup_general_mod_giant (M), M, str, res)) {
			return TRUE;
		}

/* Do the divisibility test */

		retval = isexpdiv (a, N, M, res);

	} while (retval == -1);

	if (retval) {
		if (*res)
			sprintf (buf, "%s is a Base %d Wieferich prime!! (%d decimal digits)\n", str, a, nbdg);
		else
			sprintf (buf, "%s is not a Base %d Wieferich prime. RES64: %s\n", str, a, res64);
		OutputBoth (buf);
	}



/* Clean up and return */

	gwypfree (M);
	if (retval == TRUE)
		IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
	gwypdone ();

	return (retval);
}

void TestWieferich ()
{
	char str[10];
	int n, res;

	N = newgiant (2);

	for (n=3; n<10000; n+=2) {
		if (!isPrime (n))
			continue;
		itog (n, N);
		sprintf (str, "%d", n);
		slowIsWieferich (str, &res);
	}

	gwypfree (N);
}

/* Test if N is a probable prime.  The number N can be of ANY form. */

int slowIsPRP (
	char	*str,		/* string representation of N */
	int	*res)
{
	int	retval, a;

/* Setup the gwypnum code */

	a = IniGetInt (INI_FILE, (char*)"FBase", 3);

	do {
		gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
		gwypsetmaxmulbyconst (a);
		if (!setupok (gwypsetup_general_mod_giant (N), N, str, res)) {
			return TRUE;
		}

/* Do the PRP test */

		retval = commonPRP (a, res, str);
	} while (retval == -1);

/* Clean up and return */

	if (retval == TRUE)
		IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
	gwypdone ();
	return (retval);
}

/* Test if N is the next prime in a CC1 chain.  The number N can be of ANY form. */

int slowIsCC1P (
	char	*str,		/* string representation of N */
	int	*res)
{
	int	retval, a;

/* Setup the gwypnum code */


	a = IniGetInt (INI_FILE, (char*)"FBase", 3);

	do {
		gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
		gwypsetmaxmulbyconst (a);
		if (!setupok (gwypsetup_general_mod_giant (N), N, str, res)) {
			return TRUE;
		}

/* Do the Pocklington test */

		retval = commonCC1P (a, res, str);
	} while (retval == -1);

/* Clean up and return */

	if (retval == TRUE)
		IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
	gwypdone ();
	return (retval);
}

/* Test if N is the next prime in a CC2 chain.  The number N can be of ANY form. */

int slowIsCC2P (
	char	*str,				/* string representation of N */
	int	*res)
{
	char	buf[sgkbufsize+256]; 
	int	retval, P;

/* Setup the gwypnum code */


//	Compute P for the Morrison test (Q = 1)

	P = genLucasBaseP (N, IniGetInt (INI_FILE, (char*)"PBase", 3));
	if (P < 0) {
		if (P == -1)
			sprintf (buf, "Cannot compute P to test %s...\nThis is surprising, please, let me know that!!\nMy E-mail is jpenne@free.fr\n", str);
		else
			sprintf (buf, "%s has a small factor : %d !!\n", str, abs(P));
		OutputBoth (buf);
		return (TRUE); 
	}

	do {
		gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
		gwypsetmaxmulbyconst(P);
		if (!setupok (gwypsetup_general_mod_giant (N), N, str, res)) {
			return TRUE;
		}

/* Do the Morrison test */

		retval = commonCC2P (P, res, str);
	} while (retval == -1);

/* Clean up and return */

	if (retval == TRUE)
		IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
	gwypdone ();
	return (retval);
}

/* Test if a small N is a probable prime. */
/* Compute b^(N-1) mod N */

int isProbablePrime (void)
{
	int	retval;
	giant	x;

	if (isone (N)) return (FALSE);
	x = newgiant (N->sign + 8);
	itog (IniGetInt (INI_FILE, (char*)"FBase", 3), x);
	powermodg (x, N, N);
	iaddg (-IniGetInt (INI_FILE, (char*)"FBase", 3), x);
	retval = isZero (x);
	gwypfree (x);
	return (retval);
}

int isPRPinternal (
	char *str, double dk, 
	unsigned long base,
	unsigned long n,
	int incr,
	int *res)
{
//	J.P. shadow        char buf[100];
	char	filename[20], buf[sgkbufsize+256]; // Lei - not need such long char
	unsigned long retval, fcontinue = FALSE;

	tempFileName (filename, 'L', N);	  // See if resuming a Lucas or Frobenius PRP test
	fcontinue = fileExists (filename);
	tempFileName (filename, 'F', N);
	fcontinue = fcontinue || fileExists (filename);


	if (dk >= 1.0) {
		if (fcontinue || IniGetInt(INI_FILE, (char*)"LucasPRPtest", 0)) {
			if (!fcontinue)
				gwypclear_timers ();				// Init. timers
			retval = fastIsFrobeniusPRP (dk, base, n, incr, str, res);
		}
		else {
			retval = fastIsPRP (dk, base, n, incr, str, res);
			if (retval && *res && !Fermat_only && !IniGetInt(INI_FILE, (char*)"FermatPRPtest", 0))
				retval = fastIsFrobeniusPRP (dk, base, n, incr, str, res);
			else if (retval == TRUE)// If not stopped by user...
				IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
		}
	}
	else if (Nlen < 50) {
		*res = isProbablePrime();
		if (*res)
		{
#ifndef	WIN32
			OutputStr((char*)"\033[7m");
#endif
			sprintf (buf, "%s is a probable prime.\n", str);
			OutputBoth(buf);
#ifndef	WIN32
			OutputStr((char*)"\033[0m");
#endif
		}
		else {
			sprintf (buf, "%s is not prime.\n", str);
			OutputBoth (buf);
		}
		retval = TRUE;
	}
	else {
		if (fcontinue || IniGetInt(INI_FILE, (char*)"LucasPRPtest", 0)) {
			if (!fcontinue)
				gwypclear_timers ();				// Init. timers
			retval = slowIsFrobeniusPRP (str, res);
		}
		else {
			retval = slowIsPRP (str, res);
			if (retval && *res && !Fermat_only && !IniGetInt(INI_FILE, (char*)"FermatPRPtest", 0))
				retval = slowIsFrobeniusPRP (str, res);
			else if (retval == TRUE)// If not stopped by user...
				IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
		}
	}
	return retval;
}
int gisPRPinternal (
	char *str, 
        double dk, 
	giant gb,
	unsigned long n,
	int incr,
	int *res)
{
//	J.P. shadow        char buf[100];
	char	filename[20], buf[sgkbufsize+256]; // Lei - not need such long char
	unsigned long retval, fcontinue = FALSE, smallbase = 0;


	if (abs(gb->sign) <= 2)	{	// Test if the base is a small integer
            smallbase = gb->n[0];
            if (abs(gb->sign) == 2)
                smallbase += 65536*gb->n[1];
        }

        tempFileName (filename, 'L', N); // See if resuming a Lucas or Frobenius PRP test
	fcontinue = fileExists (filename);
	tempFileName (filename, 'F', N);
	fcontinue = fcontinue || fileExists (filename);


	if (dk >= 1.0) {
		if (fcontinue || IniGetInt(INI_FILE, (char*)"LucasPRPtest", 0)) {
			if (!fcontinue)
				gwypclear_timers ();				// Init. timers
			retval = fastIsFrobeniusPRP (dk, smallbase, n, incr, str, res);
		}
		else {
			retval = fastIsPRP (dk, smallbase, n, incr, str, res);
			if (retval && *res && !Fermat_only && !IniGetInt(INI_FILE, (char*)"FermatPRPtest", 0))
				retval = fastIsFrobeniusPRP (dk, smallbase, n, incr, str, res);
			else if (retval == TRUE)// If not stopped by user...
				IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
		}
	}
	else if (Nlen < 50) {
		*res = isProbablePrime();
		if (*res)
		{
#ifndef	WIN32
			OutputStr((char*)"\033[7m");
#endif
			sprintf (buf, "%s is a probable prime.\n", str);
			OutputBoth(buf);
#ifndef	WIN32
			OutputStr((char*)"\033[0m");
#endif
		}
		else {
			sprintf (buf, "%s is not prime.\n", str);
			OutputBoth (buf);
		}
		retval = TRUE;
	}
	else {
		if (fcontinue || IniGetInt(INI_FILE, (char*)"LucasPRPtest", 0)) {
			if (!fcontinue)
				gwypclear_timers ();				// Init. timers
			retval = slowIsFrobeniusPRP (str, res);
		}
		else {
			retval = slowIsPRP (str, res);
			if (retval && *res && !Fermat_only && !IniGetInt(INI_FILE, (char*)"FermatPRPtest", 0))
				retval = slowIsFrobeniusPRP (str, res);
			else if (retval == TRUE)// If not stopped by user...
				IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
		}
	}
	return retval;
}

#define NPG	0   // NEWPGEN output format, not AP mode
#define NPGAP	1   // NEWPGEN output format, AP mode
#define	ABCC	2   // Carol ABC format
#define	ABCK	3   // Kynea ABC format
#define ABCCW	4   // Cullen/Woodall ABC format
#define ABCFF	5   // FermFact output ABC format
#define ABCGM	6   // Gaussian Mersenne ABC format
#define ABCLEI  7   // Lei ABC format
#define ABCSP	8   // (2^n+1)/3 ABC format
#define NPGCC1	9   // Cunningham chain first kind
#define NPGCC2	10  // Cunningham chain second kind


// New ABC formats for k*b^n+c

#define ABCFKGS	11  // Fixed k:  b and n specified on each input line
#define ABCFKAS	12  // Fixed k:  b, n, and c specified on each input line
#define ABCFBGS	13  // Fixed b:  k and n specified on each input line
#define ABCFBAS	14  // Fixed b:  k, n, and c specified on each input line
#define ABCFNGS	15  // Fixed n:  k and b specified on each input line
#define ABCFNAS	16  // Fixed n:  k, b, and c specified on each input line
#define ABCVARGS    17	// k, b, and n specified on each input line
#define ABCVARAS    18	// k, b, n, and c specified on each input line
#define ABCVARAQS   19	// k, b, n, c and quotient specified on each input line
#define	ABCRU	20  // (10^n-1)/9 Repunits
#define	ABCGRU	21  // (b^n-1)/(b-1) Generalized Repunits

#define ABCGF	22  // ABC format for generalized Fermat numbers
#define ABCDN	23  // b^n-b^m+c format, m < n <= 2*m
#define ABCDNG	24  // General b^n-b^m+c format, m < n <= 2*m
#define ABCWFT	25  // Format used for Wieferich test
#define ABCWFS	26  // Format used for Wieferich search
#define ABCGPT	27  // Format used for General prime test (APRCL)

int IsPRP (							// General PRP test
	unsigned long format, 
	char *sgk,
	unsigned long base,
	unsigned long n, 
	int incr,
	unsigned long shift,
	int	*res) 
{  
	char	str[sgkbufsize+256], sgk1[sgkbufsize], buf[sgkbufsize+256]; 
	unsigned long bits, retval;
	double dk;
	giant gd, gr;

	if (format == ABCRU || format == ABCGRU) {	// Repunits or Generalized Repunits
		sprintf (str, "(%lu^%lu-1)/%lu", base, n, base-1);
		gk = newgiant (1);
		itog (1, gk);
	}
	else if (!(format == ABCC || format == ABCK)) {
		gk = newgiant (strlen(sgk)/2 + 8);	// Allocate one byte per decimal digit + spares
		ctog (sgk, gk);						// Convert k string to giant
		gshiftleft (shift, gk);				// Shift k multiplier if requested
		gtoc (gk, sgk1, sgkbufsize);		// Updated k string
		if (mask & MODE_DUAL) {
			sprintf (str, "%lu^%lu%c%d", base, n, incr < 0 ? '-' : '+', abs(incr));
		}
		else if (format != NPGAP) {			// Not MODE_AP
			if (!strcmp(sgk1, "1"))
				if (format == ABCVARAQS)
					sprintf (str, "(%lu^%lu%c%d)/%s", base, n, incr < 0 ? '-' : '+', abs(incr), sgd);
				else
					sprintf (str, "%lu^%lu%c%d", base, n, incr < 0 ? '-' : '+', abs(incr));
			else
				if (format == ABCVARAQS)
					sprintf (str, "(%s*%lu^%lu%c%d)/%s", sgk1, base, n, incr < 0 ? '-' : '+', abs(incr), sgd);
				else
					sprintf (str, "%s*%lu^%lu%c%d", sgk1, base, n, incr < 0 ? '-' : '+', abs(incr));
		}
		else {								// MODE_AP
			if (!strcmp(sgk1, "1"))
				sprintf (str, "%lu^%lu+1", base, n);
			else
				sprintf (str, "%lu^%lu+2*%s-1", base, n, sgk1);
		}
	}
	else {
		gk = newgiant ((n>>4)+8);
		itog (1, gk);						// Compute k multiplier
		gshiftleft (n-2, gk);				// Warning : here, n is exponent+1 !
		if (format == ABCK) {
			iaddg (1, gk);
			sprintf (str, "%s*2^%lu%c1 = (2^%lu+1)^2 - 2", sgk, n, '-', n-1);
		}
		else {
			iaddg (-1, gk);
			sprintf (str, "%s*2^%lu%c1 = (2^%lu-1)^2 - 2", sgk, n, '-', n-1);
		}
	}

	if ((gformat == ABCDN) || (gformat == ABCDNG)) {// Compute gk = gb^(n-m)-1
		bits = ndiff*log (base)/log (2);
		gk = newgiant ((bits >> 2) + 8);
		itog (base, gk);
		power (gk, ndiff);
		iaddg (-1, gk);
		sprintf (str, "%lu^%lu-%lu^%lu%c%d", base, n+ndiff, base, n, incr < 0 ? '-' : '+', abs(incr));
	}
	
	bits = (unsigned long) ((n * log(base)) / log(2) + bitlen(gk)); 
	N =  newgiant ((bits >> 2) + 8);   // Allocate memory for N

//	Compute the number we are testing.

	itog (base, N);
	power (N, n);

	if (format == NPGAP) {	// mode AP
		addg(gk, N);
		addg(gk, N);
	}
	else {			// not mode AP
		mulg (gk, N); 
	}

	iaddg (incr, N);

	if (format == ABCRU || format == ABCGRU) {
		if (!isPrime (n)) {
			sprintf (buf, "%s is not prime because %lu is not prime!\n", str, n);
			OutputBoth (buf);
			*res = FALSE;
			gwypfree (N);
			gwypfree (gk);
			return TRUE;
		}
		uldivg (base - 1, N);
		strong = FALSE;				// Do a simple Fermat PRP test (not strong).
	}
	else if (format == ABCVARAQS) {
		gd = newgiant (strlen(sgd)/2 + 8);	// Allocate one byte per decimal digit + spares
		gr = newgiant ((bits >> 4) + 8);	// Allocate memory for the remainder
		ctog (sgd, gd);						// Convert quotient string to giant
		gtog (N, gr);
		modg (gd, gr);
		if (!isZero(gr)) {
			sprintf (buf, "%s is not an integer!\n", str);
			OutputBoth (buf);
			*res = FALSE;
			gwypfree (gr);
			gwypfree (gd);
			gwypfree (N);
			gwypfree (gk);
			return TRUE;
		}
		else {
			divg (gd, N);
			strong = FALSE;			// Do a simple Fermat PRP test (not strong).
		}
	}

		Nlen = bitlen (N); 
		klen = bitlen(gk);
                nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.

		if (klen > 53) {					// we must use generic reduction
			dk = 0.0;
		}
		else {								// we can use DWT ; compute the multiplier as a double
		dk = (double)gk->n[0];
		if (gk->sign > 1)
			dk += 65536.0*(double)gk->n[1];
		if (gk->sign > 2)
			dk += 65536.0*65536.0*(double)gk->n[2];
		if (gk->sign > 3)
			dk += 65536.0*65536.0*65536.0*(double)gk->n[3];
		}
		
		if (setupok ((nbdg < 400), N, str, res)) // Force APRCL test for small numbers...
                    retval = isPRPinternal (str, dk, base, n, incr, res);
                else
                    retval = TRUE;
//	}

	strong = TRUE;		// Restore Strong Fermat PRP test

	if (format == ABCVARAQS) {
		gwypfree (gr);
		gwypfree (gd);
	}
	gwypfree (N);
	gwypfree (gk);
	return retval;
}

int gIsPRP (							// General PRP test
	unsigned long format, 
	char *sgk,
	char *sgb,
	giant gb,
	unsigned long n,
	int incr,
	unsigned long shift,
	int	*res) 
{  
	char	str[sgkbufsize+256], sgk1[sgkbufsize], buf[sgkbufsize+256]; 
	unsigned long bits, retval, smallbase = 0;
	double dk;
	giant gd, gr;

	if (abs(gb->sign) <= 2)	{	// Test if the base is a small integer
            smallbase = gb->n[0];
            if (abs(gb->sign) == 2)
                smallbase += 65536*gb->n[1];
        }

	if (format == ABCRU || format == ABCGRU) {	// Repunits or Generalized Repunits
		if (smallbase)
			sprintf (str, "(%lu^%lu-1)/%lu", smallbase, n, smallbase-1);
		else
			sprintf (str, "(%s^%lu-1)/(%s-1)", sgb, n, sgb);
		gk = newgiant (1);
		itog (1, gk);
	}
	else if (!(format == ABCC || format == ABCK)) {
		gk = newgiant (strlen(sgk)/2 + 8);	// Allocate one byte per decimal digit + spares
		ctog (sgk, gk);						// Convert k string to giant
		gshiftleft (shift, gk);				// Shift k multiplier if requested
		gtoc (gk, sgk1, sgkbufsize);		// Updated k string
		if (mask & MODE_DUAL) {
			sprintf (str, "%s^%lu%c%d", sgb, n, incr < 0 ? '-' : '+', abs(incr));
		}
		else if (format != NPGAP) {			// Not MODE_AP
			if (!strcmp(sgk1, "1"))
				if (format == ABCVARAQS)
					sprintf (str, "(%s^%lu%c%d)/%s", sgb, n, incr < 0 ? '-' : '+', abs(incr), sgd);
				else
					sprintf (str, "%s^%lu%c%d", sgb, n, incr < 0 ? '-' : '+', abs(incr));
			else
				if (format == ABCVARAQS)
					sprintf (str, "(%s*%s^%lu%c%d)/%s", sgk1, sgb, n, incr < 0 ? '-' : '+', abs(incr), sgd);
				else
					if ((n != 0) || (incr != 0))
						sprintf (str, "%s*%s^%lu%c%d", sgk1, sgb, n, incr < 0 ? '-' : '+', abs(incr));
					else
						sprintf (str, "%s", sgk1);
		}
		else {								// MODE_AP
			if (!strcmp(sgk1, "1"))
				sprintf (str, "%s^%lu+1", sgb, n);
			else
				sprintf (str, "%s^%lu+2*%s-1", sgb, n, sgk1);
		}
	}
	else {
		gk = newgiant ((n>>4)+8);
		itog (1, gk);						// Compute k multiplier
		gshiftleft (n-2, gk);				// Warning : here, n is exponent+1 !
		if (format == ABCK) {
			iaddg (1, gk);
			sprintf (str, "%s*2^%lu%c1 = (2^%lu+1)^2 - 2", sgk, n, '-', n-1);
		}
		else {
			iaddg (-1, gk);
			sprintf (str, "%s*2^%lu%c1 = (2^%lu-1)^2 - 2", sgk, n, '-', n-1);
		}
	}

	if ((gformat == ABCDN) || (gformat == ABCDNG)) {// Compute gk = gb^(n-m)-1
		bits = ndiff*bitlen (gb);
		gk = newgiant ((bits >> 2) + 8);
		gtog (gb, gk);
		power (gk, ndiff);
		iaddg (-1, gk);
		sprintf (str, "%s^%lu-%s^%lu%c%d", sgb, n+ndiff, sgb, n, incr < 0 ? '-' : '+', abs(incr));
	}

	if (smallbase)
		bits = (unsigned long) ((n * log((double) smallbase)) / log(2.0) + bitlen(gk));
	else
		bits = n * bitlen(gb) + bitlen(gk); 
	N =  newgiant ((bits >> 2) + 8);   // Allocate memory for N

//	Compute the number we are testing.

	gtog (gb, N);
	power (N, n);

	if (format == NPGAP) {	// mode AP
		addg(gk, N);
		addg(gk, N);
	}
	else {			// not mode AP
		mulg (gk, N); 
	}

	iaddg (incr, N);

	if (format == ABCRU || format == ABCGRU) {
		if (!isPrime (n)) {
			sprintf (buf, "%s is not prime because %lu is not prime!\n", str, n);
			OutputBoth (buf);
			*res = FALSE;
			gwypfree (N);
			gwypfree (gk);
			return TRUE;
		}
		iaddg (-1, gb);
		divg (gb, N);				// Divide N by (base-1)
		iaddg (1, gb);
		quotient = TRUE;
//		strong = FALSE;				// Do a simple Fermat PRP test (not strong).
	}
	else if (format == ABCVARAQS) {
		gd = newgiant (strlen(sgd)/2 + 8);	// Allocate one byte per decimal digit + spares
		gr = newgiant ((bits >> 4) + 8);	// Allocate memory for the remainder
		ctog (sgd, gd);						// Convert divisor string to giant
		gtog (N, gr);
		modg (gd, gr);
		if (!isZero(gr)) {
			sprintf (buf, "%s is not an integer!\n", str);
			OutputBoth (buf);
			*res = FALSE;
			free (gr);
			free (gd);
			free (N);
			free (gk);
			return TRUE;
		}
		else {
			divg (gd, N);
                        quotient = TRUE;
//			strong = FALSE;			// Do a simple Fermat PRP test (not strong).
		}
	}

		Nlen = bitlen (N); 
		klen = bitlen(gk);
                nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.

		if (klen > 53) {					// we must use generic reduction
			dk = 0.0;
		}
		else {								// we can use DWT ; compute the multiplier as a double
		dk = (double)gk->n[0];
		if (gk->sign > 1)
			dk += 65536.0*(double)gk->n[1];
		if (gk->sign > 2)
			dk += 65536.0*65536.0*(double)gk->n[2];
		if (gk->sign > 3)
			dk += 65536.0*65536.0*65536.0*(double)gk->n[3];
		}
		
		if (setupok ((nbdg < 400), N, str, res)) // Force APRCL test for small numbers...
                    retval = gisPRPinternal (str, dk, gb, n, incr, res);
                else
                    retval = TRUE;
//	}

	strong = TRUE;		// Restore Strong Fermat PRP test

	if (format == ABCVARAQS) {
		gwypfree (gr);
		gwypfree (gd);
	}
	quotient = FALSE;
	gwypfree (N);
	gwypfree (gk);
	return retval;
}

int IsCCP (	// General test for the next prime in a Cunningham chain
	unsigned long format, 
	char *sgk,
	unsigned long base,
	unsigned long n, 
	int incr,
	unsigned long shift,
	int	*res) 
{  
	char	str[sgkbufsize+256], sgk1[sgkbufsize]; 
	unsigned long bits, retval;
	double dk;



	gk = newgiant (strlen(sgk)/2 + 8);	// Allocate one byte per decimal digit + spares
	ctog (sgk, gk);						// Convert k string to giant
	gshiftleft (shift, gk);				// Shift k multiplier if requested
	gtoc (gk, sgk1, sgkbufsize);		// Updated k string
		if (mask & MODE_DUAL) {
			sprintf (str, "%lu^%lu%c%d", base, n, incr < 0 ? '-' : '+', abs(incr));
		}
		else
			if (!strcmp(sgk1, "1"))
				sprintf (str, "%lu^%lu%c%d", base, n, incr < 0 ? '-' : '+', abs(incr));
			else
				sprintf (str, "%s*%lu^%lu%c%d", sgk1, base, n, incr < 0 ? '-' : '+', abs(incr));

	bits = (unsigned long) ((n * log(base)) / log(2) + bitlen(gk)); 
	N =  newgiant ((bits >> 2) + 8);   // Allocate memory for N

//	Compute the number we are testing.

	itog (base, N);
	power (N, n);
	mulg (gk, N);
	iaddg (incr, N);

		klen = bitlen(gk);
                nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.

		if (klen > 53) {					// we must use generic reduction
			dk = 0.0;
		}
		else {								// we can use DWT ; compute the multiplier as a double
		dk = (double)gk->n[0];
		if (gk->sign > 1)
			dk += 65536.0*(double)gk->n[1];
		if (gk->sign > 2)
			dk += 65536.0*65536.0*(double)gk->n[2];
		if (gk->sign > 3)
			dk += 65536.0*65536.0*65536.0*(double)gk->n[3];
		}

		if (!setupok ((nbdg < 400), N, str, res)) // Force APRCL test for small numbers...
                    retval = TRUE;
 		else if (dk >= 1.0)
			if (format == NPGCC1)
				retval = fastIsCC1P (dk, base, n, incr, str, res);
			else if (format == NPGCC2)
				retval = fastIsCC2P (dk, base, n, incr, str, res);
			else
				retval = FALSE;
		else
			if (format == NPGCC1)
				retval = slowIsCC1P (str, res);
			else if (format == NPGCC2)
				retval = slowIsCC2P (str, res);
			else
				retval = FALSE;

//	}
	gwypfree (N);
	gwypfree (gk);
	return retval;
}

int gIsCCP (	// General test for the next prime in a Cunningham chain
	unsigned long format, 
	char *sgk,
	char *sgb,
	giant gb,
	unsigned long n, 
	int incr,
	unsigned long shift,
	int	*res) 
{  
	char	str[sgkbufsize+256], sgk1[sgkbufsize]; 
	unsigned long bits, retval, smallbase = 0;
	double dk;



	if (abs(gb->sign) <= 2)	{	// Test if the base is a small integer
            smallbase = gb->n[0];
            if (abs(gb->sign) == 2)
                smallbase += 65536*gb->n[1];
        }

	gk = newgiant (strlen(sgk)/2 + 8);	// Allocate one byte per decimal digit + spares
	ctog (sgk, gk);						// Convert k string to giant
	gshiftleft (shift, gk);				// Shift k multiplier if requested
	gtoc (gk, sgk1, sgkbufsize);		// Updated k string
		if (mask & MODE_DUAL) {
			sprintf (str, "%s^%lu%c%d", sgb, n, incr < 0 ? '-' : '+', abs(incr));
		}
		else
			if (!strcmp(sgk1, "1"))
				sprintf (str, "%s^%lu%c%d", sgb, n, incr < 0 ? '-' : '+', abs(incr));
			else
				sprintf (str, "%s*%s^%lu%c%d", sgk1, sgb, n, incr < 0 ? '-' : '+', abs(incr));

	if (smallbase)	
		bits = (unsigned long) ((n * log((double) smallbase)) / log(2.0) + bitlen(gk));
	else
		bits = n * bitlen(gb) + bitlen(gk);
	N =  newgiant ((bits >> 2) + 8);   // Allocate memory for N

//	Compute the number we are testing.

	gtog (gb, N);
	power (N, n);
	mulg (gk, N);
	iaddg (incr, N);

		klen = bitlen(gk);
                nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.

		if (klen > 53 || generic || !smallbase) {					// we must use generic reduction
			dk = 0.0;
		}
		else {								// we can use DWT ; compute the multiplier as a double
		dk = (double)gk->n[0];
		if (gk->sign > 1)
			dk += 65536.0*(double)gk->n[1];
		if (gk->sign > 2)
			dk += 65536.0*65536.0*(double)gk->n[2];
		if (gk->sign > 3)
			dk += 65536.0*65536.0*65536.0*(double)gk->n[3];
		}

		if (!setupok ((nbdg < 400), N, str, res)) // Force APRCL test for small numbers...
                    retval = TRUE;
 		else if (dk >= 1.0)
			if (format == NPGCC1)
				retval = fastIsCC1P (dk, smallbase, n, incr, str, res);
			else if (format == NPGCC2)
				retval = fastIsCC2P (dk, smallbase, n, incr, str, res);
			else
				retval = FALSE;
		else
			if (format == NPGCC1)
				retval = slowIsCC1P (str, res);
			else if (format == NPGCC2)
				retval = slowIsCC2P (str, res);
			else
				retval = FALSE;

//	}
	gwypfree (N);
	gwypfree (gk);
	return retval;
}

unsigned long gcd (
	unsigned long x,
	unsigned long y)
{
	unsigned long w;

	if (!x || x==y)
		return y;
	if (!y)
		return x;

	if (x < y) {
		w = x;
		x = y;
		y = w;
	}

	while (1) {
		x = x%y;
		if (!x)
			break;
		w = x;
		x = y;
		y = w;
	}

	return (y);
}

void findbpf (unsigned long base) {// find all prime factors of the base
	unsigned long b, p;
	int i;

	for (i=0; i<30; i++)
            bpf[i] = bpc[i] = vpf[i] = 0;	// clean up

// bpf[i] : base prime factor, vpf[i] : its exponent, bpc[i] = base/bpf[i]
// A 32 bits integer can have at most 9 distinct prime factors.

	i = 0;

	b = base;  // copy of base, to be completely factored...

	if (!(base & 1)) {	// divisor two?
            bpf[i] = 2;
            while (!(b & 1)) {
                vpf[i]++;	// compute the power of two
                b >>= 1;
            }
            bpc[i++] = base/2;
            if (isPrime (b)) {	// b may be the last factor!
                bpf[i] = b;
                vpf[i] = 1;
                bpc[i] = base/b;
                return;
            }
	}

	if (!(base%3)) {		// divisor three?
            bpf[i] = 3;
            while (!(b % 3)) {
                vpf[i]++;	// compute the power of three
                b /= 3;
            }
            bpc[i++] = base/3;
            if (isPrime (b)) {	// b may be the last factor!
                bpf[i] = b;
                vpf[i] = 1;
                bpc[i] = base/b;
                return;
            }
	}

	p = 5;

	while (p*p <= base) {		// other divisors?
            if (!(base%p) && isPrime (p)) {
                bpf[i] = p;
                while (!(b % p)) {
                    vpf[i]++;   // compute the power of p
                    b /= p;
                }
                bpc[i++] = base/p;
                if (isPrime (b)) {  // b may be the last factor!
                    bpf[i] = b;
                    vpf[i] = 1;
                    bpc[i] = base/b;
                    return;
                }
            }
            p += 2;
            if (!(base%p) && isPrime (p)) {
                bpf[i] = p;
                while (!(b % p)) {
                    vpf[i]++;   // compute the power of p
                    b /= p;
                }
                bpc[i++] = base/p;
                if (isPrime (b)) {  // b may be the last factor!
                    bpf[i] = b;
                    vpf[i] = 1;
                    bpc[i] = base/b;
                    return;
                }
            }
            p += 4;
	}

	if (i == 0 || base < 4) {   // the base is prime!
            bpf[0] = base;
            vpf[0] = bpc[0] = 1;
	}
}

int findgbpf (giant gbase) {    // find all prime factors of a large integer base
	unsigned long p, pmax;
	int i;
	double db;
	giant b;
	b = newgiant (2*abs(gbase->sign) + 8);
	for (i=0; i<30; i++) {
		bpf[i]  = vpf[i] = 0;	// clean up
		if (gbpc[i] != NULL) {
			free (gbpc[i]);
			gbpc[i] = NULL;
		}
		if (gbpf[i] != NULL) {
			free (gbpf[i]);
			gbpf[i] = NULL;
		}
	}

// bpf[i] or gbpf[i] : base prime factor, vpf[i] : its exponent, gbpc[i] = gbase/bpf[i] or gbase/gbpf[i]
// We expect the base has no more than 30 distinct prime factors.

	i = 0;
	gtog (gbase, b);   // copy of base, to be completely factored...

	if (!bitval (gbase, 0)) {  // divisor two?
            bpf[i] = 2;
            while (!bitval (b, 0)) {
                vpf[i]++;  // compute the exponent of two
                gshiftright (1, b);
            }
            gbpc[i] = newgiant (2*abs(gbase->sign) + 8);
            gtog (gbase, gbpc[i]);
            gshiftright (1, gbpc[i]);
            i++;
            if ((b->sign <= 2) && isPrime (p = (unsigned long)gtoi (b))) {
                // b may be the last factor!
                bpf[i] = p;
                vpf[i] = 1;
                gbpc[i] = newgiant (2*abs(gbase->sign) + 8);
                gtog (gbase, gbpc[i]);
                uldivg (p, gbpc[i]);
                free (b);
                return TRUE;
            }
	}
	if (!gmodi (3, gbase)) {   // divisor three?
            bpf[i] = 3;
            while (!gmodi (3, b)) {
                vpf[i]++;  // compute the exponent of three
                uldivg (3, b);
            }
            gbpc[i] = newgiant (2*abs(gbase->sign) + 8);
            gtog (gbase, gbpc[i]);
            uldivg (3, gbpc[i]);
            i++;
            if ((b->sign <= 2) && isPrime (p = (unsigned long)gtoi (b))) {
                // b may be the last factor!
                bpf[i] = p;
                vpf[i] = 1;
                gbpc[i] = newgiant (2*abs(gbase->sign) + 8);
                gtog (gbase, gbpc[i]);
                uldivg (p, gbpc[i]);
                free (b);
                return TRUE;
            }
	}
	if (bitlen(b) > 53) {
            // The cofactor is still a large integer...
            pmax = 1<<20;
	}
	else {								// Compute the cofactor as a double
            db = (double)b->n[0];
            if (b->sign > 1)
                db += 65536.0*(double)b->n[1];
            if (b->sign > 2)
                db += 65536.0*65536.0*(double)b->n[2];
            if (b->sign > 3)
                db += 65536.0*65536.0*65536.0*(double)b->n[3];
		pmax = (unsigned long)floor (sqrt (db));
		if (pmax > (1<<20))
                    pmax = 1<<20;
// 2**40 > 2**32, so, a cofactor not larger than 32bit must be prime!
	}

	p = 5;

	while (p <= pmax) {	// other divisors?
            if (isPrime(p) && !gmodi (p, gbase)) {
                bpf[i] = p;
                while (!gmodi (p, b)) {
                    vpf[i]++;	// compute the exponent of p
                    uldivg (p, b);
                }
                gbpc[i] = newgiant (2*abs(gbase->sign) + 8);
                gtog (gbase, gbpc[i]);
                uldivg (p, gbpc[i]);
                i++;
                if ((b->sign <= 2) && isPrime ((unsigned long)gtoi (b))){               // b may be the last prime factor!
                    bpf[i] = (unsigned long)gtoi (b);
                    vpf[i] = 1;
                    gbpc[i] = newgiant (2*abs(gbase->sign) + 8);
                    gtog (gbase, gbpc[i]);
                    uldivg (bpf[i], gbpc[i]);
                    free (b);
                    return TRUE;
                }
            }
            p += 2;
            if (isPrime(p) && !gmodi (p, gbase)) {
                bpf[i] = p;
                while (!gmodi (p, b)) {
                    vpf[i]++;   // compute the exponent of p
                    uldivg (p, b);
                }
                gbpc[i] = newgiant (2*abs(gbase->sign) + 8);
                gtog (gbase, gbpc[i]);
                uldivg (p, gbpc[i]);
                i++;
                if ((b->sign <= 2) && isPrime ((unsigned long)gtoi (b))){               // b may be the last prime factor!
                    bpf[i] = (unsigned long)gtoi (b);
                    vpf[i] = 1;
                    gbpc[i] = newgiant (2*abs(gbase->sign) + 8);
                    gtog (gbase, gbpc[i]);
                    uldivg (bpf[i], gbpc[i]);
                    free (b);
                    return TRUE;
                }
            }
            p += 4;
	}
	if (isone(b)) {	// The factorization is complete
            free(b);
            return TRUE;
	}
	else if (abs(b->sign) <= 2) {
            // The cofactor is a small integer prime !
            bpf[i] = (unsigned long)gtoi (b);
            vpf[i] = 1;
            gbpc[i] = newgiant (2*abs(gbase->sign) + 8);
            gtog (gbase, gbpc[i]);
            uldivg (bpf[i], gbpc[i]);
            free (b);
            return TRUE;
	}
	else {
            bpf[i] = vpf[i] = 1;
            // To signal a large integer cofactor
            gbpc[i] = newgiant (2*abs(gbase->sign) + 8);
            gtog (gbase, gbpc[i]);
            divg (b, gbpc[i]);
            gbpf[i] = newgiant (2*abs(gbase->sign) + 8);
            if (bitlen(b) <= 40 || (gaprcltest (b, 0, 0) == 2)) {
                gtog (b, gbpf[i]);  // The cofactor is prime !
                free (b);
                return TRUE;
            }
            else {
                gfact (b, gbpf[i], 0, 0, debug);
                // Try to factorize using R.Crandall code
                gtog (gbase, gbpc[i]);
                divg (gbpf[i], gbpc[i]);
                i++;
                bpf[i] = vpf[i] = 1;
                // To signal a large integer cofactor
                gbpc[i] = newgiant (2*abs(gbase->sign) + 8);
                gtog (gbase, gbpc[i]);
                divg (b, gbpc[i]);
                gbpf[i] = newgiant (2*abs(gbase->sign) + 8);
                gtog (b, gbpf[i]);
                free (b);
                if ((bitlen(gbpf[i-1]) <= 40) && 
                (bitlen(gbpf[i]) <= 40))
                    return TRUE;        // The two factors are prime !
                else {
                    if ((gaprcltest (gbpf[i-1], 0, 0) == 2) && (gaprcltest (gbpf[i], 0 ,0) == 2))
                        return TRUE;    // The two factors are prime !
                    else
                        return FALSE;
                    // The factors must be factorized further...
                }
            }
	}
}

int Lucasequence (
	giant modulus,
	giant exponent,
	unsigned long P,
	unsigned long base,
	int jmin,
	int jmax,
	char *str,
	char *buf,
	int	*res)
{
	unsigned long bit, iters, D;
	unsigned long bbits, mask=0x80000000, frestart=FALSE;
	unsigned long ulone;
//	unsigned long maxrestarts, ulone;
	gwypnum x, y, gwinvD, gwinv2;
	gwypnum a11, a12, a21, a22, b11, b12, b21, b22, c11, c12, c21, c22, p, pp, s;
	giant	tmp, tmp2;
	char	filename[20], fft_desc[256];
	long	write_time = DISK_WRITE_TIME * 60;
	long	j;
	int	echk, saving, stopping;
	time_t	start_time, current_time;
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;


//	maxrestarts = IniGetInt(INI_FILE, (char*)"MaxRestarts", 10);


restart:

	/* Get the current time */

	time (&start_time);

/* Allocate memory */

	x = gwypalloc ();
	y = gwypalloc ();
	gwinvD = gwypalloc ();
	gwinv2 = gwypalloc ();

	a11 = gwypalloc ();		// allocate memory for matrix computing...
	a12 = gwypalloc ();
	a21 = gwypalloc ();
	a22 = gwypalloc ();
	b11 = gwypalloc ();
	b12 = gwypalloc ();
	b21 = gwypalloc ();
	b22 = gwypalloc ();
	c11 = gwypalloc ();
	c12 = gwypalloc ();
	c21 = gwypalloc ();
	c22 = gwypalloc ();
	p = gwypalloc ();
	pp = gwypalloc ();
	s = gwypalloc ();

	tmp = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);
	tmp2 = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);

/* Init, */


	gtog (exponent, tmp2);
	uldivg (base, tmp2);	// tmp2 = exponent/base

	Nlen = bitlen (tmp2);

/* Init filename */

	tempFileName (filename, 'L', N);

/* Optionally resume from save file and output a message */
/* indicating we are resuming a test */

	if (fileExists (filename) && readFromFileB (filename, &bit, &P, &nrestarts, bpf, x, y)) {
		char	fmt_mask[80];
		double	pct;
		pct = trunc_percent (bit * 100.0 / Nlen);
		sprintf (fmt_mask,
			 "Resuming Lucas sequence at bit %%ld [%%.%df%%%%]\n",
			 PRECISION);
		sprintf (buf, fmt_mask, bit, pct);
		OutputStr (buf);
		if (verbose)
			writeResults (buf);	
		D = P*P - 4;
	}

/* Otherwise, output a message indicating we are starting test */

	else {
		OutputStr (buf);
		if (verbose)
			writeResults (buf);	
		D = P*P - 4;
		bit = 1;
//		dbltogw (2.0, x);			// Initial values
//		dbltogw ((double)P, y);
		itogwyp (2, x);			// Initial values
		itogwyp (P, y);
	}


	itog (D, tmp);						// Compute inverse of D modulo N
	invg (modulus, tmp);
	gianttogwyp (tmp, gwinvD);	// Convert it to gwypnum
	gtog (modulus, tmp);				// Compute inverse of 2 modulo N
	iaddg (1, tmp);						// Which is (N+1)/2
	gshiftright (1, tmp);
	gianttogwyp (tmp, gwinv2);	// Convert it to gwypnum

/* Output a message about the FFT length */

	gwypfft_description (fft_desc);
#ifdef WIN32
	sprintf (buf, "%s, P = %d\n", fft_desc, (int)P);
#else
	sprintf (buf, "%s, P = %d", fft_desc, (int)P);
#endif
	OutputStr (buf);
	LineFeed ();
	if (verbose) {
#if !defined(WIN32) 
		strcat (buf, "\n");
#endif
		writeResults (buf);
	}
	ReplaceableLine (1);	/* Remember where replaceable line is */

/* Init the title */

	title ((char*)"Computing Lucas sequence...");

/* Main loop... */

	iters = 0;
	gwypstart_timer (0);
	gwypstart_timer (1);
	will_try_larger_fft = FALSE;
	while (bit <= Nlen) {

/* Error check the last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		stopping = stopCheck ();
		echk = stopping || ERRCHK || (bit <= 50) || (bit >= Nlen-50);
		if (((bit & 127) == 0) || (bit == 1) || (bit == (lasterr_point-1))) {
			echk = 1;
			time (&current_time);
			saving = ((current_time - start_time > write_time) || (bit == 1) || (bit == (lasterr_point-1)));
		} else
			saving = 0;

/* Process this bit */


		gwypsetnormroutine (0, echk, 0);

		if (bitval (tmp2, Nlen-bit)) {
			gwypsetaddin (-(int)P);
			if (/*(bit+26 < Nlen) && (bit > 26) &&*/
				((bit != lasterr_point) || (!maxerr_recovery_mode[1] && !maxerr_recovery_mode[2]))) {
                            if (cufftonly)
                                gwypmul (y, x);
                            else
                                cuda_gwypmul (y, x, 3);
                            will_try_larger_fft = FALSE;
			}
			else {
				gwypmul_carefully (y, x);
				if (bit == lasterr_point)
					maxerr_recovery_mode[1] = FALSE;
				will_try_larger_fft = TRUE;
			}
			CHECK_IF_ANY_ERROR(x, (bit), Nlen, 1)
			gwypsetaddin (-2);
			if (/*(bit+26 < Nlen) && (bit > 26) && */
				((bit != lasterr_point) || !maxerr_recovery_mode[2])) {
                            if (cufftonly)
                                gwypsquare (y);
                            else if (zp || generic)
                                cuda_gwypsquare (y, 3);
                            else if(bit==1 || it==0)  
                                {cuda_gwypsquare (y,1);it=1;}
                            else  if(bit != (lasterr_point-1)) 
                                cuda_gwypsquare (y,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (y,2);
                            will_try_larger_fft = FALSE;
			}
			else {
				gwypsquare_carefully (y);
				if (bit == lasterr_point)
					maxerr_recovery_mode[2] = FALSE;
				will_try_larger_fft = TRUE;
			}
			CHECK_IF_ANY_ERROR(y, (bit), Nlen, 2)
		}
		else {
			gwypsetaddin (-(int)P);
			if (/*(bit+26 < Nlen) && (bit > 26) && */
				((bit != lasterr_point) || (!maxerr_recovery_mode[3] && !maxerr_recovery_mode[4]))) {
                            if (cufftonly)
                                gwypmul (x, y);
                            else
                                cuda_gwypmul (x, y, 3);
                            will_try_larger_fft = FALSE;
			}
			else {
                            gwypmul_carefully (x, y);
                            if (bit == lasterr_point)
                                    maxerr_recovery_mode[3] = FALSE;
                            will_try_larger_fft = TRUE;
			}
			CHECK_IF_ANY_ERROR(y, (bit), Nlen, 3)
			gwypsetaddin (-2);
			if (/*(bit+26 < Nlen) && (bit > 26) &&*/
				((bit != lasterr_point) || !maxerr_recovery_mode[4])) {
                            if (cufftonly)
                                gwypsquare (x);
                            else if (zp || generic)
                                cuda_gwypsquare (x, 3);
                            else if(bit==1 || it==0)  
                                {cuda_gwypsquare (x,1);it=1;}
                            else  if(bit != (lasterr_point-1)) 
                                cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (x,2);
                            will_try_larger_fft = FALSE;
			}
			else {
				gwypsquare_carefully (x);
				if (bit == lasterr_point)
					maxerr_recovery_mode[4] = FALSE;
				will_try_larger_fft = TRUE;
			}
			CHECK_IF_ANY_ERROR(x, (bit), Nlen, 4)
		}

 /* That iteration succeeded, bump counters */

		if (will_try_larger_fft && (bit == lasterr_point))
			saving = 1;					// Be sure to restart after this recovery iteration!
		will_try_larger_fft = FALSE;
		bit++;
		iters++;

/* Print a message every so often */

		if (bit % ITER_OUTPUT == 0) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (bit * 100.0 / Nlen);
			sprintf (fmt_mask, "%%.%df%%%% of %%ld", PRECISION);
			sprintf (buf, fmt_mask, pct, Nlen);
			title (buf);
			ReplaceableLine (2);	/* Replace line */
			sprintf (fmt_mask,
				 "%%s, bit: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, str, bit, Nlen, pct);
			OutputStr (buf);
			if (ERRCHK && bit > 30) {
				OutputStr ((char*)".  Round off: ");
				sprintf (buf, "%10.10f", reallyminerr);
				OutputStr (buf);
				sprintf (buf, " to %10.10f", reallymaxerr);
				OutputStr (buf);
			}
			gwypend_timer (0);
			if (CUMULATIVE_TIMING) {
				OutputStr ((char*)".  Time thusfar: ");
			} else {
				OutputStr ((char*)".  Time per bit: ");
				gwypdivide_timer (0, iters);
				iters = 0;
			}
			gwypprint_timer (0, TIMER_NL | TIMER_OPT_CLR);
			gwypstart_timer (0);
		}

/* Print a results file message every so often */

		if (bit % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
			sprintf (buf, "Bit %ld / %ld\n", bit, Nlen);
			writeResults (buf);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving || stopping) {
			write_time = DISK_WRITE_TIME * 60;
			if (! writeToFileB (filename, bit, P, nrestarts, bpf, x, y)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (buf);
				if (write_time > 600) write_time = 600;
			}	
			time (&start_time);

/* If an escape key was hit, write out the results and return */

			if (stopping) {
				gwypfree (tmp);
				gwypfree (tmp2);
				gwypfree (x);				// Clean up
				gwypfree (y);
				gwypfree (gwinvD);
				gwypfree (gwinv2);
				gwypfree (a11);			// Clean up the matrix
				gwypfree (a12);
				gwypfree (a21);
				gwypfree (a22);
				gwypfree (b11);
				gwypfree (b12);
				gwypfree (b21);
				gwypfree (b22);
				gwypfree (c11);
				gwypfree (c12);
				gwypfree (c21);
				gwypfree (c22);
				gwypfree (p);
				gwypfree (pp);
				gwypfree (s);
//				gwypdone ();
				Nlen = bitlen (N);
				*res = FALSE;
				return (FALSE);
			}
		}

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next iteration so that we can compare our */
/* residues to programs that start counter at zero or one. */

		if (interimResidues && bit % interimResidues < 2) {
			gwyptogiant (x, tmp);		// The modulo reduction is done here
			modg (N, tmp);				// External modulus and gwypnums one may be different...
			if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
			else if (abs(tmp->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
			else if (abs(tmp->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
			sprintf (buf, "%s interim residue %s at bit %ld\n", str, res64, bit);
			OutputBoth (buf);
		}

/* Write a save file every "interimFiles" iterations. */

		if (interimFiles && bit % interimFiles == 0) {
			char	interimfile[20];
			sprintf (interimfile, "%.8s.%03lu",
				 filename, bit / interimFiles);
			if (! writeToFileB (interimfile, bit, P, nrestarts, bpf, x, y)) {
				sprintf (buf, WRITEFILEERR, interimfile);
				OutputBoth (buf);
			}
		}
	}

	will_try_larger_fft = TRUE;	// All following errors are considered unrecoverable...

	// Compute the matrix at (N+1)/base

	gwypsetaddin (0);		// Reset addin constant.
	gwypcopy (x, a22);	        // a22 = V((N+1)/base)
	gwypcopy (y, a12);		// a12 = V((N+1)/base+1)
	gwypcopy (y, a21);		// a21 = V((N+1)/base+1)
	gwypcopy (y, a11);		// a11 = V((N+1)/base+1)
	gwypcopy (x, y);			// Now, y = V((N+1)/base)
	gwypsetnormroutine (0, 1, 1);	// set mul. by const.
	gwypsetmulbyconst (2);
//	gwypmul (gwinvD, a21);		// a21 = D^-1*2*V((N+1)/base+1) modulo N
        if (cufftonly)
            gwypmul (gwinvD, a21);
        else
            cuda_gwypmul (gwinvD, a21, 3);
	CHECK_IF_ANY_ERROR(a21, (Nlen), Nlen, 6)
	gwypsetmulbyconst (P);
//	gwypmul (gwinvD, x);		// x =  D^-1*P*V((N+1)/base) modulo N
        if (cufftonly)
            gwypmul (gwinvD, x);
        else
            cuda_gwypmul (gwinvD, x, 3);
	CHECK_IF_ANY_ERROR(x, (Nlen), Nlen, 6)
	gwypsub (x, a21);		// a21 = D^-1*(2*V((N+1)/base+1)-P*V((N+1)/base)) = U(N+1)/base modulo N
//	gwypmul (gwinvD, a11);		// a11 = D^-1*P*V((N+1)/base+1) modulo N
        if (cufftonly)
            gwypmul (gwinvD, a11);
        else
            cuda_gwypmul (gwinvD, a11, 3);
	CHECK_IF_ANY_ERROR(a11, (Nlen), Nlen, 6)
	gwypsetmulbyconst (2);
//	gwypmul (gwinvD, y);		// xx = D^-1*2*V((N+1)/base) modulo N
        if (cufftonly)
            gwypmul (gwinvD, y);
        else
            cuda_gwypmul (gwinvD, y, 3);
	CHECK_IF_ANY_ERROR(y, (Nlen), Nlen, 6)
	gwypsub (y, a11);		// a11 = D^-1*(P*V((N+1)/base+1)-2*V((N+1)/base)) = U((N+1)/base+1) modulo N
	gwypsetnormroutine (0, 1, 0);	// reset mul by const
//	gwypmul (gwinv2, a22);		// a22 = 2^-1*V((N+1)/base)
        if (cufftonly)
            gwypmul (gwinv2, a22);
        else
            cuda_gwypmul (gwinv2, a22, 3);
	CHECK_IF_ANY_ERROR(a22, (Nlen), Nlen, 6)
//	gwypmul (gwinv2, a12);		// a12 = 2^-1*V((N+1)/base+1)
        if (cufftonly)
            gwypmul (gwinv2, a12);
        else
            cuda_gwypmul (gwinv2, a12, 3);
	CHECK_IF_ANY_ERROR(a12, (Nlen), Nlen, 6)
	gwypsetmulbyconst (P);
	gwypsetnormroutine (0, 1, 1);	// set mul. by const.
	gwypcopy (a11, x);		// x = U((N+1)/base+1)
//	gwypmul (gwinv2, x);		// x = 2^-1*P*U((N+1)/base+1)
        if (cufftonly)
            gwypmul (gwinv2, x);
        else
            cuda_gwypmul (gwinv2, x, 3);
	CHECK_IF_ANY_ERROR(x, (Nlen), Nlen, 6)
	gwypsub (x, a12);		// a12 = 2^-1(V((N+1)/base+1)-P*U((N+1)/base+1))
	gwypcopy (a21, x);		// x = U((N+1)/base)
//	gwypmul (gwinv2, x);		// x = 2^-1*P*U((N+1)/base)
        if (cufftonly)
            gwypmul (gwinv2, x);
        else
            cuda_gwypmul (gwinv2, x, 3);
	CHECK_IF_ANY_ERROR(x, (Nlen), Nlen, 6)
	gwypsub (x, a22);		// a22 = 2^-1(V((N+1)/base)-P*U((N+1)/base))
	gwypsetnormroutine (0, 1, 0);	// reset mul by const

//	gwyptogiant (a21, tmp);		// tmp = U((N+1)/base) modulo N

	gwypcopy (a11, c11);		// Save the current matrix
	gwypcopy (a12, c12);
	gwypcopy (a21, c21);
	gwypcopy (a22, c22);

	gwypcopy (a11, b11);		// Copy the current matrix
	gwypcopy (a12, b12);
	gwypcopy (a21, b21);
	gwypcopy (a22, b22);

	bbits = base;
	ulone = 1;
	while (!(bbits & mask)) {
		ulone <<= 1;
		bbits <<= 1;
	}
	bbits <<= 1;
	ulone <<= 1;

	while (ulone) {				// Finish to compute U(N+1)

// Square the matrix

		gwypcopy (a12, p);		// a12-->p
		gwypcopy (a12, pp);		// a12-->pp
		gwypadd3 (a11, a22, s);		// a11+a12-->s
//		gwypmul (a21, p);			// a21*a12-->p
                if (cufftonly)
                    gwypmul (a21, p);
                else
                    cuda_gwypmul (a21, p, 3);
		CHECK_IF_ANY_ERROR(p, (Nlen), Nlen, 6)
//		gwypsquare (a22);			// a22*a22-->a22
                if (cufftonly)
                    gwypsquare (a22);
                else
                    cuda_gwypsquare (a22, 3);
		CHECK_IF_ANY_ERROR(a22, (Nlen), Nlen, 6)
//		gwypmul (s, a21);			// (a11+a22)*a21-->a21 T
                if (cufftonly)
                    gwypmul (s, a21);
                else
                    cuda_gwypmul (s, a21, 3);
		CHECK_IF_ANY_ERROR(a21, (Nlen), Nlen, 6)
		gwypadd (p, a22);			// a21*a12+a22*a22-->a22 T
//		gwypmul (s, a12);			// (a11+a22)*a12-->a12 T
                if (cufftonly)
                    gwypmul (s, a12);
                else
                    cuda_gwypmul (s, a12, 3);
		CHECK_IF_ANY_ERROR(a12, (Nlen), Nlen, 6)
//		gwypsquare (a11);			// a11*a11-->a11
                if (cufftonly)
                    gwypsquare (a11);
                else
                    cuda_gwypsquare (a11, 3);
		CHECK_IF_ANY_ERROR(a11, (Nlen), Nlen, 6)
		gwypadd (p, a11);			// a21*a12+a11*a11-->a11 T

// Multiply it if required

		if (bbits & mask) {
			gwypcopy (a11, p);		// a11-->p
			gwypcopy (a21, pp);		// a21-->pp
//			gwypmul (b11, a11);	        // b11*a11-->a11
                        if (cufftonly)
                            gwypmul (b11, a11);
                        else
                            cuda_gwypmul (b11, a11, 3);
			CHECK_IF_ANY_ERROR(a11, (Nlen), Nlen, 6)
//			gwypmul (b12, pp);		// b12*a21-->pp
                        if (cufftonly)
                            gwypmul (b12, pp);
                        else
                            cuda_gwypmul (b12, pp, 3);
			CHECK_IF_ANY_ERROR(pp, (Nlen), Nlen, 6)
			gwypadd (pp, a11);		// b11*a11+b12*a21-->a11 T
//			gwypmul (b21, p);		// b21*a11-->p
                        if (cufftonly)
                            gwypmul (b21, p);
                        else
                            cuda_gwypmul (b21, p, 3);
			CHECK_IF_ANY_ERROR(p, (Nlen), Nlen, 6)
//			gwypmul (b22, a21);	        // b22*a21-->a21
                        if (cufftonly)
                            gwypmul (b22, a21);
                        else
                            cuda_gwypmul (b22, a21, 3);
			CHECK_IF_ANY_ERROR(a21, (Nlen), Nlen, 6)
			gwypadd (p, a21);		// b21*a11+b22*a21-->a21 T
			gwypcopy (a12, p);		// a12-->p
			gwypcopy (a22, pp);		// a22-->pp
//			gwypmul (b11, a12);	        // b11*a12-->a12
                        if (cufftonly)
                            gwypmul (b11, a12);
                        else
                            cuda_gwypmul (b11, a12, 3);
			CHECK_IF_ANY_ERROR(a12, (Nlen), Nlen, 6)
//			gwypmul (b12, pp);		// b12*a22-->pp
                        if (cufftonly)
                            gwypmul (b12, pp);
                        else
                            cuda_gwypmul (b12, pp, 3);
			CHECK_IF_ANY_ERROR(pp, (Nlen), Nlen, 6)
			gwypadd (pp, a12);		// b11*a12+b12*a22-->a12 T
//			gwypmul (b21, p);		// b21*a12-->p
                        if (cufftonly)
                            gwypmul (b21, p);
                        else
                            cuda_gwypmul (b21, p, 3);
			CHECK_IF_ANY_ERROR(p, (Nlen), Nlen, 6)
//			gwypmul (b22, a22);	        // b22*a22-->a22
                        if (cufftonly)
                            gwypmul (b22, a22);
                        else
                            cuda_gwypmul (b22, a22, 3);
			CHECK_IF_ANY_ERROR(a22, (Nlen), Nlen, 6)
			gwypadd (p, a22);		// b21*a12+b22*a22-->a22 T

		}
		bbits <<= 1;
		ulone <<= 1;
	}

	clearline (100);

	gwyptogiant (a21, tmp);			// tmp = U(N+1) modulo N
	ReplaceableLine (2);			/* Replace line */

	if (!isZero (tmp)) {
//	if (!gwiszero(a21)) {
		*res = FALSE;				/* Not a prime */
		if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
		else if (abs(tmp->sign) < 3)
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
		else if (abs(tmp->sign) < 4)
			sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
		else
			sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
		if (IniGetInt(INI_FILE, (char*)"Verify", 0))
			sprintf (buf, "%s is not prime. P = %lu, Lucas RES64: %s", str, P, res64);
		else
			sprintf (buf, "%s is not prime, although Fermat PSP! P = %lu, Lucas RES64: %s", str, P, res64);
	}
	else {
		sprintf (buf, "%s may be prime, trying to compute gcd's\n", str);
		if (verbose)
			OutputBoth (buf);
		else
			OutputStr (buf);
		for (j=jmax; j>=jmin; j--) {
			if (bpf[j] == 1)			// base prime factor already tested
				continue;
			gwypcopy (c11, a11);			// Computing U((N+1)/q)
			gwypcopy (c12, a12);
			gwypcopy (c21, a21);
			gwypcopy (c22, a22);
			bbits = bpc[j];
			ulone = 1;
			while (!(bbits & mask)) {
				bbits <<= 1;
				ulone <<= 1;
			}
			bbits <<= 1;
			ulone <<= 1;
			while (ulone) {

// Square the matrix

				gwypcopy (a12, p);		// a12-->p
				gwypcopy (a12, pp);		// a12-->pp
				gwypadd3 (a11, a22, s);		// a11+a12-->s
//				gwypmul (a21, p);		// a21*a12-->p
                                if (cufftonly)
                                    gwypmul (a21, p);
                                else
                                    cuda_gwypmul (a21, p, 3);
				CHECK_IF_ANY_ERROR(p, (Nlen), Nlen, 6)
//				gwypsquare (a22);		// a22*a22-->a22
                                if (cufftonly)
                                    gwypsquare (a22);
                                else
                                    cuda_gwypsquare (a22, 3);
				CHECK_IF_ANY_ERROR(a22, (Nlen), Nlen, 6)
//				gwypmul (s, a21);		// (a11+a22)*a21-->a21 T
                                if (cufftonly)
                                    gwypmul (s, a21);
                                else
                                    cuda_gwypmul (s, a21, 3);
				CHECK_IF_ANY_ERROR(a21, (Nlen), Nlen, 6)
				gwypadd (p, a22);		// a21*a12+a22*a22-->a22 T
//				gwypmul (s, a12);		// (a11+a22)*a12-->a12 T
                                if (cufftonly)
                                    gwypmul (s, a12);
                                else
                                    cuda_gwypmul (s, a12, 3);
				CHECK_IF_ANY_ERROR(a12, (Nlen), Nlen, 6)
//				gwypsquare (a11);		// a11*a11-->a11
                                if (cufftonly)
                                    gwypsquare (a11);
                                else
                                    cuda_gwypsquare (a11, 3);
				CHECK_IF_ANY_ERROR(a11, (Nlen), Nlen, 6)
				gwypadd (p, a11);		// a21*a12+a11*a11-->a11 T

// Multiply it if required

				if (bbits & mask) {
					gwypcopy (a11, p);	// a11-->p
					gwypcopy (a21, pp);	// a21-->pp
//					gwypmul (b11, a11);	// b11*a11-->a11
                                        if (cufftonly)
                                            gwypmul (b11, a11);
                                        else
                                            cuda_gwypmul (b11, a11, 3);
					CHECK_IF_ANY_ERROR(a11, (Nlen), Nlen, 6)
//					gwypmul (b12, pp);	// b12*a21-->pp
                                        if (cufftonly)
                                            gwypmul (b12, pp);
                                        else
                                            cuda_gwypmul (b12, pp, 3);
					CHECK_IF_ANY_ERROR(pp, (Nlen), Nlen, 6)
					gwypadd (pp, a11);	// b11*a11+b12*a21-->a11 T
//					gwypmul (b21, p);	// b21*a11-->p
                                        if (cufftonly)
                                            gwypmul (b21, p);
                                        else
                                            cuda_gwypmul (b21, p, 3);
					CHECK_IF_ANY_ERROR(p, (Nlen), Nlen, 6)
//					gwypmul (b22, a21);	// b22*a21-->a21
                                        if (cufftonly)
                                            gwypmul (b22, a21);
                                        else
                                            cuda_gwypmul (b22, a21, 3);
					CHECK_IF_ANY_ERROR(a21, (Nlen), Nlen, 6)
					gwypadd (p, a21);	// b21*a11+b22*a21-->a21 T
					gwypcopy (a12, p);	// a12-->p
					gwypcopy (a22, pp);	// a22-->pp
//					gwypmul (b11, a12);	// b11*a12-->a12
                                        if (cufftonly)
                                            gwypmul (b11, a12);
                                        else
                                            cuda_gwypmul (b11, a12, 3);
					CHECK_IF_ANY_ERROR(a12, (Nlen), Nlen, 6)
//					gwypmul (b12, pp);	// b12*a22-->pp
                                        if (cufftonly)
                                            gwypmul (b12, pp);
                                        else
                                            cuda_gwypmul (b12, pp, 3);
					CHECK_IF_ANY_ERROR(pp, (Nlen), Nlen, 6)
					gwypadd (pp, a12);	// b11*a12+b12*a22-->a12 T
//					gwypmul (b21, p);	// b21*a12-->p
                                        if (cufftonly)
                                            gwypmul (b21, p);
                                        else
                                            cuda_gwypmul (b21, p, 3);
					CHECK_IF_ANY_ERROR(p, (Nlen), Nlen, 6)
//					gwypmul (b22, a22);	// b22*a22-->a22
                                        if (cufftonly)
                                            gwypmul (b22, a22);
                                        else
                                            cuda_gwypmul (b22, a22, 3);
					CHECK_IF_ANY_ERROR(a22, (Nlen), Nlen, 6)
					gwypadd (p, a22);	// b21*a12+b22*a22-->a22 T

				}
				bbits <<= 1;
				ulone <<= 1;
			}
			gwyptogiant (a21, tmp);
			if (isZero (tmp)) {
				sprintf (buf, "%s may be prime, but N divides U((N+1)/%lu), P = %lu\n", str, bpf[j], P);
				if (verbose)
					OutputBoth (buf);
				else
					OutputStr (buf);
				frestart = TRUE;
				_unlink (filename);
				continue;
//				break;
			}
			else {
				gcdg (modulus, tmp);
				if (isone (tmp)) {
					sprintf (buf, "U((N+1)/%lu) is coprime to N!\n", bpf[j]);
					OutputStr (buf);
					if (verbose)
						writeResults (buf);	
					bpf[j] = 1;
				}
				else {
					*res = FALSE;	/* Not a prime */
					if (IniGetInt(INI_FILE, (char*)"Verify", 0))
						sprintf (buf, "%s is not prime, although Lucas PSP!! (P = %lu)", str, P);
					else
						sprintf (buf, "%s is not prime, although Fermat and Lucas PSP!! (P = %lu)", str, P);
					break;
				}
			}
		}
	}

	if (*res && !frestart)
		sprintf (buf, "%s is prime! (%d decimal digits, P = %lu)", str, nbdg, P);

	will_try_larger_fft = FALSE;// Reset the "unrecoverable" condition.
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (x);			// Clean up
	gwypfree (y);
	gwypfree (gwinvD);
	gwypfree (gwinv2);
	gwypfree (a11);		// Clean up the matrix
	gwypfree (a12);
	gwypfree (a21);
	gwypfree (a22);
	gwypfree (b11);
	gwypfree (b12);
	gwypfree (b21);
	gwypfree (b22);
	gwypfree (c11);
	gwypfree (c12);
	gwypfree (c21);
	gwypfree (c22);
	gwypfree (p);
	gwypfree (pp);
	gwypfree (s);

/* Cleanup and return */

	Nlen = bitlen (N);
//	gwypdone ();
	_unlink (filename);
	lasterr_point = 0;

	if (frestart)
		return -2;
	else {
		IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
		return TRUE;
	}

/* An error occured, sleep, then try restarting at last save point. */

error:
	Nlen = bitlen (N);
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (x);				// Clean up
	gwypfree (y);
	gwypfree (gwinvD);
	gwypfree (gwinv2);
	gwypfree (a11);			// Clean up the matrix
	gwypfree (a12);
	gwypfree (a21);
	gwypfree (a22);
	gwypfree (b11);
	gwypfree (b12);
	gwypfree (b21);
	gwypfree (b22);
	gwypfree (c11);
	gwypfree (c12);
	gwypfree (c21);
	gwypfree (c22);
	gwypfree (p);
	gwypfree (pp);
	gwypfree (s);
	*res = FALSE;

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
		sprintf (buf, ERRMSG5, checknumber, str);
		OutputBoth (buf);
		will_try_larger_fft = FALSE;
		_unlink (filename);
		return (FALSE);
	}

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
		will_try_larger_fft = FALSE;
		return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
		OutputBoth (ERRMSG8);
		IniWriteInt(INI_FILE, (char*)"FFT_Increment", IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
		_unlink (filename);
		will_try_larger_fft = FALSE;
		return (-1);
	}
	goto restart;
}


int plusminustest ( 
	char *sgk,
	unsigned long base,
	unsigned long n, 
	int incr,
	unsigned long shift,
	int	*res) 
{ 
	char	filename[20], buf[sgkbufsize+256], str[sgkbufsize+256], sgk1[sgkbufsize], fft_desc[256], oldres64[17]; 
	unsigned long bits, explen, bbits, iters, bit, mask=0x80000000, frestart=FALSE;
	unsigned long newa, maxrestarts, ulone, factorized_part = 0;
//	uint32_t hi = 0, lo = 0;
	double dk;
	giant tmp, tmp2, tmp3;
	gwypnum x, y;
	long	a, P, write_time = DISK_WRITE_TIME * 60;
	int	echk, saving, stopping, jmin, jmax, j, retval, Psample;
	time_t	start_time, current_time;
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;

	if ((gformat == ABCDN) || (gformat == ABCDNG)) {// Compute gk = gb^(n-m)-1
                bits = ndiff*(unsigned long)ceil(log ((double)base)/log (2.0));// initial gksize
		gk = newgiant ((bits >> 2) + 8);
		itog (base, gk);
		power (gk, ndiff);
		iaddg (-1, gk);
		sprintf (str, "%lu^%lu-%lu^%lu%c%d", base, n+ndiff, base, n, incr < 0 ? '-' : '+', abs(incr));
	}
        else {
            gk = newgiant (strlen(sgk)/2 + 8);	// Allocate one byte per decimal digit + spares
            ctog (sgk, gk);		   // Convert k string to giant
            gshiftleft (shift, gk);	   // Shift k multiplier if requested
            gtoc (gk, sgk1, sgkbufsize);   // Updated k string
            if (!strcmp(sgk1, "1"))
		sprintf (str, "%lu^%lu%c%d", base, n, (incr < 0) ? '-' : '+', abs(incr));
            else
		sprintf (str, "%s*%lu^%lu%c%d", sgk1, base, n, (incr < 0) ? '-' : '+', abs(incr));
        }

	klen = bitlen(gk);

//	Be sure the base does not divide the gk multiplier :

	while (!(gmodi (base, gk))) {
		uldivg (base, gk);
		n++;
	}

	if ((int)klen != bitlen(gk))	// Has k been updated ?
		gtoc (gk, sgk1, sgkbufsize);
	else
		strcpy (sgk1, sgk);

	bits = (unsigned long) ((n * log(base)) / log(2) + bitlen(gk)); 
	N =  newgiant ((bits >> 2) + 16);  // Allocate memory for N


//	Compute the number we are testing.

	itog (base, N);
	power (N, n);

	Nlen = bitlen (N);	// Bit length of base^n

	mulg (gk, N); 

	iaddg (incr, N);

	klen = bitlen(gk);

	if (klen > 53) {	// we must use generic reduction
		dk = 0.0;
	}
	else {								// we can use DWT ; compute the multiplier as a double
		dk = (double)gk->n[0];
		if (gk->sign > 1)
			dk += 65536.0*(double)gk->n[1];
		if (gk->sign > 2)
			dk += 65536.0*65536.0*(double)gk->n[2];
		if (gk->sign > 3)
			dk += 65536.0*65536.0*65536.0*(double)gk->n[3];
	}


	nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.

	if ((klen > Nlen) && (nbdg > 400)) {
            if ((gformat == ABCDN) || (gformat == ABCDNG))
                sprintf(buf, "%lu^%lu-1 > %lu^%lu, so, only a PRP test is done for %s.\n", base, ndiff, base, n, str);
            else
                sprintf(buf, "%s > %lu^%lu, so, only a PRP test is done for %s.\n", sgk1, base, n, str);
	    OutputBoth(buf);
            retval = isPRPinternal (str, dk, base, n, incr, res);
            gwypfree(gk);
            gwypfree(N);
            return retval;
	}

	Nlen = bitlen (N);					// Bit length of N
	findbpf (base);						// Factorize the base
	for (jmax=9; (jmax>=0) && !bpf[jmax]; jmax--);
	jmin = 0;							// Choose the minimum required factored part.
	if (jmax) {							// The base is composite...
            factorized_part = bitlen (gk);
		for (j=jmax; j>=0; j--) {
			factorized_part += (unsigned long)floor(n*vpf[j]*log ((double)bpf[j])/log(2.0));
			if ((2*factorized_part) > Nlen)
				break;
		}
		jmin = j;
//	}
		sprintf (buf, "Base factorized as : ");

		for (j=0; j<=jmax; j++) {
			if (j<jmax) {
                            sprintf (buf+strlen(buf), "%lu", bpf[j]);
                            if (vpf[j]>1)
                                sprintf (buf+strlen(buf), "^%lu*", vpf[j]);
                            else
                                sprintf (buf+strlen(buf), "*");
			}
			else {
                            sprintf (buf+strlen(buf), "%lu", bpf[j]);
                            if (vpf[j]>1)
                                sprintf (buf+strlen(buf), "^%lu\n", vpf[j]);
                            else
                                sprintf (buf+strlen(buf), "\n");
			}
		}
		if (!setuponly)
			if (verbose)
				OutputBoth(buf);
			else
				OutputStr(buf);
        }       // End the base is composite.

	sprintf (buf, "Base prime factor(s) taken : ");

	for (j=jmin; j<=jmax; j++) {
		if (j<jmax)
			sprintf (buf+strlen(buf), "%lu, ", bpf[j]);
		else
			sprintf (buf+strlen(buf), "%lu\n", bpf[j]);
	}
        if (verbose)
            OutputBoth(buf);
        else
            OutputStr(buf);

	maxrestarts = IniGetInt(INI_FILE, (char*)"MaxRestarts", 10);
	nrestarts = IniGetInt (INI_FILE, (char*)"NRestarts", 0);
	if (!(a = IniGetInt (INI_FILE, (char*)"FermatBase", 0)))
		a = IniGetInt (INI_FILE, (char*)"FBase", 3);// The base for the PRP and Pocklington tests
	if (!(P = IniGetInt (INI_FILE, (char*)"LucasBaseP", 0))) {
		Psample = genLucasBaseP (N, IniGetInt (INI_FILE, (char*)"PBase", 3));
		if (Psample < 0) {
			if (Psample == -1)
				sprintf (buf, "Cannot compute P to test %s...\nThis is surprising, please, let me know that!!\nMy E-mail is jpenne@free.fr\n", str);
			else
				sprintf (buf, "%s has a small factor : %d !!\n", str, abs(Psample));
			OutputBoth (buf);
			gwypfree(gk);
			gwypfree(N);
			return (TRUE); 
		}
		else
			P = Psample;
	}
											// The Discriminant for the Morrison test
//	D = P*P-4;								// D = P^2 - 4*Q with Q = 1


// restart:

/* Setup the gwypnum code */


restart:

	*res = TRUE;

	M = newgiant ((bits >> 2) + 16);	// Allocate memory for M
	tmp = newgiant ((bits >> 2) + 16);	// Allocate memory for tmp
	gtog (N, tmp);
	iaddg (-1, tmp);			// tmp = N-1

	gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
	if (incr == +1) {
		gwypsetmaxmulbyconst (a);
		uldivg (base, tmp);		// tmp = (N-1)/base
		explen = bitlen (tmp);
		if (!setupok (gwypsetup (dk, base, n, +1,N), N, str, res)) {
			gwypfree(gk);
			gwypfree(N);
			gwypfree (M);
			gwypfree (tmp);
//			*res = FALSE;		// Not proven prime...
			return TRUE; 
		}
		tmp2 =newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);	// Allocate memory for tmp2
		tmp3 =newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);	// Allocate memory for tmp3
	}
	else {
		gwypsetmaxmulbyconst (max(a, P));
		gtog (N, M);
		iaddg (1, M);
		explen = bitlen (tmp);
		if (!setupok (gwypsetup (dk, base, n, -1, N), N, str, res)) {
			gwypfree(gk);
			gwypfree(N);
			gwypfree (M);
			gwypfree (tmp);
//			*res = FALSE;		// Not proven prime...
			return TRUE; 
		}
		tmp2 =newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);	// Allocate memory for tmp2
		tmp3 =newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);	// Allocate memory for tmp3
		tempFileName (filename, 'L', N);
		if (fileExists (filename)) {	// Resuming a Lucas sequence...
			goto DoLucas;
		}
		else if (IniGetInt(INI_FILE, (char*)"Verify", 0) || IniGetInt(INI_FILE, (char*)"PRPdone", 0)) {
			gwypclear_timers ();			// Init. timers
			sprintf (buf, "Starting Lucas sequence for %s...\n", str);
			goto DoLucas;				// Starting directly a Lucas sequence...
		}
	}

/* Init filename */

	tempFileName (filename, 'z', N);

/* Allocate memory */

	x = gwypalloc ();
	y = gwypalloc ();

/* Optionally resume from save file and output a message */
/* indicating we are resuming a test */

	//cuda if (fileExists (filename) && readFromFileB (filename, &bit, &a, &nrestarts, bpf, x, y)) {
	if (fileExists (filename) && readFromFileB (filename,(unsigned long *) &bit,(unsigned long *) &a, &nrestarts, bpf, x, y)) {
		char	fmt_mask[80];
		double	pct;
		pct = trunc_percent (bit * 100.0 / explen);
		sprintf (fmt_mask,
			 "Resuming N%%c%%d prime test of %%s at bit %%ld [%%.%df%%%%]\n",
			 PRECISION);
		sprintf (buf, fmt_mask, (incr < 0) ? '+' : '-', abs(incr), str, bit, pct);
		OutputStr (buf);
		if (verbose)
			writeResults (buf);	
	}

/* Otherwise, output a message indicating we are starting test */

	else {
		if (frestart) {
			sprintf (buf, "Restarting N%c%d prime test of %s\n", (incr < 0) ? '+' : '-', abs(incr), str);
			frestart = FALSE;
		}
		else {
			gwypclear_timers ();		// Init. timers
                        if (showdigits)
                            sprintf (buf, "Starting N%c%d prime test of %s (%d decimal digits)\n", incr < 0 ? '+' : '-', abs(incr), str, nbdg);
                        else
                            sprintf (buf, "Starting N%c%d prime test of %s\n", incr < 0 ? '+' : '-', abs(incr), str);
		}
		OutputStr (buf);
		if (verbose)
			writeResults (buf);	
		bit = 1;
		itogwyp (a, x);
//		dbltogw ((double)a, x);
	}

/* Get the current time */

	gwypstart_timer (0);
	gwypstart_timer (1);
	time (&start_time);

/* Output a message about the FFT length */

	gwypfft_description (fft_desc);

#ifdef WIN32
	sprintf (buf, "%s, a = %ld\n", fft_desc, a);
#else
	sprintf (buf, "%s, a = %ld", fft_desc, a);
#endif
	OutputStr (buf);
	LineFeed ();
	if (verbose) {
#if !defined(WIN32) 
		strcat (buf, "\n");
#endif
		writeResults (buf);
	}
	ReplaceableLine (1);	/* Remember where replaceable line is */

/* Init the title */

	if (incr > 0)
		title ((char*)"Pocklington prime test in progress...");
	else
		title ((char*)"Fermat PRP test in progress...");

/* Do the PRP test */

	gwypsetmulbyconst (a);
	gwypsetaddin(0);
	iters = 0;
	
	while (bit < explen) {

/* Error check the first and last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		stopping = stopCheck ();
		echk = stopping || ERRCHK || (bit <= 50) || (bit >= explen-50);
		if (((bit & 127) == 0) || (bit == 1) || (bit == (lasterr_point-1))) {
			echk = 1;
			time (&current_time);
			saving = ((current_time - start_time > write_time) || (bit == 1) || (bit == (lasterr_point-1)));
		} else
			saving = 0;

/* Process this bit */


		if (bitval (tmp, explen-bit-1)) {
			gwypsetnormroutine (0, echk, 1);
		} else {
			gwypsetnormroutine (0, echk, 0);
		}

		if (/*(bit+25 < explen) && (bit > 25) && */((bit != lasterr_point) || !maxerr_recovery_mode[6])) {
                    if (cufftonly)
                        gwypsquare (x);
                    else if (zp || generic)
                        cuda_gwypsquare (x, 3);
                    else if(bit==1 || it==0)  
                        {cuda_gwypsquare (x,1);it=1;}
                    else  if(bit != (lasterr_point-1)) 
                        cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                    else
                        cuda_gwypsquare (x,2);
		}
		else {
			gwypsquare_carefully (x);
			will_try_larger_fft = TRUE;
			if (bit == lasterr_point)
				maxerr_recovery_mode[6] = FALSE;
		}

		CHECK_IF_ANY_ERROR (x, (bit), explen, 6);

/* That iteration succeeded, bump counters */

		if (will_try_larger_fft && (bit == lasterr_point))
			saving = 1;					// Be sure to restart after this recovery iteration!
		will_try_larger_fft = FALSE;
		bit++;
		iters++;

/* Print a message every so often */

		if (bit % ITER_OUTPUT == 0) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (bit * 100.0 / explen);
			sprintf (fmt_mask, "%%.%df%%%% of %%ld", PRECISION);
			sprintf (buf, fmt_mask, pct, explen);
			title (buf);
			ReplaceableLine (2);	/* Replace line */
			sprintf (fmt_mask,
				 "%%s, bit: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, str, bit, explen, pct);
			OutputStr (buf);
			if (ERRCHK && bit > 30) {
				OutputStr ((char*)".  Round off: ");
				sprintf (buf, "%10.10f", reallyminerr);
				OutputStr (buf);
				sprintf (buf, " to %10.10f", reallymaxerr);
				OutputStr (buf);
			}
			gwypend_timer (0);
			if (CUMULATIVE_TIMING) {
				OutputStr ((char*)".  Time thusfar: ");
			} else {
				OutputStr ((char*)".  Time per bit: ");
				gwypdivide_timer (0, iters);
				iters = 0;
			}
			gwypprint_timer (0, TIMER_NL | TIMER_OPT_CLR);
			gwypstart_timer (0);
		}

/* Print a results file message every so often */

		if (bit % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
			sprintf (buf, "Bit %ld / %ld\n", bit, explen);
			writeResults (buf);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving || stopping) {
			write_time = DISK_WRITE_TIME * 60;
			if (! writeToFileB (filename, bit, a, nrestarts, bpf, x, y)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (buf);
				if (write_time > 600) write_time = 600;
			}	
			time (&start_time);

/* If an escape key was hit, write out the results and return */

			if (stopping) {
				gwypfree (N);
				gwypfree (gk);
				gwypfree (M);
				gwypfree (tmp);
				gwypfree (tmp2);
				gwypfree (tmp3);
				gwypfree (x);
				gwypfree (y);
				gwypdone ();
				*res = FALSE;
				return (FALSE);
			}
		}

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next iteration so that we can compare our */
/* residues to programs that start counter at zero or one. */

		if (interimResidues && bit % interimResidues < 2) {
  			gwyptogiant (x, tmp2);	// The modulo reduction is done here
			if (abs(tmp2->sign) < 2)		// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp2->n[0]);
			else if (abs(tmp2->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp2->n[1], tmp2->n[0]);
			else if (abs(tmp2->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp2->n[2], tmp2->n[1], tmp2->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp2->n[3], tmp2->n[2], tmp2->n[1], tmp2->n[0]);
			sprintf (buf, "%s interim residue %s at bit %ld\n",str, res64, bit);
			OutputBoth (buf);
		}

/* Write a save file every "interimFiles" iterations. */

		if (interimFiles && bit % interimFiles == 0) {
			char	interimfile[20];
			sprintf (interimfile, "%.8s.%03lu",
				 filename, bit / interimFiles);
			if (! writeToFileB (interimfile, bit, a, nrestarts, bpf, x, y)) {
				sprintf (buf, WRITEFILEERR, interimfile);
				OutputBoth (buf);
			}
		}
	}

	clearline (100);

	if (incr == +1) {
		will_try_larger_fft = TRUE;			// All following errors are considered unrecoverable...
		bbits = base;
		ulone = 1;
		while (!(bbits & mask)) {
			bbits <<= 1;
			ulone <<= 1;
		}
		bbits <<= 1;
		ulone <<= 1;
		gwypcopy (x, y);
		gwypsetnormroutine (0, 1, 0);
		while (ulone) {
                    if (cufftonly)
                        gwypsquare (x);
                    else
                        cuda_gwypsquare (x, 3);
                    CHECK_IF_ANY_ERROR (x, (explen), explen, 6);
                    if (bbits & mask) {
                        if (cufftonly)
                            gwypmul (y, x);
                        else
                            cuda_gwypmul (y, x, 3);
                        CHECK_IF_ANY_ERROR (x, (explen), explen, 6);
                    } 
                    bbits <<= 1;
                    ulone <<= 1;
		}
	}

/* See if we've found a probable prime.  If not, format a 64-bit residue. */
/* Old versions of PRP used a non-standard 64-bit residue, computing */
/* 3^N-3 mod N rather than the more standard 3^(N-1) mod N.  Since */
/* some projects recorded these non-standard residues, output that */
/* residue too.  Note that some really old versions printed out the */
/* 32-bit chunks of the non-standard residue in reverse order. */

	ReplaceableLine (2);	/* Replace line */
	gwyptogiant (x, tmp);
	if (!isone (tmp)) {
		*res = FALSE;	/* Not a prime */
		if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
		else if (abs(tmp->sign) < 3)
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
		else if (abs(tmp->sign) < 4)
			sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
		else
			sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
		smulg ((unsigned short)a, tmp); modg (N, tmp); iaddg (-a, tmp);
		if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
			sprintf (oldres64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
		else if (abs(tmp->sign) < 3)
			sprintf (oldres64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
		else if (abs(tmp->sign) < 4)
			sprintf (oldres64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
		else
			sprintf (oldres64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
		if (IniGetInt (INI_FILE, (char*)"OldRes64", 1))
			sprintf (buf, "%s is not prime.  RES64: %s.  OLD64: %s", str, res64, oldres64);
		else
			sprintf (buf, "%s is not prime.  RES64: %s", str, res64);
	}

	if (*res) {
		gwypend_timer (1);
		_unlink (filename);
		if (incr == -1) {				// Morrison test ; start the Lucas sequence
			gwypfree (x);
			gwypfree (y);
			sprintf (buf, "%s may be prime. Starting Lucas sequence...\n", str);
			IniWriteInt(INI_FILE, (char*)"PRPdone", 1);
DoLucas:
			do {
				retval = Lucasequence (N, M, P, base, jmin, jmax, str, buf, res);
				if (retval == -2) {		// Restart required using next base
					nrestarts++;
					if (nrestarts > maxrestarts) {
						sprintf (buf, "Giving up after %lu restarts...", nrestarts);
						frestart = FALSE;
						*res = FALSE;		// Not proven prime...
						retval = TRUE;
						break;
					}
					IniWriteInt (INI_FILE, (char*)"NRestarts", nrestarts);
					P = genLucasBaseP (N, P+1);
					IniWriteInt (INI_FILE, (char*)"LucasBaseP", P);
//					D = P*P-4;
				}
				if (retval < 0)	{		// Restart required for any reason
					sprintf (buf, "Restarting Lucas sequence with P = %lu\n", P);
					gwypdone ();
								// Setup again the gwypnum code.
					gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
					gwypsetmaxmulbyconst (max(a, P));
					if (!setupok (gwypsetup (dk, base, n, -1, N), N, str, res)) {
						gwypfree(gk);
						gwypfree(N);
						gwypfree (M);
						gwypfree (tmp);
						gwypfree (tmp2);
						gwypfree (tmp3);
//						*res = FALSE;		// Not proven prime...
						return TRUE; 
					}
				}
			}	while (retval < 0);
			if (retval == FALSE) {
				gwypfree (N);
				gwypfree (gk);
				gwypfree (M);
				gwypfree (tmp);
				gwypfree (tmp2);
				gwypfree (tmp3);
				gwypdone ();
				*res = FALSE;
				return FALSE;
			}
		}
		else {	// Pocklington test ; compute the gcd's
                    sprintf (buf, "Computing GCD'S...");
                    title (buf);
                    sprintf (buf, "%s may be prime, trying to compute gcd's\n", str);
                    OutputStr (buf);
                    for (j=jmax; j>=jmin; j--) {
                        if (bpf[j] == 1)
                            // base prime factor already tested
                            continue;
                        gwypcopy (y, x);    // Computing a^((N-1)/q)
                        bbits = bpc[j];
                        ulone = 1;
                        while (!(bbits & mask)) {
                            bbits <<= 1;
                            ulone <<= 1;
                        }
                        bbits <<= 1;
                        ulone <<= 1;
                        while (ulone) {
                            if (cufftonly)
                                gwypsquare (x);
                            else if (zp || generic)
                                cuda_gwypsquare (x, 3);
                            else if(bit==1 || it==0)  
                                {cuda_gwypsquare (x,1);it=1;}
                            else  if(bit != (lasterr_point-1)) 
                                cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (x,2);
                            CHECK_IF_ANY_ERROR (x, (explen), explen, 6);
                            if (bbits & mask) {
                                if (cufftonly)
                                    gwypmul (y, x);
                                else
                                    cuda_gwypmul (y, x, 3);
                                CHECK_IF_ANY_ERROR (x, (explen), explen, 6);
                            }
                            bbits <<= 1;
                            ulone <<= 1;
                        }
                        gwyptogiant (x, tmp);
                        if (isone (tmp)) {
                            if (frestart)
                                continue;
                            if (a==2)
                                // Choose prime bases to have less restarts...
                                newa = 3;
                            else {
                                if (!(a&1))
                                    newa = a + 1;
                                else
                                    newa = a + 2;
                                while (!isPrime(newa))
                                    newa += 2;
                            }
                            nrestarts++;
                            if (nrestarts > maxrestarts) {
                                sprintf (buf, "%s may be prime, but N divides %ld^((N-1)/%lu))-1, giving up after %lu restarts...", str, a, bpf[j], maxrestarts);
                                frestart = FALSE;
                                *res = FALSE;   // Not proven prime...
                            }
                            else {
                                sprintf (buf, "%s may be prime, but N divides %ld^((N-1)/%lu))-1, restarting with a=%lu", str, a, bpf[j], newa);
                                a = newa;
                                IniWriteInt (INI_FILE, (char*)"NRestarts", nrestarts);
						IniWriteInt (INI_FILE, (char*)"FermatBase", a);
						frestart = TRUE;
                            }
                        }
                        else {
                            iaddg (-1, tmp);
                            gcdg (N, tmp);
                            if (isone (tmp)) {
                                sprintf (buf, "%ld^((N-1)/%lu)-1 is coprime to N!\n", a, bpf[j]);
                                OutputStr (buf);
                                if (verbose)
                                    writeResults (buf);
                                bpf[j] = 1;
                                // success for this prime factor of the base, continue
                            }
                            else {
                                *res = FALSE;   /* Not a prime */
                                sprintf (buf, "%s is not prime, although %ld Fermat PSP!!.", str, a);
                                break;  // No need to continue...
                            }
                        }
                    }
                    will_try_larger_fft = FALSE;	
                    if (*res && !frestart)
                        sprintf (buf, "%s is prime! (%d decimal digits)", str, nbdg);
		}
	}
	if (!frestart) {
		gwypfree (N);
		gwypfree (gk);
	}
	gwypfree (M);
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (tmp3);
	if (incr == +1) {
		gwypfree (x);
		gwypfree (y);
	}

#if defined(WIN32) && !defined(_CONSOLE)

	sprintf (buf+strlen(buf), "  Time : "); 
//	ReplaceableLine (2);	/* Replace line */ 

#else

	clearline(100);

#ifdef _CONSOLE
	OutputBoth(buf);
#else
	if (*res) {
		OutputStr((char*)"\033[7m");
		OutputBoth(buf);
		OutputStr((char*)"\033[0m");
	}
	else
		OutputBoth(buf);
#endif

	sprintf (buf, "  Time : "); 

#endif

/* Output the final timings */

	gwypend_timer (1);
//	sprintf (buf+strlen(buf)-1, "  Time: ");
	gwypwrite_timer (buf+strlen(buf), 1, TIMER_CLR | TIMER_NL); 
	if (!frestart) {
		OutputBoth (buf);
		IniWriteString (INI_FILE, (char*)"NRestarts", NULL);
		if (incr == 1)
			IniWriteString (INI_FILE, (char*)"FermatBase", NULL);
		else
			IniWriteString (INI_FILE, (char*)"LucasBaseP", NULL);
	}
	else {
		OutputStr (buf);
		if (verbose && (incr == 1))
			writeResults (buf);
	}

/* Cleanup and return */

	gwypdone ();
	_unlink (filename);
	lasterr_point = 0;
	if (frestart)
		goto restart;
//	gwypdone ();
	if (IniGetInt(INI_FILE, (char*)"PRPdone", 0))
		IniWriteString(INI_FILE, (char*)"PRPdone", NULL);
	IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
	return (TRUE);

/* An error occured, sleep, then try restarting at last save point. */

error:
	gwypfree (M);
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (tmp3);
	gwypfree (x);
	gwypfree (y);
//	gwypdone ();
	*res = FALSE;

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
		sprintf (buf, ERRMSG5, checknumber, str);
		OutputBoth (buf);
		gwypfree (N);
		gwypfree (gk);
		gwypdone ();
		_unlink (filename);
		if (IniGetInt(INI_FILE, (char*)"PRPdone", 0))
			IniWriteString(INI_FILE, (char*)"PRPdone", NULL);
		will_try_larger_fft = FALSE;
		IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
		return (TRUE);
	}

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) { 
		gwypdone ();
		will_try_larger_fft = FALSE;
		return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
		OutputBoth (ERRMSG8);
		IniWriteInt(INI_FILE, (char*)"FFT_Increment", IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
		_unlink (filename);
		will_try_larger_fft = FALSE;
	}
	goto restart;

}

int gLucasequence (
	giant modulus,
	giant exponent,
	unsigned long P,
	giant gb,
	int jmin,
	int jmax,
	char *str,
	char *buf,
	int	*res)
{
	unsigned long bit, iters, D, explen;
	unsigned long frestart=FALSE;
	gwypnum x, y, gwinvD, gwinv2;
	gwypnum a11, a12, a21, a22, b11, b12, b21, b22, c11, c12, c21, c22, p, pp, s;
	giant	tmp, tmp2;
	char	filename[20], fft_desc[256];
	long	write_time = DISK_WRITE_TIME * 60;
	long	j;
	int	echk, saving, stopping;
	time_t	start_time, current_time;
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;


//	maxrestarts = IniGetInt(INI_FILE, (char*)"MaxRestarts", 10);


restart:

	/* Get the current time */

	time (&start_time);

/* Allocate memory */

	x = gwypalloc ();
	y = gwypalloc ();
	gwinvD = gwypalloc ();
	gwinv2 = gwypalloc ();

	a11 = gwypalloc ();		// allocate memory for matrix computing...
	a12 = gwypalloc ();
	a21 = gwypalloc ();
	a22 = gwypalloc ();
	b11 = gwypalloc ();
	b12 = gwypalloc ();
	b21 = gwypalloc ();
	b22 = gwypalloc ();
	c11 = gwypalloc ();
	c12 = gwypalloc ();
	c21 = gwypalloc ();
	c22 = gwypalloc ();
	p = gwypalloc ();
	pp = gwypalloc ();
	s = gwypalloc ();

	tmp = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);
	tmp2 = newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);

/* Init, */


	gtog (exponent, tmp2);
	divg (gb, tmp2);	// tmp2 = exponent/base

	Nlen = bitlen (tmp2);

/* Init filename */

	tempFileName (filename, 'L', N);

/* Optionally resume from save file and output a message */
/* indicating we are resuming a test */

	if (fileExists (filename) && readFromFileB (filename, &bit, &P, &nrestarts, bpf, x, y)) {
		char	fmt_mask[80];
		double	pct;
		pct = trunc_percent (bit * 100.0 / Nlen);
		sprintf (fmt_mask,
			 "Resuming Lucas sequence at bit %%ld [%%.%df%%%%]\n",
			 PRECISION);
		sprintf (buf, fmt_mask, bit, pct);
		OutputStr (buf);
		if (verbose)
			writeResults (buf);	
		D = P*P - 4;
	}

/* Otherwise, output a message indicating we are starting test */

	else {
		OutputStr (buf);
		if (verbose)
			writeResults (buf);	
		D = P*P - 4;
		bit = 1;
//		dbltogw (2.0, x);	// Initial values
//		dbltogw ((double)P, y);
		itogwyp (2, x);		// Initial values
		itogwyp (P, y);
	}


	itog (D, tmp);						// Compute inverse of D modulo N
	invg (modulus, tmp);
	gianttogwyp (tmp, gwinvD);	// Convert it to gwypnum
	gtog (modulus, tmp);		// Compute inverse of 2 modulo N
	iaddg (1, tmp);			// Which is (N+1)/2
	gshiftright (1, tmp);
	gianttogwyp (tmp, gwinv2);	// Convert it to gwypnum

/* Output a message about the FFT length */

	gwypfft_description (fft_desc);
#ifdef WIN32
	sprintf (buf, "%s, P = %d\n", fft_desc, (int)P);
#else
	sprintf (buf, "%s, P = %d", fft_desc, (int)P);
#endif
	OutputStr (buf);
	LineFeed ();
	if (verbose) {
#if !defined(WIN32) 
		strcat (buf, "\n");
#endif
		writeResults (buf);
	}
	ReplaceableLine (1);	/* Remember where replaceable line is */

/* Init the title */

	title ((char*)"Computing Lucas sequence...");

/* Main loop... */

	iters = 0;
	gwypstart_timer (0);
	gwypstart_timer (1);
	will_try_larger_fft = FALSE;
	while (bit <= Nlen) {

/* Error check the last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		stopping = stopCheck ();
		echk = stopping || ERRCHK || (bit <= 50) || (bit >= Nlen-50);
		if (((bit & 127) == 0) || (bit == 1) || (bit == (lasterr_point-1))) {
			echk = 1;
			time (&current_time);
			saving = ((current_time - start_time > write_time) || (bit == 1) || (bit == (lasterr_point-1)));
		} else
			saving = 0;

/* Process this bit */


		gwypsetnormroutine (0, echk, 0);

		if (bitval (tmp2, Nlen-bit)) {
			gwypsetaddin (-(int)P);
			if (/*(bit+26 < Nlen) && (bit > 26) &&*/
				((bit != lasterr_point) || (!maxerr_recovery_mode[1] && !maxerr_recovery_mode[2]))) {
                            if (cufftonly)
                                gwypmul (y, x);
                            else
                                cuda_gwypmul (y, x, 3);
                            will_try_larger_fft = FALSE;
			}
			else {
				gwypmul_carefully (y, x);
				if (bit == lasterr_point)
					maxerr_recovery_mode[1] = FALSE;
				will_try_larger_fft = TRUE;
			}
			CHECK_IF_ANY_ERROR(x, (bit), Nlen, 1)
			gwypsetaddin (-2);
			if (/*(bit+26 < Nlen) && (bit > 26) && */
				((bit != lasterr_point) || !maxerr_recovery_mode[2])) {
                            if (cufftonly)
                                gwypsquare (y);
                            else if (zp || generic)
                                cuda_gwypsquare (y, 3);
                            else if(bit==1 || it==0)  
                                {cuda_gwypsquare (y,1);it=1;}
                            else  if(bit != (lasterr_point-1)) 
                                cuda_gwypsquare (y,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (y,2);
                            will_try_larger_fft = FALSE;
			}
			else {
				gwypsquare_carefully (y);
				if (bit == lasterr_point)
					maxerr_recovery_mode[2] = FALSE;
				will_try_larger_fft = TRUE;
			}
			CHECK_IF_ANY_ERROR(y, (bit), Nlen, 2)
		}
		else {
			gwypsetaddin (-(int)P);
			if (/*(bit+26 < Nlen) && (bit > 26) && */
				((bit != lasterr_point) || (!maxerr_recovery_mode[3] && !maxerr_recovery_mode[4]))) {
                            if (cufftonly)
                                gwypmul (x, y);
                            else
                                cuda_gwypmul (x, y, 3);
                            will_try_larger_fft = FALSE;
			}
			else {
                            gwypmul_carefully (x, y);
                            if (bit == lasterr_point)
                                    maxerr_recovery_mode[3] = FALSE;
                            will_try_larger_fft = TRUE;
			}
			CHECK_IF_ANY_ERROR(y, (bit), Nlen, 3)
			gwypsetaddin (-2);
			if (/*(bit+26 < Nlen) && (bit > 26) &&*/
				((bit != lasterr_point) || !maxerr_recovery_mode[4])) {
                            if (cufftonly)
                                gwypsquare (x);
                            else if (zp || generic)
                                cuda_gwypsquare (x, 3);
                            else if(bit==1 || it==0)  
                                {cuda_gwypsquare (x,1);it=1;}
                            else  if(bit != (lasterr_point-1)) 
                                cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (x,2);
                            will_try_larger_fft = FALSE;
			}
			else {
				gwypsquare_carefully (x);
				if (bit == lasterr_point)
					maxerr_recovery_mode[4] = FALSE;
				will_try_larger_fft = TRUE;
			}
			CHECK_IF_ANY_ERROR(x, (bit), Nlen, 4)
		}

 /* That iteration succeeded, bump counters */

		if (will_try_larger_fft && (bit == lasterr_point))
			saving = 1;
// Be sure to restart after this recovery iteration!
		will_try_larger_fft = FALSE;
		bit++;
		iters++;

/* Print a message every so often */

		if (bit % ITER_OUTPUT == 0) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (bit * 100.0 / Nlen);
			sprintf (fmt_mask, "%%.%df%%%% of %%ld", PRECISION);
			sprintf (buf, fmt_mask, pct, Nlen);
			title (buf);
			ReplaceableLine (2);	/* Replace line */
			sprintf (fmt_mask,
				 "%%s, bit: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, str, bit, Nlen, pct);
			OutputStr (buf);
			if (ERRCHK && bit > 30) {
				OutputStr ((char*)".  Round off: ");
				sprintf (buf, "%10.10f", reallyminerr);
				OutputStr (buf);
				sprintf (buf, " to %10.10f", reallymaxerr);
				OutputStr (buf);
			}
			gwypend_timer (0);
			if (CUMULATIVE_TIMING) {
				OutputStr ((char*)".  Time thusfar: ");
			} else {
				OutputStr ((char*)".  Time per bit: ");
				gwypdivide_timer (0, iters);
				iters = 0;
			}
			gwypprint_timer (0, TIMER_NL | TIMER_OPT_CLR);
			gwypstart_timer (0);
		}

/* Print a results file message every so often */

		if (bit % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
			sprintf (buf, "Bit %ld / %ld\n", bit, Nlen);
			writeResults (buf);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving || stopping) {
			write_time = DISK_WRITE_TIME * 60;
			if (! writeToFileB (filename, bit, P, nrestarts, bpf, x, y)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (buf);
				if (write_time > 600) write_time = 600;
			}	
			time (&start_time);

/* If an escape key was hit, write out the results and return */

			if (stopping) {
				gwypfree (tmp);
				gwypfree (tmp2);
				gwypfree (x);				// Clean up
				gwypfree (y);
				gwypfree (gwinvD);
				gwypfree (gwinv2);
				gwypfree (a11);			// Clean up the matrix
				gwypfree (a12);
				gwypfree (a21);
				gwypfree (a22);
				gwypfree (b11);
				gwypfree (b12);
				gwypfree (b21);
				gwypfree (b22);
				gwypfree (c11);
				gwypfree (c12);
				gwypfree (c21);
				gwypfree (c22);
				gwypfree (p);
				gwypfree (pp);
				gwypfree (s);
//				gwypdone ();
				Nlen = bitlen (N);
				*res = FALSE;
				return (FALSE);
			}
		}

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next iteration so that we can compare our */
/* residues to programs that start counter at zero or one. */

		if (interimResidues && bit % interimResidues < 2) {
			gwyptogiant (x, tmp);		// The modulo reduction is done here
			modg (N, tmp);				// External modulus and gwypnums one may be different...
			if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
			else if (abs(tmp->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
			else if (abs(tmp->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
			sprintf (buf, "%s interim residue %s at bit %ld\n", str, res64, bit);
			OutputBoth (buf);
		}

/* Write a save file every "interimFiles" iterations. */

		if (interimFiles && bit % interimFiles == 0) {
			char	interimfile[20];
			sprintf (interimfile, "%.8s.%03lu",
				 filename, bit / interimFiles);
			if (! writeToFileB (interimfile, bit, P, nrestarts, bpf, x, y)) {
				sprintf (buf, WRITEFILEERR, interimfile);
				OutputBoth (buf);
			}
		}
	}

	will_try_larger_fft = TRUE;	// All following errors are considered unrecoverable...

	// Compute the matrix at (N+1)/base

	gwypsetaddin (0);		// Reset addin constant.
	gwypcopy (x, a22);	        // a22 = V((N+1)/base)
	gwypcopy (y, a12);		// a12 = V((N+1)/base+1)
	gwypcopy (y, a21);		// a21 = V((N+1)/base+1)
	gwypcopy (y, a11);		// a11 = V((N+1)/base+1)
	gwypcopy (x, y);		// Now, y = V((N+1)/base)
	gwypsetnormroutine (0, 1, 1);	// set mul. by const.
	gwypsetmulbyconst (2);
//	gwypmul (gwinvD, a21);		
// a21 = D^-1*2*V((N+1)/base+1) modulo N
        if (cufftonly)
            gwypmul (gwinvD, a21);
        else
            cuda_gwypmul (gwinvD, a21, 3);
	CHECK_IF_ANY_ERROR(a21, (Nlen), Nlen, 6)
	gwypsetmulbyconst (P);
//	gwypmul (gwinvD, x);
// x =  D^-1*P*V((N+1)/base) modulo N
        if (cufftonly)
            gwypmul (gwinvD, x);
        else
            cuda_gwypmul (gwinvD, x, 3);
	CHECK_IF_ANY_ERROR(x, (Nlen), Nlen, 6)
	gwypsub (x, a21);
// a21 = D^-1*(2*V((N+1)/base+1)-P*V((N+1)/base)) = U(N+1)/base modulo N
//	gwypmul (gwinvD, a11); a11 = D^-1*P*V((N+1)/base+1) modulo N
        if (cufftonly)
            gwypmul (gwinvD, a11);
        else
            cuda_gwypmul (gwinvD, a11, 3);
	CHECK_IF_ANY_ERROR(a11, (Nlen), Nlen, 6)
	gwypsetmulbyconst (2);
//	gwypmul (gwinvD, y);// xx = D^-1*2*V((N+1)/base) modulo N
        if (cufftonly)
            gwypmul (gwinvD, y);
        else
            cuda_gwypmul (gwinvD, y, 3);
	CHECK_IF_ANY_ERROR(y, (Nlen), Nlen, 6)
	gwypsub (y, a11);
// a11 = D^-1*(P*V((N+1)/base+1)-2*V((N+1)/base)) = U((N+1)/base+1) modulo N
	gwypsetnormroutine (0, 1, 0);	// reset mul by const
//	gwypmul (gwinv2, a22); a22 = 2^-1*V((N+1)/base)
        if (cufftonly)
            gwypmul (gwinv2, a22);
        else
            cuda_gwypmul (gwinv2, a22, 3);
	CHECK_IF_ANY_ERROR(a22, (Nlen), Nlen, 6)
//	gwypmul (gwinv2, a12); a12 = 2^-1*V((N+1)/base+1)
        if (cufftonly)
            gwypmul (gwinv2, a12);
        else
            cuda_gwypmul (gwinv2, a12, 3);
	CHECK_IF_ANY_ERROR(a12, (Nlen), Nlen, 6)
	gwypsetmulbyconst (P);
	gwypsetnormroutine (0, 1, 1);	// set mul. by const.
	gwypcopy (a11, x);		// x = U((N+1)/base+1)
//	gwypmul (gwinv2, x);		// x = 2^-1*P*U((N+1)/base+1)
        if (cufftonly)
            gwypmul (gwinv2, x);
        else
            cuda_gwypmul (gwinv2, x, 3);
	CHECK_IF_ANY_ERROR(x, (Nlen), Nlen, 6)
	gwypsub (x, a12);
// a12 = 2^-1(V((N+1)/base+1)-P*U((N+1)/base+1))
	gwypcopy (a21, x);		// x = U((N+1)/base)
//	gwypmul (gwinv2, x);		// x = 2^-1*P*U((N+1)/base)
        if (cufftonly)
            gwypmul (gwinv2, x);
        else
            cuda_gwypmul (gwinv2, x, 3);
	CHECK_IF_ANY_ERROR(x, (Nlen), Nlen, 6)
	gwypsub (x, a22);  // a22 = 2^-1(V((N+1)/base)-P*U((N+1)/base))
	gwypsetnormroutine (0, 1, 0);	// reset mul by const

//	gwyptogiant (a21, tmp);		// tmp = U((N+1)/base) modulo N

	gwypcopy (a11, c11);	// Save the current matrix
	gwypcopy (a12, c12);
	gwypcopy (a21, c21);
	gwypcopy (a22, c22);

	gwypcopy (a11, b11);	// Copy the current matrix
	gwypcopy (a12, b12);
	gwypcopy (a21, b21);
	gwypcopy (a22, b22);

        explen = bitlen (gb);
        bit = 1;
	while (bit < explen) {	// Finish to compute U(N+1)

// Square the matrix

		gwypcopy (a12, p);		// a12-->p
		gwypcopy (a12, pp);		// a12-->pp
		gwypadd3 (a11, a22, s);		// a11+a12-->s
//		gwypmul (a21, p);		// a21*a12-->p
                if (cufftonly)
                    gwypmul (a21, p);
                else
                    cuda_gwypmul (a21, p, 3);
		CHECK_IF_ANY_ERROR(p, (Nlen), Nlen, 6)
//		gwypsquare (a22);	        // a22*a22-->a22
                if (cufftonly)
                    gwypsquare (a22);
                else
                    cuda_gwypsquare (a22, 3);
		CHECK_IF_ANY_ERROR(a22, (Nlen), Nlen, 6)
//		gwypmul (s, a21);	        // (a11+a22)*a21-->a21 T
                if (cufftonly)
                    gwypmul (s, a21);
                else
                    cuda_gwypmul (s, a21, 3);
		CHECK_IF_ANY_ERROR(a21, (Nlen), Nlen, 6)
		gwypadd (p, a22);     // a21*a12+a22*a22-->a22 T
//		gwypmul (s, a12);     // (a11+a22)*a12-->a12 T
                if (cufftonly)
                    gwypmul (s, a12);
                else
                    cuda_gwypmul (s, a12, 3);
		CHECK_IF_ANY_ERROR(a12, (Nlen), Nlen, 6)
//		gwypsquare (a11);     // a11*a11-->a11
                if (cufftonly)
                    gwypsquare (a11);
                else
                    cuda_gwypsquare (a11, 3);
		CHECK_IF_ANY_ERROR(a11, (Nlen), Nlen, 6)
		gwypadd (p, a11);     // a21*a12+a11*a11-->a11 T

// Multiply it if required

		if (bitval (gb, explen-bit-1)) {
			gwypcopy (a11, p);		// a11-->p
			gwypcopy (a21, pp);		// a21-->pp
//			gwypmul (b11, a11);	        // b11*a11-->a11
                        if (cufftonly)
                            gwypmul (b11, a11);
                        else
                            cuda_gwypmul (b11, a11, 3);
			CHECK_IF_ANY_ERROR(a11, (Nlen), Nlen, 6)
//			gwypmul (b12, pp);		// b12*a21-->pp
                        if (cufftonly)
                            gwypmul (b12, pp);
                        else
                            cuda_gwypmul (b12, pp, 3);
			CHECK_IF_ANY_ERROR(pp, (Nlen), Nlen, 6)
			gwypadd (pp, a11);   // b11*a11+b12*a21-->a11 T
//			gwypmul (b21, p);    // b21*a11-->p
                        if (cufftonly)
                            gwypmul (b21, p);
                        else
                            cuda_gwypmul (b21, p, 3);
			CHECK_IF_ANY_ERROR(p, (Nlen), Nlen, 6)
//			gwypmul (b22, a21);  // b22*a21-->a21
                        if (cufftonly)
                            gwypmul (b22, a21);
                        else
                            cuda_gwypmul (b22, a21, 3);
			CHECK_IF_ANY_ERROR(a21, (Nlen), Nlen, 6)
			gwypadd (p, a21);		// b21*a11+b22*a21-->a21 T
			gwypcopy (a12, p);		// a12-->p
			gwypcopy (a22, pp);		// a22-->pp
//			gwypmul (b11, a12);	        // b11*a12-->a12
                        if (cufftonly)
                            gwypmul (b11, a12);
                        else
                            cuda_gwypmul (b11, a12, 3);
			CHECK_IF_ANY_ERROR(a12, (Nlen), Nlen, 6)
//			gwypmul (b12, pp);		// b12*a22-->pp
                        if (cufftonly)
                            gwypmul (b12, pp);
                        else
                            cuda_gwypmul (b12, pp, 3);
			CHECK_IF_ANY_ERROR(pp, (Nlen), Nlen, 6)
			gwypadd (pp, a12);		// b11*a12+b12*a22-->a12 T
//			gwypmul (b21, p);		// b21*a12-->p
                        if (cufftonly)
                            gwypmul (b21, p);
                        else
                            cuda_gwypmul (b21, p, 3);
			CHECK_IF_ANY_ERROR(p, (Nlen), Nlen, 6)
//			gwypmul (b22, a22);	        // b22*a22-->a22
                        if (cufftonly)
                            gwypmul (b22, a22);
                        else
                            cuda_gwypmul (b22, a22, 3);
			CHECK_IF_ANY_ERROR(a22, (Nlen), Nlen, 6)
			gwypadd (p, a22);    // b21*a12+b22*a22-->a22 T

		}
		bit++;
	}

	clearline (100);

	gwyptogiant (a21, tmp);	// tmp = U(N+1) modulo N
	ReplaceableLine (2);	/* Replace line */

	if (!isZero (tmp)) {
//	if (!gwiszero(a21)) {
		*res = FALSE;	/* Not a prime */
		if (abs(tmp->sign) < 2)	// make a 64 bit residue correct !!
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
		else if (abs(tmp->sign) < 3)
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
		else if (abs(tmp->sign) < 4)
			sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
		else
			sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
		if (IniGetInt(INI_FILE, (char*)"Verify", 0))
			sprintf (buf, "%s is not prime. P = %lu, Lucas RES64: %s", str, P, res64);
		else
			sprintf (buf, "%s is not prime, although Fermat PSP! P = %lu, Lucas RES64: %s", str, P, res64);
	}
	else {
		sprintf (buf, "%s may be prime, trying to compute gcd's\n", str);
		if (verbose)
			OutputBoth (buf);
		else
			OutputStr (buf);
		for (j=jmax; j>=jmin; j--) {
			if (bpf[j] == 1)
                            // base prime factor already tested
				continue;
			gwypcopy (c11, a11);			// Computing U((N+1)/q)
			gwypcopy (c12, a12);
			gwypcopy (c21, a21);
			gwypcopy (c22, a22);
			explen = bitlen (gbpc[j]);
                        bit = 1;
			while (bit < explen) {

// Square the matrix

				gwypcopy (a12, p);	// a12-->p
				gwypcopy (a12, pp);	// a12-->pp
				gwypadd3 (a11, a22, s);	// a11+a12-->s
//				gwypmul (a21, p);	// a21*a12-->p
                                if (cufftonly)
                                    gwypmul (a21, p);
                                else
                                    cuda_gwypmul (a21, p, 3);
				CHECK_IF_ANY_ERROR(p, (Nlen), Nlen, 6)
//				gwypsquare (a22);	// a22*a22-->a22
                                if (cufftonly)
                                    gwypsquare (a22);
                                else
                                    cuda_gwypsquare (a22, 3);
				CHECK_IF_ANY_ERROR(a22, (Nlen), Nlen, 6)
//				gwypmul (s, a21); (a11+a22)*a21-->a21 T
                                if (cufftonly)
                                    gwypmul (s, a21);
                                else
                                    cuda_gwypmul (s, a21, 3);
				CHECK_IF_ANY_ERROR(a21, (Nlen), Nlen, 6)
				gwypadd (p, a22);
                                // a21*a12+a22*a22-->a22 T
//				gwypmul (s, a12); (a11+a22)*a12-->a12 T
                                if (cufftonly)
                                    gwypmul (s, a12);
                                else
                                    cuda_gwypmul (s, a12, 3);
				CHECK_IF_ANY_ERROR(a12, (Nlen), Nlen, 6)
//				gwypsquare (a11); a11*a11-->a11
                                if (cufftonly)
                                    gwypsquare (a11);
                                else
                                    cuda_gwypsquare (a11, 3);
				CHECK_IF_ANY_ERROR(a11, (Nlen), Nlen, 6)
				gwypadd (p, a11);
                                // a21*a12+a11*a11-->a11 T

// Multiply it if required

				if (bitval (gbpc[j], explen-bit-1)) {
					gwypcopy (a11, p);// a11-->p
					gwypcopy (a21, pp);// a21-->pp
//				gwypmul (b11, a11); b11*a11-->a11
                                        if (cufftonly)
                                            gwypmul (b11, a11);
                                        else
                                            cuda_gwypmul (b11, a11, 3);
					CHECK_IF_ANY_ERROR(a11, (Nlen), Nlen, 6)
//					gwypmul (b12, pp); b12*a21-->pp
                                        if (cufftonly)
                                            gwypmul (b12, pp);
                                        else
                                            cuda_gwypmul (b12, pp, 3);
					CHECK_IF_ANY_ERROR(pp, (Nlen), Nlen, 6)
					gwypadd (pp, a11);
                                        // b11*a11+b12*a21-->a11 T
//					gwypmul (b21, p); b21*a11-->p
                                        if (cufftonly)
                                            gwypmul (b21, p);
                                        else
                                            cuda_gwypmul (b21, p, 3);
					CHECK_IF_ANY_ERROR(p, (Nlen), Nlen, 6)
//				gwypmul (b22, a21); b22*a21-->a21
                                        if (cufftonly)
                                            gwypmul (b22, a21);
                                        else
                                            cuda_gwypmul (b22, a21, 3);
					CHECK_IF_ANY_ERROR(a21, (Nlen), Nlen, 6)
					gwypadd (p, a21);
                                        // b21*a11+b22*a21-->a21 T
					gwypcopy (a12, p);// a12-->p
					gwypcopy (a22, pp);// a22-->pp
//				gwypmul (b11, a12); b11*a12-->a12
                                        if (cufftonly)
                                            gwypmul (b11, a12);
                                        else
                                            cuda_gwypmul (b11, a12, 3);
					CHECK_IF_ANY_ERROR(a12, (Nlen), Nlen, 6)
//					gwypmul (b12, pp); b12*a22-->pp
                                        if (cufftonly)
                                            gwypmul (b12, pp);
                                        else
                                            cuda_gwypmul (b12, pp, 3);
					CHECK_IF_ANY_ERROR(pp, (Nlen), Nlen, 6)
					gwypadd (pp, a12);
                                        // b11*a12+b12*a22-->a12 T
//					gwypmul (b21, p); b21*a12-->p
                                        if (cufftonly)
                                            gwypmul (b21, p);
                                        else
                                            cuda_gwypmul (b21, p, 3);
					CHECK_IF_ANY_ERROR(p, (Nlen), Nlen, 6)
//				gwypmul (b22, a22); b22*a22-->a22
                                        if (cufftonly)
                                            gwypmul (b22, a22);
                                        else
                                            cuda_gwypmul (b22, a22, 3);
					CHECK_IF_ANY_ERROR(a22, (Nlen), Nlen, 6)
					gwypadd (p, a22);
                                        // b21*a12+b22*a22-->a22 T

				}
				bit++;
			}
			gwyptogiant (a21, tmp);
			if (isZero (tmp)) {
                            gtoc(gbpf[j], bpfstring, strlen(sgb));
				sprintf (buf, "%s may be prime, but N divides U((N+1)/%s), P = %lu\n", str, bpfstring, P);
				if (verbose)
					OutputBoth (buf);
				else
					OutputStr (buf);
				frestart = TRUE;
				_unlink (filename);
				continue;
//				break;
			}
			else {
				gcdg (modulus, tmp);
				if (isone (tmp)) {
                                    gtoc(gbpf[j], bpfstring, strlen(sgb));
                                    sprintf (buf, "U((N+1)/%s) is coprime to N!\n", bpfstring);
                                    OutputStr (buf);
                                    if (verbose)
                                        writeResults (buf);	
                                    bpf[j] = 1;
				}
				else {
                                    *res = FALSE; /* Not a prime */
                                    if (IniGetInt(INI_FILE, (char*)"Verify", 0))
                                        sprintf (buf, "%s is not prime, although Lucas PSP!! (P = %lu)", str, P);
                                    else
                                        sprintf (buf, "%s is not prime, although Fermat and Lucas PSP!! (P = %lu)", str, P);
					break;
				}
			}
		}
	}

	if (*res && !frestart)
		sprintf (buf, "%s is prime! (%d decimal digits, P = %lu)", str, nbdg, P);

	will_try_larger_fft = FALSE;// Reset the "unrecoverable" condition.
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (x);			// Clean up
	gwypfree (y);
	gwypfree (gwinvD);
	gwypfree (gwinv2);
	gwypfree (a11);		// Clean up the matrix
	gwypfree (a12);
	gwypfree (a21);
	gwypfree (a22);
	gwypfree (b11);
	gwypfree (b12);
	gwypfree (b21);
	gwypfree (b22);
	gwypfree (c11);
	gwypfree (c12);
	gwypfree (c21);
	gwypfree (c22);
	gwypfree (p);
	gwypfree (pp);
	gwypfree (s);

/* Cleanup and return */

	Nlen = bitlen (N);
//	gwypdone ();
	_unlink (filename);
	lasterr_point = 0;

	if (frestart)
		return -2;
	else {
		IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
		return TRUE;
	}

/* An error occured, sleep, then try restarting at last save point. */

error:
	Nlen = bitlen (N);
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (x);				// Clean up
	gwypfree (y);
	gwypfree (gwinvD);
	gwypfree (gwinv2);
	gwypfree (a11);			// Clean up the matrix
	gwypfree (a12);
	gwypfree (a21);
	gwypfree (a22);
	gwypfree (b11);
	gwypfree (b12);
	gwypfree (b21);
	gwypfree (b22);
	gwypfree (c11);
	gwypfree (c12);
	gwypfree (c21);
	gwypfree (c22);
	gwypfree (p);
	gwypfree (pp);
	gwypfree (s);
	*res = FALSE;

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
		sprintf (buf, ERRMSG5, checknumber, str);
		OutputBoth (buf);
		will_try_larger_fft = FALSE;
		_unlink (filename);
		return (FALSE);
	}

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
		will_try_larger_fft = FALSE;
		return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
		OutputBoth (ERRMSG8);
		IniWriteInt(INI_FILE, (char*)"FFT_Increment", IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
		_unlink (filename);
		will_try_larger_fft = FALSE;
		return (-1);
	}
	goto restart;
}


int gplusminustest ( 
	char *sgk,
	char *sgb,
	giant gb,
	unsigned long n, 
	int incr,
	unsigned long shift,
	int	*res) 
{ 
	char	filename[20], buf[sgkbufsize+256], str[sgkbufsize+256], sgk1[sgkbufsize], fft_desc[256], oldres64[17]; 
	unsigned long bits, explen, iters, bit, frestart=FALSE;
	unsigned long newa, maxrestarts;
//	uint32_t hi = 0, lo = 0;
	double dk;
	giant grem, tmp, tmp2, tmp3;
	gwypnum x, y;
	long	a, P, write_time = DISK_WRITE_TIME * 60;
	int	echk, saving, stopping, jmin, jmax, j, retval, Psample;
        int     factorized, factorized_part = 0;
	time_t	start_time, current_time;
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;

/*	gk = newgiant (strlen(sgk)/2 + 8);
// Allocate one byte per decimal digit + spares
	ctog (sgk, gk);		// Convert k string to giant
	gshiftleft (shift, gk);	// Shift k multiplier if requested
	gtoc (gk, sgk1, sgkbufsize);   // Updated k string
	if (!strcmp(sgk1, "1"))
		sprintf (str, "%lu^%lu%c%d", base, n, (incr < 0) ? '-' : '+', abs(incr));
	else
		sprintf (str, "%s*%lu^%lu%c%d", sgk1, base, n, (incr < 0) ? '-' : '+', abs(incr));

	klen = bitlen(gk);
*/
	if ((gformat == ABCDN) || (gformat == ABCDNG)) {
            // Compute gk = gb^(n-m)-1
            bits = ndiff*bitlen (gb);
            gk = newgiant ((bits >> 2) + 8);
            gtog (gb, gk);
            power (gk, ndiff);
            iaddg (-1, gk);
            sprintf (str, "%s^%lu-%s^%lu%c%d", sgb, n+ndiff, sgb, n, incr < 0 ? '-' : '+', abs(incr));
	}
	else {
            gk = newgiant (strlen(sgk)/2 + 8);
            // Allocate one byte per decimal digit + spares
            ctog (sgk, gk);	// Convert k string to giant
            grem = newgiant (2*abs(gk->sign) + 8);
            // place for mod (gk, gb)
            gshiftleft (shift, gk); // Shift k multiplier if requested
            gtoc (gk, sgk1, sgkbufsize);	// Updated k string
            if (!strcmp(sgk1, "1"))
                sprintf (str, "%s^%lu%c%d", sgb, n, incr < 0 ? '-' : '+', abs(incr));
            else
                sprintf (str, "%s*%s^%lu%c%d", sgk1, sgb, n, incr < 0 ? '-' : '+', abs(incr));
	}

	bits = n * bitlen(gb) + bitlen(gk); 
	N =  newgiant ((bits >> 2) + 8);   // Allocate memory for N

//	Be sure the base does not divide the gk multiplier :

/*	while (!(gmodi (base, gk))) {
		uldivg (base, gk);
		n++;
	}
*/
	if ((gformat != ABCDN) && (gformat != ABCDNG)) {
		while (!isone(gk)) {
			gtog (gk,grem);
			modg (gb, grem);
			if (!isZero(grem))
				break;
			divg (gb, gk);
			n++;
		}
		free (grem);
	}

	if ((int)klen != bitlen(gk))	// Has k been updated ?
		gtoc (gk, sgk1, sgkbufsize);
	else
		strcpy (sgk1, sgk);

/*	bits = (unsigned long) ((n * log(base)) / log(2) + bitlen(gk)); 
	N =  newgiant ((bits >> 2) + 16);  // Allocate memory for N
*/

//	Compute the number we are testing.
        
	gtog (gb, N);
	power (N, n);

	Nlen = bitlen (N);	// Bit length of base^n

	mulg (gk, N); 

	iaddg (incr, N);

	klen = bitlen(gk);

	if (klen > 53) {					// we must use generic reduction
		dk = 0.0;
	}
	else {								// we can use DWT ; compute the multiplier as a double
		dk = (double)gk->n[0];
		if (gk->sign > 1)
			dk += 65536.0*(double)gk->n[1];
		if (gk->sign > 2)
			dk += 65536.0*65536.0*(double)gk->n[2];
		if (gk->sign > 3)
			dk += 65536.0*65536.0*65536.0*(double)gk->n[3];
	}


	nbdg = gnbdg (N, 10);	
        // Compute the number of decimal digits of the tested number.

	if ((klen > Nlen) && (nbdg > 400)) {
//	    sprintf(buf, "%s > %lu^%lu, so, only a Strong PRP test is done for %s.\n", sgk1, base, n, str);
            if ((gformat == ABCDN) || (gformat == ABCDNG))
                sprintf(buf, "%s^%lu-1 > %s^%lu, so, only a PRP test is done for %s.\n", sgb, ndiff, sgb, n, str);
            else
                sprintf(buf, "%s > %s^%lu, so, only a PRP test is done for %s.\n", sgk, sgb, n, str);
            OutputBoth(buf);
            retval = gisPRPinternal (str, dk, gb, n, incr, res);
            gwypfree(gk);
            gwypfree(N);
            return retval;
	}

	Nlen = bitlen (N); // Bit length of N
//	findbpf (base);	   // Factorize the base
	factorized = findgbpf (gb); // Factorize the base if possible...

	for (jmax=29; (jmax>=0) && !bpf[jmax]; jmax--);
	jmin = 0;							// Choose the minimum required factored part.
	if (jmax) {							// The base is composite...
            factorized_part = bitlen (gk);
            for (j=jmax; j>0; j--) {
                if (bpf[j] == 1)
                    factorized_part += n*vpf[j]*bitlen(gbpf[j]);
                else
                    factorized_part += (unsigned long)floor(n*vpf[j]*log ((double)bpf[j])/log(2.0));
                if ((2*factorized_part) > Nlen)
                    break;
            }
            jmin = j;
            sprintf (buf, "Base factorized as : ");

            for (j=0; j<=jmax; j++) {
                if (j<jmax) {
                    if (bpf[j] == 1)
                        gtoc(gbpf[j], buf+strlen(buf), strlen(sgb));
                    else
                        sprintf (buf+strlen(buf), "%lu", bpf[j]);
                    if (vpf[j]>1)
                        sprintf (buf+strlen(buf), "^%lu*", vpf[j]);
                    else
                        sprintf (buf+strlen(buf), "*");
                }
                else {
                    if (bpf[j] == 1)
                        gtoc(gbpf[j], buf+strlen(buf), strlen(sgb));
                    else
                        sprintf (buf+strlen(buf), "%lu", bpf[j]);
                    if (vpf[j]>1)
                        sprintf (buf+strlen(buf), "^%lu\n", vpf[j]);
                    else
                        sprintf (buf+strlen(buf), "\n");
                }
            }
            if (!setuponly)
                if (verbose)
                    OutputBoth(buf);
                else
                    OutputStr(buf);
                        
	}

	sprintf (buf, "Base prime factor(s) taken : ");

	for (j=jmin; j<=jmax; j++) {
            if (j<jmax)
                if (bpf[j] == 1) {
                    gtoc(gbpf[j], buf+strlen(buf), strlen(sgb));
                    sprintf (buf+strlen(buf), ", ");
                }
                else
                    sprintf (buf+strlen(buf), "%lu, ", bpf[j]);
            else
                if (bpf[j] == 1) {
                    gtoc(gbpf[j], buf+strlen(buf), strlen(sgb));
                    if (factorized)
                        sprintf (buf+strlen(buf), "\n");
                    else
                        sprintf (buf+strlen(buf), " (Must be proven prime or factorized externally)\n");
                }
                else
                    sprintf (buf+strlen(buf), "%lu\n", bpf[j]);
	}

        if (verbose)
            OutputBoth(buf);
        else
            OutputStr(buf);

	maxrestarts = IniGetInt(INI_FILE, (char*)"MaxRestarts", 10);
	nrestarts = IniGetInt (INI_FILE, (char*)"NRestarts", 0);
	if (!(a = IniGetInt (INI_FILE, (char*)"FermatBase", 0)))
		a = IniGetInt (INI_FILE, (char*)"FBase", 3);// The base for the PRP and Pocklington tests
	if (!(P = IniGetInt (INI_FILE, (char*)"LucasBaseP", 0))) {
		Psample = genLucasBaseP (N, IniGetInt (INI_FILE, (char*)"PBase", 3));
		if (Psample < 0) {
			if (Psample == -1)
				sprintf (buf, "Cannot compute P to test %s...\nThis is surprising, please, let me know that!!\nMy E-mail is jpenne@free.fr\n", str);
			else
				sprintf (buf, "%s has a small factor : %d !!\n", str, abs(Psample));
			OutputBoth (buf);
			gwypfree(gk);
			gwypfree(N);
			return (TRUE); 
		}
		else
			P = Psample;
	}
											// The Discriminant for the Morrison test
//	D = P*P-4;								// D = P^2 - 4*Q with Q = 1


// restart:

/* Setup the gwypnum code */


restart:

	*res = TRUE;

	M = newgiant ((bits >> 2) + 16);	// Allocate memory for M
	tmp = newgiant ((bits >> 2) + 16);	// Allocate memory for tmp
	gtog (N, tmp);
	iaddg (-1, tmp);			// tmp = N-1

	gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
	if (incr == +1) {
            gwypsetmaxmulbyconst (a);
            divg (gb, tmp);		// tmp = (N-1)/base
            explen = bitlen (tmp);
            if (!setupok (gwypsetup (dk, 1, n, +1, N), N, str, res))
// set base to 1 to force generic mode...
            {
                gwypfree(gk);
                gwypfree(N);
                gwypfree (M);
                gwypfree (tmp);
//		*res = FALSE;		// Not proven prime...
                return TRUE; 
            }
            tmp2 =newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);	// Allocate memory for tmp2
            tmp3 =newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);	// Allocate memory for tmp3
	}
	else {
            gwypsetmaxmulbyconst (max(a, P));
            gtog (N, M);
            iaddg (1, M);
            explen = bitlen (tmp);
            if (!setupok (gwypsetup (dk, 1, n, -1, N), N, str, res)) {
// set base to 1 to force generic mode...
                gwypfree(gk);
                gwypfree(N);
                gwypfree (M);
                gwypfree (tmp);
//		*res = FALSE;		// Not proven prime...
                return TRUE; 
            }
            tmp2 =newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);	// Allocate memory for tmp2
            tmp3 =newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16);	// Allocate memory for tmp3
            tempFileName (filename, 'L', N);
            if (fileExists (filename)) {// Resuming a Lucas sequence...
                goto DoLucas;
            }
            else if (IniGetInt(INI_FILE, (char*)"Verify", 0) ||         IniGetInt(INI_FILE, (char*)"PRPdone", 0)) {
                gwypclear_timers ();	// Init. timers
                sprintf (buf, "Starting Lucas sequence for %s...\n", str);
                goto DoLucas;	// Starting directly a Lucas sequence...
            }
	}

/* Init filename */

	tempFileName (filename, 'z', N);

/* Allocate memory */

	x = gwypalloc ();
	y = gwypalloc ();

/* Optionally resume from save file and output a message */
/* indicating we are resuming a test */

	//cuda if (fileExists (filename) && readFromFileB (filename, &bit, &a, &nrestarts, bpf, x, y)) {
	if (fileExists (filename) && readFromFileB (filename,(unsigned long *) &bit,(unsigned long *) &a, &nrestarts, bpf, x, y)) {
		char	fmt_mask[80];
		double	pct;
		pct = trunc_percent (bit * 100.0 / explen);
		sprintf (fmt_mask,
			 "Resuming N%%c%%d prime test of %%s at bit %%ld [%%.%df%%%%]\n",
			 PRECISION);
		sprintf (buf, fmt_mask, (incr < 0) ? '+' : '-', abs(incr), str, bit, pct);
		OutputStr (buf);
		if (verbose)
			writeResults (buf);	
	}

/* Otherwise, output a message indicating we are starting test */

	else {
		if (frestart) {
			sprintf (buf, "Restarting N%c%d prime test of %s\n", (incr < 0) ? '+' : '-', abs(incr), str);
			frestart = FALSE;
		}
		else {
			gwypclear_timers ();		// Init. timers
                        if (showdigits)
                            sprintf (buf, "Starting N%c%d prime test of %s (%d decimal digits)\n", incr < 0 ? '+' : '-', abs(incr), str, nbdg);
                        else
                            sprintf (buf, "Starting N%c%d prime test of %s\n", incr < 0 ? '+' : '-', abs(incr), str);
		}
		OutputStr (buf);
		if (verbose)
			writeResults (buf);	
		bit = 1;
		itogwyp (a, x);
//		dbltogw ((double)a, x);
	}

/* Get the current time */

	gwypstart_timer (0);
	gwypstart_timer (1);
	time (&start_time);

/* Output a message about the FFT length */

	gwypfft_description (fft_desc);

#ifdef WIN32
	sprintf (buf, "%s, a = %ld\n", fft_desc, a);
#else
	sprintf (buf, "%s, a = %ld", fft_desc, a);
#endif
	OutputStr (buf);
	LineFeed ();
	if (verbose) {
#if !defined(WIN32) 
		strcat (buf, "\n");
#endif
		writeResults (buf);
	}
	ReplaceableLine (1);	/* Remember where replaceable line is */

/* Init the title */

	if (incr > 0)
		title ((char*)"Pocklington prime test in progress...");
	else
		title ((char*)"Fermat PRP test in progress...");

/* Do the PRP test */

	gwypsetmulbyconst (a);
	gwypsetaddin(0);
	iters = 0;
	
	while (bit < explen) {

/* Error check the first and last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		stopping = stopCheck ();
		echk = stopping || ERRCHK || (bit <= 50) || (bit >= explen-50);
		if (((bit & 127) == 0) || (bit == 1) || (bit == (lasterr_point-1))) {
			echk = 1;
			time (&current_time);
			saving = ((current_time - start_time > write_time) || (bit == 1) || (bit == (lasterr_point-1)));
		} else
			saving = 0;

/* Process this bit */


		if (bitval (tmp, explen-bit-1)) {
			gwypsetnormroutine (0, echk, 1);
		} else {
			gwypsetnormroutine (0, echk, 0);
		}

		if (/*(bit+25 < explen) && (bit > 25) && */((bit != lasterr_point) || !maxerr_recovery_mode[6])) {
                    if (cufftonly)
                        gwypsquare (x);
                    else if (zp || generic)
                        cuda_gwypsquare (x, 3);
                    else if(bit==1 || it==0)  
                        {cuda_gwypsquare (x,1);it=1;}
                    else  if(bit != (lasterr_point-1)) 
                        cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                    else
                        cuda_gwypsquare (x,2);
		}
		else {
			gwypsquare_carefully (x);
			will_try_larger_fft = TRUE;
			if (bit == lasterr_point)
				maxerr_recovery_mode[6] = FALSE;
		}

		CHECK_IF_ANY_ERROR (x, (bit), explen, 6);

/* That iteration succeeded, bump counters */

		if (will_try_larger_fft && (bit == lasterr_point))
			saving = 1;					// Be sure to restart after this recovery iteration!
		will_try_larger_fft = FALSE;
		bit++;
		iters++;

/* Print a message every so often */

		if (bit % ITER_OUTPUT == 0) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (bit * 100.0 / explen);
			sprintf (fmt_mask, "%%.%df%%%% of %%ld", PRECISION);
			sprintf (buf, fmt_mask, pct, explen);
			title (buf);
			ReplaceableLine (2);	/* Replace line */
			sprintf (fmt_mask,
				 "%%s, bit: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, str, bit, explen, pct);
			OutputStr (buf);
			if (ERRCHK && bit > 30) {
				OutputStr ((char*)".  Round off: ");
				sprintf (buf, "%10.10f", reallyminerr);
				OutputStr (buf);
				sprintf (buf, " to %10.10f", reallymaxerr);
				OutputStr (buf);
			}
			gwypend_timer (0);
			if (CUMULATIVE_TIMING) {
				OutputStr ((char*)".  Time thusfar: ");
			} else {
				OutputStr ((char*)".  Time per bit: ");
				gwypdivide_timer (0, iters);
				iters = 0;
			}
			gwypprint_timer (0, TIMER_NL | TIMER_OPT_CLR);
			gwypstart_timer (0);
		}

/* Print a results file message every so often */

		if (bit % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
			sprintf (buf, "Bit %ld / %ld\n", bit, explen);
			writeResults (buf);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving || stopping) {
			write_time = DISK_WRITE_TIME * 60;
			if (! writeToFileB (filename, bit, a, nrestarts, bpf, x, y)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (buf);
				if (write_time > 600) write_time = 600;
			}	
			time (&start_time);

/* If an escape key was hit, write out the results and return */

			if (stopping) {
				gwypfree (N);
				gwypfree (gk);
				gwypfree (M);
				gwypfree (tmp);
				gwypfree (tmp2);
				gwypfree (tmp3);
				gwypfree (x);
				gwypfree (y);
				gwypdone ();
				*res = FALSE;
				return (FALSE);
			}
		}

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next iteration so that we can compare our */
/* residues to programs that start counter at zero or one. */

		if (interimResidues && bit % interimResidues < 2) {
  			gwyptogiant (x, tmp2);	// The modulo reduction is done here
			if (abs(tmp2->sign) < 2)		// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp2->n[0]);
			else if (abs(tmp2->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp2->n[1], tmp2->n[0]);
			else if (abs(tmp2->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp2->n[2], tmp2->n[1], tmp2->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp2->n[3], tmp2->n[2], tmp2->n[1], tmp2->n[0]);
			sprintf (buf, "%s interim residue %s at bit %ld\n",str, res64, bit);
			OutputBoth (buf);
		}

/* Write a save file every "interimFiles" iterations. */

		if (interimFiles && bit % interimFiles == 0) {
			char	interimfile[20];
			sprintf (interimfile, "%.8s.%03lu",
				 filename, bit / interimFiles);
			if (! writeToFileB (interimfile, bit, a, nrestarts, bpf, x, y)) {
				sprintf (buf, WRITEFILEERR, interimfile);
				OutputBoth (buf);
			}
		}
	}

	clearline (100);

	if (incr == +1) {
		will_try_larger_fft = TRUE;			// All following errors are considered unrecoverable...
		bit = 1;
		explen = bitlen (gb);
/*		bbits = base;
		ulone = 1;
		while (!(bbits & mask)) {
			bbits <<= 1;
			ulone <<= 1;
		}
		bbits <<= 1;
		ulone <<= 1;*/
		gwypcopy (x, y);
		gwypsetnormroutine (0, 1, 0);
//		while (ulone) {
		while (bit < explen) {
                    if (cufftonly)
                        gwypsquare (x);
                    else
                        cuda_gwypsquare (x, 3);
//                    CHECK_IF_ANY_ERROR (x, (explen), explen, 6);
//                    if (bbits & mask) {
                    CHECK_IF_ANY_ERROR (x, (bit), explen, 6);
                    if (bitval (gb, explen-bit-1)) {
                        if (cufftonly)
                            gwypmul (y, x);
                        else
                            cuda_gwypmul (y, x, 3);
//                        CHECK_IF_ANY_ERROR (x, (explen), explen, 6);
                        CHECK_IF_ANY_ERROR (x, (bit), explen, 6);
                    } 
//                    bbits <<= 1;
//                    ulone <<= 1;
                    bit++;
		}
	}

/* See if we've found a probable prime.  If not, format a 64-bit residue. */
/* Old versions of PRP used a non-standard 64-bit residue, computing */
/* 3^N-3 mod N rather than the more standard 3^(N-1) mod N.  Since */
/* some projects recorded these non-standard residues, output that */
/* residue too.  Note that some really old versions printed out the */
/* 32-bit chunks of the non-standard residue in reverse order. */

	ReplaceableLine (2);	/* Replace line */
	gwyptogiant (x, tmp);
	if (!isone (tmp)) {
		*res = FALSE;	/* Not a prime */
		if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
		else if (abs(tmp->sign) < 3)
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
		else if (abs(tmp->sign) < 4)
			sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
		else
			sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
		smulg ((unsigned short)a, tmp); modg (N, tmp); iaddg (-a, tmp);
		if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
			sprintf (oldres64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
		else if (abs(tmp->sign) < 3)
			sprintf (oldres64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
		else if (abs(tmp->sign) < 4)
			sprintf (oldres64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
		else
			sprintf (oldres64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
		if (IniGetInt (INI_FILE, (char*)"OldRes64", 1))
			sprintf (buf, "%s is not prime.  RES64: %s.  OLD64: %s", str, res64, oldres64);
		else
			sprintf (buf, "%s is not prime.  RES64: %s", str, res64);
	}

	if (*res) {
		gwypend_timer (1);
		_unlink (filename);
		if (!factorized) {
			gwypfree (x);
			gwypfree (y);
//			gwdone (gwdata);
                        gwypdone();
			sprintf (buf, "%s is a Probable Prime (Base incompletely factorized).", str);
		}
		if (incr == -1) {				// Morrison test ; start the Lucas sequence
			gwypfree (x);
			gwypfree (y);
			sprintf (buf, "%s may be prime. Starting Lucas sequence...\n", str);
			IniWriteInt(INI_FILE, (char*)"PRPdone", 1);
DoLucas:
			do {
				retval = gLucasequence (N, M, P, gb, jmin, jmax, str, buf, res);
				if (retval == -2) {		// Restart required using next base
					nrestarts++;
					if (nrestarts > maxrestarts) {
						sprintf (buf, "Giving up after %lu restarts...", nrestarts);
						frestart = FALSE;
						*res = FALSE;		// Not proven prime...
						retval = TRUE;
						break;
					}
					IniWriteInt (INI_FILE, (char*)"NRestarts", nrestarts);
					P = genLucasBaseP (N, P+1);
					IniWriteInt (INI_FILE, (char*)"LucasBaseP", P);
//					D = P*P-4;
				}
				if (retval < 0)	{
                                // Restart required for any reason
					sprintf (buf, "Restarting Lucas sequence with P = %lu\n", P);
					gwypdone ();
                                // Setup again the gwypnum code.
					gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
					gwypsetmaxmulbyconst (max(a, P));
					if (!setupok (gwypsetup (dk, 1, n, -1, N), N, str, res)) {
						gwypfree(gk);
						gwypfree(N);
						gwypfree (M);
						gwypfree (tmp);
						gwypfree (tmp2);
						gwypfree (tmp3);
//						*res = FALSE;		// Not proven prime...
						return TRUE; 
					}
				}
			}	while (retval < 0);
			if (retval == FALSE) {
				gwypfree (N);
				gwypfree (gk);
				gwypfree (M);
				gwypfree (tmp);
				gwypfree (tmp2);
				gwypfree (tmp3);
				gwypdone ();
				*res = FALSE;
				return FALSE;
			}
		}
		else {	// Pocklington test ; compute the gcd's
                    sprintf (buf, "Computing GCD'S...");
                    title (buf);
                    sprintf (buf, "%s may be prime, trying to compute gcd's\n", str);
                    OutputStr (buf);
                    for (j=jmax; j>=jmin; j--) {
                        if (bpf[j] == 0)
                            // base prime factor already tested
                            continue;
                        gwypcopy (y, x);    // Computing a^((N-1)/q)
/*                        bbits = bpc[j];
                        ulone = 1;
                        while (!(bbits & mask)) {
                            bbits <<= 1;
                            ulone <<= 1;
                        }
                        bbits <<= 1;
                        ulone <<= 1;
                        while (ulone) {*/
                        bit = 1;
                        explen = bitlen (gbpc[j]);
                        while (bit < explen) {
                           if (cufftonly)
                                gwypsquare (x);
                            else if (zp || generic)
                                cuda_gwypsquare (x, 3);
                            else if(bit==1 || it==0)  
                                {cuda_gwypsquare (x,1);it=1;}
                            else  if(bit != (lasterr_point-1)) 
                                cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (x,2);
//                          CHECK_IF_ANY_ERROR (x, (explen), explen, 6);
//                          if (bbits & mask) {
                            CHECK_IF_ANY_ERROR (x, (bit), explen, 6);
                            if (bitval (gbpc[j], explen-bit-1)) {
                                if (cufftonly)
                                    gwypmul (y, x);
                                else
                                    cuda_gwypmul (y, x, 3);
                                CHECK_IF_ANY_ERROR (x, (explen), explen, 6);
                            }
//                          bbits <<= 1;
//                          ulone <<= 1;
                            bit++;
                        }
                        gwyptogiant (x, tmp);
                        if (isone (tmp)) {
                            if (frestart)
                                continue;
                            if (a==2)
                                // Choose prime bases to have less restarts...
                                newa = 3;
                            else {
                                if (!(a&1))
                                    newa = a + 1;
                                else
                                    newa = a + 2;
                                while (!isPrime(newa))
                                    newa += 2;
                            }
                            nrestarts++;
                            if (nrestarts > maxrestarts) {
/* sprintf (buf, "%s may be prime, but N divides %ld^((N-1)/%lu))-1, giving up after %lu restarts...", str, a, bpf[j], maxrestarts);*/
                                if (bpf[j] == 1) {
                                    gtoc(gbpf[j], bpfstring, strlen(sgb));
                                    sprintf (buf, "%s may be prime, but N divides %lu^((N-1)/%s))-1, giving up after %lu restarts...", str, a, bpfstring, maxrestarts);
                                }
                                else
                                    sprintf (buf, "%s may be prime, but N divides %lu^((N-1)/%lu))-1, giving up after %lu restarts...", str, a, bpf[j], maxrestarts);
                                frestart = FALSE;
                                *res = FALSE;   // Not proven prime...
                            }
                            else {
/*                                sprintf (buf, "%s may be prime, but N divides %ld^((N-1)/%lu))-1, restarting with a=%lu", str, a, bpf[j], newa);*/
                                if (bpf[j] == 1) {
                                    gtoc(gbpf[j], bpfstring, strlen(sgb));
                                    sprintf (buf, "%s may be prime, but N divides %lu^((N-1)/%s))-1, restarting with a=%lu", str, a, bpfstring, newa);
                                }
                                else
                                    sprintf (buf, "%s may be prime, but N divides %lu^((N-1)/%lu))-1, restarting with a=%lu", str, a, bpf[j], newa);
                                a = newa;
                                IniWriteInt (INI_FILE, (char*)"NRestarts", nrestarts);
						IniWriteInt (INI_FILE, (char*)"FermatBase", a);
						frestart = TRUE;
                            }
                        }
                        else {
                            iaddg (-1, tmp);
                            gcdg (N, tmp);
                            if (isone (tmp)) {
//                                sprintf (buf, "%ld^((N-1)/%lu)-1 is coprime to N!\n", a, bpf[j]);
                                if (bpf[j] == 1) {
                                    gtoc(gbpf[j], bpfstring, strlen(sgb));
                                    sprintf (buf, "%lu^((N-1)/%s)-1 is coprime to N!\n", a, bpfstring);
                                }
                                else
                                    sprintf (buf, "%lu^((N-1)/%lu)-1 is coprime to N!\n", a, bpf[j]);
                                OutputStr (buf);
                                if (verbose)
                                    writeResults (buf);
                                bpf[j] = 1;
                                // success for this prime factor of the base, continue
                            }
                            else {
                                *res = FALSE;   /* Not a prime */
                                sprintf (buf, "%s is not prime, although %ld Fermat PSP!!.", str, a);
                                break;  // No need to continue...
                            }
                        }
                    }
                    will_try_larger_fft = FALSE;	
                    if (*res && !frestart)
                        sprintf (buf, "%s is prime! (%d decimal digits)", str, nbdg);
		}
	}
	if (!frestart) {
		gwypfree (N);
		gwypfree (gk);
	}
	gwypfree (M);
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (tmp3);
	if (incr == +1) {
		gwypfree (x);
		gwypfree (y);
	}

#if defined(WIN32) && !defined(_CONSOLE)

	sprintf (buf+strlen(buf), "  Time : "); 
//	ReplaceableLine (2);	/* Replace line */ 

#else

	clearline(100);

#ifdef _CONSOLE
	OutputBoth(buf);
#else
	if (*res) {
		OutputStr((char*)"\033[7m");
		OutputBoth(buf);
		OutputStr((char*)"\033[0m");
	}
	else
		OutputBoth(buf);
#endif

	sprintf (buf, "  Time : "); 

#endif

/* Output the final timings */

	gwypend_timer (1);
//	sprintf (buf+strlen(buf)-1, "  Time: ");
	gwypwrite_timer (buf+strlen(buf), 1, TIMER_CLR | TIMER_NL); 
	if (!frestart) {
		OutputBoth (buf);
		IniWriteString (INI_FILE, (char*)"NRestarts", NULL);
		if (incr == 1)
			IniWriteString (INI_FILE, (char*)"FermatBase", NULL);
		else
			IniWriteString (INI_FILE, (char*)"LucasBaseP", NULL);
	}
	else {
		OutputStr (buf);
		if (verbose && (incr == 1))
			writeResults (buf);
	}

/* Cleanup and return */

	gwypdone ();
	_unlink (filename);
	lasterr_point = 0;
	if (frestart)
		goto restart;
//	gwypdone ();
	if (IniGetInt(INI_FILE, (char*)"PRPdone", 0))
		IniWriteString(INI_FILE, (char*)"PRPdone", NULL);
	IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
	return (TRUE);

/* An error occured, sleep, then try restarting at last save point. */

error:
	gwypfree (M);
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (tmp3);
	gwypfree (x);
	gwypfree (y);
//	gwypdone ();
	*res = FALSE;

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
		sprintf (buf, ERRMSG5, checknumber, str);
		OutputBoth (buf);
		gwypfree (N);
		gwypfree (gk);
		gwypdone ();
		_unlink (filename);
		if (IniGetInt(INI_FILE, (char*)"PRPdone", 0))
			IniWriteString(INI_FILE, (char*)"PRPdone", NULL);
		will_try_larger_fft = FALSE;
		IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
		return (TRUE);
	}

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) { 
		gwypdone ();
		will_try_larger_fft = FALSE;
		return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
		OutputBoth (ERRMSG8);
		IniWriteInt(INI_FILE, (char*)"FFT_Increment", IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
		_unlink (filename);
		will_try_larger_fft = FALSE;
	}
	goto restart;

}

/*
 	Primality testing of k*2^n-1 numbers with the Lucas Lehmer Riesel
	algorithm, using the gwypnums for fast multiplications and squarings.
	Second attempt for a full IBDWT version, using Colin Percival method improved
	by George Woltman.
	Jean Penn?May 2004.
*/

int isLLRP ( 
	unsigned long format, 
	char *sgk,
	unsigned long b_else,	// Lei
	unsigned long n, 
	unsigned long binput,		// Lei
	unsigned long ninput,		// Lei
	unsigned long shift,
	int	*res) 
{ 
	unsigned long iters, index; 
	unsigned long gksize, j, k; 
	unsigned long mask, last, bit, bits; 
	long	vindex, retval;
	gwypnum	x, y; 
	giant	tmp; 
	char	filename[20], buf[sgkbufsize+256], str[sgkbufsize+256],
			sgk1[sgkbufsize], fft_desc[256]; 
	long	write_time = DISK_WRITE_TIME * 60; 
	int		echk, saving, stopping, v1; 
	time_t	start_time, current_time; 
	double	reallyminerr = 1.0; 
	double	reallymaxerr = 0.0; 
	double	dk;

// Lei
	double ddk;
	unsigned long idk = 0;
	giant gk1;
// Lei end
	if (!(format == ABCC || format == ABCK)) {
		gksize = strlen(sgk);				// J. P. Initial gksize

// Lei
		if (b_else != 1) {					// Compute the length of b_else^ninput
			ddk = (double) b_else;
			ddk = ninput * log10 (ddk);
			idk = (long) ddk + 1;
			gksize += idk;					// J. P. Add it to gksize
		}
// Lei end
		else
			idk = 0;
// Lei end
		if ((format == ABCDN) || (format == ABCDNG)) {	// Compute gk = gb^(n-m)-1
			gksize = ndiff*(unsigned long)ceil(log ((double)binput)/log (2.0))+idk;// initial gksize
			gk = newgiant ((gksize >> 2) + 8);		// Allocate space for gk
			itog (binput, gk);
			power (gk, ndiff);
			iaddg (-1, gk);
			sprintf (str, "%lu^%lu-%lu^%lu-1", binput, ninput+ndiff, binput, ninput);
		}
		else {
			gksize = 8*strlen(sgk) + idk;		// J. P. Initial gksize
			gk = newgiant ((gksize >> 2) + 8);	// Allocate space for gk
			ctog (sgk, gk);						// Convert k string to giant
		}

		klen = bitlen(gk);			// Bit length ok initial k multiplier

		if (klen > 53) {		// we must use generic reduction
			dk = 0.0;
		}
		else {					// we can use DWT ; compute the multiplier as a double
			dk = (double)gk->n[0];
			if (gk->sign > 1)
				dk += 65536.0*(double)gk->n[1];
			if (gk->sign > 2)
				dk += 65536.0*65536.0*(double)gk->n[2];
			if (gk->sign > 3)
				dk += 65536.0*65536.0*65536.0*(double)gk->n[3];
		}
// Lei
		if (b_else != 1) {					// Compute the big multiplier
			gk1 = newgiant ((gksize>>2) + 8);
			itog (b_else, gk1);		
			power (gk1, ninput);
			mulg (gk1, gk);
			gwypfree (gk1);
// J.P. shadow   gtoc (gk, sgk1, sgkbufsize);    // Updated k string
		}
// Lei end
		if (shift > 0) {
			gshiftleft (shift, gk);			// Shift k multiplier if requested
			if (b_else != 1)
				strcpy (sgk1, sgk);			// Lei, J.P.
			else
				gtoc (gk, sgk1, sgkbufsize);// Updated k string
		}
		else {
			strcpy (sgk1, sgk);
//	J.P. shadow		if (b_else == 1) strcpy (sgk1, sgk);	// Lei
		}
		if ((format != ABCDN) && (format != ABCDNG))
                    if (b_else != 1)	// Lei, J.P.
			sprintf (str, "%s*%lu^%lu%c1", sgk, binput, ninput, '-');// Number N to test, as a string
                    else
			sprintf (str, "%s*2^%lu%c1", sgk1, n, '-');	// Number N to test, as a string

//	gk must be odd for the LLR test, so, adjust gk and n if necessary.

		while (!bitval(gk, 0)) {
			gshiftright (1, gk);	// update k as a giant
			n++;
		}
	}
	else {
		gk = newgiant ((n>>4)+8);
		itog (1, gk);						// Compute k multiplier
		gshiftleft (n-2, gk);				// Warning : here, n is exponent+1 !
		if (format == ABCK) {
			iaddg (1, gk);
			sprintf (str, "%s*2^%lu%c1 = (2^%lu+1)^2 - 2", sgk, n, '-', n-1);
		}
		else {
			iaddg (-1, gk);
			sprintf (str, "%s*2^%lu%c1 = (2^%lu-1)^2 - 2", sgk, n, '-', n-1);
		}
	}

	klen = bitlen(gk);					// Bit length ok k multiplier
	bits = n + klen;					// Bit length of N
	N =  newgiant ((bits >> 2) + 8);	// Allocate memory for N

//	Compute the number we are testing.

	gtog (gk, N);
	gshiftleft (n, N);
//	mulg (gk, N); 
	iaddg (-1, N);

	Nlen = bitlen (N); 
	nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.

	if ((klen > n) && (nbdg > 400)) {
		if ((format == ABCDN) || (format == ABCDNG))
			sprintf(buf, "2^%lu-1 > 2^%lu, so we can only do a PRP test for %s.\n", ndiff, n, str);
		else
			sprintf(buf, "%s > 2^%lu, so we can only do a PRP test for %s.\n", sgk, n, str);
	    OutputBoth(buf);
// Lei
// Lei shadow   retval = isPRPinternal (str, dk, 2, n, -1, res);
		if ((format == ABCDN) || (format == ABCDNG))
			sprintf (str, "%lu^%lu-%lu^%lu-1", binput, ninput+ndiff, binput, ninput);
		else
			sprintf (str, "%s*%lu^%lu%c1", sgk, binput, ninput, '-');     // Number N to test, as a string
		retval = isPRPinternal (str, dk, binput, ninput, -1, res);
// Lei end

		gwypfree(gk);
		gwypfree(N);
		return retval;
	}
// Lei
//	J.P. shadow sprintf (buf, "Should try prp?\n");
//	J.P. shadow OutputStr (buf);
// Lei end

	if (!IniGetInt(INI_FILE, (char*)"Verify", 0) && !IniGetInt(INI_FILE, (char*)"PRPdone", 0) && (Nlen/klen < 10.0) && (nbdg > 400)) {
								// We have better to do ad first a PRP test.
// Lei
// Lei shadow   retval = isPRPinternal (str, dk, 2, n, -1, res);
		strcpy (buf, str);
		if ((format == ABCDN) || (format == ABCDNG))
			sprintf (str, "%lu^%lu-%lu^%lu-1", binput, ninput+ndiff, binput, ninput);
		else
			sprintf (str, "%s*%lu^%lu%c1", sgk, binput, ninput, '-');     // Number N to test, as a string
				Fermat_only = TRUE;
                retval = isPRPinternal (str, dk, binput, ninput, -1, res);
				Fermat_only = FALSE;
// Lei end

		if (!*res) {
			gwypfree(gk);
			gwypfree(N);
			return retval;
		}
		IniWriteInt(INI_FILE, (char*)"PRPdone", 1);
		strcpy (str, buf);	// Lei
	}

// Lei
//	J.P. shadow sprintf (buf, "Can I get here?\n");
//	J.P. shadow OutputStr (buf);
// Lei end


	k = gk->n[0];
	if(abs(gk->sign) == 2) {	// k is a "small" integer
		k += 65536*gk->n[1];
	}
	else if (abs(gk->sign) > 2)
		k = 0;					// to indicate that k is a big integer.

//restart: 

restart: 

	vindex = 1;					// First attempt

//	p = Nlen; 

	*res = TRUE;		/* Assume it is prime */ 

	tempFileName (filename, 'z', N); //cuda
	if (fileExists (filename) && readFFTLENFromFile (filename, &j, x, NULL));

	gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
	if (!setupok (gwypsetup (dk, binput, ninput, -1, N), N, str, res)) { 
		gwypfree(gk);
		gwypfree(N);
//		*res = FALSE;
		return TRUE;
	}

	x = gwypalloc (); 
	y = gwypalloc ();
 
	tmp =  newgiant (FFTLEN*sizeof(double)/sizeof(short) + 16); 

	last = n-1;

 	gwypsetnormroutine (0, ERRCHK, 0); 

/* Init filename */ 

	//cuda tempFileName (filename, 'z', N); 
 
/* Init the title */ 
 
	title ((char*)"L.L.R. prime test in progress...");
 
/* Optionally resume from save file and output a message */ 
/* indicating we are resuming a test */ 
 
	//cuda if ((g_roundoff==0) && fileExists (filename) && readFromFile (filename, &j, x, NULL)) { 
	if (fileExists (filename) && readFromFile (filename, &j, x, NULL)) { 
		char	fmt_mask[80]; 
		double	pct; 
		pct = trunc_percent (j * 100.0 / n); 
		sprintf (fmt_mask, 
			"Resuming LLR test of %%s at iteration %%ld [%%.%df%%%%]\n", 
			PRECISION); 
		sprintf (buf, fmt_mask, str, j, pct); 
		OutputStr (buf); 
		gwypstart_timer (0); 
		gwypstart_timer (1); 
		time (&start_time); 
	} 
 

/* Otherwise, output a message indicating we are starting test, */ 
/* or resuming the computing of U0. */
 
	else { 
	    if (k==1) {
			if (!isPrime (n)) {
				sprintf (buf, "The Mersenne number %s is not prime because %lu is not prime.\n", str, n); 
				OutputBoth (buf); 
				gwypfree (tmp);
				gwypfree(gk);
				gwypfree(N);
				gwypfree (x); 
				gwypfree (y);
				gwypdone();
				*res = FALSE;
				gwypend_timer (1); 
				return(TRUE);
			}
/*			sprintf (buf, 
				"Prime95 or Mprime are much better to test this Mersenne number !!\n");
			if (verbose)
				OutputBoth(buf);
			else
				OutputStr(buf); */
			v1 = 4;
//			dbltogw ((double) v1, x);
			itogwyp (v1, x);
			gwypclear_timers ();		// Init. timers
			gwypstart_timer (0); 
			gwypstart_timer (1); 
			time (&start_time); 
			goto MERSENNE;
	    }

	    filename[0] = 'u';
	    if ((v1 = gen_v1(gk, n, 0, vindex)) < 0) {
			if (v1 == -1)
				sprintf (buf, "Cannot compute V1 to test %s...\nThis is surprising, please, let me know that!!\nMy E-mail is jpenne@free.fr\n", str);
			else
				sprintf (buf, "%s has a small factor : %d !!\n", str, abs(v1));
			OutputBoth (buf); 
			gwypfree (tmp);
			gwypfree(gk);
			gwypfree(N);
			gwypfree (x); 
			gwypfree (y);
			gwypdone();
			*res = FALSE;
			gwypend_timer (1); 
			return(TRUE);
	    }

	    if (fileExists (filename) && readFromFile (filename, &j, x, y)) { 
			char	fmt_mask[80]; 
			double	pct; 
			pct = trunc_percent (100.0 - j * 100.0 / klen); 
			sprintf (fmt_mask, 
			 "Resuming test of %%s (computing U0) at iteration %%ld [%%.%df%%%%]", 
			 PRECISION); 
			sprintf (buf, fmt_mask, str, klen - j, pct); 
			OutputStr (buf); 
			LineFeed();
			ReplaceableLine (1);	/* Remember where replacable line is */ 
//			will_try_larger_fft = TRUE;
			gwypstart_timer (0); 
			gwypstart_timer (1); 
			time (&start_time); 
	    } 
	    else {
			gwypclear_timers ();		// Init. timers
			gwypstart_timer (0); 
			gwypstart_timer (1); 
			time (&start_time); 
			if (setuponly) {
				if (FFTLEN != OLDFFTLEN) {
					OutputBoth (str); 
					OutputBoth ((char*)" : "); 
				}
			}
			else {
				if (showdigits)
					sprintf (buf, "Starting Lucas Lehmer Riesel prime test of %s (%d decimal digits)\n", str, nbdg);
				else
					sprintf (buf, "Starting Lucas Lehmer Riesel prime test of %s\n", str);
				if (verbose)
					OutputBoth(buf);
				else
					OutputStr (buf); 
			}
			gwypfft_description (fft_desc);
			sprintf (buf, "%s\n", fft_desc);
			if (setuponly) {
				if (FFTLEN != OLDFFTLEN) {
					OutputBoth(buf);
					OLDFFTLEN = FFTLEN;
				}
			}
			else if (verbose)
				OutputBoth(buf);
			else {
				OutputStr(buf);
			}
			if (setuponly) {
				stopping = stopCheck (); 
				gwypfree (tmp);
				gwypfree(gk);
				gwypfree(N);
				gwypfree (x); 
				gwypfree (y);
				gwypdone();
				*res = FALSE;
				gwypend_timer (1); 
				return(!stopping);
			}
			sprintf (buf, "V1 = %d ; Computing U0...", v1);
			OutputStr (buf); 
			LineFeed();
			ReplaceableLine (1);	/* Remember where replacable line is */ 
			itogwyp (v1, x);
			gwypcopy (x, y);
                        it = 0; //cuda
			will_try_larger_fft = FALSE;
			gwypsetnormroutine (0, 1, 0);
			gwypsetaddin (-2);
			if ((1 != lasterr_point) || !maxerr_recovery_mode[0]) {
                            if (cufftonly)
				gwypsquare (y);
                            else
                                cuda_gwypsquare (y,3);
			}
			else {
				gwypsquare_carefully (y);
//				will_try_larger_fft = TRUE;
				if (1 == lasterr_point)
					maxerr_recovery_mode[0] = FALSE;
			}
			CHECK_IF_ANY_ERROR(y, 1, klen, 0)
			if (will_try_larger_fft && (1 == lasterr_point))
				saving = 1;	// Be sure to restart after this recovery iteration!
			will_try_larger_fft = FALSE;
			j = klen - 2;
 	    }
					/* Computing u0 (cf Hans Riesel) */
	    iters = 0; 
	    while (j>0) {

/* Process this iteration */ 

			mask = 1<<j;

			if (k)
				bit = (k&mask);
			else
				bit = bitval (gk, j);
 
			index = klen-j--;
			iters++;

/* Error check the first 50 iterations, before writing an */ 
/* intermediate file (either user-requested stop or a */ 
/* 30 minute interval expired), and every 128th iteration. */ 
 
			stopping = stopCheck (); 
			echk = stopping || ERRCHK || (index <= 50); 
			if (((index & 127) == 0) || (index == 2) || (index == (lasterr_point-1))) {
				echk = 1;
				time (&current_time);
				saving = ((current_time - start_time > write_time) || (index == 2) || (index == (lasterr_point-1)));
			} else
				saving = 0;

			gwypsetnormroutine (0, echk, 0);

			if (bit) {
				gwypsetaddin (-v1);
				if ((index != lasterr_point) || (!maxerr_recovery_mode[1] && !maxerr_recovery_mode[2])) {
                                    if (cufftonly)
					gwypmul (y, x);
                                    else
					cuda_gwypmul (y, x, 3);
					will_try_larger_fft = FALSE;
				}
				else {
					gwypmul_carefully (y, x);
//					will_try_larger_fft = TRUE;
					if (index == lasterr_point)
						maxerr_recovery_mode[1] = FALSE;
				}
				CHECK_IF_ANY_ERROR(x, (index), klen, 1)
				gwypsetaddin (-2);
				if ((index != lasterr_point) || !maxerr_recovery_mode[2]) {
                                    if (cufftonly)
					gwypsquare (y);
                                    else
                                        cuda_gwypsquare (y,3);
                                    will_try_larger_fft = FALSE;
				}
				else {
					gwypsquare_carefully (y);
//					will_try_larger_fft = TRUE;
					if (index == lasterr_point)
						maxerr_recovery_mode[2] = FALSE;
				}
				CHECK_IF_ANY_ERROR(y, (index), klen, 2)
			}
			else {
				gwypsetaddin (-v1);
				if ((index != lasterr_point) || (!maxerr_recovery_mode[3] && !maxerr_recovery_mode[4])) {
                                    if (cufftonly)
					gwypmul (x, y);
                                    else
					cuda_gwypmul (x, y, 3);
					will_try_larger_fft = FALSE;
				}
				else {
					gwypmul_carefully (x, y);
//					will_try_larger_fft = TRUE;
					if (index == lasterr_point)
						maxerr_recovery_mode[3] = FALSE;
				}
				CHECK_IF_ANY_ERROR(y, (index), klen, 3)
				gwypsetaddin (-2);
				if ((index != lasterr_point) || !maxerr_recovery_mode[4]) {
                                    if (cufftonly)
					gwypsquare (x);
                                    else
                                        cuda_gwypsquare (x,3);
                                    will_try_larger_fft = FALSE;
				}
				else {
					gwypsquare_carefully (x);
//					will_try_larger_fft = TRUE;
					if (index == lasterr_point)
						maxerr_recovery_mode[4] = FALSE;
				}
				CHECK_IF_ANY_ERROR(x, (index), klen, 4)
			}

			if (will_try_larger_fft && (index == lasterr_point))
				saving = 1;					// Be sure to restart after this recovery iteration!
			will_try_larger_fft = FALSE;

/* Print a message every so often */ 
 
			if (index % ITER_OUTPUT == 0) { 
				char	fmt_mask[80]; 
				double	pct; 
				pct = trunc_percent (100.0 - j * 100.0 / klen); 
				if (strlen (str) < 40) {
					sprintf (fmt_mask, "%%s, %%.%df%%%% of %%ld", PRECISION); 
					sprintf (buf, fmt_mask, str, pct, klen); 
				}
				else {
					sprintf (fmt_mask, "%%.%df%%%% of %%ld", PRECISION); 
					sprintf (buf, fmt_mask, pct, klen); 
				}
				title (buf); 
				ReplaceableLine (2);	/* Replace line */ 
				sprintf (fmt_mask, 
				 "%%s, iteration : %%ld / %%ld [%%.%df%%%%]", 
				 PRECISION); 
				sprintf (buf, fmt_mask, str, index, klen, pct); 
				OutputStr (buf); 
				if (ERRCHK && index > 30) { 
					OutputStr ((char*)".  Round off: "); 
					sprintf (buf, "%10.10f", reallyminerr); 
					OutputStr (buf); 
					sprintf (buf, " to %10.10f", reallymaxerr); 
					OutputStr (buf); 
				} 
				gwypend_timer (0); 
				if (CUMULATIVE_TIMING) { 
					OutputStr ((char*)".  Time thusfar : "); 
				} 
				else { 
					OutputStr ((char*)".  Time per iteration : "); 
					gwypdivide_timer (0, iters); 
					iters = 0; 
				} 
				gwypprint_timer (0, TIMER_NL | TIMER_OPT_CLR); 
				gwypstart_timer (0); 
			} 
 
/* Print a results file message every so often */ 
 
			if (index % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) { 
				sprintf (buf, "Iteration %ld / %ld\n", index, klen); 
				writeResults (buf); 
			} 
 
/* Write results to a file every DISK_WRITE_TIME minutes */ 
/* On error, retry in 10 minutes (it could be a temporary */ 
/* disk-full situation) */ 
 
			if (saving || stopping) { 
				write_time = DISK_WRITE_TIME * 60; 
				if (! writeToFile (filename, j, x, y)) { 
					sprintf (buf, WRITEFILEERR, filename); 
					OutputBoth (buf); 
					if (write_time > 600) write_time = 600; 
				} 
				time (&start_time); 
 
/* If an escape key was hit, write out the results and return */ 
 
				if (stopping) {
					gwypfree (tmp);
					gwypfree(gk);
					gwypfree(N);
					gwypfree (x); 
					gwypfree (y);
					gwypdone();
					return (FALSE); 
				}
			} 
	    }

		gwypsetaddin (-v1);
		if ((klen != lasterr_point) || !maxerr_recovery_mode[5])
                    if (cufftonly)
                        gwypmul (y, x);
                    else
                        cuda_gwypmul (y, x, 3);
		else {
			gwypmul_carefully (y, x);
//			will_try_larger_fft = TRUE;
			if (klen == lasterr_point)
				maxerr_recovery_mode[5] = FALSE;
		}
		CHECK_IF_ANY_ERROR(x, klen, klen, 5)
		will_try_larger_fft = FALSE;

		ReplaceableLine (2);	/* Replace line */ 
		sprintf (buf, "V1 = %d ; Computing U0...done.\n", v1);
		OutputStr(buf);
		if (verbose) {
			sprintf (buf, "V1 = %d ; Computing U0...done.\n", v1);
			writeResults (buf); 
		}
							/* End of x = u0 computing */
		_unlink (filename);	/* Remove the save file */
	    filename[0] = 'z';	/* restore filename which was modified... */

MERSENNE:
	    sprintf (buf, "Starting Lucas-Lehmer loop..."); 
	    OutputStr (buf); 
		LineFeed();
		j = 1;
	} 


/* Do the Lucas Lehmer Riesel Prime test */ 

	will_try_larger_fft = FALSE;
	ReplaceableLine (1);	/* Remember where replacable line is */  
	iters = 0; 
	gwypsetaddin (-2);
	it = 0; //cuda
//last = 1000; //cuda
	while (j<last) { 

/* Error check the first and last 50 iterations, before writing an */ 
/* intermediate file (either user-requested stop or a */ 
/* 30 minute interval expired), and every 128th iteration. */ 
		stopping = stopCheck (); 
		echk = stopping || ERRCHK || (j <= 50) || (j >= last - 50); 
		if (((j & 127) == 0) || (j == 1) || (j == (lasterr_point-1))) {
			echk = 1;
			time (&current_time);
			saving = ((current_time - start_time > write_time) || (j == 1) || (j == (lasterr_point-1)));
		} else
			saving = 0;

/* Process this iteration */ 

		gwypsetnormroutine (0, echk, 0);

                if (/*(j > 30) && (j < last - 30) && */((j != lasterr_point) || !maxerr_recovery_mode[6])) {
                    if (cufftonly)
                        gwypsquare (x);
                    else if (zp || generic)
                        cuda_gwypsquare (x, 3);
                    else if(j==1 || it==0)
                        {cuda_gwypsquare (x,1);it=1;}
                    else if(j != (lasterr_point-1))
                        cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                    else
                        cuda_gwypsquare (x,2);
		}
		else {
			gwypsquare_carefully (x);
//			will_try_larger_fft = TRUE;
			if (j == lasterr_point)
				maxerr_recovery_mode[6] = FALSE;
		}
		CHECK_IF_ANY_ERROR(x, j, last, 6)
		if (will_try_larger_fft && (j == lasterr_point))
			saving = 1;					// Be sure to restart after this recovery iteration!
		will_try_larger_fft = FALSE;
		j++; 
		iters++; 

/* Print a message every so often */ 
 
		if (j % ITER_OUTPUT == 0) { 
			char	fmt_mask[80]; 
			double	pct; 
			pct = trunc_percent (j * 100.0 / n); 
			if (strlen (str) < 40) {
				sprintf (fmt_mask, "%%.%df%%%% of %%s", PRECISION); 
				sprintf (buf, fmt_mask, pct, str); 
			}
			else {
				sprintf (fmt_mask, "%%.%df%%%% of %%ld", PRECISION); 
				sprintf (buf, fmt_mask, pct, n); 
			}
			title (buf); 
			ReplaceableLine (2);	/* Replace line */ 
			sprintf (fmt_mask, 
				 "%%s, iteration : %%ld / %%ld [%%.%df%%%%]", 
				 PRECISION); 
			sprintf (buf, fmt_mask, str, j, n, pct); 
			OutputStr (buf); 
			if (ERRCHK && j > 30) { 
				OutputStr ((char*)".  Round off: "); 
				sprintf (buf, "%10.10f", reallyminerr); 
				OutputStr (buf); 
				sprintf (buf, " to %10.10f", reallymaxerr); 
				OutputStr (buf); 
			} 
			gwypend_timer (0); 
			if (CUMULATIVE_TIMING) { 
				OutputStr ((char*)".  Time thusfar : "); 
			} 
			else { 
				OutputStr ((char*)".  Time per iteration : "); 
				gwypdivide_timer (0, iters); 
				iters = 0; 
			} 
			gwypprint_timer (0, TIMER_NL | TIMER_OPT_CLR); 
			gwypstart_timer (0); 
		} 
 
/* Print a results file message every so often */ 
 
		if (j % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) { 
			sprintf (buf, "Iteration %ld / %ld\n", j, n); 
			writeResults (buf); 
		} 
 
/* Write results to a file every DISK_WRITE_TIME minutes */ 
/* On error, retry in 10 minutes (it could be a temporary */ 
/* disk-full situation) */ 
 
		if (saving || stopping) { 
			write_time = DISK_WRITE_TIME * 60; 

			if (! writeToFile (filename, j, x, NULL)) { 
				sprintf (buf, WRITEFILEERR, filename); 
				OutputBoth (buf); 
				if (write_time > 600) write_time = 600; 
			} 
			time (&start_time); 

 
/* If an escape key was hit, write out the results and return */ 
 
			if (stopping) {
				gwypfree (tmp);
				gwypfree(gk);
				gwypfree(N);
				gwypfree (x); 
				gwypfree (y);
				gwypdone();
				return (FALSE); 
			}
		} 

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next iteration so that we can compare our */
/* residues to programs that start counter at zero or one. */

		if (interimResidues && j % interimResidues < 2) {
				gwyptogiant (x, tmp);	// The modulo reduction is done here
			if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
			else if (abs(tmp->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
			else if (abs(tmp->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
			sprintf (buf, "%s interim residue %s at iteration %ld\n", str, res64, j);
			OutputBoth (buf);
		}

/* Write a save file every "interimFiles" iterations. */

		if (interimFiles && j % interimFiles == 0) {
			char	interimfile[20];
			sprintf (interimfile, "%.8s.%03lu",
				 filename, j / interimFiles);
			if (! writeToFile (interimfile, j, x, NULL)) {
				sprintf (buf, WRITEFILEERR, interimfile);
				OutputBoth (buf);
			}
		}
	} 
	clearline (100);

	gwyptogiant (x, tmp); 
	if (!isZero (tmp)) { 
		*res = FALSE;				/* Not a prime */ 
		if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
		else if (abs(tmp->sign) < 3)
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
		else if (abs(tmp->sign) < 4)
			sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
		else
			sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
	} 

/* Print results and cleanup */ 

	if (*res) 
		sprintf (buf, "%s is prime! (%d decimal digits)", str, nbdg); 
	else
		sprintf (buf, "%s is not prime.  LLR Res64: %s", str, res64); //msft

#if defined(WIN32) && !defined(_CONSOLE)

	sprintf (buf+strlen(buf), "  Time : "); 
	ReplaceableLine (2);	/* Replace line */ 

#else

	clearline(100);

#ifdef _CONSOLE
	OutputBoth(buf);
#else
	if (*res) {
		OutputStr((char*)"\033[7m");
		OutputBoth(buf);
		OutputStr((char*)"\033[0m");
	}
	else
		OutputBoth(buf);
#endif

	sprintf (buf, "  Time : "); 

#endif

	gwypend_timer (1); 
	gwypwrite_timer (buf+strlen(buf), 1, TIMER_CLR | TIMER_NL); 
	OutputBoth (buf); 
	gwypfree (tmp);
	gwypfree(gk);
	gwypfree(N);
	gwypfree (x); 
	gwypfree (y);
	gwypdone (); 
	filename[0] = 'z';
	_unlink (filename); 
	if (IniGetInt(INI_FILE, (char*)"PRPdone", 0))
		IniWriteString(INI_FILE, (char*)"PRPdone", NULL);
	IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
	lasterr_point = 0;
        g_fftlen = 0;
	return (TRUE); 
 
/* An error occured, sleep, then try restarting at last save point. */ 

error:
	gwypfree (tmp);
	gwypfree (x); 
	gwypfree (y); 
//	gwypdone (); 
	*res = FALSE;

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
		sprintf (buf, ERRMSG5, checknumber, str);
		OutputBoth (buf);
		gwypfree(gk);
		gwypfree(N);
		filename[0] = 'u';
		_unlink (filename);
		filename[0] = 'z';
		gwypdone ();
		_unlink (filename); 
		if (IniGetInt(INI_FILE, (char*)"PRPdone", 0))
			IniWriteString(INI_FILE, (char*)"PRPdone", NULL);
		will_try_larger_fft = FALSE;
                g_fftlen = 0;
		return (TRUE);
	}

/* Output a message saying we are restarting */ 
 
	if (sleep5) OutputBoth (ERRMSG2); 
	OutputBoth (ERRMSG3); 
 
/* Sleep five minutes before restarting */ 
 
	if (sleep5 && ! SleepFive ()) {
		gwypdone ();
		will_try_larger_fft = FALSE;
		return (FALSE); 
	}

/* Restart */ 
 
	if (will_try_larger_fft) {
		OutputBoth (ERRMSG8);
		IniWriteInt(INI_FILE, (char*)"FFT_Increment", IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
		_unlink (filename);
		will_try_larger_fft = FALSE;
	}
	goto restart; 

} 


int isLLRW ( 
	unsigned long format, 
	char *sgk,
	unsigned long n,
	unsigned long shift,
	int	*res) 
{ 
	unsigned long bits, gksize; 
	long retval;
	char str[sgkbufsize+256], sgk1[sgkbufsize]; 

	gksize = strlen(sgk);
	gk = newgiant ((gksize>>1) + 8);	// Allocate one byte per decimal digit + spares
	ctog (sgk, gk);						// Convert k string to giant

	if (shift > 0) {
		gshiftleft (shift, gk);			// Shift k multiplier if requested
		gtoc (gk, sgk1, sgkbufsize);	// Updated k string
	}
	else
		strcpy (sgk1, sgk);

	sprintf (str, "%s*2^%lu%c1", sgk1, n, '-');	// Number N to test, as a string

	bits = n + bitlen(gk);				// Bit length of N
	N =  newgiant ((bits>>4) + 8);		// Allocate memory for N

//	Compute the number we are testing.

	itog (1, N);
	gshiftleft (n, N);
	mulg (gk, N); 
	iaddg (-1, N);
	retval = slowIsWieferich (str, res);
	gwypfree (gk);
	gwypfree (N);
	return retval;
} 


int isProthW ( 
	unsigned long format, 
	char *sgk,
	unsigned long n,
	unsigned long shift,
	int	*res) 
{ 
	unsigned long bits, gksize;
	long retval;
	char	str[sgkbufsize+256], sgk1[sgkbufsize]; 

	gksize = strlen(sgk);
	gk = newgiant ((gksize>>1) + 8);	// Allocate one byte per decimal digit + spares
	ctog (sgk, gk);						// Convert k string to giant

	if (shift > 0) {
		gshiftleft (shift, gk);			// Shift k multiplier if requested
		gtoc (gk, sgk1, sgkbufsize);	// Updated k string
	}
	else
		strcpy (sgk1, sgk);

	sprintf (str, "%s*2^%lu%c1", sgk1, n, '+');	// Number N to test, as a string

	bits = n + bitlen(gk);				// Bit length of N
	N =  newgiant ((bits>>4) + 8);		// Allocate memory for N

//	Compute the number we are testing.

	itog (1, N);
	gshiftleft (n, N);
	mulg (gk, N); 
	iaddg (1, N);

	retval =  slowIsWieferich (str, res);

	gwypfree (gk);
	gwypfree (N);
	return retval;
} 


int isProthP ( 
	unsigned long format, 
	char *sgk,
        unsigned long b_else,	// Lei
	unsigned long n,
	unsigned long binput,		// Lei
	unsigned long ninput,		// Lei
	unsigned long shift,
	int	*res) 
{ 
	unsigned long iters, gksize; 
//	unsigned long p; 
	unsigned long bit, bits; 
	long	a, retval;
	gwypnum	x; 
	giant	tmp, tmp2; 
	char	filename[20], buf[sgkbufsize+256], 
		str[sgkbufsize+256], fft_desc[256], sgk1[sgkbufsize]; 
	long	write_time = DISK_WRITE_TIME * 60; 
	int	echk, saving, stopping; 
	time_t	start_time, current_time; 
	double	reallyminerr = 1.0; 
	double	reallymaxerr = 0.0; 
	double dk;

// Lei
	double ddk;
	unsigned long idk = 0;
	giant gk1;
// Lei end

	gksize = strlen(sgk);				// J.P. Initial gksize

// Lei
	if (b_else != 1) {					// Compute the length of b_else^ninput
		ddk = (double) b_else;
		ddk = ninput * log10 (ddk);
	    idk = (long) ddk + 1;
		gksize += idk;					// J.P. Add it to gksize
	}
// Lei end

	else
		idk = 0;
// Lei end
	if ((format == ABCDN) || (format == ABCDNG)) {	// Compute gk = gb^(n-m)-1
		gksize = ndiff*(unsigned long)ceil(log ((double)binput)/log (2.0))+idk;// initial gksize
		gk = newgiant ((gksize >> 2) + 8);	// Allocate space for gk
		itog (binput, gk);
		power (gk, ndiff);
		iaddg (-1, gk);
		sprintf (str, "%lu^%lu-%lu^%lu+1", binput, ninput+ndiff, binput, ninput);
	}
	else {
		gksize = 8*strlen(sgk) + idk;		// J. P. Initial gksize
		gk = newgiant ((gksize >> 2)  + 8);	// Allocate space for gk
		ctog (sgk, gk);						// Convert k string to giant
	}

	klen = bitlen(gk);					// Length of initial k multiplier

	if (klen > 53) {					// we must use generic reduction
		dk = 0.0;
	}
	else {								// we can use DWT, compute k as a double
		dk = (double)gk->n[0];
		if (gk->sign > 1)
			dk += 65536.0*(double)gk->n[1];
		if (gk->sign > 2)
			dk += 65536.0*65536.0*(double)gk->n[2];
		if (gk->sign > 3)
			dk += 65536.0*65536.0*65536.0*(double)gk->n[3];
	}

// Lei
	if (b_else != 1) {					// Compute the big multiplier
		gk1 = newgiant ((gksize>>2) + 8);
		itog (b_else, gk1);
                power (gk1, ninput);
		mulg (gk1, gk);
		gwypfree (gk1);
	}
// Lei end

	if (shift > 0) {
		gshiftleft (shift, gk);			// Shift k multiplier if requested
		if (b_else != 1)
			strcpy (sgk1, sgk);			// Lei, J.P.
		else
			gtoc (gk, sgk1, sgkbufsize);// Updated k string
	}
	else
		strcpy (sgk1, sgk);

	if ((format != ABCDN) && (format != ABCDNG))
		if (b_else != 1)	// Lei, J.P.
			sprintf (str, "%s*%lu^%lu%c1", sgk, binput, ninput, '+');// Number N to test, as a string
		else
			sprintf (str, "%s*2^%lu%c1", sgk1, n, '+');	// Number N to test, as a string


	bits = n + bitlen(gk);				// Bit length of N
	N =  newgiant ((bits>>4) + 8);		// Allocate memory for N

//	Compute the number we are testing.

	itog (1, N);
	gshiftleft (n, N);
	mulg (gk, N); 
	iaddg (1, N);

//	gk must be odd for the Proth test, so, adjust gk and n if necessary.

	while (bitval(gk, 0) == 0) {
	    gshiftright (1, gk);			// update k as a giant
	    n++;							// update the exponent
	}

	Nlen = bitlen (N); 
	klen = bitlen(gk);
	nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.


	if ((klen > n) && (nbdg > 400)) {
		if ((format == ABCDN) || (format == ABCDNG))
		    sprintf(buf, "2^%lu-1 > 2^%lu, so we can only do a PRP test for %s.\n", ndiff, n, str);
		else
		    sprintf(buf, "%s > 2^%lu, so we can only do a PRP test for %s.\n", sgk, n, str);
	    OutputBoth(buf);

// Lei
// Lei shadow   retval = isPRPinternal (str, dk, 2, n, 1, res);
		if ((format == ABCDN) || (format == ABCDNG))
			sprintf (str, "%lu^%lu-%lu^%lu+1", binput, ninput+ndiff, binput, ninput);
		else
			sprintf (str, "%s*%lu^%lu%c1", sgk, binput, ninput, '+');     // Number N to test, as a string
                retval = isPRPinternal (str, dk, binput, ninput, 1, res);
// Lei end

		gwypfree(gk);
		gwypfree(N);
		return (retval);
	}

/* Init the title */ 
 
	title ((char*)"Proth prime test in progress...");
 
/* Compute the base for the Proth algorithm. */
 
if ((a = genProthBase(gk, n)) < 0) {
	if (a == -1)
		sprintf (buf, "Cannot compute a to test %s...\nThis is surprising, please, let me know that!!\nMy E-mail is jpenne@free.fr\n", str);
	else
		sprintf (buf, "%s has a small factor : %ld !!\n", str, abs(a));
	OutputBoth (buf); 
	*res = FALSE;
	gwypfree(gk);
	gwypfree(N);
	return(TRUE);
}

//restart:

	gwypsetmaxmulbyconst (a);

restart:

//	p = Nlen; 

	*res = TRUE;						/* Assume it is a prime */ 
	
	tempFileName (filename, 'z', N);//cuda
	if (fileExists (filename) && readFFTLENFromFile (filename, &bit, x, NULL)) ;//cuda

	gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
	if (!setupok (gwypsetup (dk, binput, ninput, +1, N), N, str, res)) {
		gwypfree(gk);
		gwypfree(N);
//		*res = FALSE;
		return TRUE;
	}

/* Init tmp = (N-1)/2 to compute a^(N-1)/2 mod N */

	tmp = newgiant(FFTLEN*sizeof(double)/sizeof(short) + 16);
	tmp2 = newgiant(FFTLEN*sizeof(double)/sizeof(short) + 16);
	gtog (N, tmp);
	iaddg (-1, tmp);
	gshiftright (1, tmp);
	Nlen = bitlen (tmp);

/* Init filename */

	//cudatempFileName (filename, 'z', N);

/* Get the current time */
/* Allocate memory */

	x = gwypalloc ();

/* Optionally resume from save file and output a message */
/* indicating we are resuming a test */

	if (fileExists (filename) && readFromFile (filename, &bit, x, NULL)) {
		char	fmt_mask[80];
		double	pct;
		pct = trunc_percent (bit * 100.0 / Nlen);
		sprintf (fmt_mask,
			 "Resuming Proth prime test of %%s at bit %%ld [%%.%df%%%%]\n",
			 PRECISION);
		sprintf (buf, fmt_mask, str, bit, pct);
		OutputStr (buf);
		if (verbose)
			writeResults (buf);
	}

/* Otherwise, output a message indicating we are starting test */

	else {
		gwypclear_timers ();	// Make all timers clean...
		if (setuponly) {
			if (FFTLEN != OLDFFTLEN) {
				OutputBoth (str); 
				OutputBoth ((char*)" : "); 
			}
		}
		else {
			if (showdigits)
				sprintf (buf, "Starting Proth prime test of %s (%d decimal digits)\n", str, nbdg);
			else
				sprintf (buf, "Starting Proth prime test of %s\n", str);
			OutputStr (buf);
			if (verbose)
				writeResults (buf);
		}
		bit = 1;
		itogwyp (a, x);
	}

	gwypstart_timer (0);	// Start loop timer.
	gwypstart_timer (1);	// Start global timer
	time (&start_time);	// Start intermediate file saving time

/* Output a message about the FFT length and the Proth base. */

	gwypfft_description (fft_desc);
#ifdef WIN32
	sprintf (buf, "%s, a = %ld\n", fft_desc, a);
#else
	sprintf (buf, "%s, a = %ld", fft_desc, a);
#endif
	if (!setuponly || (FFTLEN != OLDFFTLEN)) {
		OutputStr (buf);
		if (!setuponly)
			LineFeed();
	}
	sprintf (buf, "%s, a = %ld\n", fft_desc, a);
	if (setuponly) {
		stopping = stopCheck (); 
		if (FFTLEN != OLDFFTLEN) {
			writeResults (buf);
			OLDFFTLEN = FFTLEN;
		}
		gwypfree (tmp);
		gwypfree (tmp2);
		gwypfree(gk);
		gwypfree(N);
		gwypfree (x);
		gwypdone ();
		*res = FALSE;
		return (!stopping);
	}
	else if (verbose) {
#if !defined(WIN32) 
		strcat (buf, "\n");
#endif
		writeResults (buf);
	}
	ReplaceableLine (1);	/* Remember where replaceable line is */

/* Do the Proth test */

	gwypsetmulbyconst (a);
	iters = 0;
//Nlen = 1000; //cuda
	while (bit < Nlen) {

/* Error check the first and last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		stopping = stopCheck ();
		echk = stopping || ERRCHK || (bit <= 50) || (bit >= Nlen-50);
		if (((bit & 127) == 0) || (bit == 1) || (bit == (lasterr_point-1))) {
			echk = 1;
			time (&current_time);
			saving = ((current_time - start_time > write_time) || (bit == 1) || (bit == (lasterr_point-1)));
		} else
			saving = 0;

/* Process this bit */

		if (!(interimResidues && ((bit+1) % interimResidues < 2)) && 
			(bit >= 30) && (bit < Nlen-31) && !maxerr_recovery_mode[6]);

		if (bitval (tmp, Nlen-bit-1)) {
			gwypsetnormroutine (0, echk, 1);
		} else {
			gwypsetnormroutine (0, echk, 0);
		}
                if (/*(j > 30) && (j < last - 30) && */((bit != lasterr_point) || !maxerr_recovery_mode[6])) {
                    if (cufftonly)
                        gwypsquare (x);
                    else if (zp || generic)
                        cuda_gwypsquare (x, 3);
                    else if(bit==1 || it==0)  
                        {cuda_gwypsquare (x,1);it=1;}
                    else  if(bit != (lasterr_point-1)) 
                        cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                    else
                        cuda_gwypsquare (x,2);
		}
		else {
			gwypsquare_carefully (x);
			will_try_larger_fft = TRUE;
			if (bit == lasterr_point)
				maxerr_recovery_mode[6] = FALSE;
		}

		CHECK_IF_ANY_ERROR (x, (bit), Nlen, 6);

/* That iteration succeeded, bump counters */

		if (will_try_larger_fft && (bit == lasterr_point))
			saving = 1;					// Be sure to restart after this recovery iteration!
		will_try_larger_fft = FALSE;
		bit++;
		iters++;

/* Print a message every so often */

		if (bit % ITER_OUTPUT == 0) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (bit * 100.0 / Nlen);
			if (strlen (str) < 40) {
				sprintf (fmt_mask, "%%.%df%%%% of %%s", PRECISION);
				sprintf (buf, fmt_mask, pct, str);
			}
			else {
				sprintf (fmt_mask, "%%.%df%%%% of %%ld", PRECISION);
				sprintf (buf, fmt_mask, pct, Nlen);
			}
			title (buf);
			ReplaceableLine (2);	/* Replace line */
			sprintf (fmt_mask,
				 "%%s, bit: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, str, bit, Nlen, pct);
			OutputStr (buf);
			if (ERRCHK && bit > 30) {
				OutputStr ((char*)".  Round off: ");
				sprintf (buf, "%10.10f", reallyminerr);
				OutputStr (buf);
				sprintf (buf, " to %10.10f", reallymaxerr);
				OutputStr (buf);
			}
			gwypend_timer (0);
			if (CUMULATIVE_TIMING) {
				OutputStr ((char*)".  Time thusfar: ");
			} else {
				OutputStr ((char*)".  Time per bit: ");
				gwypdivide_timer (0, iters);
				iters = 0;
			}
			gwypprint_timer (0, TIMER_NL | TIMER_OPT_CLR);
			gwypstart_timer (0);
		}

/* Print a results file message every so often */

		if (bit % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
			sprintf (buf, "Bit %ld / %ld\n", bit, Nlen);
			writeResults (buf);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving || stopping) {
			write_time = DISK_WRITE_TIME * 60;
			if (! writeToFile (filename, bit, x, NULL)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (buf);
				if (write_time > 600) write_time = 600;
			}
			time (&start_time);

/* If an escape key was hit, write out the results and return */

			if (stopping) {
				gwypfree (tmp);
				gwypfree (tmp2);
				gwypfree(gk);
				gwypfree(N);
				gwypfree (x);
				gwypdone ();
				*res = FALSE;		// To avoid credit message !
				return (FALSE);
			}
		}

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next iteration so that we can compare our */
/* residues to programs that start counter at zero or one. */

		if (interimResidues && bit % interimResidues < 2) {
			gwyptogiant (x, tmp2);		// The modulo reduction is done here
//			iaddg (1, tmp2);			// Compute the (unnormalized) residue
			if (abs(tmp2->sign) < 2)	// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp2->n[0]);
			else if (abs(tmp2->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp2->n[1], tmp2->n[0]);
			else if (abs(tmp2->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp2->n[2], tmp2->n[1], tmp2->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp2->n[3], tmp2->n[2], tmp2->n[1], tmp2->n[0]);
			sprintf (buf, "%s interim residue %s at bit %ld\n", str, res64, bit);
			OutputBoth (buf);
		}

/* Write a save file every "interimFiles" iterations. */

		if (interimFiles && bit % interimFiles == 0) {
			char	interimfile[20];
			sprintf (interimfile, "%.8s.%03lu",
				 filename, bit / interimFiles);
			if (! writeToFile (interimfile, bit, x, NULL)) {
				sprintf (buf, WRITEFILEERR, interimfile);
				OutputBoth (buf);
			}
		}
	}

/* See if we've found a Proth prime.  If not, format a 64-bit residue. */

	clearline (100);

	gwyptogiant (x, tmp2);			// The modulo reduction is done here
	iaddg (1, tmp2);				// Compute the (unnormalized) residue
	if (gcompg (N, tmp2)) {
		*res = FALSE;				/* Not a prime */
		if (abs(tmp2->sign) < 2)	// make a 64 bit residue correct !!
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp2->n[0]);
		else if (abs(tmp2->sign) < 3)
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp2->n[1], tmp2->n[0]);
		else if (abs(tmp2->sign) < 4)
			sprintf (res64, "%04X%04X%04X%04X", 0, tmp2->n[2], tmp2->n[1], tmp2->n[0]);
		else
			sprintf (res64, "%04X%04X%04X%04X", tmp2->n[3], tmp2->n[2], tmp2->n[1], tmp2->n[0]);
	}
	
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree(gk);
	gwypfree(N);
	gwypfree (x);


/* Print results.  Do not change the format of this line as Jim Fougeron of */
/* PFGW fame automates his QA scripts by parsing this line. */

	if (*res)
		sprintf (buf, "%s is prime! (%d decimal digits)", str, nbdg); 
	else
		sprintf (buf, "%s is not prime.  Proth RES64: %s", str, res64);

#if defined(WIN32) && !defined(_CONSOLE)

	sprintf (buf+strlen(buf), "  Time : "); 
	ReplaceableLine (2);	/* Replace line */ 

#else

	clearline(100);

#ifdef _CONSOLE
	OutputBoth(buf);
#else
	if (*res) {
		OutputStr((char*)"\033[7m");
		OutputBoth(buf);
		OutputStr((char*)"\033[0m");
	}
	else
		OutputBoth(buf);
#endif

	sprintf (buf, "  Time : "); 

#endif

/* Output the final timings */

	gwypend_timer (1);
	gwypwrite_timer (buf+strlen(buf), 1, TIMER_CLR | TIMER_NL); 
	OutputBoth(buf);

/* Cleanup and return */

	gwypdone ();
	_unlink (filename);
	IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
	lasterr_point = 0;
        g_fftlen = 0;
	return (TRUE);

/* An error occured, sleep, then try restarting at last save point. */

error:
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (x);
//	gwypdone ();
	*res = FALSE;

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
		sprintf (buf, ERRMSG5, checknumber, str);
		OutputBoth (buf);
		gwypfree(gk);
		gwypfree(N);
		gwypdone ();
		_unlink (filename);
		will_try_larger_fft = FALSE;
                g_fftlen = 0;
		return (TRUE);
	}

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
		gwypdone ();
		will_try_larger_fft = FALSE;
		return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
		OutputBoth (ERRMSG8);
		IniWriteInt(INI_FILE, (char*)"FFT_Increment", IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
		_unlink (filename);
		will_try_larger_fft = FALSE;
	}
	goto restart;
} 

#define LOWFACTORLIMIT 10000	// To factor lower exponent candidates is not useful...
int res1, res2;
long	a;

/********************************* GAUSSIAN-MERSENNE PRIME SEARCH CODE *****************************************

The main purpose of this code is the seach for non-real Gaussian-Mersenne primes.

A Gaussian-Mersenne number can be written as : GM(p) = (1+/-i)^p - 1
It is easy to prove that such a number may be prime in the Gauss ring, only if p is prime,
so we are only considering odd prime values for the exponent p.

More generally, a non-real Gaussian integer x+i*y is prime iff its norm x^2+y^2 is a rational prime.

The norm of GM(p) is N(p) = 2^p + sign*2^((p+1)/2) + 1
where sign = +1 if p = 3 or 5 modulo 8, -1 if p = 1 or 7 modulo 8.

This norm can be rewritten as N(p) = [2^((p-1)/2) + sign]*2^((p+1)/2) + 1, which shows that
N(p) is a Proth number, so the Proth algorithm may be used to prove its primality.
A drawback is that the multiplier in the [] increases so rapidly that it should require using generic mode...

The algorithm implemented below avoids this drawback, and has been suggested to me by Harsh Aggarwal.

The starting point is the Aurifeuillian factorization of M(p) = 4^p+1 :

M(p) = 4^p+1 = (2^p + 2^((p+1)/2) + 1)(2^p - 2^((p+1)/2) + 1)

One of these two factors is the norm N(p) of GM(p) while the other, N'(p) is the norm of another
Gaussian integer : GF(p) = (1+/-i)^p + 1
p beeing odd, such a number is always divisible by 2+/-i, so, N'(p) (like M(p)) has always the trivial factor 5.
But, p beeing prime, GQ(p) = GF(p)/(2+/-i) may be a Gaussian prime, iff N'/5 is a rational prime...

Now, the idea is to run the Proth algorithm, but doing the squarings modulo M(p), and then doing the modulo N
reduction only on the final result. Then, the performances for a given p may be approximatively the same as
for a Lucas-Lehmer test with exponent 1.4*p.

Moreover, using an interim result of the main loop, we have all the info to get the PRP test result of N'/5!

Here is the algorithm in pseudo - C :

GMtest() {

Compute M(p), N(p), N'(p), N'/5;
Compute a such Jacobi(a, N) = -1; // (the Proth base)

x = a;

for (i=1; i<=p-1; i++) {		// This main loop implies only squarings modulo 2^(2*p) + 1, which is optimal!
	x = x*x modulo M;
	if (x == (p-1)/2) y = x;
}

// We have now x = a^(2^(n-1)) modulo M and y = a^(2^((n-1)/2)) modulo M.
// To do the Proth test, we need now to compute R = a^(N-1)/2 nodulo N;
// But, (N-1)/2 = 2^(p-1) + (sign)*2^((p-1)/2), which shows us how to complete :

if (sign == 1)
	R = x*y modulo N;
else if (sign == -1)
	R = x*(y^-1 modulo N) modulo N;

	if (R == -1 modulo N)
		printf("N is prime!");
	else
		printf("N is not prime.");

// The PRP test on N'/5 is slightly more complicated :
// We need to test if a^(N'/5) == a modulo N'/5, which implies a^N' = a^5 modulo N'/5 and finally :
// R' = a^(N'-1] = a^4 modulo N'/5 (a and N'/5 beeing co-prime).
// But, N'-1 = 2*[2^(p-1) - (sign)*2^((p-1)/2)], so :

if (sign == 1)
	R' = x*x*y^-1*y^-1; // all computed modulo N'/5
else if (sign == -1)
	R' = x*x*y*y;		// all computed modulo N'/5
		
	if (R' == a^4 modulo N'/5)
		printf("N'/5 is a-PRP!");
	else
		printf("N'/5 is not prime.");

} // End

Finally, to be efficient, this prime search needs eliminating as many candidates as possible by prefactoring.
To do that I adapted the George Woltman's "factor32.asm" code to find the factors of 4^p+1, and then, of
N and/or N'/5.
This feature is included here, and allows to factor up to 2^86, if really needed.
Also, there is an option to do factoring only jobs.

Jean Penn?March 30 2006

****************************************************************************************************************/

int isGMNP ( 
	char *sgk,
	unsigned long n,
	int	*res) 
{ 
	unsigned long iters; 
	unsigned long ubx, uby, atemp, abits = 0; 
	unsigned long bit, bits, explen, expx, expy, loopshift; 
	gwypnum	x, y; 
	giant	tmp, tmp2, tmp3, apow4; 
	char	filename[20], buf[sgkbufsize+256], 
		str[sgkbufsize+256], strp[sgkbufsize+256], fft_desc[256]; 
	long	write_time = DISK_WRITE_TIME * 60; 
	int	echk, saving, stopping, sign; 
	time_t	start_time, current_time; 
	double	reallyminerr = 1.0; 
	double	reallymaxerr = 0.0; 
	double dk;

	if (!isPrime (n) || n == 2) {
		sprintf (buf, "Gaussian-Mersenne prime test not done because %lu is not an odd prime.\n", n); 
		OutputBoth (buf); 
		*res = FALSE;
		return(TRUE);
	}

	sign = (((n&7) == 3) || ((n&7) == 5))? 1 : 0;	// 1 if positive, 0 if negative
	sprintf (str, "2^%lu%c%s+1",  n, (sign) ? '+' : '-', sgk);	// Number N to test, as a string
	sprintf (strp, "(2^%lu%c%s+1)/5",  n, (sign) ? '-' : '+', sgk);	// Number N' to test, as a string

	bits = 2*n;							// Bit length of M = N*N'
	M = newgiant ((bits>>3) + 8);		// Allocate memory for M = N*N'
	N = newgiant ((bits>>4) + 8);		// Allocate memory for N
	NP = newgiant ((bits>>4) + 8);		// Allocate memory for N'
	gk = newgiant ((bits>>4) + 8);		// Allocate memory for gk
	testn =  newgiant ((bits>>2) + 16);	// For factoring
	testnp = newgiant ((bits>>2) + 16);

//	gk is the multiplier when N is written as gk*2^exponent + 1
//	N = 2^n + s*2^((n+1)/2) + 1 = (2^((n-1)/2) + s)*2^((n+1)/2) + 1
//	So, gk = 2^((n-1)/2) + s and exponent = (n+1)/2, where s is +1 or -1
//	It is only used to compute the Proth base.

//	Compute the numbers we are testing or using.

	itog (1, N);
	gshiftleft (n, N);	    // N  = 2^n
	iaddg (1, N);		    // N  = 2^n + 1
	gtog (N, NP);		    // N' = 2^n + 1
	itog (1, M);
	gshiftleft ((n+1)/2, M);    // M  = 2^((n+1)/2)
	gtog (M, gk);
	gshiftright (1, gk);	    // gk  = 2^((n-1)/2)

	if (sign) {		    // n = 3 or 5 modulo 8
		addg (M, N);	    // N  = 2^n + 2^((n+1)/2) + 1
		iaddg (1, gk);	    // gk  = 2^((n-1)/2) + 1
		subg (M, NP);	    // N' = 2^n - 2^((n+1)/2) + 1
	}
	else {			    // n = 1 or 7 modulo 8
		subg (M, N);	    // N  = 2^n - 2^((n+1)/2) + 1
		iaddg (-1, gk);	    // gk  = 2^((n-1)/2) - 1
		addg (M, NP);	    // N' = 2^n + 2^((n+1)/2) + 1
	}

	uldivg (5, NP);		    // NP = N'/5
	itog (1, M);
	gshiftleft (2*n, M);	    // M  = 2^(2*n)
	iaddg (1, M);		    // M  = N*N' = 2^(2*n) + 1

	Nlen = 2*n+1; 
	nbdg1 = gnbdg (N, 10);	    // Size of the two candidates
	nbdg2 = gnbdg (NP, 10);

	res1 = res2 = 1;	    // Assume N and NP are prime...
        if (!setupok ((nbdg1 < 400), N, str, &res1)&&!setupok ((nbdg2 < 400), NP, strp, &res2)) {   // Force APRCL test for small numbers...
		*res = (res1 || res2);
		gwypfree(N);
		gwypfree(NP);
		gwypfree(M);
		gwypfree(gk);
		gwypfree(testn);
		gwypfree(testnp);
		return (TRUE); 
	}

 	dk = 1.0;	        // k == 1 for the modulo N*N'

/* Compute the base for the Proth algorithm. */
 
	if ((a = genProthBase(gk, (n+1)/2)) < 0) {
		if (a == -1)
			sprintf (buf, "Cannot compute a to test %s...\nThis is surprising, please, let me know that!!\nMy E-mail is jpenne@free.fr\n", str);
		else
			sprintf (buf, "%s has a small factor : %ld !!\n", str, abs(a));
		OutputBoth (buf); 
		*res = res1 = res2 = FALSE;
		gwypfree(gk);
		gwypfree(N);
		gwypfree(NP);
		gwypfree(M);
		gwypfree(testn);
		gwypfree(testnp);
		return(TRUE);
	}


//restart:

	gwypsetmaxmulbyconst (a);

restart:

/* Assume intermediate results of the length of N*N'. */ 

	*res = TRUE;						/* Assume it is a prime */ 

	gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));

        if (!setupok (gwypsetup (dk, 2, 2*n, +1, M), M, str, res)) { 	// Setup the DWT mode
//		*res = res1 = res2 = FALSE;
		gwypfree(gk);
		gwypfree(N);
		gwypfree(NP);
		gwypfree(M);
		gwypfree(testn);
		gwypfree(testnp);
		return TRUE;
	}

	expx = n-1;
	expy = expx/2;

/* More initializations... */

	tmp = newgiant((Nlen >> 2) + 8);
	tmp2 = newgiant((Nlen >> 2) + 8);
	tmp3 = newgiant((Nlen >> 2) + 8);
	apow4 = newgiant(32);
	itog (a, apow4);
	smulg ((unsigned short)a, apow4);
	smulg ((unsigned short)a, apow4);
	smulg ((unsigned short)a, apow4);

/* Init filename */

	tempFileName (filename, 'z', N);

/* Get the current time */
/* Allocate memory */

	x = gwypalloc ();
	y = gwypalloc ();

/* Optionally resume from save file and output a message */
/* indicating we are resuming a test */

/* ubx and uby are the current units bit positions, in x and y, respectively.*/

	if (fileExists (filename) && gmreadFromFile (filename, &bit, &ubx, &uby, x, y)) {
		char	fmt_mask[80];
		double	pct;
		pct = trunc_percent (bit * 100.0 / expx);
		sprintf (fmt_mask,
			 "Resuming Proth prime test of %%s at bit %%ld [%%.%df%%%%]\n",
			 PRECISION);
		sprintf (buf, fmt_mask, str, bit, pct);
		OutputStr (buf);
		if (verbose)
			writeResults (buf);
	}

/* Otherwise, output a message indicating we are starting test */

	else {
		gwypclear_timers ();		// Init. timers
		if (setuponly) {
			if (FFTLEN != OLDFFTLEN) {
				OutputBoth (str); 
				OutputBoth ((char*)(char*)" : "); 
			}
		}
		else {
			if (showdigits)
				sprintf (buf, "Starting Proth prime test of %s (%d decimal digits)\n", str, nbdg1);
			else
				sprintf (buf, "Starting Proth prime test of %s\n", str);
			OutputStr (buf);
			if (verbose)
				writeResults (buf);
		}

		bit = 1;

/* Compute a random shift for the initial value */

		srand ((unsigned int) time (NULL));
		ubx = (rand() << 16) + rand();
		atemp = a;
		while (atemp) {						// Compute the bit length of the Proth base a
			atemp >>= 1;
			abits++;
		}
		ubx = ubx % (bits-abits);			// Be sure that the shift is not too large...
		uby = 0;


/* Compute the left shifted initial value */

		itog (a, tmp3);
		gshiftleft (ubx, tmp3);

		gianttogwyp (tmp3, x);
		gianttogwyp (M, y);
	}

	gwypstart_timer (0);
	gwypstart_timer (1);
	time (&start_time);		// Get current time

/* Output a message about the FFT length and the Proth base. */

	gwypfft_description (fft_desc);
#ifdef WIN32
	sprintf (buf, "%s, a = %ld\n", fft_desc, a);
#else
	sprintf (buf, "%s, a = %ld", fft_desc, a);
#endif
	if (!setuponly || (FFTLEN != OLDFFTLEN)) {
		OutputStr (buf);
		if (!setuponly)
			LineFeed();
	}
	sprintf (buf, "%s, a = %ld\n", fft_desc, a);
	if (setuponly) {
		stopping = stopCheck (); 
		if (FFTLEN != OLDFFTLEN) {
			writeResults (buf);
			OLDFFTLEN = FFTLEN;
		}
		gwypfree (tmp);
		gwypfree (tmp2);
		gwypfree (tmp3);
		gwypfree (apow4);
		gwypfree(gk);
		gwypfree(N);
		gwypfree(NP);
		gwypfree(M);
		gwypfree(testn);
		gwypfree(testnp);
		gwypfree (x);
		gwypfree (y);
		gwypdone ();
		*res = res1 = res2 = FALSE;
		return (!stopping);
	}
	else if (verbose) {
#if !defined(WIN32) 
		strcat (buf, "\n");
#endif
		writeResults (buf);
	}
	ReplaceableLine (1);	/* Remember where replaceable line is */

/* Init the title */

	title ((char*)"G.M.N. prime test in progress...");

/* Do the Proth test */

	iters = 0;
	loopshift = (bit >= expy) ? expy : 0;
	explen = (bit >= expy) ? expx : expy;
	while (bit <= expx) {

/* Error check the first and last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		stopping = stopCheck ();
		echk = stopping || ERRCHK || (bit <= (50+loopshift)) || (bit >= explen-50);
		if (((bit & 127) == 0) || (bit == 1) || (bit == (lasterr_point-1))) {
			echk = 1;
			time (&current_time);
			saving = ((current_time - start_time > write_time) || (bit == 1) || (bit == (lasterr_point-1)));
		} else
			saving = 0;

/* Process this bit */

		if (!(interimResidues && ((bit+1) % interimResidues < 2)) && 
			(bit >= (30+loopshift)) && (bit < explen-31) && !maxerr_recovery_mode[6]);


		gwypsetnormroutine (0, echk, 0);
		if (/*(bit > (30+loopshift)) && (bit < explen-30) && */((bit != lasterr_point) || !maxerr_recovery_mode[6])) {
                    if (cufftonly)
                        gwypsquare (x);
                    else if (zp || generic)
                        cuda_gwypsquare (x, 3);
                    else if(bit==1 || it==0)  
                        {cuda_gwypsquare (x,1);it=1;}
                    else  if(bit != (lasterr_point-1)) 
                        cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                    else
                        cuda_gwypsquare (x,2);
                }
		else {
			gwypsquare_carefully (x);
			will_try_larger_fft = TRUE;
			if (bit == lasterr_point)
				maxerr_recovery_mode[6] = FALSE;
		}

		ubx <<= 1;
		if (ubx >= bits) ubx -= bits;		// Compute the doubled shift modulo 2*n

		if (bit == expy) {
			gwypcopy (x, y);
			uby = ubx;
			loopshift = expy;
			explen = expx;
		}

		CHECK_IF_ANY_ERROR (x, (bit), explen, 6);

/* That iteration succeeded, bump counters */

		if (will_try_larger_fft && (bit == lasterr_point))
			saving = 1;					// Be sure to restart after this recovery iteration!
		will_try_larger_fft = FALSE;
		bit++;
		iters++;

/* Print a message every so often */

		if (bit % ITER_OUTPUT == 0) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (bit * 100.0 / expx);
			if (strlen (str) < 40) {
				sprintf (fmt_mask, "%%.%df%%%% of %%s", PRECISION);
				sprintf (buf, fmt_mask, pct, str);
			}
			else {
				sprintf (fmt_mask, "%%.%df%%%% of %%ld", PRECISION);
				sprintf (buf, fmt_mask, pct, explen);
			}
			title (buf);
			ReplaceableLine (2);	/* Replace line */
			sprintf (fmt_mask,
				 "%%s, bit: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, str, bit, expx, pct);
			OutputStr (buf);
			if (ERRCHK && bit > 30) {
				OutputStr ((char*)".  Round off: ");
				sprintf (buf, "%10.10f", reallyminerr);
				OutputStr (buf);
				sprintf (buf, " to %10.10f", reallymaxerr);
				OutputStr (buf);
			}
			gwypend_timer (0);
			if (CUMULATIVE_TIMING) {
				OutputStr ((char*)".  Time thusfar: ");
			} else {
				OutputStr ((char*)".  Time per bit: ");
				gwypdivide_timer (0, iters);
				iters = 0;
			}
			gwypprint_timer (0, TIMER_NL | TIMER_OPT_CLR);
			gwypstart_timer (0);
		}

/* Print a results file message every so often */

		if (bit % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
			sprintf (buf, "Bit %ld / %ld\n", bit, explen);
			writeResults (buf);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving || stopping) {
			write_time = DISK_WRITE_TIME * 60;
			if (! gmwriteToFile (filename, bit, ubx, uby, x, y)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (buf);
				if (write_time > 600) write_time = 600;
			}
			time (&start_time);


/* If an escape key was hit, write out the results and return */

			if (stopping) {
				gwypfree (tmp);
				gwypfree (tmp2);
				gwypfree (tmp3);
				gwypfree (apow4);
				gwypfree(gk);
				gwypfree(N);
				gwypfree(NP);
				gwypfree(M);
				gwypfree(testn);
				gwypfree(testnp);
				gwypfree (x);
				gwypfree (y);
				gwypdone ();
				*res = res1 = res2 = FALSE;
				return (FALSE);
			}
		}

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next iteration so that we can compare our */
/* residues to programs that start counter at zero or one. */

		if (interimResidues && bit >= expy && bit % interimResidues < 2) {

			itog (1, tmp3);					// Restore the value of x from the shifted one.
			gshiftleft (ubx, tmp3);
			invg (M,tmp3);
			gtog (M, testn);
			if (ubx&2)						// View if a sign change on x is necessary.
				subg (tmp3, testn);
			else
				gtog (tmp3, testn);

			gwyptogiant (x, tmp3);	// The modulo reduction is done here
			mulg (tmp3, testn);
			modg (M, testn);

			itog (1, tmp3);			// Restore the value of y from the shifted one.
			gshiftleft (uby, tmp3);
			invg (M,tmp3);
			gtog (M, testnp);
			if (uby&2)				// View if a sign change on y is necessary.
				subg (tmp3, testnp);
			else
				gtog (tmp3, testnp);
			gwyptogiant (y, tmp3);	// The modulo reduction is done here
			mulg (tmp3, testnp);
			modg (M, testnp);
			gtog (testn, tmp);
			gtog (testnp, tmp2);

			if (sign) {
				mulg (tmp2, tmp);
				modg (N, tmp);
				iaddg (1, tmp);		// Compute the (unnormalized) residue
			}
			else {
				invg (N, tmp2);
				mulg (tmp2, tmp);
				modg (N, tmp);
				iaddg (1, tmp);
			}
			if (abs(tmp->sign) < 2)	// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
			else if (abs(tmp->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
			else if (abs(tmp->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
			sprintf (buf, "GM%ld interim residue %s at iteration %ld\n", n, res64, bit);
			OutputBoth (buf);
		}

/* Write a save file every "interimFiles" iterations. */

		if (interimFiles && bit % interimFiles == 0) {
			char	interimfile[20];
			sprintf (interimfile, "%.8s.%03lu",
				 filename, bit / interimFiles);
			if (! gmwriteToFile (interimfile, bit, ubx, uby, x, y)) {
				sprintf (buf, WRITEFILEERR, interimfile);
				OutputBoth (buf);
			}
		}
	}

	clearline (100);

	itog (1, tmp3);			// Restore the value of x from the shifted one.
	gshiftleft (ubx, tmp3);
	invg (M,tmp3);
	gtog (M, testn);
	if (ubx&2)				// View if a sign change on x is necessary.
		subg (tmp3, testn);
	else
		gtog (tmp3, testn);

	gwyptogiant (x, tmp3);	// The modulo reduction is done here
	mulg (tmp3, testn);
	modg (M, testn);

	itog (1, tmp3);			// Restore the value of y from the shifted one.
	gshiftleft (uby, tmp3);
	invg (M,tmp3);
	gtog (M, testnp);
	if (uby&2)				// View if a sign change on y is necessary.
		subg (tmp3, testnp);
	else
		gtog (tmp3, testnp);
	gwyptogiant (y, tmp3);	// The modulo reduction is done here
	mulg (tmp3, testnp);
	modg (M, testnp);
	gtog (testn, tmp);
	gtog (testnp, tmp2);

	if (sign) {
		mulg (tmp2, tmp);
		modg (N, tmp);
		iaddg (1, tmp);		// Compute the (unnormalized) residue
	}
	else {
		invg (N, tmp2);
		mulg (tmp2, tmp);
		modg (N, tmp);
		iaddg (1, tmp);
	}

/* See if we've found a Proth prime.  If not, format a 64-bit residue. */

	if (gcompg (N, tmp) != 0) {
		res1 = FALSE;
		if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
		else if (abs(tmp->sign) < 3)
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
		else if (abs(tmp->sign) < 4)
			sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
		else
			sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
	}


/* Print results.  Do not change the format of this line as Jim Fougeron of */
/* PFGW fame automates his QA scripts by parsing this line. */

	if (res1)
		sprintf (buf, "%s is prime! (%d decimal digits)\n", str, nbdg1);
	else
		sprintf (buf, "%s is not prime.  Proth RES64: %s\n", str, res64);

#if defined(WIN32) && !defined(_CONSOLE)

	ReplaceableLine (2);	/* Replace line */ 
	OutputBoth (buf);

#else

	clearline(100);

#ifdef _CONSOLE
	OutputBoth(buf);
#else
	if (res1) {
		OutputStr((char*)"\033[7m");
		OutputBoth(buf);
		OutputStr((char*)"\033[0m");
	}
	else
		OutputBoth(buf);
#endif

	sprintf (buf, "  Time : "); 

#endif

	gtog (testn, tmp);
	gtog (testnp, tmp2);

	if (sign) {
		mulg (tmp, tmp);
		modg (NP, tmp);
		mulg (tmp2, tmp2);
		mulg (apow4, tmp2);
		modg (NP, tmp2);
	}
	else {
		mulg (tmp, tmp);
		modg (NP, tmp);
		mulg (tmp2, tmp);
		modg (NP, tmp);
		mulg (tmp2, tmp);
		modg (NP, tmp);
		gtog (apow4, tmp2);
		modg (NP, tmp2);
	}

	if (gcompg (tmp2, tmp) != 0) {
		subg (tmp2, tmp);
		res2 = FALSE;				/* Not a prime */
		if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
		else if (abs(tmp->sign) < 3)
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
		else if (abs(tmp->sign) < 4)
			sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
		else
			sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
	}


/* Print results.  Do not change the format of this line as Jim Fougeron of */
/* PFGW fame automates his QA scripts by parsing this line. */

	if (res2)
		sprintf (buf, "%s is %ld-PRP! (%d decimal digits)", strp, a, nbdg2);
	else
		sprintf (buf, "%s is not prime.  RES64: %s", strp, res64);

#ifdef WIN32

	sprintf (buf+strlen(buf), "  Time: ");

#else

	if (res2) {
		OutputStr((char*)"\033[7m");
		OutputBoth (buf);
		OutputStr((char*)"\033[0m");
	}
	else
		OutputBoth (buf);
	sprintf (buf, "  Time: ");

#endif

/* Output the final timings */

	gwypend_timer (1);
	gwypwrite_timer (buf+strlen(buf), 1, TIMER_CLR | TIMER_NL); 
	OutputBoth (buf);

	*res = (res1 || res2);

/* Cleanup and return */

        gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (tmp3);
	gwypfree (apow4);
	gwypfree(gk);
	gwypfree(N);
	gwypfree(NP);
	gwypfree(M);
	gwypfree(testn);
	gwypfree(testnp);
	gwypfree (x);
	gwypfree (y);

	gwypdone ();
	_unlink (filename);
	IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
	lasterr_point = 0;
	return (TRUE);

/* An error occured, sleep, then try restarting at last save point. */

error:
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (tmp3);
	gwypfree (apow4);
	gwypfree (x);
	gwypfree (y);
//	gwypdone ();

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
		sprintf (buf, ERRMSG5, checknumber, str);
		OutputBoth (buf);
		*res = res1 = res2 = FALSE;
		gwypfree(gk);
		gwypfree(N);
		gwypfree(NP);
		gwypfree(M);
		gwypfree(testn);
		gwypfree(testnp);
		gwypdone ();
		_unlink (filename);
		will_try_larger_fft = FALSE;
		return (TRUE);
	}

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
		gwypdone ();
		will_try_larger_fft = FALSE;
		return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
		OutputBoth (ERRMSG8);
		IniWriteInt(INI_FILE, (char*)"FFT_Increment", IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
		_unlink (filename);
		will_try_larger_fft = FALSE;
	}
	goto restart;
} 

/************************************** Strong Fermat PRP test code for Wagstaff numbers *************************************

Wagstaff numbers are numbers of the form W(n) = (2^n+1)/3 where n is an odd integer.

/****************************************************************************************************************************/


int isWSPRP ( 
	char *sgk,
	unsigned long n,
	int	*res) 
{ 
	unsigned long iters; 
	unsigned long ubx, uby, atemp, abits = 0;
	unsigned long bit, bits, expx, dovrbareix = vrbareix; 
	gwypnum	x, y; 
	giant	tmp, tmp2, gx0; 
	char	filename[20], buf[sgkbufsize+256], 
		fft_desc[256], oldres64[17]; 
	long	a, write_time = DISK_WRITE_TIME * 60; 
	int	echk, saving, stopping; 
	time_t	start_time, current_time; 
	double	reallyminerr = 1.0; 
	double	reallymaxerr = 0.0; 

	if (!isPrime (n) || n == 2) {
		sprintf (buf, "(2^%lu+1)/3 SPRP test not done because %lu is not an odd prime.\n", n, n); 
		OutputBoth (buf); 
		*res = FALSE;
		return(TRUE);
	}

	bits = n;							// Bit length of NP
	Nlen = bits + 1;					// for read/write intermediate files
	M = newgiant ((bits>>3) + 8);		// Allocate memory for M
	NP = newgiant ((bits>>3) + 8);		// Allocate memory for NP
	testn =  newgiant ((bits>>3) + 16);	// For factoring

//	Compute the numbers we are testing or using.

	itog (1, M);
	gshiftleft (n, M);	// M  = 2^n
	iaddg (1, M);		// M  = 2^n + 1
	gtog (M, NP);		// NP  = 2^n + 1
	uldivg (3, NP);		// NP  = (2^n + 1)/3
        nbdg = gnbdg (NP, 10);

        if (!setupok ((nbdg < 400), NP, sgk, res)) {
		gwypfree(NP);   // Force APRCL test for small numbers...
		gwypfree(M);
		gwypfree(testn);
		return (TRUE); 
	}

// Test if we are resuming a PRP test.

	tempFileName (filename, 's', NP);
	if (fileExists (filename)) {				// Resuming a Fermat SPRP test
		dovrbareix = FALSE;
		goto restart;
	}

	tempFileName (filename, 'z', NP);
	if (fileExists (filename)) {				// Resuming a Vrba-Reix test
		dovrbareix = TRUE;
	}

restart:

	if (dovrbareix) {						// Compute the seed for the Vrba-Reix test
		gx0 =  newgiant ((bits >> 4) + 8);	// Allocate memory for gx0
		gtog (NP, gx0);						// gx0 = NP
		iaddg (3, gx0);						// gx0 = N+3
		gshiftright (1, gx0);				// gx0 = (N+3)/2 = 3/2 mod N = 3*2^(-1) mod N
		expx = n-1;
		tempFileName (filename, 'z', NP);	// Set the filename to zxxxxxxx
	}
	else {									// Set he base for the SPRP test
		a = IniGetInt (INI_FILE, (char*)"FBase", 3);
		gwypsetmaxmulbyconst (a);
		expx = n;
		tempFileName (filename, 's', NP);	// Set the filename to sxxxxxxx
	}

	*res = TRUE;						/* Assume it is a prime */ 

	gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
	if (!setupok (gwypsetup (1.0, 2, n, +1, NP), NP, sgk, res)) { 	// Setup the DWT mode
//		*res = FALSE;
		gwypfree(NP);
		gwypfree(M);
		gwypfree(testn);
		return (TRUE); 
	}

/* More initializations... */

	tmp = newgiant((bits >> 2) + 8);
	tmp2 = newgiant((bits >> 2) + 8);

/* Get the current time */
/* Allocate memory */

	x = gwypalloc ();
	y = gwypalloc ();

/* Optionally resume from save file and output a message */
/* indicating we are resuming a test */

/* ubx and uby are the current units bit positions, in x and y, respectively.*/

	if (fileExists (filename) && gmreadFromFile (filename, &bit, &ubx, &uby, x, y)) {
		char	fmt_mask[80];
		double	pct;
		pct = trunc_percent (bit * 100.0 / expx);
		if (dovrbareix)
			sprintf (fmt_mask,
				"Resuming Vrba-Reix test of %%s at bit %%ld [%%.%df%%%%]\n",
				 PRECISION);
		else
			sprintf (fmt_mask,
				"Resuming SPRP test of %%s at bit %%ld [%%.%df%%%%]\n",
				 PRECISION);
		sprintf (buf, fmt_mask, sgk, bit, pct);
		OutputStr (buf);
		if (verbose)
			writeResults (buf);
	}

/* Otherwise, output a message indicating we are starting test */

	else {
		gwypclear_timers ();		// Init. timers
		if (showdigits) {
			if (dovrbareix)
				sprintf (buf, "Starting Vrba-Reix test of %s (%d decimal digits)\n", sgk, nbdg);
			else
				sprintf (buf, "Starting SPRP test of %s (%d decimal digits)\n", sgk, nbdg);
		}
		else {
			if (dovrbareix)
				sprintf (buf, "Starting Vrba-Reix test of %s\n", sgk);
			else
				sprintf (buf, "Starting SPRP test of %s\n", sgk);
		}
		OutputStr (buf);
		if (verbose)
			writeResults (buf);
		bit = 1;

/* Compute a random shift for the initial value */

		srand ((unsigned int) time (NULL));
		ubx = (rand() << 16) + rand();
		if (dovrbareix) {
			ubx = ubx % (bits);			// To be sure that the shift is not too large...
		}
		else {
			atemp = a;
			while (atemp) {				// Compute the bit length of the Fermat base a
				atemp >>= 1;
				abits++;
			}
			ubx = ubx % (bits-abits);	// To be sure that the shift is not too large...
		}
		uby = 0;


/* Compute the left shifted initial value */

		if (dovrbareix) {
			gtog (gx0, tmp2);
			gshiftleft (ubx, tmp2);
			modg (M, tmp2);
		}
		else {
			itog (a, tmp2);
			gshiftleft (ubx, tmp2);
		}
		gianttogwyp (tmp2, x);
		gwypcopy (x, y);
	}

	gwypstart_timer (0);
	gwypstart_timer (1);
	time (&start_time);				// Get current time

/* Output a message about the FFT length. */

	gwypfft_description (fft_desc);
#ifdef WIN32
	sprintf (buf, "%s\n", fft_desc);
#else
	sprintf (buf, "%s", fft_desc);
#endif
	OutputStr (buf);
	LineFeed ();
	if (verbose) {
#if !defined(WIN32) 
		strcat (buf, "\n");
#endif
		writeResults (buf);
	}
	ReplaceableLine (1);	/* Remember where replaceable line is */

/* Init the title */

	if (dovrbareix)
		title ((char*)"Wagstaff numbers Vrba-Reix test in progress...");
	else
		title ((char*)"Wagstaff numbers SPRP test in progress...");

/* Do the PRP test */

	iters = 0;
	while (bit <= expx) {

/* Error check the first and last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		stopping = stopCheck ();
		echk = stopping || ERRCHK || (bit <= 50) || (bit >= expx-50);
		if (((bit & 127) == 0) || (bit == 1) || (bit == (lasterr_point-1))) {
			echk = 1;
			time (&current_time);
			saving = ((current_time - start_time > write_time) || (bit == 1) || (bit == (lasterr_point-1)));
		} else
			saving = 0;

/* Process this bit */

		if (!(interimResidues && ((bit+1) % interimResidues < 2)) && 
			(bit >= 30) && (bit < expx-31) && !maxerr_recovery_mode[6]);


		gwypsetnormroutine (0, echk, 0);


		ubx <<= 1;
		if (ubx >= bits) ubx -= bits;				// Compute the doubled shift modulo n

		if (dovrbareix)								// Fix-up the addin constant
			if (ubx&1)								// See if a change of sign is needed
				gwypsetaddinatpowerofb (2, ubx);
			else
				gwypsetaddinatpowerofb (-2, ubx);

		if (/*(bit > 30) && (bit < expx-30) && */((bit != lasterr_point) || !maxerr_recovery_mode[6])) {
                    if (cufftonly)
                        gwypsquare (x);
                    else if (zp || generic)
                        cuda_gwypsquare (x, 3);
                    else if(bit==1 || it==0)  
                        {cuda_gwypsquare (x,1);it=1;}
                    else  if(bit != (lasterr_point-1)) 
                        cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                    else
                        cuda_gwypsquare (x,2);
                }
		else {
			gwypsquare_carefully (x);
			will_try_larger_fft = TRUE;
			if (bit == lasterr_point)
				maxerr_recovery_mode[6] = FALSE;
		}
		if (!dovrbareix && bit == (expx - 1)) {
			gwypcopy (x, y);
			uby = ubx;
		}


		CHECK_IF_ANY_ERROR (x, (bit), expx, 6);

/* That iteration succeeded, bump counters */

		if (will_try_larger_fft && (bit == lasterr_point))
			saving = 1;					// Be sure to restart after this recovery iteration!
		will_try_larger_fft = FALSE;
		bit++;
		iters++;

/* Print a message every so often */

		if (bit % ITER_OUTPUT == 0) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (bit * 100.0 / expx);
			sprintf (fmt_mask, "%%.%df%%%% of %%ld", PRECISION);
			sprintf (buf, fmt_mask, pct, expx);
			title (buf);
			ReplaceableLine (2);	/* Replace line */
			sprintf (fmt_mask,
				 "%%s, bit: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, sgk, bit, expx, pct);
			OutputStr (buf);
			if (ERRCHK && bit > 30) {
				OutputStr ((char*)".  Round off: ");
				sprintf (buf, "%10.10f", reallyminerr);
				OutputStr (buf);
				sprintf (buf, " to %10.10f", reallymaxerr);
				OutputStr (buf);
			}
			gwypend_timer (0);
			if (CUMULATIVE_TIMING) {
				OutputStr ((char*)".  Time thusfar: ");
			} else {
				OutputStr ((char*)".  Time per bit: ");
				gwypdivide_timer (0, iters);
				iters = 0;
			}
			gwypprint_timer (0, TIMER_NL | TIMER_OPT_CLR);
			gwypstart_timer (0);
		}

/* Print a results file message every so often */

		if (bit % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
			sprintf (buf, "Bit %ld / %ld\n", bit, expx);
			writeResults (buf);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving || stopping) {
			write_time = DISK_WRITE_TIME * 60;
			if (! gmwriteToFile (filename, bit, ubx, uby, x, y)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (buf);
				if (write_time > 600) write_time = 600;
			}
			time (&start_time);


/* If an escape key was hit, write out the results and return */

			if (stopping) {
				gwypfree (tmp);
				gwypfree (tmp2);
				gwypfree(NP);
				gwypfree(M);
				gwypfree(testn);
				if (dovrbareix)
					gwypfree (gx0);
				gwypfree (x);
				gwypfree (y);
				gwypdone ();
				*res = FALSE;
				return (FALSE);
			}
		}

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next iteration so that we can compare our */
/* residues to programs that start counter at zero or one. */

		if (interimResidues && bit % interimResidues < 2) {

			itog (1, tmp2);					// Restore the value of x from the shifted one.
			gshiftleft (ubx, tmp2);
			invg (M,tmp2);
			gtog (M, tmp);
			if (ubx&1)						// View if a sign change on x is necessary.
				subg (tmp2, tmp);
			else
				gtog (tmp2, tmp);

			gwyptogiant (x, tmp2);		// The modulo M reduction is done here
			mulg (tmp2, tmp);
			modg (M, tmp);

			modg (NP, tmp);
			if (!dovrbareix)
				iaddg (-a*a, tmp);		// Compute the (unnormalized) residue modulo NP
			if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
			else if (abs(tmp->sign) < 3)
				sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
			else if (abs(tmp->sign) < 4)
				sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
			else
				sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
			sprintf (buf, "%s interim residue %s at iteration %ld\n", sgk, res64, bit);
			OutputBoth (buf);
		}

/* Write a save file every "interimFiles" iterations. */

		if (interimFiles && bit % interimFiles == 0) {
			char	interimfile[20];
			sprintf (interimfile, "%.8s.%03lu",
				 filename, bit / interimFiles);
			if (! gmwriteToFile (interimfile, bit, ubx, uby, x, y)) {
				sprintf (buf, WRITEFILEERR, interimfile);
				OutputBoth (buf);
			}
		}
	}

	clearline (100);

	itog (1, tmp2);					// Restore the value of x from the shifted one.
	gshiftleft (ubx, tmp2);
	invg (M,tmp2);
	gtog (M, tmp);
	if (ubx&1)						// View if a sign change on x is necessary.
		subg (tmp2, tmp);
	else
		gtog (tmp2, tmp);

	gwyptogiant (x, tmp2);			// The modulo M reduction is done here
	mulg (tmp2, tmp);
	modg (M, tmp);					// Result modulo M


	if (dovrbareix) {
		modg (NP, tmp);
		subg (gx0, tmp);
	}
	else {
		itog (a*a, tmp2);
		invg (NP, tmp2);		// a^(-2) modulo NP
		mulg (tmp2, tmp);		// tmp = a^(2^n-2) = a^(3*(NP-1)) --> the very base is a^3 !!
		modg (NP, tmp);			// Compute the (unnormalized) residue modulo NP
	}

/* Do the Strong PRP test. If the number is proven composite, format a 64-bit residue. */

	if ((!dovrbareix && !isone (tmp)) || (dovrbareix && !isZero (tmp))) {
		*res = FALSE;				/* Not a prime */
		if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
		else if (abs(tmp->sign) < 3)
			sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
		else if (abs(tmp->sign) < 4)
			sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
		else
			sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
		smulg ((unsigned short)a, tmp);
		smulg ((unsigned short)a, tmp);
		smulg ((unsigned short)a, tmp);
		modg (NP, tmp);
		iaddg (-a*a*a, tmp);
		if (abs(tmp->sign) < 2)		// make a 64 bit residue correct !!
			sprintf (oldres64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
		else if (abs(tmp->sign) < 3)
			sprintf (oldres64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
		else if (abs(tmp->sign) < 4)
			sprintf (oldres64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
		else
			sprintf (oldres64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
		if (vrbareix && !dovrbareix)
			if (IniGetInt (INI_FILE, (char*)"OldRes64", 1))
				sprintf (buf, "%s is not prime, although Vrba-Reix PSP!!  RES64: %s.  OLD64: %s", sgk, res64, oldres64);
			else
				sprintf (buf, "%s is not prime, although Vrba-Reix PSP!!  RES64: %s", sgk, res64);
		else if (!vrbareix && !dovrbareix)
			if (IniGetInt (INI_FILE, (char*)"OldRes64", 1))
				sprintf (buf, "%s is not prime.  RES64: %s.  OLD64: %s", sgk, res64, oldres64);
			else
				sprintf (buf, "%s is not prime.  RES64: %s", sgk, res64);
		else if (!vrbareix && dovrbareix)
			sprintf (buf, "%s is not prime, although Strong Fermat PSP!!  Vrba-Reix RES64: %s", sgk, res64);
		else
			sprintf (buf, "%s is not prime.  Vrba-Reix RES64: %s", sgk, res64);
	}
	else if (!dovrbareix) {			// May be a prime, continue the SPRP test
		itog (1, tmp2);				// Restore the value of y from the shifted one.
		gshiftleft (uby, tmp2);
		invg (M,tmp2);
		gtog (M, tmp);
		if (uby&1)					// View if a sign change on y is necessary.
			subg (tmp2, tmp);
		else
			gtog (tmp2, tmp);

		gwyptogiant (y, tmp2);		// The modulo M reduction is done here
		mulg (tmp2, tmp);
		modg (M, tmp);

		modg (NP, tmp);
		iaddg (a, tmp);
		if (gcompg (NP, tmp) != 0 && (tmp->sign != 1 || tmp->n[0] != 2*a)) {
			*res = FALSE;			/* Not a prime */
			if (vrbareix)
				sprintf (buf, "%s is not prime, although Vrba-Reix and Base %lu - Fermat PSP!!", sgk, a*a*a);
			else
				sprintf (buf, "%s is not prime, although Base %lu - Fermat PSP!!", sgk, a*a*a);
		}
		else {
			sprintf (buf, "%s is Base %lu - Strong Fermat PRP! (%d decimal digits)", sgk, a*a*a, nbdg);
		}
	}
	else
		sprintf (buf, "%s is Vrba-Reix PRP! (%d decimal digits)", sgk, nbdg);




#if defined(WIN32) && !defined(_CONSOLE)

	sprintf (buf+strlen(buf), "  Time : "); 
	ReplaceableLine (2);	/* Replace line */ 

#else

	clearline(100);

#ifdef _CONSOLE
	OutputBoth(buf);
#else
	if (*res) {
		OutputStr((char*)"\033[7m");
		OutputBoth(buf);
		OutputStr((char*)"\033[0m");
	}
	else
		OutputBoth(buf);
#endif

	sprintf (buf, "  Time : "); 

#endif

/* Output the final timings */

	gwypend_timer (1);
	gwypwrite_timer (buf+strlen(buf), 1, TIMER_CLR | TIMER_NL); 
	OutputBoth (buf);

/* Cleanup and return */

	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (x);
	gwypfree (y);
	IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
	lasterr_point = 0;
	gwypdone ();
	_unlink (filename);
	if (dualtest && *res)				// If dual test required and positive result
		if (vrbareix && dovrbareix) {
			gwypfree (gx0);
			dovrbareix = FALSE;
			goto restart;				// Do now a Fermat SPRP test
		}
		else if (!vrbareix && !dovrbareix) {
			dovrbareix = TRUE;
			goto restart;				// Do now a Vrba-Reix test
		}
	gwypfree(NP);
	gwypfree(M);
	gwypfree(testn);
	if (dovrbareix)
		gwypfree (gx0);
	return (TRUE);

/* An error occured, sleep, then try restarting at last save point. */

error:
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (x);
	gwypfree (y);
//	gwypdone ();

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
		sprintf (buf, ERRMSG5, checknumber, sgk);
		OutputBoth (buf);
		*res = FALSE;
		gwypfree(NP);
		gwypfree(M);
		gwypfree(testn);
		if (dovrbareix)
			gwypfree (gx0);
		gwypdone ();
		_unlink (filename);
		will_try_larger_fft = FALSE;
		return (TRUE);
	}

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
		gwypdone ();
		will_try_larger_fft = FALSE;
		return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
		OutputBoth (ERRMSG8);
		IniWriteInt(INI_FILE, (char*)"FFT_Increment", IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
		_unlink (filename);
		will_try_larger_fft = FALSE;
	}
	goto restart;

}

static unsigned __int64 li, smallbase, smallk, lastfactor;

int ispoweroftwo (
	unsigned long n)
{
	if (!n)
		return (FALSE);
	while (!(n&1))
		n = n >> 1;
	return (n == 1);
}

int process_num (
	unsigned long format,
	char *sgk,
        char *sgb,
	unsigned long n,
	int	incr,
	unsigned long shift,
	int	*res)
{
	int	retval;
	char outbuf[sgkbufsize+256];

// Lei -remove a line and replace
//	unsigned long ninput = n, binput = base;
	unsigned long ninput = n, binput, base, b_2up = 1, b_else = 1, superPRP = 1;
	long mult;
// Lei end

	gformat = format; // save format in a global.

	if(mult = IniGetInt(INI_FILE, (char*)"StopOnPrimedK", 0)) {
		sprintf (outbuf, "ks%s", sgk);
		if(IniGetInt(INI_FILE, outbuf, 0) >= mult) {// is the count for this k value reached ?
			*res = FALSE;							// then, skip this test
			return TRUE;
		}
	}
	else if(mult = IniGetInt(INI_FILE, (char*)"StopOnPrimedN", 0)) {
		sprintf (outbuf, "ns%lu", n);
		if(IniGetInt(INI_FILE, outbuf, 0) >= mult) {// is the count for this n value reached ?
			*res = FALSE;							// then, skip this test
			return TRUE;
		}
	}
	else if(mult = IniGetInt(INI_FILE, (char*)"StopOnPrimedB", 0)) {
		sprintf (outbuf, "bs%s", sgb);
		if(IniGetInt(INI_FILE, outbuf, 0) >= mult) {// is the count for this base value reached ?
			*res = FALSE;							// then, skip this test
			return TRUE;
		}
	}

	if (format == ABCGM)
		return (isGMNP (sgk, n, res));	// Do the primality test of a Gaussian Mersenne norm

	if (format == ABCSP)   // Do the PRP test of a Wagstaff number
		return (isWSPRP (sgk, n, res));

	gb = newgiant (strlen(sgb)/2 + 8); // Allocate one byte per decimal digit + spares
	ctog (sgb, gb);	        // Convert b string to giant
	if (gb->sign <= 2) {	// Test if the base is a small integer...
            base = gb->n[0];// Then, get the base in an unsigned long
            if (gb->sign == 2)
                base += 65536*gb->n[1];
            binput = base;
            while (!(base&1) && base > 2) { // Divide the base by two as far as possible
		base >>= 1;
		n += ninput;
            }
            

            if (base != 2) {	// Test if the base was a power of two

// Lei
		n -= ninput;
		b_else = base;	// Is odd...
                b_2up = binput / b_else;// binput = b_else*b_2up
                if ((b_2up > b_else) && (!((format == ABCC) || (format == ABCK) || (format == ABCRU) || (format == ABCGRU) || (format == ABCVARAQS)))) {
                    superPRP = 0;   // Then b_2up^n > b_else^n
		}
		else {
// Lei end

                    base = binput;  // Do not modify because PRP will be forced...
                    n = ninput;
		}
            }
//	}

		globalb = base;	    // Keep the base of the candidate in a global

//	Replaced by Lei :
//	if (base == 2 && !IniGetInt (INI_FILE, (char*)"ForcePRP", 0) && ((incr == -1) || (incr == +1))) {
//		if (incr == -1)
//			retval = isLLRP (format, sgk, n, shift, res);
//		else
//			retval = isProthP (format, sgk, n, shift, res);
//	}

// Lei mod
	if (((base == 2) || (superPRP == 0)) && !IniGetInt (INI_FILE, (char*)"ForcePRP", 0) && ((incr == -1) || (incr == +1)) && (format != ABCVARAQS)) {
		if (incr == -1)
			retval = IniGetInt(INI_FILE, (char*)"TestW", 0) ? isLLRW (format, sgk, n, shift, res) : isLLRP (format, sgk, b_else, n, binput, ninput, shift, res);
		else
			retval = IniGetInt(INI_FILE, (char*)"TestW", 0) ? isProthW (format, sgk, n, shift, res) : isProthP (format, sgk, b_else, n, binput, ninput, shift, res);
	}
// end Lei mod

	else if ((format == NPGCC1 || format == NPGCC2) && !IniGetInt (INI_FILE, (char*)"ForcePRP", 0)) {
		retval = IsCCP (format, sgk, base, n, incr, shift, res);
	}
	else if (!IniGetInt (INI_FILE, (char*)"ForcePRP", 0) && (incr == +1 || incr == -1) && (format != ABCVARAQS) && 
		(format != ABCRU) && (format != ABCGRU))
		retval = plusminustest (sgk, base, n, incr, shift, res);
	else  {
		retval = IsPRP (format, sgk, base, n, incr, shift, res);
	}
//	return (retval);
        }   // End gb is a small integer.
	else if ((format == NPGCC1 || format == NPGCC2) && !IniGetInt (INI_FILE, (char*)"ForcePRP", 0)) {
		retval = gIsCCP (format, sgk, sgb, gb, n, incr, shift, res);
	}
	else if (!IniGetInt (INI_FILE, (char*)"ForcePRP", 0) && (incr == +1 || incr == -1) && (format != ABCVARAQS) && 
	(format != ABCRU) && (format != ABCGRU))
		retval = gplusminustest (sgk, sgb, gb, n, incr, shift, res);
	else  {
//		Fermat_only = TRUE;     // JP 30/01/17
		retval = gIsPRP (format, sgk, sgb, gb, n, incr, shift, res);
//		Fermat_only = FALSE;    // JP 30/01/17
        }
	free (gb);
	return (retval);
}

char	outpf[] = "gqplus.res", outmf[] = "gqminus.res";
char	gqpstring[] = "ABC (2^$a+2^(($a+1)/2)+1)/5\n";
char	gqmstring[] = "ABC (2^$a-2^(($a+1)/2)+1)/5\n";


int primeContinue ()
{

	int	work, nargs, hiline, completed = FALSE;
	unsigned long format, shift, begline, rising_ns, rising_ks, last_processed_n;
	char *pinput;

/* Set appropriate priority */

	SetPriority ();

/* Case off the work type */

	work = IniGetInt (INI_FILE, (char*)"Work", 0);

/* Handle a sieving program output file */

	if (work == 0) {
	    char	inputfile[80], outputfile[80], oldinputfile[80], cmaxroundoff[10], cpcfftlim[10], sgk[sgkbufsize], buff[sgkbufsize+256];
		char	hbuff[sgkbufsize+256], outbuf[sgkbufsize+256], last_processed_k[sgkbufsize+256];
	    FILE *fd;
	    unsigned long i, chainlen, m, n, base, nfudge, nn;
	    int	firstline, line, hline, resultline,
			outfd, outfdp, outfdm, res, incr, sign, argcnt, validheader = FALSE;
	    char c;

#ifdef	WIN32
		giant initgiants = newgiant (1<<19);	// Create a giant of maximal size
#else
		giant initgiants = newgiant (-1);	// Create a giant of maximal size
#endif
		gwypfree (initgiants);					// And free it, to initialize the popg / pushg routines


	    IniGetString (INI_FILE, (char*)"PgenInputFile", inputfile, IBSIZE, NULL);
	    IniGetString (INI_FILE, (char*)"OldInputFile", oldinputfile, IBSIZE, inputfile);	// Default it to PgenInputFile! JP 26/02/17
	    IniGetString (INI_FILE, (char*)"PgenOutputFile", outputfile, IBSIZE, NULL);
	    IniGetString (INI_FILE, (char*)"MaxRoundOff", cmaxroundoff, 5, (char*)"0.40");
		maxroundoff = atof (cmaxroundoff);
	    IniGetString (INI_FILE, (char*)"PercentFFTLimit", cpcfftlim, 5, (char*)"0.50");
            pcfftlim = atof (cpcfftlim);
            if (!strcmp (inputfile, oldinputfile))
                firstline = IniGetInt (INI_FILE, (char*)"PgenLine", 1);		// Continuing on the same file
            else
                firstline = 1;		
                // Processing a new file
            last_processed_n = (unsigned long)IniGetInt(INI_FILE, (char*)"Last_Processed_n", 0);
            IniGetString(INI_FILE, (char*)"Last_Processed_k",last_processed_k, sgkbufsize, NULL);
//	    firstline = IniGetInt (INI_FILE, (char*)"PgenLine", 1);
	    hline = IniGetInt (INI_FILE, (char*)"HeaderLine", 0);
	    verbose = IniGetInt (INI_FILE, (char*)"Verbose", 0);

// Transmit the pointers to user output fuctions to the gwypnum system.

		gwypsetoutputs (OutputStr, OutputBoth);

	    setuponly = IniGetInt (INI_FILE, (char*)"SetupOnly", 0);
/*	    fd = fopen (inputfile, "r");

	    if (fd == NULL) {
			IniWriteInt (INI_FILE, (char*)"Workdone", 1);
			return;
	    } */

		begline = IniGetInt(INI_FILE, (char*)"BegLine", 0);
		testgm  = IniGetInt(INI_FILE, (char*)"TestGM", 1);
		testgq  = IniGetInt(INI_FILE, (char*)"TestGQ", 0);
		testfac  = IniGetInt(INI_FILE, (char*)"TestFac", 0);
		facfrom =  IniGetInt(INI_FILE, (char*)"FacFrom", 0);
		facto =  IniGetInt(INI_FILE, (char*)"FacTo", 0);
		debug =  IniGetInt(INI_FILE, (char*)"Debug", 0);
		zcomplex =  IniGetInt(INI_FILE, (char*)"Zcomplex", 1);
		vrbareix  = IniGetInt(INI_FILE, (char*)"VrbaReixTest", 0);
		dualtest = IniGetInt(INI_FILE, (char*)"DualTest", 0);
		hiline =  IniGetInt(INI_FILE, (char*)"HiLine", 0);
		rising_ns =  IniGetInt(INI_FILE, (char*)"Rising_ns", 0);
		rising_ks =  IniGetInt(INI_FILE, (char*)"Rising_ks", 0);
		nofac =  IniGetInt(INI_FILE, (char*)"NoPrefactoring", 0);
		showdigits =  IniGetInt(INI_FILE, (char*)"ShowDigits", 0);

/* A new option to create interim save files every N iterations. */
/* This allows two machines to simultanously work on the same exponent */
/* and compare results along the way. */

		interimFiles = IniGetInt (INI_FILE, (char*)"InterimFiles", 0);
		interimResidues = IniGetInt (INI_FILE, (char*)"InterimResidues", interimFiles);

/* Option to slow down the program by sleeping after every iteration.  You */
/* might use this on a laptop or a computer running in a hot room to keep */
/* temperatures down and thus reduce the chance of a hardware error.  */

		throttle = IniGetInt (INI_FILE, (char*)"Throttle", 0);

OPENFILE :
                fd = fopen (inputfile, "r");

                if (fd == NULL) {
			IniWriteInt (INI_FILE, (char*)"WorkDone", 1);
			return (FALSE);
                }

		if (!strncmp (buff, (char*)"TestWieferichcode", 17)) {	// Very particular test code...
			TestWieferich ();
			IniWriteInt (INI_FILE, (char*)"Workdone", 1);
			return (FALSE);
                }

		sprintf (SVINI_FILE, "save_%s", INI_FILE);		// set the name of the backup Ini File

// Process each line in the output file

		for (line=0; ; line++) {

// Set termination on error conditions

		abonillsum = IniGetInt(INI_FILE, (char*)"Abortonillsum", 0);
		abonmismatch = IniGetInt(INI_FILE, (char*)"Abortonmismatch", 0);
		abonroundoff = IniGetInt(INI_FILE, (char*)"Abortonroundoff", 0);
		if (IniGetInt(INI_FILE, (char*)"Abortonerror", 0))
			abonillsum = abonmismatch = abonroundoff = 1;

                generic =  IniGetInt(INI_FILE, (char*)"ForceGeneric", 0);
                cufftonly = IniGetInt(INI_FILE, (char*)"Cufftonly", 0);
                    
// Blank the input line

			for (i=0; i<(sgkbufsize+256); i++)
				buff[i] = ' ';
			buff[sgkbufsize+255] = '\n';

// Read the line, break at EOF

			if (fgets (buff, sgkbufsize+256, fd) == NULL) {
				IniWriteInt (INI_FILE, (char*)"Workdone", 1);
				rising_ns = rising_ks = FALSE;
				break;
			}
			else
				IniWriteInt (INI_FILE, (char*)"Workdone", 0);

// Skip this line if requested (we processed it on an earlier run)
// (but don't ignore last header line found!)

			if (hiline && line > hiline) {
				IniWriteInt (INI_FILE, (char*)"Workdone", 1);
				break;
			}

			if (!strncmp (buff, "ABC", 3)) {	// ABC format header found

				sprintf (hbuff, "%s", buff);	// Save the header
				IniWriteInt (INI_FILE, (char*)"HeaderLine", line);	// Save the header line number
				hline = line;
				validheader = TRUE;				// Assume it is valid...

				for (pinput=buff+3; *pinput && isspace(*pinput); pinput++);

				for (i=0;i<strlen(pinput);i++)
					if (isspace(pinput[i]))
						pinput[i] = '\0';		// Suppress the EOL characters if necessary.

				if (!strcmp (pinput, grepustring)) {
					format = ABCGRU;
				}
				else if (!strcmp (pinput, abcadstring)) {
					format = ABCVARAQS;
				}
				else if (!strcmp (pinput, diffnumstring))
					format = ABCDNG;
				else if (!strcmp (pinput, ckstring)) {
					format = ABCK;
				}
				else if (!strcmp (pinput, repustring)) {
					format = ABCRU;
				}
				else if (!strcmp (pinput, cwstring)) {
					format = ABCCW;
				}
				else if (!strcmp (pinput, abcastring)) {
					format = ABCVARAS;
				}
				else if (!strcmp (pinput, ffstring)) {
					format = ABCFF;
				}
				else if (!strcmp (pinput, gmstring)) {
					format = ABCGM;
				}
				else if (!strcmp (pinput, spstring)) {
					format = ABCSP;
				}
				else if (!strcmp (pinput, gfstring)) {
					format = ABCGF;
				}
				else if (!strcmp (pinput, wfsstring)) {
					format = ABCWFS;
				}
				else if (!strcmp (pinput, wftstring)) {
					format = ABCWFT;
				}
				else if (!strcmp (pinput, numberstring)) {
					format = ABCGPT;
				}
				else if (sscanf(pinput, fnpstring, &n, &incr) == 2) {
					format = ABCFNGS;
				}
				else if (sscanf(pinput, fnmstring, &n, &incr) == 2) {
					format = ABCFNGS;
					incr = - incr;
				}
				else if (sscanf(pinput, fnpstring, &n) == 1) { 
					format = ABCFNAS;
				}
				else if (sscanf(pinput, abcpstring, &incr) == 1) {
					format = ABCVARGS;
				}
				else if (sscanf(pinput, abcmstring, &incr) == 1) {
					format = ABCVARGS;
					incr = - incr;
				}
				else if (sscanf(pinput, diffnumpstring, &incr) == 1) {
					format = ABCDN;
				}
				else if (sscanf(pinput, diffnummstring, &incr) == 1) {
					format = ABCDN;
					incr = - incr;
				}
				else if (sscanf(pinput, fkpstring, &smallk, &incr) == 2) {
					sprintf(sgk, $LLF, smallk);
                                        // unsigned fixed k...	
					format = ABCFKGS;
				}
				else if (sscanf(pinput, fkmstring, &smallk, &incr) == 2) {
					sprintf(sgk, $LLF, smallk);
                                        // unsigned fixed k...	
					format = ABCFKGS;
					incr = - incr;
				}
				else if (sscanf(pinput, fkpstring, &smallk) == 1) { 
					sprintf(sgk, $LLF, smallk);
                                        // unsigned fixed k...	
					format = ABCFKAS;
				}
				else if (sscanf(pinput, fbpstring, &smallbase, &incr) == 2) {
					sprintf (sgb, $LLF, smallbase);	// unsigned fixed base...	
					format = ABCFBGS;
				}
				else if (sscanf(pinput, fbmstring, &smallbase, &incr) == 2) {
					sprintf (sgb, $LLF, smallbase);	// unsigned fixed base...	
					format = ABCFBGS;
					incr = - incr;
				}
				else if (sscanf(pinput, fbastring, &smallbase) == 1) { 
					sprintf (sgb, $LLF, smallbase);	// unsigned fixed base...	
					format = ABCFBAS;
				}
				else {
					OutputBoth ((char*)"Invalid ABC format, next data lines will be flushed...\n");
					validheader = FALSE;		// Invalid header found...
				}

				if (format == ABCGM) {
					if (!facto)
						sprintf (pinput+strlen (gmstring),
							" // Let GM(p) = (1+/-i)^p-1, GQ(p) = ((1+/-i)^p+1)/(2+/-i) if p>3, (1+/-i)^p+1 if p<=3\n");
					if (!facto && !fileExists (outpf)) {
						outfdp = _open (outpf, _O_TEXT | _O_RDWR | _O_CREAT, 0666);
						if (outfdp) {
							_write (outfdp, gqpstring, strlen (gqpstring));
							_close (outfdp);
						}	
					}
					if (!facto && !fileExists (outmf)) {
						outfdm = _open (outmf, _O_TEXT | _O_RDWR | _O_CREAT, 0666);
						if (outfdm) {
							_write (outfdm, gqmstring, strlen (gqmstring));
							_close (outfdm);
						}
					}
				}
				continue;				// Read next line, but do not change PgenLine!
			}							// End ABC format header found

			else if (((argcnt = sscanf (buff, $LLF":%c:%lu:"$LLF":%lu\n", &li, &c, &chainlen, &smallbase, &mask)) > 1) || !line) {
				if (argcnt < 4) {
					OutputBoth ((char*)"Missing or invalid NewPGen header, next data lines will be flushed...\n");
					validheader = FALSE;			// Invalid NewPGen header...
				}
				else {
					sprintf (sgb, $LLF, smallbase);	// Newpgen format admits only unsigned base...	
					validheader = TRUE;
					if (argcnt == 4)
						mask = 0;
					strcpy (hbuff, buff);			// Save the header
					IniWriteInt (INI_FILE, (char*)"HeaderLine", line);	// Save the header line number
					hline = line;
					format = NPG;
					if (mask & 0x40) {
						OutputStr ((char*)"Primorial NewPgen files are not supported...\n");
						validheader = FALSE;
					}
					if (chainlen == 0) chainlen = 1;
				}
				continue;				// Read next line, but do not change PgenLine!
			}							// End NewPGen header found

			else {						// Processing a data line
				if (((!rising_ns && !rising_ks) || (rising_ns && rising_ks)) && (line < firstline))
					continue;			// Skip this line if requested (we processed it on an earlier run)

				if (!validheader)
					continue;			// Flush data until a valid header is found...

				shift = 0;				// Only one value for the k multiplier

				if (format == NPG) {	// NEWPGEN output

// THIS SECTION IS FOR BACKWARDS COMPATIBILITY WITH PREVIOUS PRP.EXE
// That version used the one character code to determine what to do.
// The new version uses the mask field.

					if (mask == 0 || IniGetInt (INI_FILE, (char*)"UseCharCode", 0)) {

// The variable c is a one character code as follows:
//  P : k.b^n+1 (Plus)
//  M : k.b^n-1 (Minus)
//  T: k.b^n+-1 (Twin)
//  S: k.b^n-1; k.b^(n+1)-1 (SG (CC 1st kind len 2))
//  C: k.b^n+1; k.b^(n+1)+1 (CC 2nd kind len 2)
//  B: k.b^n+-1; k.b^(n+1)+-1 (BiTwin)
//  J: k.b^n+-1; k.b^(n+1)-1 (Twin/SG)
//  K: k.b^n+-1; k.b^(n+1)+1 (Twin/CC)
//  Y : k.b^n+1 + others (Lucky Plus)
//  Z : k.b^n-1 + others (Lucky Minus)
//  1: CC 1st kind chain
//  2: CC 2nd kind chain
//  3: BiTwin chain
// Undo the increment of n that newpgen did on types 1, 2, 3
// Map P, M, Y, Z, T, S, C, B to their more generic counterparts

						nfudge = 0;
						if (c == '1') nfudge = 1;
						if (c == '2') nfudge = 1;
						if (c == '3') nfudge = 1;
						if (c == 'P') c = '2', chainlen = 1;
						if (c == 'M') c = '1', chainlen = 1;
//						if (c == 'Y') c = '2', chainlen = 1;
//						if (c == 'Z') c = '1', chainlen = 1;
						if (c == 'T') c = '3', chainlen = 1;
						if (c == 'S') c = '1', chainlen = 2;
						if (c == 'C') c = '2', chainlen = 2;
						if (c == 'B') c = '3', chainlen = 2;


// Process each line in the newpgen output file

// allow k to be a big integer
						if (sscanf (buff+begline, "%s %lu", sgk, &n) != 2)
							continue;				// Skip invalid line

						if (!isDigitString(sgk))
							continue;				// Skip invalid line

						if (rising_ns && !rising_ks && (n <= last_processed_n))
							continue;				// Skip already processed n's


						if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
							continue;				// Skip already processed k's

						if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
							fclose (fd);			// Unlock the file during the test...
                                                        
// Test numbers according to the c variable

						nn = n;
						if (c == 'Y') {
							nn--;
						}
						if (c == 'Z') {
							nn--;
						}

						for (i = 0; i < chainlen; i++) {
							if (c == '1' || c == '3') {
								if (! process_num (format, sgk, sgb, n - nfudge + i, -1, shift, &res))
									goto done;
								if (!res)
									break;
								if (c == '1')
									format = NPGCC1;
							}
							if (c == '2' || c == '3') {
								if (! process_num (format, sgk, sgb, n - nfudge + i, +1, shift, &res))
									goto done;
								if (!res)
									break;
								if (c == '2')
									format = NPGCC2;
							}
							if (c == 'J') {	// Twin/SG
								int	res2;
								if (! process_num (format, sgk, sgb, n, -1, shift, &res))
									goto done;
								if (!res)
									break;
								if (! process_num (format, sgk, sgb, n+1, -1, shift, &res))
									goto done;
								if (! process_num (format, sgk, sgb, n, +1, shift, &res2))
									goto done;
								res |= res2;
								format = NPG;
								break;
							}
							if (c == 'K') {	// Twin/CC
								int	res2;
								if (! process_num (format, sgk, sgb, n, +1, shift, &res))
									goto done;
								if (!res)
									break;
								if (! process_num (format, sgk, sgb, n, -1, shift, &res))
									goto done;
								if (! process_num (format, sgk, sgb, n+1, +1, shift, &res2))
									goto done;
								res |= res2;
								format = NPG;
								break;
							}
							if (c == 'Y') {	// Lucky Plus
								int	res2;
								if (! process_num (format, sgk, sgb, nn+1, +1, shift, &res))
									goto done;
								if (!res)
									break;
								if (! process_num (format, sgk, sgb, nn+1, -1, shift, &res))
									goto done;
								if (! process_num (format, sgk, sgb, nn, +1, shift, &res2))
									goto done;
								res |= res2;
								if (! process_num (format, sgk, sgb, nn+2, +1, shift, &res2))
									goto done;
								res |= res2;
								format = NPG;
								break;
							}
							if (c == 'Z') {	// Lucky Minus
								int	res2;
								if (! process_num (format, sgk, sgb, nn+1, -1, shift, &res))
									goto done;
								if (!res)
									break;
								if (! process_num (format, sgk, sgb, nn+1, +1, shift, &res))
									goto done;
								if (! process_num (format, sgk, sgb, nn, -1, shift, &res2))
									goto done;
								res |= res2;
								if (! process_num (format, sgk, sgb, nn+2, -1, shift, &res2))
									goto done;
								res |= res2;
								format = NPG;
								break;
							}
							if (c == 'A') {	// AP mode
								format = NPGAP;
								if (! process_num (format, sgk, sgb, n, -1, shift, &res))
									goto done;
								format = NPG;
								if (!res) break;
							}
						}		// End loop on chain length
						format = NPG;
					}			// End of old section


// THIS IS THE NEW SECTION.  It uses both the mask field and the
// character code to determine what to do

					else {

// NEWPGEN output files use the mask as defined below:
// #define MODE_PLUS    0x01	/* k.b^n+1
// #define MODE_MINUS   0x02	/* k.b^n-1
// #define MODE_2PLUS   0x04	/* k.b^(n+1)+1 (*)
// #define MODE_2MINUS  0x08	/* k.b^(n+1)-1 (*)
// #define MODE_4PLUS   0x10	/* k.b^(n+2)+1 (*)
// #define MODE_4MINUS  0x20	/* k.b^(n+2)-1 (*)
// #define MODE_PRIMORIAL 0x40	/* PRIMORIAL - can't handle this
// #define MODE_PLUS5  0x80	/* k.b^n+5
// #define MODE_2MINUS3 0x100	/* 2k.b^n-3 JP 23/08/17 */
// #define MODE_AP	    0x200	/* 2^n+2k-1
// #define MODE_PLUS7  0x800	/* k.b^n+7
// #define MODE_2PLUS3 0x1000	/* 2k.b^n+3
// #define MODE_DUAL 0x8000
// #define MODE_PLUS_DUAL 0x8001	/* b^n+k
// #define MODE_MINUS_DUAL 0x8002	/* b^n-k
// #define MODE_NOTGENERALISED 0x400
// Those entries that have a (*) next to them are modified if the
// MODE_NOTGENERALISED flag is set.  If it is set, they are changed
// as follows
// MODE_2PLUS      2k.b^n+1
// MODE_2MINUS     2k.b^n-1
// MODE_4PLUS      4k.b^n+1
// MODE_4MINUS     4k.b^n-1
// Similarly, longer chains are affected in the same way (so if the base
// is 3 and we are after a CC of the 1st kind of length 4, rather that
// looking at k.3^n-1 & k.3^(n+1)-1 & k.3^(n+2)-1 & k.3^(n+3)-1 we look
// at k.3^n-1 & 2k.3^n-1 & 4k.3^n-1 & 8k.3^n-1).

// allow k to be a big integer

						if (sscanf (buff+begline, "%s %lu", sgk, &n) != 2)
							continue;	// Skip invalid line

						if (!isDigitString(sgk))
							continue;	// Skip invalid line


                                                if (rising_ns && !rising_ks && (n <= last_processed_n))
							continue;				// Skip already processed n's

						if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
							continue;				// Skip already processed k's

						if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
							fclose (fd);			// Unlock the file during the test...

// Undo the increment of n that newpgen did on types 1, 2, 3

						nn = n;
//						if (c == '1' || c == '2' || c == '3')
//							nn--;

						if (c == 'S')
							chainlen = 2;
						if (c == 'C')
							chainlen = 2;
						if (c == 'B')
							chainlen = 2;

						if ((mask & MODE_PLUS) && (mask & MODE_2MINUS) &&
								(mask & MODE_2PLUS) && (mask & MODE_4PLUS)) {
							nn--;
						}
						if ((mask & MODE_MINUS) && (mask & MODE_2MINUS) &&
								(mask & MODE_2PLUS) && (mask & MODE_4MINUS)) {
							nn--;
						}

// Test numbers according to the mask variable
// The J and K types (Twin/CC and Twin/SG) are special in that they
// are output if either a Twin OR a CC/SG is found

						shift = 0;

						for (i = 0; i < chainlen; i++) {
							if ((mask & MODE_MINUS) && (mask & MODE_PLUS) &&
								(mask & MODE_2MINUS)) {	// Twin/SG
								int	res2;
								if (! process_num (format, sgk, sgb, nn, -1, shift, &res))
									goto done;
								if (!res)
									break;
								if (! process_num (format, sgk, sgb, nn+1, -1, shift, &res))
									goto done;
								if (! process_num (format, sgk, sgb, nn, +1, shift, &res2))
									goto done;
								res |= res2;
								break;
							}
							if ((mask & MODE_MINUS) && (mask & MODE_PLUS) &&
								(mask & MODE_2PLUS)) {	// Twin/CC
								int	res2;
								if (! process_num (format, sgk, sgb, nn, +1, shift, &res))
									goto done;
								if (!res)
									break;
								if (! process_num (format, sgk, sgb, nn, -1, shift, &res))
									goto done;
								if (! process_num (format, sgk, sgb, nn+1, +1, shift, &res2))
									goto done;
								res |= res2;
								break;
							}
							if ((mask & MODE_PLUS) && (mask & MODE_2MINUS) &&
								(mask & MODE_2PLUS) && (mask & MODE_4PLUS)) {	// Lucky Plus
								int	res2;
								if (! process_num (format, sgk, sgb, nn+1, +1, shift, &res))
									goto done;
								if (!res)
									break;
								if (! process_num (format, sgk, sgb, nn+1, -1, shift, &res))
									goto done;
								if (! process_num (format, sgk, sgb, nn, +1, shift, &res2))
									goto done;
								res |= res2;
								if (! process_num (format, sgk, sgb, nn+2, +1, shift, &res2))
									goto done;
								res |= res2;
								break;
							}
							if ((mask & MODE_MINUS) && (mask & MODE_2MINUS) &&
								(mask & MODE_2PLUS) && (mask & MODE_4MINUS)) {	// Lucky Minus
								int	res2;
								if (! process_num (format, sgk, sgb, nn+1, -1, shift, &res))
									goto done;
								if (!res)
									break;
								if (! process_num (format, sgk, sgb, nn+1, +1, shift, &res))
									goto done;
								if (! process_num (format, sgk, sgb, nn, -1, shift, &res2))
									goto done;
								res |= res2;
								if (! process_num (format, sgk, sgb, nn+2, -1, shift, &res2))
									goto done;
								res |= res2;
								break;
							}
							if (mask & MODE_MINUS) {
								if (mask & MODE_DUAL) {
									if (! process_num (format, (char*)"1", sgb, nn, -atoi(sgk), shift, &res))
										goto done;
								}
								else
									if (! process_num (format, sgk, sgb, nn, -1, shift, &res))
										goto done;
								if (!res)
									break;
							}
							if (mask & MODE_PLUS) {
								if (mask & MODE_DUAL) {
									if (! process_num (format, (char*)"1", sgb, nn, atoi(sgk), shift, &res))
										goto done;
								}
								else
									if (! process_num (format, sgk, sgb, nn, +1, shift, &res))
										goto done;
								if (!res)
									break;
							}
							if (mask & MODE_PLUS5) {
								if (! process_num (format, sgk, sgb, nn, +5, shift, &res))
									goto done;
								if (!res)
									break;
							}
							if (mask & MODE_PLUS7) {
								if (! process_num (format, sgk, sgb, nn, +7, shift, &res))
									goto done;
								if (!res)
									break;
							}
							if (mask & MODE_2PLUS3) {
								shift = 1;
								format = NPGCC1;
								if (! process_num (format, sgk, sgb, nn, +3, shift, &res))
									goto done;
								shift = 0;
								format = NPG;
								if (!res)
									break;
							}
							if (mask & MODE_2MINUS3) {
								shift = 1;
								format = NPGCC2;
								if (! process_num (format, sgk, sgb, nn, -3, shift, &res))
									goto done;
								shift = 0;
								format = NPG;
								if (!res)
									break;
							}
							if (mask & MODE_AP) {
								format = NPGAP;
								if (! process_num (format, sgk, sgb, nn, -1, shift, &res))
									goto done;
								format = NPG;
								if (!res)
									break;
							}

// Bump k or n for the next iteration or for the MODE_2PLUS and
// MODE_2MINUS flags

							if (mask & MODE_NOTGENERALISED)
								shift += 1; 
							else
								nn += 1;

// If chainlength is more than 1, then we let the for loop do the work
// rather than the MODE_2PLUS, etc. flags

							if (chainlen > 1) {
								if ((mask & MODE_2MINUS) || (mask & MODE_4MINUS))
									format = NPGCC1;
								else if ((mask & MODE_2PLUS) || (mask & MODE_4PLUS))
									format = NPGCC2;
								else
									format = NPG;
								continue;
							}
							if (mask & MODE_2MINUS) {
								if (! process_num (format, sgk, sgb, nn, -1, shift, &res))
									goto done;
								if (!res)
									break;
							}
							if (mask & MODE_2PLUS) {
								if (! process_num (format, sgk, sgb, nn, +1, shift, &res))
									goto done;
								if (!res)
									break;
							}

// Bump k or n for the MODE_4PLUS and MODE_4MINUS flags

							if (mask & MODE_NOTGENERALISED)
								shift += 1; 
							else
								nn += 1;

							if (mask & MODE_4MINUS) {
								if (! process_num (format, sgk, sgb, nn, -1, shift, &res))
									goto done;
								if (!res)
									break;
							}
							if (mask & MODE_4PLUS) {
								if (! process_num (format, sgk, sgb, nn, +1, shift, &res))
									goto done;
								if (!res)
									break;
							}
						}	// End of loop on chain length
						format = NPG;
					}		// End of new section

// If all numbers tested were primes or PRPs, copy the line to the output file

					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%s %lu\n", sgk, n);	// write the result
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}			// End of NewPGen format processing

				else if (format == ABCCW) {			// Cullen/Woodall
                                        if (sscanf (buff+begline, "%lu %s %d", &n, sgb, &incr) != 3)
						continue;				// Skip invalid line
					if (!isDigitString (sgb))
						continue;				// Skip invalid line
					if (rising_ns && !rising_ks  && (n <= last_processed_n))
						continue;				// Skip already processed n's
					sprintf (sgk, "%lu", n);
					if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
						continue;				// Skip already processed k's
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					sprintf (sgk, "%lu", n);
					if (! process_num (format, sgk, sgb, n, incr, shift, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%lu %lu %d\n", n, base, incr); 
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCFF)	{	// FermFact output
												// allow k to be a big integer
					if (sscanf (buff+begline, "%s %lu", sgk, &n) != 2)
						continue;				// Skip invalid line
					if (!isDigitString(sgk))
						continue;				// Skip invalid line
					if (rising_ns && !rising_ks && (n <= last_processed_n))
						continue;				// Skip already processed n's
					if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
						continue;				// Skip already processed k's
					sprintf (sgb, "2");
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					if (! process_num (format, sgk, (char*)"2", n, +1, shift, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%s %lu\n", sgk, n); 
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCLEI)       {	// Lei output
											// allow k to be a big integer
					if (sscanf (buff+begline, "%s %lu", sgk, &n) != 2)
						continue;			// Skip invalid line
					if (!isDigitString(sgk))
						continue;			// Skip invalid line
					if (rising_ns && !rising_ks && (n <= last_processed_n))
						continue;				// Skip already processed n's
					if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
						continue;				// Skip already processed k's
					sprintf (sgb, "2");
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					if (! process_num (format, sgk, (char*)"2", n, -1, shift, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%s %lu\n", sgk, n);
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCFKGS)	{	// Fixed k:  b and n specified on each input line
					if (sscanf (buff+begline, "%s %lu", sgb, &n) != 2)
						continue;	
                                                        // Skip invalid line
					if (!isDigitString (sgb))
						continue;				// Skip invalid line
					if (rising_ns && (n <= last_processed_n))
						continue;	
                                                // Skip already processed n's
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					if (! process_num (format, sgk, sgb, n, incr, shift, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%s %lu\n", sgb, n); 
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCFKAS)	{	// Fixed k:  b, n, and c specified on each input line
					if (sscanf (buff+begline, "%s %lu %d", sgb, &n, &incr) != 3)
						continue;				// Skip invalid line
					if (!isDigitString(sgk))
						continue;				// Skip invalid line
					if (!isDigitString (sgb))
						continue;				// Skip invalid line
					if (rising_ns && (n <= last_processed_n))
						continue;				// Skip already processed n's
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					if (! process_num (format, sgk, sgb, n, incr, shift, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%s %lu %d\n", sgb, n, incr); 
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCFBGS)	{	// Fixed b:  k and n specified on each input line
					if (sscanf (buff+begline, "%s %lu", sgk, &n) != 2)
						continue;				// Skip invalid line
					if (!isDigitString(sgk))
						continue;				// Skip invalid line
					if (rising_ns && !rising_ks && (n <= last_processed_n))
						continue;				// Skip already processed n's
					if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
						continue;				// Skip already processed k's
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					if (! process_num (format, sgk, sgb, n, incr, shift, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%s %lu\n", sgk, n); 
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCFBAS)	{	// Fixed b:  k, n, and c specified on each input line
					if (sscanf (buff+begline, "%s %lu %d", sgk, &n, &incr) != 3)
						continue;				// Skip invalid line
					if (!isDigitString(sgk))
						continue;				// Skip invalid line
					if (rising_ns && !rising_ks && (n <= last_processed_n))
						continue;				// Skip already processed n's
					if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
						continue;				// Skip already processed k's
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					if (! process_num (format, sgk, sgb, n, incr, shift, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%s %lu %d\n", sgk, n, incr); 
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCFNGS)	{	// Fixed n:  k and b specified on each input line
					if (sscanf (buff+begline, "%s %s", sgk, sgb) != 2)
						continue;				// Skip invalid line
					if (!isDigitString(sgk))
						continue;				// Skip invalid line
					if (!isDigitString (sgb))
						continue;				// Skip invalid line
					if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
						continue;				// Skip already processed k's
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					if (! process_num (format, sgk, sgb, n, incr, shift, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%s %s\n", sgk, sgb); 
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCFNAS)	{	// Fixed n:  k, b, and c specified on each input line
					if (sscanf (buff+begline, "%s %s %d", sgk, sgb, &incr) != 3)
						continue;				// Skip invalid line
					if (!isDigitString(sgk))
						continue;				// Skip invalid line
					if (!isDigitString (sgb))
						continue;				// Skip invalid line
					if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
						continue;				// Skip already processed k's
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					if (! process_num (format, sgk, sgb, n, incr, shift, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%s %s %d\n", sgk, sgb, incr); 
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCVARGS)	{	// k, b, and n specified on each input line
					if (sscanf (buff+begline, "%s %s %lu", sgk, sgb, &n) != 3)
						continue;				// Skip invalid line
					if (!isDigitString(sgk))
						continue;				// Skip invalid line
					if (!isDigitString (sgb))
						continue;				// Skip invalid line
					if (rising_ns && !rising_ks && (n <= last_processed_n))
						continue;				// Skip already processed n's
					if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
						continue;				// Skip already processed k's
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					if (! process_num (format, sgk, sgb, n, incr, shift, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%s %s %lu\n", sgk, sgb, n); 
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCVARAS)	{	// k, b, n, and c specified on each input line
					if (sscanf (buff+begline, "%s %s %lu %d", sgk, sgb, &n, &incr) != 4)
						continue;				// Skip invalid line
					if (!isDigitString(sgk))
						continue;				// Skip invalid line
					if (!isDigitString (sgb))
						continue;				// Skip invalid line
					if (rising_ns && !rising_ks && (n <= last_processed_n))
						continue;				// Skip already processed n's
					if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
						continue;				// Skip already processed k's
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					if (! process_num (format, sgk, sgb, n, incr, shift, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%s %s %lu %d\n", sgk, sgb, n, incr); 
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCRU)	{	// Repunits, n is the only parameter.
					if (sscanf (buff+begline, "%lu", &n) != 1)
						continue;				// Skip invalid line
					sprintf (sgb, "10");
					sprintf (sgk, "1");
					if (rising_ns && !rising_ks && (n <= last_processed_n))
						continue;				// Skip already processed n's
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					if (! process_num (format, (char*)"1", (char*)"10", n, -1, 0, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%lu\n", n); 
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCGRU)	{	// Generalized Repunits, b, n, are the two parameters
					if (sscanf (buff+begline, "%s %lu", sgb, &n) != 2)
						continue;				// Skip invalid line
					if (!isDigitString (sgb))
						continue;				// Skip invalid line
					sprintf (sgk, "1");
					if (rising_ns && !rising_ks && (n <= last_processed_n))
						continue;				// Skip already processed n's
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					if (! process_num (format, (char*)"1", sgb, n, -1, 0, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%s %lu\n", sgb, n); 
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCGF)	{	// Generalized Fermat, sgb, n, are the two parameters
					if (sscanf (buff+begline, "%s %lu", sgb, &n) != 2)
						continue;				// Skip invalid line
					if (!isDigitString(sgb))
						continue;				// Skip invalid line
					if (!ispoweroftwo(n))
						continue;				// Skip invalid line
					sprintf (sgk, "1");
					if (! process_num (format, (char*)"1", sgb, n, 1, 0, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
                                                        _write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%s %lu\n", sgb, n);	// write the result
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCDN)	{	// b^n-b^m+c numbers ; sgb, n, m are the three parameters
					if (sscanf (buff+begline, "%s %lu %lu", sgb, &n, &m) != 3)
						continue;				// Skip invalid line
					if (!isDigitString(sgb))
						continue;				// Skip invalid line
					if (n <= m)
						continue;				// Skip invalid line
					ndiff = n-m;				// Save difference of exponents in a global
					sprintf (sgk, "1");
					if (! process_num (format, (char *)"1", sgb, m, incr, 0, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char *)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%s %lu %lu\n", sgb, n, m);	// write the result
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char  *)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCDNG)	{	// b^n-b^m+c numbers ; sgb, n, m, c are the four parameters
					if (sscanf (buff+begline, "%s %lu %lu %d", sgb, &n, &m, &incr) != 4)
						continue;				// Skip invalid line
					if (!isDigitString(sgb))
						continue;				// Skip invalid line
					if (n <= m)
						continue;				// Skip invalid line
					ndiff = n-m;				// Save difference of exponents in a global
					sprintf (sgk, "1");
					if (! process_num (format, (char *)"1", sgb, m, incr, 0, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char *)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%s %lu %lu %d\n", sgb, n, m, incr);	// write the result
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char *)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCVARAQS)	{	// k, b, n, c and d specified on each input line
					if (sscanf (buff+begline, "%s %s %lu %d %s", sgk, sgb, &n, &incr, sgd) != 5)
						continue;				// Skip invalid line
					if (!isDigitString(sgk))
						continue;				// Skip invalid line
					if (!isDigitString (sgb))
						continue;				// Skip invalid line
					if (!isDigitString(sgd))
						continue;				// Skip invalid line
					if (rising_ns && !rising_ks && (n <= last_processed_n))
						continue;				// Skip already processed n's
					if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
						continue;				// Skip already processed k's
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					if (! process_num (format, sgk, sgb, n, incr, shift, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
                                                        _write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%s %s %lu %d %s\n", sgk, sgb, n, incr, sgd); 
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCGM)	{	// Gaussian Mersenne
					if ((nargs = sscanf (buff+begline, "%lu %lu %lu", &n, &facn, &facnp)) < 1)
						continue;				// Skip invalid line
					else if (nargs == 1)		// Not prefactored.
						facn = facnp = 0;
					else if (nargs == 2) {		// Second argument is how far already factored, in bits)
						if (!facfrom)
							facfrom = facn;
						facn = facnp = 0;
					}
					if (rising_ns && !rising_ks && (n <= last_processed_n))
						continue;				// Skip already processed n's
					sprintf (sgk, "2^%lu", (n+1)/2);
					sprintf (sgb, "2");
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					if (! process_num (format, sgk, (char*)"2", n, +1, shift, &res))
						goto done;
#ifndef X86_64
					if (facto) {				// If factoring, print a job progress message every so often
						if (n/pdivisor-pquotient == 1) {
							sprintf (outbuf, "%lu candidates factored, %lu factors found, %lu remaining\n"
								, factored, eliminated, factored - eliminated);
							OutputBoth (outbuf);
							pquotient = n/pdivisor;
						}
						else if (n/pdivisor-pquotient > 1)
							pquotient = n/pdivisor;
					}
#endif
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						sign = (((n&7) == 3) || ((n&7) == 5))? 1 : 0;	// 1 if positive, 0 if negative
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
#ifndef X86_64
							if (facto)
								if (n >= LOWFACTORLIMIT)
									sprintf (outbuf, "%lu %lu\n", n, facto); 
								else if (facfrom)
									sprintf (outbuf, "%lu %lu\n", n, facfrom); 
								else
									sprintf (outbuf, "%lu\n", n); 
							else if (res1 && res2)
#else
							if (res1 && res2)
#endif
								if (a)
									sprintf (outbuf, "%lu (GM(%lu) is Prime in Z+iZ and the norm of GQ(%lu) is %ld-PRP.)\n", n, n, n, a); 
								else
									sprintf (outbuf, "%lu (GM(%lu) and GQ(%lu) are Prime in Z+iZ.)\n", n, n, n); 
							else if (res1)
								sprintf (outbuf, "%lu (GM(%lu) is Prime in Z+iZ.)\n", n, n); 
							else
								if (a)
									sprintf (outbuf, "%lu (The norm of GQ(%lu) is %ld-PRP.)\n", n, n, a); 
								else
									sprintf (outbuf, "%lu (GQ(%lu) is Prime in Z+iZ.)\n", n, n); 
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						if (res2) {
							if (sign) {
								outfdm = _open (outmf, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
								if (outfdm) {
									sprintf (outbuf, "%lu\n", n); 
									_write (outfdm, outbuf, strlen (outbuf));
									_close (outfdm);
								}
							}
							else {
								outfdp = _open (outpf, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
								if (outfdp) {
									sprintf (outbuf, "%lu\n", n); 
									_write (outfdp, outbuf, strlen (outbuf));
									_close (outfdp);
								}
							}
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCSP)	{	// SPRP test of (2^n+1)/3 numbers
					if ((nargs = sscanf (buff+begline, "%lu %lu", &n, &facn)) < 1)
						continue;				// Skip invalid line
					else if (nargs == 1)		// Not prefactored.
						facn = facnp = 0;
					else if (nargs == 2) {		// Second argument is how far already factored, in bits)
						if (!facfrom)
							facfrom = facn;
						facn = facnp = 0;
					}
					if (rising_ns && !rising_ks  && (n <= last_processed_n))
						continue;				// Skip already processed n's
					sprintf (sgk, "(2^%lu+1)/3", n);
					sprintf (sgb, "2");
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					if (! process_num (format, sgk, (char*)"2", n, +1, shift, &res))
						goto done;
#ifndef X86_64
					if (facto) {				// If factoring, print a job progress message every so often
						if (n/pdivisor-pquotient == 1) {
								sprintf (outbuf, "%lu candidates factored, %lu factors found, %lu remaining\n"
									, factored, eliminated, factored - eliminated);
								OutputBoth (outbuf);
							pquotient = n/pdivisor;
						}
						else if (n/pdivisor-pquotient > 1)
							pquotient = n/pdivisor;
					}
#endif
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
#ifndef X86_64
						if (facto)
							if (n >= LOWFACTORLIMIT)
								sprintf (outbuf, "%lu %lu\n", n, facto); 
							else if (facfrom)
								sprintf (outbuf, "%lu %lu\n", n, facfrom); 
							else
								sprintf (outbuf, "%lu\n", n); 
						else
							sprintf (outbuf, "%lu\n", n); 
#else
						sprintf (outbuf, "%lu\n", n); 
#endif
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
				else if (format == ABCK) {							// Carol/Kynea
					if (sscanf (buff+begline, "%lu %d", &n, &incr) != 2)
						continue;						// Skip invalid line
					if (rising_ns && !rising_ks  && (n <= last_processed_n))
						continue;				// Skip already processed n's
					if (incr == 1) {
						format = ABCK;					// Kynea number
						sprintf (sgk, "(2^%lu+1)", n-1);
					}
					else if (incr == -1) {
						format = ABCC;					// Carol number
						sprintf (sgk, "(2^%lu-1)", n-1);
					}
					else
						continue;
					sprintf (sgb, "2");
					if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
						fclose (fd);			// Unlock the file during the test...
					if (! process_num (format, sgk, (char*)"2", n+1, -1, shift, &res))
						goto done;
					if (res) {
						resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
						outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
						if (outfd) {
							if (hline >= resultline) {	// write the relevant header
								_write (outfd, hbuff, strlen (hbuff));
							}
							sprintf (outbuf, "%lu %d\n", n, incr); 
							_write (outfd, outbuf, strlen (outbuf));
							_close (outfd);
						}
						IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
					}
				}
			}				// End processing a data line

			if ((!rising_ns && !rising_ks) || (rising_ns && rising_ks))
				IniWriteInt (INI_FILE, (char*)"PgenLine", line + 1);		// Point on the next line
			if (rising_ns && !rising_ks) {
				IniWriteInt (INI_FILE, (char*)"Last_Processed_n", n);		// Point on the next n
				last_processed_n = n;
			}
			if (rising_ks && !rising_ns) {
				IniWriteString (INI_FILE, (char*)"Last_Processed_k", sgk); // Point on the next k
				strcpy (last_processed_k, sgk);
			}
			if(n>=(unsigned long)IniGetInt(INI_FILE, (char*)"MaxN", 2147483647)) {
				break;
			}
			if (res) {
				if(IniGetInt(INI_FILE, (char*)"BeepOnSuccess", 0)) {
					do {	// If stopping on this success, beep infinitely!
#if !defined(WIN32) || defined(_CONSOLE)
						//cuda flashWindowAndBeep (20);
#else
						flashWindowAndBeep (50);
#endif
					} while (!stopCheck () && IniGetInt(INI_FILE, (char*)"StopOnSuccess", 0));
				}
				if(IniGetInt(INI_FILE, (char*)"StopOnSuccess", 0)) {
					goto done;
				}
				else if(IniGetInt(INI_FILE, (char*)"StopOnPrimedK", 0)) {
					sprintf (outbuf, "ks%s", sgk);
					IniWriteInt (INI_FILE, outbuf, 1+IniGetInt(INI_FILE, outbuf, 0));
							// Increment this k success count
					save_IniFile (INI_FILE, SVINI_FILE);	// make a backup of INI_FILE
				}
				else if(IniGetInt(INI_FILE, (char*)"StopOnPrimedN", 0)) {
					sprintf (outbuf, "ns%lu", n);
					IniWriteInt (INI_FILE, outbuf, 1+IniGetInt(INI_FILE, outbuf, 0));
							// Increment this n success count
					save_IniFile (INI_FILE, SVINI_FILE);	// make a backup of INI_FILE
				}
				else if(IniGetInt(INI_FILE, (char*)"StopOnPrimedB", 0)) {
					sprintf (outbuf, "bs%s", sgb);
					IniWriteInt (INI_FILE, outbuf, 1+IniGetInt(INI_FILE, outbuf, 0));
							// Increment this base success count
					save_IniFile (INI_FILE, SVINI_FILE);	// make a backup of INI_FILE
				}
			}
			if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
				goto OPENFILE;
		}					// End of loop on input lines
		IniWriteString (INI_FILE, (char*)"ResultLine", NULL);		// delete the result line
		_unlink (SVINI_FILE);								// delete the backup of INI_FILE
		completed = TRUE;
done:
		if(IniGetInt(INI_FILE, (char*)"StopOnSuccess", 0) && res) {
			if ((!rising_ns && !rising_ks) || (rising_ns && rising_ks))
				IniWriteInt (INI_FILE, (char*)"PgenLine", line + 1);		// Point on the next line
			if (rising_ns && !rising_ks)
				IniWriteInt (INI_FILE, (char*)"Last_Processed_n", n);		// Point on the next n
			if (rising_ks && !rising_ns)
				IniWriteString (INI_FILE, (char*)"Last_Processed_k", sgk); // Point on the next k
		}
		else if (!aborted && ((!rising_ns && !rising_ks) || (rising_ns && rising_ks)))
			IniWriteInt (INI_FILE, (char*)"PgenLine", line);		// Point again on the current line...
                IniWriteString (INI_FILE, (char*)"MaxRoundOff", NULL);
		if (facto) {
			sprintf (outbuf, "%lu candidates factored, %lu factors found, %lu remaining\n"
				, factored, eliminated, factored - eliminated);
			OutputBoth (outbuf);
		}
		if ((!rising_ns && !rising_ks) || (rising_ns && rising_ks))
			fclose (fd);
		IniWriteString(INI_FILE, (char*)"OldInputFile", inputfile);		// Save the just processed input file name.
	}						// End Work == 0

// Handle an expr

	else {					// Work != 0
		OutputStr ((char*)"Expression testing not yet implemented.\n");
		IniWriteInt (INI_FILE, (char*)"Workdone", 1);
	}
	aborted = FALSE;
	return (completed);
}
