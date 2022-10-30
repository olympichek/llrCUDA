/*----------------------------------------------------------------------
| This file contains routines and global variables that are common for
| all operating systems the program has been ported to.  It is included
| in one of the source code files of each port.  See Llr.h for the
| common #defines and common routine definitions.
+---------------------------------------------------------------------*/
 
#define CHECK_IF_ANY_ERROR(X,J,N,K) \
		checknumber = K;\
		restarting = FALSE;\
		will_try_larger_fft = FALSE;\
\
/* Check for excessive roundoff error  */\
\
		if (MAXERR > maxroundoff) {\
			lasterr_point = J;\
			if (J == last_bit[K] &&\
			    MAXERR == last_maxerr[K] && !abonroundoff && !care/*will_try_larger_fft*/) {\
				clearline(100);\
				OutputBoth (ERROK);\
				GWERROR = 0;\
				clearline(100);\
				IniWriteInt(INI_FILE, (char*)"Error_Count", error_count = IniGetInt(INI_FILE, (char*)"Error_Count", 0) + 1);\
				if (error_count > MAX_ERROR_COUNT) {\
					OutputBoth (ERRMSG9);\
					IniWriteString(INI_FILE, (char*)"Error_Count", NULL);\
					will_try_larger_fft = TRUE;\
					last_bit[K]  = 0;\
					last_maxerr[K]  = 0.0;\
					maxerr_recovery_mode[K] = FALSE;\
				}\
				else {\
                                    OutputBoth (ERRMSG6);\
                                    maxerr_recovery_mode[K] = TRUE;\
				}\
				sleep5 = FALSE;\
				restarting = TRUE;\
				goto error;\
			} else {\
				char	msg[80];\
				sprintf (msg, ERRMSG1C, MAXERR, maxroundoff);\
				sprintf (buf, ERRMSG0L, J, N, msg);\
				clearline(100);\
				OutputBoth (buf);\
				if (J == last_bit[K])\
					will_try_larger_fft = TRUE;\
				if (will_try_larger_fft) {\
					OutputBoth (ERRMSG8);\
					IniWriteString(INI_FILE, (char*)"Error_Count", NULL);\
					last_bit[K]  = 0;\
					last_maxerr[K]  = 0.0;\
				}\
				else {\
					last_bit[K] = J;\
					last_maxerr[K] = MAXERR;\
				}\
				maxerr_recovery_mode[K] = FALSE;\
				sleep5 = FALSE;\
				restarting = TRUE;\
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

extern void print (double*, int);   // debugging tool...

int it = 0;         //cuda ; init to zero, JP 22/06/17
ssize_t wc = 0;     // default return value for write function. JP 29/11/20

unsigned long format;

// Some ABC format strings

char ckstring[] = "(2^$a$b)^2-2";	// Carol/Kynea
char cwstring[] = "$a*$b^$a$c";		// Cullen/Woodall
char ffstring[] = "$a*2^$b+1";		// FermFact output
char ffmstring[] = "$a*2^$b-1";		// Lei output
char gmstring[] = "4^$a+1";		// Gaussian Mersenne norms

char gfstring[] = "$a^$b+1";		// Special primality test for generalized Fermat numbers
char spstring[] = "(2^$a+1)/3";		// Special SPRP test for Wagstaff numbers
char dpstring[] = "DivPhi($a*$b^$c+1)"; // like-PRP test for DivPhi(...,2)
char repustring[] = "(10^$a-1)/9";	// PRP test for repunits numbers
char grepustring[] = "($a^$b-1)/($a-1)";// PRP test for generalized repunits numbers
char diffnumpstring[] = "$a^$b-$a^$c+%d";// If $b>$c, it is [$a^($b-$c)-1]*$a^$c+%d so, form K*B^N+C
char diffnummstring[] = "$a^$b-$a^$c-%d";// If $b>$c, it is [$a^($b-$c)-1]*$a^$c-%d so, form K*B^N+C
char diffnumstring[] = "$a^$b-$a^$c$d";	// General diffnum format
// Fixed k and c forms for k*b^n+c

char fkpstring[] = ""$LLF"*$a^$b+%d";
char fkmstring[] = ""$LLF"*$a^$b-%d";

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

// k*b^n+c format with k, b, c fixed

char abcnstring[] = "%lu*%lu^n%d";

// General (k*b^n+c)/d format 

char abcadstring[] = "($a*$b^$c$d)/$e";

// (k*b^n+c)/d format with k, b, c, d fixed

char abcndstring[] = "(%lu*%lu^n%d)/%lu";

// Test the primality of a number given as a string

char numberstring[] = "$a";

// Test if $a is a base $b Wieferich prime

char wftstring[] = "$a$b";

// Search for base $c Wieferich primes in the range $a to $b

char wfsstring[] = "$a$b$c";

#define ABCVARAQS   19 // k, b, n, c and divisor specified on each input line or header
#define ABCVARAS    18 // k, b, n, and c specified on each input line or header
#define ABCDP	    28  // Format used for DivPhi()


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

#define _log2(x) log(x)/log(2)

char	greatbuf[10001] = {0};
char	INI_FILE[80] = {0};
char	SVINI_FILE[80] = {0};
char	RESFILE[80] = {0};
char	LOGFILE[80] = {0};
char	EXTENSION[8] = {0};

// int fftlen = 0;
int ERRCHK = 0;
unsigned int PRIORITY = 1;
unsigned int CPU_AFFINITY = 99;
unsigned int CPU_MASK = 0;
unsigned int CPU_TYPE = 0;
unsigned long volatile ITER_OUTPUT = 0;
unsigned long volatile ITER_OUTPUT_RES = 99999999;
unsigned long volatile DISK_WRITE_TIME = 30;
unsigned long INTERIM_FILES = 0;
unsigned long INTERIM_RESIDUES = 0;
int	CLASSIC_OUTPUT = 0;
int	OUTPUT_ROUNDOFF = 0;
int	CUMULATIVE_ROUNDOFF = 1;
int	NUM_BACKUP_FILES = 3;   // was 3
int	NUM_JACOBI_BACKUP_FILES = 2;    // was 2
int	RUN_ON_BATTERY = 1;
int	TRAY_ICON = TRUE;
int	HIDE_ICON = FALSE;
int     mul_final = 0;
double UNOFFICIAL_CPU_SPEED = 0.0;
unsigned int WORKTODO_COUNT = 0;/* Count of valid work lines */
unsigned int ROLLING_AVERAGE = 0;
unsigned int PRECISION = 2;
unsigned long NUM_CPUS = 1;	/* Number of CPUs/Cores in the computer */
unsigned long NUM_WORKER_THREADS = 1;
int	TWO_BACKUP_FILES = 1;
int	HIGH_RES_TIMER = 0; 
int	usingDivPhi_m = 0;
int	string_rep_truncated = FALSE;
int     recovering = FALSE;
unsigned long CPU_SPEED = 25;
int ZERO_PADDED_FFT;
int GWERROR = 0;

double maxdiffmult = 1.0;

/* PRP and LLR global variables */

#define	sgkbufsize 20000

giant	N = NULL;	/* Number being tested */
giant	NP = NULL;	/* Number being tested */
giant	M = NULL;	/* Gaussian Mersenne modulus = N*NP */
giant	gk = NULL;	/* k multiplier */
giant	gb = NULL;	/* Generalized Fermat base may be a large integer... */

unsigned long Nlen = 0;	/* Bit length of number being LLRed or PRPed */
unsigned long klen = 0;	/* Number of bits of k multiplier */
unsigned long n_orig;	/* exponent associated to initial base sgb */
long OLDFFTLEN = 0; /* previous value of FFTLEN, used by setuponly option */
unsigned long ndiff = 0;/* used for b^n-b^m+c number processing */
unsigned long gformat;	/* used for b^n-b^m+c number processing */
							
struct work_unit *w = NULL;	// Work to do data as described in Prime95

/* Other gwypnum globals */

giant testn, testnp;
unsigned long facn = 0, facnp = 0;
int resn = 0, resnp = 0;
char facnstr[80], facnpstr[80];
// char m_pgen_input[IBSIZE], m_pgen_output[IBSIZE], oldm_pgen_input[IBSIZE];
// char keywords[10][IBSIZE], values[10][IBSIZE];
 char multiplier[IBSIZE], base[IBSIZE], exponent[IBSIZE], exponent2[IBSIZE], addin[IBSIZE], divisors[IBSIZE];
// char inifilebuf[IBSIZE];
char    string_rep[80];
char    sgd[sgkbufsize];
char    sgb[sgkbufsize];
char    bpfstring[sgkbufsize];
char    cMAXBPD[10]={0};
/*
#ifndef X86_64

void setupf();
int factor64();
void psetupf();
int pfactor64();

#endif

void* aligned_malloc (unsigned long, unsigned long);
void  aligned_free (void *);
*/
static unsigned long last_bit[10] = {0};
static double last_suminp[10] = {0.0};
static double last_sumout[10] = {0.0};
static double last_maxerr[10] = {0.0};
static double maxroundoff = 0.40;
double MAXBPD = 37.0;
// static double fpart = 0.0;

static unsigned long mask;

extern int zcomplex;    // CUDA
extern int cufftonly;   // CUDA

unsigned int fermat_only = FALSE;

unsigned int strong = TRUE;
unsigned int quotient = FALSE;
unsigned int vrbareix = FALSE;
unsigned int dualtest = FALSE;
unsigned int bpsw = FALSE;
unsigned int setuponly = FALSE;
unsigned int nosaving = FALSE;
unsigned int abonillsum = FALSE;
unsigned int abonmismatch = FALSE;
unsigned int testgm = FALSE;
unsigned int testgq = FALSE;
unsigned int testfac = FALSE;
unsigned int nofac = FALSE;
unsigned int general = FALSE;
unsigned int eps2 = FALSE;
unsigned int vdebug = FALSE;
unsigned int restarting = FALSE;
// Actually set by Roundoff error condition ; reset by the macro.
unsigned int care = FALSE;  // Set after a careful iteration ; reset after a normal one.
unsigned int nbfftinc = 0;  // Number of required FFT increments.
unsigned int maxfftinc = 3; // Maximum accepted FFT increments.
unsigned int aborted = FALSE;
unsigned int abonroundoff = FALSE;
unsigned int will_try_larger_fft = FALSE;
unsigned int checknumber = 0;
unsigned int error_count = 0;
// Number of reproducible errors that previously occured.
unsigned int sleep5 = FALSE, showdigits = FALSE;
unsigned int maxerr_recovery_mode [10] = {0};
unsigned int lasterr_point = 0;
unsigned long primolimit = 30000;
unsigned long nextifnear = FALSE;
unsigned long maxaprcl = 400;
unsigned long interimFiles, interimResidues, throttle, facfrom, facto;
unsigned long factored = 0, eliminated = 0;
unsigned long pdivisor = 1000000, pquotient = 1;
unsigned long bpf[30], bpc[30], vpf[30];
// Base prime factors, cofactors, power of p.f.
unsigned long *t = NULL,*ta = NULL;// Precrible arrays, dynamically allocated.
unsigned long smallprime[168] =	// Small primes < 1000
{2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,
179,    181,    191,    193,    197,    199,    211,    223,    227,    229,
233,    239,    241,    251,    257,    263,    269,    271,    277,    281,
283,    293,    307,    311,    313,    317,    331,    337,    347,    349,
353,    359,    367,    373,    379,    383,    389,    397,    401,    409,
419,    421,    431,    433,    439,    443,    449,    457,    461,    463,
467,    479,    487,    491,    499,    503,    509,    521,    523,    541,
547,    557,    563,    569,    571,    577,    587,    593,    599,    601,
607,    613,    617,    619,    631,    641,    643,    647,    653,    659,
661,    673,    677,    683,    691,    701,    709,    719,    727,    733,
739,    743,    751,    757,    761,    769,    773,    787,    797,    809,
811,    821,    823,    827,    829,    839,    853,    857,    859,    863,
877,    881,    883,    887,    907,    911,    919,    929,    937,    941,
947,    953,    967,    971,    977,    983,    991,    997};
unsigned long smallphi[10] = {1,2,4,6,10,12,16,18,22,28};
giantstruct*	gbpf[30] = {NULL}; // Large integer prime factors
giantstruct*	gbpc[30] = {NULL}; // Large integer cofactors
unsigned long nrestarts = 0;
// Nb. of restarts for an N+1 or N-1 prime test
int      nbdg = 0, nbdg1, nbdg2;
unsigned long globalb = 2;	// base of the candidate in a global
double	globalk = 1.0;   // k value of the candidate in a global

double	pcfftlim = 0.5;	
int readlg, writelg;	// Values returned by low level I/O

// double timers[10] = {0.0};	/* Up to five separate timers */

double smargin = 0.0;

// int is_valid_double(double);
int genProthBase(giant, uint32_t);
long generalLucasBase(giant , uint32_t *, uint32_t *);
unsigned long gcd (
    unsigned long x,
    unsigned long y);

// int ispower (unsigned long,
//	unsigned long);

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

#define MAX_ERROR_COUNT 5

/* Routines missing from GMP */

#define mpz_add_si(d,s,addin)	if (addin >= 0) mpz_add_ui(d,s,(unsigned int)addin); else mpz_sub_ui(d,s,(unsigned int)(-addin));
#define mpz_mul_d(d,s,flt)	{ mpz_t t; mpz_init_set_d(t,flt); mpz_mul(d,s,t); mpz_clear(t); }
#define mpz_eq(a,b)		mpz_cmp(a,b) == 0


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
char ERRMSG9[] = "Too much errors ; Restarting with next larger FFT length...\n";
char ERRMSG60[] = "ERROR: Comparing double-check values failed.  Rolling back to iteration %lu.\n";
char ERRMSG70[] = "ERROR: Comparing Gerbicz checksum values failed.  Rolling back to iteration %lu.\n";
char ERRMSG80[] = "ERROR: Invalid FFT data.  Restarting from last save file.\n";
char ERRMSG90[] = "ERROR: Invalid PRP state.  Restarting from last save file.\n";
char WRITEFILEERR[] = "Error writing intermediate file: %s\n";

void	trace(int n) {			// Debugging tool...
    char buf[100];
    
    sprintf(buf, "OK until number %d\n", n);
    if (verbose)
        OutputBoth (buf);
    else
        OutputStr (buf);
}

//  giant gcd and invert functions implementations using GNU MP library functions

unsigned long gwypgcdui (giant g, unsigned long x) {
    mpz_t rop, op1;
    unsigned long result;
    mpz_init (rop);             // allocate mpz memory
    mpz_init (op1);
    gtompz (g, op1);            // g ==> op1
    result = mpz_gcd_ui (rop, op1, x);  // returns zero if x is zero, but then, gcd is g...
    mpz_clear (rop);            // free mpz memory
    mpz_clear (op1);
    return (result);
}

void gwypgcdg (giant g1, giant g2) {
    mpz_t rop, op1, op2;
    mpz_init (rop);             // allocate mpz memory
    mpz_init (op1);
    mpz_init (op2);
    gtompz (g1, op1);           // g1 ==> op1
    gtompz (g2, op2);           // g2 ==> op2
    mpz_gcd (rop, op1, op2);    // gcd (op1, op2) ==> rop
    mpztog (rop, g2);           // rop ==> g2
    mpz_clear (rop);            // free mpz memory
    mpz_clear (op1);
    mpz_clear (op2);
}

int gwypinvg (giant g1, giant g2) {
    int result;
    mpz_t rop, op1, op2;
    mpz_init (rop);             // allocate mpz memory
    mpz_init (op1);
    mpz_init (op2);
    gtompz (g1, op1);           // g1 ==> op1
    gtompz (g2, op2);           // g2 ==> op2
    result = mpz_invert (rop, op2, op1);    // invert of op2 mod op1 ==> rop
    if (result == 0) {        
        mpz_gcd (rop, op1, op2);
        OutputBoth ((char*)"inverse does not exist, so, gcd is returned in second operand...\n");
    }
    mpztog (rop, g2);           // rop ==> g2
    mpz_clear (rop);            // free mpz memory
    mpz_clear (op1);
    mpz_clear (op2);
    return (result);
}

//*******************************************************************************
    
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

void clearbuf (char *b) {
    int i;
    for (i=0;i<strlen (b);i++)
        b[i] = ' ';
    b[i] = '\0';
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
    BlinkIcon (10);	    /* Blink icon for 10 seconds */
    Sleep (10000);
    ChangeIcon (IDLE_ICON); /* Idle icon for rest of 5 minutes */
    for (i = 0; i < 290; i++) {
        Sleep (1000);
        if (escapeCheck ()) return (FALSE);
    }
    ChangeIcon (WORKING_ICON);
                            /* And back to the working icon */
    return (TRUE);
}

/* Generate the scaling factors for ITER_OUTPUT in the rare cases where the user */
/* has used some undoc.txt settings to change how often the title is output or to */
/* make the frequency roughly the same in all windows even if using different FFT sizes. */

void calc_output_frequencies (
	double	*output_frequency,	/* Calculated adjustment to ITER_OUTPUT */
	double	*output_title_frequency)/* Calculated adjustment to ITER_OUTPUT for title */
{
	int	title_freq;
//	double	exp, temp;

        *output_frequency = 1.0;
	/* Calculate the title frequency as a fraction of the output frequency */
	title_freq = (int) IniGetInt (INI_FILE, (char*)"TitleOutputFrequency", 1);
	if (title_freq < 1) title_freq = 1;
	*output_title_frequency = *output_frequency / (double) title_freq;
}

/* Truncate a percentage to the requested number of digits. */
/* Truncating prevents 99.5% from showing up as 100% complete. */

double trunc_percent (
	double	percent)
{
    if (percent > 100.0) percent = 100.0;
    percent -= 0.5 * pow (10.0, - (double) PRECISION);
    if (percent < 0.0)
        return (0.0);
    return (percent);
}

 
//  Test if a string contains only valid digits. 
 
int isDigitString(char *s) { 
    while (*s) { 
	if (!isdigit(*s)) 
            return (FALSE); 
	s++; 
    } 
    return (TRUE); 
} 

// The eight timer routines are now in the gwypnum.cu  file

void OutputTimeStamp ()
{
    time_t	this_time;
    char	tbuf[40], buf[40];

    time (&this_time);
    strcpy (tbuf, ctime (&this_time)+4);
    tbuf[12] = 0;
    sprintf (buf, "[%s] ", tbuf);
    OutputStr (buf);
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
        if (p % i == 0) 
            return (FALSE);
    return (TRUE);
}

/* Determine the names of the INI files */

void nameIniFiles (
    int	named_ini_files)
{
    char   buf[120];

    if (named_ini_files < 0) {
        strcpy (INI_FILE, "llr.ini");
        strcpy (RESFILE, "lresults.txt");
        strcpy (LOGFILE, "lprime.log");
        strcpy (EXTENSION, "");
    } 
    else {
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
        if (_chdir (buf))
            printf ("Could not change the working directory...\n");
;
        IniFileOpen (INI_FILE, 0);
    }
}

/* Read the INI files */

void readIniFiles () 
{ 
    int	temp; 
    const char *p;
    char cpustring[100];
 
//  getCpuInfo (); 
 
    PRECISION = (unsigned int) IniGetInt (INI_FILE, (char*)"PercentPrecision", 2); 
    if (PRECISION > 6) 
        PRECISION = 6; 
 
    ITER_OUTPUT = IniGetInt (INI_FILE, (char*)"OutputIterations", 10000); 
    if (ITER_OUTPUT <= 0) 
        ITER_OUTPUT = 1; 
    ITER_OUTPUT_RES = IniGetInt (INI_FILE,
        (char*)"ResultsFileIterations", 99999999); 
    if (ITER_OUTPUT_RES < 1000) 
        ITER_OUTPUT_RES = 1000; 
    DISK_WRITE_TIME = IniGetInt (INI_FILE, (char*)"DiskWriteTime", 30); 
    TWO_BACKUP_FILES = (int) IniGetInt (INI_FILE,
        (char*)"TwoBackupFiles", 1); 
    RUN_ON_BATTERY = (int) IniGetInt (INI_FILE, (char*)"RunOnBattery",
        1); 
 
    temp = (int) IniGetInt (INI_FILE, (char*)"ErrorCheck", 0); 
    ERRCHK = (temp != 0); 
    PRIORITY = (unsigned int) IniGetInt (INI_FILE, (char*)"Priority", 
        1); 
    CPU_AFFINITY = (unsigned int) IniGetInt (INI_FILE,
        (char*)"Affinity", 99); 
        if (CPU_AFFINITY!=99) {
            CPU_MASK = 1<<CPU_AFFINITY;
//            p = IniSectionGetNthStringRaw (INI_FILE, NULL, "Affinity", 1);
            IniGetString (INI_FILE, (char*)"Affinity", cpustring, 100, (char*)"99");
            p = strchr (cpustring, ',');
            while (p != NULL) {
                p++;
                sscanf (p, "%d", &CPU_AFFINITY);
                CPU_MASK |= 1<<CPU_AFFINITY;
                p = strchr (p, ',');
            }
            printf ("CPU_MASK = %d\n", CPU_MASK);
        }
    HIDE_ICON = (int) IniGetInt (INI_FILE, (char*)"HideIcon", 0); 
    TRAY_ICON = (int) IniGetInt (INI_FILE, (char*)"TrayIcon", 1); 
 
/* Guess the CPU type if it isn't known.  Otherwise, validate it. */ 
 
//  getCpuInfo (); 
 
/* Other oddball options */ 
 
    CUMULATIVE_TIMING = IniGetInt (INI_FILE, (char*)"CumulativeTiming",
        0); 
//  HIGH_RES_TIMER = isHighResTimerAvailable (); 
} 
 
/*----------------------------------------------------------------------
| Portable routines to read and write ini files!  NOTE:  These only
| work if you open no more than 5 ini files.  Also you must not
| change the working directory at any time during program execution.
+---------------------------------------------------------------------*/

struct IniLine {
    char    *keyword;
    char    *value;
    int	    active;
};

struct IniCache {
    char    *filename;
    int	    immediate_writes;
    int	    dirty;
    unsigned int num_lines;
    unsigned int array_size;
    struct   IniLine **lines;
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
    FILE    *fd;
    unsigned int i;
    char    line[80];
    char    *val;

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
    if (fd == NULL) 
        return (p);

    while (fgets (line, 80, fd)) {
        if (line[strlen(line)-1] == '\n') 
            line[strlen(line)-1] = 0;
        if (line[0] == 0) 
            continue;
        if (line[strlen(line)-1] == '\r')
            line[strlen(line)-1] = 0;
        if (line[0] == 0) 
            continue;
        val = strchr (line, '=');
        if (val == NULL) {
            char    buf[130];
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
    char    buf[100];

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
        wc = _write (fd, buf, strlen (buf));
    }
    p->dirty = 0;
    _close (fd);
}

void save_IniFile (char *filename, char *savedfilename) {
    struct IniCache *p;	
    p = openIniFile (filename, 1);  // open and read the source IniFile.
    p->filename = savedfilename;
                            // put the target filename in the structure.
    writeIniFile (p);	    // Write the target.
    p->filename = filename; // Restore the structure in memory.
}

void truncated_strcpy (
    char    *buf,
    unsigned int bufsize,
    char    *val)
{
    if (strlen (val) >= bufsize) {
        memcpy (buf, val, bufsize-1);
        buf[bufsize-1] = 0;
    } 
    else {
        strcpy (buf, val);
    }
}

void IniGetString (
    char    *filename,
    char    *keyword,
    char    *val,
    unsigned int val_bufsize,
    char    *default_val)
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
            } 
            else {
                truncated_strcpy (val, val_bufsize, default_val);
            }
            return;
        }
        if (p->lines[i]->active &&
            stricmp (keyword, p->lines[i]->keyword) == 0) 
            break;
    }

/* Copy info from the line structure to the user buffers */

    truncated_strcpy (val, val_bufsize, p->lines[i]->value);
}

long IniGetInt (
    char    *filename,
    char    *keyword,
    long    default_val)
{
    char    buf[20], defval[20];
    
    sprintf (defval, "%ld", default_val);
    IniGetString (filename, keyword, buf, 20, defval);
    return (atol (buf));
}

void IniWriteString (
    char    *filename,
    char    *keyword,
    char    *val)
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
            if (val != NULL && strcmp (val, p->lines[i]->value) == 0)   
                return;
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
    char    *filename,
    char    *keyword,
    long    val)
{
    char    buf[40];
    
    sprintf (buf, "%ld", val);
    IniWriteString (filename, keyword, buf);
}

void IniFileOpen (
    char    *filename,
    int	immediate_writes)
{
    struct IniCache *p;
    
    p = openIniFile (filename, 1);
    p->immediate_writes = immediate_writes;
}

void IniFileClose (
    char    *filename)
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
    char    *filename)
{
    struct IniCache *p;
    int	fd;
    unsigned int j;
    char    buf[100];

/* Create and write out the INI file */

    p = openIniFile (filename, 0);
    fd = _open (p->filename, _O_CREAT | _O_TRUNC | _O_WRONLY | _O_TEXT, 0666);
    if (fd < 0) 
        return (FALSE);
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
    char    *filename)
{
    struct IniCache *p;

    p = openIniFile (filename, 0);
    return (p->num_lines);
}

void IniGetLineAsString (
    char    *filename,
    unsigned int line,
    char    *keyword,
    unsigned int keyword_bufsize,
    char    *val,
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
    char    *filename,
    unsigned int line,
    char    *keyword,
    unsigned int keyword_bufsize,
    long    *val)
{
    char    buf[20];
    
    IniGetLineAsString (filename, line, keyword, keyword_bufsize, buf, 20);
    *val = atol (buf);
}

void IniReplaceLineAsString (
    char    *filename,
    unsigned int line,
    char    *keyword,
    char    *val)
{
    IniDeleteLine (filename, line);
    IniInsertLineAsString (filename, line, keyword, val);
}

void IniReplaceLineAsInt (
    char    *filename,
    unsigned int line,
    char    *keyword,
    long    val)
{
    char    buf[20];
    
    sprintf (buf, "%ld", val);
    IniReplaceLineAsString (filename, line, keyword, buf);
}

void IniInsertLineAsString (
    char    *filename,
    unsigned int line,
    char    *keyword,
    char    *val)
{
    struct IniCache *p;
    unsigned int i;

/* Open ini file, do not reread it as that could change the line numbers! */

    p = openIniFile (filename, 0);

/* Adjust line number if it doesn't make sense */

    if (line == 0) 
        line = 1;
    if (line > p->num_lines+1) 
        line = p->num_lines+1;

/* Make sure the line array has room for the new line */

    growIniLineArray (p);

/* Shuffle lines down in the array to make room for the new line */

    for (i = p->num_lines; i >= line; i--) 
        p->lines[i] = p->lines[i-1];
    p->num_lines++;

/* Allocate and fill in a new line structure */

    p->lines[line-1] = (struct IniLine *) malloc (sizeof (struct   
        IniLine));
    p->lines[line-1]->keyword = (char *) malloc (strlen (keyword) + 1);
    p->lines[line-1]->value = (char *) malloc (strlen (val) + 1);
    p->lines[line-1]->active = TRUE;
    strcpy (p->lines[line-1]->keyword, keyword);
    strcpy (p->lines[line-1]->value, val);

/* Write the INI file back to disk */

    writeIniFile (p);
}

void IniInsertLineAsInt (
    char    *filename,
    unsigned int line,
    char    *keyword,
    long    val)
{
    char    buf[20];
    
    sprintf (buf, "%ld", val);
    IniInsertLineAsString (filename, line, keyword, buf);
}

void IniAppendLineAsString (
    char    *filename,
    char    *keyword,
    char    *val)
{
    struct IniCache *p;
    
    p = openIniFile (filename, 0);
    IniInsertLineAsString (filename, p->num_lines+1, keyword, val);
}

void IniAppendLineAsInt (
    char    *filename,
    char    *keyword,
    long    val)
{
    char    buf[20];
    
    sprintf (buf, "%ld", val);
    IniAppendLineAsString (filename, keyword, buf);
}

void IniDeleteLine (
    char    *filename,
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

    for (i = line; i < p->num_lines; i++) 
        p->lines[i-1] = p->lines[i];
    p->num_lines--;

/* Write the INI file back to disk */

    writeIniFile (p);
}

void IniDeleteAllLines (
    char    *filename)
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

/****************************************************************************/
/*                Routines to read and write float values		    */
/****************************************************************************/

float IniGetFloat (			/* Get a floating point value from the global section of the INI file */
	const char *filename,
	const char *keyword,
	float	default_val)
{
	char	buf[20], defval[20];
	sprintf (defval, "%f", default_val);
	IniGetString ((char*)filename,(char*) keyword, buf, 20, defval);
	return ((float) atof (buf));
}

void IniWriteFloat (			/* Write a floating point value to the global section of the INI file */
	const char *filename,
	const char *keyword,
	float	val)
{
	/* Assume FLT_MAX is 3.40282e+038, the maximum significant digits that */
	/* can be stored in this buf is 12. ((sizeof(buf))-sizeof("-.E+038")) */
 	char	buf[20];
	sprintf (buf, "%11g", val);
 	IniWriteString ((char*)filename, (char*)keyword, buf);

}


/* Output string to screen or results file */

void OutputSomewhere (
    char    *buf)
{
    if (NO_GUI) writeResults (buf);
    else OutputStr (buf);
}

/* Output string to both the screen and results file */

void OutputBoth (
    char    *buf)
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
    char    *str)
{
    int	fd;
    unsigned long filelen;
    static  time_t  last_time = 0;
    time_t  this_time;

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
            char    buf[48];
            last_time = this_time;
            buf[0] = '[';
            strcpy (buf+1, ctime (&this_time));
            sprintf (buf+25, " - ver %s]\n", VERSION);
            wc = _write (fd, buf, strlen (buf));
        }

/* Output the message */

        wc = _write (fd, str, strlen (str));
    }

/* Display message about full log file */
	
    else {
        char    *fullmsg = (char*)"Prime.log file full.  Please delete it.\n";
        OutputStr (fullmsg);
        if (filelen < 251000)
            wc = _write (fd, fullmsg, strlen (fullmsg));
    }
    _close (fd);
}

int gmodi (uint32_t, giant);


/* Generate temporary file name */
 
void tempFileName ( 
    char *buf, 
    char c,
    giant NN) 
{ 
    int remainder;
 
    remainder = gmodi(19999981, NN);
    sprintf (buf, "%1c%07i", c, remainder % 10000000); 
} 

/* See if the given file exists */

int fileExists (
    char    *filename)
{
    int	fd;
    
    fd = _open (filename, _O_RDONLY | _O_BINARY);
    if (fd < 0) 
        return (0);
    _close (fd);
    return (1);
}

/* Open the results file and write a line to the end of it. */

int writeResults (
    char    *msg)
{
    static  time_t  last_time = 0;
    time_t  this_time;
    int	fd;

    if (IniGetInt (INI_FILE, (char*)"NoLresultFile", 0))
        return (TRUE);

/* Open file, position to end */

    fd = _open (RESFILE, _O_TEXT | _O_RDWR |
        _O_CREAT | _O_APPEND, 0666);
    if (fd < 0) {
        LogMsg ((char*)"Error opening the results file ; see result below :\n");
        LogMsg (msg);// Log the unwrited message (Darren Bedwell'request)
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
            wc = _write (fd, buf, 27);
    }

/* Output the message */

    if (_write (fd, msg, strlen (msg)) < 0) 
        goto fail;
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
	if (_read (fd, &len, sizeof (long)) != sizeof (long)) {
            gwypfree (tmp);
            return (FALSE);
        }
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
    int fd,
    gwypnum  g,
    long    *sum)
{
//  giant   tmp;
    long    len, bytes;

    if (_read (fd, &len, sizeof (long)) != sizeof (long)) 
        return (FALSE);
    bytes = len * sizeof (double);
    if (_read (fd, g, bytes) != bytes) 
        return (FALSE);
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
        if (_write (fd, &len, sizeof (long)) != sizeof (long)) {
            gwypfree (tmp);
            return (FALSE);
        }
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
    gwypnum g,
    long    *sum)
{
//  giant   tmp;
    long    len, bytes;

    len = FFTLEN;
    if (_write (fd, &len, sizeof (long)) != sizeof (long)) 
        return (FALSE);
    bytes = len * sizeof (double);
    if (_write (fd, g, bytes) != bytes) 
        return (FALSE);
    *sum = 0;
    return (TRUE);
}

int read_long (
    int	fd,
    unsigned long *val,
    long    *sum)
{
    if (_read (fd, val, sizeof (long)) != sizeof (long)) 
        return (FALSE);
    *sum += *val;
    return (TRUE);
}

int write_long (
    int	fd,
    unsigned long val,
    long    *sum)
{
    if (_write (fd, &val, sizeof (long)) != sizeof (long)) 
        return (FALSE);
    *sum += val;
    return (TRUE);
}

int read_double (
    int	fd,
    double *val,
    long    *sum)
{
    if (_read (fd, val, sizeof (double)) != sizeof (double)) 
        return (FALSE);
    *sum += (long)floor(*val);
    return (TRUE);
}

int write_double (
    int	fd,
    double val,
    long    *sum)
{
    if (_write (fd, &val, sizeof (double)) != sizeof (double)) 
        return (FALSE);
    *sum += (long)floor(val);
    return (TRUE);
}

int read_long298 (
	int	fd,
	unsigned long *val,
	unsigned long *sum)
{
	uint32_t tmp;

	if (_read (fd, &tmp, sizeof (uint32_t)) != sizeof (uint32_t))
		return (FALSE);
	if (sum != NULL) *sum = (uint32_t) (*sum + tmp);
	*val = tmp;
	return (TRUE);
}

int write_long298 (
	int	fd,
	unsigned long val,
	unsigned long *sum)
{
	uint32_t tmp;

	tmp = (uint32_t) val;
	if (_write (fd, &tmp, sizeof (uint32_t)) != sizeof (uint32_t))
		return (FALSE);
	if (sum != NULL) *sum = (uint32_t) (*sum + tmp);
	return (TRUE);
}


int read_slong (
	int	fd,
	long	*val,
	unsigned long *sum)
{
	int32_t tmp;

	if (_read (fd, &tmp, sizeof (int32_t)) != sizeof (int32_t))
		return (FALSE);
	if (sum != NULL) *sum = (uint32_t) (*sum + (uint32_t) tmp);
	*val = tmp;
	return (TRUE);
}

int write_slong (
	int	fd,
	long	val,
	unsigned long *sum)
{
	int32_t tmp;

	tmp = (int32_t) val;
	if (_write (fd, &tmp, sizeof (int32_t)) != sizeof (int32_t))
		return (FALSE);
	if (sum != NULL) *sum = (uint32_t) (*sum + (uint32_t) tmp);
	return (TRUE);
}

int read_double298 (
	int	fd,
	double	*val,
	unsigned long *sum)
{
	if (_read (fd, val, sizeof (double)) != sizeof (double))
		return (FALSE);
	if (sum != NULL) *sum = (uint32_t) (*sum + (uint32_t) *val);
	return (TRUE);
}

int write_double298 (
	int	fd,
	double	val,
	unsigned long *sum)
{
	if (_write (fd, &val, sizeof (double)) != sizeof (double))
		return (FALSE);
	if (sum != NULL) *sum += (uint32_t) (*sum + (uint32_t) val);
	return (TRUE);
}

int read_gwypnumjp (
    int fd,
    gwypnum  g,
    unsigned long    *sum)
{
    long    len, i;

    if (_read (fd, &len, sizeof (long)) != sizeof (long)) 
        return (FALSE);
    *sum += len;
    for (i=0;i < len;i++) {
        if (!read_double298 (fd, &g[i], 0))
            return (FALSE);
        *sum += (uint32_t) (*sum + (uint32_t) g[i]);
    }
    return (TRUE);
}

int write_gwypnumjp (
    int	fd,
    gwypnum g,
    unsigned long    *sum)
{
    long    len, i;

    len = FFTLEN;
    if (_write (fd, &len, sizeof (long)) != sizeof (long)) 
        return (FALSE);
    *sum += len;
    for (i=0;i < len;i++) {
        if (!write_double298 (fd, g[i], 0))
            return (FALSE);
        *sum += (uint32_t) (*sum + (uint32_t) g[i]);
    }
    return (TRUE);
}

/* Routines to read and write the common header portion of all save files */
/* The save file format is: */
/*	u32		magic number  (different for ll, p-1, prp, tf, ecm) */
/*	u32		version number */
/*	double		k in k*b^n+c */
/*	u32		b in k*b^n+c */
/*	u32		n in k*b^n+c */
/*	s32		c in k*b^n+c */
/*	double		pct complete */
/*	char(11)	stage */
/*	char(1)		pad */
/*	u32		checksum of all following data */

/* Read and test the "magic number" in the first 4 bytes of the file. */
/* We use this to detect save files created by this program prior to the */
/* common header format. */

/* Routines to read and write a byte array from and to a save file */

int read_array (
	int	fd,
	char	*buf,
	unsigned long len,
	unsigned long *sum)
{
	unsigned long i;
	unsigned char *ubuf;

	if ((unsigned long)_read (fd, buf, len) != len) return (FALSE);
	ubuf = (unsigned char *) buf;
	if (sum != NULL)
		for (i = 0; i < len; i++)
			*sum = (uint32_t) (*sum + ubuf[i]);
	return (TRUE);
}

int write_array (
	int	fd,
	const char *buf,
	unsigned long len,
	unsigned long *sum)
{
	unsigned long i;
	unsigned char *ubuf;

	if (len == 0) return (TRUE);
	if ((unsigned long)_write (fd, buf, len) != len) return (FALSE);
	ubuf = (unsigned char *) buf;
	if (sum != NULL)
		for (i = 0; i < len; i++)
			*sum = (uint32_t) (*sum + ubuf[i]);
	return (TRUE);
}

#define CHECKSUM_OFFSET	72

int read_checksum (
	int	fd,
	unsigned long *sum)
{
	_lseek (fd, CHECKSUM_OFFSET, SEEK_SET);
	if (!read_long298 (fd, sum, NULL)) return (FALSE);
	return (TRUE);
}

int write_checksum (
	int	fd,
	unsigned long sum)
{
	_lseek (fd, CHECKSUM_OFFSET, SEEK_SET);
	if (!write_long298 (fd, sum, NULL)) return (FALSE);
	return (TRUE);
}

int read_magicnum (
	int	fd,
	unsigned long magicnum)
{
	unsigned long filenum;

/* Read the magic number from the first 4 bytes */

	_lseek (fd, 0, SEEK_SET);
	if (!read_long298 (fd, &filenum, NULL)) return (FALSE);

/* Return TRUE if the magic number matches the caller's desired magic number */

	return (filenum == magicnum);
}

/* Read the rest of the common header */

int read_header (
	int	fd,
	unsigned long *version,
	struct work_unit *w,
	unsigned long *sum)
{
	double	k;
	unsigned long b, n;
	long	c;
	char	pad[5];
	char	stage[11];
	double	pct_complete;
	unsigned long trash_sum;
;

/* Skip past the magic number in the first 4 bytes */

	_lseek (fd, sizeof (uint32_t), SEEK_SET);

/* Read the header */

	if (!read_long298 (fd, version, NULL)) return (FALSE);
	if (!read_double298 (fd, &k, NULL)) return (FALSE);
	if (!read_long298 (fd, &b, NULL)) return (FALSE);
	if (!read_long298 (fd, &n, NULL)) return (FALSE);
	if (!read_slong (fd, &c, NULL)) return (FALSE);
	if (!read_array (fd, stage, 11, NULL)) return (FALSE);
	if (!read_array (fd, pad, 5, NULL)) return (FALSE);
	if (!read_double298 (fd, &pct_complete, NULL)) return (FALSE);
	if (sum == NULL) sum = &trash_sum;
	_lseek (fd, CHECKSUM_OFFSET, SEEK_SET);
        if (_read (fd, sum, sizeof (unsigned long)) != sizeof (unsigned long)) return (FALSE);
        
/* Validate the k,b,n,c values */

	if (k != w->k || b != w->b || n != w->n || c != w->c) return (FALSE);

/* Set the work unit's stage and pct_complete fields */

	stage[10] = 0;
	strcpy (w->stage, stage);
	if (pct_complete < 0.0) pct_complete = 0.0;
	if (pct_complete > 1.0) pct_complete = 1.0;
	w->pct_complete = pct_complete;

/* Return success */

	return (TRUE);
}

int write_header (
	int	fd,
	unsigned long magicnum,
	unsigned long version,
	struct work_unit *w)
{
	char   pad[5] = {0};
        // to have a header size multiple of 8
	unsigned long sum = 0;

	if (!write_long298 (fd, magicnum, NULL)) return (FALSE);
	if (!write_long298 (fd, version, NULL)) return (FALSE);
	if (!write_double298 (fd, w->k, NULL)) return (FALSE);
	if (!write_long298 (fd, w->b, NULL)) return (FALSE);
	if (!write_long298 (fd, w->n, NULL)) return (FALSE);
	if (!write_slong (fd, w->c, NULL)) return (FALSE);
	if (!write_array (fd, w->stage, 11, NULL)) return (FALSE);
	if (!write_array (fd, pad, 5, NULL)) return (FALSE);
	if (!write_double298 (fd, w->pct_complete, NULL)) return (FALSE);
//	if (!write_long298 (fd, sum, NULL)) return (FALSE);
	_lseek (fd, CHECKSUM_OFFSET, SEEK_SET);
	if (_write (fd, &sum, sizeof (unsigned long)) != sizeof (unsigned long)) return (FALSE);
	return (TRUE);
}


int writeToFile (
    char    *filename,
    unsigned long j,
    gwypnum x,
    gwypnum y)
{
    char    newfilename[16], errmsg[100];
    int	fd;
    unsigned long magicnum, version;
    long    sum = 0, i;


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
    if (fd < 0) 
        return (FALSE);
    
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

    if (!write_long (fd, j, &sum)) 
        goto writeerr;

/* Write the data values */

    if (!write_gwypnum (fd, x, &sum)) 
        goto writeerr;
    if (y != NULL && !write_gwypnum (fd, y, &sum))
        goto writeerr; 

/* Write the checksum */

    if (_write (fd, &sum, sizeof (long)) != sizeof (long)) 
        goto writeerr;

/* Save the timers */

    for (i=0; i<NBTIMERS; i++) {
        if (gwyptimers[i+NBTIMERS] != 0.0) {
            // if the timer was running
            gwypend_timer (i);// update and save it
            if (! write_double (fd, gwyptimers[i], &sum)) {
                gwypstart_timer (i);
            // and then, restart it, even if write is in error!
                goto writeerr;
            }
            if (! write_double (fd, gwyptimers[i+NBTIMERS], &sum)) { 
                // save the timer status
                gwypstart_timer (i);
                // and then, restart it, even if write is in error!
                goto writeerr;
            }
            gwypstart_timer (i);
                // and then, restart it!
        }
        else {
            if (! write_double (fd, gwyptimers[i], &sum)) 
                goto writeerr;  // save the timer
            if (! write_double (fd, gwyptimers[i+NBTIMERS], &sum)) 
                goto writeerr;	// save its status
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
    char   *filename,
    unsigned long *j,
    gwypnum  x,
    gwypnum  y)
{
    int  fd;
    unsigned long magicnum, version;
//       long    sum = 0, i;
    char errmsg[100];

/* Open the intermediate file */

    fd = _open (filename, _O_BINARY | _O_RDONLY);
    if (fd < 0) goto error;

/* Read the file header */

    if (_read (fd, &magicnum, sizeof (long)) != sizeof (long))
        goto readerr;
    if (magicnum != 0x9f2b3cd4) 
        goto readerr;

    if (_read (fd, &version, sizeof (long)) != sizeof (long)) 
        goto readerr;
    if (version != 1 && version != 2) 
        goto readerr;
    if (_read (fd, &s_FFTLEN, sizeof (int)) != sizeof (int)) 
        goto readerr; //cuda
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
    char    *filename,
    unsigned long *j,
    gwypnum	x,
    gwypnum	y)
{
    int	fd;
    unsigned long magicnum, version;
    long    sum = 0, i;
    char    errmsg[100];

/* Open the intermediate file */

    fd = _open (filename, _O_BINARY | _O_RDONLY);
    if (fd < 0)
        goto error;

/* Read the file header */

    if (_read (fd, &magicnum, sizeof (long)) != sizeof (long))
        goto readerr;
    if (magicnum != 0x9f2b3cd4)
        goto readerr;
    if (_read (fd, &version, sizeof (long)) != sizeof (long))
        goto readerr;
    if (version != 1 && version != 2)
        goto readerr;
    if (_read (fd, &s_FFTLEN, sizeof (int)) != sizeof (int))
        goto readerr; //cuda

/* Read the file data */

    if (!read_long (fd, j, &sum))
        goto readerr;

/* Read the values */

    if (!read_gwypnum (fd, x, &sum))
        goto readerr;
    if (y != NULL && !read_gwypnum (fd, y, &sum))
        goto readerr; 

/* Read and compare the checksum */

    if (_read (fd, &i, sizeof (long)) != sizeof (long))
        goto readerr;
    if (i != sum)
        goto readerr;

/* Read the timers and their status */

    for (i=0; i<NBTIMERS; i++) {
        if (!read_double (fd, &gwyptimers[i], &sum))
            goto readerr;
        if (!read_double (fd, &gwyptimers[i+NBTIMERS], &sum))
            goto readerr;
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
    char    *filename,
    unsigned long j,
    unsigned long B,
    unsigned long nr,
    unsigned long *bpf,
    gwypnum x,
    gwypnum y)
{
    char    newfilename[16],errmsg[100];
    int	    fd;
    unsigned long magicnum, version;
    long    sum = 0, i;


    if (debug) {
        sprintf (errmsg, "Writing %s intermediate file, j = %lu\n",filename, j);
        OutputBoth (errmsg);
    }

/* If we are allowed to create multiple intermediate files, then */
/* write to a file called yNNNNNNN. */

    strcpy (newfilename, filename);
    if (TWO_BACKUP_FILES)
        newfilename[0] = 'y';

/* Create the intermediate file */

    fd = _open (newfilename, _O_BINARY|_O_WRONLY|_O_TRUNC|_O_CREAT, 0666);
    if (fd < 0)
        return (FALSE);

/* Write the file header. */

    magicnum = 0x9f2b3cd4;
    if (_write (fd, &magicnum, sizeof (long)) != sizeof (long))
        goto writeerr;
    version = 1;
    if (_write (fd, &version, sizeof (long)) != sizeof (long))
        goto writeerr;

/* Write the file data */

    if (! write_long (fd, j, &sum))
        goto writeerr;
    if (! write_long (fd, B, &sum))
        goto writeerr;
    if (! write_long (fd, nr, &sum))
        goto writeerr;
    for (i=0; i<10; i++) {
        if (! write_long (fd, bpf[i], &sum))
            goto writeerr;
    }

/* Write the data values */

    if (! write_gwypnum (fd, x, &sum))
        goto writeerr;
    if (y != NULL && ! write_gwypnum (fd, y, &sum)) 
        goto writeerr; 

/* Write the checksum */

    if (_write (fd, &sum, sizeof (long)) != sizeof (long))
        goto writeerr;

/* Save the timers */

    for (i=0; i<NBTIMERS; i++) {
        if (gwyptimers[i+NBTIMERS] != 0.0) {
                    // if the timer was running
            gwypend_timer (i);// update and save it
            if (! write_double (fd, gwyptimers[i], &sum)) {
                gwypstart_timer (i);
                // and then, restart it, even if write is in error!
                goto writeerr;
            }
            if (! write_double (fd, gwyptimers[i+NBTIMERS], &sum)) {
                // save the timer status
                gwypstart_timer (i);
                // and then, restart it, even if write is in error!
                goto writeerr;
            }
            gwypstart_timer (i);
                // and then, restart it!
        }
        else {
            if (! write_double (fd, gwyptimers[i], &sum))
                goto writeerr;	// save the timer
            if (! write_double (fd, gwyptimers[i+NBTIMERS], &sum))
                goto writeerr;	// save its status
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
    char    *filename,
    unsigned long *j,
    unsigned long *B,
    unsigned long *nr,
    unsigned long *bpf,
    gwypnum x,
    gwypnum y)
{
    int	    fd;
    unsigned long magicnum, version;
    long    sum = 0, i;
    char    errmsg[100];

/* Open the intermediate file */

    fd = _open (filename, _O_BINARY | _O_RDONLY);
    if (fd < 0)
        goto error;

/* Read the file header */

    if (_read (fd, &magicnum, sizeof (long)) != sizeof (long))
        goto readerr;
    if (magicnum != 0x9f2b3cd4)
        goto readerr;
    if (_read (fd, &version, sizeof (long)) != sizeof (long))
        goto readerr;
    if (version != 1 && version != 2)
        goto readerr;

/* Read the file data */

    if (! read_long (fd, j, &sum))
        goto readerr;
    if (! read_long (fd, B, &sum))
        goto readerr;
    if (! read_long (fd, nr, &sum))
        goto readerr;
    for (i=0; i<10; i++) {
        if (! read_long (fd, &bpf[i], &sum))
            goto readerr;
    }

/* Read the values */

    if (! read_gwypnum (fd, x, &sum))
        goto readerr;
    if (y != NULL && ! read_gwypnum (fd, y, &sum))
        goto readerr; 

/* Read and compare the checksum */

    if (_read (fd, &i, sizeof (long)) != sizeof (long))
        goto readerr;
    if (i != sum)
        goto readerr;

/* Read the timers and their status */

    for (i=0; i<NBTIMERS; i++) {
        if (! read_double (fd, &gwyptimers[i], &sum))
            goto readerr;
        if (! read_double (fd, &gwyptimers[i+NBTIMERS], &sum))
            goto readerr;
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
    char    *filename,
    unsigned long j,
    unsigned long ubx,
    unsigned long uby,
    gwypnum x,
    gwypnum y)
{
    char    newfilename[16], errmsg[100];
    int	    fd;
    unsigned long magicnum, version;
    long    sum = 0, i;


    if (debug) {
        sprintf (errmsg, "Writing %s intermediate file, j = %lu\n",filename, j);
        OutputBoth (errmsg);
    }

/* If we are allowed to create multiple intermediate files, then */
/* write to a file called yNNNNNNN. */

    strcpy (newfilename, filename);
    if (TWO_BACKUP_FILES)
        newfilename[0] = 'y';

/* Create the intermediate file */

    fd = _open (newfilename, _O_BINARY|_O_WRONLY|_O_TRUNC|_O_CREAT, 0666);
    if (fd < 0)
        return (FALSE);

/* Write the file header. */

    magicnum = 0x9f2b3cd4;
    if (_write (fd, &magicnum, sizeof (long)) != sizeof (long))
        goto writeerr;
    version = 1;
    if (_write (fd, &version, sizeof (long)) != sizeof (long))
        goto writeerr;

/* Write the file data */

    if (! write_long (fd, j, &sum))
        goto writeerr;
    if (! write_long (fd, ubx, &sum))
        goto writeerr;
    if (! write_long (fd, uby, &sum))
        goto writeerr;

/* Write the data values */

    if (! write_gwypnum (fd, x, &sum))
        goto writeerr;
    if (y != NULL && ! write_gwypnum (fd, y, &sum))
        goto writeerr; 

/* Write the checksum */

    if (_write (fd, &sum, sizeof (long)) != sizeof (long))
        goto writeerr;

/* Save the timers */

    for (i=0; i<NBTIMERS; i++) {
        if (gwyptimers[i+NBTIMERS] != 0.0) {
                // if the timer was running
            gwypend_timer (i);
                // update and save it
            if (! write_double (fd, gwyptimers[i], &sum)) {
                gwypstart_timer (i);
                // and then, restart it, even if write is in error!
                goto writeerr;
            }
            if (! write_double (fd, gwyptimers[i+NBTIMERS], &sum)) {
                // save the timer status
                gwypstart_timer (i);
                // and then, restart it, even if write is in error!
                goto writeerr;
            }
            gwypstart_timer (i);
            // and then, restart it!
        }
        else {
            if (! write_double (fd, gwyptimers[i], &sum))
                goto writeerr;// save the timer
            if (! write_double (fd, gwyptimers[i+NBTIMERS], &sum))
                goto writeerr;// save its status
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
    char    *filename,
    unsigned long *j,
    unsigned long *ubx,
    unsigned long *uby,
    gwypnum x,
    gwypnum y)
{
    int	fd;
    unsigned long magicnum, version;
    long    sum = 0, i;
    char    errmsg[100];

/* Open the intermediate file */

    fd = _open (filename, _O_BINARY | _O_RDONLY);
    if (fd < 0)
        goto error;

/* Read the file header */

    if (_read (fd, &magicnum, sizeof (long)) != sizeof (long))
        goto readerr;
    if (magicnum != 0x9f2b3cd4)
        goto readerr;
    if (_read (fd, &version, sizeof (long)) != sizeof (long))
        goto readerr;
    if (version != 1 && version != 2)
        goto readerr;

/* Read the file data */

    if (! read_long (fd, j, &sum))
        goto readerr;
    if (! read_long (fd, ubx, &sum))
        goto readerr;
    if (! read_long (fd, uby, &sum))
        goto readerr;

/* Read the values */

    if (! read_gwypnum (fd, x, &sum))
        goto readerr;
    if (y != NULL && ! read_gwypnum (fd, y, &sum))
        goto readerr; 

/* Read and compare the checksum */

    if (_read (fd, &i, sizeof (long)) != sizeof (long))
        goto readerr;
    if (i != sum)
        goto readerr;

/* Read the timers and their status */

    for (i=0; i<NBTIMERS; i++) {
        if (! read_double (fd, &gwyptimers[i], &sum))
            goto readerr;
        if (! read_double (fd, &gwyptimers[i+NBTIMERS], &sum))
            goto readerr;
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
    char    *filename,
    unsigned long j,
    unsigned long D,
    unsigned long nr,
    unsigned long *bpf,
    gwypnum x,
    gwypnum y,
    gwypnum z,
    gwypnum t)
{
    char   newfilename[16], errmsg[100];
    int    fd;
    unsigned long magicnum, version;
    long   sum = 0, i;


    if (debug) {
        sprintf (errmsg, "Writing %s intermediate file, j = %lu\n",filename, j);
        OutputBoth (errmsg);
    }

/* If we are allowed to create multiple intermediate files, then */
/* write to a file called yNNNNNNN. */

    strcpy (newfilename, filename);
    if (TWO_BACKUP_FILES)
        newfilename[0] = 'y';

/* Create the intermediate file */

    fd = _open (newfilename, _O_BINARY|_O_WRONLY|_O_TRUNC|_O_CREAT, 0666);
	if (fd < 0)
            return (FALSE);

/* Write the file header. */

    magicnum = 0x9f2b3cd4;
    if (_write (fd, &magicnum, sizeof (long)) != sizeof (long))
        goto writeerr;
    version = 1;
    if (_write (fd, &version, sizeof (long)) != sizeof (long))
        goto writeerr;

/* Write the file data */

    if (! write_long (fd, j, &sum))
        goto writeerr;
    if (! write_long (fd, D, &sum))
        goto writeerr;
    if (! write_long (fd, nr, &sum))
        goto writeerr;
    for (i=0; i<10; i++) {
        if (! write_long (fd, bpf[i], &sum))
            goto writeerr;
    }

/* Write the data values */

    if (! write_gwypnum (fd, x, &sum))
        goto writeerr;
    if (y != NULL && ! write_gwypnum (fd, y, &sum))
        goto writeerr; 
    if (z != NULL && ! write_gwypnum (fd, z, &sum))
        goto writeerr; 
    if (t != NULL && ! write_gwypnum (fd, t, &sum))
        goto writeerr; 

/* Write the checksum */

    if (_write (fd, &sum, sizeof (long)) != sizeof (long))
        goto writeerr;

/* Save the five gwyptimers */

    for (i=0; i<NBTIMERS; i++) {
        if (gwyptimers[i+NBTIMERS] != 0.0) {
                    // if the timer was running
            gwypend_timer (i);// update and save it
            if (! write_double (fd, gwyptimers[i], &sum)) {
                gwypstart_timer (i);
                    // and then, restart it, even if write is in error!
                goto writeerr;
            }
            if (! write_double (fd, gwyptimers[i+NBTIMERS], &sum)) {
                    // save the timer status
                gwypstart_timer (i);
                    // and then, restart it, even if write is in error!
                goto writeerr;
            }
            gwypstart_timer (i);
                    // and then, restart it!
        }
        else {
            if (! write_double (fd, gwyptimers[i], &sum))
                goto writeerr;	// save the timer
            if (! write_double (fd, gwyptimers[i+NBTIMERS], &sum))
                goto writeerr;	// save its status
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

int LreadFromFile (					// To restore a Lucas sequence matrix
    char    *filename,
    unsigned long *j,
    unsigned long *D,
    unsigned long *nr,
    unsigned long *bpf,
    gwypnum x,
    gwypnum y,
    gwypnum z,
    gwypnum t)
{
    int	fd;
    unsigned long magicnum, version;
    long    sum = 0, i;
    char    errmsg[100];

/* Open the intermediate file */

    fd = _open (filename, _O_BINARY | _O_RDONLY);
    if (fd < 0)
        goto error;

/* Read the file header */

    if (_read (fd, &magicnum, sizeof (long)) != sizeof (long))
        goto readerr;
    if (magicnum != 0x9f2b3cd4)
        goto readerr;
    if (_read (fd, &version, sizeof (long)) != sizeof (long))
        goto readerr;
    if (version != 1 && version != 2)
        goto readerr;

/* Read the file data */

    if (! read_long (fd, j, &sum))
        goto readerr;
    if (! read_long (fd, D, &sum))
        goto readerr;
    if (! read_long (fd, nr, &sum))
        goto readerr;
    for (i=0; i<10; i++) {
        if (! read_long (fd, &bpf[i], &sum))
            goto readerr;
    }

/* Read the values */

    if (! read_gwypnum (fd, x, &sum))
        goto readerr;
    if (y != NULL && ! read_gwypnum (fd, y, &sum))
        goto readerr; 
    if (z != NULL && ! read_gwypnum (fd, z, &sum))
        goto readerr; 
    if (t != NULL && ! read_gwypnum (fd, t, &sum))
        goto readerr; 

/* Read and compare the checksum */

    if (_read (fd, &i, sizeof (long)) != sizeof (long))
        goto readerr;
    if (i != sum)
        goto readerr;

/* Read the timers and their status */

    for (i=0; i<NBTIMERS; i++) {
        if (! read_double (fd, &gwyptimers[i], &sum))
            goto readerr;
        if (! read_double (fd, &gwyptimers[i+NBTIMERS], &sum))
            goto readerr;
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
    gwypnum gg,
    int N
)
{
    int j;
    long    val;
    char buf[20];

    OutputStr ((char*)"\n");
    for(j=0; j<N; ++j)  {
        val = (long)gg[j];
        if (val) {
            sprintf (buf, "%ld ", val);
            OutputBoth (buf);
        }
    }
    OutputBoth ((char*)"\n");
    return 0;
}

// There are (at least) 8 PRP or Prime residue types for testing N=(k*b^n+c)/d:

#define	PRP_TYPE_FERMAT		1	// Fermat PRP.  Calculate a^(N-1) mod N.  PRP if result = 1
#define	PRP_TYPE_SPRP		2	// SPRP variant.  Calculate a^((N-1)/2) mod N.  PRP if result = +/-1
#define	PRP_TYPE_FERMAT_VAR	3	// Type 1 variant,b=2,d=1. Calculate a^(N-c) mod N.  PRP if result = a^-(c-1)
#define	PRP_TYPE_SPRP_VAR	4	// Type 2 variant,b=2,d=1. Calculate a^((N-c)/2) mod N.  PRP if result = +/-a^-((c-1)/2)
#define	PRP_TYPE_COFACTOR	5	// Cofactor variant.  Calculate a^(N*d-1) mod N*d.  PRP if result = a^(d-1) mod N
#define PROTH_TYPE			6	// Proth prime variant. Calculate a^((N-1)/2) mod N with Jacobi (a,N) = -1. Prime iff result = -1
#define	GMN_TYPE			7	// A special variant of Proth prime test
#define	WAGSTAFF_TYPE		8	// A special variant of Cofactor SPRP test 

// Primenet encourages programs to return type 1 PRP residues as that has been the standard for prime95, PFGW, LLR for many years.

#define PRP_MAGICNUM		0x87f2a91b
#define PRP_VERSION		4

#ifndef SEC1
#define SEC1(p)			0UL
#define SEC2(p,hi,lo,u,e)	0UL
#define SEC3(p)			0UL
#define SEC4(p)			0UL
#define SEC5(p,b1,b2)		0UL
#endif

char res64[40];
    /* VP : This variable has been made global */
char alt_res64[40]; /* JP 01/07/20 : This variable has been made global */

/* Data structure used in reading save files and their backups as well as */
/* renaming bad save files. */

typedef struct read_save_file_state {
	int	thread_num;
	int	read_attempt;
	int	a_save_file_existed;
	int	a_non_bad_save_file_existed;
	int	num_original_bad_files;
	int	num_save_files_renamed;
	char	base_filename[80];
	char	current_filename[80];
} readSaveFileState;

/* Data structure used in writing save files and their backups */

typedef struct write_save_file_state {
	char	base_filename[80];
	int	num_ordinary_save_files;
	int	num_special_save_files;		/* Example: Number of save files to keep that passed the Jacobi error check */
	uint64_t special;			/* Bit array for which ordinary save files are special */
} writeSaveFileState;

/* Increment the error counter.  The error counter is one 32-bit field containing 5 values.  Prior to version 29.3, this was */
/* a one-bit flag if this is a continuation from a save file that did not track error counts, a 7-bit count of errors that were */
/* reproducible, a 8-bit count of ILLEGAL SUMOUTs or zeroed FFT data or corrupt units_bit, a 8-bit count of convolution errors */
/* above 0.4, and a 8-bit count of SUMOUTs not close enough to SUMINPs. */
/* NOTE:  The server considers an LL run clean if the error code is XXaaYY00 and XX = YY and aa is ignored.  That is, repeatable */
/* round off errors and all ILLEGAL SUMOUTS are ignored. */
/* In version 29.3, a.k.a. Wf in result lines, the 32-bit field changed.  See comments in the code below. */

void inc_error_count (
	int	type,
	unsigned long *error_count)
{
	unsigned long addin, orin, maxval;

	addin = orin = 0;
	if (type == 0) addin = 1, maxval = 0xF;				// SUMINP != SUMOUT
	else if (type == 4) addin = 1 << 4, maxval = 0x0F << 4;		// Jacobi error check
	else if (type == 1) addin = 1 << 8, maxval = 0x3F << 8;		// Roundoff > 0.4
	else if (type == 5) orin = 1 << 14;				// Zeroed FFT data
	else if (type == 6) orin = 1 << 15;				// Units bit, counter, or other value corrupted
	else if (type == 2) addin = 1 << 16, maxval = 0xF << 16;	// ILLEGAL SUMOUT
	else if (type == 7) addin = 1 << 20, maxval = 0xF << 20;	// High reliability (Gerbicz or dblchk) PRP error
	else if (type == 3) addin = 1 << 24, maxval = 0x3F << 24;	// Repeatable error

	if (addin && (*error_count & maxval) != maxval) *error_count += addin;
	*error_count |= orin;
}

/* Create a message if the non-repeatable error count is more than zero */
/* Returns TRUE if the non-repeatable error count is more than zero. */

int make_error_count_message (
	unsigned long error_count,
	int	message_type,		/* 1 = very small, 2 = one line, 3 = multi-line, 0x8000 = confidence always excellent */
	char	*buf,
	int	buflen)
{
	int	count_repeatable, count_suminp, count_roundoff, count_illegal_sumout, count_total;
	int	count_jacobi, count_gerbicz, count_bad_errors;
	int	force_high_confidence;
	char	local_buf[400], counts_buf[200], confidence[25];

/* Massage input argument */

	force_high_confidence = message_type & 0x8000;
	message_type = message_type & 0x7FFF;

/* Parse the error counts variable */

	count_repeatable = (error_count >> 24) & 0x3F;
	count_illegal_sumout = (error_count >> 16) & 0xF;
	count_roundoff = (error_count >> 8) & 0x3F;
	count_jacobi = (error_count >> 4) & 0xF;
	count_gerbicz = (error_count >> 20) & 0xF;
	count_suminp = error_count & 0xF;

/* Return if no hardware errors have occurred */

	count_total = count_illegal_sumout + count_suminp + count_roundoff + count_jacobi + count_gerbicz;
	if (count_total - count_repeatable == 0) return (FALSE);

/* Format the error counts */

	counts_buf[0] = 0;

	if (message_type == 1) {
		sprintf (counts_buf, ", errors: %d", count_total);
	}

	if (message_type == 2) {
		sprintf (counts_buf, "%d error%s", count_total, count_total > 1 ? "s": "");
		if (count_repeatable == 1) {
			strcpy (local_buf, ", 1 was repeatable (not an error)");
			strcat (counts_buf, local_buf);
		} else if (count_repeatable > 1) {
			sprintf (local_buf, ", %d%s were repeatable (not errors)",
				 count_repeatable, count_repeatable == 0x3F ? " or more" : "");
			strcat (counts_buf, local_buf);
		}
	}

	if (message_type == 3) {
		if (count_jacobi >= 1) {
			sprintf (local_buf, "%d%s Jacobi error%s, ",
				 count_jacobi, count_jacobi == 0xF ? " or more" : "", count_jacobi > 1 ? "s" : "");
			strcat (counts_buf, local_buf);
		}
		if (count_gerbicz >= 1) {
			sprintf (local_buf, "%d%s Gerbicz/double-check error%s, ",
				 count_gerbicz, count_gerbicz == 0xF ? " or more" : "", count_gerbicz > 1 ? "s" : "");
			strcat (counts_buf, local_buf);
		}
		if (count_roundoff >= 1) {
			sprintf (local_buf, "%d%s ROUNDOFF > 0.4, ",
				 count_roundoff, count_roundoff == 0x3F ? " or more" : "");
			strcat (counts_buf, local_buf);
		}
		if (count_suminp >= 1) {
			sprintf (local_buf, "%d%s SUM(INPUTS) != SUM(OUTPUTS), ",
				 count_suminp, count_suminp == 0xF ? " or more" : "");
			strcat (counts_buf, local_buf);
		}
		if (count_illegal_sumout >= 1) {
			sprintf (local_buf, "%d%s ILLEGAL SUMOUT/bad FFT data, ",
				 count_illegal_sumout, count_illegal_sumout == 0xF ? " or more" : "");
			strcat (counts_buf, local_buf);
		}
		counts_buf[strlen(counts_buf)-2] = 0;
		if (count_repeatable >= 1) {
			if (count_repeatable == 1)
				strcpy (local_buf, "of which 1 was repeatable (not a hardware error)");
			else
				sprintf (local_buf, "of which %d were repeatable (not hardware errors)", count_repeatable);
			if (strlen (counts_buf) <= 40) strcat (counts_buf, " ");
			else strcat (counts_buf, "\n");
			strcat (counts_buf, local_buf);
		}
		strcat (counts_buf, ".\n");
	}

/* Guess our confidence in the end result */

	count_bad_errors = count_jacobi + count_suminp + count_roundoff - count_repeatable;
	if (force_high_confidence) count_bad_errors = 0;
	strcpy (confidence, count_bad_errors == 0 ? "excellent" :
			    count_bad_errors <= 3 ? "fair" :
			    count_bad_errors <= 6 ? "poor" : "very poor");

/* Put it all together to form our full message */

	if (message_type == 1) {
		sprintf (local_buf, ", %s, confidence: %s", counts_buf, confidence);
	}
	if (message_type == 2) {
		if (count_jacobi || count_gerbicz) sprintf (local_buf, "Hardware errors!  %s.  Confidence in end result is %s.\n", counts_buf, confidence);
		else sprintf (local_buf, "Possible hardware errors!  %s.  Confidence in end result is %s.\n", counts_buf, confidence);
	}
	if (message_type == 3) {
		if (count_jacobi || count_gerbicz) strcpy (local_buf, "Hardware errors have occurred during the test!");
		else strcpy (local_buf, "Possible hardware errors have occurred during the test!");
		if (strlen (counts_buf) <= 25) strcat (local_buf, " ");
		else strcat (local_buf, "\n");
		strcat (local_buf, counts_buf);
		sprintf (local_buf+strlen(local_buf), "Confidence in final result is %s.\n", confidence);
	}

/* Copy as much of our result as possible to the caller's buffer */

	if ((int) strlen (local_buf) >= buflen) local_buf[buflen-1] = 0;
	strcpy (buf, local_buf);
	return (TRUE);
}

int openWriteSaveFile (
	writeSaveFileState *state)
{
	char	output_filename[32];
	int	fd;

/* If we are allowed to create multiple intermediate files, then use a .write extension */
/* The value 99, not accessible via the GUI, is a special value meaning overwrite the */
/* existing save file -- a very dangerous choice.  You might use this for a floppy or */
/* small USB stick installation where there is no room for two save files. */
/* NOTE: This behavior is different than v24 where when the user selected one save */
/* file, then he got the dangerous overwrite option. */

	if (state->num_ordinary_save_files == 99)
		strcpy (output_filename, state->base_filename);
	else
		sprintf (output_filename, "%s.write", state->base_filename);

/* Now save to the intermediate file */

	fd = _open (output_filename, _O_BINARY | _O_WRONLY | _O_TRUNC | _O_CREAT, CREATE_FILE_ACCESS);
	return (fd);
}

/* Close the save file we finished writing.  If necessary, delete old */
/* save file, and rename just written save file. */

void closeWriteSaveFile (
	writeSaveFileState *state,
	int	fd)
{
	char	dest_filename[32], src_filename[32];
	int	rename_count;

/* Flush data to disk and close the save file. */

	_commit (fd);
	_close (fd);

/* If no renaming is needed, we're done */

	if (state->num_ordinary_save_files == 99) return;

/* Save files that are special will be one step further down the chain after renaming */

	state->special <<= 1;

/* Decide how many save files need renaming (does the last ordinary file deserve to move into the special save files?) */

	rename_count = bittst (&state->special, state->num_ordinary_save_files) ?
			       state->num_ordinary_save_files + state->num_special_save_files : state->num_ordinary_save_files;

/* Delete the last file in the rename chain */

	if (rename_count == 1) strcpy (dest_filename, state->base_filename);
	else if (rename_count == 2) sprintf (dest_filename, "%s.bu", state->base_filename);
	else sprintf (dest_filename, "%s.bu%d", state->base_filename, rename_count-1);
	_unlink (dest_filename);

/* Perform the proper number of renames */

	while (rename_count--) {
		if (rename_count == 0) sprintf (src_filename, "%s.write", state->base_filename);
		else if (rename_count == 1) strcpy (src_filename, state->base_filename);
		else if (rename_count == 2) sprintf (src_filename, "%s.bu", state->base_filename);
		else sprintf (src_filename, "%s.bu%d", state->base_filename, rename_count-1);
		rename (src_filename, dest_filename);
		strcpy (dest_filename, src_filename);
	}
}

/* Mark the current save file as special (a super good save file -- Jacobi or Gerbicz checked) */

void setWriteSaveFileSpecial (
	writeSaveFileState *state)
{
	state->special |= 1;
}

/* Close and delete the save file we were writing.  This is done */
/* when an error occurs while writing the save file. */

void deleteWriteSaveFile (
	writeSaveFileState *state,
	int	fd)
{
	char	output_filename[32];

/* Close and delete the save file */

	_close (fd);
	if (state->num_ordinary_save_files == 99)
		strcpy (output_filename, state->base_filename);
	else
		sprintf (output_filename, "%s.write", state->base_filename);
	_unlink (output_filename);
}

/* Delete save files when work unit completes. */

void unlinkSaveFiles (
	writeSaveFileState *state)
{
	int	i, maxbad;
	char	unlink_filename[80];

	maxbad = IniGetInt (INI_FILE, (char*)"MaxBadSaveFiles", 10);
	for (i = 1; i <= maxbad; i++) {
		sprintf (unlink_filename, "%s.bad%d", state->base_filename, i);
		_unlink (unlink_filename);
	}
	if (state->num_ordinary_save_files != 99) {
		for (i = 1; i < state->num_ordinary_save_files + state->num_special_save_files; i++) {
			if (i == 1) sprintf (unlink_filename, "%s.bu", state->base_filename);
			else sprintf (unlink_filename, "%s.bu%d", state->base_filename, i);
			_unlink (unlink_filename);
		}
	}
	if (state->base_filename[0] == 'p') {
		sprintf (unlink_filename, "q%s", state->base_filename+1);
		_unlink (unlink_filename);
		sprintf (unlink_filename, "r%s", state->base_filename+1);
		_unlink (unlink_filename);
	}
	sprintf (unlink_filename, "%s.write", state->base_filename);
	_unlink (unlink_filename);
	_unlink (state->base_filename);
}

/* Format a message for writing to the results file and sending to the */
/* server.  We prepend the assignment ID to the message.  If operating */
/* offline, then prepend other information. */

void formatMsgForResultsFile (
	char	*buf,		/* Msg to prepend to, 200 characters max */
	struct work_unit *w)
{
	char	newbuf[2000];

/* Output a USERID/COMPID prefix for result messages 

	if (!USERID[0])
		strcpy (newbuf, buf);
	else if (!COMPID[0])
		sprintf (newbuf, "UID: %s, %s", USERID, buf);
	else
		sprintf (newbuf, "UID: %s/%s, %s", USERID, COMPID, buf);

/* Output the assignment ID too.  That alone should be enough info to */
/* credit the correct userID.  However, we still output the user ID as it */
/* is far more human friendly than an assignment ID. */

	if (w->assignment_uid[0])
		sprintf (newbuf + strlen (newbuf) - 1,
			 ", AID: %s\n", w->assignment_uid);

/* Now truncate the message to 200 characters */

	newbuf[200] = 0;
	strcpy (buf, newbuf);
}

/* Rotate a p-bit giant right shift_count bits.  Used to undo shifted FFT data. */

int rotategp (				/* Return false if failed due to memory allocation error */
	giant	v,			/* Giant to rotate right */
	unsigned long p,		/* Exponent */
	unsigned long shift_count       /* Number of bits to rotate right */
	)
{
	giant	vlo, modulus, mvlo;


/* If rotate count is zero, no work needed */

	if (shift_count >= p)		// only to be careful...
		shift_count = shift_count%p;

	if (shift_count == 0) return (1);
        
/* Convert current iteration to binary */

//	vlo = 		// JP 24/04/20
        vlo = newgiant(2*FFTLEN*sizeof(double)/sizeof(short) + 16);
	if (vlo == NULL) return (0);
//	mvlo = newgiant ((p >> 3) + 5);		// JP 06/05/20
        mvlo = newgiant(2*FFTLEN*sizeof(double)/sizeof(short) + 16);
	if (mvlo == NULL) return (0);
//	modulus = newgiant ((p >> 3) + 5);	// JP 24/04/20
        modulus = newgiant(2*FFTLEN*sizeof(double)/sizeof(short) + 16);
	if (modulus == NULL) return (0);
	itog ((unsigned long)w->k, modulus);					// modulus is 2^p-1 (Mersenne),2^p+1 (Wagstaff or Gaussian Mersenne norm),k*2^p+1 Proth or k*2^p-1 Riesel
	gshiftleft (p, modulus);
	iaddg (w->c, modulus);
	itog (1, vlo);
	gshiftleft (shift_count, vlo);	// inverse of 2^shift_count modulo k*2^p+c
//        invg (modulus,vlo);
        gwypinvg (modulus,vlo);
        gtog (modulus, mvlo);
	if ((w->c>0) && ((w->prp_residue_type==GMN_TYPE)?shift_count&2:shift_count&1))
		subg (vlo, mvlo);					// change the sign if necessary
	else
		gtog (vlo, mvlo);
	mulg (mvlo, v);							// restore unshifted v
	modg (modulus, v);
	free (vlo);
	free (modulus);
	free (mvlo);
	return (1);
}

/* Format the ETA for output to the worker window */

void formatETA (
	double	howlong,		/* how long to complete (in seconds) */
	char	*buf)
{
	double days, hours, minutes, seconds;
	days = floor (howlong / 86400.0);  howlong -= days * 86400.0;
	hours = floor (howlong / 3600.0);  howlong -= hours * 3600.0;
	minutes = floor (howlong / 60.0);  howlong -= minutes * 60.0;
	seconds = floor (howlong);
	if (days >= 3.0)
		sprintf (buf, ", ETA: %dd %02d:%02d", (int) days, (int) hours, (int) minutes);
	else
		sprintf (buf, ", ETA: %02d:%02d:%02d", (int) (days * 24.0 + hours), (int) minutes, (int) seconds);
}

/* Compare two (possibly shifted) gwypnums for equality.  Used in PRP error-checking. */

int areTwoPRPValsEqual (
	unsigned long p,		/* Mersenne exponent (for shifting) */		
	gwypnum	val1,			/* Value #1 */
	unsigned long units_bit1,	/* Shift count #1 */
	gwypnum	val2,			/* Value #2 */
	unsigned long units_bit2	/* Shift count #2 */
)
{
	giant	tmp1, tmp2;
	int	diff, err_code1, err_code2, tmpsize;
        tmpsize = 2*FFTLEN*sizeof(double)/sizeof(short) + 16;
	tmp1 = newgiant(tmpsize);
	err_code1 = gwyptogiant (val1, tmp1) || isZero (tmp1);
	rotategp (tmp1, p, units_bit1);
	tmp2 = newgiant (tmpsize);
	err_code2 = gwyptogiant (val2, tmp2) || isZero (tmp2);
	rotategp (tmp2, p, units_bit2);
	diff = gcompg (tmp1, tmp2);
        free (tmp1);
        free (tmp2);
	return (err_code1 == 0 && err_code2 == 0 && !diff);
}

/* multiplication or squaring  state to be written to / read from save files */

struct program_state {
	gwypnum	x;			/* The current value in our left-to-right exponentiation */
	gwypnum	y;					/* used when some interim value is needed. */
	unsigned long counter;		/* Current "iteration" counter */
	unsigned long units_bit;	/* For shifting FFT data -- allows more robust double-checking */
	unsigned long units_bit2;	/* Keep the shift count for y */
	unsigned long error_count;	/* Count of any errors that have occurred during test */
	unsigned int prp_base;		/* Fermat PRP base, default is 3 */
	int	two_power_opt;			/* TRUE is power of two optimizations are enabled (N adjusted to create long run of squarings) */
	int	residue_type;			/* Type of residue to generate (5 different residue types are supported) */
	int	error_check_type;		/* 0=none, 1=Gerbicz, 2=double-checking */
	int	state;					/* State variable see definitions below */
	gwypnum	alt_x;				/* When doing ERRCHK_DBLCHK, this is the alternate x value */
								/* When doing ERRCHK_GERBICZ, this is the comparison checksum value being calculated */
	unsigned long alt_units_bit;/* When doing ERRCHK_DBLCHK, this is the alternate shift count */
	gwypnum	u0;					/* Saved first value of the Gerbicz checksum function */
	gwypnum	d;					/* Last computed value of the Gerbicz checksum function */
	unsigned long L;			/* Iterations between multiplies in computing Gerbicz checksum */
	unsigned long start_counter;	/* Counter at start of current Gerbicz or double-check block */
	unsigned long next_mul_counter;	/* Counter when next Gerbicz multiply takes place */
	unsigned long end_counter;	/* Counter when current Gerbicz or double-check block ends */
	giant gx;					/* x result of the test in giant form */
	giant gy;					/* y result of the test in giant form */
} ps;

#define ERRCHK_NONE		0	/* No high-reliability error-checking */
#define ERRCHK_GERBICZ	1	/* Gerbicz high-reliability error-checking -- very low overhead */
#define ERRCHK_DBLCHK	2	/* Run test twice, comparing residues along the way. Highly-reliable, very expensive. */

#define STATE_NORMAL				 0	/* Normal left-to-right exponentiation */
#define STATE_DCHK_PASS1			10	/* Do squarings for a while, then rollback counter and switch to pass 2 */
#define STATE_DCHK_PASS2			11	/* Do squarings for a while, compare vals, then switch back to pass 1 */
#define STATE_GERB_START_BLOCK		22	/* Determine how many iters are in the next Gerbicz block */
#define STATE_GERB_MID_BLOCK		23	/* Do squarings for L iterations */
#define STATE_GERB_MID_BLOCK_MULT	24	/* Do checksum multiply after the L-th squaring (except last block) */
#define STATE_GERB_END_BLOCK		25	/* After L^2 squarings, do alt squarings to compute 2nd Gerbicz compare value */
#define STATE_GERB_END_BLOCK_MULT	26	/* After L alt squarings, do one last mul to compute checksum #2 */
#define STATE_GERB_FINAL_MULT		27	/* Do one last mul to compute checksum #1 and then compare checksum values */
/* Write intermediate PRP results to a file */
/* The PRP save file format is: */
/*	u32		magic number  (different for ll, p-1, prp, tf, ecm) */
/*	u32		version number */
/*	double		pct complete */
/*	char(11)	stage */
/*	char(1)		pad */
/*	u32		checksum of following data */
/*	u32		error_count */
/*	u32		iteration counter */
/*	u32		prp base (version number >= 2) */
/*	u32		shift count (version number >= 2) */
/*	u32		power-of-two optimization was used (version number >= 2) */
/*	u32		residue type (version number >= 3) */
/*	u32		error-check type (version number >= 3) */
/*	u32		state (version number >= 3) */
/*	u32		alternate shift count (version number >= 3) */
/*	u32		L - iterations between Gerbicz multiplies (version number >= 3) */
/*	u32		error-checking start counter (version number >= 3) */
/*	u32		error-checking next Gerbicz multiply counter (version number >= 3) */
/*	u32		error-checking end counter (version number >= 3) */
/*	gwypnum		FFT data for x (u32 len, array u32s) */
/*	gwypnum		FFT data for alt_x (u32 len, array u32s) (version number >= 3) */
/*	gwypnum		FFT data for u0 (u32 len, array u32s) (version number >= 3) */
/*	gwypnum		FFT data for d (u32 len, array u32s) (version number >= 3) */

/* Prepare for writing save files */

void writeSaveFileStateInit (
	writeSaveFileState *state,
	char	*filename,
	int	num_special_save_files)
{
	strcpy (state->base_filename, filename);
	state->num_ordinary_save_files = NUM_BACKUP_FILES;
	state->num_special_save_files = num_special_save_files;
	state->special = 0;			/* Init with "no ordinary files are special" */
}

int writeSaveFile (
	writeSaveFileState *write_save_file_state,
	struct work_unit *w,
	struct program_state *ps)
{
	int	fd, i;
	unsigned long sum = 0;

/* Now save to the intermediate file */
	fd = openWriteSaveFile (write_save_file_state);
	if (fd < 0) return (FALSE);

	if (!write_header (fd, PRP_MAGICNUM, PRP_VERSION, w)) goto err;

	if (!write_long298 (fd, ps->error_count, &sum)) goto err;
	if (!write_long298 (fd, ps->counter, &sum)) goto err;
	if (!write_long298 (fd, ps->prp_base, &sum)) goto err;
	if (!write_long298 (fd, ps->units_bit, &sum)) goto err;
	if (!write_long298 (fd, ps->units_bit2, &sum)) goto err;
	if (!write_long298 (fd, ps->two_power_opt, &sum)) goto err;
	if (!write_long298 (fd, ps->residue_type, &sum)) goto err;
	if (!write_long298 (fd, ps->error_check_type, &sum)) goto err;
	if (!write_long298 (fd, ps->state, &sum)) goto err;
	if (!write_long298 (fd, ps->alt_units_bit, &sum)) goto err;
	if (!write_long298 (fd, ps->L, &sum)) goto err;
	if (!write_long298 (fd, ps->start_counter, &sum)) goto err;
	if (!write_long298 (fd, ps->next_mul_counter, &sum)) goto err;
	if (!write_long298 (fd, ps->end_counter, &sum)) goto err;
	if (!write_gwypnumjp (fd, ps->x, &sum)) goto err;
	if (!write_gwypnumjp (fd, ps->y, &sum)) goto err;

	if (ps->state != STATE_NORMAL && ps->state != STATE_GERB_MID_BLOCK && ps->state != STATE_GERB_MID_BLOCK_MULT) {
            if (!write_gwypnumjp (fd, ps->alt_x, &sum)) goto err;
	}

	if (ps->state != STATE_NORMAL && ps->state != STATE_DCHK_PASS1 && ps->state != STATE_DCHK_PASS2 &&
	    ps->state != STATE_GERB_START_BLOCK && ps->state != STATE_GERB_FINAL_MULT) {
            if (!write_gwypnumjp (fd, ps->u0, &sum)) goto err;
	}

	if (ps->state != STATE_NORMAL && ps->state != STATE_DCHK_PASS1 && ps->state != STATE_DCHK_PASS2 &&
	    ps->state != STATE_GERB_START_BLOCK) {
            if (!write_gwypnumjp (fd, ps->d, &sum)) goto err;
	}
	
/* Save the five timers */

	for (i=0; i<5; i++) {
		if (gwyptimers[i+5] != 0.0) {// if the timer was running
			gwypend_timer (i);			// update and save it
			if (! write_double (fd, gwyptimers[i], (long*)&sum)) {
				gwypstart_timer (i);	// and then, restart it, even if write is in error!
				goto err;
			}
			if (! write_double (fd, gwyptimers[i+5], (long*)&sum)) { // save the timer status
				gwypstart_timer (i);	// and then, restart it, even if write is in error!
				goto err;
			}
			gwypstart_timer (i);	// and then, restart it!
		}
		else {
			if (! write_double (fd, gwyptimers[i], (long*)&sum)) goto err;	// save the timer
			if (! write_double (fd, gwyptimers[i+5], (long*)&sum)) goto err;	// save its status
		}
	}
	
	_lseek (fd, CHECKSUM_OFFSET, SEEK_SET);
	if (_write (fd, &sum, sizeof (unsigned long)) != sizeof (unsigned long)) return (FALSE);

	closeWriteSaveFile (write_save_file_state, fd);
	return (TRUE);

/* An error occured.  Delete the current file. */

err:	deleteWriteSaveFile (write_save_file_state, fd);
	return (FALSE);
}

/* Read the data portion of an intermediate PRP results file */

int readPRPSaveFile (
	char	*filename,		/* Save file name */
	struct work_unit *w,		/* Work unit */
	struct program_state *ps)		/* PRP state structure to read and fill in */
{
	int	fd, i, ner = 0;
	unsigned long savefile_prp_base, sum, filesum, version;
	fd = _open (filename, _O_BINARY | _O_RDONLY);
	if (fd <= 0) return (FALSE);
	if (!read_magicnum (fd, PRP_MAGICNUM)) goto err;
        ner=1;
	if (!read_header (fd, &version, w, &filesum)) goto err;
        ner=2;
	if (version == 0 || version > PRP_VERSION) goto err;
        ner=3;

	sum = 0;
	if (!read_long298 (fd, &ps->error_count, &sum)) goto err;
        ner=4;
	if (!read_long298 (fd, &ps->counter, &sum)) goto err;
        ner=5;

	if (version <= 1) {
		savefile_prp_base = 3;
		ps->units_bit = 0;
		ps->two_power_opt = FALSE;
	} else {
		unsigned long savefile_two_power_opt;
		if (!read_long298 (fd, &savefile_prp_base, &sum)) goto err;
        ner=6;
		if (!read_long298 (fd, &ps->units_bit, &sum)) goto err;
        ner=7;
		if (!read_long298 (fd, &ps->units_bit2, &sum)) goto err;
        ner=8;
		if (!read_long298 (fd, &savefile_two_power_opt, &sum)) goto err;
        ner=9;
		ps->two_power_opt = savefile_two_power_opt;
	}
	if (savefile_prp_base != ps->prp_base) goto err;
        ner=10;

	if (version <= 2) {
		ps->state = STATE_NORMAL;
		ps->error_check_type = ERRCHK_NONE;
	} else {
		unsigned long savefile_state, savefile_residue_type, savefile_error_check_type;
		// We might be able to handle some mismatched residue types, for now don't since this should never happen
		if (!read_long298 (fd, &savefile_residue_type, &sum)) goto err;
        ner=11;
		if (savefile_residue_type != (unsigned long)ps->residue_type) goto err;
        ner=12;
		// We could handle some mismatched error check types by looking at the state
		// variable, for now don't since this should never happen
		if (!read_long298 (fd, &savefile_error_check_type, &sum)) goto err;
        ner=13;
		if (savefile_error_check_type != (unsigned long)ps->error_check_type) goto err;
        ner=14;
		if (!read_long298 (fd, &savefile_state, &sum)) goto err;
        ner=15;
		ps->state = savefile_state;
		if (!read_long298 (fd, &ps->alt_units_bit, &sum)) goto err;
        ner=16;
		if (!read_long298 (fd, &ps->L, &sum)) goto err;
        ner=17;
		if (!read_long298 (fd, &ps->start_counter, &sum)) goto err;
        ner=18;
		if (!read_long298 (fd, &ps->next_mul_counter, &sum)) goto err;
        ner=19;
		if (!read_long298 (fd, &ps->end_counter, &sum)) goto err;
        ner=20;
	}

	// In version 3, we did not delay the final multiply in calculation of checksum #1.
	// We must ignore some save files because the version 3 and version 4 states are subtly different.
	if (version == 3 && (ps->state == STATE_GERB_MID_BLOCK_MULT ||
			     ps->state == STATE_GERB_END_BLOCK ||
			     ps->state == STATE_GERB_END_BLOCK_MULT)) goto err;
        ner=21;

	// All PRP states wrote an x value
	if (!read_gwypnumjp (fd, ps->x, &sum)) goto err;
        ner=22;
	if (!read_gwypnumjp (fd, ps->y, &sum)) goto err;
        ner=23;

	// In version 3, we only wrote x to the save file at Gerbicz start block.  In version 4, we write x and the
	// identical alt_x.  There is added error protection by always having at least two gwypnum values in memory.
	if (version == 3 && ps->state == STATE_GERB_START_BLOCK) {
		gwypcopy (ps->x, ps->alt_x);
		ps->alt_units_bit = ps->units_bit;
	}
	else if (ps->state != STATE_NORMAL && ps->state != STATE_GERB_MID_BLOCK && ps->state != STATE_GERB_MID_BLOCK_MULT) {
            if (!read_gwypnumjp (fd, ps->alt_x, &sum)) goto err;
            ner=24;
	}

	// Most PRP Gerbicz states wrote a u0 value
	if (ps->state != STATE_NORMAL && ps->state != STATE_DCHK_PASS1 && ps->state != STATE_DCHK_PASS2 &&
	    ps->state != STATE_GERB_START_BLOCK && ps->state != STATE_GERB_FINAL_MULT) {
            if (!read_gwypnumjp (fd, ps->u0, &sum)) goto err;
            ner=25;
	}

	// Most PRP Gerbicz states wrote a d value
	if (ps->state != STATE_NORMAL && ps->state != STATE_DCHK_PASS1 && ps->state != STATE_DCHK_PASS2 &&
	    ps->state != STATE_GERB_START_BLOCK) {
            if (!read_gwypnumjp (fd, ps->d, &sum)) goto err;
            ner=26;
	}
	
/* Read the five timers and their status, and restart them! */

	for (i=0; i<5; i++) {
		if (! read_double (fd, &gwyptimers[i], (long*)&sum)) goto err;
        ner=27;
		if (! read_double (fd, &gwyptimers[i+5], (long*)&sum)) goto err;
		gwypstart_timer (i);
        ner=28;
	}

	// Validate checksum and return
	
	if (filesum != sum) goto err;
	_close (fd);
	return (TRUE);
err:	_close (fd);
        trace (ner);
	return (FALSE);
}

#define BIT 1
#define ITER 0

void writeresidue (
    gwypnum s,	// The gwypnum data
    giant m,	// The current external modulus
    giant t,	// A temporary giant file	
    char *b,	// The output buffer
    const char *str,// The tested number as a string
    const int bit,// The iteration or bit number
    const int kind  // Bit or Iteration
)
{
    char restr[20];

    gwyptogiant (s, t);
    // The modulo reduction is done here
    modg (m, t);
    // External modulus and gwypnum's one may be different...
    if (abs(t->sign) < 1)
    // make a 32 bit residue correct !!
        sprintf (restr, "%08lX%08lX", (unsigned long)0, (unsigned long)0);
    else if (abs(t->sign) < 2)
        sprintf (restr, "%08lX%08lX", (unsigned long)0, (unsigned long)t->n[0]);
    else
        sprintf (restr, "%08lX%08lX", (unsigned long)t->n[1], (unsigned long)t->n[0]);
    sprintf (b, "%s interim residue %s at %s %d\n", str, restr, kind? "bit" : "iteration", bit);
    if (verbose)
        OutputBoth (b);
    else
        OutputStr (b);
}

void formatresidue (giant tmp, char *resbuf) {
    if (abs(tmp->sign) < 2)
    // make a 64 bit residue correct !!
        sprintf (resbuf, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
    else if (abs(tmp->sign) < 3)
        sprintf (resbuf, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
    else if (abs(tmp->sign) < 4)
        sprintf (resbuf, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
    else
        sprintf (resbuf, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
}

// Compute the number of digits of a large integer.

int gnbdg (giant nb, int digitbase) {
    giant gnbmax;
    int templ;
    if (digitbase==2)
        return (bitlen (nb));
    else {
        templ = (int)(floor((double)bitlen (nb) * log (2.0) /log ((double)digitbase)));
            // default value
        gnbmax = newgiant (abs (nb->sign) + 2);
        itog (digitbase, gnbmax);
        power (gnbmax, templ);
        while (gcompg (gnbmax, nb) < 0) {
            templ++;                                      // adjust the value if it is too small
            ulmulg ((unsigned long) digitbase, gnbmax);  // and adjust the comparand
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
        return (6);					// Invalid numeric string...
    if (prptest)
        retcode = mpz_strongbpsw_prp (n1);		// Strong BPSW PRP test
    else
        retcode = mpz_aprtcle (n1, verbose);
        // Prime test possibly with print out.
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
        return (8);// Could not create the file...

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
        &aprprinfo)) {
            errcode = GetLastError ();	// Echec...
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
        return (-1);    // Number too large...
    if (N->sign == 1) {
        // Trial divisions test is sufficient for this small number...
        return (isPrime (N->n[0]) ? 12 : 10);
    }
    gtoc (N, greatbuf, 10000);
    return (aprcltest(prptest, verbose));
}

int setupok (int errcode, giant bignumber, char *string, int *resultat)    // Test if the call to gwypsetup is successful
{
    char buff[256];
    int resaprcl;

    if (!errcode)           // The setup is OK
        return TRUE;
    else if (errcode == 1) {// Number too small...
        
/* Format the string representation of the test number */

	if (w->known_factors == NULL) {
		strcpy (string_rep, string); // 02/07/20
		string_rep_truncated = FALSE;
	} else {
		if (strchr (string, '^') == NULL)
			strcpy (string_rep, string);
		else
			sprintf (string_rep, "(%s)", string);
		if (strlen (w->known_factors) < 40) {
			char	*p;
			strcat (string_rep, "/");
			strcat (string_rep, w->known_factors);
			while ((p = strchr (string_rep, ',')) != NULL) *p = '/';
			string_rep_truncated = FALSE;
		} else {
			strcat (string_rep, "/known_factors");
			string_rep_truncated = TRUE;
                }
		if (string_rep_truncated) {
			char	*bigbuf;
			bigbuf = (char *) malloc (strlen (w->known_factors) + 100);
			if (bigbuf != NULL) {
				sprintf (bigbuf, "Known factors used for the test were: %s\n", w->known_factors);
				OutputStr (bigbuf);
                                free (bigbuf);
                        }
		}
	}
        nbdg = gnbdg (bignumber, 10);
        gwypstart_timer(1);
        if (nbdg > maxaprcl)
            // Make only a Strong BPSW PRP test
            resaprcl = gaprcltest (bignumber, 1, 0);
        else if (debug)
            // Primality test while showing progress 
            resaprcl = gaprcltest (bignumber, 0, 2);
        else// Primality test silently done
            resaprcl = gaprcltest (bignumber, 0, 0);
        gwypend_timer (1);
        if (resaprcl == 10) {
            sprintf (buff,"%s is not prime. (Trial divisions)", string_rep);
        }
        else if (resaprcl == 12)
            sprintf (buff,"%s is prime! (%d decimal digits, Trial divisions)", string_rep, nbdg);
        else if (resaprcl == 0)
            sprintf (buff,"%s is not prime. (APRCL test) ", string_rep);
        else if (resaprcl == 1)
            sprintf (buff,"%s is a probable BPSW prime! (%d decimal digits, APRCL test) ", string_rep, nbdg);
        else if (resaprcl == 2)
            sprintf (buff,"%s is prime! (%d decimal digits, APRCL test)", string_rep, nbdg);
        else if (resaprcl == 6)
            sprintf (buff,"Invalid numerical string in %s\n", string_rep);
        else if (resaprcl == 7)
            sprintf (buff,"APRCL error while testing %s...\n", string_rep);
        else {
            if (resaprcl == 9)
                sprintf (buff, "APRCL primality test not available for %s\n", string_rep);
            else
                sprintf (buff,"Unexpected return value : %d, while APRCL testing %s...\n", resaprcl, string_rep);
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
        return FALSE;
    }
    else {
        sprintf (buff, "%s : Fatal error at setup, code = %d", string, errcode);
        strcpy(buff+strlen(buff), "\n");
        OutputBoth (buff);
        return FALSE;
    }
}

/* Compare the final giant in a PRP run.  Different PRP residue types check for different final values */

int isPRPg (
	giant	v,			/* Final result of PRP powering */
	giant	N,			/* Number we are PRPing, (k*b^n+c)/known_factors */
	struct work_unit *w,		/* Number being tested */
	unsigned int prp_base,		/* PRP base */
	int	prp_residue_type)	/* Type of PRP test performed */
{

	int	result = FALSE;
	giant kbnc, known_factors, reduced_v, compare_val;
        
/* Standard Fermat PRP, test for one */

            if (prp_residue_type == PRP_TYPE_FERMAT) return (isone (v));

/* SPRP test is PRP if result is one or minus one */

            if (prp_residue_type == PRP_TYPE_SPRP || prp_residue_type == PROTH_TYPE) {
		if (prp_residue_type == PRP_TYPE_SPRP && isone (v))
                    return (TRUE);
                iaddg (1, v);
                result = (gcompg (v, N) == 0);
                iaddg (-1, v);
		return (result);
            }

/* Handle the cofactor case.  We calculated v = a^(N*KF-1) mod (N*KF).  We have a PRP if (v mod N) = (a^(KF-1)) mod N */

            if (prp_residue_type == PRP_TYPE_COFACTOR) {
		// Calculate k*b^n+c
		if (w->k > 0.0) {
			kbnc = newgiant ((w->n + (unsigned long)_log2(w->k)) + 5);
			known_factors = newgiant ((w->n + (unsigned long)_log2(w->k)) + 5);
			reduced_v = newgiant ((w->n + (unsigned long)_log2(w->k)) + 5);
			compare_val = newgiant ((w->n + (unsigned long)_log2(w->k)) + 5);
			itog ((unsigned long)w->k, kbnc);
		}
		else {
			kbnc = newgiant ((w->n + bitlen(gk)) + 5);
			known_factors = newgiant ((w->n + bitlen(gk)) + 5);
			reduced_v = newgiant ((w->n + bitlen(gk)) + 5);
			compare_val = newgiant ((w->n + bitlen(gk)) + 5);
			gtog (gk, kbnc);
		}
		gshiftleft (w->n, kbnc);
		iaddg (w->c, kbnc);
		gtog (kbnc, known_factors);
		divg (N, known_factors);
		iaddg (-1, known_factors);
		itog (prp_base, compare_val);
		powermodg (compare_val, known_factors, kbnc);
		modg (N, compare_val);
		modg (kbnc, v);
		gtog (v, reduced_v);
		modg (N, reduced_v);
		result = (gcompg (reduced_v, compare_val) == 0);
		free (kbnc);
		free (known_factors);
		free (reduced_v);
		free (compare_val);
            }

	return (result);
}

/* Mul giant by a power of the PRP base.  Used in optimizing PRPs for k*2^n+c numbers. */

void basemulg (                         // 15/04/21 using now giant <-> GMP conversions 
	giant	v,			/* Giant to multiply by base^power */
	struct work_unit *w,		/* Number being tested */
	unsigned int prp_base,		/* PRP base */
	int	power)			/* Desired power of the PRP base */
{
	mpz_t	modulus, prp_base_power, tmp;

/* If power is zero, then multiply by base^0 is a no-op */

	if (power == 0) return;

/* Generate the modulus (k*b^n+c), b is known to be 2 */

//	mpz_init_set_d (modulus, w->k);
        mpz_init (modulus);
        gtompz (gk, modulus);           // 17/04/21 k may be a large integer...
	mpz_mul_2exp (modulus, modulus, w->n);
	mpz_add_si (modulus, modulus, w->c);

/* Calculate prp_base^power mod k*b^n+c */

	mpz_init_set_si (tmp, power);
	mpz_init_set_ui (prp_base_power, prp_base);
	mpz_powm (prp_base_power, prp_base_power, tmp, modulus);

/* Copy the giant value to tmp.  Multiply by prp_base_power to get the final result */

	gtompz (v, tmp);
	mpz_mul (tmp, tmp, prp_base_power);
	mpz_mod (tmp, tmp, modulus);
	mpztog (tmp, v);

/* Cleanup and return */

	mpz_clear (tmp);
	mpz_clear (prp_base_power);
	mpz_clear (modulus);
}

/* Test if M divides a^(N-1) - 1 -- gwypsetup has already been called. */

int isexpdiv (
    long a,
    giant N,
    giant M,
    int	*res)
{
    unsigned long bit, iters;
    gwypnum x;
    giant   tmp;
    char    filename[20], buf[sgkbufsize+256],
        fft_desc[256], oldres64[17];
    long    write_time = DISK_WRITE_TIME * 60;
    int	echk, saving, stopping;
    time_t  start_time, current_time;
    double  reallyminerr = 1.0;
    double  reallymaxerr = 0.0;


/* Init, subtract 1 from N to compute a^(N-1) mod M */

    iaddg (-1, N);
    Nlen = bitlen (N);
    nbdg = gnbdg (N, 10);
    // Compute the number of decimal digits of the tested number.
    *res = TRUE;/* Assume the residue is one */

/* Init filename */

    tempFileName (filename, 'z', N);

/* Allocate memory */

    x = gwypalloc ();
    nbllr_mallocs++;
    tmp = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);
    nbllr_mallocs++;

/* Optionally resume from save file and output a message */
/* indicating we are resuming a test */

    if (fileExists (filename) && readFromFile (filename, &bit, x, NULL)) {
        char    fmt_mask[80];
        double  pct;
        pct = trunc_percent (bit * 100.0 / Nlen);
        sprintf (fmt_mask,
        "Resuming divisibility test of %%d^(N-1)-1 at bit %%ld [%%.%df%%%%]\n", PRECISION);
        sprintf (buf, fmt_mask, a, bit, pct);
//		OutputStr (buf);
//		if (verbose)
//			writeResults (buf);
    }

/* Otherwise, output a message indicating we are starting test */

    else {
        gwypclear_timers ();// Init. the timers
        sprintf (buf, "Starting divisibility test of %ld^(N-1)-1\n", a);
//		OutputStr (buf);
//		if (verbose)
//			writeResults (buf);
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
//	OutputStr (buf);
//	LineFeed ();
//	if (verbose) {
#if !defined(WIN32) 
//		strcat (buf, "\n");
#endif
//		writeResults (buf);
//	}
    ReplaceableLine (1);
        /* Remember where replaceable line is */

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
        }
        else {
            gwypsetnormroutine (0, echk, 0);
        }

        if (/*(bit+25 < Nlen) && (bit > 25) && */((bit != lasterr_point) || !maxerr_recovery_mode[6])) {
            if (cufftonly)
                gwypsquare (x);
            else if (zp || generic)
                cuda_gwypsquare (x, 3);
            else if(bit==1 || it==0) {
                cuda_gwypsquare (x,1);
                it=1;
            }
            else  if(bit != (lasterr_point-1)&&(bit+25 < Nlen) && (bit > 25)) 
                cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
            else
                cuda_gwypsquare (x,2);
            care = FALSE;
        }
        else {
            gwypsquare_carefully (x);
            care = TRUE;
        }

        CHECK_IF_ANY_ERROR (x, (bit), Nlen, 6);

/* That iteration succeeded, bump counters */

        if (bit == lasterr_point)
            saving = 1;
// Be sure to restart after this recovery iteration!
        bit++;
        iters++;

/* Print a message every so often */

        if (bit % ITER_OUTPUT == 0) {
            char    fmt_mask[80];
            double  pct;
            pct = trunc_percent (bit * 100.0 /  
                Nlen);
            sprintf (fmt_mask, "%%.%df%%%% of %%ld",
                PRECISION);
            sprintf (buf, fmt_mask, pct, Nlen);
            title (buf);
            ReplaceableLine (2);/* Replace line */
            sprintf (fmt_mask,
                "%%d^(N-1)-1, bit: %%ld / %%ld [%%.%df%%%%]", PRECISION);
            sprintf (buf, fmt_mask, a, bit, Nlen, 
                pct);
            OutputStr (buf);
            if (ERRCHK && bit > 30) {
                OutputStr ((char*)".  Round off: ");
                sprintf(buf,"%10.10f",reallyminerr);
                OutputStr (buf);
                sprintf (buf," to %10.10f",
                    reallymaxerr);
                OutputStr (buf);
            }
            gwypend_timer (0);
            if (CUMULATIVE_TIMING) {
                OutputStr ((char*)".  Time thusfar: ");
            }
            else {
                OutputStr ((char*)".  Time per bit: ");
                gwypdivide_timer (0, iters);
                iters = 0;
            }
            gwypprint_timer (0, TIMER_NL | TIMER_OPT_CLR);
            gwypstart_timer (0);
        }

/* Print a results file message every so often */

        if (bit % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
            sprintf (buf, "Bit %ld / %ld\n", bit,
                Nlen);
            writeResults (buf);
        }

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

        if (saving || stopping) {
            write_time = DISK_WRITE_TIME * 60;
            saving = FALSE;
            if (! writeToFile (filename, bit, x, NULL)) {
                sprintf (buf, WRITEFILEERR,
                    filename);
                OutputBoth (buf);
                if (write_time > 600)
                    write_time = 600;
            }	
            time (&start_time);

/* If an escape key was hit, write out the results and return */

            if (stopping) {
                iaddg (1, N);// Restore the modulus
                gwypfree (tmp);
                nbllr_frees++;
                gwypfree (x);
                nbllr_frees++;
//		gwypdone ();
                *res = FALSE;
                return (FALSE);
            }
        }

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next iteration so that we can compare our */
/* residues to programs that start counter at zero or one. */

        if (interimResidues && bit % interimResidues < 2) {
            gwyptogiant (x, tmp);
            // The modulo reduction is done here
            if (abs(tmp->sign) < 2)
            // make a 64 bit residue correct !!
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
            char    interimfile[20];
            sprintf (interimfile, "%.8s.%03lu",
                filename, bit / interimFiles);
            if (! writeToFile (interimfile, bit, x,
                NULL)) {
                sprintf (buf, WRITEFILEERR,    
                    interimfile);
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
        if (abs(tmp->sign) < 2)
            // make a 64 bit residue correct !!
            sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
        else if (abs(tmp->sign) < 3)
            sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
        else if (abs(tmp->sign) < 4)
            sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
        else
            sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
        smulg ((unsigned short)a, tmp);
        modg (N, tmp);
        iaddg (-a, tmp);
        if (abs(tmp->sign) < 2)
            // make a 64 bit residue correct !!
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
    nbllr_frees++;
    gwypfree (x);
    nbllr_frees++;

/* Output the final timings */

    gwypend_timer (1);
    sprintf (buf+strlen(buf)-1, "  Time: ");
    ReplaceableLine (2);	/* Replace line */
    gwypwrite_timer (buf+strlen(buf), 1, TIMER_CLR | TIMER_NL); 
//  if (verbose)
//	OutputBoth (buf);
//  else
//	OutputStr (buf);

/* Cleanup and return */

    iaddg (1, N);   // Restore the modulus
    Nlen = bitlen (N);
//	gwypdone ();
    _unlink (filename);
    lasterr_point = 0;
    return (TRUE);

/* An error occured, sleep, then try restarting at last save point. */

error:
    iaddg (1, N);   // Restore the value of N
    Nlen = bitlen (N);
    gwypfree (tmp);
    nbllr_frees++;
    gwypfree (x);
    nbllr_frees++;
//	gwypdone ();
    *res = FALSE;   // To avoid credit mesage...

    if (abonroundoff && MAXERR > maxroundoff) {
            // Abort...
        sprintf (buf,ERRMSG7,a);
        aborted = TRUE;
        OutputBoth (buf);
//      gwypdone ();
        _unlink (filename);
        if(IniGetInt(INI_FILE, (char*)"StopOnAbort", 0)) {
            IniWriteInt (INI_FILE, (char*)"PgenLine", IniGetInt(INI_FILE, (char*)"PgenLine", 0) + 1);
                // Point on the next line
            return (FALSE);
        }
        else
            return (TRUE);
    }

//	gwypdone ();

/* Output a message saying we are restarting */

    if (sleep5)
        OutputBoth (ERRMSG2);
    OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

    if (sleep5 && ! SleepFive ()) {
        return (FALSE);
    }

/* Restart */

    if (will_try_larger_fft) {
        IniWriteInt(INI_FILE, (char*)"FFT_Increment", nbfftinc =  IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
        if (nbfftinc == maxfftinc)
            abonroundoff = TRUE;
            // Don't accept any more Roundoff error.
        _unlink (filename);
    }
    return (-1);
}

/* This code is adapted from George Woltman's prp one in Version 29.8 of Prime95 */
/* Thanks to George's work, it was not difficult to implement it in LLR such as  */
/* the Gerbicz error checking is now available by default for all base two tests */


int GerbiczTest (
	unsigned long a,
	int	*res,
	char *str,
	struct work_unit *w)
{
	unsigned long explen, final_counter, iters, errchk, bit, init_units_bit = 0;
        //	unsigned long restart_error_count = 0;	/* On a restart, use this error count rather than the one from a save file */
	int	have_res2048;
	unsigned long final_counter_y, max_counter, gerb_counter, loopshift = 0, last_counter = 0xFFFFFFFF;	/* Iteration of last error */
	long	/*stop_counter = 0,*/ restart_counter = -1;	/* On a restart, this specifies how far back to rollback save files */
	giant	exp, tmp, tmp2, tmp3; 
//	struct	program_state ps;
	int	interim_counter_off_one, interim_mul, PositiveResult;
	int	first_iter_msg, echk, gpuflag = 3, iters_left = 0; 
	double	inverse_explen;
	double	reallyminerr = 1.0; 
	double	reallymaxerr = 0.0; 
	double	allowable_maxerr, output_frequency, output_title_frequency;
//	char	string_rep[80];
	char	filename[20], buf[sgkbufsize+256], 
		fft_desc[256], res2048[513]; 
	long	write_time = DISK_WRITE_TIME * 60; 
//	int	string_rep_truncated;
	int	error_count_messages/*, zeroed_d = FALSE*/;
//	readSaveFileState read_save_file_state; /* Manage savefile names during reading */
	writeSaveFileState write_save_file_state; /* Manage savefile names during writing */
	
/* Init program state */
        memset (&ps, 0, sizeof (ps));
        if (w->prp_base)
		ps.prp_base = w->prp_base;
	else
		ps.prp_base = a;

	if (w->prp_residue_type)
		ps.residue_type = w->prp_residue_type;
	else
		ps.residue_type = IniGetInt (INI_FILE, (char*)"PRPResidueType", PRP_TYPE_COFACTOR);
	if (w->known_factors == NULL && ps.residue_type == PRP_TYPE_COFACTOR)
		ps.residue_type = PRP_TYPE_FERMAT;

/* Set flag if we will perform power-of-two optimizations.  These optimizations reduce the number of mul-by-small constants */
/* by computing a^(k*2^n) which gives us a long run of simple squarings.  These squarings let us do Gerbicz error checking. */

	ps.two_power_opt = !IniGetInt (INI_FILE, (char*)"PRPStraightForward", 0) && (w->b == 2) &&
			    ((w->known_factors == NULL) || (ps.residue_type == PRP_TYPE_COFACTOR));

	interim_counter_off_one = (ps.two_power_opt && w->k == 1.0 && w->c < 0);
	interim_mul = (ps.two_power_opt && w->c < 0);

/* Flag the PRP tests that require multiplying the final a^exp to account for c */

	if (ps.two_power_opt && (ps.residue_type == PRP_TYPE_FERMAT || ps.residue_type == PRP_TYPE_COFACTOR))
		mul_final = w->c - 1;
	else if (ps.two_power_opt && ps.residue_type == PRP_TYPE_SPRP)
		mul_final = (w->c - 1) / 2;
	else
		mul_final = 0;

/* Determine what highly-reliable error-checking will be done (if any) */

	errchk = IniGetInt (INI_FILE, (char*)"ErrorChecking", 1);
	if (errchk == 0) ps.error_check_type = ERRCHK_NONE;
	if (errchk == 1) ps.error_check_type = ps.two_power_opt ? ERRCHK_GERBICZ : ERRCHK_NONE;
	if (errchk == 2) ps.error_check_type = ps.two_power_opt ? ERRCHK_GERBICZ : ERRCHK_DBLCHK;
	if (errchk == 3) ps.error_check_type = ERRCHK_DBLCHK;

/* Calculate the exponent we will use to do our left-to-right binary exponentiation */

	if (ps.residue_type == GMN_TYPE) {
		exp = newgiant (4*FFTLEN*sizeof(double)/sizeof(short) + 16);
		tmp = newgiant (4*FFTLEN*sizeof(double)/sizeof(short) + 16);
		tmp2 = newgiant (4*FFTLEN*sizeof(double)/sizeof(short) + 16);
	}
	else {
		exp = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);
		tmp = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);
		tmp2 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);
	}
	nbllr_mallocs+=3;

/* As a small optimization, base 2 numbers are computed as a^(k*2^n) or a^(k*2^(n-1)) mod N with the final result */
/* multiplied by a^(c-1).  This eliminates tons of mul-by-consts at the expense of lots of bookkeepping headaches */
/* and one squaring if k=1 and c<0. */

	if (ps.two_power_opt) {
		int	gerbicz_squarings;
		if (ps.residue_type == PRP_TYPE_FERMAT ||
		    ps.residue_type == PRP_TYPE_FERMAT_VAR ||
		    ps.residue_type == PRP_TYPE_COFACTOR)
			gerbicz_squarings = w->n;
		else
			gerbicz_squarings = w->n - 1;
		itog (2, exp);
		power (exp, gerbicz_squarings);
		if (w->k != 1.0)
                    mulg (gk, exp);
	}

/* PRP co-factor test.  Do a PRP test on k*b^n+c rather than (k*b^n+c)/known_factors. */

	else if (ps.residue_type == PRP_TYPE_COFACTOR) {
		itog (w->b, exp);
		power (exp, w->n);
		if (w->k != 1.0)
                    mulg (gk, exp);
		iaddg (w->c, exp);
		iaddg (-1, exp);
	}

/* Standard PRP test.  Subtract 1 from N to compute a^(N-1) mod N */

	else {
		gtog (N, exp);
		iaddg (-1, exp);
	}

/* Get the exact bit length of the binary exponent.  We will perform bitlen(exp)-1 squarings for the PRP test. */


        explen = Nlen = bitlen (exp);
	final_counter = explen - 1;
	final_counter_y = (final_counter/2)-1;	// Warning : ps.counter starts from 0 !
	if (ps.error_check_type == ERRCHK_GERBICZ) // 14/04/21
            gerb_counter = final_counter-1049;  // max of L + min of iters_left 13/04/21
        else
            gerb_counter = final_counter;

/* Allocate memory */

	ps.x = gwypalloc ();
        nbllr_mallocs++;
	ps.y = gwypalloc ();
        nbllr_mallocs++;
	if (ps.error_check_type == ERRCHK_GERBICZ || ps.error_check_type == ERRCHK_DBLCHK) {
            ps.alt_x = gwypalloc ();
            nbllr_mallocs++;
        }
	if (ps.error_check_type == ERRCHK_GERBICZ) {
            ps.u0 = gwypalloc ();
            nbllr_mallocs++;
            ps.d = gwypalloc ();
            nbllr_mallocs++;
	}

/* Set the proper starting value and state if no save file was present */

	if (ps.counter == 0) {
            if (ps.residue_type == GMN_TYPE)
                tmp3 = newgiant (4*FFTLEN*sizeof(double)/sizeof(short) + 16);
            else
                tmp3 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);
            nbllr_mallocs++;
            itog (ps.prp_base, tmp3);
            /* For base==2 numbers, k==1 and abs(c)==1 we support FFT data shifting */
            if (w->b==2 && (w->k==1.0) && (abs(w->c)==1) && w->n>1000 && IniGetInt
                (INI_FILE, (char*)"Shifting", 1)) {
                // Generate a random initial shift count
                srand ((unsigned) time (NULL));
                init_units_bit = (rand () << 16) + rand ();
                // Let user override random initial shift count
                init_units_bit = IniGetInt (INI_FILE, (char*)"InitialShiftCount", init_units_bit);
                // Initial shift count can't be larger than n-64 (the -64 avoids wraparound in setting intial value)
                init_units_bit = init_units_bit % (maxbitsinfftlen/2 - 64);
                gshiftleft (init_units_bit, tmp3);
                ps.units_bit = ps.alt_units_bit = init_units_bit;
            } else {
                ps.units_bit = ps.alt_units_bit = init_units_bit = 0;
            }
            gianttogwyp (tmp3, ps.x);
            free (tmp3);
            nbllr_frees++;
            gwypcopy (ps.x, ps.y);	// temporary...
            ps.units_bit2 = 0;

/* The easy state case is no high-reliability error-checking */

		if (ps.error_check_type == ERRCHK_NONE) {
			ps.state = STATE_NORMAL;
			ps.start_counter = 0;			// Value not used
			ps.end_counter = final_counter;	// JP		// Value not used
			ps.alt_units_bit = ps.units_bit; // JP
		}

/* The next easiest case is double-the-work error-checking comparing residues at specified intervals */

		else if (ps.error_check_type == ERRCHK_DBLCHK) {
			ps.state = STATE_DCHK_PASS1;
			ps.start_counter = 0;
			ps.end_counter = IniGetInt (INI_FILE, (char*)"DoublecheckCompareInterval", 100000);
			if (ps.end_counter > final_counter) ps.end_counter = final_counter;
			if (ps.units_bit == 0) {
				gwypcopy (ps.x, ps.alt_x);
				ps.alt_units_bit = 0;
			} else {
				gwypadd3 (ps.x, ps.x, ps.alt_x);
				ps.alt_units_bit = ps.units_bit + 1;
				if (ps.alt_units_bit > ((ps.residue_type == GMN_TYPE)?2*(w->n):w->n)) {
					ps.alt_units_bit -= (ps.residue_type == GMN_TYPE)?2*(w->n):w->n;
				}
			}
		}

/* The final case of high-reliability error-checking is Gerbicz error checking, described below. */
/* See http://mersenneforum.org/showthread.php?t=22471 for background information. */
/* In a nutshell, if PRPing k*2^n+c we calculate (prp_base^k)^(2^n) which ends with a long string of squarings. */
/* Let u[0] = (prp_base^k), our n squarings are defined as:
	u[i]=u[0]^(2^i) mod mp, for i=0..n
   We define a "checksum" function below, which is updated every L-th squaring:
	d[t]=prod(i=0,t,u[L*i]) mod mp
   The key idea is that the checksum function can be calculated two different ways, thus the two can be compared to detect an error:
	d[t]=d[t-1]*u[L*t] mod mp		(checksum #1)
	d[t]=u[0]*d[t-1]^(2^L) mod mp		(checksum #2)
   The larger L we choose, the lower the error-checking overhead cost.  For L of 1000, we catch errors every 1 million iterations
   with an overhead of just 0.2%. */
/* For extra protection, we always keep two copies of the gwypnum in memory.  If either one "goes bad" error checking will catch this. */

		else if (ps.error_check_type == ERRCHK_GERBICZ) {
			// Both STATE_DCHK_PASS1 and STATE_GERB_START_BLOCK expect alt_x to be a copy of x
			gwypcopy (ps.x, ps.alt_x);
			ps.alt_units_bit = ps.units_bit;
			// We first compute (prp_base^k) by double-checking 
			if (w->k != 1.0) {
                            ps.state = STATE_DCHK_PASS1;
                            ps.start_counter = 0;
                            if (w->k != 0.0)
                                ps.end_counter = (int)ceil(_log2(w->k));
                            else
                                ps.end_counter = bitlen (gk);
			} else {
                            ps.state = STATE_GERB_START_BLOCK;
			}
		}
	}              // end ps.counter == 0
	
/* Init vars for Test/Status and CommunicateWithServer */

	if (ps.residue_type == PROTH_TYPE || ps.residue_type == GMN_TYPE)
		strcpy (w->stage, "Proth");
	else
		strcpy (w->stage, "PRP");
	inverse_explen = 1.0 / (double) final_counter;
	w->pct_complete = (double) ps.counter * inverse_explen;
	calc_output_frequencies (&output_frequency, &output_title_frequency);

/* If we are near the maximum exponent this fft length can test, then we */
/* will error check all iterations */

//	near_fft_limit = gwypnear_fft_limit (pcfftlim);


/* Figure out the maximum round-off error we will allow.  By default this is 27/64 when near the FFT limit and 26/64 otherwise. */
/* We've found that this default catches errors without raising too many spurious error messages.  We let the user override */
/* this default for user "Never Odd Or Even" who tests exponents well beyond an FFT's limit.  He does his error checking by */
/* running the first-test and double-check simultaneously. */

	allowable_maxerr = IniGetFloat (INI_FILE, "MaxRoundoffError", (float) (/*near_fft_limit*/0 ? 0.421875 : 0.40625));

/* Get setting for verbosity of hardware error messages.  Force output of "confidence is excellent" when error checking. */

	error_count_messages = IniGetInt (INI_FILE, (char*)"ErrorCountMessages", 3);
	if (ps.error_check_type != ERRCHK_NONE) error_count_messages |= 0x8000;

/* Init filename */

	tempFileName (filename, 'z', N);
        if (!quotient)
            strcpy (string_rep, str);
        
/* Init the write save file state.  This remembers which save files are Gerbicz-checked.  Do this initialization */
/* before the restart for roundoff errors so that error recovery does not destroy thw write save file state. */

	writeSaveFileStateInit (&write_save_file_state, filename, NUM_JACOBI_BACKUP_FILES);

/* Optionally resume from save file and output a message */
/* indicating we are resuming a test */
	if (fileExists (filename) && readPRPSaveFile (filename, w, &ps)) {
		char	fmt_mask[80];
		double	pct;
                restarting = TRUE;
//                stop_counter = ps.counter;
		gwypstart_timer (2);
		if (restart_counter < 0 || ps.counter <= (unsigned long) restart_counter)
			first_iter_msg = TRUE;
		bit = ps.counter+1;
		pct = trunc_percent (bit * 100.0 / Nlen);
		if (ps.residue_type == PROTH_TYPE || ps.residue_type == GMN_TYPE)
			sprintf (fmt_mask,
				"Resuming Proth prime test of %%s at bit %%ld [%%.%df%%%%]\n",
				PRECISION);
		else if (ps.residue_type == PRP_TYPE_SPRP)
			sprintf (fmt_mask,
				"Resuming Strong PRP test of %%s at bit %%ld [%%.%df%%%%]\n",
				PRECISION);
		else
			sprintf (fmt_mask,
				"Resuming Fermat PRP test of %%s at bit %%ld [%%.%df%%%%]\n",
				PRECISION);
		sprintf (buf, fmt_mask, string_rep, bit, pct);
		OutputStr (buf);
		if (verbose || restarting)
			writeResults (buf);
	}

/* Otherwise, output a message indicating we are starting test */

	else {
		gwypclear_timers ();	// Make all timers clean...
		gwypstart_timer (2);
		if (ps.residue_type == PROTH_TYPE || ps.residue_type == GMN_TYPE)
			if (showdigits)
				sprintf (buf, "Starting Proth prime test of %s (%d decimal digits)\n", string_rep, nbdg);
			else
				sprintf (buf, "Starting Proth prime test of %s\n", string_rep);
		else if (ps.residue_type == PRP_TYPE_SPRP)
				sprintf (buf, "Starting Strong PRP test of %s\n", string_rep);
		else
				sprintf (buf, "Starting Fermat PRP test of %s\n", string_rep);
		OutputStr (buf);
		if (verbose || restarting)
			writeResults (buf);
		ps.counter = 0;
		ps.error_count = 0;
		first_iter_msg = TRUE;
		bit = 1;
	}

	gwypstart_timer (1);	// Start loop timer

/* Output a message about the FFT length and the Proth base. */

//	topology_print_children (hwloc_get_root_obj (hwloc_topology), 0);
	gwypfft_description (fft_desc);
	sprintf (buf, "%s, a = %lu\n", fft_desc, a);

	OutputStr (buf);
        
	if (verbose || restarting)
		writeResults (buf);

	ReplaceableLine (1);	/* Remember where replaceable line is */

/* Get setting for verbosity of hardware error messages.  Force output of "confidence is excellent" when error checking. 

	error_count_messages = IniGetInt (INI_FILE, (char*)"ErrorCountMessages", 3);
	if (ps.error_check_type != ERRCHK_NONE) error_count_messages |= 0x8000;

/* Do the Proth or PRP test */
	gwypsetmulbyconst (a);
	iters = 0;
	if (ps.residue_type == GMN_TYPE) {
		loopshift = (ps.counter >= final_counter_y) ? final_counter_y : 0;
		max_counter = (ps.counter >= final_counter_y) ? final_counter : final_counter_y;
	}
	else {
		max_counter = final_counter;
		loopshift = 0;
	}
        
        it = 0; // cuda
        
	while (ps.counter < final_counter) {
		gwypnum	x;  /* Pointer to number to square */
		unsigned long *units_bit;
                            /* Pointer to units_bit to update */
		int	saving, saving_highly_reliable, /*sending_residue, */interim_residue, interim_file, stop_reason;
		int	actual_frequency;

/* If this is the first iteration of a Gerbicz error-checking block, then */
/* determine "L" -- the number of squarings between each Gerbicz multiplication */
/* We end this Gerbicz block after L^2 iterations.  */
/* If there aren't many iterations left, revert to simple double-checking. */

		if (ps.state == STATE_GERB_START_BLOCK) {
 
                    	iters_left = final_counter - ps.counter;
			if (iters_left < 49) {
				ps.state = STATE_DCHK_PASS1;
				ps.start_counter = ps.counter;
				ps.end_counter = final_counter;
				if (ps.alt_units_bit) {
					gwypadd3 (ps.alt_x, ps.alt_x, ps.alt_x);
					ps.alt_units_bit = ps.alt_units_bit + 1;
					if (ps.alt_units_bit > ((ps.residue_type == GMN_TYPE)?2*(w->n):w->n)) {
						ps.alt_units_bit -= (ps.residue_type == GMN_TYPE)?2*(w->n):w->n;
					}
				}
			} else {
				int	gerbicz_block_size;
				double	adjustment;
				ps.state = STATE_GERB_MID_BLOCK;
				adjustment = IniGetFloat (INI_FILE, "PRPGerbiczCompareIntervalAdj", 1.0);
				if (adjustment < 0.001 || adjustment > 1.0) adjustment = 0.5;
				gerbicz_block_size = (int) (adjustment * IniGetInt (INI_FILE, (char*)"GerbiczCompareInterval", 1000000));
				if (gerbicz_block_size < 25) gerbicz_block_size = 25;
				if (gerbicz_block_size > iters_left) gerbicz_block_size = iters_left;
				ps.L = (unsigned long) sqrt ((double) gerbicz_block_size);
				ps.start_counter = ps.counter;
				ps.next_mul_counter = ps.counter + ps.L;
				ps.end_counter = ps.counter + ps.L * ps.L;
				gwypswap (ps.alt_x, ps.u0);		// Set u0 to a copy of x
//				gwypcopy (ps.alt_x, ps.u0);		// Set u0 to a copy of x
				gwypcopy (ps.x, ps.d);	// Set d[0] to a copy of x
				if (IniGetInt (INI_FILE, (char*)"GerbiczVerbosity", 1) > 1) {
					sprintf (buf, "Start Gerbicz block of size %ld at iteration %ld.\n", ps.L * ps.L, ps.start_counter+1);
					OutputBoth (buf);
				}
			}
		}

/* Save if we are stopping, right after we pass an errored iteration, several iterations before retesting */
/* an errored iteration so that we don't have to backtrack very far to do a gwsquare_carefully iteration */
/* (we don't do the iteration immediately before because a save operation may change the FFT data and make */
/* the error non-reproducible), and finally save if the save file timer has gone off. */
		stop_reason = stopCheck ();
		saving = stop_reason || (ps.counter == last_counter-8) || (ps.counter == last_counter);
		saving_highly_reliable = FALSE;

/* Round off error check the first and last 50 iterations, before writing a save file, near an FFT size's limit, */
/* or check every iteration option is set, and every 128th iteration. */

		echk = ERRCHK || ps.counter < 50 || ps.counter >= final_counter-50 || saving ||
		       (ps.error_check_type == ERRCHK_NONE && (/*near_fft_limit*/0 || ((ps.counter & 127) == 0)));
		gwyp_clear_maxerr ();

/* Check if we should output residue to screen, or create an intermediate save file */

/*		sending_residue = (ps.state == STATE_NORMAL || ps.state == STATE_DCHK_PASS1 || ps.state == STATE_GERB_MID_BLOCK) &&
				  ps.counter > 0 &&
				  ((ps.counter+1-interim_counter_off_one) == 500000 ||
				   ((ps.counter+1-interim_counter_off_one) % 5000000 == 0 && IniGetInt (INI_FILE, (char*)"SendInterimResidues", 1)));*/
		interim_residue = INTERIM_RESIDUES &&
				  (ps.state == STATE_NORMAL || ps.state == STATE_DCHK_PASS1 || ps.state == STATE_GERB_MID_BLOCK) &&
				  (ps.counter > 0 && (ps.counter+1-interim_counter_off_one) % INTERIM_RESIDUES == 0);
		interim_file = INTERIM_FILES &&
			       (ps.state == STATE_NORMAL || ps.state == STATE_DCHK_PASS2 || ps.state == STATE_GERB_MID_BLOCK) &&
			       (ps.counter > 0 && (ps.counter+1) % INTERIM_FILES == 0);

/* Do one PRP or Proth iteration */

		gwyptimers[1] = 0.0;
		gwypstart_timer (1);

/* If we are doing one of the Gerbicz multiplies (not a squaring), then handle that here */

		if (ps.state == STATE_GERB_MID_BLOCK_MULT) {
			gwypsetnormroutine (0, 1, 0);	/* Always roundoff error check multiplies */
                        if (cufftonly)
                            gwypmul (ps.x, ps.d);   /* "Safe" multiply that does not change ps.x */
                        else
                            cuda_gwypmul (ps.x, ps.d, 3);/* "Safe" multiply that does not change ps.x */
			x = ps.d;	/* Set pointer for checking roundoff errors, sumouts, etc. */
		} else if (ps.state == STATE_GERB_END_BLOCK_MULT) {
			gwypsetnormroutine (0, 1, 0);	/* Always roundoff error check multiplies */
                        if (cufftonly)
                            gwypmul (ps.u0, ps.alt_x);
                                    /* Multiply to calc checksum #2.  u0 value can be destroyed. */
                        else
                            cuda_gwypmul (ps.u0, ps.alt_x, 3);	/* Multiply to calc checksum #2. */
			x = ps.alt_x;	/* Set pointer for checking roundoff errors, sumouts, etc. */
		} else if (ps.state == STATE_GERB_FINAL_MULT) {
			gwypcopy (ps.x, ps.u0);	// Copy x (before using it) for next Gerbicz block
			gwypsetnormroutine (0, 1, 0);	/* Always roundoff error check multiplies */
                        if (cufftonly)
                            gwypmul (ps.x, ps.d);
                                    /* "Safe" multiply to compute final d[t] value (checksum #1) */
                        else
                            cuda_gwypmul (ps.x, ps.d, 3);
                                    /* "Safe" multiply to compute final d[t] value (checksum #1) */
			x = ps.d;   /* Set pointer for checking roundoff errors, sumouts, etc. */
		}

/* Otherwise, do a squaring iteration */

		else {

/* Use state to decide which number we are squaring */

			if (ps.state == STATE_NORMAL || ps.state == STATE_DCHK_PASS1 || ps.state == STATE_GERB_MID_BLOCK) {
				x = ps.x;
				units_bit = &ps.units_bit;
			} else {			// (ps.state == STATE_DCHK_PASS2 || ps.state == STATE_GERB_END_BLOCK) {
				x = ps.alt_x;
				units_bit = &ps.alt_units_bit;
			}

/* Process this bit.  Use square carefully the first and last 30 iterations. */
/* This should avoid any pathological non-random bit pattterns.  Also square */
/* carefully during an error recovery. This will protect us from roundoff */
/* errors up to (1.0 - 0.40625). */

			if (bitval (exp, final_counter-ps.counter-1)) {
				gwypsetnormroutine (0, echk, 1);
			} else {
				gwypsetnormroutine (0, echk, 0);
			}
                        if (maxerr_recovery_mode[6] && ps.counter == last_counter) {
                            gwypsquare_carefully (x);
//                            gwypsquare (x);
                            maxerr_recovery_mode[6] = 0;
                            last_counter = 0xFFFFFFFF;
                            echk = 0;
                            care = TRUE;
			} else if (ps.counter <= (50+loopshift) || ps.counter > (/*max (gerb_counter,*/ max_counter-40)) {             // 10/04/21
                            gwypsquare_carefully (x);   // 40 in place of 50 06/04/21
//                            gwypsquare (x);
                            care = TRUE;
                        }
                        else if (cufftonly)
                            {
                            gwypsquare (x);
                            care = TRUE;
                        }
                        else if ((ps.error_check_type == ERRCHK_GERBICZ) && (ps.counter > gerb_counter))
                            {
                            gwypsquare (x);
                            care = TRUE;
                        }
                        else {
                            if (zp || generic || (ps.counter <= (51+loopshift)))
                                gpuflag = 3;
                            else if  (ps.counter >= max_counter-41)   // 41 in place of 51 06/04/21
                                gpuflag = 2;
                            else if (saving)
                                gpuflag = 2;
                            else if (restarting) {
                                gpuflag = 1;
                                restarting = FALSE;
                            }
                            else if ((ps.error_check_type == ERRCHK_DBLCHK) && (ps.counter == (ps.end_counter-1)))
                                gpuflag = 2;
                            else if (/*(ps.error_check_type == ERRCHK_GERBICZ) &&*/ (ps.counter==(ps.next_mul_counter-1)))
                                gpuflag = 2;
                            else if ((ps.error_check_type == ERRCHK_GERBICZ) && (ps.counter>(gerb_counter-1)))
                                gpuflag = 2;
                            else if ((ps.error_check_type == ERRCHK_GERBICZ)
                                && (x == ps.alt_x))
                                gpuflag = 3;
                            else if ((ps.error_check_type == ERRCHK_GERBICZ) && (x == ps.d))
                                gpuflag = 3;
                            else if ((ps.start_counter < ps.counter)&&(ps.counter < ps.end_counter))
                                gpuflag = 0;
                            else
                                gpuflag = 3;
                            cuda_gwypsquare (x, gpuflag);
                            care = FALSE;
                        }
			*units_bit = *units_bit*2;
			if (*units_bit > ((ps.residue_type == GMN_TYPE)?2*(w->n):w->n)) {
				*units_bit -= (ps.residue_type == GMN_TYPE)?2*(w->n):w->n;
			}
			if ((ps.residue_type == GMN_TYPE) && (ps.counter == final_counter_y)) {
				gwypcopy (x, ps.y);
				ps.units_bit2 = *units_bit;
				loopshift = final_counter_y;
				max_counter = final_counter;
			}

		}
		
		if (balerr) {
                    OutputBoth ((char*)"Normalization error, restarting from the beginning using fewer GPU computing...\n");
                    recovering = TRUE;
                    unlinkSaveFiles (&write_save_file_state);
                    goto restart;
                }

// introduce an error every random # iterations when debugging highly reliable error checking

//#define INTRODUCE_ERRORS
#ifdef INTRODUCE_ERRORS
		if ((rand () & 0x7FFF) == 134)  // one out of 32768
			*x += 5.0;
#endif
		
/* End iteration timing and increase count of iterations completed */

		gwypend_timer (1);
		gwyptimers[0] += gwyptimers[1];
		iters++;

/* Update min/max round-off error */

		if (echk) {
			if (ps.counter > 30 && gwyp_get_maxerr() < reallyminerr) reallyminerr = gwyp_get_maxerr ();
			if (gwyp_get_maxerr () > reallymaxerr) reallymaxerr = gwyp_get_maxerr ();
		}

/* Check for excessive roundoff error.  If round off is too large, repeat the iteration to see if this was */
/* a hardware error.  If it was repeatable then repeat the iteration using a safer, slower method.  This can */
/* happen when operating near the limit of an FFT.  NOTE: with the introduction of Gerbicz error-checking we */
/* ignore some of these errors as the Gerbicz check will catch any problems later.  However, if the round off */
/* error is really large, then results are certainly corrupt and we roll back immmediately. */

		if (echk && gwyp_get_maxerr () > allowable_maxerr) {
			if (ps.counter == last_counter && gwyp_get_maxerr () == last_maxerr[6]) {
				OutputBoth (ERROK);
				inc_error_count (3, &ps.error_count);
//				gwyp_clear_error ();
				OutputBoth (ERRMSG5);
				maxerr_recovery_mode[6] = 1;
				restart_counter = ps.counter;		/* rollback to this iteration or earlier */
				sleep5 = FALSE;
				goto restart;
			} else {
				char	msg[100];
				sprintf (msg, ERRMSG1C, gwyp_get_maxerr (), allowable_maxerr);
				sprintf (buf, ERRMSG0, ps.counter+1, final_counter, msg);
				OutputBoth (buf);
				inc_error_count (1, &ps.error_count);
				if (ps.error_check_type == ERRCHK_NONE ||
				    gwyp_get_maxerr () > IniGetFloat (INI_FILE, "RoundoffRollbackError", (float) 0.475)) {
					last_counter = ps.counter;
					last_maxerr[6] = gwyp_get_maxerr ();
					restart_counter = ps.counter;		/* rollback to this iteration or earlier */
					sleep5 = FALSE;
                                        recovering = TRUE;   // 01/02/21
					goto restart;
				}
			}
		}

/* Update counter, percentage complete */

		ps.counter++;
		w->pct_complete = (double) ps.counter * inverse_explen;
		if (ps.error_check_type == ERRCHK_DBLCHK) {
			unsigned long true_counter;
			true_counter = ps.start_counter + ((ps.counter - ps.start_counter) >> 1);
			if (ps.state == STATE_DCHK_PASS2) true_counter += ((ps.end_counter - ps.start_counter) >> 1);
			w->pct_complete = (double) true_counter * inverse_explen;
		}

/* Output the title every so often */

		actual_frequency = (int) (ITER_OUTPUT * output_title_frequency);
		if (actual_frequency < 1) actual_frequency = 1;
		if (ps.counter % actual_frequency == 0 || first_iter_msg) {
			if ((ps.residue_type == PROTH_TYPE)||(ps.residue_type == GMN_TYPE))
				sprintf (buf, "%.*f%% of Proth %s", (int) PRECISION, 100.0*trunc_percent (w->pct_complete), string_rep);
			else if (ps.residue_type == PRP_TYPE_SPRP)
				sprintf (buf, "%.*f%% of Strong PRP %s", (int) PRECISION, 100.0*trunc_percent (w->pct_complete), string_rep);
			else
				sprintf (buf, "%.*f%% of Fermat PRP %s", (int) PRECISION, 100.0*trunc_percent (w->pct_complete), string_rep);
			title (buf);
		}

//		ReplaceableLine (2);	/* Replace line */ 

/* Print a message every so often */

		actual_frequency = (int) (ITER_OUTPUT * output_frequency);
		if (actual_frequency < 1) actual_frequency = 1;
		if ((ps.counter % actual_frequency == 0 && ps.state != STATE_GERB_MID_BLOCK_MULT &&
		     ps.state != STATE_GERB_END_BLOCK && ps.state != STATE_GERB_END_BLOCK_MULT &&
		     ps.state != STATE_GERB_FINAL_MULT) || first_iter_msg) {
			sprintf (buf, "Iteration: %ld / %ld [%.*f%%]",
				 ps.counter, final_counter, (int) PRECISION, 100.0*trunc_percent (w->pct_complete));
			/* Append a short form total errors message */
			if ((error_count_messages & 0xFF) == 1)
				make_error_count_message (ps.error_count, error_count_messages, buf + strlen (buf),
							  (int) (sizeof (buf) - strlen (buf)));
			/* Truncate first message */
			if (first_iter_msg) {
				strcat (buf, ".\r");	// JP 17/04/20
                                gwypclear_timer (0);
				first_iter_msg = FALSE;
			}
			/* In v28.5 and later, format a consise message including the ETA */
			else if (!CLASSIC_OUTPUT) {
				double speed;
				/* Append roundoff error */
				if ((OUTPUT_ROUNDOFF || ERRCHK) && reallymaxerr >= 0.001) {
					sprintf (buf+strlen(buf), ", roundoff: %5.3f", reallymaxerr);
					if (!CUMULATIVE_ROUNDOFF) reallyminerr = 1.0, reallymaxerr = 0.0;
				}
				/* Append ms/iter */
				speed = gwyptimer_value (0) / (double) iters;
				sprintf (buf+strlen(buf), ", ms/iter: %6.3f", speed * 1000.0);
				gwypclear_timer (0);
				iters = 0;
				/* Append ETA */
				formatETA ((final_counter - ps.counter) * speed, buf+strlen(buf));
#if defined (_CONSOLE) || defined (__linux__) || defined (__FreeBSD__) || defined (__APPLE__)
				strcat (buf, "\r");
#else
				strcat (buf, "\n");
#endif
			}
			/* Format the classic (pre-v28.5) message */
			else {
				/* Append optional roundoff message */
				if (ERRCHK && ps.counter > 30) {
					sprintf (buf+strlen(buf), ".  Round off: %10.10f to %10.10f", reallyminerr, reallymaxerr);
					if (!CUMULATIVE_ROUNDOFF) reallyminerr = 1.0, reallymaxerr = 0.0;
				}
				if (CUMULATIVE_TIMING) {
					strcat (buf, ".  Total time: ");
					gwypprint_timer (0, TIMER_NL);
				} else {
					strcat (buf, ".  Per iteration time: ");
					gwypdivide_timer (0, iters);
					gwypprint_timer (0, TIMER_NL | TIMER_CLR);
					iters = 0;
				}
			}
			
			OutputStr (buf);


/* Output a verbose message showing the error counts.  This way a user is likely to */
/* notice a problem without reading the results.txt file. */

			if ((error_count_messages & 0xFF) >= 2 &&
			    make_error_count_message (ps.error_count, error_count_messages, buf, sizeof (buf)))
				OutputStr (buf);
		}

		ReplaceableLine (2);	/* Replace line */ 

/* Print a results file message every so often */

		if ((ps.counter % ITER_OUTPUT_RES == 0 && ps.state != STATE_GERB_MID_BLOCK_MULT &&
		     ps.state != STATE_GERB_END_BLOCK && ps.state != STATE_GERB_END_BLOCK_MULT &&
		     ps.state != STATE_GERB_FINAL_MULT) || (NO_GUI && stop_reason)) {
			sprintf (buf, "Iteration %ld / %ld\n", ps.counter, final_counter);
			writeResults (buf);
		}

/* If double-checking, at end of pass 1 rollback counter and start computing alt_x. */
/* If double-checking, at end of pass 2 compare values and move onto next block. */

		if (ps.state == STATE_DCHK_PASS1) {
			if (ps.counter < ps.end_counter)
				;   // Do next iteration
			else if (ps.counter == ps.end_counter) {
                                    // Switch to alt_x computations
				ps.state = STATE_DCHK_PASS2;
				ps.counter = ps.start_counter;
			} else {					// Can't happen
				OutputBoth (ERRMSG9);
				inc_error_count (6, &ps.error_count);
				restart_counter = -1;			/* rollback to any save file */
				sleep5 = FALSE;
				goto restart;
			}
		}

		if (ps.state == STATE_DCHK_PASS2) {
			if (ps.counter < ps.end_counter)
                                ;   // Do next iteration
			else if (ps.counter == ps.end_counter) {
 				if (!areTwoPRPValsEqual ( (ps.residue_type == GMN_TYPE)?2*(w->n):w->n, ps.x, ps.units_bit, ps.alt_x, ps.alt_units_bit)) {
					sprintf (buf, ERRMSG60, ps.start_counter);
					OutputBoth (buf);
					inc_error_count (7, &ps.error_count);
					restart_counter = ps.start_counter;		/* rollback to this iteration */
					sleep5 = FALSE;
					goto restart;
				}
				/* If doing a full double-check, start next block of iterations */
				if (ps.error_check_type == ERRCHK_DBLCHK) {
					ps.state = STATE_DCHK_PASS1;
					ps.start_counter = ps.counter;
					ps.end_counter = ps.start_counter + IniGetInt (INI_FILE, (char*)"PRPDoublecheckCompareInterval", 100000);
					if (ps.end_counter > final_counter) {
                                            ps.end_counter = final_counter;
                                        }
				}
				/* Otherwise, we're doing the first or last few iterations of a Gerbicz error-check. */
				/* Set state to start next Gerbicz block (in case we just computed prp_base^k). */
				else {
                                    ps.state = STATE_GERB_START_BLOCK;
				}
				/* We've reached a verified iteration, create a save file and mark it highly reliable. */
				/* But if there are less than 1000 iterations left on a reliable machine */
				/* don't bother creating the save file. */
				if (final_counter - ps.counter >= 1000 || ps.error_count > 0) {
					saving = TRUE;
					saving_highly_reliable = TRUE;
				}
			} else {					// Can't happen
				OutputBoth (ERRMSG9);
				inc_error_count (6, &ps.error_count);
				restart_counter = -1;			/* rollback to any save file */
				sleep5 = FALSE;
				goto restart;
			}
		}

/* If Gerbicz error-checking, handle all the possible Gerbicz states.  See if this is an L-th iteration that needs */
/* to do a checksum multiply.  Also check if this is the L^2-th iteration where we switch to an alternate method to */
/* compute checksum #2. */
/* NOTE: We do not need to worry about shift_counts in the checksum as both checksums will end up with the same shift count. */

		// Just did a normal PRP squaring

		if (ps.state == STATE_GERB_MID_BLOCK) {
			if (ps.counter < ps.next_mul_counter)
                            ;
						// Do next iteration
			else if (ps.counter == ps.end_counter) {	// Delay last checksum #1 multiply, start checksum #2 calculation
				if (IniGetInt (INI_FILE, (char*)"GerbiczVerbosity", 1) > 1) OutputStr ((char*)"Start Gerbicz error check.\n");
				// At end of Gerbicz block, switch to "L" squarings of alt_x to create Gerbicz checksum #2 value
				gwypcopy (ps.d, ps.alt_x);	// Copy d[t-1] to alt_x
				ps.state = STATE_GERB_END_BLOCK;	// Squaring alt_x state
				ps.counter -= ps.L;			// L squarings
			} else if (ps.counter == ps.next_mul_counter) {	// Do a checksum #1 multiply next
				// Back counter up by one and do one multiply in the computation of Gerbicz checksum #1 value
				ps.state = STATE_GERB_MID_BLOCK_MULT;
				ps.counter -= 1;
			} else {					// Can't happen
				OutputBoth (ERRMSG9);
				inc_error_count (6, &ps.error_count);
				restart_counter = -1;			/* rollback to any save file */
				sleep5 = FALSE;
				goto restart;
			}
		}

		// Just did a a checksum #1 multiply at the end of a block of L normal squarings


		else if (ps.state == STATE_GERB_MID_BLOCK_MULT) {
			if (ps.counter < ps.end_counter) {		// In middle of Gerbicz block, do another "L" squarings
				ps.state = STATE_GERB_MID_BLOCK;
				ps.next_mul_counter += ps.L;
			} else {					// Can't happen
				OutputBoth (ERRMSG9);
				inc_error_count (6, &ps.error_count);
				restart_counter = -1;			/* rollback to any save file */
				sleep5 = FALSE;
				goto restart;
			}
		}

		// Just did a checksum #2 squaring
		else if (ps.state == STATE_GERB_END_BLOCK) {
			if (ps.counter < ps.end_counter)
				;		// Do next iteration in computing checksum #2
			else if (ps.counter == ps.end_counter) {	// Next do final multiply in computing checksum #2
				ps.state = STATE_GERB_END_BLOCK_MULT;
				ps.counter -= 1;
			} else {					// Can't happen
				OutputBoth (ERRMSG9);
				inc_error_count (6, &ps.error_count);
				restart_counter = -1;			/* rollback to any save file */
				sleep5 = FALSE;
				goto restart;
			}
		}

		// Just did the checksum #2 multiply at the end of L checksum #2 squarings
		else if (ps.state == STATE_GERB_END_BLOCK_MULT) {
			if (ps.counter == ps.end_counter) {		// Next do final multiply in computing checksum #1
				ps.state = STATE_GERB_FINAL_MULT;
				ps.counter -= 1;
			} else {					// Can't happen
				OutputBoth (ERRMSG9);
				inc_error_count (6, &ps.error_count);
				restart_counter = -1;			/* rollback to any save file */
				sleep5 = FALSE;
				goto restart;
			}
		}

		// Just did the final checksum #1 multiply which we delayed until checksum #2 completed
		else if (ps.state == STATE_GERB_FINAL_MULT) {
			double	gerbicz_block_size_adjustment;
			// We adjust the compare interval size downward when errors occur and upwards when they dont.
			// That way buggy machines will lose fewer iterations when rolling back.
			gerbicz_block_size_adjustment = IniGetFloat (INI_FILE, "PRPGerbiczCompareIntervalAdj", 1.0);
			if (gerbicz_block_size_adjustment < 0.001 || gerbicz_block_size_adjustment > 1.0) gerbicz_block_size_adjustment = 0.5;
			// Compare alt_x, d (the two Gerbicz checksum values that must match)
			if (!areTwoPRPValsEqual ((ps.residue_type == GMN_TYPE)?2*(w->n):w->n, ps.alt_x, 0, ps.d, 0)) {
				sprintf (buf, ERRMSG70, ps.start_counter);
				OutputBoth (buf);
				gerbicz_block_size_adjustment *= 0.25;		/* This will halve next L */
				if (gerbicz_block_size_adjustment < 0.001) gerbicz_block_size_adjustment = 0.001;
				IniWriteFloat (INI_FILE, "PRPGerbiczCompareIntervalAdj", (float) gerbicz_block_size_adjustment);
				inc_error_count (7, &ps.error_count);
				restart_counter = ps.start_counter;		/* rollback to this iteration */
				sleep5 = FALSE;
                                        if (!generic) {
                                            MAXBPD = 35.0;
                                            OutputBoth ((char*)"Restarting from the beginning...\n");
                                            unlinkSaveFiles (&write_save_file_state);
                                        }
				goto restart;
			}
			if (IniGetInt (INI_FILE, (char*)"GerbiczVerbosity", 1)) {
				clearline(100);
                                if (verbose) {
                                    sprintf (buf, "Gerbicz error check passed at iteration %ld.\n", ps.counter);
                                    OutputBoth (buf);
                                }
                                else {
                                    sprintf (buf, "Gerbicz error check passed at iteration %ld.\r", ps.counter);
                                    OutputStr (buf);
                                }
			}
			gerbicz_block_size_adjustment *= 1.0473;		/* 30 good blocks to double L */
			if (gerbicz_block_size_adjustment > 1.0) gerbicz_block_size_adjustment = 1.0;
			IniWriteFloat (INI_FILE, "PRPGerbiczCompareIntervalAdj", (float) gerbicz_block_size_adjustment);
			/* Start next Gerbicz block.  Both x and alt_x must be identical at start of next block. */
			ps.state = STATE_GERB_START_BLOCK;
			gwypswap (ps.alt_x, ps.u0);
			ps.alt_units_bit = ps.units_bit;
			/* We've reached a verified iteration, create a save file and mark it highly reliable. */
			/* But if there are less than 1000 iterations left on a reliable machine, don't bother creating the save file. */
			if (1/*final_counter - ps.counter >= 1000 || gerbicz_block_size_adjustment < 1.0*/) {
				saving = TRUE;
				saving_highly_reliable = TRUE;
			}
		}

/* Write results to a file every DISK_WRITE_TIME minutes */

		if (saving) {
                    if (stop_reason)
                        gwypend_timer (2);
                    if (! writeSaveFile ( &write_save_file_state, w, &ps)) {
                        sprintf (buf, WRITEFILEERR, filename);
                        OutputBoth (buf);
                    }
                    // Mark save files that contain verified computations.  This will keep the save file
                    // for a longer period of time (i.e. will not be replaced by a save file that does
                    // not also contain verified computations).
                    if (saving_highly_reliable)
                        setWriteSaveFileSpecial (&write_save_file_state);
                    saving = FALSE;
                    saving_highly_reliable = FALSE;
		}

/* If an escape key was hit, write out the results and return */

		if (stop_reason) {
			stop_reason = FALSE;
//                        stop_counter = ps.counter;
			if (ps.residue_type == PROTH_TYPE || ps.residue_type == GMN_TYPE)
				sprintf (buf, "Stopping Proth prime test of %s at iteration %ld [%.*f%%]\n",
					string_rep, ps.counter, (int) PRECISION, trunc_percent (100.0*w->pct_complete));
			else if (ps.residue_type == PRP_TYPE_SPRP)
				sprintf (buf, "Stopping Strong PRP test of %s at iteration %ld [%.*f%%]\n",
					string_rep, ps.counter, (int) PRECISION, trunc_percent (100.0*w->pct_complete));
				else
				sprintf (buf, "Stopping Fermat PRP test of %s at iteration %ld [%.*f%%]\n",
					string_rep, ps.counter, (int) PRECISION, trunc_percent (100.0*w->pct_complete));
			clearline(100);
			OutputStr (buf);
                        if (ps.error_check_type == ERRCHK_GERBICZ || ps.error_check_type == ERRCHK_DBLCHK) {
                            gwypfree (ps.alt_x);
                            nbllr_frees++;
                        }
                        if (ps.error_check_type == ERRCHK_GERBICZ) {
                            gwypfree(ps.u0);
                            nbllr_frees++;
                            gwypfree(ps.d);
                            nbllr_frees++;
                        }
			goto exit;
		}

/* Output the 64-bit residue at specified interims. */

		if (interim_residue) {
                        if (ps.residue_type == GMN_TYPE)
                            tmp3 = newgiant (4*FFTLEN*sizeof(double)/sizeof(short) + 16);
                        else
                            tmp3 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);
                        nbllr_mallocs++;
			if (gwyptogiant (x, tmp3)) {
				free (tmp3);
                                nbllr_frees++;
				OutputBoth (ERRMSG80);
				inc_error_count (2, &ps.error_count);
				last_counter = ps.counter;		/* create save files before and after this iteration */
				restart_counter = -1;			/* rollback to any save file */
				sleep5 = TRUE;
				goto restart;
			}
			rotategp (tmp3, (ps.residue_type == GMN_TYPE)?2*(w->n):w->n, *units_bit);
			if (interim_mul)
				basemulg (tmp3, w, ps.prp_base, -1);
			if (w->known_factors && ps.residue_type != PRP_TYPE_COFACTOR)
				modg (N, tmp3);
			if ((ps.residue_type == PROTH_TYPE)||(ps.residue_type == GMN_TYPE))
				sprintf (buf, "%s interim Proth residue %08lX%08lX at iteration %ld\n",
					string_rep, (unsigned long) tmp3->n[1], (unsigned long) tmp3->n[0],
					ps.counter - interim_counter_off_one);
			else if (ps.residue_type == PRP_TYPE_SPRP)
				sprintf (buf, "%s interim Strong PRP residue %08lX%08lX at iteration %ld\n",
					string_rep, (unsigned long) tmp3->n[1], (unsigned long) tmp3->n[0],
					ps.counter - interim_counter_off_one);
			else
				sprintf (buf, "%s interim Fermat PRP residue %08lX%08lX at iteration %ld\n",
					string_rep, (unsigned long) tmp3->n[1], (unsigned long) tmp3->n[0],
					ps.counter - interim_counter_off_one);
			OutputBoth (buf);
                        free (tmp3);
                        nbllr_frees++;
		}

/* Write a save file every INTERIM_FILES iterations. */

		if (interim_file) {
			char	interimfile[32];
			writeSaveFileState state;
			sprintf (interimfile, "%s.%03ld", filename, ps.counter / INTERIM_FILES);
			writeSaveFileStateInit (&state, interimfile, 0);
			state.num_ordinary_save_files = 99;
			writeSaveFile (&state, w, &ps);
		}
		bit++;
	}
        
/* Restaure the default adjustment, and free up some memory */

	if (ps.error_check_type == ERRCHK_GERBICZ) {
            IniWriteFloat (INI_FILE, "PRPGerbiczCompareIntervalAdj", 1.0);
            gwypfree (ps.u0);
            nbllr_frees++;
            gwypfree (ps.d);
            nbllr_frees++;
        }

        
/* Make sure PRP state is valid.  We cannot be in the middle of a double-check or in the middle of a Gerbicz block */

	if (ps.state != STATE_NORMAL && ps.state != STATE_DCHK_PASS1 && ps.state != STATE_GERB_START_BLOCK) {
		OutputBoth (ERRMSG90);
		inc_error_count (6, &ps.error_count);
		restart_counter = -1;			/* rollback to any save file */
		sleep5 = FALSE;
		goto restart;
	}
	

/* See if we've found a probable prime.  If not, format a 64-bit residue. */

	if (gwyptogiant (ps.x, tmp)) {
		OutputBoth (ERRMSG8);
		inc_error_count (2, &ps.error_count);
		restart_counter = -1;			/* rollback to any save file */
		sleep5 = TRUE;
		goto restart;
	}
	
	rotategp (tmp, (ps.residue_type == GMN_TYPE)?2*(w->n):w->n, ps.units_bit);

        if (mul_final)
		basemulg (tmp, w, ps.prp_base, mul_final);
         
                
	if (w->known_factors && ps.residue_type != PRP_TYPE_COFACTOR)
		modg (N, tmp);
                
	if (ps.residue_type != GMN_TYPE) {
		PositiveResult = isPRPg (tmp, N, w, ps.prp_base, ps.residue_type);
		if (!PositiveResult) {
			if (ps.residue_type == PROTH_TYPE)
				iaddg (1, tmp);	// to match with V3.8.23 Proth residue
//			sprintf (res64, "%08lX%08lX", (unsigned long) ((tmp->sign > 1) ? tmp->n[1] : 0), (unsigned long) tmp->n[0]);
                        formatresidue (tmp, res64);
			have_res2048 = (tmp->sign > 64);
			if (have_res2048) {
				int i;
				for (i = 63; i >= 0; i--) sprintf (res2048+504-i*8, "%08lX", (unsigned long) tmp->n[i]);
			}
		}
	}

/* If we are doing highly reliable error checking, then make sure the calculation of the final residue was error free! */
/* Perform the same calculations above but use alt_x. */

	if (ps.state == STATE_DCHK_PASS1 || ps.state == STATE_GERB_START_BLOCK) {
		int	alt_match, alt_isProbablePrime;
		if (gwyptogiant (ps.alt_x, tmp2)) {
			OutputBoth (ERRMSG8);
			inc_error_count (2, &ps.error_count);
			restart_counter = -1;			/* rollback to any save file */
			sleep5 = TRUE;
			goto restart;
		}
		
		rotategp (tmp2, (ps.residue_type == GMN_TYPE)?2*(w->n):w->n, ps.alt_units_bit);
                
		if (mul_final)
			basemulg (tmp2, w, ps.prp_base, mul_final);
                
		if (w->known_factors && ps.residue_type != PRP_TYPE_COFACTOR)
			modg (N, tmp2);
                
		if (ps.residue_type != GMN_TYPE) {
			alt_isProbablePrime = isPRPg (tmp2, N, w, ps.prp_base, ps.residue_type);
			alt_match = (PositiveResult == alt_isProbablePrime);
			if (alt_match && !alt_isProbablePrime) {
				if (ps.residue_type == PROTH_TYPE)
					iaddg (1, tmp2);	// to match with V3.8.23 Proth residue
//				sprintf (alt_res64, "%08lX%08lX", (unsigned long) ((tmp2->sign > 1) ? tmp2->n[1] : 0), (unsigned long) tmp2->n[0]);
                                formatresidue (tmp2, alt_res64);
				alt_match = !strcmp (res64, alt_res64);
			}
		}
		else {
			alt_match = TRUE;// Test only x == alt-x
		}
		
		if (!alt_match || (gcompg (tmp, tmp2)!=0)) {
			OutputBoth (ERRMSG80);
			inc_error_count (2, &ps.error_count);
			restart_counter = -1;
                            /* rollback to any save file */
			sleep5 = TRUE;
			goto restart;
		}
		
		gwypfree (ps.alt_x);
                nbllr_frees++;
	}

	if (ps.residue_type == GMN_TYPE) {
		ps.gx = tmp;
                
		if (gwyptogiant (ps.y, tmp2)) {
			OutputBoth (ERRMSG8);
			inc_error_count (2, &ps.error_count);
			restart_counter = -1;
                            /* rollback to any save file */
			sleep5 = TRUE;
			goto restart;
		}
		
		rotategp (tmp2, (ps.residue_type == GMN_TYPE)?2*(w->n):w->n, ps.units_bit2);
		ps.gy = tmp2;

		gwypfree (ps.x);
                nbllr_frees++;
		gwypfree (ps.y);
                nbllr_frees++;
                 
	/* Delete the continuation files. */

		unlinkSaveFiles (&write_save_file_state);
	
	}
	else {

		/* Print results */

		if (PositiveResult) {
			if (ps.residue_type == PROTH_TYPE)
				sprintf (buf, "%s is prime", string_rep);
			else	if (ps.residue_type == PRP_TYPE_SPRP) {
				sprintf (buf, "%s is a Strong Probable prime", string_rep);
				if (ps.prp_base != 3)
					sprintf (buf+strlen(buf), " (%u-PRP)", ps.prp_base);
			}
			else {
				sprintf (buf, "%s is a Fermat Probable prime", string_rep);
				if (ps.prp_base != 3)
					sprintf (buf+strlen(buf), " (%u-PRP)", ps.prp_base);
			}

			strcat (buf, "!");
			sprintf (buf+strlen(buf), " (%d decimal digits)", nbdg); 
			*res = TRUE;
		}
		else {
			sprintf (buf, "%s is not prime.  ", string_rep);
			if (ps.prp_base != 3)
				sprintf (buf+strlen(buf), "Base-%u ", ps.prp_base);
			if (ps.residue_type == PROTH_TYPE)
				sprintf (buf+strlen(buf), "Proth ");
			else if (ps.residue_type != PRP_TYPE_FERMAT)
				sprintf (buf+strlen(buf), "Type-%d ", ps.residue_type);
			sprintf (buf+strlen(buf), "RES64: %s.", res64);
			*res = FALSE;
		}

/* Print known factors */

		if (string_rep_truncated) {
			char	*bigbuf;
			bigbuf = (char *) malloc (strlen (w->known_factors) + 100);
			if (bigbuf != NULL) {
				sprintf (bigbuf, "Known factors used for PRP test were: %s\n", w->known_factors);
				OutputStr (bigbuf);
				if ((PositiveResult && IniGetInt (INI_FILE, (char*)"OutputPrimes", 0)) ||
			    (!PositiveResult && IniGetInt (INI_FILE, (char*)"OutputComposites", 0)))
					writeResults (bigbuf);
				free (bigbuf);
			}
		}

	/* Delete the continuation files. */

		unlinkSaveFiles (&write_save_file_state);
	

		clearline (100);

		free (tmp);
                nbllr_frees++;
		free (tmp2);
                nbllr_frees++;
		gwypfree (ps.x);
                nbllr_frees++;
		gwypfree (ps.y);
                nbllr_frees++;
                free (exp);
                nbllr_frees++;


#if defined(WIN32) && !defined(_CONSOLE)

		sprintf (buf+strlen(buf), "  Time : "); 
		ReplaceableLine (2);	/* Replace line */ 

#else

		clearline(100);

		if (*res) {
#if defined (_CONSOLE)
			hConsole = GetStdHandle(STD_OUTPUT_HANDLE);	// Access to Console attributes
			SetConsoleTextAttribute(hConsole, BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_RED);
			OutputBoth(buf);
			SetConsoleTextAttribute(hConsole, FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED);
#else
			OutputStr((char*)"\033[7m");
			OutputBoth(buf);
			OutputStr((char*)"\033[0m");
#endif
		}
		else
			OutputBoth(buf);

		sprintf (buf, (char*)"  Time : "); 

#endif

/* Output the final timings */

		gwypend_timer (2);
		gwypwrite_timer (buf+strlen(buf), 2, TIMER_CLR | TIMER_NL); 
		OutputBoth(buf);

/* Cleanup and return */

	}
	
	lasterr_point = 0;
	return (TRUE);

/* An error occured, output a message saying we are restarting, sleep, */
/* then try restarting at last save point. */

restart:
	if (sleep5) 
		OutputBoth (ERRMSG2);
        if (generic)    // 25/05/22
            OutputBoth (ERRMSG3);

/* Save the incremented error count to be used in the restart rather than the error count read from a save file */

//	restart_error_count = ps.error_count;


/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
		return (FALSE);
	}

/* Return so that last continuation file is read in */

        IniWriteInt(INI_FILE, (char*)"FFT_Increment", nbfftinc =  IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
        if (recovering || (nbfftinc == maxfftinc)) {
            if (nbfftinc >= maxfftinc) {
                OutputBoth ((char*)"Too much retries, restarting from the beginning using fewer GPU computing...\n");
                unlinkSaveFiles (&write_save_file_state);
            } 
            else
                OutputBoth ((char*)"Restarting using fewer GPU computing...\n");
            recovering = TRUE;
//            unlinkSaveFiles (&write_save_file_state);
        }
        
	free (tmp);
        nbllr_frees++;
	free (tmp2);
        nbllr_frees++;
	gwypfree (ps.x);
        nbllr_frees++;
	gwypfree (ps.y);
        nbllr_frees++;
	free (exp);
        nbllr_frees++;
	return (-1);;

exit:
	free (tmp);
        nbllr_frees++;
	free (tmp2);
        nbllr_frees++;
	gwypfree (ps.x);
        nbllr_frees++;
	gwypfree (ps.y);
        nbllr_frees++;
	free (exp);
        nbllr_frees++;
//        gwypdone();
	*res = FALSE;		// To avoid credit message !
	return (FALSE);
}

/* Test for an N+1 Frobenius probable prime -- gwypsetup has already been called. */

int commonFrobeniusPRP (
	long P,
        long Q,
	int	*res, char *str)
{
//	unsigned long bit, firstone = 0, iters, D, bitv;
	unsigned long bit, iters, D;
	gwypnum x, y, gwA, gw2, gwQ;
	giant	tmp, tmp2, tmp3, A;
	char	filename[20], buf[sgkbufsize+256], fft_desc[256];
	long	write_time = DISK_WRITE_TIME * 60;
	int	echk, saving, stopping, gpuflag = 3;
	time_t	start_time, current_time;
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;


/* Allocate memory */
	x = gwypalloc ();
	y = gwypalloc ();
	gwA = gwypalloc ();
	gw2 = gwypalloc ();
	gwQ = gwypalloc ();

	tmp = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);
	tmp2 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);
	tmp3 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);
	A = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);

	D = P*P - 4*Q;
//	Compute needed large integer and gwypnum constants

	itogwyp (2, gw2);
	itog (Q, tmp);	        // Compute gwQ for Frobenius test
	if (Q < 0)
            addg (N, tmp);	// only a positive giant can be converted to gwypnum
	gianttogwyp (tmp, gwQ);
	gtog (tmp, A);		// Compute A = P*P*Q^-1 - 2 mod N
//	invg (N, A);
	gwypinvg (N, A);
	ulmulg (P*P, A);
	iaddg (-2, A);
	modg (N, A);
	gianttogwyp (A, gwA);
        nbdg = gnbdg (N, 10);
        
	*res = TRUE;		/* Assume it is a probable prime */
	
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
                restarting = TRUE;
		goto Frobeniusresume;
	}
	else {

/* Init, compute (N+1)/2 to compute x and y mod N */

		gtog (N, tmp3);
		iaddg (1, tmp3);
		gshiftright (1, tmp3);
		Nlen = bitlen (tmp3);

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
			itogwyp (2, x);		// Initial values
//                        itog (2, tmp);
//                        gianttogwyp (tmp, x);
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
	sprintf (buf, "%s, P = %ld, Q = %ld\n", fft_desc, P, Q);
#else
	sprintf (buf, "%s, P = %ld, Q = %ld", fft_desc, P, Q);
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

	iters = 0;
        it = 0;                 // JP 24/01/21
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
//			if (abs(inc)==1)
                            gwypsetaddin (0);
			if ((bit+26 < Nlen) && (bit > 26) &&
				((bit != lasterr_point) || (!maxerr_recovery_mode[1] && !maxerr_recovery_mode[2]))) {
                            gpuflag = 3;
                            if (cufftonly)
                                gwypmul (y, x);
                            else
                                cuda_gwypmul (y, x, gpuflag);
                            care = FALSE;
			}
			else {
                            gwypmul_carefully (y, x);
                            care = TRUE;
			}
			CHECK_IF_ANY_ERROR(x, (bit), Nlen, 1)
			gwypsub3 (x, gwA, x);
			if (abs(inc)==1)
				gwypsetaddin (-2);
			if ((bit+26 < Nlen) && (bit > 26) &&
				((bit != lasterr_point) || !maxerr_recovery_mode[2])) {
                            if (cufftonly)
                                gwypsquare (y);
                            else
                                cuda_gwypsquare (y,gpuflag);
                            care = FALSE;
			}
			else {
                            gwypsquare_carefully (y);
                            care = TRUE;
			}
			CHECK_IF_ANY_ERROR(y, (bit), Nlen, 2)
			if (abs(inc)!=1)
				gwypsubquick (gw2, y);
		}
		else {
//			if (abs(inc)==1)
                            gwypsetaddin (0);
			if ((bit+26 < Nlen) && (bit > 26) &&
				((bit != lasterr_point) || (!maxerr_recovery_mode[3] && !maxerr_recovery_mode[4]))) {
                            gpuflag = 3;
                            if (cufftonly)
                                gwypmul (x, y);
                            else
                                cuda_gwypmul (x, y, gpuflag);
                            care = FALSE;
			}
			else {
                            gwypmul_carefully (x, y);
                            care = TRUE;
			}
			CHECK_IF_ANY_ERROR(y, (bit), Nlen, 3)
			gwypsub3 (y, gwA, y);
			if (abs(inc)==1)
                            gwypsetaddin (-2);
			if ((bit+26 < Nlen) && (bit > 26) &&
				((bit != lasterr_point) || !maxerr_recovery_mode[4])) {
                            if (cufftonly)
                                gwypsquare (x);
                            else
                                cuda_gwypsquare (x,gpuflag);
                            care = FALSE;
			}
			else {
                            gwypsquare_carefully (x);
                            care = TRUE;
			}
			CHECK_IF_ANY_ERROR(x, (bit), Nlen, 4)
			if (abs(inc)!=1)
				gwypsubquick (gw2, x);
		}

 /* That iteration succeeded, bump counters */

		if (bit == lasterr_point)
			saving = 1;	// Be sure to restart after this recovery iteration!
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
			sprintf (buf, fmt_mask, quotient?string_rep:str, bit, Nlen, pct);
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
                                gwypfree (gwQ);
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
			sprintf (buf, "%s interim residue %s at bit %ld\n", quotient?string_rep:str, res64, bit);
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
	
/* See if we've found a Lucas probable prime.  If not, format a 64-bit residue. */

	clearline (100);

	gwyptogiant (x, tmp);			// V(m)
	gwyptogiant (y, tmp2);			// V(m+1)
	mulg (A, tmp);				// A*V(m)
	gshiftleft (1, tmp2);			// 2*V(m+1)
	subg (tmp2, tmp);			// A*V(m)-2*V(m+1)
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
			if (bpsw)
				sprintf (buf, "%s is not prime(P = %ld, Q = %ld), BPSW RES64: %s", quotient?string_rep:str, P, Q, res64);
			else
				sprintf (buf, "%s is not prime(P = %ld, Q = %ld), Lucas RES64: %s", quotient?string_rep:str, P, Q, res64);
		else
			if (bpsw)
				sprintf (buf, "%s is Fermat PSP, but composite!! (P = %ld, Q = %ld), BPSW RES64: %s", quotient?string_rep:str, P, Q, res64);
			else
				sprintf (buf, "%s is Fermat PSP, but composite!! (P = %ld, Q = %ld), Lucas RES64: %s", quotient?string_rep:str, P, Q, res64);
	}
	if (*res) {	// N may be prime ; do now the Frobenius PRP test
		_unlink (filename);	          // Remove Lucas save file

		tempFileName (filename, 'F', N);  // Reinit file name

		bit = 1;
                gwypcopy (gwQ, y);
//		lasterr_point = 0;			// Reset a possible Lucas roundoff error point
		if (IniGetInt(INI_FILE, (char*)"LucasPRPtest", 0))
			if (bpsw)
				sprintf (buf, "%s is BPSW PRP, Starting Frobenius test sequence\n", quotient?string_rep:str);
			else
				sprintf (buf, "%s is Lucas PRP, Starting Frobenius test sequence\n", quotient?string_rep:str);
		else
			if (bpsw)
				sprintf (buf, "%s is Fermat and BPSW PRP, Starting Frobenius test sequence\n", quotient?string_rep:str);
			else
				sprintf (buf, "%s is Fermat and Lucas PRP, Starting Frobenius test sequence\n", quotient?string_rep:str);
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
		sprintf (buf, "%s, Q = %ld\n", fft_desc, Q);
#else
		sprintf (buf, "%s, Q = %ld", fft_desc, Q);
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
                it = 0;                 // JP 24/01/21
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

			if ((bit+25 < Nlen) && (bit > 25) && ((bit != lasterr_point) || !maxerr_recovery_mode[6])) {
                            if (zp || generic)
                                gpuflag = 3;
                            else if(bit==1 || it==0)  
                                {gpuflag = 3;it = 1;}
                            else if (restarting)  {
                                gpuflag = 3;
                                restarting = FALSE;
                            }
                            else  if(bit != (lasterr_point-1)) 
                                gpuflag = ((bit==Nlen-26) || saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0;
                            else
                                gpuflag = 3;
                            if (cufftonly)
                                gwypsquare (y);
                            else
                                cuda_gwypsquare (y, gpuflag);
                            care = FALSE;
                        }
			else {
                            gwypsquare_carefully (y);
                            care = TRUE;
			}

			CHECK_IF_ANY_ERROR (y, (bit), Nlen, 6);

/* That iteration succeeded, bump counters */

			if (bit == lasterr_point)
                            saving = 1;					// Be sure to restart after this recovery iteration!
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
				sprintf (buf, fmt_mask,quotient?string_rep:str , bit, Nlen, pct);
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
                                        gwypfree (gwQ);
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
			sprintf (buf, "%s interim residue %s at bit %ld\n", quotient?string_rep:str, res64, bit);
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
//		gwypmul_carefully (x, y);	// y = B*V(m)-2
                if (cufftonly)
                    gwypmul (x, y);
                else
                    cuda_gwypmul (x, y, 3);
		CHECK_IF_ANY_ERROR (y, (Nlen), Nlen, 6);
//		care = FALSE;
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
                            if (bpsw)
                                sprintf (buf, "%s is BPSW PSP (P = %ld, Q = %ld), but composite!!. Frobenius RES64: %s", quotient?string_rep:str, P, Q, res64);
                            else
                                sprintf (buf, "%s is Lucas PSP (P = %ld, Q = %ld), but composite!!. Frobenius RES64: %s", quotient?string_rep:str, P, Q, res64);
			else
                            if (bpsw)
                                sprintf (buf, "%s is Fermat and BPSW PSP (P = %ld, Q = %ld), but composite!!. Frobenius RES64: %s", quotient?string_rep:str, P, Q, res64);
                            else
                                sprintf (buf, "%s is Fermat and Lucas PSP (P = %ld, Q = %ld), but composite!!. Frobenius RES64: %s", quotient?string_rep:str, P, Q, res64);
		}
	}

/* Print results.  */

	clearline (100);

	if (*res)
            if (IniGetInt(INI_FILE, (char*)"LucasPRPtest", 0))
                if (bpsw)
                    sprintf (buf, "%s is BPSW and Frobenius PRP! (P = %ld, Q = %ld, D = %ld)", quotient?string_rep:str, P, Q, D);
                else
                    sprintf (buf, "%s is Lucas and Frobenius PRP! (P = %ld, Q = %ld, D = %ld)", quotient?string_rep:str, P, Q, D);
            else
                if (bpsw)
                    sprintf (buf, "%s is Fermat, BPSW and Frobenius PRP! (P = %ld, Q = %ld, D = %ld)",quotient?string_rep:str , P, Q, D);
                else
                    sprintf (buf, "%s is Fermat, Lucas and Frobenius PRP! (P = %ld, Q = %ld, D = %ld)", quotient?string_rep:str, P, Q, D);

	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (tmp3);
	gwypfree (A);
	gwypfree (x);				// Clean up
	gwypfree (y);
	gwypfree (gwA);
	gwypfree (gw2);
        gwypfree (gwQ);

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
        gwypfree (gwQ);
//	gwypdone ();
	*res = FALSE;

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
            aborted = TRUE;
            sprintf (buf, ERRMSG5, checknumber, str);
            OutputBoth (buf);
//	    gwypdone ();
            _unlink (filename);
            if(IniGetInt(INI_FILE, (char*)"StopOnAbort", 0)) {
                IniWriteInt (INI_FILE, (char*)"PgenLine", IniGetInt(INI_FILE, (char*)"PgenLine", 0) + 1);	// Point on the next line
                return (FALSE);
            }
            else
                return (TRUE);
	}

//	gwypdone ();

/* Output a message saying we are restarting */
	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
		return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
            IniWriteInt(INI_FILE, (char*)"FFT_Increment", nbfftinc =  IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
            if (nbfftinc == maxfftinc)
                abonroundoff = TRUE;	// Don't accept any more Roundoff error.
            _unlink (filename);
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
	giant	tmp, tmp2;
	char	filename[20], buf[sgkbufsize+256], fft_desc[256], oldres64[17];
	long	write_time = DISK_WRITE_TIME * 60;
	int	echk, saving, stopping, zres;
	time_t	start_time, current_time;
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;

	Nlen = bitlen (N);
	tmp = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);
	tmp2 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);

/* Init, subtract 1 from N to compute a^(N-1) mod N */

	gtog (N, tmp2);
	iaddg (-1, tmp2);
        
// tmp2 == N-1 == t*2^s (for Proth numbers, t == k, s == n)
// t == (N-1)>>s ; s == firstone
        
	while (bitval (tmp2, firstone) == 0)	// position of first one bit in N-1
		firstone++;   // firstone == s
	nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.
	*res = TRUE;		/* Assume it is a probable prime */

/* Init filename */

	tempFileName (filename, 'z', N);

/* Allocate memory */

	x = gwypalloc ();
	y = gwypalloc ();
	gwypminusone = gwypalloc ();
	gwypone = gwypalloc ();

	itogwyp (1, gwypone);
        gianttogwyp (tmp2, gwypminusone);

/* Optionally resume from save file and output a message */
/* indicating we are resuming a test */

	if (fileExists (filename) && readFromFile (filename, &bit, x, NULL)) {
		char	fmt_mask[80];
		double	pct;
		pct = trunc_percent (bit * 100.0 / Nlen);
		sprintf (fmt_mask,
			 "Resuming probable prime test of %%s at bit %%ld [%%.%df%%%%]\n",
			 PRECISION);
		sprintf (buf, fmt_mask, quotient?string_rep:str, bit, pct);
		OutputStr (buf);
		if (verbose)
			writeResults (buf);
	}

/* Otherwise, output a message indicating we are starting test */

	else {
		gwypclear_timers ();	// Init. timers
 		if(!usingDivPhi_m) {
                    if (showdigits)
			sprintf (buf, "Starting probable prime test of %s (%d decimal digits)\n", quotient?string_rep:str, nbdg);
                    else
			sprintf (buf, "Starting probable prime test of %s\n", quotient?string_rep:str);
                    OutputStr (buf);
                }
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
//        strong = 0;
/* Init the title */
        if (strong)
            title ((char*)"Strong Fermat PRP test in progress...");
        else
            title ((char*)"Fermat PRP test in progress...");

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

		bitpos = Nlen-bit-1;

		if (bitval (tmp2, bitpos)) {
			gwypsetnormroutine (0, echk, 1);
		} else {
			gwypsetnormroutine (0, echk, 0);
		}

		if ((bit != lasterr_point) || !maxerr_recovery_mode[6]) {
                    if (cufftonly) {
                        gwypsquare (x);
                    }
                    else if (zp || generic) {
                        cuda_gwypsquare (x, 3);
                    }
                    else if((bit==1) || (it==0))  
                        {cuda_gwypsquare (x,1);it=1;}
                    else  if(bit != (lasterr_point-1) && (bit+25 < Nlen) && (bit > 25)) {
                        cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                    }
                    else {
                        cuda_gwypsquare (x,2);
                    }
                    care = FALSE;
                }
		else {
                    gwypsquare (x);
                    care = TRUE;
		}

		CHECK_IF_ANY_ERROR (x, (bit), Nlen, 6);

/* That iteration succeeded, bump counters */

		if (bit == lasterr_point)
                    saving = 1;	// Be sure to restart after this recovery iteration!
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
			sprintf (buf, fmt_mask, quotient?string_rep:str, bit, Nlen, pct);
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
			sprintf (buf, "%s interim residue %s at bit %ld\n", quotient?string_rep:str, res64, bit);
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
		if (strong) {
			if (bit == (Nlen-1)) {
				gwypsub3 (gwypone, x, y);
                                        // Test for x == 1 mod N
				zres = gwypiszero (y);
				if (zres == 1)	// success.
					break;
				gwypadd3 (gwypone, x, y);
                                        // Test for x == -1 mod N
				zres = gwypiszero (y);
				if (zres == 1)	// success.
					break;
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
                
	if (strong && (bit < Nlen)) {
		sprintf (buf, "%s is base %lu-Strong Fermat PRP! (%d decimal digits)", quotient?string_rep:str, a, nbdg);
	}
	else {

            gwyptogiant (x, tmp);
            modg ((quotient||(format == ABCDP))?M:N, tmp);

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
                if (usingDivPhi_m)
                        sprintf (buf, "%s does not divide %lu^%s^%d-1", quotient?string_rep:str, a, sgb, usingDivPhi_m);
                else
		if (IniGetInt (INI_FILE, (char*)"OldRes64", 1))
			sprintf (buf, "%s is not prime.  RES64: %s.  OLD64: %s", quotient?string_rep:str, res64, oldres64);
		else
			sprintf (buf, "%s is not prime.  RES64: %s", quotient?string_rep:str, res64);
            }
            else {
                if (strong) {
                    *res = FALSE;	/* Not a prime */
                    sprintf (buf, "%s is not prime, although base %lu-Fermat PSP!!", quotient?string_rep:str, a);
                }
                else if (usingDivPhi_m)
                    sprintf (buf, "%s Divides %lu^%s^%d-1", quotient?string_rep:str, a, sgb, usingDivPhi_m);
                else
                    sprintf (buf, "%s is base %lu-Fermat PRP! (%d decimal digits)", quotient?string_rep:str, a, nbdg);
            }
        }
        
/* Print known factors */

        if (string_rep_truncated) {
            char    *bigbuf;
            bigbuf = (char *) malloc (strlen (w->known_factors) + 100);
            if (bigbuf != NULL) {
                sprintf (bigbuf, "Known factors used for PRP test were: %s\n", w->known_factors);
                OutputStr (bigbuf);
                free (bigbuf);
            }
        }

/* Print results.  Do not change the format of this line as Jim Fougeron of */
/* PFGW fame automates his QA scripts by parsing this line. */

	gwypfree (tmp);
        gwypfree (tmp2);
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

	_unlink (filename);
	lasterr_point = 0;
	return (TRUE);

/* An error occured, sleep, then try restarting at last save point. */

error:
	gwypfree (tmp);
        gwypfree (tmp2);
	gwypfree (x);
	gwypfree (y);
	gwypfree (gwypone);
	gwypfree (gwypminusone);
	*res = FALSE;					// To avoid credit mesage...

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
            aborted = TRUE;
            sprintf (buf, ERRMSG5, checknumber, str);
            OutputBoth (buf);
            _unlink (filename);
            if(IniGetInt(INI_FILE, (char*)"StopOnAbort", 0)) {
                IniWriteInt (INI_FILE, (char*)"PgenLine", IniGetInt(INI_FILE, (char*)"PgenLine", 0) + 1);	// Point on the next line
                return (FALSE);
            }
            else
                return (TRUE);
	}

//	gwypdone ();

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
            return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
            IniWriteInt(INI_FILE, (char*)"FFT_Increment", nbfftinc =  IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
            if (nbfftinc == maxfftinc)
                abonroundoff = TRUE;	// Don't accept any more Roundoff error.
            _unlink (filename);
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
		*res = FALSE;
		return (TRUE);
	}
	
	Nlen = bitlen (N);
	nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.

/* Allocate memory */

	x = gwypalloc ();

	exponent = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);
	tmp = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);
        
/* Init, subtract 1 from N to compute a^(N-1) mod N */

        gtog (N, exponent);
	iaddg (-1, exponent);
	*res = TRUE;		/* Assume it is a prime */

/* Init filename */

	tempFileName (filename, 'z', N);


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
                    else  if(bit != (lasterr_point-1)&&(bit+25 < Nlen) && (bit > 25)) 
                        cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                    else
                        cuda_gwypsquare (x,2);
                    care = FALSE;
                }
		else {
                    gwypsquare_carefully (x);
                    care = TRUE;
		}

		CHECK_IF_ANY_ERROR (x, (bit), Nlen, 6);

/* That iteration succeeded, bump counters */

		if (bit == lasterr_point)
                    saving = 1;					// Be sure to restart after this recovery iteration!
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
                                gwypfree (exponent);
				gwypfree (tmp);
				gwypfree (x);
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
            aborted = TRUE;
            sprintf (buf, ERRMSG5, checknumber, str);
            OutputBoth (buf);
//		gwypdone ();
            _unlink (filename);
            if(IniGetInt(INI_FILE, (char*)"StopOnAbort", 0)) {
                IniWriteInt (INI_FILE, (char*)"PgenLine", IniGetInt(INI_FILE, (char*)"PgenLine", 0) + 1);	// Point on the next line
                return (FALSE);
            }
            else
                return (TRUE);
	}

//	gwypdone ();

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
            return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
            IniWriteInt(INI_FILE, (char*)"FFT_Increment", nbfftinc =  IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
            if (nbfftinc == maxfftinc)
                abonroundoff = TRUE;	// Don't accept any more Roundoff error.
            _unlink (filename);
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

	exponent = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);// Allocate memory for exponent
	tmp = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);	// Allocate memory for tmp
	tmp2 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);	// Allocate memory for tmp2
	tmp3 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);	// Allocate memory for tmp3

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
			 "Resuming Morrison prime test of %s at bit %%ld [%%.%df%%%%]\n",
			 str, PRECISION);
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
		if (showdigits)
			sprintf (buf, "Starting Morrison prime test of %s (%d decimal digits)\n", str, nbdg);
		else
			sprintf (buf, "Starting Morrison prime test of %s\n", str);
		OutputStr (buf);
		if (verbose)
			writeResults (buf);
		bit = 1;
		itogwyp (2, x);
		itogwyp (P, y);

	}

	itog (D, tmp3);    // Compute the inverse of D modulo N
//	invg (N, tmp3);
	gwypinvg (N, tmp3);
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
                            care = FALSE;
			}
			else {
                            gwypmul_carefully (y, x);
                            care = TRUE;
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
                            else  if(bit != (lasterr_point-1)&&(bit+26 < explen) && (bit > 26)) 
                                cuda_gwypsquare (y,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (y,2);
                            care = FALSE;
 			}
			else {
                            gwypsquare_carefully (y);
                            care = TRUE;
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
                            care = FALSE;
			}
			else {
                            gwypmul_carefully (x, y);
                            care = TRUE;
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
                            else  if(bit != (lasterr_point-1)&&(bit+26 < explen) && (bit > 26)) 
                                cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (x,2);
                            care = FALSE;
			}
			else {
                            gwypsquare_carefully (x);
                            care = TRUE;
			}
			CHECK_IF_ANY_ERROR(x, (bit), explen, 4)
		}

 /* That iteration succeeded, bump counters */

		if (bit == lasterr_point)
			saving = 1;					// Be sure to restart after this recovery iteration!
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

	care = TRUE; // All errors are now flagged as unrecoverable...

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
            aborted = TRUE;
            sprintf (buf, ERRMSG5, checknumber, str);
            OutputBoth (buf);
//		gwypdone ();
            _unlink (filename);
            if(IniGetInt(INI_FILE, (char*)"StopOnAbort", 0)) {
                IniWriteInt (INI_FILE, (char*)"PgenLine", IniGetInt(INI_FILE, (char*)"PgenLine", 0) + 1);	// Point on the next line
                return (FALSE);
            }
            else
                return (TRUE);
	}

//	gwypdone ();

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
            return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
            IniWriteInt(INI_FILE, (char*)"FFT_Increment", nbfftinc =  IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
            if (nbfftinc == maxfftinc)
                abonroundoff = TRUE;	// Don't accept any more Roundoff error.
            _unlink (filename);
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
	int	retval, a, gcdNa;
        char    buf[1024];

	if (bpsw)
		a = IniGetInt (INI_FILE, (char*)"FBase", 2);
	else
		a = IniGetInt (INI_FILE, (char*)"FBase", 3);
 	if (usingDivPhi_m) 
		a = 2;

//      Test if the candidate is coprime to the PRP base...
        
        gcdNa = gwypgcdui (N, a);
        if ((format != ABCDP) && (gcdNa != 1)) {
            sprintf (buf, "%s has a small factor : %d\n", str, gcdNa);
            OutputStr (buf);
            return (TRUE);
        }

//      Init work_unit used by the Gerbicz code
        
	w->k = k;
	w->b = b;
	w->n = n;
	w->c = c;
	w->prp_base = a;
	w->prp_residue_type = strong? PRP_TYPE_SPRP : PRP_TYPE_FERMAT;
	if (quotient) 
		w->prp_residue_type = PRP_TYPE_COFACTOR;
	else
		w->known_factors = NULL;
	do {
		gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
		gwypsetmaxmulbyconst (abs(a));
		if (!setupok (gwypsetup (k, b, n, c,(/*quotient||*/(format == ABCDP))?M:N), N, str, res)) {
			return TRUE;
		}

/* Do the PRP test */

		if (b == 2)
			retval = GerbiczTest (a, res, str, w);
		else
			retval = commonPRP (a, res, str);
                gwypdone ();
                if (recovering)
                    cufftonly = TRUE;   // restart while using fewer GPU code...
	} while (retval == -1);
        
/* Clean up and return */

//	gwypdone ();
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
		gwypsetmaxmulbyconst (abs(a));
		if (!setupok (gwypsetup (k, b, n, c, N), N, str, res)) {
			return TRUE;
		}

/* Do the Pocklington test */

		retval = commonCC1P (a, res, str);
                gwypdone ();
	} while (retval == -1);

/* Clean up and return */

//	gwypdone ();
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
		gwypsetmaxmulbyconst(abs(P));
		if (!setupok (gwypsetup (k, b, n, c, N), N, str, res)) {
			return TRUE;
		}

/* Do the Morrison test */

		retval = commonCC2P (P, res, str);
                gwypdone ();
	} while (retval == -1);

/* Clean up and return */

//	gwypdone ();
	return (retval);
}

/* Test if k*b^n+c (or a factor of it) is a Frobenius probable prime. */

int fastIsFrobeniusPRP (
	double	k,			/* k in k*b^n+c */
	unsigned long b,		/* b in k*b^n+c */
	unsigned long n,		/* n in k*b^n+c */
	signed long c,			/* c in k*b^n+c */
	char *str,
	int	*res)
{
	char	buf[sgkbufsize+256]; 
	int	retval;
	int32_t P = 3, Q = 0;
	long D, sign1, sign2, abs_d, sqrtabs_d;

/* Init. the relevant constants */

	if (bpsw) {
            P = 1;
            sign1 = N->n[0]&2? -1 : 1;	    // sign of kronecker(-1, N)
            for (abs_d = 5;abs_d < 2147483647;abs_d += 2) { // JP 27/10/22
                sign2 = abs_d&2 ? -1 : 1;   // requested sign of D
                if ((D = isLucasBaseQ (N, abs_d, ((sign2 == -1) && (sign1 == -1))? 1 : -1)) == TRUE) {
                    sqrtabs_d = (int)floor(sqrt((double)abs_d));
                    D = sign2*abs_d;
                    Q = (1-D)/4;
                    if ((Q == 1) || ((globalk == 1.0) && ispower(Q, globalb)) || (sqrtabs_d * sqrtabs_d == abs_d)) // abs(Q) was useless here... JP 25/03/22
                        continue;	// Avoid abs(Q) == 1 , k == 1 and abs(Q) power of b, or D perfect square...
                    else
                        break;
                }
                else if (D != FALSE) {
                    sprintf (buf, "%s has a small factor : %ld !!\n", str, abs(D));
                    OutputBoth (buf);
                    return (TRUE); 
                }
            }
            if (D == FALSE) {
                sprintf (buf, "Unable to further test %s ; it may be a perfect square !!\n", str);
                OutputBoth (buf);
                return (TRUE); 
            }
	}
	else {

		P = IniGetInt (INI_FILE, (char*)"PBase", 3);

		D = generalLucasBase (N, (uint32_t*)&P, (uint32_t*)&Q);
		if (D < 0) {
			if (D == -1)
				sprintf (buf, "Cannot compute D to test %s...\nThis is surprising, please, let me know that!!\nMy E-mail is jpenne@free.fr\n", str);
			else
				sprintf (buf, "%s has a small factor : %ld !!\n", str, abs(D));
			OutputBoth (buf);
			return (TRUE); 
		}
	}

	do {
		gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
		gwypsetmaxmulbyconst (bpsw?2:3); // JP 30/10/22
                if (recovering)
                    cufftonly = TRUE;   // 19/04/21
		if (!setupok (gwypsetup (k, b, n, c, N), N, str, res)) {
			return TRUE;
		}
		
/* Do the Frobenius PRP test */

		retval = commonFrobeniusPRP (P, Q, res, str);
                gwypdone ();

	} while (retval == -1);

/* Clean up and return */
	
//	gwypdone ();
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
	long D, sign1, sign2, abs_d, sqrtabs_d;

/* Setup the relevant constants. */

	if (bpsw) {
		P = 1;
		sign1 = N->n[0]&2? -1 : 1;			// sign of kronecker (-1, N)
		for (abs_d = 5;abs_d < 2147483647;abs_d += 2) { // JP 27/10/22
			sign2 = abs_d&2 ? -1 : 1;		// requested sign of D
			if ((D = isLucasBaseQ (N, abs_d, ((sign2 == -1) && (sign1 == -1))? 1 : -1)) == TRUE) {
				sqrtabs_d = (int)floor(sqrt((double)abs_d));
				D = sign2*abs_d;
				Q = (1-D)/4;
				if ((Q == 1) || ((globalk == 1.0) && ispower(Q, globalb)) || (sqrtabs_d * sqrtabs_d == abs_d)) // abs(Q) was useless here... JP 25/03/22
					continue;		// Avoid abs(Q) == 1 , k == 1 and abs(Q) power of b, or D perfect square...
				else
					break;
			}
			else if (D != FALSE) {
				sprintf (buf, "%s has a small factor : %ld !!\n", str, abs(D));
				OutputBoth (buf);
				return (TRUE); 
			}
		}
		if (D == FALSE) {
			sprintf (buf, "Unable to further test %s ; it may be a perfect square !!\n", str);
			OutputBoth (buf);
			return (TRUE); 
		}
	}
	else {
		P = IniGetInt (INI_FILE, (char*)"PBase", 3);

		D = generalLucasBase (N, (uint32_t*)&P, (uint32_t*)&Q);
		if (D < 0) {
			if (D == -1)
				sprintf (buf, "Cannot compute D to test %s...\nThis is surprising, please, let me know that!!\nMy E-mail is jpenne@free.fr\n", str);
			else
				sprintf (buf, "%s has a small factor : %ld !!\n", str, abs(D));
			OutputBoth (buf);
			return (TRUE); 
		}
	}

	do {
		gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
		gwypsetmaxmulbyconst (bpsw?2:3); // JP 30/10/22
                if (recovering)
                    cufftonly = TRUE;   // 19/04/21
		if (!setupok (gwypsetup_general_mod_giant (N), N, str, res)) {
			return TRUE;
		}

/* Do the Frobenius PRP test */

		retval = commonFrobeniusPRP (P, Q, res, str);
                gwypdone ();

	} while (retval == -1);

/* Clean up and return */

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
		gwypsetmaxmulbyconst (abs(a));
		if (!setupok (gwypsetup_general_mod_giant (M), M, str, res)) {
			return TRUE;
		}

/* Do the divisibility test */

		retval = isexpdiv (a, N, M, res);
                gwypdone ();

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
//	gwypdone ();

	return (retval);
}

void TestWieferich ()
{
	char str[10];
	int n, res;

	N = newgiant (4);

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
	int	retval, a, gcdNa;
        char    buf[1024];

/* Setup the gwypnum code */

	if (bpsw)
		a = IniGetInt (INI_FILE, (char*)"FBase", 2);
	else
		a = IniGetInt (INI_FILE, (char*)"FBase", 3);
 	if (usingDivPhi_m) 
		a = 2;
        
//      Test if the candidate is coprime to the PRP base...

        gcdNa = gwypgcdui (N, a);
        if ((format != ABCDP) && (gcdNa != 1)) {
            sprintf (buf, "%s has a small factor : %d\n", str, gcdNa);
            OutputStr (buf);
            return (TRUE);
        }
        

//      Init work_unit used by the Gerbicz code

	w->prp_base = a;
	w->prp_residue_type = strong? PRP_TYPE_SPRP : PRP_TYPE_FERMAT;
	if (quotient) 
		w->prp_residue_type = PRP_TYPE_COFACTOR;
	else
		w->known_factors = NULL;
	do {
		gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
		gwypsetmaxmulbyconst (abs(a));
		if (!setupok (gwypsetup_general_mod_giant (quotient||(format == ABCDP)?M:N), N, str, res)) {
			return TRUE;
		}

/* Do the PRP test */

		if (w->b == 2) {
                    retval = GerbiczTest (a, res, str, w);
                } 
		else {
                    retval = commonPRP (a, res, str);
                } 
                gwypdone ();
                if (recovering)
                    cufftonly = TRUE;   // restart while using fewer GPU code...
	} while (retval == -1);

/* Clean up and return */

//	gwypdone ();
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
		gwypsetmaxmulbyconst (abs(a));
		if (!setupok (gwypsetup_general_mod_giant (N), N, str, res)) {
			return TRUE;
		}

/* Do the Pocklington test */

		retval = commonCC1P (a, res, str);
                gwypdone ();
	} while (retval == -1);

/* Clean up and return */

//	gwypdone ();
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
		gwypsetmaxmulbyconst(abs(P));
		if (!setupok (gwypsetup_general_mod_giant (N), N, str, res)) {
			return TRUE;
		}

/* Do the Morrison test */

		retval = commonCC2P (P, res, str);
                gwypdone ();
	} while (retval == -1);

/* Clean up and return */

//	gwypdone ();
	return (retval);
}

/* Test if a small N is a probable prime. */
/* Compute b^(N-1) mod N */

int isProbablePrime (void)
{
	int	retval;
	giant	x;

	if (isone (N)) return (FALSE);
	x = newgiant (2*N->sign + 8);
	itog (IniGetInt (INI_FILE, (char*)"FBase", 3), x);
	powermodg (x, N, N);
	iaddg (-IniGetInt (INI_FILE, (char*)"FBase", 3), x);
	retval = isZero (x);
	gwypfree (x);
	return (retval);
}

int isPRPinternal (
	char *str,
        double dk, 
	unsigned long base,
	unsigned long n,
	int incr,
	int *res)
{
//	J.P. shadow        char buf[100];
	char	filename[20], buf[sgkbufsize+256]; // Lei - not need such long char
	unsigned long retval, fcontinue = FALSE;
	w->k = dk;
	w->b = base;
	w->n = n;
	w->c = incr;	// Transmit c to Gerbicz code!
	
	tempFileName (filename, 'L', N);
            // See if resuming a Lucas or Frobenius PRP test
	fcontinue = fileExists (filename);
	tempFileName (filename, 'F', N);
	fcontinue = fcontinue || fileExists (filename);
	if (dk >= 1.0) {
            if (fcontinue || IniGetInt(INI_FILE, (char*)"LucasPRPtest", 0)) {
                if (!fcontinue)
                    gwypclear_timers ();    // Init. timers
                retval = fastIsFrobeniusPRP (dk, base, n, incr, str, res);
            }
            else {
                retval = fastIsPRP (dk, base, n, incr, str, res);
                if (retval && *res && !fermat_only && !IniGetInt(INI_FILE, (char*)"FermatPRPtest", 0))
                    retval = fastIsFrobeniusPRP (dk, base, n, incr, str, res);
            }
	}
	else if (Nlen < 50) {
            *res = isProbablePrime();
            if (*res) {
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
                    gwypclear_timers ();    // Init. timers
			retval = slowIsFrobeniusPRP (str, res);
            }
            else {
                retval = slowIsPRP (str, res);
                if (retval && *res && !fermat_only && !IniGetInt(INI_FILE, (char*)"FermatPRPtest", 0))
                    retval = slowIsFrobeniusPRP (str, res);
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
			if (retval && *res && !fermat_only && !IniGetInt(INI_FILE, (char*)"FermatPRPtest", 0))
				retval = fastIsFrobeniusPRP (dk, smallbase, n, incr, str, res);
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
			if (retval && *res && !fermat_only && !IniGetInt(INI_FILE, (char*)"FermatPRPtest", 0))
				retval = slowIsFrobeniusPRP (str, res);
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
//#define ABCVARAS  18  // k, b, n, and c specified on each input line
//#define ABCVARAQS 19  // k, b, n, c and quotient specified on each input line
#define	ABCRU	20  // (10^n-1)/9 Repunits
#define	ABCGRU	21  // (b^n-1)/(b-1) Generalized Repunits

#define ABCGF	22  // ABC format for generalized Fermat numbers
#define ABCDN	23  // b^n-b^m+c format, m < n <= 2*m
#define ABCDNG	24  // General b^n-b^m+c format, m < n <= 2*m
#define ABCWFT	25  // Format used for Wieferich test
#define ABCWFS	26  // Format used for Wieferich search
#define ABCGPT	27  // Format used for General prime test (APRCL)
//#define ABCDP	28  // Format used for DivPhi()

int IsPRP (	    // General PRP test
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
		sprintf (str, "(%lu^%lu-1)/%lu", base, n_orig, base-1);
		gk = newgiant (1);
		itog (1, gk);
	}
	else if (!(format == ABCC || format == ABCK)) {
            gk = newgiant (strlen(sgk)/2 + 8);	// Allocate one byte per decimal digit + spares
            ctog (sgk, gk);			// Convert k string to giant
            gshiftleft (shift, gk);		// Shift k multiplier if requested
            gtoc (gk, sgk1, sgkbufsize);	// Updated k string
            if (mask & MODE_DUAL) {
                sprintf (str, "%lu^%lu%c%d", base, n_orig, incr < 0 ? '-' : '+', abs(incr));
            }
            else if (format != NPGAP) {   // Not MODE_AP
                if (!strcmp(sgk1, "1")) {
                    if (format == ABCVARAQS) {
                        char *p;
                        while ((p = strchr (sgd,'.'))!=NULL)
                            *p = '/';
                        sprintf (str, "%s^%lu%c%d", sgb, n_orig, incr < 0 ? '-' : '+', abs(incr));
                    }
                    else
                        sprintf (str, "%s^%lu%c%d", sgb, n_orig, incr < 0 ? '-' : '+', abs(incr));
                }
                else {
                    if (format == ABCVARAQS) {
                        char *p;
                        while ((p = strchr (sgd,'.'))!=NULL)
                            *p = '/';
                        sprintf (str, "%s*%s^%lu%c%d", sgk1, sgb, n_orig, incr < 0 ? '-' : '+', abs(incr));
                    }
                    else {
                        if ((n != 0) || (incr != 0))
                            sprintf (str, "%s*%s^%lu%c%d", sgk1, sgb, n_orig, incr < 0 ? '-' : '+', abs(incr));
                        else
                            sprintf (str, "%s", sgk1);
                    }
                }
            }
	}
	else {	// MODE_AP
		gk = newgiant ((n>>3)+8);
		itog (1, gk);		// Compute k multiplier
		gshiftleft (n-2, gk);	// Warning : here, n is exponent+1 !
		if (format == ABCK) {
			iaddg (1, gk);
			sprintf (str, "%s*2^%lu%c1 = (2^%lu+1)^2 - 2", sgk, n_orig, '-', n-1);
		}
		else {
			iaddg (-1, gk);
			sprintf (str, "%s*2^%lu%c1 = (2^%lu-1)^2 - 2", sgk, n_orig, '-', n-1);
		}
	}

	if ((gformat == ABCDN) || (gformat == ABCDNG)) {// Compute gk = gb^(n-m)-1
		bits = (unsigned long)(ndiff*log (base)/log (2));
		gk = newgiant ((bits >> 2) + 8);
		itog (base, gk);
		power (gk, ndiff);
		iaddg (-1, gk);
		sprintf (str, "%lu^%lu-%lu^%lu%c%d", base, n+ndiff, base, n_orig, incr < 0 ? '-' : '+', abs(incr));
	}
	bits = (unsigned long) ((n * log(base)) / log(2) + bitlen(gk)); 
	N =  newgiant ((bits >> 2) + 8);   // Allocate memory for N

 	if (format == ABCDP) {
            M =  newgiant ((bits >> 2) + 8);   // Allocate memory for M
 		globalk = dk = (double)gk->n[0]; 
 		fermat_only = TRUE;  // 25/04/21
 		gtog (gb, M);        // 26/04/21 compute the correct modulus!
                power (M, n);
                addg (M, M);
                iaddg (1, M);
 		for(usingDivPhi_m = n; usingDivPhi_m > 0; usingDivPhi_m--) {
 			sprintf (buf, "Starting test if %s divides Phi(%s^%d,2)\n", str, sgb, usingDivPhi_m);
			OutputStr (buf);
 			gtog (gb, N);
 			power (N, usingDivPhi_m);
 			iaddg (1, N);
 			retval = isPRPinternal (str, dk, base, n, incr, res);
 			if(!*res) break;
 		}
 		if(usingDivPhi_m<(int)n) {
			sprintf(buf, "Conclusion: %s Divides Phi(%s^%d,2)\n", str, sgb, usingDivPhi_m+1);
			OutputBoth (buf);
		}
 		fermat_only = usingDivPhi_m = 0; //  0 in place of FALSE 24/04/21 */
 		free (N);
                free (M);
 		free (gk);
 		return retval;
 	}
 	
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
		char factor[sgkbufsize+256], *p, *p2;
		gd = newgiant (strlen(sgd) + 8);	// Allocate one byte per decimal digit + spares
		gr = newgiant ((bits >> 2) + 8);	// Allocate memory for the remainder
		p = sgd;
                M =  newgiant ((bits >> 2) + 8);        // Allocate memory for M
                gtog (N, M);                            // keep M = N*known factors
                while (p != NULL) {
			strcpy (factor, p);             // copy the tail of the chain.
			if ((p2 = strchr (p,'/'))!=NULL) {// search for next factor.
                            factor[p2-p] = '\0';        // terminate the present factor.
                            p = &p2[1];                 // Point on next factor or end.
			}
			else
                            p = NULL;
			if (!isDigitString(factor)) {
				sprintf (buf,"invalid digit string : %s\n", factor);
				OutputBoth (buf);
				*res = FALSE;
				free (gr);
				free (gd);
				free (N);
				free (gk);
				return TRUE;
			}
			ctog (factor, gd);   // Convert divisor string to giant
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
                                quotient = FALSE;
				return TRUE;
			}
			else {
                            divg (gd, N);
			}
			if (p == NULL)
                            break;
		}
                w->prp_residue_type = PRP_TYPE_COFACTOR;
		w->known_factors = sgd;
		quotient = TRUE;
		strong = FALSE;
                    // Do a simple Fermat PRP test (not strong).
	}

                    /* Format the string representation of the test number */
                if (w->known_factors == NULL) {
                    strcpy (string_rep, str); // 02/07/20
                    string_rep_truncated = FALSE;
                } else {
                    if (strchr (str, '^') == NULL)
                        strcpy (string_rep, str);
                    else
			sprintf (string_rep, "(%s)", str);
                    if (strlen (w->known_factors) < 40) {
			char	*p;
			strcat (string_rep, "/");
			strcat (string_rep, w->known_factors);
			while ((p = strchr (string_rep, ',')) != NULL) *p = '/';
			string_rep_truncated = FALSE;
                    } else {
			strcat (string_rep, "/known_factors");
			string_rep_truncated = TRUE;
                    }
                }

		Nlen = bitlen (N); 
		klen = bitlen(gk);
                nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.

		if (klen > 53) {	// we must use generic reduction
			dk = 0.0;
		}
		else {			// we can use DWT ; compute the multiplier as a double
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

	strong = TRUE;		// Restore Strong Fermat PRP test

	if (format == ABCVARAQS) {
		gwypfree (gr);
		gwypfree (gd);
                gwypfree (M);
	}
	gwypfree (N);
	gwypfree (gk);
	return retval;
}

int gIsPRP (			// General PRP test
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
		gk = newgiant ((n>>3)+8);
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
                        quotient = FALSE;
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
		ctog (sgd, gd);				// Convert divisor string to giant
                M =  newgiant ((bits >> 2) + 8);        // Allocate memory for M
                gtog (N, M);                            // keep M = N*known factors
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
                        quotient = FALSE;
			return TRUE;
		}
		else {
			divg (gd, N);
		}
		w->prp_residue_type = PRP_TYPE_COFACTOR;
		w->known_factors = sgd;
		quotient = TRUE;
		strong = FALSE;
                    // Do a simple Fermat PRP test (not strong).
	}

                    /* Format the string representation of the test number */
                if (w->known_factors == NULL) {
                    strcpy (string_rep, str); // 02/07/20
                    string_rep_truncated = FALSE;
                } else {
                    if (strchr (str, '^') == NULL)
                        strcpy (string_rep, str);
                    else
			sprintf (string_rep, "(%s)", str);
                    if (strlen (w->known_factors) < 40) {
			char	*p;
			strcat (string_rep, "/");
			strcat (string_rep, w->known_factors);
			while ((p = strchr (string_rep, ',')) != NULL) *p = '/';
			string_rep_truncated = FALSE;
                    } else {
			strcat (string_rep, "/known_factors");
			string_rep_truncated = TRUE;
                    }
                }

		Nlen = bitlen (N); 
		klen = bitlen(gk);
                nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.

		if (klen > 53) {	// we must use generic reduction
			dk = 0.0;
		}
		else {			// we can use DWT ; compute the multiplier as a double
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
                gwypfree (M);
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
	ctog (sgk, gk);				// Convert k string to giant
	gshiftleft (shift, gk);			// Shift k multiplier if requested
	gtoc (gk, sgk1, sgkbufsize);		// Updated k string
		if (mask & MODE_DUAL) {
			sprintf (str, "%lu^%lu%c%d", base, n_orig, incr < 0 ? '-' : '+', abs(incr));
		}
		else
			if (!strcmp(sgk1, "1"))
				sprintf (str, "%lu^%lu%c%d", base, n_orig, incr < 0 ? '-' : '+', abs(incr));
			else
				sprintf (str, "%s*%lu^%lu%c%d", sgk1, base, n_orig, incr < 0 ? '-' : '+', abs(incr));

	bits = (unsigned long) ((n * log(base)) / log(2) + bitlen(gk)); 
	N =  newgiant ((bits >> 2) + 8);   // Allocate memory for N

//	Compute the number we are testing.

	itog (base, N);
	power (N, n);
	mulg (gk, N);
	iaddg (incr, N);

		klen = bitlen(gk);
                nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.

		if (klen > 53) {	// we must use generic reduction
			dk = 0.0;
		}
		else {			// we can use DWT ; compute the multiplier as a double
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

		if (klen > 53 || generic || !smallbase) {	// we must use generic reduction
			dk = 0.0;
		}
		else {						// we can use DWT ; compute the multiplier as a double
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

	tmp = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);
	tmp2 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);

/* Init, */


	gtog (exponent, tmp2);
//	uldivg (base, tmp2);	// tmp2 = exponent/base

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
//	invg (modulus, tmp);
	gwypinvg (modulus, tmp);
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
                            care = FALSE;
			}
			else {
                            gwypmul_carefully (y, x);
                            care = TRUE;
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
                            else  if(bit != (lasterr_point-1)&&(bit+26 < Nlen) && (bit > 26)) 
                                cuda_gwypsquare (y,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (y,2);
                            care = FALSE;
			}
			else {
                            gwypsquare_carefully (y);
                            care = TRUE;
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
                            care = FALSE;
			}
			else {
                            gwypmul_carefully (x, y);
                            care = TRUE;
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
                            else  if(bit != (lasterr_point-1)&&(bit+26 < Nlen) && (bit > 26)) 
                                cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (x,2);
                            care = FALSE;
			}
			else {
                            gwypsquare_carefully (x);
                            care = TRUE;
			}
			CHECK_IF_ANY_ERROR(x, (bit), Nlen, 4)
		}

 /* That iteration succeeded, bump counters */

		if (bit == lasterr_point)
                    saving = 1;	// Be sure to restart after this recovery iteration!
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

	care = TRUE;	// All following errors are considered unrecoverable...

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
//				gcdg (modulus, tmp);
				gwypgcdg (modulus, tmp);
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

	care = FALSE;// Reset the "unrecoverable" condition.
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
            aborted = TRUE;
            sprintf (buf, ERRMSG5, checknumber, str);
            OutputBoth (buf);
            _unlink (filename);
            if(IniGetInt(INI_FILE, (char*)"StopOnAbort", 0)) {
                IniWriteInt (INI_FILE, (char*)"PgenLine", IniGetInt(INI_FILE, (char*)"PgenLine", 0) + 1);	// Point on the next line
                return (FALSE);
            }
            else
                return (TRUE);
	}

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
            return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
            IniWriteInt(INI_FILE, (char*)"FFT_Increment", nbfftinc =  IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
            if (nbfftinc == maxfftinc)
                abonroundoff = TRUE;	// Don't accept any more Roundoff error.
            _unlink (filename);
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
		sprintf (str, "%lu^%lu-%lu^%lu%c%d", base, n_orig+ndiff, base, n_orig, incr < 0 ? '-' : '+', abs(incr));
	}
        else {
            gk = newgiant (strlen(sgk)/2 + 8);	// Allocate one byte per decimal digit + spares
            ctog (sgk, gk);		   // Convert k string to giant
            gshiftleft (shift, gk);	   // Shift k multiplier if requested
            gtoc (gk, sgk1, sgkbufsize);   // Updated k string
            if (!strcmp(sgk1, "1"))
		sprintf (str, "%lu^%lu%c%d", base, n_orig, (incr < 0) ? '-' : '+', abs(incr));
            else
		sprintf (str, "%s*%lu^%lu%c%d", sgk1, base, n_orig, (incr < 0) ? '-' : '+', abs(incr));
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
	M = newgiant ((bits >> 2) + 16);   // Allocate memory for M


//	Compute the number we are testing.

	itog (base, M);
	power (M, n-1);
        gtog (M, N);
        ulmulg (base, N);

	Nlen = bitlen (N); // Bit length of base^n

	mulg (gk, N);      // N is now gk*base^n
        mulg (gk, M);      // M = (N-incr)/base         

	iaddg (incr, N);   // N = gk*base^n+incr

	klen = bitlen(gk);

	if (klen > 53) {	// we must use generic reduction
		dk = 0.0;
	}
	else {			// we can use DWT ; compute the multiplier as a double
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
            gwypfree(M);
            return retval;
	}

	Nlen = bitlen (N); // Bit length of N
	findbpf (base);	   // Factorize the base
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
                        gwypfree(M);
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
//	M = newgiant ((bits >> 2) + 16);	// Allocate memory for M
	tmp = newgiant ((bits >> 2) + 16);	// Allocate memory for tmp
	gtog (N, tmp);
	iaddg (-1, tmp);			// tmp = N-1

	gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
	if (incr == +1) {
		gwypsetmaxmulbyconst (abs(a));
                gtog (M,tmp);                   // tmp = (N-1)/base
		explen = bitlen (tmp);
		if (!setupok (gwypsetup (dk, base, n, +1,N), N, str, res)) {
			gwypfree(gk);
			gwypfree(N);
			gwypfree (M);
			gwypfree (tmp);
//			*res = FALSE;		// Not proven prime...
			return TRUE; 
		}
		tmp2 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);// Allocate memory for tmp2
		tmp3 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);// Allocate memory for tmp3
	}
	else {
		gwypsetmaxmulbyconst (max(abs(a), abs(P)));
//		gtog (N, M);
//		iaddg (1, M);
		explen = bitlen (tmp);
		if (!setupok (gwypsetup (dk, base, n, -1, N), N, str, res)) {
			gwypfree(gk);
			gwypfree(N);
			gwypfree (M);
			gwypfree (tmp);
//			*res = FALSE;		// Not proven prime...
			return TRUE; 
		}
		tmp2 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);// Allocate memory for tmp2
		tmp3 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);// Allocate memory for tmp3
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
                    else  if(bit != (lasterr_point-1)&&(bit+25 < explen) && (bit > 25)) 
                        cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                    else
                        cuda_gwypsquare (x,2);
                    care = FALSE;
		}
		else {
                    gwypsquare_carefully (x);
                    care = TRUE;
		}

		CHECK_IF_ANY_ERROR (x, (bit), explen, 6);

/* That iteration succeeded, bump counters */

		if (bit == lasterr_point)
                    saving = 1;	// Be sure to restart after this recovery iteration!
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
		care = TRUE;	// All following errors are considered unrecoverable...
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
					gwypsetmaxmulbyconst (max(abs(a), abs(P)));
					if (!setupok (gwypsetup (dk, base, n, -1, N), N, str, res)) {
						gwypfree(gk);
						gwypfree(N);
						gwypfree (M);
						gwypfree (tmp);
						gwypfree (tmp2);
						gwypfree (tmp3);
                                                gwypdone ();    // 07/02/21
//						*res = FALSE;	// Not proven prime...
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
                            else  if(bit != (lasterr_point-1)&&bit != (explen-1)) 
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
                                sprintf (buf, "%s may be prime, but N divides %ld^((N-1)/%lu))-1, giving up after %lu restarts...\n", str, a, bpf[j], maxrestarts);
                                frestart = FALSE;
                                *res = FALSE;   // Not proven prime...
                            }
                            else {
                                sprintf (buf, "%s may be prime, but N divides %ld^((N-1)/%lu))-1, restarting with a=%lu\n", str, a, bpf[j], newa);
                                a = newa;
                                IniWriteInt (INI_FILE, (char*)"NRestarts", nrestarts);
						IniWriteInt (INI_FILE, (char*)"FermatBase", a);
						frestart = TRUE;
                            }
                        }
                        else {
                            iaddg (-1, tmp);
//                            gcdg (N, tmp);
                            gwypgcdg (N, tmp);
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
                    if (*res && !frestart)
                        sprintf (buf, "%s is prime! (%d decimal digits)", str, nbdg);
		}
	}
	if (!frestart) {
		gwypfree (N);
                gwypfree (M);
		gwypfree (gk);
	}
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
	return (TRUE);

/* An error occured, sleep, then try restarting at last save point. */

error:
//	gwypfree (M);
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (tmp3);
	gwypfree (x);
	gwypfree (y);
	gwypdone ();
	*res = FALSE;

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
            aborted = TRUE;
            sprintf (buf, ERRMSG5, checknumber, str);
            OutputBoth (buf);
            gwypfree (N);
            gwypfree (M);
            gwypfree (gk);
//            gwypdone ();
            _unlink (filename);
            if (IniGetInt(INI_FILE, (char*)"PRPdone", 0))
                IniWriteString(INI_FILE, (char*)"PRPdone", NULL);
            if(IniGetInt(INI_FILE, (char*)"StopOnAbort", 0)) {
                IniWriteInt (INI_FILE, (char*)"PgenLine", IniGetInt(INI_FILE, (char*)"PgenLine", 0) + 1);	// Point on the next line
                return (FALSE);
            }
            else
                return (TRUE);
	}

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) { 
//            gwypdone ();
            gwypfree (M);
            return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
            IniWriteInt(INI_FILE, (char*)"FFT_Increment", nbfftinc =  IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
            if (nbfftinc == maxfftinc)
                abonroundoff = TRUE;	// Don't accept any more Roundoff error.
            _unlink (filename);
	}
//        gwypdone ();    // JP 26/11/18
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

	tmp = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);
	tmp2 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);

/* Init, */


	gtog (exponent, tmp2);
//	divg (gb, tmp2);	// tmp2 = exponent/base

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
//	invg (modulus, tmp);
	gwypinvg (modulus, tmp);
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
                            care = FALSE;
			}
			else {
                            gwypmul_carefully (y, x);
                            care = TRUE;
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
                            else  if(bit != (lasterr_point-1)&&(bit+26 < Nlen) && (bit > 26)) 
                                cuda_gwypsquare (y,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (y,2);
                            care = FALSE;
			}
			else {
                            gwypsquare_carefully (y);
                            care = TRUE;
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
                            care = FALSE;
			}
			else {
                            gwypmul_carefully (x, y);
                            care = TRUE;
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
                            else  if(bit != (lasterr_point-1)&&(bit+26 < Nlen) && (bit > 26)) 
                                cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                            else
                                cuda_gwypsquare (x,2);
                            care = FALSE;
			}
			else {
                            gwypsquare_carefully (x);
                            care = TRUE;
			}
			CHECK_IF_ANY_ERROR(x, (bit), Nlen, 4)
		}

 /* That iteration succeeded, bump counters */

		if (will_try_larger_fft && (bit == lasterr_point))
                    saving = 1;
// Be sure to restart after this recovery iteration!
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

	care = TRUE;	// All following errors are considered unrecoverable...

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
//				gcdg (modulus, tmp);
				gwypgcdg (modulus, tmp);
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

	care = FALSE;// Reset the "unrecoverable" condition.
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
            aborted = TRUE;
            sprintf (buf, ERRMSG5, checknumber, str);
            OutputBoth (buf);
            _unlink (filename);
            if(IniGetInt(INI_FILE, (char*)"StopOnAbort", 0)) {
                IniWriteInt (INI_FILE, (char*)"PgenLine", IniGetInt(INI_FILE, (char*)"PgenLine", 0) + 1);	// Point on the next line
                return (FALSE);
            }
            else
                return (TRUE);
	}

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
            return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
            IniWriteInt(INI_FILE, (char*)"FFT_Increment", nbfftinc =  IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
            if (nbfftinc == maxfftinc)
                abonroundoff = TRUE;	// Don't accept any more Roundoff error.
            _unlink (filename);
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
	N =  newgiant ((bits >> 2) + 16);   // Allocate memory for N
	M = newgiant ((bits >> 2) + 16);   // Allocate memory for M

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
        
	gtog (gb, M);
	power (M, n-1);
        gtog (M, N);
        mulg (gb, N);

	Nlen = bitlen (N); // Bit length of base^n

	mulg (gk, N); 
        mulg (gk, M);      // M = (N-incr)/base         

	iaddg (incr, N);

	klen = bitlen(gk);

	if (klen > 53) {	// we must use generic reduction
		dk = 0.0;
	}
	else {			// we can use DWT ; compute the multiplier as a double
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
            gwypfree(M);
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
//	M = newgiant ((bits >> 2) + 16);	// Allocate memory for M
	tmp = newgiant ((bits >> 2) + 16);	// Allocate memory for tmp
	gtog (N, tmp);
	iaddg (-1, tmp);			// tmp = N-1

	gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
	if (incr == +1) {
            gwypsetmaxmulbyconst (abs(a));
//            divg (gb, tmp);		// tmp = (N-1)/base
//            gtog (gb, tmp);
//            power (tmp, n-1);
//            mulg (gk, tmp);             // tmp = (N-1)/base
            gtog (M,tmp);                 // tmp = (N-1)/base
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
            tmp2 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);// Allocate memory for tmp2
            tmp3 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);// Allocate memory for tmp3
	}
	else {
            gwypsetmaxmulbyconst (max(abs(a), abs(P)));
//            gtog (N, M);
//            iaddg (1, M);
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
            tmp2 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);// Allocate memory for tmp2
            tmp3 = newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16);// Allocate memory for tmp3
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

		if ((bit+25 < explen) && (bit > 25) && ((bit != lasterr_point) || !maxerr_recovery_mode[6])) {
                    if (cufftonly)
                        gwypsquare (x);
                    else if (zp || generic)
                        cuda_gwypsquare (x, 3);
                    else if(bit==1 || it==0)  
                        {cuda_gwypsquare (x,1);it=1;}
                    else  if(bit != (lasterr_point-1)&&(bit+25 < explen) && (bit > 25)) 
                        cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                    else
                        cuda_gwypsquare (x,2);
                    care = FALSE;
		}
		else {
                    gwypsquare_carefully (x);
                    care = TRUE;
		}

		CHECK_IF_ANY_ERROR (x, (bit), explen, 6);

/* That iteration succeeded, bump counters */

		if (will_try_larger_fft && (bit == lasterr_point))
                    saving = 1;	// Be sure to restart after this recovery iteration!
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
            care = TRUE;    // All following errors are considered unrecoverable...
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
					gwypsetmaxmulbyconst (max(abs(a),abs(P)));
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
                            else  if(bit != (lasterr_point-1)&&bit != (explen-1)) 
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
                                if (bpf[j] == 1) {
                                    gtoc(gbpf[j], bpfstring, strlen(sgb));
                                    sprintf (buf, "%s may be prime, but N divides %lu^((N-1)/%s))-1, giving up after %lu restarts...\n", str, a, bpfstring, maxrestarts);
                                }
                                else
                                    sprintf (buf, "%s may be prime, but N divides %lu^((N-1)/%lu))-1, giving up after %lu restarts...\n", str, a, bpf[j], maxrestarts);
                                frestart = FALSE;
                                *res = FALSE;   // Not proven prime...
                            }
                            else {
                                if (bpf[j] == 1) {
                                    gtoc(gbpf[j], bpfstring, strlen(sgb));
                                    sprintf (buf, "%s may be prime, but N divides %lu^((N-1)/%s))-1, restarting with a=%lu\n", str, a, bpfstring, newa);
                                }
                                else
                                    sprintf (buf, "%s may be prime, but N divides %lu^((N-1)/%lu))-1, restarting with a=%lu\n", str, a, bpf[j], newa);
                                a = newa;
                                IniWriteInt (INI_FILE, (char*)"NRestarts", nrestarts);
						IniWriteInt (INI_FILE, (char*)"FermatBase", a);
						frestart = TRUE;
                            }
                        }
                        else {
                            iaddg (-1, tmp);
//                            gcdg (N, tmp);
                            gwypgcdg (N, tmp);
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
                    care = FALSE;	
                    if (*res && !frestart)
                        sprintf (buf, "%s is prime! (%d decimal digits)", str, nbdg);
		}
	}
	if (!frestart) {
		gwypfree (N);
                gwypfree (M);
		gwypfree (gk);
	}
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
	return (TRUE);

/* An error occured, sleep, then try restarting at last save point. */

error:
//	gwypfree (M);
	gwypfree (tmp);
	gwypfree (tmp2);
	gwypfree (tmp3);
	gwypfree (x);
	gwypfree (y);
	gwypdone ();
	*res = FALSE;

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
            aborted = TRUE;
            sprintf (buf, ERRMSG5, checknumber, str);
            OutputBoth (buf);
            gwypfree (N);
            gwypfree (M);
            gwypfree (gk);
//            gwypdone ();
            _unlink (filename);
            if (IniGetInt(INI_FILE, (char*)"PRPdone", 0))
                IniWriteString(INI_FILE, (char*)"PRPdone", NULL);
            if(IniGetInt(INI_FILE, (char*)"StopOnAbort", 0)) {
                IniWriteInt (INI_FILE, (char*)"PgenLine", IniGetInt(INI_FILE, (char*)"PgenLine", 0) + 1);	// Point on the next line
                return (FALSE);
            }
            else
                return (TRUE);
	}

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) { 
//            gwypdone ();
            gwypfree (M);
            return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
            IniWriteInt(INI_FILE, (char*)"FFT_Increment", nbfftinc =  IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
            if (nbfftinc == maxfftinc)
                abonroundoff = TRUE;	// Don't accept any more Roundoff error.
            _unlink (filename);
	}
//        gwypdone ();
	goto restart;

}

/*
 	Primality testing of k*2^n-1 numbers with the Lucas Lehmer Riesel
	algorithm, using the gwypnums for fast multiplications and squarings.
	Second attempt for a full IBDWT version, using Colin Percival method improved
	by George Woltman.
	Jean Penne May 2004.
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
	long	retval;
	gwypnum	x, y; 
	giant	tmp; 
	char	filename[20], buf[sgkbufsize+256], str[sgkbufsize+256],
			sgk1[sgkbufsize], fft_desc[256]; 
	long	write_time = DISK_WRITE_TIME * 60; 
	int		echk, saving, stopping, v1, prp_res = 0, fftfermat = 0; 
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
			ddk = ninput /** log10 (ddk)*/;
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
                        nbllr_mallocs++;
                        itog (binput, gk);
			power (gk, ndiff);
			iaddg (-1, gk);
			sprintf (str, "%lu^%lu-%lu^%lu-1", binput, ninput+ndiff, binput, ninput);
		}
		else {
			gksize = 8*strlen(sgk) + idk;		// J. P. Initial gksize
			gk = newgiant ((gksize >> 2) + 8);	// Allocate space for gk
                        nbllr_mallocs++;
			ctog (sgk, gk);				// Convert k string to giant
		}
// Lei
		if (b_else != 1) {				// Compute the big multiplier
			gk1 = newgiant ((gksize>>2) + 8);
                        nbllr_mallocs++;
			itog (b_else, gk1);	
			power (gk1, ninput);
			mulg (gk1, gk);
			gwypfree (gk1);
                        nbllr_frees++;
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
                nbllr_mallocs++;
		itog (1, gk);					// Compute k multiplier
		gshiftleft (n-2, gk);				// Warning : here, n is exponent+1 !
		if (format == ABCK) {
			iaddg (1, gk);
			sprintf (str, "%s*2^%lu%c1 = (2^%lu+1)^2 - 2", sgk, n_orig, '-', n_orig-1);
		}
		else {
			iaddg (-1, gk);
			sprintf (str, "%s*2^%lu%c1 = (2^%lu-1)^2 - 2", sgk, n_orig, '-', n_orig-1);
		}
	}

	klen = bitlen(gk);			// Bit length ok k multiplier
	bits = n + klen;			// Bit length of N
	N =  newgiant ((bits >> 2) + 8);	// Allocate memory for N
//	N =  newgiant ((bits>>3) + 8);	        // Allocate memory for N incr 11/11/20
        nbllr_mallocs++;


//	Compute the number we are testing.

	gtog (gk, N);
	gshiftleft (n, N);
	iaddg (-1, N);

	Nlen = bitlen (N); 
	nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.

	if (klen > 53 || generic) {// we must use generic reduction. 12/06/20
		dk = 0.0;
	}
	else {	// we can use DWT ; compute the multiplier as a double
		dk = (double)gk->n[0];
		if (gk->sign > 1)
			dk += 65536.0*(double)gk->n[1];
		if (gk->sign > 2)
			dk += 65536.0*65536.0*(double)gk->n[2];
		if (gk->sign > 3)
			dk += 65536.0*65536.0*65536.0*(double)gk->n[3];
	}
	if ((klen > n) && (nbdg > 400)) {
		if ((format == ABCDN) || (format == ABCDNG))
			sprintf(buf, "2^%lu-1 > 2^%lu, so we can only do a PRP test for %s.\n", ndiff, n, str);
		else
			sprintf(buf, "%s > 2^%lu, so we can only do a PRP test for %s.\n", sgk, n, str);
                OutputBoth(buf);
		if ((format == ABCDN) || (format == ABCDNG))
			sprintf (str, "%lu^%lu-%lu^%lu-1", binput, ninput+ndiff, binput, ninput);
		else
			sprintf (str, "%s*%lu^%lu%c1", sgk, binput, ninput, '-');
                    // Number N to test, as a string
//BUG!		retval = isPRPinternal (str, dk, binput, ninput, -1, res);
		retval = isPRPinternal (str, dk, binput, n, -1, res);

		gwypfree(gk);
                nbllr_frees++;
		gwypfree(N);
                nbllr_frees++;
//                recovering = FALSE;       // 20/04/21
		return retval;
	}
	else if (klen > 53)        // k is a big integer
            dk = 0.0;

	if (!IniGetInt(INI_FILE, (char*)"Verify", 0) && !IniGetInt(INI_FILE, (char*)"PRPdone", 0) && (Nlen/klen < 10.0 || IniGetInt (INI_FILE, (char*)"ErrorChecking", 1)) && (nbdg > 400)) {
// We have better to do ad first a Gerbicz PRP test.
		strcpy (buf, str);
		if ((format == ABCDN) || (format == ABCDNG))
			sprintf (str, "%lu^%lu-%lu^%lu-1", binput, ninput+ndiff, binput, ninput);
		else
			sprintf (str, "%s*%lu^%lu%c1", sgk, binput, ninput, '-');     // Number N to test, as a string
				fermat_only = TRUE;
//              strong = FALSE;
                retval = isPRPinternal (str, dk, 2, n, -1, res);
				fermat_only = FALSE;

		if (!*res) {
			gwypfree(gk);
                        nbllr_frees++;
			gwypfree(N);
                        nbllr_frees++;
//                        recovering = FALSE;       // 20/04/21
			return retval;
		}
		IniWriteInt(INI_FILE, (char*)"PRPdone", 1);
                prp_res = *res;
		strcpy (str, buf);	// Lei
                fftfermat = FFTLEN;
	}

        
	k = gk->n[0];
	if(abs(gk->sign) == 2) {	// k is a "small" integer
		k += 65536*gk->n[1];
	}
	else if (abs(gk->sign) > 2)    {
		k = 0;			// to indicate that k is a big integer.
		dk = 0.0;
        }

restart: 

	*res = TRUE;		/* Assume it is prime */ 

	tempFileName (filename, 'z', N); //cuda
	if (fileExists (filename) && readFFTLENFromFile (filename, &j, x, NULL));        
	gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
        if (dk == 0.0) {    // generic modular reduction needed...
            if (!setupok (gwypsetup_general_mod_giant (N), N, str, res)) {
//                gwypdone ();
                gwypfree(gk);
                nbllr_frees++;
                gwypfree(N);
                nbllr_frees++;
//	        *res = FALSE;
                return TRUE;
            }
        }	// end dk == 0.0
//	else if (!setupok (gwypsetup (dk, binput, ninput, -1, N), N, str, res)) {    BUG!
	else if (!setupok (gwypsetup (dk, 2, n, -1, N), N, str, res)) {
            // JP 18/10/22
		gwypfree(gk);
                nbllr_frees++;
		gwypfree(N);
                nbllr_frees++;
//		*res = FALSE;
                return TRUE;
	}

	x = gwypalloc (); 
        nbllr_mallocs++;
	y = gwypalloc ();
        nbllr_mallocs++;
 
	tmp =  newgiant (2*FFTLEN*sizeof(double)/sizeof(short) + 16); 
        nbllr_mallocs++;
        
	if (fftfermat && (FFTLEN < fftfermat)) { // J.P. 22/05/22
            will_try_larger_fft = TRUE;
            goto error;
        }
        
	last = n-1;

 	gwypsetnormroutine (0, ERRCHK, 0); 

/* Init filename */ 

	//cuda tempFileName (filename, 'z', N); 
 
/* Init the title */ 
 
	title ((char*)"L.L.R. prime test in progress...");
 
/* Optionally resume from save file and output a message */ 
/* indicating we are resuming a test */ 
 
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
                                nbllr_frees++;
				gwypfree(gk);
                                nbllr_frees++;
				gwypfree(N);
                                nbllr_frees++;
				gwypfree (x); 
                                nbllr_frees++;
				gwypfree (y);
                                nbllr_frees++;
				gwypdone();
				*res = FALSE;
				gwypend_timer (1); 
				return(TRUE);
			}
                        if (showdigits)
                                sprintf (buf, "Starting Lucas Lehmer prime test of %s (%d decimal digits)\n", str, nbdg);
                        else
                                sprintf (buf, "Starting Lucas Lehmer prime test of %s\n", str);
                        if (verbose)
                            OutputBoth(buf);
                        else
                            OutputStr (buf); 
			gwypfft_description (fft_desc);
			sprintf (buf, "%s\n", fft_desc);
                        if (verbose)
                            OutputBoth(buf);
                        else
                            OutputStr (buf); 
			v1 = 4;
			itogwyp (v1, x);
			gwypclear_timers ();		// Init. timers
			gwypstart_timer (0); 
			gwypstart_timer (1); 
			time (&start_time); 
			goto MERSENNE;
	    }

	    filename[0] = 'u';
	    if ((v1 = gen_v1(gk, n, general, eps2, vdebug)) < 0) {
			if (v1 == -1)
				sprintf (buf, "Cannot compute V1 to test %s...\nThis is surprising, please, let me know that!!\nMy E-mail is jpenne@free.fr\n", str);
			else
				sprintf (buf, "%s has a small factor : %d !!\n", str, abs(v1));
			OutputBoth (buf); 
			gwypfree (tmp);
                        nbllr_frees++;
			gwypfree(gk);
                        nbllr_frees++;
			gwypfree(N);
                        nbllr_frees++;
			gwypfree (x); 
                        nbllr_frees++;
			gwypfree (y);
                        nbllr_frees++;
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
                            nbllr_frees++;
                            gwypfree(gk);
                            nbllr_frees++;
                            gwypfree(N);
                            nbllr_frees++;
                            gwypfree (x); 
                            nbllr_frees++;
                            gwypfree (y);
                            nbllr_frees++;
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
			gwypsetnormroutine (0, 1, 0);
			gwypsetaddin (-2);
			if ((1 != lasterr_point) || !maxerr_recovery_mode[0]) {
                            if (cufftonly)
				gwypsquare (y);
                            else
                                cuda_gwypsquare (y,3);
                            care = FALSE;
			}
			else {
                            gwypsquare_carefully (y);
			    care = TRUE;
			}
			CHECK_IF_ANY_ERROR(y, 1, klen, 0)
			if (1 == lasterr_point)
                            saving = 1;	// Be sure to restart after this recovery iteration!
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
                                    care = FALSE;
				}
				else {
                                    gwypmul_carefully (y, x);
                                    care = TRUE;
				}
				CHECK_IF_ANY_ERROR(x, (index), klen, 1)
				gwypsetaddin (-2);
				if ((index != lasterr_point) || !maxerr_recovery_mode[2]) {
                                    if (cufftonly)
					gwypsquare (y);
                                    else
                                        cuda_gwypsquare (y,3);
                                    care = FALSE;
				}
				else {
                                    gwypsquare_carefully (y);
                                    care = TRUE;
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
                                    care = FALSE;
				}
				else {
                                    gwypmul_carefully (x, y);
                                    care = TRUE;
				}
				CHECK_IF_ANY_ERROR(y, (index), klen, 3)
 				gwypsetaddin (-2);
				if ((index != lasterr_point) || !maxerr_recovery_mode[4]) {
                                    if (cufftonly)
					gwypsquare (x);
                                    else
                                        cuda_gwypsquare (x,3);
                                    care = FALSE;
				}
				else {
                                    gwypsquare_carefully (x);
                                    care = TRUE;
				}
				CHECK_IF_ANY_ERROR(x, (index), klen, 4)
			}

			if (will_try_larger_fft && (index == lasterr_point))
                            saving = 1;					// Be sure to restart after this recovery iteration!

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
                                    nbllr_frees++;
                                    gwypfree(gk);
                                    nbllr_frees++;
                                    gwypfree(N);
                                    nbllr_frees++;
                                    gwypfree (x); 
                                    nbllr_frees++;
                                    gwypfree (y);
                                    nbllr_frees++;
                                    gwypdone();
                                    return (FALSE); 
				}
			} 
	    }

		gwypsetaddin (-v1);
		if ((klen != lasterr_point) || !maxerr_recovery_mode[5]) {
                    if (cufftonly)
                        gwypmul (y, x);
                    else
                        cuda_gwypmul (y, x, 3);
                    care = FALSE;
                }
		else {
                    gwypmul_carefully (y, x);
                    care = TRUE;
		}
		CHECK_IF_ANY_ERROR(x, klen, klen, 5)

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

	    sprintf (buf, "Starting Lucas-Lehmer loop..."); 
	    OutputStr (buf); 
            LineFeed();
MERSENNE:
            j = 1;
	} 

/* Do the Lucas Lehmer Riesel Prime test */ 

	ReplaceableLine (1);	/* Remember where replacable line is */  
	iters = 0; 
	gwypsetaddin (-2);
	it = 0; //cuda
	while (j<last) { 

/* Error check the first and last 50 iterations, before writing an */ 
/* intermediate file (either user-requested stop or a */ 
/* 30 minute interval expired), and every 128th iteration. */ 
		stopping = stopCheck (); 
		echk = stopping || ERRCHK || (j <= 30) || (j >= last - 30); 
		if (((j & 127) == 0) || (j == 1) || (j == (lasterr_point-1))) {
			echk = 1;
			time (&current_time);
			saving = ((current_time - start_time > write_time) || (j == 1) || (j == (lasterr_point-1)));
		} else
			saving = 0;

/* Process this iteration */ 

		gwypsetnormroutine (0, echk, 0);

                if ((j > 30) && (j < last - 30) && ((j != lasterr_point) || !maxerr_recovery_mode[6])) {
                    if (cufftonly)
                        gwypsquare (x);
                    else if (zp || generic)
                        cuda_gwypsquare (x, 3);
                    else if(j==1 || it==0)
                        {cuda_gwypsquare (x,1);it=1;}
                    else if(j != (lasterr_point-1)&&(j > 31) && (j < last - 31))
                        cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                    else
                        cuda_gwypsquare (x,2);
                    care = FALSE;
		}
		else {
                    gwypsquare_carefully (x);
                    care = TRUE;
		}
		CHECK_IF_ANY_ERROR(x, j, last, 6)
                
		if (will_try_larger_fft && (j == lasterr_point))
                    saving = 1;	// Be sure to restart after this recovery iteration!
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
                            nbllr_frees++;
                            gwypfree(gk);
                            nbllr_frees++;
                            gwypfree(N);
                            nbllr_frees++;
                            gwypfree (x); 
                            nbllr_frees++;
                            gwypfree (y);
                            nbllr_frees++;
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
        if ((prp_res == TRUE) && (*res == FALSE)) {      // Probably a false negative result...
            sprintf (buf, "**** IT IS PROBABLY A FALSE NEGATIVE, RESTARTING USING FEWER GPU COMPUTING... ****\n");
            OutputBoth (buf); 
            prp_res = FALSE;
            gwypfree (tmp);
            nbllr_frees++;
            gwypfree (x); 
            nbllr_frees++;
            gwypfree (y); 
            nbllr_frees++;
            gwypdone (); 
            IniWriteInt(INI_FILE, (char*)"FFT_Increment", nbfftinc = IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
            cufftonly = TRUE;
            filename[0] = 'z';
            _unlink (filename); 
            goto restart;
        }
	gwypfree (tmp);
        nbllr_frees++;
	gwypfree(gk);
        nbllr_frees++;
	gwypfree(N);
        nbllr_frees++;
	gwypfree (x); 
        nbllr_frees++;
	gwypfree (y);
        nbllr_frees++;
	gwypdone (); 
	filename[0] = 'z';
	_unlink (filename); 
	if (IniGetInt(INI_FILE, (char*)"PRPdone", 0))
		IniWriteString(INI_FILE, (char*)"PRPdone", NULL);
	lasterr_point = 0;
        g_fftlen = 0;
	return (TRUE); 
 
/* An error occured, sleep, then try restarting at last save point. */ 

error:
	gwypfree (tmp);
        nbllr_frees++;
	gwypfree (x); 
        nbllr_frees++;
	gwypfree (y); 
        nbllr_frees++;
	gwypdone (); 
	*res = FALSE;

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
		aborted = TRUE;
		sprintf (buf, ERRMSG5, checknumber, str);
		OutputBoth (buf);
		gwypfree(gk);
                nbllr_frees++;
		gwypfree(N);
                nbllr_frees++;
		filename[0] = 'u';
		_unlink (filename);
		filename[0] = 'z';
//		gwypdone ();
		_unlink (filename); 
		if (IniGetInt(INI_FILE, (char*)"PRPdone", 0))
			IniWriteString(INI_FILE, (char*)"PRPdone", NULL);
                g_fftlen = 0;
		if(IniGetInt(INI_FILE, (char*)"StopOnAbort", 0)) {
			IniWriteInt (INI_FILE, (char*)"PgenLine", IniGetInt(INI_FILE, (char*)"PgenLine", 0) + 1);	// Point on the next line
			return (FALSE);
		}
		else
                    return (TRUE);
	}

/* Output a message saying we are restarting */ 
 
	if (sleep5) OutputBoth (ERRMSG2);
        if (!fftfermat)
            OutputBoth (ERRMSG3); 
 
/* Sleep five minutes before restarting */ 
 
	if (sleep5 && ! SleepFive ()) {
//            gwypdone ();
            return (FALSE); 
	}

/* Restart */ 
 
	if (will_try_larger_fft) {
            IniWriteInt(INI_FILE, (char*)"FFT_Increment", nbfftinc = IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
            if (nbfftinc == maxfftinc)
                abonroundoff = TRUE; // Don't accept any more Roundoff error.
            _unlink (filename);
	}
	will_try_larger_fft = FALSE;
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
        nbllr_mallocs++;

	ctog (sgk, gk);				// Convert k string to giant

	if (shift > 0) {
		gshiftleft (shift, gk);		// Shift k multiplier if requested
		gtoc (gk, sgk1, sgkbufsize);	// Updated k string
	}
	else
		strcpy (sgk1, sgk);

	sprintf (str, "%s*2^%lu%c1", sgk1, n_orig, '-');	// Number N to test, as a string

	bits = n + bitlen(gk);				// Bit length of N
	N =  newgiant ((bits>>3) + 8);		// Allocate memory for N
        nbllr_mallocs++;

//	Compute the number we are testing.

	itog (1, N);
	gshiftleft (n, N);
	mulg (gk, N); 
	iaddg (-1, N);
	retval = slowIsWieferich (str, res);
	gwypfree (gk);
        nbllr_frees++;
	gwypfree (N);
        nbllr_frees++;
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

	sprintf (str, "%s*2^%lu%c1", sgk1, n_orig, '+');	// Number N to test, as a string

	bits = n + bitlen(gk);				// Bit length of N
	N =  newgiant ((bits>>3) + 8);		// Allocate memory for N

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
	unsigned long gksize; 
	unsigned long bits; 
	long	a, retval;
	char	buf[sgkbufsize+256], 
		str[sgkbufsize+256], sgk1[sgkbufsize]; 
	long	write_time = DISK_WRITE_TIME * 60; 
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
		gksize = 8*strlen(sgk) + idk;	// J. P. Initial gksize
		gk = newgiant ((gksize)  + 8);  // Allocate space for gk incr 11/11/20
		ctog (sgk, gk);			// Convert k string to giant
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
			sprintf (str, "%s*2^%lu%c1", sgk1, n_orig, '+');	// Number N to test, as a string


	bits = n + bitlen(gk);	        // Bit length of N
	N =  newgiant ((bits>>3) + 8);	// Allocate memory for N incr 11/11/20

//	Compute the number we are testing.

	itog (1, N);
	gshiftleft (n, N);
	mulg (gk, N); 
	iaddg (1, N);

//	gk must be odd for the Proth test, so, adjust gk and n if necessary.

	while (bitval(gk, 0) == 0) {
	    gshiftright (1, gk);	// update k as a giant
	    n++;							// update the exponent
	}

	Nlen = bitlen (N); 
	klen = bitlen(gk);
	nbdg = gnbdg (N, 10);	// Compute the number of decimal digits of the tested number.
        
	if (klen > 53 || generic) {    // we must use generic reduction
		dk = 0.0;
	}
	else {	// we can use DWT ; compute the multiplier as a double
		dk = (double)gk->n[0];
		if (gk->sign > 1)
			dk += 65536.0*(double)gk->n[1];
		if (gk->sign > 2)
			dk += 65536.0*65536.0*(double)gk->n[2];
		if (gk->sign > 3)
			dk += 65536.0*65536.0*65536.0*(double)gk->n[3];
	}
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
// BUG!               retval = isPRPinternal (str, dk, binput, ninput, 1, res);
                retval = isPRPinternal (str, dk, 2, n, 1, res);
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
        
	w->k = dk;		// init work_unit for the Gerbicz test.
	w->n = n;
	w->b = 2;
	w->c = 1;
	w->prp_base = a;
	w->prp_residue_type = PROTH_TYPE;

	do {
//restart:
                gwypsetmaxmulbyconst (abs(a));
		if (dk == 0.0) {  // generic modular reduction needed...
			gwypset_larger_fftlen_count((char)IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
			if (!setupok (gwypsetup_general_mod_giant (N), N, str, res)) {
                                gwypdone ();
				free(gk);
				free(N);
//				*res = FALSE;
				return TRUE;
			}
		}			// end dk == 0.0
		else {
//restart2:
			gwypset_larger_fftlen_count((char)IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));
//			if (!setupok (gwypsetup (dk, binput, ninput, +1, N), N, str, res)) {   BUG!
			if (!setupok (gwypsetup (dk, 2, n, +1, N), N, str, res)) {// JP 18/10/22
                                gwypdone ();
				free(gk);
				free(N);
//				*res = FALSE;
				return TRUE;
			}
                }

/* Do the Proth test */

			retval = GerbiczTest (a, res, str, w);
                        gwypdone ();
                        if (recovering)
                            cufftonly = TRUE;   // restart while using fewer GPU code...
	} while (retval == -1);

/* Clean up and return */
//                gwypdone ();
		free (gk);
		free (N);
		return (retval);
} 

#define LOWFACTORLIMIT 10000	// To factor lower exponent candidates is not useful...
int res1, res2;
long	a;

void normal_addg (giant, giant);

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

// We have now x = a^(2^(p-1)) modulo M and y = a^(2^((p-1)/2)) modulo M.
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

Jean Penne March 30 2006

****************************************************************************************************************/

int isGMNP ( 
	char *sgk,
	unsigned long n,
	int	*res) 
{ 
	unsigned long bits; 
	giant	tmp, tmp2, apow4; 
	char	buf[sgkbufsize+256], 
		str[sgkbufsize+256], strp[sgkbufsize+256]; 
	int	retval, sign; 
	double dk;

	if (!isPrime (n) || n == 2) {
		sprintf (buf, "Gaussian-Mersenne prime test not done because %lu is not an odd prime.\n", n); 
		OutputBoth (buf); 
		*res = FALSE;
		return(TRUE);
	}

	sign = (((n&7) == 3) || ((n&7) == 5))? 1 : 0;	// 1 if positive, 0 if negative
	sprintf (str, "2^%lu%c%s+1",  n_orig, (sign) ? '+' : '-', sgk);	// Number N to test, as a string
	sprintf (strp, "(2^%lu%c%s+1)/5",  n_orig, (sign) ? '-' : '+', sgk);	// Number N' to test, as a string

	bits = 2*n+1;							// Bit length of M = N*N'
	M = newgiant ((bits>>1) + 8);	// Allocate memory for M = N*N'
	N = newgiant ((bits>>2) + 8);	// Allocate memory for N
	NP = newgiant ((bits>>2) + 8);	// Allocate memory for N'
	gk = newgiant ((bits>>2) + 8);	// Allocate memory for gk

// gk is the multiplier when N is written as gk*2^exponent + 1
// N = 2^n + s*2^((n+1)/2) + 1 = (2^((n-1)/2) + s)*2^((n+1)/2) + 1
// So, gk = 2^((n-1)/2)+s and exponent = (n+1)/2, with s = +1 or -1
// It is only used to compute the Proth base.

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
		return (TRUE); 
	}

 	dk = 1.0;	// k == 1 for the modulo M = N*NP

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
		return(TRUE);
	}


//restart:

// 	gwypsetmaxmulbyconst (abs(a));

//restart:

/* Assume intermediate results of the length of N*N'. */ 

	*res = TRUE;						/* Assume it is a prime */ 

	gwypset_larger_fftlen_count(IniGetInt(INI_FILE, (char*)"FFT_Increment", 0));

        if (!setupok (gwypsetup (dk, 2, 2*n, +1, M), M, str, res)) { 	// Setup the DWT mode
//		*res = res1 = res2 = FALSE;
		gwypfree(gk);
		gwypfree(N);
		gwypfree(NP);
		gwypfree(M);
		return TRUE;
	}


/* More initializations... */

	tmp = newgiant(4*FFTLEN*sizeof(double)/sizeof(short) + 16);
	tmp2 = newgiant(4*FFTLEN*sizeof(double)/sizeof(short) + 16);
//	tmp3 = newgiant(2*FFTLEN*sizeof(double)/sizeof(short) + 16);
	apow4 = newgiant(32);
	itog (a, apow4);
	smulg ((unsigned short)a, apow4);
	smulg ((unsigned short)a, apow4);
	smulg ((unsigned short)a, apow4);

	w->prp_residue_type = GMN_TYPE;
	w->known_factors = NULL;
	w->k = dk;
	w->n = n;
	w->b = 2;
	w->c = 1;
        
	do {
		retval = GerbiczTest (a, res, str, w);
                gwypdone();
                if (recovering)
                    cufftonly = TRUE;   // restart while using fewer GPU code...
	}	while (retval == -1);
        
//        gwypdone();
	if (retval == FALSE)
		goto EXIT;

	clearline (100);

	gtog (ps.gx, tmp);
	modg (M, tmp);
	gtog (ps.gy, tmp2);
	modg (M, tmp2);
	if (sign) {
		mulg (tmp2, tmp);
		modg (N, tmp);
		iaddg (1, tmp);	// Compute the (unnormalized) residue
	}
	else {
		gwypinvg (N, tmp2);   // 19/04/21 residue must match with LLR 3.8.24
		mulg (tmp2, tmp);
		modg (N, tmp);
		iaddg (1, tmp);
	}
	
/* See if we've found a Proth prime.  If not, format a 64-bit residue. */

//	if (resaprcl1 != 2)	{
// Exclude this output if previous APRCL positive result
        if (gcompg (N, tmp) != 0) {
                res1 = FALSE;	/* Not a prime */
            if (abs(tmp->sign) < 2)
                // make a 64 bit residue correct !!
                sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
            else if (abs(tmp->sign) < 3)
                sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
            else if (abs(tmp->sign) < 4)
                sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
            else
                sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
            }
            else {
                res1 = TRUE;
            }


/* Print results.  Do not change the format of this line as Jim Fougeron of */
/* PFGW fame automates his QA scripts by parsing this line. */

            if (res1) {
#if defined (_CONSOLE)
                sprintf (buf, "%s is prime! (%d decimal digits)", str, nbdg1);
#else
//			if ((resaprcl2 != 1) && (resaprcl2 != 2))
                sprintf (buf, "%s is prime! (%d decimal digits)\n", str, nbdg1);
//			else
//				sprintf (buf, "%s is prime! (%lu decimal digits)", str, nbdg1);
#endif
            }
            else {
//			if ((resaprcl2 != 1) && (resaprcl2 != 2))
                sprintf (buf, "%s is not prime.  Proth RES64: %s\n", str, res64);
//			else
//				sprintf (buf, "%s is not prime.  Proth RES64: %s", str, res64);
            }

#if defined(WIN32) && !defined(_CONSOLE)

//		if ((resaprcl2 != 1) && (resaprcl2 != 2))
            OutputBoth (buf);				// Avoid a double display...

#else

            clearline(100);

            if (res1) {
//			if ((resaprcl2 != 1) && (resaprcl2 != 2)) {
#if defined (_CONSOLE)
                hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
                // Access to Console attributes
                SetConsoleTextAttribute(hConsole, BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_RED);
                OutputBoth(buf);    // Avoid a double display...
                SetConsoleTextAttribute(hConsole, FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED);
                OutputBoth("\n");
#else
                OutputStr((char*)"\033[7m");
                OutputBoth(buf);    // Avoid a double display...
                OutputStr((char*)"\033[0m");
#endif
//			}
            }
            else {
//			if ((resaprcl2 != 1) && (resaprcl2 != 2))
                OutputBoth(buf);    //Avoid a double display...
            }

//		if ((resaprcl2 != 1) && (resaprcl2 != 2))
            sprintf (buf, "  Time : "); // Avoid a double display...

#endif
//	}									// End Exclude...
	gtog (ps.gx, tmp);
	modg (M, tmp);
	gtog (ps.gy, tmp2);
	modg (M, tmp2);

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

//	if ((resaprcl2 != 1) && (resaprcl2 != 2)) {
// Exclude this output if previous APRCL positive result

        if (gcompg (tmp2, tmp) != 0) {
            subg (tmp2, tmp);
            res2 = FALSE;   /* Not a prime */
            if (abs(tmp->sign) < 2)
                // make a 64 bit residue correct !!
                sprintf (res64, "%04X%04X%04X%04X", 0, 0, 0, tmp->n[0]);
            else if (abs(tmp->sign) < 3)
                sprintf (res64, "%04X%04X%04X%04X", 0, 0, tmp->n[1], tmp->n[0]);
            else if (abs(tmp->sign) < 4)
                sprintf (res64, "%04X%04X%04X%04X", 0, tmp->n[2], tmp->n[1], tmp->n[0]);
            else
                sprintf (res64, "%04X%04X%04X%04X", tmp->n[3], tmp->n[2], tmp->n[1], tmp->n[0]);
        }
        else
            res2 = TRUE;


/* Print results.  Do not change the format of this line as Jim Fougeron of */
/* PFGW fame automates his QA scripts by parsing this line. */


        if (res2)
            sprintf (buf, "%s is %lu-PRP! (%d decimal digits)", strp, a, nbdg2);
        else
            sprintf (buf, "%s is not prime.  RES64: %s", strp, res64);

#ifdef WIN32

        sprintf (buf+strlen(buf), "  Time: ");

#else

        if (res2) {
#if defined (_CONSOLE)
            hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
            // Access to Console attributes
            SetConsoleTextAttribute(hConsole, BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_RED);
            OutputBoth(buf);
            SetConsoleTextAttribute(hConsole, FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED);
#else
            OutputStr((char*)"\033[7m");
            OutputBoth(buf);
            OutputStr((char*)"\033[0m");
#endif
		}
            else
                OutputBoth (buf);
            sprintf (buf, "  Time: ");

#endif

//	}								// End Exclude...
//	else
//		sprintf (buf+strlen(buf), "  Time: ");

/* Output the final timings */

//	if (!(resaprcl1 == 2) || !((resaprcl2 == 2) || (resaprcl2 == 1))) {
            gwypend_timer (2);
            gwypwrite_timer (buf+strlen(buf), 2, TIMER_CLR | TIMER_NL); 
            OutputBoth (buf);
//		if (res2 /*&& (resaprcl2 != 1) && (resaprcl2 != 2)*/ && (nbdg2 < primolimit))
//			MakePrimoInput (NP, strp);
//	}

	*res = (res1 || res2);

/* Cleanup and return */


//	gwypfree(ps.gx);
//	gwypfree(ps.gy);

EXIT:

	gwypfree(ps.gx);
	gwypfree(ps.gy);
	gwypfree(tmp);
	gwypfree(tmp2);
	gwypfree(apow4);
	gwypfree(gk);
	gwypfree(N);
	gwypfree(NP);
	gwypfree(M);
//        gwypdone();
	lasterr_point = 0;
        

	return (retval);

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
	if (fileExists (filename)) {			// Resuming a Vrba-Reix test
		dovrbareix = TRUE;
	}

restart:

	if (dovrbareix) {				// Compute the seed for the Vrba-Reix test
		gx0 =  newgiant ((bits >> 3) + 8);	// Allocate memory for gx0
		gtog (NP, gx0);				// gx0 = NP
		iaddg (3, gx0);				// gx0 = N+3
		gshiftright (1, gx0);			// gx0 = (N+3)/2 = 3/2 mod N = 3*2^(-1) mod N
		expx = n-1;
		tempFileName (filename, 'z', NP);	// Set the filename to zxxxxxxx
	}
	else {									// Set he base for the SPRP test
		a = IniGetInt (INI_FILE, (char*)"FBase", 3);
		gwypsetmaxmulbyconst (abs(a));
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
                    else  if(bit != (lasterr_point-1)&&(bit > 30) && (bit < expx-30)) 
                        cuda_gwypsquare (x,(saving || stopping || (interimResidues && (bit+1) % interimResidues < 2)) ? 2:0);
                    else
                        cuda_gwypsquare (x,2);
                    care = FALSE;
                }
		else {
                    gwypsquare_carefully (x);
                    care = TRUE;
		}
		if (!dovrbareix && bit == (expx - 1)) {
                    gwypcopy (x, y);
                    uby = ubx;
		}


		CHECK_IF_ANY_ERROR (x, (bit), expx, 6);

/* That iteration succeeded, bump counters */

		if (will_try_larger_fft && (bit == lasterr_point))
                    saving = 1;	// Be sure to restart after this recovery iteration!
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
//			invg (M,tmp2);
			gwypinvg (M,tmp2);
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
//	invg (M,tmp2);
	gwypinvg (M,tmp2);
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
//		invg (NP, tmp2);		// a^(-2) modulo NP
		gwypinvg (NP, tmp2);		// a^(-2) modulo NP
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
//		invg (M,tmp2);
		gwypinvg (M,tmp2);
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
	gwypdone ();

	if (abonroundoff && MAXERR > maxroundoff) {	// Abort...
            aborted = TRUE;
            sprintf (buf, ERRMSG5, checknumber, sgk);
            OutputBoth (buf);
            *res = FALSE;
            gwypfree(NP);
            gwypfree(M);
            gwypfree(testn);
            if (dovrbareix)
                gwypfree (gx0);
//            gwypdone ();
            _unlink (filename);
            if(IniGetInt(INI_FILE, (char*)"StopOnAbort", 0)) {
                IniWriteInt (INI_FILE, (char*)"PgenLine", IniGetInt(INI_FILE, (char*)"PgenLine", 0) + 1);	// Point on the next line
                return (FALSE);
            }
            else
                return (TRUE);
	}

/* Output a message saying we are restarting */

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) {
//            gwypdone ();
            *res = FALSE;
            return (FALSE);
	}

/* Restart */

	if (will_try_larger_fft) {
            IniWriteInt(INI_FILE, (char*)"FFT_Increment", nbfftinc =  IniGetInt(INI_FILE, (char*)"FFT_Increment", 0) + 1);
            if (nbfftinc == maxfftinc)
                abonroundoff = TRUE;	// Don't accept any more Roundoff error.
            _unlink (filename);
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
    unsigned long ninput = n, base, binput, b_2up = 1, b_else = 1, superPRP = 1;
    long mult;
// Lei end
    n_orig = n;

    gformat = format; // save format in a global.

    if(mult = IniGetInt(INI_FILE, (char*)"StopOnPrimedK", 0)) {
        sprintf (outbuf, "ks%s", sgk);
        if(IniGetInt(INI_FILE, outbuf, 0) >= mult) {
                    // is the count for this k value reached ?
            *res = FALSE;   // then, skip this test
//            return TRUE;
            retval = TRUE;
            goto EXITPR;
        }
    }
    else if(mult = IniGetInt(INI_FILE, (char*)"StopOnPrimedN", 0)) {
        sprintf (outbuf, "ns%lu", n);
        if(IniGetInt(INI_FILE, outbuf, 0) >= mult) {
                    // is the count for this n value reached ?
            *res = FALSE;   // then, skip this test
//            return TRUE;
            retval = TRUE;
            goto EXITPR;
        }
    }
    else if(mult = IniGetInt(INI_FILE, (char*)"StopOnPrimedB", 0)) {
        sprintf (outbuf, "bs%s", sgb);
        if(IniGetInt(INI_FILE, outbuf, 0) >= mult) {
                    // is the count for this base value reached ?
            *res = FALSE;   // then, skip this test
//            return TRUE;
            retval = TRUE;
            goto EXITPR;
        }
    }

    if (format == ABCGM) {  // Do the primality test of a Gaussian Mersenne norm
//        return (isGMNP (sgk, n, res));
        retval = isGMNP (sgk, n, res);
        goto EXITPR;

    }

    if (format == ABCSP) {  // Do the PRP test of a Wagstaff number
//        return (isWSPRP (sgk, n, res));
        retval = isWSPRP (sgk, n, res);
        goto EXITPR;
    }

    gb = newgiant (strlen(sgb)/2 + 8);
    nbllr_mallocs++;
                // Allocate one byte per decimal digit + spares
    ctog (sgb, gb); // Convert b string to giant
    if (gb->sign <= 2) {// Test if the base is a small integer...
        base = gb->n[0];// Then, get the base in an unsigned long
        if (gb->sign == 2)
            base += 65536*gb->n[1];
        binput = base;
        while (!(base&1) && base > 2) {
            // Divide the base by two as far as possible
            base >>= 1;
            gshiftright (1, gb); // update also gb JP 21/10/22
            n += ninput;
        }
            

        if (base != 2) {    // Test if the base was a power of two
// Lei
            n -= ninput;
            b_else = base;  // Is odd...
            b_2up = binput / b_else;// binput = b_else*b_2up
            if ((b_2up > b_else) && (!((format == ABCC) || (format == ABCK) || (format == ABCRU) || (format == ABCGRU) || (format == ABCVARAQS)))) {
                superPRP = 0;   // Then b_2up^n > b_else^n
            }
            else {
// Lei end
                base = binput;
                itog (binput, gb); // Restore also gb JP 21/10/22
                    // Do not modify because PRP will be forced...
                n = ninput;
            }
        }
        globalb = base;
            // Keep the base of the candidate in a global

// Lei mod
        if (format == ABCDP) {
            retval = IsPRP (format, sgk, base, n, incr, shift, res);
            free (gb);
            nbllr_frees++;
//            return (retval);
            goto EXITPR;
        }
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
        else if (!IniGetInt (INI_FILE, (char*)"ForcePRP", 0) && (incr ==   +1 || incr == -1) && (format != ABCVARAQS) && 
            (format != ABCRU) && (format != ABCGRU)) {
            retval = plusminustest (sgk, base, n, incr, shift, res);
        }
        else  {
            retval = IsPRP (format, sgk, base, n, incr, shift, res);
        }
    }   // End gb is a small integer.
    else if ((format == NPGCC1 || format == NPGCC2) && !IniGetInt (INI_FILE, (char*)"ForcePRP", 0)) {
        retval = gIsCCP (format, sgk, sgb, gb, n, incr, shift, res);
    }
    else if (!IniGetInt (INI_FILE, (char*)"ForcePRP", 0) && (incr == +1 || incr == -1) && (format != ABCVARAQS) && 
	(format != ABCRU) && (format != ABCGRU)) {
        retval = gplusminustest (sgk, sgb, gb, n, incr, shift, res);
    }
    else  {
//		fermat_only = TRUE;     // JP 30/01/17
        retval = gIsPRP (format, sgk, sgb, gb, n, incr, shift, res);
//		fermat_only = FALSE;    // JP 30/01/17
    }
    free (gb);
    nbllr_frees++;
EXITPR :
    if (nbllr_mallocs!=nbllr_frees)
        printf ("Number of LLR mallocs = %d although Number of LLR frees = %d\n", nbllr_mallocs, nbllr_frees);
    nbllr_mallocs = nbllr_frees = 0;
    return (retval);
}

char	outpf[] = "gqplus.res", outmf[] = "gqminus.res";
char	gqpstring[] = "ABC (2^$a+2^(($a+1)/2)+1)/5\n";
char	gqmstring[] = "ABC (2^$a-2^(($a+1)/2)+1)/5\n";


int primeContinue ()
{

    int	work, nargs, hiline, completed = FALSE;
    unsigned long /*format, */shift, begline, rising_ns, rising_ks, last_processed_n;
    char *pinput;

/* Set appropriate priority */

    SetPriority ();

/* Case off the work type */

    work = IniGetInt (INI_FILE, (char*)"Work", 0);

/* Handle a sieving program output file */

    if (work == 0) {
        char    inputfile[80], outputfile[80], oldinputfile[80], cmaxroundoff[10], cpcfftlim[10], sgk[sgkbufsize], buff[sgkbufsize+256];
        char	hbuff[sgkbufsize+256], outbuf[sgkbufsize+256], last_processed_k[sgkbufsize+256];
        FILE *fd;
        unsigned long i, chainlen, m, n, base, nfudge, nn, k, b, d;
        int firstline, line, hline, resultline,
        outfd, outfdp, outfdm, res, incr, sign, argcnt, validheader = FALSE;
        char c;
#ifdef	WIN32
        giant initgiants = newgiant (1<<19);
            // Create a giant of maximal size
#else
        giant initgiants = newgiant (-1);
            // Create a giant of maximal size
#endif
        gwypfree (initgiants);
            // And free it, to initialize the popg / pushg routines

        IniGetString (INI_FILE, (char*)"PgenInputFile", inputfile,    IBSIZE, NULL);
        IniGetString (INI_FILE, (char*)"OldInputFile", oldinputfile, IBSIZE, inputfile);// Default it to PgenInputFile! JP 26/02/17
        IniGetString (INI_FILE, (char*)"PgenOutputFile", outputfile, IBSIZE, NULL);
        IniGetString (INI_FILE, (char*)"MaxRoundOff", cmaxroundoff, 5, (char*)"0.40");
        maxroundoff = atof (cmaxroundoff);
        IniGetString (INI_FILE, (char*)"MAXBPD", cMAXBPD, 5, (char*)"37.0");
        MAXBPD = atof (cMAXBPD);
        IniGetString (INI_FILE, (char*)"PercentFFTLimit", cpcfftlim, 5, (char*)"0.50");
        pcfftlim = atof (cpcfftlim);
        if (!strcmp (inputfile, oldinputfile))
            firstline = IniGetInt (INI_FILE, (char*)"PgenLine", 1);		// Continuing on the same file
        else
            firstline = 1;		
            // Processing a new file
        last_processed_n = (unsigned long)IniGetInt(INI_FILE, (char*)"Last_Processed_n", 0);
        IniGetString(INI_FILE, (char*)"Last_Processed_k",last_processed_k, sgkbufsize, NULL);
        bpsw = IniGetInt(INI_FILE, (char*)"BPSW", 0);
        hline = IniGetInt (INI_FILE, (char*)"HeaderLine", 0);
        verbose = IniGetInt (INI_FILE, (char*)"Verbose", 0);

// Transmit the pointers to user output fuctions to the gwypnum system.

        gwypsetoutputs (OutputStr, OutputBoth);

        setuponly = IniGetInt (INI_FILE, (char*)"SetupOnly", 0);
        begline = IniGetInt(INI_FILE, (char*)"BegLine", 0);
        fermat_only  = IniGetInt(INI_FILE, (char*)"Fermat_only", 0);
	strong = IniGetInt (INI_FILE, (char*)"StrongFermat", 0);
        testgm  = IniGetInt(INI_FILE, (char*)"TestGM", 1);
        testgq  = IniGetInt(INI_FILE, (char*)"TestGQ", 0);
        testfac  = IniGetInt(INI_FILE, (char*)"TestFac", 0);
        facfrom =  IniGetInt(INI_FILE, (char*)"FacFrom", 0);
        facto =  IniGetInt(INI_FILE, (char*)"FacTo", 0);
        general =  IniGetInt(INI_FILE, (char*)"Vgeneral", 0);
        eps2 =  IniGetInt(INI_FILE, (char*)"Veps2", 0);
        vdebug =  IniGetInt(INI_FILE, (char*)"Vdebug", 0);
        debug =  IniGetInt(INI_FILE, (char*)"Debug", 0);
        tdebug =  IniGetInt(INI_FILE, (char*)"TimeDebug", 0);
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

        if (!strncmp (buff, (char*)"TestWieferichcode", 17)) {
                // Very particular test code...
            TestWieferich ();
            IniWriteInt (INI_FILE, (char*)"Workdone", 1);
            return (FALSE);
        }
        sprintf (SVINI_FILE, "save_%s", INI_FILE);
                    // set the name of the backup Ini File

// Process each line in the input file

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
            
            w->known_factors = NULL;
            quotient = FALSE;
            
// Read the line, break at EOF

            if (fgets (buff, sgkbufsize+256, fd) == NULL) {
                IniWriteInt (INI_FILE, (char*)"Workdone", 1);
                rising_ns = rising_ks = FALSE;
                break;
            }
            else {
                IniWriteInt (INI_FILE, (char*)"Workdone", 0);
            }

// Skip this line if requested (we processed it on an earlier run)
// (but don't ignore last header line found!)

            if (hiline && line > hiline) {
                IniWriteInt (INI_FILE, (char*)"Workdone", 1);
                break;
            }

            if (!strncmp (buff, "ABC", 3)) {
                    // ABC format header found
                sprintf (hbuff, "%s", buff);	// Save the header
                IniWriteInt (INI_FILE, (char*)"HeaderLine", line);
                    // Save the header line number
                hline = line;
                validheader = TRUE; // Assume it is valid...
                for (pinput=buff+3; *pinput && isspace(*pinput); pinput++);
                for (i=0;i<strlen(pinput);i++)
                    if (isspace(pinput[i]))
                        pinput[i] = '\0';
                            // Suppress the EOL characters if necessary.

                if (!strcmp (pinput, grepustring)) {
                    format = ABCGRU;
                }
                else if (!strcmp (pinput, abcadstring)) {
                    format = ABCVARAQS;
                    nargs = 5;
                }
                else if (nargs = sscanf (pinput, abcndstring, &k, &b, &incr, &d)==4) {
                    format = ABCVARAQS;
                    sprintf (sgk, "%lu", k); // convert to string
                    sprintf (sgb, "%lu", b); // convert to string
                    sprintf (sgd, "%lu", d); // convert to string
                    nargs = 1;
                }
                else if (nargs = sscanf(pinput, abcnstring, &k, &b, &incr) == 3) {
                    format = ABCVARAS;
                    sprintf(sgk, "%lu", k); // convert to string
                    sprintf(sgb, "%lu", b); // convert to string
                    nargs = 1;
                }
                else if (!strcmp (pinput, diffnumstring)) {
                    format = ABCDNG;
                }
                else if (!strcmp (pinput, ckstring)) {
                    format = ABCK;
                }
                else if (!strcmp (pinput, dpstring)) {
                    format = ABCDP;
                }
                else if (!strcmp (pinput, repustring)) {
                    format = ABCRU;
                }
                else if (!strcmp (pinput, cwstring)) {
                    format = ABCCW;
                }
                else if (!strcmp (pinput, abcastring)) {
                    format = ABCVARAS;
                    nargs = 4;
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
                else if (sscanf(pinput, fkpstring, &smallk, &incr) == 2){
                    sprintf(sgk, $LLF, smallk); // unsigned fixed k...	
                    format = ABCFKGS;
                }
                else if (sscanf(pinput, fkmstring, &smallk, &incr) == 2){
                    sprintf(sgk, $LLF, smallk); // unsigned fixed k...	
                    format = ABCFKGS;
                    incr = - incr;
                }
                else if (sscanf(pinput, fkpstring, &smallk) == 1) { 
                    sprintf(sgk, $LLF, smallk); // unsigned fixed k...	
					format = ABCFKAS;
                }
                else if (sscanf(pinput, fbpstring, &smallbase, &incr) == 2){
                    sprintf (sgb, $LLF, smallbase);// unsigned fixed b...	
					format = ABCFBGS;
                }
                else if (sscanf(pinput, fbmstring, &smallbase, &incr) == 2){
                    sprintf (sgb, $LLF, smallbase);// unsigned fixed b...	
                    format = ABCFBGS;
                    incr = - incr;
                }
                else if (sscanf(pinput, fbastring, &smallbase) == 1) { 
                    sprintf (sgb, $LLF, smallbase);// unsigned fixed b...	
                    format = ABCFBAS;
                }
                else {
                    OutputBoth ((char*)"Invalid ABC format, next data lines will be flushed...\n");
                    validheader = FALSE;    // Invalid header found...
                }

                if (format == ABCGM) {
                    if (!facto)
                        sprintf (pinput+strlen (gmstring),
                        " // Let GM(p) = (1+/-i)^p-1, GQ(p) = ((1+/-i)^p+1)/(2+/-i) if p>3, (1+/-i)^p+1 if p<=3\n");
                    if (!facto && !fileExists (outpf)) {
                        outfdp = _open (outpf, _O_TEXT | _O_RDWR | _O_CREAT, 0666);
                        if (outfdp) {
                            wc = _write (outfdp, gqpstring, strlen (gqpstring));
                            _close (outfdp);
                        }	
                    }
                    if (!facto && !fileExists (outmf)) {
                        outfdm = _open (outmf, _O_TEXT | _O_RDWR | _O_CREAT, 0666);
                        if (outfdm) {
                            wc = _write (outfdm, gqmstring, strlen (gqmstring));
                            _close (outfdm);
                        }
                    }
                }
                continue;    // Read next line, but do not change PgenLine!
            }		     // End ABC format header found
            else if (((argcnt = sscanf (buff, $LLF":%c:%lu:"$LLF":%lu\n", &li, &c, &chainlen, &smallbase, &mask)) > 1) || !line) {
                if (argcnt < 4) {
                    OutputBoth ((char*)"Missing or invalid NewPGen header, next data lines will be flushed...\n");
                    validheader = FALSE;    // Invalid NewPGen header...
                }
                else {      // Valid NewPGen header
                    sprintf (sgb, $LLF, smallbase);
                        // Newpgen format admits only unsigned long base...	
                    validheader = TRUE;
                    if (argcnt == 4)
                        mask = 0;
                    strcpy (hbuff, buff);   // Save the header
                    IniWriteInt (INI_FILE, (char*)"HeaderLine", line);
                                            // Save the header line number
                    hline = line;
                    format = NPG;
                    if (mask & 0x40) {
                        OutputStr ((char*)"Primorial NewPgen files are not supported...\n");
                        validheader = FALSE;
                    }
                    if (chainlen == 0) chainlen = 1;
                }
                continue;   // Read next line, but do not change PgenLine!
            }		    // End NewPGen header found

            else {	    // Processing a data line
                if (((!rising_ns && !rising_ks) || (rising_ns && rising_ks)) && (line < firstline))
                    continue;
        // Skip this line if requested (we processed it on an earlier run)
                if (!validheader)
                    continue;
        // Flush data until a valid header is found...
                shift = 0;  // Only one value for the k multiplier

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
                        if (c == 'P') /*c = '2',*/ chainlen = 1; // Reentrance ?? 29/04/21
                        if (c == 'M') /*c = '1',*/ chainlen = 1;
//		        if (c == 'Y') /*c = '2',*/ chainlen = 1;
//			if (c == 'Z') /*c = '1',*/ chainlen = 1;
                        if (c == 'T') /*c = '3',*/ chainlen = 1;
                        if (c == 'S') /*c = '1',*/ chainlen = 2;
                        if (c == 'C') /*c = '2',*/ chainlen = 2;
                        if (c == 'B') /*c = '3',*/ chainlen = 2;


// Process each line in the newpgen output file
// allow k to be a big integer
                        
                        if (sscanf (buff+begline, "%s %lu", sgk, &n) != 2)
                            continue;	// Skip invalid line
                        if (!isDigitString(sgk))
                            continue;	// Skip invalid line
                        if (rising_ns && !rising_ks && (n <= last_processed_n))
                            continue;	// Skip already processed n's
                        if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
                            continue;	// Skip already processed k's

                        if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                            fclose (fd);
                            // Unlock the file during the test...
                                                        
// Test numbers according to the c variable

                        nn = n;
                        if (c == 'Y') {
                            nn--;
                        }
                        if (c == 'Z') {
                            nn--;
                        }
                        for (i = 0; i < chainlen; i++) {
                            if (c == '1' || c == '3' || c == 'M' || c == 'S' || c == 'T' || c == 'B') {
                                if (! process_num (format, sgk, sgb, n - nfudge + i, -1, shift, &res))
                                    goto done;
                                if (!res)
                                    break;
                                if (c == '1' || c == 'M' || c == 'S')
                                    format = NPGCC1;
                            }
                            if (c == '2' || c == '3' || c == 'P' || c == 'C' || c == 'T' || c == 'B') {
                                if (! process_num (format, sgk, sgb, n - nfudge + i, +1, shift, &res))
                                    goto done;
                                if (!res)
                                    break;
                                if (c == '2' || c == 'P' || c == 'C')
                                    format = NPGCC2;
                            }
                            if (c == 'J') {	// Twin/SG
                                int res2;
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
                                int res2;
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
                                int res2;
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
                            if (c == 'A') {// Arithmetic Progression mode 
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
                            continue;	// Skip already processed n's
                        if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
                            continue;	// Skip already processed k's
                        if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                            fclose (fd);
                                    // Unlock the file during the test...

// Undo the increment of n that newpgen did on types 1, 2, 3

                        nn = n;
//			if (c == '1' || c == '2' || c == '3')
//			   nn--;
                        if (c == 'S' && !(mask & MODE_PLUS))    // 30/04/21
                            chainlen = 2;
                        if (c == 'C' && !(mask & MODE_PLUS))   // 01/05/21
                            chainlen = 2;
                        if (c == 'B')
                            chainlen = 2;
                        if ((mask & MODE_PLUS) && (mask & MODE_2MINUS) && (mask & MODE_2PLUS) && (mask & MODE_4PLUS))
                            nn--;
                        if ((mask & MODE_MINUS) && (mask & MODE_2MINUS) && (mask & MODE_2PLUS) && (mask & MODE_4MINUS))
                            nn--;

// Test numbers according to the mask variable
// The J and K types (Twin/CC and Twin/SG) are special in that they
// are output if either a Twin OR a CC/SG is found

                        shift = 0;
                        for (i = 0; i < chainlen; i++) {
                            if ((mask & MODE_MINUS) && (mask & MODE_PLUS) && (mask & MODE_2MINUS)) {	// Twin/SG
                                int res2;
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
                            if ((mask & MODE_MINUS) && (mask & MODE_PLUS) && (mask & MODE_2PLUS)) {	// Twin/CC
                                int res2;
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
                            if ((mask & MODE_PLUS) && (mask & MODE_2MINUS)&& (mask & MODE_2PLUS) && (mask & MODE_4PLUS)) {
                                // Lucky Plus
                                int res2;
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
                            if ((mask & MODE_MINUS) && (mask & MODE_2MINUS) &&(mask & MODE_2PLUS) && (mask & MODE_4MINUS)) {	       // Lucky Minus
                                int res2;
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
                                else {
                                    if (! process_num (format, sgk, sgb, nn, -1, shift, &res))
                                        goto done;
                                }
                                if (!res)
                                    break;
                            }
                            if (mask & MODE_PLUS) {
                                if (mask & MODE_DUAL) {
                                    if (! process_num (format, (char*)"1", sgb, nn, atoi(sgk), shift, &res))
                                        goto done;
                                }
                                else {
                                    if (! process_num (format, sgk, sgb, nn, +1, shift, &res))
                                        goto done;
                                }
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
                            if (hline >= resultline) {
                                // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%s %lu\n", sgk, n);
                                // write the result
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	        // update the result line
                    }
                }	// End of NewPGen format processing

                else if (format == ABCCW) {	// Cullen/Woodall
                    if (sscanf (buff+begline, "%lu %s %d", &n, sgb, &incr) != 3)
                        continue;	// Skip invalid line
                    if (!isDigitString (sgb))
                        continue;	// Skip invalid line
                    if (rising_ns && !rising_ks  && (n <= last_processed_n))
                        continue;	// Skip already processed n's
                    sprintf (sgk, "%lu", n);
                    if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
                        continue;	// Skip already processed k's
                    if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);// Unlock the file during the test...
                    sprintf (sgk, "%lu", n);
                    if (! process_num (format, sgk, sgb, n, incr, shift, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {
                                // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%lu %lu %d\n", n, base, incr); 
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
                    }
                }
                else if (format == ABCFF)   {   // FermFact output                                                                                                        
                        // allow k to be a big integer
                    if (sscanf (buff+begline, "%s %lu", sgk, &n) != 2)
                        continue;	// Skip invalid line
                    if (!isDigitString(sgk))
                        continue;	// Skip invalid line
                    if (rising_ns && !rising_ks && (n <= last_processed_n))
                        continue;	// Skip already processed n's
                    if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
                        continue;	// Skip already processed k's
                    sprintf (sgb, "2");
                    if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);// Unlock the file during the test...
                    if (! process_num (format, sgk, (char*)"2", n, +1, shift, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {
                                    // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%s %lu\n", sgk, n); 
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
                    }
                }
                else if (format == ABCLEI)  {   // Lei output
                            // allow k to be a big integer
                    if (sscanf (buff+begline, "%s %lu", sgk, &n) != 2)
                        continue;	// Skip invalid line
                    if (!isDigitString(sgk))
                        continue;	// Skip invalid line
                    if (rising_ns && !rising_ks && (n <= last_processed_n))
                        continue;	// Skip already processed n's
                    if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
                        continue;	// Skip already processed k's
                    sprintf (sgb, "2");
                    if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);// Unlock the file during the test...
                    if (! process_num (format, sgk, (char*)"2", n, -1, shift, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {
                                    // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%s %lu\n", sgk, n);
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);
                                    // update the result line
                    }
                }
                else if (format == ABCFKGS)	{
                    // Fixed k:  b and n specified on each input line
                    if (sscanf (buff+begline, "%s %lu", sgb, &n) != 2)
                        continue;	// Skip invalid line
                    if (!isDigitString (sgb))
                        continue;	// Skip invalid line
                    if (rising_ns && (n <= last_processed_n))
                        continue;	// Skip already processed n's
                    if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);// Unlock the file during the test...
                    if (! process_num (format, sgk, sgb, n, incr, shift, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {
                                // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%s %lu\n", sgb, n); 
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
                    }
                }
                else if (format == ABCFKAS)	{
                    // Fixed k: b, n, and c specified on each input line
                    if (sscanf (buff+begline, "%s %lu %d", sgb, &n, &incr) != 3)
                        continue;   // Skip invalid line
                    if (!isDigitString(sgk))
                        continue;	// Skip invalid line
                    if (!isDigitString (sgb))
                        continue;	// Skip invalid line
                    if (rising_ns && (n <= last_processed_n))
                        continue;	// Skip already processed n's
                    if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);
                            // Unlock the file during the test...
                    if (! process_num (format, sgk, sgb, n, incr, shift, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {
                                // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%s %lu %d\n", sgb, n, incr); 
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
                    }
                }
                else if (format == ABCFBGS)	{
                    // Fixed b:  k and n specified on each input line
                    if (sscanf (buff+begline, "%s %lu", sgk, &n) != 2)
                        continue;	// Skip invalid line
                    if (!isDigitString(sgk))
                        continue;	// Skip invalid line
                    if (rising_ns && !rising_ks && (n <= last_processed_n))
                        continue;	// Skip already processed n's
                    if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
                        continue;	// Skip already processed k's
                    if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);// Unlock the file during the test...
                    if (! process_num (format, sgk, sgb, n, incr, shift, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {
                                // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%s %lu\n", sgk, n); 
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
                    }
                }
                else if (format == ABCFBAS)	{
                    // Fixed b: k, n, and c specified on each input line
                    if (sscanf (buff+begline, "%s %lu %d", sgk, &n, &incr) != 3)
                        continue;	// Skip invalid line
                    if (!isDigitString(sgk))
                        continue;	// Skip invalid line
                    if (rising_ns && !rising_ks && (n <= last_processed_n))
                        continue;	// Skip already processed n's
                    if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
                        continue;	// Skip already processed k's
                    if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);
                            // Unlock the file during the test...
                    if (! process_num (format, sgk, sgb, n, incr, shift, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {
                                // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%s %lu %d\n", sgk, n, incr); 
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	         // update the result line
                    }
                }
                else if (format == ABCFNGS)	{
                    // Fixed n:  k and b specified on each input line
                    if (sscanf (buff+begline, "%s %s", sgk, sgb) != 2)
                        continue;	// Skip invalid line
                    if (!isDigitString(sgk))
                        continue;	// Skip invalid line
                    if (!isDigitString (sgb))
                        continue;	// Skip invalid line
                    if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
                        continue;	// Skip already processed k's
                    if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);// Unlock the file during the test...
                    if (! process_num (format, sgk, sgb, n, incr, shift, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {
                                // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%s %s\n", sgk, sgb); 
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	     // update the result line
                    }
                }
                else if (format == ABCFNAS)	{
                    // Fixed n:  k, b, and c specified on each input line
                    if (sscanf (buff+begline, "%s %s %d", sgk, sgb, &incr) != 3)
                        continue;	// Skip invalid line
                    if (!isDigitString(sgk))
                        continue;	// Skip invalid line
                    if (!isDigitString (sgb))
                        continue;	// Skip invalid line
                    if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
                        continue;	// Skip already processed k's
                    if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);// Unlock the file during the test...
                    if (! process_num (format, sgk, sgb, n, incr, shift, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {
                                // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%s %s %d\n", sgk, sgb, incr); 
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	     // update the result line
                    }
                }
                else if (format == ABCVARGS)    {
                    // k, b, and n specified on each input line
                    if (sscanf (buff+begline, "%s %s %lu", sgk, sgb, &n) != 3)
                        continue;	// Skip invalid line
                    if (!isDigitString(sgk))
                        continue;	// Skip invalid line
                    if (!isDigitString (sgb))
                        continue;	// Skip invalid line
                    if (rising_ns && !rising_ks && (n <= last_processed_n))
                        continue;	// Skip already processed n's
                    if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
                        continue;	// Skip already processed k's
                    if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);// Unlock the file during the test...
                    if (! process_num (format, sgk, sgb, n, incr, shift, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {
                                // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%s %s %lu\n", sgk, sgb, n); 
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	     // update the result line
                    }
                }
                else if ((format == ABCVARAS) && (nargs == 4))	{// k, b, n, and c specified on each input line
                    if (sscanf(buff + begline, "%s %s %lu %d", sgk, sgb, &n, &incr)!=4)
                        continue;   // Skip invalid line
                    if (!isDigitString(sgk))
                        continue;   // Skip invalid line
                    if (!isDigitString (sgb))
                        continue;   // Skip invalid line
                    if (rising_ns && !rising_ks && (n <= last_processed_n))
                        continue;   // Skip already processed n's
                    if (rising_ks && !rising_ns && (digitstrcmp (sgk,
                        last_processed_k) <= 0))
                        continue;   // Skip already processed k's
                    if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);// Unlock the file during the test...
                    if (! process_num (format, sgk, sgb, n, incr, shift, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {	// write the relevant header
                                writelg = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%s %s %lu %d\n", sgk, sgb, n, incr); 
                            writelg = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);// update the result line
                    }
                }
                else if ((format == ABCVARAS) && (nargs == 1)) {// Only n specified on each input line
                    if (sscanf(buff + begline, "%lu", &n)!=1)
                        continue;   // Skip invalid line
                    if (!isDigitString(sgk))
                        continue;   // Skip invalid line
                    if (!isDigitString(sgb))
                        continue;   // Skip invalid line
                    if (rising_ns && (n <= last_processed_n))
                        continue;   // Skip already processed n's
                    if (rising_ns)
                        fclose(fd); // Unlock the file during the test...
                    if (!process_num(format, sgk, sgb, n, incr, shift, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open(outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {// write the relevant header
                                writelg = _write(outfd, hbuff, strlen(hbuff));
                            }
                            sprintf(outbuf, "%lu\n", n);
                            writelg = _write(outfd, outbuf, strlen(outbuf));
                            _close(outfd);
                        }
                        IniWriteInt(INI_FILE, (char*)"ResultLine", line);// update the result line
                    }
                }
                else if (format == ABCRU)	{
                    // Repunits, n is the only parameter.
                    if (sscanf (buff+begline, "%lu", &n) != 1)
                        continue;	// Skip invalid line
                    sprintf (sgb, "10");
                    sprintf (sgk, "1");
                    if (rising_ns && !rising_ks && (n <= last_processed_n))
                        continue;	// Skip already processed n's
                    if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);// Unlock the file during the test...
                    if (! process_num (format, (char*)"1", (char*)"10", n, -1, 0, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {
                                // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%lu\n", n); 
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	     // update the result line
                    }
                }
                else if (format == ABCDP)   {
                    // DivPhi, k, b, n specified on each input line
                    if (sscanf (buff+begline, "%s %s %lu", sgk, sgb, &n) != 3)
                        continue;   // Skip invalid line
                    if (!isDigitString(sgk))
                        continue;   // Skip invalid line
                    if (!isDigitString (sgb))
                        continue;   // Skip invalid line
                    incr = 1;
                    if (! process_num (format, sgk, sgb, n, 1, 0, &res))
                        goto done;
                }
                else if (format == ABCGRU)	{
                    // Generalized Repunits, b, n, are the two parameters
                    if (sscanf (buff+begline, "%s %lu", sgb, &n) != 2)
                        continue;	// Skip invalid line
                    if (!isDigitString (sgb))
                        continue;	// Skip invalid line
                    sprintf (sgk, "1");
                    if (rising_ns && !rising_ks && (n <= last_processed_n))
                        continue;	// Skip already processed n's
                        if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);// Unlock the file during the test...
                    if (! process_num (format, (char*)"1", sgb, n, -1, 0, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {
                                // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%s %lu\n", sgb, n); 
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                            }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	     // update the result line
                    }
                }
                else if (format == ABCGF)	{
                    // Generalized Fermat, sgb, n, are the two parameters
                    if (sscanf (buff+begline, "%s %lu", sgb, &n) != 2)
                        continue;	// Skip invalid line
                    if (!isDigitString(sgb))
                        continue;   // Skip invalid line
                    if (!ispoweroftwo(n))
                        continue;	// Skip invalid line
                    sprintf (sgk, "1");
                    if (! process_num (format, (char*)"1", sgb, n, 1, 0, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {
                                    // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%s %lu\n", sgb, n);
                                // write the result
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	     // update the result line
                    }
                }
                else if (format == ABCDN)	{
                    // b^n-b^m+c numbers ; sgb, n, m are the 3 parameters
                    if (sscanf (buff+begline, "%s %lu %lu", sgb, &n, &m) != 3)
                        continue;	// Skip invalid line
                    if (!isDigitString(sgb))
                        continue;	// Skip invalid line
                    if (n <= m)
                        continue;	// Skip invalid line
                    ndiff = n-m;
                        // Save difference of exponents in a global
                    sprintf (sgk, "1");
                    if (! process_num (format, (char *)"1", sgb, m, incr, 0, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char *)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {
                                // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%s %lu %lu\n", sgb, n, m);	   // write the result
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char  *)"ResultLine", line);	       // update the result line
                    }
                }
                else if (format == ABCDNG)	{
                    // b^n-b^m+c numbers ; sgb, n, m, c are the 4 parameters
                    if (sscanf (buff+begline, "%s %lu %lu %d", sgb, &n, &m, &incr) != 4)
                        continue;	// Skip invalid line
                    if (!isDigitString(sgb))
                        continue;	// Skip invalid line
                    if (n <= m)
                        continue;	// Skip invalid line
                    ndiff = n-m;
                        // Save difference of exponents in a global
                    sprintf (sgk, "1");
                    if (! process_num (format, (char *)"1", sgb, m, incr, 0, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char *)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {
                                // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%s %lu %lu %d\n", sgb, n, m, incr);	// write the result
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char *)"ResultLine", line);	        // update the result line
                    }
                }
                else if (format == ABCVARAQS && nargs == 5) {	// k, b, n, c and d specified on each input line
                    if (sscanf (buff+begline, "%s %s %lu %d %s", sgk, sgb, &n, &incr, sgd) != 5)
                        continue;   // Skip invalid line
                    if (!isDigitString(sgk))
                        continue;   // Skip invalid line
                    if (!isDigitString (sgb))
                        continue;   // Skip invalid line
//		    if (!isDigitString(sgd))
//			continue;   // Skip invalid line
                    if (rising_ns && !rising_ks && (n <= last_processed_n))
                        continue;   // Skip already processed n's
                    if (rising_ks && !rising_ns && (digitstrcmp (sgk,
                        last_processed_k) <= 0))
                        continue;   // Skip already processed k's
                    if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);    // Unlock the file during the test...
                    if (! process_num (format, sgk, sgb, n, incr, shift, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {// write the relevant header
                            writelg = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%s %s %lu %d %s\n", sgk, sgb, n, incr, sgd); 
                            writelg = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
                    }
                }
                else if (format == ABCVARAQS && nargs == 1) {// Only n specified on each input line.
                    if (sscanf (buff+begline, "%lu", &n) != 1)
                        continue;   // Skip invalid line
                    if (!isDigitString(sgk))
                        continue;   // Skip invalid line
                    if (!isDigitString (sgb))
                        continue;   // Skip invalid line
//		    if (!isDigitString(sgd))
//			continue;   // Skip invalid line
                    if (rising_ns && (n <= last_processed_n))
                        continue;   // Skip already processed n's
//		    if (rising_ks && !rising_ns && (digitstrcmp (sgk, last_processed_k) <= 0))
//			continue;   // Skip already processed k's
                    if (rising_ns)
                        fclose (fd);// Unlock the file during the test...
                    if (! process_num (format, sgk, sgb, n, incr, shift, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {	// write the relevant header
                                writelg = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%lu\n", n); 
                            writelg = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);// update the result line
                    }
                }
                else if (format == ABCGM)	{// Gaussian Mersenne
                    if ((nargs = sscanf (buff+begline, "%lu %lu %lu", &n, &facn, &facnp)) < 1)
                        continue;	// Skip invalid line
                    else if (nargs == 1)
                        // Not prefactored.
                        facn = facnp = 0;
                    else if (nargs == 2) {
                        // Second argument is how far already factored, in bits)
                        if (!facfrom)
                            facfrom = facn;
                        facn = facnp = 0;
                    }
                    if (rising_ns && !rising_ks && (n <= last_processed_n))
                        continue;	// Skip already processed n's
                    sprintf (sgk, "2^%lu", (n+1)/2);
                    sprintf (sgb, "2");
                    if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);// Unlock the file during the test...
                    if (! process_num (format, sgk, (char*)"2", n, +1, shift, &res))
                        goto done;
#ifndef X86_64
                    if (facto) {// If factoring, print a job progress message every so often
                        if (n/pdivisor-pquotient == 1) {
                            sprintf (outbuf, "%lu candidates factored, %lu factors found, %lu remaining\n", factored, eliminated, factored - eliminated);
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
                            if (hline >= resultline) {
                                // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        if (res2) {
                            if (sign) {
                                outfdm = _open (outmf, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                                if (outfdm) {
                                    sprintf (outbuf, "%lu\n", n); 
                                    wc = _write (outfdm, outbuf, strlen (outbuf));
                                    _close (outfdm);
                                }
                            }
                            else {
                                outfdp = _open (outpf, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                                if (outfdp) {
                                    sprintf (outbuf, "%lu\n", n); 
                                    wc = _write (outfdp, outbuf, strlen (outbuf));
                                    _close (outfdp);
                                }
                            }
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	     // update the result line
                    }
                }
                else if (format == ABCSP)   {
                    // SPRP test of (2^n+1)/3 numbers
                    if ((nargs = sscanf (buff+begline, "%lu %lu", &n, &facn)) < 1)
                        continue;	// Skip invalid line
                    else if (nargs == 1)
                                        // Not prefactored.
                        facn = facnp = 0;
                    else if (nargs == 2) {
                        // Second argument is how far already factored, in bits)
                        if (!facfrom)
                            facfrom = facn;
                        facn = facnp = 0;
                    }
                    if (rising_ns && !rising_ks  && (n <= last_processed_n))
                        continue;	// Skip already processed n's
                    sprintf (sgk, "(2^%lu+1)/3", n);
                    sprintf (sgb, "2");
                    if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);// Unlock the file during the test...
                    if (! process_num (format, sgk, (char*)"2", n, +1, shift, &res))
                        goto done;
#ifndef X86_64
                    if (facto) {// If factoring, print a job progress message every so often
                        if (n/pdivisor-pquotient == 1) {
                            sprintf (outbuf, "%lu candidates factored, %lu factors found, %lu remaining\n", factored, eliminated, factored - eliminated);
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
                            if (hline >= resultline) {
                                // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	     // update the result line
                    }
                }
                else if (format == ABCK) {	// Carol/Kynea
                    if (sscanf (buff+begline, "%lu %d", &n, &incr) != 2)
                        continue;	// Skip invalid line
                    if (rising_ns && !rising_ks  && (n <= last_processed_n))
                        continue;	// Skip already processed n's
                    if (incr == 1) {
                        format = ABCK;  // Kynea number
                        sprintf (sgk, "(2^%lu+1)", n-1);
                    }
                    else if (incr == -1) {
                        format = ABCC;  // Carol number
                        sprintf (sgk, "(2^%lu-1)", n-1);
                    }
                    else
                        continue;
                    sprintf (sgb, "2");
                    if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                        fclose (fd);// Unlock the file during the test...
                    if (! process_num (format, sgk, (char*)"2", n+1, -1, shift, &res))
                        goto done;
                    if (res) {
                        resultline = IniGetInt(INI_FILE, (char*)"ResultLine", 0);
                        outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
                        if (outfd) {
                            if (hline >= resultline) {
                                // write the relevant header
                                wc = _write (outfd, hbuff, strlen (hbuff));
                            }
                            sprintf (outbuf, "%lu %d\n", n, incr); 
                            wc = _write (outfd, outbuf, strlen (outbuf));
                            _close (outfd);
                        }
                        IniWriteInt (INI_FILE, (char*)"ResultLine", line);	// update the result line
                    }
                }
            }	// End processing a data line
            
//            if (nbllr_mallocs!=nbllr_frees)
//                printf ("Number of LLR mallocs = %d although Number of LLR frees = %d\n", nbllr_mallocs, nbllr_frees);

            IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
            recovering = FALSE; // 21/04/21
            
            if ((!rising_ns && !rising_ks) || (rising_ns && rising_ks))
                IniWriteInt (INI_FILE, (char*)"PgenLine", line + 1);		// Point on the next line
            if (rising_ns && !rising_ks) {
                IniWriteInt (INI_FILE, (char*)"Last_Processed_n", n);		// Point on the next n
                last_processed_n = n;
            }
            if (rising_ks && !rising_ns) {
                IniWriteString (INI_FILE, (char*)"Last_Processed_k", sgk);  // Point on the next k
                strcpy (last_processed_k, sgk);
            }
            if(n>=(unsigned long)IniGetInt(INI_FILE, (char*)"MaxN", 2147483647)) {
                break;
            }
            if (res) {
                if(IniGetInt(INI_FILE, (char*)"BeepOnSuccess", 0)) {
                    do {// If stopping on this success, beep infinitely!
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
                else if(IniGetInt(INI_FILE, (char*)"StopOnPrimedK", 0)){
                    sprintf (outbuf, "ks%s", sgk);
                    IniWriteInt (INI_FILE, outbuf, 1+IniGetInt(INI_FILE, outbuf, 0));
                        // Increment this k success count
                    save_IniFile (INI_FILE, SVINI_FILE);
                        // make a backup of INI_FILE
                }
                else if(IniGetInt(INI_FILE, (char*)"StopOnPrimedN", 0)){
                    sprintf (outbuf, "ns%lu", n);
                    IniWriteInt (INI_FILE, outbuf, 1+IniGetInt(INI_FILE, outbuf, 0));
                        // Increment this n success count
                    save_IniFile (INI_FILE, SVINI_FILE);
                        // make a backup of INI_FILE
                }
                else if(IniGetInt(INI_FILE, (char*)"StopOnPrimedB", 0)){
                    sprintf (outbuf, "bs%s", sgb);
                    IniWriteInt (INI_FILE, outbuf, 1+IniGetInt(INI_FILE, outbuf, 0));
                        // Increment this base success count
                    save_IniFile (INI_FILE, SVINI_FILE);
                        // make a backup of INI_FILE
                }
            }
            if ((rising_ns && !rising_ks) || (!rising_ns && rising_ks))
                goto OPENFILE;
        }       // End of loop on input lines
        IniWriteString (INI_FILE, (char*)"ResultLine", NULL);
                // delete the result line
        _unlink (SVINI_FILE);					// delete the backup of INI_FILE
        completed = TRUE;
done:
//        if (nbllr_mallocs!=nbllr_frees)
//            printf ("Number of LLR mallocs = %d although Number of LLR frees = %d\n", nbllr_mallocs, nbllr_frees);
        IniWriteString(INI_FILE, (char*)"FFT_Increment", NULL);
        recovering = FALSE; // 21/04/21
        
        if(IniGetInt(INI_FILE, (char*)"StopOnSuccess", 0) && res) {
            if ((!rising_ns && !rising_ks) || (rising_ns && rising_ks))
                IniWriteInt (INI_FILE, (char*)"PgenLine", line + 1);	// Point on the next line
            if (rising_ns && !rising_ks)
                IniWriteInt (INI_FILE, (char*)"Last_Processed_n", n);	// Point on the next n
            if (rising_ks && !rising_ns)
                IniWriteString (INI_FILE, (char*)"Last_Processed_k", sgk);   // Point on the next k
        }
        else if (!aborted && ((!rising_ns && !rising_ks) || (rising_ns && rising_ks)))
            IniWriteInt (INI_FILE, (char*)"PgenLine", line);	// Point again on the current line...
        IniWriteString (INI_FILE, (char*)"MaxRoundOff", NULL);
        if (facto) {
            sprintf (outbuf, "%lu candidates factored, %lu factors found, %lu remaining\n", factored, eliminated, factored - eliminated);
                OutputBoth (outbuf);
        }
        if ((!rising_ns && !rising_ks) || (rising_ns && rising_ks))
            fclose (fd);
        IniWriteString(INI_FILE, (char*)"OldInputFile", inputfile);		       // Save the just processed input file name.
    }  // End Work == 0

// Handle an expr

    else {	// Work != 0
        OutputStr ((char*)"Expression testing not yet implemented.\n");
        IniWriteInt (INI_FILE, (char*)"Workdone", 1);
    }
    aborted = FALSE;
    return (completed);
}
