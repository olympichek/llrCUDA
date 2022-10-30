/* Copyright 1995-2000 Just For Fun Software, Inc. */
/* Author:  George Woltman */
/* Email: woltman@alum.mit.edu */
/* Adapted for LLR program by Jean Penné */
/* Email : jpenne@free.fr */

/* Include files */

//#include <string.h>
//#include <stdio.h>

#include <ctype.h>
#include <fcntl.h>
#include <math.h>
#include <memory.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#if defined (__linux__) || defined (__FreeBSD__) || defined (__APPLE__)
#include <dirent.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#define $LLF "%qi"
#define __int64 long long
#else
#include <direct.h>
#include <dos.h>
#include <io.h>
#include <time.h>
#include <process.h>
#include <windows.h>
#define $LLF "%I64d"
#endif
#include <sys/timeb.h>

#if !defined (__linux__) && !defined (__FreeBSD__) && !defined (__APPLE__)

#include <io.h>
#include <stdlib.h>

#endif

#include "lprime.h"
//#include "gwypini.c"

/* Routine definitions */

void options_cpu ();

/* Get line from the user (stdin) */

void get_line (
	char	*buf)
{
	char	c;
	int	i;
	for (i = 0; ; i++) {
		if (_read (0, &c, 1) != 1) break;
		if (c == '\n' || c == 0) break;
		if (i < 80) *buf++ = c;
	}
	*buf++ = 0;
}

/* Get a number from the user */

unsigned long get_number (
	unsigned long dflt)
{
	char	line[80];
	get_line (line);
	if (line[0] == 0) return (dflt);
	return (atol (line));
}

/* Ask a Yes/No question */

void askYN (
	char	*str,
	int	*val)
{
	char	buf[80];
	printf ("%s (%s): ", str, *val ? "Y" : "N");
	get_line (buf);
	if (buf[0] == 0) return;
	*val = (buf[0] == 'Y' || buf[0] == 'y');
}

/* Ask a number question */

void askNum (
	char	*str,
	unsigned long *val,
	unsigned long min,
	unsigned long max)
{
	char	buf[80];
	unsigned long newval;
	printf ("%s (%ld): ", str, *val);
loop:	get_line (buf);
	if (buf[0] == 0) return;
	newval = atol (buf);
	if (min || max) {
		if (newval < min || newval > max) {
			printf ("Please enter a value between %ld and %ld. ",
				min, max);
			goto loop;
		}
	}
	*val = newval;
}

/* Ask a number question */

void askNumNoDflt (
	char	*str,
	unsigned long *val,
	unsigned long min,
	unsigned long max)
{
	char	buf[80];
	unsigned long newval;
	printf ("%s: ", str);
loop:	get_line (buf);
	if (buf[0] == 0) goto loop;
	newval = atol (buf);
	if (min || max) {
		if (newval < min || newval > max) {
			printf ("Please enter a value between %ld and %ld. ",
				min, max);
			goto loop;
		}
	}
	*val = newval;
}

/* Ask a string question */

void askStr (
	char	*str,
	char	*val,
	unsigned long maxlen)
{
	char	buf[80];
	if (val[0])
		printf ("%s (%s): ", str, val);
	else
		printf ("%s: ", str);
loop:	get_line (buf);
	if (buf[0] == 0) return;
	if (strlen (buf) > maxlen) {
		printf ("Maximum string length is %ld characters. ", maxlen);
		goto loop;
	}
	strcpy (val, buf);
}

/* Wait for user input - gives the user time to read the screen */

void askOK ()
{
	char	str[80];
	if (THREAD_KILL) return;
	printf ("\nHit enter to continue: ");
	get_line (str);
}

/* Ask user if he is satisfied with his dialog responses */

int askOkCancel ()
{
	char	buf[80];
	if (THREAD_KILL) return (FALSE);
	printf ("\nAccept the answers above? (Y): ");
	get_line (buf);
	return (buf[0] == 0 || buf[0] == 'Y' || buf[0] == 'y');
}

/* Ask user if he is satisfied with his dialog responses */

int askYesNo (
	char	dflt)
{
	char	buf[80];
	if (THREAD_KILL) return (FALSE);
	printf (" (%c): ", dflt);
	get_line (buf);
	if (buf[0] == 0) buf[0] = dflt;
	return (buf[0] == 'Y' || buf[0] == 'y');
}

/* Output a long string with a max of 75 characters to a line */

void outputLongLine (
	char	*buf)
{
	char	line[80];
	char	*p;
	int	i, j;

	for (p = buf; ; ) {
		for (i = 0; i < 75; i++) {
			line[i] = p[i];
			if (p[i] == 0 || p[i] == '\n') { j = i; break; }
			if (p[i] == ' ' || p[i] == '.' || p[i] == ',') j = i;
		}
		line[j+1] = 0;
		printf ("%s", line);
		if (p[j] == 0) break;
		if (p[j] != '\n') printf ("\n");
		p += j + 1;
		while (*p == ' ') p++;
	}
}

/* Test/InputData dialog */

void test_inputdata ()
{
	char	m_pgen_input[80], m_pgen_output[80];
	unsigned long m_pgen_line;

	IniGetString (INI_FILE, (char*)"PgenInputFile", m_pgen_input, 80, NULL);
	IniGetString (INI_FILE, (char*)"PgenOutputFile", m_pgen_output, 80, NULL);
	m_pgen_line = IniGetInt (INI_FILE, (char*)"PgenLine", 1);

	askStr ((char*)"Input file (from NewPgen): ", m_pgen_input, 76);
	askStr ((char*)"Output file (Results): ", m_pgen_output, 76);
	askNum ((char*)"Line number", &m_pgen_line, 1, 999999999);

	if (askOkCancel ()) {
		IniWriteInt (INI_FILE, (char*)"Work", 0);
		IniWriteString (INI_FILE, (char*)"PgenInputFile", m_pgen_input);
		IniWriteString (INI_FILE, (char*)"PgenOutputFile", m_pgen_output);
		IniWriteInt (INI_FILE, (char*)"PgenLine", m_pgen_line);
		IniWriteInt (INI_FILE, (char*)"WorkDone", 0);
	}
}

/* Advanced/Priority dialog */

void advanced_priority ()
{
	unsigned long m_priority;

	m_priority = PRIORITY;

	outputLongLine ((char*)"Pick a priority between 1 and 10 where 1 is the lowest priority and 10 is the highest.\n");
	outputLongLine ((char*)"It is strongly recommended that you use the default priority of 1.  Your throughput will probably not improve by using a higher priority.  The only time you should raise the priority is when another process, such as a screen saver, is stealing CPU cycles from this program.\n");
	askNum ((char*)"Priority", &m_priority, 1, 10);

	if (askOkCancel ()) {
		PRIORITY = m_priority;
		IniWriteInt (INI_FILE, (char*)"Priority", PRIORITY);
	}
}

/* Options/CPU dialog */

void options_cpu () 
{ 
//	char buf[512]; 
 
//	getCpuDescription (buf, 0); 
//	printf ("CPU Information:\n%s\n", buf); 
	printf ("No CPU Information for now...\n"); 
 
	askOK (); 
} 

/* Options/Preferences dialog */

void options_preferences ()
{
	unsigned long m_iter, m_r_iter, m_disk_write_time;
	int	m_backup;

	m_iter = ITER_OUTPUT;
	m_r_iter = ITER_OUTPUT_RES;
	m_disk_write_time = DISK_WRITE_TIME;
	m_backup = TWO_BACKUP_FILES;

	askNum ((char*)"Iterations between screen outputs", &m_iter, 1, 999999999);
	askNum ((char*)"Iterations between results file outputs",
		&m_r_iter, 10000, 999999999);
	askNum ((char*)"Minutes between disk writes", &m_disk_write_time, 10, 999999);
	askYN ((char*)"Create Two Backup Files", &m_backup);

	if (askOkCancel ()) {
		ITER_OUTPUT = m_iter;
		ITER_OUTPUT_RES = m_r_iter;
		DISK_WRITE_TIME = m_disk_write_time;
		TWO_BACKUP_FILES = m_backup;
		IniWriteInt (INI_FILE, (char*)"OutputIterations", ITER_OUTPUT);
		IniWriteInt (INI_FILE, (char*)"ResultsFileIterations", ITER_OUTPUT_RES);
		IniWriteInt (INI_FILE, (char*)"DiskWriteTime", DISK_WRITE_TIME);
		IniWriteInt (INI_FILE, (char*)"TwoBackupFiles", TWO_BACKUP_FILES);
	}
}

/* Help/About */

void help_about ()
{
	printf ("Primality Testing of k*b^n+/-1 Program - GPU Version 3.8.7\n");
	printf ("Using George Woltman's 1/k method IBDWT for k's up to 19 bits\n"); 
	printf (" << Rational Base DWT for larger k's up to maximal double\n");
	printf (" << Rational base DWT and generic reduction for larger k's or base != 2\n");
        printf ("Using GPU code written by Shoichiro Yamada for IBDWT,\n << extended further to Rational Bases and Generic reduction.\n");
        printf ("By default, all tests on base two numbers use Gerbicz error checking.\nThis is the case for PRP Fermat and SPRP tests as in Prime95 or Mprime,\nbut also for the deterministic prime tests of Proth numbers.\nLLR tests on Riesel numbers are only done after a positive Fermat PRP result.\n");
	printf ("Written : Jan. 20011 - Oct. 2022 by Jean Penne\n");
	printf ("Email : jpenne@free.fr\n");
	askOK ();
}

/* Display the main menu */

void main_menu ()
{
	unsigned long choice;

mloop:	if (THREAD_KILL) return;
	printf ("\t     Main Menu\n");
loop:	printf ("\n");
	printf ("\t 1.  Test/Input Data\n");
	printf ("\t 2.  Test/Continue\n");
	printf ("\t 3.  Test/Exit\n");
	printf ("\n");
	printf ("\t 4.  Options/CPU\n");
	printf ("\t 5.  Options/Preferences\n");
	printf ("\t 6.  Advanced/Priority\n");
	printf ("\n");
	printf ("\t 7.  Help/About\n");
	printf ("Your choice: ");
	choice = get_number (0);
	if (choice <= 0 || choice >= 8) {
		printf ("\t     Invalid choice\n");
		goto loop;
	}

/* Display the main menu and switch off the users choice */

	printf ("\n");
	switch (choice) {

/* Test/Primenet dialog */

	case 1:
		test_inputdata ();
		break;

/* Test/Continue */

	case 2:
		linuxContinue ((char*)"Another llrCUDA is running.\n");
		askOK ();
		break;

/* Test/Exit */

	case 3:
		exit(0);  // JP 08/06/17
 //             return;

/* Options/CPU dialog */

	case 4:
		options_cpu ();
		break;

/* Options/Preferences dialog */

	case 5:
		options_preferences ();
		break;

/* Advanced/Priority dialog */

	case 6:
		advanced_priority ();
		break;

/* Help/About */

	case 7:
		help_about ();
		break;

	}
	goto mloop;
}
