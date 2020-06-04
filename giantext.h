/**************************************************************
 *
 *	giantext.h
 *
 *	Header file for extension of library giants.c.
 *
 **************************************************************/

#ifndef _GIANTEXT_HEADER_
#define _GIANTEXT_HEADER_


/**************************************************************
 *
 * Function Prototypes
 *
 **************************************************************/

/* stack handling functions. */
giant	popg(void);
void		pushg(int);


/**************************************************************
 *
 * Conversion to or from strings Routines
 *
 **************************************************************/

/* Convert a positive giant in a decimal string. */

void gtoc (		/* The giant g is converted to string s. */
	giant	g,
	char	*s,
	int	sbufsize
);

/* Convert a decimal string into a giant. */

void ctog (		/* The giant g is set from string s. */
	char	*s,
	giant	g
);

/**************************************************************
 *
 * Math Functions
 *
 **************************************************************/

void ulmulg (	/* num *= mult */
	unsigned long mult,
	giant num
);


int gmodi (		/* Returns g%i as an integer */ 
	unsigned long i,
	giant g
);

void uldivg (unsigned long divisor, giant theg);

void power (giant x, int n);

#endif
