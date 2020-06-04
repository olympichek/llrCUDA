/************************** Giants.c extensions ******* J.P. 13/01/11 **************/

/* Include Files. */

#include <ctype.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "giants.h"
#include "giantext.h"

#define BITSINUSHORT 16
#define USHORTM (1<<BITSINUSHORT)

int gmodi (	/* Returns g%i as an integer */
	unsigned long i, giant g) { 
	unsigned long wordweight, j, k, size, value;
	int sign;
	if (i==1 || g->sign==0) return 0;
	if (g->sign < 0) {
	    sign = -1;
	    size = -g->sign;
	}
	else {
	    sign = 1;
	    size = g->sign;
	}
	wordweight = 1;
	value = 0;
	for (j=0; j<size; j++) {
	    value += ((unsigned long)(g->n[j]%i)*wordweight)%i;
	    if (value >= i) value -= i;
	    for (k=1; k<=BITSINUSHORT; k++) {
		wordweight <<=1;
		if (wordweight >= i) wordweight -= i;
	    }
	}
	return sign*value;
}

void ulmulg (unsigned long mult, giant num) {
	giant gmult;
	gmult = newgiant (2);
        if (mult == 0) {
            num->sign = 0;
            return;
        }
	gmult->n[0] = (unsigned short)(mult % USHORTM);
	gmult->n[1] = (unsigned short)(mult >> BITSINUSHORT);
	gmult->sign = (gmult->n[1] == 0) ? 1 : 2;
	mulg (gmult, num);
	free (gmult);
}


void ctog (		/* The giant g is set to string s. */
	char	*s,
	giant	g)
{
	for (g->sign = 0; isdigit (*s); s++) {
		smulg (10, g);
		iaddg (*s - '0', g);
	}
}

void gtoc (		/* The giant g is converted to string s. */
	giant	g,
	char	*s,
	int	sbufsize)
{
	giant	x, y, ten;
	int	i, len;
	char	c;

	assert (g->sign >= 0);

	x = newgiant (g->sign); gtog (g, x);
	y = newgiant (g->sign);
	ten = newgiant (1); itog (10, ten);
	sbufsize--;
	for (len = 0; len < sbufsize && x->sign; len++) {
		gtog (x, y);
		modg (ten, y);
		s[len] = isZero (y) ? '0' : (char) (y->n[0] + '0');	// J.P. 24/01/2011
		divg (ten, x);
	}
	for (i = 0; i < len / 2; i++) {
		c = s[i];
		s[i] = s[len-1-i];
		s[len-1-i] = c;
	}
	s[len] = 0;
	free(x);
	free(y);
	free(ten);
}

void uldivg (unsigned long divisor, giant theg) {
    giant gdivisor = newgiant (2);
    if (divisor == 0)
        gdivisor->sign = 0;
    else {
        gdivisor->n[0] = (unsigned short)(divisor % USHORTM);
        gdivisor->n[1] = (unsigned short)(divisor >> BITSINUSHORT);
        gdivisor->sign = (gdivisor->n[1] == 0) ? 1 : 2;
    }
    divg (gdivisor, theg);
    free (gdivisor);
}

void
power(
	giant		x,
	int 		n
)
/* x becomes x^n. */
{
	giant scratch2 = popg();
	gtog(x, scratch2);
	itog(1, x);
	while (n)
	{
            if (n & 1)
                mulg(scratch2, x);
            n >>= 1;
            if (n)
                squareg(scratch2);
	}
	pushg(1);
}


/************************************************************************************/
