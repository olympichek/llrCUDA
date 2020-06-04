#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "./giants.h"

#define  ULONG unsigned long
#define MAXULONG 0xFFFFFF00

void trace(int);
EXTERNC int kronecker (uint32_t, uint32_t);
EXTERNC ULONG gcd (ULONG, ULONG);

void Reduce (uint32_t x, uint32_t *d, uint32_t *b) {

// Reduce a Discriminant to a square free integer

	uint32_t div, sq;

	*d = x;
	*b = 1;

	while (!((*d)%4)) {
		*d /= 4;
		*b *= 2;
	}

	for (div = 3; (sq = div*div) <= *d; div += 2)
		while (!((*d)%sq)) {
			*d /= sq;
			*b *= div;
		}
}

uint32_t issquare (uint32_t n) {
	uint32_t s;
	s = (uint32_t)floor(sqrt(n));
	if (s*s == n)
		return s;
	else
		return 0;
}

uint32_t twopownmodm (uint32_t n, uint32_t m, uint32_t *order, uint32_t *nmodorder) {
	uint32_t tpnmodm, tp, i;
	tpnmodm = 0;
	*order = tp = 1;
	for (i=1; i<m; i++) {
	    tp <<= 1;
	    if (tp >= m) tp -= m;
	    if (i==n) tpnmodm = tp;
	    if (tp==1) {
		*order = i;
		break;
	    }
	}
	*nmodorder = n%*order;
	if (!tpnmodm) {
	    for (i=1, tp=1; i<=*nmodorder; i++) {
		tp <<= 1;
		if (tp >= m) tp -= m;
	    }
	    tpnmodm = tp;
	}
	return tpnmodm;
}

uint32_t Bachet(uint32_t u, uint32_t v, long *a, long *b) {

//  Computes a and b such as a*u+b*v = gcd(u,v),
//  returns gcd(u,v).

    uint32_t n=0, m11=1, m12=0, m21=0, m22=1;
    uint32_t q, newm11, newm21, newu;

    while (v!=0) {
	q = u/v;
	newm11 = m11*q+m12;
	m12 = m11;
	m11 = newm11;
	newm21 = m21*q+m22;
	m22 = m21;
	m21 = newm21;
	newu = v;
	v = u-q*v;
	u = newu;
	n++;
    }
    if (n&1) {
	*a = -(int)m22;
	*b = (int)m12;
    }
    else {
	*a = (int)m22;
	*b = -(int)m12;
    }
    return u;
}

int gen_v1(giant k, uint32_t n, int debug, uint32_t maxres) {
    long sign, j2N, jNd, jNa, v, vcand;
    long idred, iared, iorderd, iordera, kdiff, ndiff;
    uint32_t nbres, kmod3, rawd, d, dred, kmodd, tpnmd, i, orderd, Nmodd;
    uint32_t V, W, dsqW, aplus, aminus, b, bcand, rplus, rminus, found, iseps2;
    uint32_t nmodorderd, ared, kmoda, tpnma, ordera, nmodordera, Nmoda;
    uint32_t dk, dn, kmod, nmod, modulk, moduln;
    nbres = 0;
    if (n>2)
	j2N = 1;
    else
	j2N = -1;
    kmod3 = gmodi (3, k);
    if (maxres == 1 && kmod3)		// Consider only the simple case !
	if ((kmod3 == 1 && !(n&1)) || (kmod3 == 2 && (n&1))) {
	    if (debug)
		printf ("d = 3 divides N !\n");
	    return (-3);
	}
	else {
	    if (debug) {
		printf ("epsilon = 2+sqrt(3)\n");
		printf ("k = %lu (mod 3), n = %lu (mod 2)\n", kmod3, n&1);
		printf ("v1 = 4, d = 3, a = 6, b = 2, r = 24, +1\n");
		printf ("v1 = 4, d = 3, a = 2, b = 2, r = 8, -1\n");
	    }
	    return (4);
	}
    for (v=3; (rawd = v*v-4)<MAXULONG; v++) {
	Reduce (rawd, &d, &b);
	if (!(d&1)) {
	    sign = j2N;
	    dred = d>>1;
	}
	else {
	    sign = 1;
	    dred = d;
	}
	if (d==2) {
	    if (n > 2) continue;
	}
	else {
	  kmodd = gmodi (dred, k);
	  if (!kmodd) continue;
	  if (n>1 && (((dred-1)/2) & 1)) sign = - sign;
	  tpnmd = twopownmodm (n, dred, &orderd, &nmodorderd);
	  Nmodd = (kmodd*tpnmd-1)%dred;
	  if (!Nmodd) {
	    if (debug)
		if (d&1)
		    printf ("d = %lu divides N !\n", dred);
		else
		    printf ("d/2 = %lu divides N !\n", dred);
	    return (-(int)dred);
	  }
	  if ((jNd = jacobi(Nmodd,dred)) > 1) {
	    if (debug)
		printf ("%ld divides d = %lu and N !\n", jNd, d);
	    return (-jNd);
	  }
	  if ((sign*jNd) != -1) continue;
	}
	found = 0;
	for (W=1; (dsqW=d*W*W) < MAXULONG; W++)
	    if (dsqW>4 && (V=issquare(dsqW-4))) {
		vcand = dsqW-2;
		bcand = V*W;
		iseps2 = 1;
		found = 1;
		break;
	    }
	    else if (V=issquare(dsqW+4)) {
		vcand = V;
		bcand = W;
		iseps2 = 0;
		found = 1;
		break;
	    }
	if (found && v==vcand && b==bcand) {
	    aplus = v+2;
	    aminus = v-2;
	    rplus = 4*aplus;
	    rminus = 4*aminus;
	    if (iseps2) {
		nbres++;
		if (debug) {
		    if (V&1 || W&1)
			if (W != 1)
			    printf ("epsilon = [%lu+%lu*sqrt(%lu)]/2\n", V, W, d);
			else
			    printf ("epsilon = [%lu+sqrt(%lu)]/2\n", V, d);
		    else
			if (W/2 != 1)
			    printf ("epsilon = %lu+%lu*sqrt(%lu)\n", V/2, W/2, d);
			else
			    printf ("epsilon = %lu+sqrt(%lu)\n", V/2, d);
		    printf ("k = %lu (mod %lu), n = %lu (mod %lu)\n", kmodd, dred, nmodorderd, orderd);
		    printf ("v1 = %ld, d = %lu, a = %lu, b = %lu, r = %lu, +1,eps2\n",
				v, d, aplus, b, rplus);
		    printf ("v1 = %ld, d = %lu, a = %lu, b = %lu, r = %lu, -1,eps2\n",
				v, d, aminus, b, rminus);
		}
		if (nbres == maxres) return v;
	    }
	    else {
		for (ared=aplus, i=0; !(ared & 1); ared >>= 1) i++;
		if (i&1)
		    sign = j2N;
		else
		    sign = 1;
		kmoda = gmodi (ared, k);
		if (!kmoda) continue;
		if (n>1 && (((ared-1)/2) & 1)) sign = - sign;
		tpnma = twopownmodm (n, ared, &ordera, &nmodordera);
		Nmoda = (kmoda*tpnma-1)%ared;
		if (!Nmoda) {
	 	   if (debug)
			printf ("a/%lu = %lu divides N !\n", aplus/ared, ared);
		    return (-(int)ared);
		}
		if ((jNa = jacobi(Nmoda,ared)) > 1) {
		    if (debug)
			printf ("%ld divides a = %lu and N !\n", jNa, aplus);
		    return (-jNa);
		}
		if ((sign*jNa) != -1) continue;
		nbres++;
		if (debug) {
		    if (V&1 || W&1)
			if (W != 1)
			    printf ("epsilon = [%lu+%lu*sqrt(%lu)]/2\n", V, W, d);
			else
			    printf ("epsilon = [%lu+sqrt(%lu)]/2\n", V, d);
		    else
			if (W/2 != 1)
			    printf ("epsilon = %lu+%lu*sqrt(%lu)\n", V/2, W/2, d);
			else
			    printf ("epsilon = %lu+sqrt(%lu)\n", V/2, d);
		    dk = Bachet (dred, ared, &idred, &iared);
		    modulk = dred*(ared/dk);		// lcm(dred, ared).
		    kdiff = kmoda-kmodd;
		    if (kdiff >=0) {
			while (idred<0)
			    idred += ared;
			kmod = (kmodd + idred*(kdiff/dk)*dred)%modulk;
		    }
		    else {
			kdiff = -kdiff;
			while (iared<0)
			    iared += dred;
			kmod = (kmoda + iared*(kdiff/dk)*ared)%modulk;
		    }
		    dn = Bachet (orderd, ordera, &iorderd, &iordera);
		    moduln = orderd*(ordera/dn);	// lcm(orderd, ordera).
		    ndiff = nmodordera-nmodorderd;
		    if (ndiff >=0) {
			while (iorderd<0)
			    iorderd += ordera;
			nmod = (nmodorderd + iorderd*(ndiff/dn)*orderd)%moduln;
		    }
		    else {
			ndiff = -ndiff;
			while (iordera<0)
			    iordera += orderd;
			nmod = (nmodordera + iordera*(ndiff/dn)*ordera)%moduln;
		    }
		    printf ("k = %lu (mod %lu), n = %lu (mod %lu)\n", kmod, modulk, nmod, moduln);
		    printf ("v1 = %ld, d = %lu, a = %lu, b = %lu, r = %lu, +1\n",
				v, d, aplus, b, rplus);
		    printf ("v1 = %ld, d = %lu, a = %lu, b = %lu, r = %lu, -1\n",
				v, d, aminus, b, rminus);
		}
		if (nbres == maxres) return v;
	    }
	}
    }
    return -1;		// Unable to find a value for v...
}

int genProthBase(giant k, uint32_t n) {
	uint32_t Nmodp, kmodp, p, tpnmp, orderp, nmodorderp, kw;
	int jNp;

//	Return the least prime p such Jacobi (N, p) = -1

	if (k->sign == 1 && n < 3) {	//	Eliminate some trivial cases
		kw = k->n[0];
		if (n == 1 && kw == 1)
			return (2);
		else if (n == 2)
			return (2);
		else
			return (-1);
	}
	else {							// General case
		for (p = 3; p<=2147483647; p+=2) {
			if (!isPrime(p))
				continue;
			kmodp = gmodi (p, k);
			if (!kmodp)
				continue;
			tpnmp = twopownmodm (n, p, &orderp, &nmodorderp);
			Nmodp = (kmodp*tpnmp+1)%p;
			if (!Nmodp) {
				return (-(int)p);
			}
			if ((jNp = jacobi(Nmodp, p)) > 1) {
				return (-jNp);
			}
			if (jNp != -1)
				continue;
			return (p);
		}
		return (-1);
	}
}

int genProthBase1(giant N)
{
	uint32_t NmodD, D, dred, Nmod8;
	int jNp, chgsign;

//	Return the least D such as kronecker (D, N) = -1

//	Rem : for n>1, (N-1)/2 = k*2^(n-1) is even, so (D / N) = (N / D)

//  Then, if D = 2^s*t, (N / D) = (N / 2)^s * (N / t) = (N / 2)^s * (Nmodt / t)

	Nmod8 = N->n[0] & 7;

	for (D = 2; D<=2147483647; D++) {
		dred = D;
		chgsign = 1;
		while (!(dred&1)) {
			dred >>= 1;					// Compute the odd part of D
			if (Nmod8 == 3 || Nmod8 == 5)
				chgsign = -chgsign;
		}
		if (dred == 1)
			jNp = 1;
		else {
			NmodD = gmodi (dred, N);
			if (!NmodD)
				return (-(int)dred);
			if ((jNp = kronecker(NmodD, dred)) > 1)
				return (-jNp);
		}
		if ((jNp*chgsign) != -1)
			continue;
		return ((int)D);
		}
		return (-1);
}

int genLucasBaseQ(giant N, uint32_t DARG) {
	uint32_t NmodD, dred, Nmod8, D;
	int jNp, chgsign;

//	Return the least D = 1+4*Q such as kronecker (D, N) = -1

//  if D = 2^s*t, (N / D) = (N / 2)^s * (N / t) = (N / 2)^s * (Nmodt / t)

	Nmod8 = N->n[0] & 7;

	for (D = DARG; D<=2147483647; D+=4) {
		dred = D;
		chgsign = 1;
		while (!(dred&1)) {
			dred >>= 1;					// Compute the odd part of D
			if (Nmod8 == 3 || Nmod8 == 5)
				chgsign = -chgsign;
		}
		if (dred == 1)
			jNp = 1;
		else {
			NmodD = gmodi (dred, N);
			if (!NmodD)
				return (-(int)dred);
			if ((jNp = kronecker(NmodD, dred)) > 1)
				return (-jNp);
		}
		iaddg (-1, N);					// Compute N-1
		if (((dred-1) & 2) && (N->n[0] & 2))	// Quadratic reciprocity
			chgsign = -chgsign;
		iaddg (1, N);					// Restore N
		if ((jNp*chgsign) != -1)
			continue;
		return ((int)D);
	}
	return (-1);
}

int genLucasBaseP(giant N, uint32_t PARG) {
	uint32_t NmodD, D, dred, Nmod8, P;
	int jNp, chgsign;

//	Return the least P such as D = P^2-4 and kronecker (D, N) = -1

//  if D = 2^s*t, (N / D) = (N / 2)^s * (N / t) = (N / 2)^s * (Nmodt / t)

	Nmod8 = N->n[0] & 7;

	for (P = PARG; P*P<=2147483647; P++) {
		D = P*P-4;
		dred = D;
		chgsign = 1;
		while (!(dred&1)) {
			dred >>= 1;					// Compute the odd part of D
			if (Nmod8 == 3 || Nmod8 == 5)
				chgsign = -chgsign;
		}
		if (dred == 1)
			jNp = 1;
		else {
			NmodD = gmodi (dred, N);
			if (!NmodD)
				return (-(int)dred);
			if ((jNp = kronecker(NmodD, dred)) > 1)
				return (-jNp);
		}
		iaddg (-1, N);					// Compute N-1
		if (((dred-1) & 2) && (N->n[0] & 2))	// Quadratic reciprocity
			chgsign = -chgsign;
		iaddg (1, N);					// Restore N
		if ((jNp*chgsign) != -1)
			continue;
		return ((int)P);
	}
	return (-1);
}

long generalLucasBase(giant N, uint32_t *PARG, uint32_t *Q) {
	uint32_t NmodD, D, dred, Nmod8, NmodPQD, gcdNPQD, *P;
	int jNp, chgsign;

//	Return the least D = P^2-4*Q such as kronecker (D, N) = -1

//  if D = 2^s*t, (N / D) = (N / 2)^s * (N / t) = (N / 2)^s * (Nmodt / t)

	Nmod8 = N->n[0] & 7;

	for (P = PARG; (*P)*(*P)<=2147483647; (*P)++) {
		for ((*Q)=2; 4*(*Q)<(*P)*(*P); (*Q)++) {
			D = (*P)*(*P)-4*(*Q);
			if ((uint32_t)(floor(sqrt ((double)D)) * floor(sqrt ((double)D))) == D) {
				continue;
			}
			dred = D;
			chgsign = 1;
			while (!(dred&1)) {
				dred >>= 1;					// Compute the odd part of D
				if (Nmod8 == 3 || Nmod8 == 5)
					chgsign = -chgsign;
			}
			if (dred == 1)
				jNp = 1;
			else {
				NmodD = gmodi (dred, N);
				if (!NmodD)
					return (-(long)dred);
				if ((jNp = kronecker(NmodD, dred)) > 1)
					return (-jNp);
			}
			iaddg (-1, N);					// Compute N-1
			if (((dred-1) & 2) && (N->n[0] & 2))	// Quadratic reciprocity
				chgsign = -chgsign;
			iaddg (1, N);					// Restore N
			if ((jNp*chgsign) != -1)
				continue;
			NmodPQD = gmodi ((*P)*(*Q)*D, N);
			gcdNPQD = gcd (NmodPQD, (*P)*(*Q)*D);
			if (gcdNPQD != 1)
				return (-(long)gcdNPQD);
			return (D);
		}
	}
	return (-1);
}
