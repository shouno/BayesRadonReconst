/*
 * $Id: normaltest.c,v 1.1 2003/10/21 06:47:13 shouno Exp shouno $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mtrand.h"
#include "gaussrand.h"


#define EPSILON 1.0e-15

#ifdef _GAUSSRAND_TEST_
int main( void )
{
    int i;
    unsigned long seed;
    double m, v, r;

#define NMAX (100*100*100)
#define NMAX_1 (100*100)
#define SQR(x)  ((x)*(x))

    seed = (unsigned long)time( NULL );

    init_genrand( seed );

    m = v = 0.0;
    for( i = 0 ; i < NMAX ; i++ ){

	r = sqrt(6.0*0.5/ ((double)NMAX_1)) * gausdev() + ((double)NMAX_1) * 3;

//	printf( "r : %f\n", r );
	m += r;
	v += SQR( r );
    }
    m /= NMAX;
    v = v/NMAX - SQR(m);

    printf( "mean   = %e\n", m );
    printf( "stddev = %e\n", sqrt( v ) );

    return 0;
}
#endif

//
// Modified Numerical Recipies in C P.217
//
double gausdev( void )
{
    // generate random variables yields to N(0,1)
    static int iset = 0;
    static float gset;
    double fac, rsq, v1, v2;

    if( ! iset ){
		do{
			v1 = 2.0 * genrand_real1() - 1.0;
	    	v2 = 2.0 * genrand_real1() - 1.0;
		    rsq = v1*v1 + v2*v2;
		} while( rsq >= 1.0 || rsq < EPSILON );
		fac =  sqrt( -2.0 * log(rsq) / rsq );
		gset = v1 * fac;
		iset = 1;
		return( v2 * fac );
    }
    else{
		iset = 0;
		return( gset );
    }
}


//
//
//
double NormalRand( double m, double sd )
{
	return( m + sd * gausdev() );
}
