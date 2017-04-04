/*
 * from Numerical Recipies in C. (Japanese) in p.222 
 * Poisson Random Generator
 */
 
#include "poissonrand.h"
 
 // Log Gamma Function 
 // Numerical Recipies in C (Japanese) p.170
 double gammln( double xx )
 {
 	double x, y, tmp, ser;
 	static double cof[6] = {76.18009172947146, -86.50532032941677,
 							24.01409824083091, -1.231739572450155,
 							0.1208650973866179e-2, -0.5395239384953e-5};
	int j;
	
	y = x = xx;
	tmp = x + 5.5;
	tmp -= ( x + 0.5 ) * log( tmp );
	ser = 1.000000000190015;
	for( j = 0 ; j < 5 ; j++ )	ser += cof[j] / ++y;
	return( -tmp+log(2.5066282746310005*ser/x) );
 }
 
 double PoissonRand( double xm ) // 母数 xm のポアソン乱数生成
 {
 	static double sq, alxm, g, oldm = (-1.0);
 	double em, t, y;
 	
 	if( xm < 12.0 ){
 		if( xm != oldm ){
 			oldm = xm;
 			g = exp( -xm );
 		}
 		em = -1;
 		t = 1.0;
 		do{
 			++em;
 			t *= genrand_real1();
 		} while( t > g );
 	}
 	else{
 		if( xm != oldm ){
 			oldm = xm;
 			sq = sqrt( 2.0 * xm );
 			alxm = log( xm );
 			g = xm * alxm - gammln( xm + 1.0 );
 		}
 		do{
 			do{
 				y = tan( PI * genrand_real1() );
 				em = sq * y + xm;
 			} while( em < 0.0 );
 			em = floor( em );
 			t = 0.9 * ( 1.0 + y*y ) * exp( em * alxm - gammln( em + 1.0 ) - g );
 		} while( genrand_real1() > t );
 	}
	return em;
}


// #define _POISSONRAND_TEST_ 1

#ifdef _POISSONRAND_TEST_
#include <stdio.h>
#include <time.h>


int main( void )
{
    int i;
    unsigned long seed;
    double m, v, r;
	double lambda;
	
#define NMAX (100*100)
#define SQR( x )  ((x)*(x))
    seed = (unsigned long)time( NULL );
    init_genrand( seed );

    m = v = 0.0;
    lambda = 5.0;
    for( i = 0 ; i < NMAX ; i++ ){
		r = PoissonRand( lambda );   
		printf( "%f\n", r );
		m += r;
		v += SQR( r );
    }
    m /= NMAX;
    v = v/NMAX - SQR(m);
    
    fprintf( stderr, "mean   = %e\n", m );
    fprintf( stderr, "stddev = %e\n", sqrt( v ) );
	return 0;
}
#endif
