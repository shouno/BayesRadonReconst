#include <stdio.h>
#include <math.h>

int main( void )
{
    double x, xx;
    double xT, lft, rgt, dlt;
    int N;

    N = 4;
    lft = 2.0;
    rgt = -1.0;
    dlt = (rgt-lft)/(N-1);
    xT = (rgt-lft) + dlt;

    printf( "# (lft, rgt) = (%lf, %lf), dlt=%lf, xT =%lf\n", lft, rgt, dlt, xT );
    
    for( x = -8.0 ; x <= 8.0 ; x += 0.1 ){
	xx = lft + fmod( x-lft, xT );
	if( (x-lft)*xT < 0 ){
	    xx += xT;
	}
	printf( "%lf, %lf\n", x, xx );
    }

    return 0;
}


