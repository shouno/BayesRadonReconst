#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "Field2D.h"
#include "RadonTrans.h"
#include "gaussrand.h"
#include "poissonrand.h"

#ifndef PI 
#define PI 3.141592653589793
#endif


void ReconstructImage(  BackParams p )
{
    int i, j;
    double dx, dy, x, y;
    dx = (p.r - p.l)/(p.Nx-1);
    dy = (p.b - p.t)/(p.Ny-1);

    y = p.t;
    for( i = 0 ; i < p.Ny ; i++ ){
	x = p.l;
	for( j = 0 ; j < p.Nx ; j++ ){
	    p.out[i*p.Nx+j] = ReconstuctIntegral( x, y, p );
	    x += dx;
	}
	y += dy;
    }
}


double ReconstuctIntegral( double x, double y, BackParams p )
{
    int i, k;
    double dth = PI/p.Nth;
    double Smin = -(p.Ns/2)*p.ds;
    double Smax = ((p.Ns-1)/2)*p.ds;
    double sk, q;
    double ans = 0.0;

    for( i = 0 ; i < p.Nth ; i++ ){
	double th = i * dth;
	double s;
	s = x * cos(th) + y * sin(th);
	if( ( Smin <= s ) && ( s < Smax ) ){
	    k = (int)( (s-Smin)/p.ds );
	    sk = Smin + k * p.ds;
	    q = (s - sk) / p.ds;
	    ans += (1-q) * p.g[i*p.Ns+k] + q * p.g[i*p.Ns+k+1];
	}
    }
    return ans;
}


//
// Noise settings:
// The 'noiseflag' changes Callback functions
// The callback function should have void arguments, 
// thus, the following parameters are used for setting.
// In Normal distribution, the 1st parameter rparam1 
// is used as mean, the 2nd as a standard deviation.
// In the Delta, and Poission distribution, 
// only 1st parameter is used.
//
static double rparam1;
static double rparam2;

void SetRandParams( ForwParams p )
{
    rparam1 = rparam2 = 0.0;
    if( p.noiseflag == 1 ){ // Normal Random
	rparam2 = p.psd;
    }
}


double Delta( void )
{
    return rparam1;
}

double Poisson( void )
{
    return( PoissonRand( rparam1 ) );
}

double Normal( void )
{
    return( NormalRand( rparam1, rparam2 ) );
}

double LineIntegral( double s, double theta, Field2D F, 
		     double dh, double Lim,
		     double (*rfunc)(void) )
{
    double plus, minus;
    double x, y, t;

    plus = 0.0;
    for( t = 0.0 ; t <= Lim ; t += dh ){
	x = s * cos( theta ) - t * sin( theta );
	y = s * sin( theta ) + t * cos( theta );
	if( inArea( x, y, F ) ){
	    rparam1 = FieldValue( x, y, F );
	    // rparam2 is set by SetRandomFunc()
	    plus += (*rfunc)() * dh;
	}
	else{
	    break;
	}
    }
    minus = 0.0;
    for( t = 0.0 ; t >= -Lim ; t -= dh ){
	x = s * cos( theta ) - t * sin( theta );
	y = s * sin( theta ) + t * cos( theta );
	if( inArea( x, y, F ) ){
	    rparam1 = FieldValue( x, y, F );
	    // rparam2 is set by SetRandomFunc()
	    minus += (*rfunc)() * dh;
	}
	else{
	    break;
	}
    }
    return plus + minus;
}




void RadonTrans( ForwParams p )
{
    int i, j;
    Field2D F;
    double *Tau;
    double ds =  p.ds;
    double dtheta = (PI-0) / p.Nth;
    double theta = 0.0;
    double dh = p.Lim / p.Nh; 
    double (*rfunc)(void);

    SetupField( &F, p.l, p.r, p.t, p.b, p.in, p.Ny, p.Nx, p.periodic );
    Tau = p.out;

    SetRandParams( p );
    if(p.noiseflag == 0){
	rfunc = Delta;
    }
    else if( p.noiseflag == 1 ){
	rfunc = Normal;
    }
    else if( p.noiseflag == 2 ){
	rfunc = Poisson;
    }
    else{
	fprintf( stderr, "Illegal Noise flag\n" );
	return;
    }

    for( i = 0 ; i < p.Nth ; i++ ){
	// double s = p.Rmin + ds/2; // j = 0 での s 座標
	/* 下のようにしないと k = 0 での座標が x-y 座標の原点に対応しない */
	double s = -p.Ns/2 * ds; // j = 0 での s 座標 

	// Radon Transform
	for( j = 0 ; j < p.Ns ; j++ ){
	    Tau[i*p.Ns+j] = LineIntegral( s, theta, F, dh, p.Lim, rfunc );
	    s += ds;
	}
	theta += dtheta;
    }
    DestroyField( F );
}



void ReconstructInterface( double out[], int *ny, int *nx,
			   double *l, double *r, double *t, double *b,
			   double in[], int *nth, int *ns, double *ds )
{
    BackParams p;

    p.Ny = *ny;
    p.Nx = *nx;
    p.out = (double*)malloc( sizeof(double)*p.Ny*p.Nx);
    p.l = *l;
    p.r = *r;
    p.t = *t;
    p.b = *b;

    p.Nth = *nth;
    p.Ns = *ns;
    p.g = (double*)malloc( sizeof(double)*p.Nth*p.Ns );
    p.ds = *ds;


//	printf( "Ny, Nx = %d, %d\n", p.Ny, p.Nx );
//	printf( "rect (%f, %f) - (%f, %f)\n", p.l, p.t, p.r, p.b );
//	printf( "Nth, Ns = %d, %d, ds = %f\n", p.Nth, p.Ns, p.ds );
//	printf( "in: %x ----> out: %x\n", in, out );

    memmove( p.g, in, sizeof(double)*p.Nth*p.Ns );
	
    ReconstructImage( p );

    memmove( out, p.out, sizeof(double)*p.Ny*p.Nx );
    free( p.out );
    free( p.g );
}

			   
void RadonInterface( double out[], int *nth, int *ns, double *ds,
		     double in[], int *ny, int *nx,
		     double *l, double *r, double *t, double *b,
		     double *pmu, double *psd, int *nh,
		     int *pflag, int *periodic )
{
    ForwParams p;
	
    // Setup parameters
    p.Nth = *nth;
    p.Ns = *ns;
    p.ds = *ds;
    p.out = (double*)malloc( p.Nth * p.Ns * sizeof(double) );

    p.Nx = *nx;
    p.Ny = *ny;
    p.in = (double*)malloc( p.Nx * p.Ny * sizeof(double) );
    memmove( p.in, in, p.Nx*p.Ny * sizeof(double) );
    p.l = *l;
    p.r = *r;
    p.t = *t;
    p.b = *b;

//	printf( "Nth, Ns = %d, %d, ds = %f\n", p.Nth, p.Ns, p.ds );
//	printf( "Ny, Nx = %d, %d\n", p.Ny, p.Nx );
//	printf( "rect (%f, %f) - (%f, %f)\n", p.l, p.t, p.r, p.b );

    p.noiseflag = *pflag;
    p.pmu = *pmu;
    p.psd = *psd;

    // Default parameters
    // p.Nh = 4096;   
    p.Nh = *nh;
    p.Lim = (p.Ns*p.ds)*sqrt(2)/2;  // 積分の上限値  t in [ -Lim, Lim ] で
    p.periodic = *periodic;

    RadonTrans( p );

    memmove( out, p.out, p.Ns*p.Nth*sizeof(double) );

    free( p.out );
    free( p.in );
}





#ifdef _TEST_RADONTRANS_
int main( void )
{
    int i, j;
    int Nx, Ny;
    int Nth, Ns;
    int Nh;
    double l, r, t, b;
    double Rmax, Rmin, ds;
    double *inary;
    double *outary;
    double pmu, psd;
    int pflag = 0;
    int periodic = 1;

    Nh = 4096;
    Nx = Ny = 64;
    inary = (double*)malloc(sizeof(double)*Nx*Ny);
    memset( inary, 0, Nx*Ny*sizeof(double) );
    for( i = 10 ; i < 20 ; i++ ){
	for( j = 15 ; j < 25 ; j++ ){
	    inary[i*Nx+j] = 10.0;
	}
    }
    for( i = 28 ; i < 32 ; i++ ){
	for( j = 12 ; j < 18 ; j++ ){
	    inary[i*Nx+j] = 5.0;
	}
    }
    l = b = -1.0;
    t = r = 1.0;

    Nth = Ns = 64;
    Rmax = 1.0;
    Rmin = -1.0;
    ds = (Rmax-Rmin)/Ns;
    outary = (double*)malloc(sizeof(double)*Ns*Nth);
    memset( outary, 0, Ns*Nth * sizeof(double) );

    pmu = 0.0;
    psd = 1.0;
    pflag = 0;  // normal N( val, psd );
    RadonInterface( outary, &Nth, &Ns, &ds,
		    inary, &Ny, &Nx, &l, &r, &t, &b,
		    &pmu, &psd, &Nh,
		    &pflag, &periodic );

    for( i = 0 ; i < Nth ; i++ ){
	for( j = 0 ; j < Ns ; j++ ){
	    printf( "%f, %f, %f\n", 
		    (double)i, (double)j, outary[i*Ns+j] );
	}
	printf( "\n" );
    }

    return 0;
}

#endif /* _TEST_ */
