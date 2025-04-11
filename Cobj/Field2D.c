#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "Field2D.h"


void SetupField( Field2D *f, double l, double r, double t, double b,
		 double data[], int nrow, int ncol, int periodic )
{
    int i;
    double dlt, cur;

    f->left = l;
    f->right = r;
    f->top = t;
    f->bottom = b;
	
    f->ary = (double*)malloc(nrow*ncol*sizeof(double));
    f->xpos = (double*)malloc(ncol*sizeof(double));
    f->ypos = (double*)malloc(nrow*sizeof(double));
    memmove( f->ary, data, nrow*ncol*sizeof(double) );
    f->ncol = ncol;
    f->nrow = nrow;

    if( periodic == 1 )		f->isPeriodic = 1;
    else			f->isPeriodic = 0;

    if( ncol == 1 )	f->xpos[0] = 0.0; // 原点のみ
    else{
	cur = l;
	dlt = ( r - l ) / (ncol-1);
	for( i = 0 ; i < ncol ; i++ ){
	    f->xpos[i] = cur;
	    cur += dlt;
	}
	f->xT = ncol * dlt;
    }
    if( nrow == 1 )	f->ypos[0] = 0.0;
    else{
	cur = t;
	dlt = ( b - t ) / (nrow-1);
	for( i = 0 ; i < nrow ; i++ ){
	    f->ypos[i] = cur;
	    cur += dlt;
	}
	f->yT = - nrow * dlt;
    }
}


void DestroyField( Field2D f )
{
    free( f.ary );
    free( f.xpos );
    free( f.ypos );
}


//
// These functions are private functions:
// field position x, y should be in original area,
// even if the isPeriodic flag is true.
//
static int getRowPos( double y, Field2D f );
static int getColPos( double x, Field2D f );

static int getRowPos( double y, Field2D f )
{
    //
    // The tics vector should be arranged in reducing.  i.e. ypos[0] > ypos[nr-1]
    //
    int low, mid, high;
    int nr = f.nrow;

    if( f.ypos[0] > f.ypos[nr-1] ){
	high = 0;
	low = nr-1;
	if( y <  f.ypos[nr-1] )		return( nr );
	if( y >= f.ypos[0] )		return( 0 );
    }
    else{
	high = nr-1;
	low = 0;
	if( y <  f.ypos[0] )		return( -1 );
	if( y >= f.ypos[nr-1] )		return( nr-1 );
    }
    while( abs( low - high ) > 1 ){
	mid = ( low + high ) / 2;
	if( f.ypos[mid] <= y )		low = mid;
	else				high = mid;
    }
    return( low );

}


static int getColPos( double x, Field2D f )
{
    int low, mid, high;
    int nc = f.ncol;

    if( f.xpos[0] > f.xpos[nc-1] ){
	high = 0;
	low = nc-1;
	if( x <  f.xpos[nc-1] )		return( nc );
	if( x >= f.xpos[0] )		return( 0 );
    }
    else{
	high = nc-1;
	low = 0;
	if( x <  f.xpos[0] )		return( -1 );
	if( x >= f.xpos[nc-1] )		return( nc-1 );
    }
    while( abs( high - low ) > 1 ){
	mid = ( low + high ) / 2;
	if( f.xpos[mid] <= x )		low = mid;
	else				high = mid;
    }
    return( low );
}


int inArea( double x, double y, Field2D f )
{
    double l, r, t, b;

    l = f.left;
    r = f.right;
    t = f.top;
    b = f.bottom;

    if( f.isPeriodic ){
	return 1;
    }
    if( ( l <= x && x <= r ) &&
	( b <= y && y <= t ) ){
	return 1;
    }
    else
	return 0;
}


double FieldValue( double x, double y, Field2D f )
{
    int c, r, c1, r1;
    int nr = f.nrow;
    int nc = f.ncol;
    double dx, dy, s, t;
    double f1, f2, f3, f4;

    if( f.isPeriodic ){
	x = fmod( x-f.left, f.xT ) + f.left;
	if( (x-f.left)*f.xT < 0.0 )	x += f.xT;
	c = getColPos( x, f );

//	printf( "x = %lf, c = %d, ", x, c );

	if( f.xpos[0] < f.xpos[nc-1] ){
	    dx = f.xpos[1] - f.xpos[0];
	    if( c < 0 ){
		c1 = 0;		c = nc-1;
		s = 1 - (f.xpos[c1] - x) / dx;
	    }
	    else if( c >= nc-1 ){
		c = nc-1;	c1 = 0;
		s = (x - f.xpos[c]) / dx;
	    }
	    else{
		c1 = c + 1;
		s = ( x - f.xpos[c] ) / dx;
	    }
	}
	else{
	    dx = f.xpos[0] - f.xpos[1];
	    if( c > nc-1 ){
		c1 = nc-1;	c = 0;
		s = 1 - (f.xpos[c1] - x) / dx;
	    }
	    else if( c <= 0 ){
		c = nc-1;	c1 = 0;
		s = ( x - f.xpos[c] ) / dx;
	    }
	    else{
		c1 = c - 1;
		s = ( x - f.xpos[c] ) / dx;
	    }
	}

	y = fmod( y-f.bottom, f.yT ) + f.bottom;
	if( (y-f.bottom)*f.yT < 0.0  )	y += f.yT;
	r = getRowPos( y, f );

//	printf( "f.yT = %f, bottom = %f, y = %lf, r = %d, ", f.yT, f.bottom, y, r );

	if( f.ypos[0] < f.ypos[nr-1] ){
	    dy = f.ypos[1] - f.ypos[0];
	    if( r < 0 ){
		r = 0;		r1 = nr-1;
		t = 1 - (f.ypos[r1] - y) / dy;
	    }
	    else if( r >= nr-1 ){
		r1 = 0;		r = nr-1;
		t = (y - f.ypos[r])/dy;
	    }
	    else{
		r1 = r + 1;
		t = ( y - f.ypos[r] ) / dy;
	    }
	}
	else{
	    dy = f.ypos[0] - f.ypos[1];
	    if( r > nr-1 ){
		r1 = nr-1;	r = 0;
		t = 1 - (f.ypos[r1] - y) / dy;
	    }
	    else if( r <= 0 ){
		r = nr-1;	r1 = 0;
		t = ( y - f.ypos[r] ) / dy;
	    }
	    else{
		r1 = r - 1;
		t = ( y - f.ypos[r] ) / dy;
	    }
	}
    }
    else{
	c = getColPos( x, f );
	if( f.xpos[0] < f.xpos[nc-1] ){
	    s = 0.0;
	    if( c < 0 )			c1 = c = 0;
	    else if( c >= nc-1 )	c1 = c = nc-1;
	    else{
		c1 = c + 1;
		dx = f.xpos[c1] - f.xpos[c];
		s = ( x - f.xpos[c] ) / dx;
	    }
	}
	else{
	    s = 1.0;
	    if( c > nc-1 )		c1 = c = nc-1;
	    else if( c <= 0 )		c1 = c = 0;
	    else{
		c1 = c - 1;
		dx = f.xpos[c1] - f.xpos[c];
		s = ( x - f.xpos[c] ) / dx;
	    }
	}
	
	r = getRowPos( y, f );
	if( f.ypos[0] < f.ypos[nr-1] ){
	    t = 0.0;
	    if( r < 0 )			r1 = r = 0;
	    else if( r >= nr-1 )	r1 = r = nr-1;
	    else{
		r1 = r + 1;
		dy = f.ypos[r1] - f.ypos[r];
		t = ( y - f.ypos[r] ) / dy;
	    }
	}
	else{
	    t = 1.0;
	    if( r > nr-1 )		r1 = r = nr-1;
	    else if( r <= 0 )		r1 = r = 0;
	    else{
		r1 = r - 1;
		dy = f.ypos[r1] - f.ypos[r];
		t = ( y - f.ypos[r] ) / dy;
	    }
	}
	// cout << "( r, c ) = " << r << ", " << c << std::endl;
	// cout << "( s, t ) = " << s << ", " << t << std::endl;
    }

    f1 = f.ary[r  * nc + c ];
    f2 = f.ary[r1 * nc + c ];
    f3 = f.ary[r1 * nc + c1];
    f4 = f.ary[r  * nc + c1];

    return( (1-s) * (1-t) * f1 +
	    (1-s) *    t  * f2 +
	    s     *    t  * f3 +
	    s     * (1-t) * f4 );
}



#ifdef _TEST_FIELD2D_
int main( void )
{
    int i, j;
    int Nx, Ny;
    double l, r, t, b;
    double *dtary, x, y;
    Field2D f;
	
    Nx = Ny = 64;
    l = -1.0;
    r = 1.0;
    t = 1.0;
    b = -1.0;

    dtary = (double*)malloc(sizeof(double)*Nx*Ny);
    memset( dtary, 0, Nx*Ny*sizeof(double) );
    for( i = 10 ; i < 20 ; i++ ){
	for( j = 15 ; j < 25 ; j++ ){
	    dtary[i*Nx+j] = 1.0;
	}
    }
    for( i = 28 ; i < 32 ; i++ ){
	for( j = 12 ; j < 18 ; j++ ){
	    dtary[i*Nx+j] = 1.0;
	}
    }

    SetupField( &f, l, r, t, b, dtary, Ny, Nx, 1 );

    for( y = -3.0 ; y <= 3.0 ; y += 0.1 ){
	for( x = -3.0 ; x <= 3.0 ; x += 0.1 ){
	    printf( "%f, %f, %f\n", x, y, FieldValue( x, y, f ) );
	}
	printf( "\n" );
    }




    DestroyField( f );
    free( dtary );
	
    return 0;
}

#endif /* _TEST_ */

