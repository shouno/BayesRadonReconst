#ifndef __RADONTRANS_H__
#define __RADONTRANS_H__ 1

typedef struct _ForwParams{
    double *out;    // size should be Nth * Ns
    int Nth;
    int Ns;
    double ds;	//s の刻み. -Ns/2 が Rmin に対応

    double *in;     // size shouuld be Nx * Ny
    int Nx;
    int Ny;
    double l;	// assinged to (l, t) - (r, b) rectangle
    double r;
    double t;
    double b;

    int Nh;		// 積分きざみ
    double Lim;		// 有効積分範囲
    int periodic;	// 境界条件

    int noiseflag;	// noise adding flag
    double pmu;		// position noise mean
    double psd;		// position noise stddev
} ForwParams;


typedef struct _BackParams{
    double *out;	// size shouuld be Nx * Ny
    int Nx;
    int Ny;
    double l;	// assinged to (l, t) - (r, b) rectangle
    double r;
    double t;
    double b;

    double *g;	// size should be Nth * Ns
    int Nth;
    int Ns;
    double ds;	// s の刻み. -Ns/2 が Rmin に対応
} BackParams;


void RadonTrans( ForwParams p );
double LineIntegral( double s, double theta, Field2D F, 
		     double dh, double Lim,
		     double (*rfunc)(void) );

double Delta( void );
double Poisson( void );
double Normal( void );

void ReconstructImage( BackParams p );
double ReconstuctIntegral( double x, double y, BackParams p );


// interface for R
void RadonInterface( double out[], int *nth, int *ns, double *ds,
		     double in[], int *ny, int *nx,
		     double *l, double *r, double *t, double *b,
		     double *pmu, double *psd, int *nh,
		     int *pflag, int *periodic );

void ReconstructInterface( double out[], int *ny, int *nx,
			   double *l, double *r, double *t, double *b,
			   double in[], int *nth, int *ns, double *ds );

#endif /* __RADONTRANS_H__ */
