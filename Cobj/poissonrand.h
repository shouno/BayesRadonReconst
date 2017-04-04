#ifndef POISSONRAND_H_
#define POISSONRAND_H_

#include <math.h>
#include "mtrand.h"
 
#ifndef PI
#define PI 3.141592653589793
#endif

#ifdef __cplusplus
extern "C" {
#endif

double gammln( double xx );
double PoissonRand( double xm );

#ifdef __cplusplus
}
#endif

#endif /*POISSONRAND_H_*/
