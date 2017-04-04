#ifndef __FIELD2D_H_
#define __FIELD2D_H_ 1
typedef struct Field2D {
    double left;
    double right;
    double top;
    double bottom;

    double *xpos;
    double *ypos;
    double *ary;
    int nrow;
    int ncol;

    int isPeriodic;
    double xT;
    double yT;
} Field2D;

void SetupField( Field2D *f, double l, double r, double t, double b,
		 double data[], int nrow, int ncol,
		 int periodic );
void DestroyField( Field2D f );
double FieldValue( double x, double y, Field2D f );
int inArea( double x, double y, Field2D f );

#endif /* __FIELD2D_H_ */
