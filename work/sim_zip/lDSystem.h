#include <stdio.h> 
//#ifndef __lDSysstem__
//#define __lDSysstem__


void DSReadLn (FILE *stream, char str[], int);

void DSReadGeometryStr (FILE *stream, int *num, double *gx, double *gy, double *gz);

void DSReadFieldStr (FILE *stream, int *num, double *gx, double *gy, double *gz, double *gsum);

int DSReadAnsysData(void);

void DSfield(double U, double rx, double ry, double rz, double *Ex, double *Ey, double *Ez);

//#endif
