#include "mex.h"

/* The computational routine */
void dist2ab(double *x,double *y, double *z,double *s2,
        double *Z2, double *S2,double *Lag,  mwSize n)
{
    mwSize i,j,l;
    l = 0;
    /* distance calculation */
    for (i=0; i<n; i++) {
       for (j=0; j<n; j++) {
        Z2[l]=(z[i]-z[j])*(z[i]-z[j]); 
        S2[l]=s2[i]+s2[j];
        Lag[l]=sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])); 
        l++;
       }
    }
}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double *x, *y, *z, *s2 ;
    double *Z2, *S2, *Lag; /* pointers to input & output matrices*/
    size_t n;      /* matrix dimensions */

    x = mxGetPr(prhs[0]); /* first input matrix */
    y = mxGetPr(prhs[1]); /* first input matrix */
    z = mxGetPr(prhs[2]); /* first input matrix */
    s2 = mxGetPr(prhs[3]); /* first input matrix */
    
    /* dimensions of input matrices */
    n = mxGetM(prhs[0]);;  
 
     /* create output matrix C */
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    Z2 = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL);
    S2 = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(n, n, mxREAL);
    Lag = mxGetPr(plhs[2]);

    /* Pass arguments to Fortran by reference */
    dist2ab(x,y,z,s2,Z2,S2,Lag, n);
}
