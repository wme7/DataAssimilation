#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mex.h"

static void L3(double y[], double y1[])
{
    y1[0] = 10.0 * (y[1] - y[0]);
    y1[1] = (28.0 - y[2]) * y[0] - y[1];
    y1[2] = y[0] * y[1] - (8.0 / 3.0) * y[2];
}

static void rk4step(double x, double xend, int n, double y[], double yout[])
{
    double* k = malloc(n * 5 * sizeof(double));
    double* k1 = k;
    double* k2 = &k[n];
    double* k3 = &k[n * 2];
    double* k4 = &k[n * 3];
    double* yy = &k[n * 4];
    double h = xend - x;
    double h2 = h / 2.0;
    int i;

    L3(y, k1);
    for (i = 0; i < 3; ++i)
	yy[i] = y[i] + h2 * k1[i];
    L3(yy, k2);
    for (i = 0; i < 3; ++i)
	yy[i] = y[i] + h2 * k2[i];
    L3(yy, k3);
    for (i = 0; i < 3; ++i)
	yy[i] = y[i] + h * k3[i];
    L3(yy, k4);

    for (i = 0; i < 3; ++i)
	yout[i] = y[i] + h * (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) / 6.0;

    free(k);
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    double* dt;
    double* x;
    double* xout;
    int N, M, n;
    int i;

    if (nlhs != 1 || nrhs != 1)
	mexErrMsgTxt("L3_step_c: wrong number of input or output arguments");

    x = mxGetPr(prhs[0]);
    N = mxGetN(prhs[0]);
    M = mxGetM(prhs[0]);
    n = N * M;

    plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    xout = mxGetPr(plhs[0]);
    for (i = 0; i < n; ++i)
	xout[i] = x[i];

    rk4step(0.0, 0.01, n, xout, xout);
}
