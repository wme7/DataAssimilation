#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mex.h"

static void L40(double x, int n, double y[], double y1[])
{
    int i;

    for (i = 0; i < n; ++i) {
        int ip1 = (i + 1) % n;
        int im1 = (i + n - 1) % n;
        int im2 = (i + n - 2) % n;

        y1[i] = (y[ip1] - y[im2]) * y[im1] - y[i] + 8.0;
    }
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

    L40(x, n, y, k1);
    for (i = 0; i < n; ++i)
	yy[i] = y[i] + h2 * k1[i];
    x += 0.5 * h;
    L40(x, n, yy, k2);
    for (i = 0; i < n; ++i)
	yy[i] = y[i] + h2 * k2[i];
    L40(x, n, yy, k3);
    for (i = 0; i < n; ++i)
	yy[i] = y[i] + h * k3[i];
    x += 0.5 * h;
    L40(x, n, yy, k4);

    for (i = 0; i < n; ++i)
	yout[i] = y[i] + h * (k1[i] + 2 * (k2[i] + k3[i]) + k4[i]) / 6.0;

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
	mexErrMsgTxt("L40_step_c: wrong number of input or output arguments");

    x = mxGetPr(prhs[0]);
    N = mxGetN(prhs[0]);
    M = mxGetM(prhs[0]);
    n = N * M;

    plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    xout = mxGetPr(plhs[0]);
    for (i = 0; i < n; ++i)
	xout[i] = x[i];

    rk4step(0.0, 0.05, n, xout, xout);
}
