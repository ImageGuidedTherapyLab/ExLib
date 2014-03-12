/* Copyright 1999-2006 The MathWorks, Inc. */

/*
 *Doc example  Chapter 5.
 *$Revision: 1.1.6.4 $
 */

#include "collect_one_external.h"
#include <math.h>

extern double measure_from_device(void);

bool collect_one(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(plhs[0])) = measure_from_device();
}

double measure_from_device(void)
{
    static double t = 0.0;
    t = t + 0.05;
    return sin(t);
}

