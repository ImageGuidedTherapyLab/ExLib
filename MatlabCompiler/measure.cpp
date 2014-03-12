/* Copyright 1999-2007 The MathWorks, Inc. */

/*
 *Doc example  Chapter 5.
 *$Revision: 1.1.6.2 $
 */

#include "collect_one_external.h"
#include <math.h>

extern double measure_from_device(void);

#ifdef __cplusplus
extern "C" {
#endif

bool collect_one(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
    printf("Hello from c");
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(plhs[0])) = measure_from_device();
    return 0;
}
#ifdef __cplusplus
}
#endif
double measure_from_device(void)
{
    static double t = 0.0;
    t = t + 0.05;
    return sin(t);
}
