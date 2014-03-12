/* Copyright 1997-2007 The MathWorks, Inc. */

/*
 * MRANKP.C
 * "Posix" C main program illustrating the use of the MATLAB Math Library.
 * Calls mlfMrank, obtained by using MCC to compile mrank.m.
 *
 * $Revision: 1.3.6.9 $
 *
 */

#include <stdio.h>
#include <math.h>
#include "libPkg.h"
#include "main_for_lib.h"



int run_main(int ac, char **av)
{
    mxArray *N;    /* Matrix containing n. */
    mxArray *R = NULL;    /* Result matrix. */
    int      n;    /* Integer parameter from command line. */

    /* Get any command line parameter. */
    
    if (ac >= 2) {
        n = atoi(av[1]);
    } else {
        n = 12;
    }

    /* Call the mclInitializeApplication routine. Make sure that the application
     * was initialized properly by checking the return status. This initialization
     * has to be done before calling any MATLAB API's or MATLAB Compiler generated
     * shared library functions. */
    if( !mclInitializeApplication(NULL,0) )
    {
        fprintf(stderr, "Could not initialize the application.\n");
        return -2;
    }
    /* Call the library intialization routine and make sure that the
     * library was initialized properly */
    if (!libPkgInitialize())
    {
      fprintf(stderr,"Could not initialize the library.\n");
      return -3;
    }
    else
    {
	/* Create a 1-by-1 matrix containing n. */
        N = mxCreateDoubleScalar(n);
      
	/* Call mlfMrank, the compiled version of mrank.m. */
	mlfMrank(1, &R, N);
	
	/* Print the results. */
	mlfPrintmatrix(R);
	
	/* Free the matrices allocated during this computation. */
	mxDestroyArray(N);
	mxDestroyArray(R);
	
	libPkgTerminate();    /* Terminate the library of M-functions */
    }
/* Note that you should call mclTerminate application in the end of
 * your application. mclTerminateApplication terminates the entire 
 * application.
 */
    mclTerminateApplication();
    return 0;
}
