/* Copyright 2005-2007 The MathWorks, Inc. */

/* catcherror.c
 * MATLAB Compiler Example.
 *
 * This example demonstrates how to use print and error handling functions
 * to control the output from a MATLAB Compiler-generated shared library.
 * 
 * This example relies on the M-file functions realacos and reveal, found in
 * realacos.m and reveal.m, respectively. It calls realacos twice, once in a 
 * way that should not fail, and once in a way that should.
 *
 * To build this example:
 *     mcc -W lib:libracos -T link:lib realacos.m reveal.m
 *     (Windows) mbuild catcherror.c libracos.lib
 *     (Unix) mbuild catcherror.c libracos.so (or .dylib on the Mac) 
 */

#include <stdio.h>

/* The MCC command above will generate an include file, libracos.h, which
 * contains a complete definition of the public interface of the libracos
 * shared library. Include this Compiler generated header file first.
 */
#include "libracos.h"

/* Print and error handling functions. The print handler will be called for
 * all non-error-related output, and the error handler will be called to 
 * process error message strings.
 * 
 * The print handler wraps a banner around the printed output, like so:
 *
 * ******** MATLAB output *******
 * *
 * * <Output message>
 * *
 * ******** MATLAB output *******
 *
 * The error handler catches the most recent error into a static buffer, 
 * for processing after the function returns, much like MATLAB's LASTERROR
 * function.
 */

static bool emitPrefix = false;

static int PrintHandler(const char *s)
{
    /* Declare and initialize all variables */
    int len = 0;
    static char * prefix = "* ";
    int written = 0;

    if (s == NULL) return 0;

    len = strlen(s);

    /* DISP adds two carriage returns. Suppress the last one. */
    if (len >= 2 && s[len-1] == '\n' && s[len-2] == '\n')
	len = len-1;

    if (emitPrefix)
	fwrite(prefix, sizeof(char), strlen(prefix), stdout);

    written = fwrite(s, sizeof(char), len, stdout);

    if (s[len-1] == '\n')
	emitPrefix = true;
    else
	emitPrefix = false;

    return written;
}

static void StartMatlabOutput()
{
    static char *startBanner = "******** Start MATLAB output *******\n";
    emitPrefix = false;
    PrintHandler(startBanner);
}

static void EndMatlabOutput()
{
    static char *endBanner = "********* End MATLAB output ********\n";
    emitPrefix = false;
    PrintHandler(endBanner);
}

static char LastError[2048];

static int ErrorHandler(const char *s)
{
    int len = 0;
    len = strlen(s);
    LastError[0] = '\0';
    strcpy(LastError, s);
    return len+1;
}

/* The main function. 
 */

int run_main(int argc, char **argv)
{
    /* Declare everything here, since C does not allow the mixing of 
     * declarations and code, as C++ does.
     */
    double in1val = 0.1625;
    double in2val = 17;
    mxArray *in1 = NULL;
    mxArray *in2 = NULL;
    
    /* Declare a variable to hold the output */
    mxArray *out = NULL;

    /* Call the mclInitializeApplication routine. This initializes the MCR;
     * if this fails, the application cannot run. Make sure that the 
     * application was initialized properly by checking the return status. 
     * This initialization has to be done before calling any MATLAB API's or 
     * MATLAB Compiler generated shared library functions.
     */

    if( !mclInitializeApplication(NULL,0) )
    {
        fprintf(stderr, "Could not initialize the application.\n");
	    return -1;
    }
    
    /* Create the input data */
    in1 = mxCreateDoubleScalar(in1val);
    in2 = mxCreateDoubleScalar(in2val);

    /* Initialize the MATLAB Compiler generated shared library. Pass in
     * pointers to our print and error handling functions.
     */
    if (!libracosInitializeWithHandlers(ErrorHandler, PrintHandler)){
        fprintf(stderr,"Could not initialize the library.\n");
        return -2;
    }
    else
    {
        /* Call the library function. This call should succeed.
	 * Note: all MATLAB Compiler generated functions return true when
	 * they succeed (complete without error).
	 */
        if (mlfRealacos(1, &out, in1))
	{
	    /* Display the return value of the library function */
	    printf("realacos(%6.4f) = \n", in1val);
	    StartMatlabOutput();
	    mlfReveal(out);
	    EndMatlabOutput();
	}
	else
	{
	    printf("Disaster! An unexpected error occurred.\n");
	    printf("Error:\n%s\n", LastError);
	}

	/* Destroy the return value since this varaible will be resued in
	 * the next function call. Since we are going to reuse the variable,
	 * we have to set it to NULL. Refer to MATLAB Compiler documentation
	 * for more information on this.
	 */
	mxDestroyArray(out); out=0;

        /* Call the library function. This call should fail, because 
	 * the input value is greater than 1.0.
	 *
	 * Note: All MATLAB Compiler generated functions return false when
	 * they fail.
	 */
        if (mlfRealacos(1, &out, in2))
	{
	    /* Display the return value of the library function.
	     * We expect this code will never execute!
	     */
	    printf("Error! mlfRealacos should have returned false.\n");
	    printf("realacos(%6.4f) = \n", in2val);
	    StartMatlabOutput();
	    mlfReveal(out);
	    EndMatlabOutput();
	}
	else
	{
	    static char *errBanner = "\n######## MATLAB Error ########\n";
	    printf("%s %s %s", errBanner, LastError, errBanner);
	}

	/* Call the library termination routine */
	libracosTerminate();
	
	/* Free the memory we allocated for the input variables. */
	mxDestroyArray(in1); in1=0;
	mxDestroyArray(in2); in2=0;
    }

    /* mclTerminateApplication simply shuts down the MCR. 
     * (You cannot restart it by calling mclInitializeApplication.
     * Call mclTerminateApplication once and only once in your application.)
     */

    mclTerminateApplication();
    return 0;
}

/* The main routine. On the Mac, the system code runs in the main thread, and
 * user code must be processed by a secondary thread.
 *
 * On other platforms, the main thread runs the user code.
 */
int main()
{
    mclmcrInitialize();
    return mclRunMain((mclMainFcnType)run_main,0,NULL);
}
