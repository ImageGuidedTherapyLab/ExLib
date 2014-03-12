/* Copyright 2005-2007 The MathWorks, Inc. */

// catcherror.cpp
// MATLAB Compiler Example.
//
// This example demonstrates how to use print and error handling functions
// to control the output from a MATLAB Compiler-generated shared library.
// 
// This example relies on the M-file functions realacos and reveal, found in
// realacos.m and reveal.m, respectively. It calls realacos twice, once in a 
// way that should not fail, and once in a way that should.
//
// To build this example:
//     mcc -W cpplib:libracos -T link:lib realacos.m reveal.m
//     (Windows) mbuild catcherror.cpp libracos.lib
//     (Unix) mbuild catcherror.cpp libracos.so (or .dylib on the Mac) 
//

#include <stdio.h>

// The MCC command above will generate an include file, libracos.h, which
// contains a complete definition of the public interface of the libracos
// shared library. Include this Compiler generated header file first.

#include "libracos.h"

// Print and error handling functions. The print handler will be called for
// all non-error-related output, and the error handler will be called to 
// process error message strings.
// 
// The print handler wraps a banner around the printed output, like so:
//
// ******** MATLAB output *******
// *
// * <Output message>
// *
// ******** MATLAB output *******
//
// The error handler catches the most recent error into a static buffer, 
// for processing after the function returns, much like MATLAB's LASTERROR
// function.

static bool emitPrefix = false;

static int PrintHandler(const char *s)
{
    // Declare and initialize all variables 
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

// The main function. 

int run_main(int argc, char **argv)
{
    // Declare everything here, since C does not allow the mixing of 
    // declarations and code, as C++ does.

    double in1val = 0.1625;
    double in2val = 17;
        
    // Call the mclInitializeApplication routine. This initializes the MCR;
    // if this fails, the application cannot run. Make sure that the 
    // application was initialized properly by checking the return status. 
    // This initialization has to be done before calling any MATLAB API's or 
    // MATLAB Compiler generated shared library functions.

    if( !mclInitializeApplication(NULL,0) )
    {
        fprintf(stderr, "Could not initialize the application.\n");
        return -1;
    }
    
    // Declare input variables. Note that any attempt to manipulate
    // mwArrays before calling mclInitializeApplication will create a seg.
    // fault (the function pointers in the MCLMCRRT library will not be
    // initialized).
    mwArray in1(in1val);
    mwArray in2(in2val);

    // Declare a variable to hold the output
    mwArray out;

    // Initialize the MATLAB Compiler generated shared library. Pass in
    // pointers to our print and error handling functions.

    if (!libracosInitializeWithHandlers(ErrorHandler, PrintHandler)){
        fprintf(stderr,"Could not initialize the library.\n");
        return -2;
    }
    else
    {
        // Call the library function. This call should succeed.
	// Note: all MATLAB Compiler generated functions return true when
	// they succeed (complete without error).
	try
	{
	    realacos(1, out, in1);
	    std::cout << "realacos(" << in1val << ") = " << std::endl;
	    StartMatlabOutput();
	    reveal(out);
	    EndMatlabOutput();
	}
	catch(const mwException& e)
	{
	    std::cout << "Disaster! An unexpected error occurred."
		      << std::endl;
	    std::cout << "Error:" << std::endl << LastError << std::endl;
	    std::cerr << e.what() << std::endl;
	}

        // Call the library function. This call should fail, because 
	// the input value is greater than 1.0.
	//
	// Note: All MATLAB Compiler generated functions return false when
	// they fail.

	try
	{
	    realacos(1, out, in2);
	    // Display the return value of the library function.
	    // We expect this code will never execute!

	    std::cout << "Error! realacos should have returned false."
		      << std::endl;
	    std::cout << "realacos(" << in2val << ") = " << std::endl;
	    StartMatlabOutput();
	    reveal(out);
	    EndMatlabOutput();
	}
	catch(mwException &e)
	{
	    static char *errBanner = "\n######## MATLAB Error ########\n";
	    std::cerr << errBanner << LastError << errBanner;
	    std::cerr << "Exception: " << std::endl << e.what() << std::endl;
	}
        catch (...)
        {
          std::cerr << "Unexpected error thrown" << std::endl;
          return -3;
        }     

	/* Call the library termination routine */
	libracosTerminate();
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
