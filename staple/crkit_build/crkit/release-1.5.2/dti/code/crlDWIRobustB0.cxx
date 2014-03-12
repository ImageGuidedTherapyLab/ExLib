/*
 * Copyright (c) 2008-2011 Children's Hospital Boston.
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 * http://www.crl.med.harvard.edu
 * For more information contact: simon.warfield@childrens.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#include "configuration.h"
#include "tclap/CmdLine.h"
#include <iostream>

using namespace std;

//===============================================================
// Generalized Instantiator stuff
//===============================================================

// I. Include used files
#include "crlDWIRobustB0.h"

// II. Define some constants for crlGeneralizedInstantiator.h
#undef		CRL_GA_SCALAR_CLASS			
#define		CRL_GA_VECTORIMAGE_CLASS	crl::DWI::DWIRobustB0
#undef		CRL_GA_VECTOR_CLASS
#undef		CRL_GA_DIFFTENSOR_CLASS		
typedef		crl::DWI::DWIRobustB0_Base AlgoBaseType;

// ONLY 3-D images as input
#define MIN_IMAGE_DIM 3
#define MAX_IMAGE_DIM 3

// III. Now include the file crlGeneralizedInstantiator.h
#include	"crlGeneralizedInstantiator.h"


/**********************************************************************************************//**
 * \fn	int main(int argc, char **argv)
 *
 * \brief	Main entry-point for this application. 
 *
 * \author	Benoit Scherrer
 * \date	December 2009
 *
 * \param	argc	Number of command-line arguments. 
 * \param	argv	Array of command-line argument strings. 
 *
 * \return	Exit-code for the process - 0 for success, else an error code. 
*************************************************************************************************/

int main(int argc, char **argv)
{
	string fileInput, fileOutput, fileMean, fileDiff, fileNewRaw;
	int nbProcs;
	float b0threshold;

	/*-------------------------------------------------
	 First parse the command line (with TCLAP)
	-------------------------------------------------*/
	try {
		TCLAP::CmdLine cmd("IMPORTANT REMARK: This tool keeps the input data type \n (c) Computational Radiology Laboratory\n       (Benoit Scherrer)", ' ', CRKIT_VERSION_STRING );

		TCLAP::ValueArg<unsigned int>	argProc		("p","proc",		"Number of processors",	false, 1, "Nb Proc", cmd);
		TCLAP::ValueArg<float>			argB0Threshold ("","b0threshold",		"Define the maximum b-value to be considered a b=0 image (default: 10)",	false, 10, "b-value", cmd);
		TCLAP::ValueArg<string>			argNewRaw	("","newraw",		"If set, output a new NHDR DWI raw file the only one B=0 robust image.", false, "", "NHDR output file", cmd);
		TCLAP::ValueArg<string>			argDiff		("d","diff",		"If set, output an image of the difference between the mean and the robust mean",	false, "", "Diff image file", cmd);
		TCLAP::ValueArg<string>			argMean		("m","mean",		"If set, output the image of the non-robust mean",	false, "", "Mean image file", cmd);
		TCLAP::ValueArg<string>			argOutput	("o","output",		"Output file [default: [inputfile]_B0 ]",	false, "", "Output file", cmd);
		TCLAP::ValueArg<string>			argInput	("i","input",		"Input NHDR file",	true, "", "Input file", cmd);
	
		// Parse !
		cmd.parse(argc,argv);

		// Get parameters
		fileInput	= argInput.getValue();
		fileOutput	= argOutput.getValue();
		nbProcs		= argProc.getValue();
		fileMean	= argMean.getValue();
		fileDiff	= argDiff.getValue();
		fileNewRaw	= argNewRaw.getValue();
		b0threshold = argB0Threshold.getValue();
	}
	/*-------------------------------------------------
	 If parsing error
	-------------------------------------------------*/
	catch ( TCLAP::ArgException& e)
	{
		std::cerr << "Error: "<<e.error()<< " for argument " <<e.argId() << std::endl;
		exit(1);
	}


	/*-------------------------------------
	Beginning of the program!
	-------------------------------------*/
	try {

		// Instantiate the algo with the generalized Instantiator
		crl::Instantiator Instantiator;
		AlgoBaseType::Pointer convertAlgo = Instantiator.Instantiate(fileInput);
		
		// Set the parameters for the algo
		convertAlgo->SetInputFileName(fileInput);
		convertAlgo->SetOutputFileName(fileOutput);
		convertAlgo->SetNumberOfThreads(nbProcs);
		convertAlgo->SetMeanFileName(fileMean);
		convertAlgo->SetDiffFileName(fileDiff);
		convertAlgo->SetNewRAWFileName(fileNewRaw);
		convertAlgo->SetB0Threshold(b0threshold);

		// Execute!
		convertAlgo->Execute();
		cout << " Done" << endl << endl;


	}

	/*-------------------------------------
	If catched an exception, show a message
	-------------------------------------*/
	catch (itk::ExceptionObject& e)
	{
		std::cout<< "ITK EXCEPTION ERROR CAUGHT"<<std::endl<< e.GetDescription() << std::endl << "Cannot continue." << std::endl ;
		return 1;
	}
	catch ( ... )
	{
		std::cout<<"UNKNOWN EXCEPTION ERROR." << std::endl << "Cannot continue."<< std::endl;
		return 1;
	}

  return 0;
}
