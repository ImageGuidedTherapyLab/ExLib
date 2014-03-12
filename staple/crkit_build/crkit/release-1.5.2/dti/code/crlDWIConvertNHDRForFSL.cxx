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
#include "crlDocumentedToolOutput.h"
#include <iostream>
#include <fstream>

using namespace std;

//===============================================================
// Generalized Instantiator stuff
//===============================================================

// I. Include used files
#include "crlDWIConvertNHDRForFSL.h"

// II. Define some constants for crlGeneralizedInstantiator.h
#undef		CRL_GA_SCALAR_CLASS			
#define		CRL_GA_VECTORIMAGE_CLASS	crlDWIConvertNHDRForFSL_Algo
#undef		CRL_GA_VECTOR_CLASS
#undef		CRL_GA_DIFFTENSOR_CLASS		
typedef		crl::crlDWIConvertNHDRForFSL_Base AlgoBaseType;

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
 * \date	January 2010
 *
 * \exception	itk	Thrown when . 
 *
 * \param	argc	Number of command-line arguments. 
 * \param	argv	Array of command-line argument strings. 
 *
 * \return	Exit-code for the process - 0 for success, else an error code. 
*************************************************************************************************/
int main(int argc, char **argv)
{
	string fileInput;
	string fnBVecs, fnBVals, fnData, fnB0 ;
	float B0max;

	/*-------------------------------------------------
	 First parse the command line (with TCLAP)
	-------------------------------------------------*/
	try {
		TCLAP::CmdLine cmd("", ' ', CRKIT_VERSION_STRING );
		crl::DocumentedToolOutput crlToolOutput(
			"CRL Diffusion Weighted Imaging ToolBox.\n"
			"Converts a NHDR RAW data to the FSL RAW data format (4D data file\n"
			"+ b-vectors text file + b-values text file). \n"
			"WARNING: By default, all b-values smaller or equal to 10s/mm2 are forced to B=0s/mm2.",
			"Benoit Scherrer", 
			"2010" );
		crlToolOutput.SetKeepOriginalDataType(true);
		cmd.setOutput( &crlToolOutput );

		TCLAP::ValueArg<float>	argB0max ("","b0max", "Set the maximum b-value to be considered as a B=0 image (by default: 10)", false, 10, "value", cmd);
		TCLAP::ValueArg<string> argBVecs  ("","bvecs",	"Set the filename for b-vectors [default: bvecs]", false, "bvecs", "filename", cmd);
		TCLAP::ValueArg<string> argBVals  ("","bvals",	"Set the filename for b-values [default: bvals]", false, "bvals", "filename", cmd);
		TCLAP::ValueArg<string> argData   ("","data",	"Set the filename for the data [default: data.nii.gz]", false, "data.nii.gz", "filename", cmd);
		TCLAP::ValueArg<string> argB0     ("","b0",		"If set, extract a B=0 volume in that file", false, "", "filename", cmd);

		TCLAP::ValueArg<string> argInput     ("i","input",	"Input file",	true, "", "Input file", cmd);
	
		// Parse !
		cmd.parse(argc,argv);

		// Get parameters
		fileInput = argInput.getValue();
		fnBVecs = argBVecs.getValue();
		fnBVals = argBVals.getValue();
		fnData = argData.getValue();
		fnB0  = argB0.getValue();
		B0max = argB0max.getValue();
	}
	catch ( TCLAP::ArgException& e)
	{
		std::cerr << "Error: "<<e.error()<< " for argument " <<e.argId() << std::endl;
		exit(1);
	}

	/*-------------------------------------
	Beginning of the program!
	-------------------------------------*/
	try {
		// Instantiate the mask algo with MyInstantiator
		crl::Instantiator Instantiator;
		AlgoBaseType::Pointer convertAlgo = Instantiator.Instantiate(fileInput );
		
		// Set the parameters for the algo
		convertAlgo->SetInputFileName(fileInput);
		convertAlgo->SetOutputFileName(fnData);
		convertAlgo->SetOutputBVecsFileName(fnBVecs);
		convertAlgo->SetOutputBValsFileName(fnBVals);
		convertAlgo->SetOutputB0FileName(fnB0);
		convertAlgo->SetB0max(B0max);

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
	cout << endl;
  	return 0;
}
