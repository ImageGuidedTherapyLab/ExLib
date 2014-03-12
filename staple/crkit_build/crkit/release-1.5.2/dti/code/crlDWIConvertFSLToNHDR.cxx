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
#include <fstream>

#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include "crlFileName.h"
#include "crlDocumentedToolOutput.h"
#include "crlVectorImageUtils.h"

using namespace std;

//===============================================================
// Generalized Instantiator stuff
//===============================================================

// I. Include used files
#include "crlDWIConvertFSLToNHDR.h"

// II. Define some constants for crlGeneralizedInstantiator.h
#define		CRL_GA_SCALAR_CLASS			crlDWIConvertFSLToNHDR_Algo
#undef		CRL_GA_VECTORIMAGE_CLASS	
#undef		CRL_GA_VECTOR_CLASS
#undef		CRL_GA_DIFFTENSOR_CLASS		
typedef		crl::crlDWIConvertFSLToNHDR_Base AlgoBaseType;

// ONLY 4-D images as input
#define MIN_IMAGE_DIM 4
#define MAX_IMAGE_DIM 4

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
	string fileInput, fileOutput;
	string fnBVecs, fnBVals, fnRef ;
	float B0max;

	/*-------------------------------------------------
	 First parse the command line (with TCLAP)
	-------------------------------------------------*/
	try {
		TCLAP::CmdLine cmd("IMPORTANT REMARK: This tool keeps the input data type (unless --ref is specified). \n (c) Computational Radiology Laboratory\n       (Benoit Scherrer)", ' ', CRKIT_VERSION_STRING );
		crl::DocumentedToolOutput crlToolOutput(
			"CRL Diffusion Weighted Imaging ToolBox.\n"
			"Converts the FSL RAW data (4D data file + b-vectors text file + "
			"b-values text file) to NHDR RAW data.\n"
			"WARNING: By default, all b-values smaller or equal to 10s/mm2 are forced to B=0s/mm2.",
			"Benoit Scherrer", 
			"2010" );
		crlToolOutput.SetKeepOriginalDataType(true);
		cmd.setOutput( &crlToolOutput );

		TCLAP::ValueArg<float>	argB0max ("","b0max", "Set the maximum b-value to be considered as a B=0 image (by default: 10)", false, 10, "value", cmd);
		TCLAP::ValueArg<string> argRef  ("","ref",	"Set a reference file for output data type. FSL indeed do not keep the input data type and convert all files in int (s32)", false, "", "filename", cmd);
		TCLAP::ValueArg<string> argBVecs  ("","bvecs",	"Set the input filename for b-vectors [default: bvecs]", false, "", "filename", cmd);
		TCLAP::ValueArg<string> argBVals  ("","bvals",	"Set the input filename for b-values [default: bvals]", false, "", "filename", cmd);
		TCLAP::ValueArg<string> argOutput ("o","output", "Output NHDR file", true, "", "Output filename", cmd);
		TCLAP::ValueArg<string> argInput  ("i","input",	"Input 4D FSL file",	true, "", "Input filename", cmd);
		
		// Parse !
		cmd.parse(argc,argv);

		// Get parameters
		fileInput = argInput.getValue();
		fileOutput = argOutput.getValue();
		fnBVecs = argBVecs.getValue();
		fnBVals = argBVals.getValue();
		fnRef = argRef.getValue();
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
		crl::FileName fnOutput(fileOutput);
		crl::FileName fnInput(fileInput);

		if ( fnOutput.getExtension()!="nhdr" )
			throw itk::ExceptionObject(__FILE__,__LINE__,"Error: the output file must have the extension NHDR", "");

		// default value: bvecs and bvals, in the same directory as the input
		if (fnBVecs=="") fnBVecs = fnInput.getPath()+"bvecs";
		if (fnBVals=="") fnBVals = fnInput.getPath()+"bvals";

		// Instantiate the algo 
		crl::Instantiator Instantiator;

		// Take into account a reference data type
		crl::Instantiator::ImageTypeInfos infos = Instantiator.GetImageType( fileInput );
		if ( fnRef!="" )
		{
			try {
				crl::Instantiator::ImageTypeInfos infos2 = Instantiator.GetImageType( fnRef );
				infos.componentType = infos2.componentType;
			}
			catch (...) {}
		}

		AlgoBaseType::Pointer convertAlgo = Instantiator.Instantiate(infos);
		if ( convertAlgo.GetPointer()==NULL )
			throw itk::ExceptionObject(__FILE__,__LINE__,"The input file type was not recognized.", "");

		// Set the parameters for the algo
		convertAlgo->SetInputDataFileName(fileInput);
		convertAlgo->SetInputBVecsFileName(fnBVecs);
		convertAlgo->SetInputBValsFileName(fnBVals);
		convertAlgo->SetOutputFileName(fileOutput);
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
