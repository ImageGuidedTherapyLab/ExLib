/*
 * Copyright (c) 2008-2009 Children's Hospital Boston.
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

//===============================================================
// Generalized Instantiator stuff
//===============================================================

// I. Include used files
#include "crlMFMTensorResampleClass.h"

// II. Define some constants for crlGeneralizedInstantiator.h
#undef		CRL_GA_SCALAR_CLASS
#undef		CRL_GA_VECTORIMAGE_CLASS	
#undef		CRL_GA_VECTOR_CLASS
#define		CRL_GA_DIFFTENSOR_CLASS		crlMFMTensorResample
typedef		crl::crlMFMTensorResample_Base AlgoBaseType;

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
 * \param	argc	Number of command-line arguments. 
 * \param	argv	Array of command-line argument strings. 
 *
 * \return	Exit-code for the process - 0 for success, else an error code. 
*************************************************************************************************/
int main(int argc, char **argv)
{
	std::vector<std::string> inputFiles;
	std::string voxelSize, interpMode, inputTransform, outputBaseFileName, inputFixedImage,  fractionsFileName;

	/*-------------------------------------------------
	 First parse the command line (with TCLAP)
	-------------------------------------------------*/
	try {
		TCLAP::CmdLine cmd("", ' ', CRKIT_VERSION_STRING );
		crl::DocumentedToolOutput crlToolOutput(
			"Resample tensors. Support multi-tensors.",
			"Benoit Scherrer", 
			"2010" );
		crlToolOutput.SetKeepOriginalDataType(true);		
		cmd.setOutput( &crlToolOutput );

		//-------------------------------------------
		// Set up TCLAP
		//-------------------------------------------
		std::vector<std::string> strInterpModes;
		strInterpModes.push_back("none");
		strInterpModes.push_back("euclidean");
		strInterpModes.push_back("logeuclidean");
		TCLAP::ValuesConstraint<std::string> interpModesArg(strInterpModes);

		TCLAP::ValueArg<std::string> argInterpolation("","interp","Tensor interpolation mode : none, euclidean interpolation or log-euclidean interpolation (default: euclidean)",false,"euclidean", &interpModesArg, cmd);
	    TCLAP::ValueArg<std::string> argFractionsFile("", "fractions", "Input 4D Fraction of Occupancy File Name", false, "", "input 4D file name", cmd);

		TCLAP::ValueArg<std::string> argTransformFile ("t","transform", "Input transform file (if not defined, use identity)", false, "", "Transform filename", cmd);

		TCLAP::ValueArg<std::string> argVoxelSize ("","voxelsize", "Convert the voxel size to a new size. Example: --voxelsize 0.8,0.8,1 ", false, "", "x,y,z");
		TCLAP::ValueArg<std::string> argFixedImage("g","geometry","Set the fixed image to get the geometry from (size, spacing, orientation, ...)", false, "", "Geometry image filename");
		std::vector<TCLAP::Arg *> modes;
		modes.push_back(&argFixedImage);
		modes.push_back(&argVoxelSize);
		cmd.xorAdd(modes);

		TCLAP::ValueArg<std::string> argOutputBase("o","output","Output base file name. If several tensors, the created files will be [base]_t0, [base]_t1, ...", true, "", "base filename", cmd);
		TCLAP::MultiArg<std::string> argInputs("i", "inputFile", "Input (moving) tensor file name", true, "input tensor file name", cmd);

		//-------------------------------------------
		// Parse !
		//-------------------------------------------
		cmd.parse(argc,argv);

		//-------------------------------------------
		// Get parameters
		//-------------------------------------------
		inputFiles = argInputs.getValue();
		outputBaseFileName = argOutputBase.getValue();
		interpMode = argInterpolation.getValue();
		inputTransform = argTransformFile.getValue();
		inputFixedImage = argFixedImage.getValue();
		voxelSize = argVoxelSize.getValue();
		fractionsFileName = argFractionsFile.getValue();

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
		//-------------------------------------------
		// Instantiate the algo : takes into account the input tensor type 
		//-------------------------------------------
		crl::Instantiator Instantiator;
		AlgoBaseType::Pointer resampleAlgo = Instantiator.Instantiate(inputFiles[0]);
		if ( resampleAlgo.GetPointer()==NULL )
			throw itk::ExceptionObject(__FILE__,__LINE__,"The input file type was not recognized.", "");

		//-------------------------------------------
		// Set the parameters for the algo
		//-------------------------------------------
		resampleAlgo->SetInputFileNames(inputFiles);
		resampleAlgo->SetOutputBaseFileName(outputBaseFileName);
		resampleAlgo->SetTransformFile(inputTransform);
		resampleAlgo->SetGeometryFixedImage(inputFixedImage);
		resampleAlgo->SetVoxelSizeMode(voxelSize);
		resampleAlgo->SetInputFractionsFileName(fractionsFileName);

		if ( interpMode=="euclidean" ) resampleAlgo->SetTensorInterpolationMode(AlgoBaseType::INTERP_EUCLIDEAN);
		else if ( interpMode=="logeuclidean" ) resampleAlgo->SetTensorInterpolationMode(AlgoBaseType::INTERP_LOGEUCLIDEAN);
		else resampleAlgo->SetTensorInterpolationMode(AlgoBaseType::INTERP_NONE);

		//-------------------------------------------
		// Execute!
		//-------------------------------------------
		resampleAlgo->Execute();		
		std::cout << " Done" << std::endl << std::endl;
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
	std::cout << std::endl;
  	return 0;
}
