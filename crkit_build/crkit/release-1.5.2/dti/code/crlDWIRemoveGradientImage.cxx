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

using namespace std;

//===============================================================
// Generalized Instantiator stuff
//===============================================================

// I. Include used files
#include "crlDWIRemoveGradientImage.h"

// II. Define some constants for crlGeneralizedInstantiator.h 
#undef		CRL_GA_SCALAR_CLASS			
#define		CRL_GA_VECTORIMAGE_CLASS	crl::DWI::RemoveGradientImage
#undef		CRL_GA_VECTOR_CLASS
#undef		CRL_GA_DIFFTENSOR_CLASS		
typedef		crl::DWI::RemoveGradientImage_Base AlgoBaseType;

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
	string fileInput, fileOutput, removeList;
	double bMax;

	/*-------------------------------------------------
	 First parse the command line (with TCLAP)
	-------------------------------------------------*/
	try {
		TCLAP::CmdLine cmd("", ' ', CRKIT_VERSION_STRING );
		crl::DocumentedToolOutput crlToolOutput(
			"CRL Diffusion Weighted Imaging ToolBox.\n"
			"Removes a set of gradient images in a NHDR RAW file." ,
			"Benoit Scherrer", 
			"2010" );
		crlToolOutput.SetKeepOriginalDataType(true);
		cmd.setOutput( &crlToolOutput );

		TCLAP::ValueArg<double>			argBMax	("","bmax",		"Keep only the images with b-value inferior or egal to the given maximum b-value.",	false, -1, "value", cmd);
		TCLAP::ValueArg<string>			argRmList	("r","remove",		"Comma-separated list of gradient images indexes to remove (zero-based indexes).",	false, "", "[i0],[i1],[i2] ...", cmd);
		TCLAP::ValueArg<string>			argOutput	("o","output",		"Output NHDR file",	true, "", "Output file", cmd);
		TCLAP::ValueArg<string>			argInput	("i","input",		"Input NHDR file",	true, "", "Input file", cmd);
	
		// Parse !
		cmd.parse(argc,argv);

		// Get parameters
		fileInput	= argInput.getValue();
		fileOutput	= argOutput.getValue();
		removeList		= argRmList.getValue();
		bMax = argBMax.getValue();

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
		convertAlgo->SetBMax(bMax);

		convertAlgo->ClearListOfIndexes();

		if ( removeList!="" )
		{
			char szBuffer[256];
			int Index, nbIndex=0;
			const char *p = removeList.c_str();
			const char *p_end = p + removeList.length();
			while (p<p_end)  
			{
				if (sscanf(p, "%d", &Index )!=1)
					throw itk::ExceptionObject(__FILE__, __LINE__, "Invalid list of image indexes to remove.");

				convertAlgo->AddIndexToRemove(Index);
				nbIndex++;
				sprintf(szBuffer, "%d", Index);
				p += strlen(szBuffer) + 1;
			}

			if ( nbIndex==0 )
				throw itk::ExceptionObject(__FILE__, __LINE__, "Invalid list of image indexes to remove.");
		}

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
