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

#include <stdlib.h>
#include "configuration.h"
#include "tclap/CmdLine.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "crlDocumentedToolOutput.h"
#include "crlCCSegmentFilter.h"

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
	std::string fileTensor, fileCC; 

	/*-------------------------------------------------
	 First parse the command line (with TCLAP)
	-------------------------------------------------*/
	try {
		TCLAP::CmdLine cmd("", ' ', CRKIT_VERSION_STRING );
		crl::DocumentedToolOutput crlToolOutput(
			"Automatically compute the CC segmentation from a tensor file.",
			"Benoit Scherrer", 
			"2011" );
		crlToolOutput.SetKeepOriginalDataType(false);
		cmd.setOutput( &crlToolOutput );

		TCLAP::ValueArg<std::string>	argOutputCC  ("o","output",	"Output CC segmentation filename",	true, "", "Output CC segmentation filename", cmd);
		TCLAP::ValueArg<std::string>	argInputTensor  ("i","input",	"Input tensor file",	true, "", "Input tensor filename", cmd);

		// Parse !
		cmd.parse(argc,argv);

		// Get parameters
		fileTensor = argInputTensor.getValue();
		fileCC = argOutputCC.getValue();
	}
	catch ( TCLAP::ArgException& e)
	{
		std::cerr << "Error: "<< e.error() << " for argument " << e.argId() << std::endl;
		exit(1);
	}
	/*-------------------------------------
	Beginning of the program!
	-------------------------------------*/
	try {
		//---------------------------------------
		// Load the tensors
		//---------------------------------------
		typedef itk::DiffusionTensor3D< float >			TensorPixelType;
		typedef itk::Image< TensorPixelType, 3 >		TensorImageType;
		typedef itk::Image< float, 3 >					FloatImageType;
		typedef itk::ImageFileReader<TensorImageType>	TensorImageReaderType;
		std::cout<<"- Load the tensor file <" << fileTensor.c_str() << ">..." << std::endl;
		TensorImageReaderType::Pointer tensorReader = TensorImageReaderType::New();
		tensorReader->SetFileName(fileTensor);
		tensorReader->Update();

		//---------------------------------------
		// Automatic CC segmentation
		//---------------------------------------
		std::cout<<"- Compute the CC..." << std::endl;
		typedef itk::Image<unsigned short, 3> LabelImageType;
		crl::CCSegmentFilter::Pointer segmFilter = crl::CCSegmentFilter::New();
		segmFilter->SetInput(tensorReader->GetOutput());
		segmFilter->Update();

		//---------------------------------------
		// Save the result
		//---------------------------------------
		std::cout<<"- Write the CC..."<<std::endl;
		typedef itk::ImageFileWriter<LabelImageType>	LabelImageWriterType;
		LabelImageWriterType::Pointer labelWriter = LabelImageWriterType::New();
		labelWriter->SetFileName( fileCC );
		labelWriter->SetInput(segmFilter->GetOutput());
		labelWriter->Update();

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
  	return 0;
}
