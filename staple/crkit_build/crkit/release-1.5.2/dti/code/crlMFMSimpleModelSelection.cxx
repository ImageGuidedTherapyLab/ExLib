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

#include "itkLabelStatisticsImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

#include "crlDocumentedToolOutput.h"
#include "crlDWIStudy.h"
#include "crlCCSegmentFilter.h"
#include "crlDiffusionTensorMeasuresFilter.h"

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
	std::string fileTensor, fileOutputMask, segmCC, fileMSE; 

	float sigmaFactor;
	int nbProc;

	/*-------------------------------------------------
	 First parse the command line (with TCLAP)
	-------------------------------------------------*/
	try {
		TCLAP::CmdLine cmd("", ' ', CRKIT_VERSION_STRING );
		crl::DocumentedToolOutput crlToolOutput(
			"Compute a two-tensor mask from a DWI acquisition, based on the shape of the one-tensor solution.",
			"Benoit Scherrer", 
			"2011" );
		crlToolOutput.SetKeepOriginalDataType(false);
		cmd.setOutput( &crlToolOutput );

		TCLAP::ValueArg<unsigned int>	argProc			("p","proc",		"Set the number of processors",	false, 1, "Nb Proc", cmd);

		TCLAP::ValueArg<std::string>	argSegmCC("","cc","Automatic segmentation of the CC to determine the threshold",false,"", "CC Output Filename", cmd);

		TCLAP::ValueArg<float>			argSigmaFactor ("","sigmafactor", "Set the sigma factor for the threshold. threshold=mean+sigmafactor*sigma (default: 2.5)", false, 2.5, "float", cmd);

		TCLAP::ValueArg<std::string>	argOutputMSE  ("","mse",	"Output MSE",	false, "", "Output MSE", cmd);
		TCLAP::ValueArg<std::string>	argOutputMask  ("o","output",	"Output 2T mask",	true, "", "Output mask", cmd);
		TCLAP::ValueArg<std::string>	argInputTensor  ("t","tensor",	"Input tensor file",	true, "", "Input tensor", cmd);

		// Parse !
		cmd.parse(argc,argv);

		// Get parameters
		fileTensor = argInputTensor.getValue();
		fileOutputMask = argOutputMask.getValue();
		nbProc = argProc.getValue();
		segmCC=argSegmCC.getValue();
		fileMSE = argOutputMSE.getValue();
		sigmaFactor = argSigmaFactor.getValue();
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
		std::cout<<"- Read "<<fileTensor.c_str()<<"..."<<std::endl;
		typedef itk::DiffusionTensor3D< float >			TensorPixelType;
		typedef itk::Image< TensorPixelType, 3 >		TensorImageType;
		typedef itk::Image< float, 3 >					FloatImageType;
		typedef itk::ImageFileReader<TensorImageType>	TensorImageReaderType;
		TensorImageReaderType::Pointer tensorReader = TensorImageReaderType::New();
		tensorReader->SetFileName(fileTensor);
		tensorReader->Update();

		//---------------------------------------
		// Automatic CC segmentation
		//---------------------------------------
		typedef itk::Image<unsigned short, 3> LabelImageType;
		LabelImageType::Pointer ccMask;
		if ( segmCC!="" ) {
			crl::CCSegmentFilter::Pointer segmFilter = crl::CCSegmentFilter::New();
			segmFilter->SetInput(tensorReader->GetOutput());
			segmFilter->Update();

			std::cout<<"- Write the CC..."<<std::endl;
			typedef itk::ImageFileWriter<LabelImageType>	LabelImageWriterType;
			LabelImageWriterType::Pointer floatWriter = LabelImageWriterType::New();
			floatWriter->SetFileName(segmCC );
			floatWriter->SetInput(segmFilter->GetOutput());
			floatWriter->Update();

			ccMask = segmFilter->GetOutput();
		}

	
		std::cout<<"- Simple MOSE mode: estimate the linear component of each tensor."<<std::endl;
		crl::DiffusionTensorMeasuresFilter::Pointer Filter = crl::DiffusionTensorMeasuresFilter::New();
		Filter->SetInput(tensorReader->GetOutput());
		Filter->SetNumberOfThreads(nbProc);
		Filter->Update();
		
		FloatImageType::Pointer measureImage = Filter->GetOutput();

		//------------------------------------------
		// Write the MSE if needed
		//------------------------------------------
		if ( fileMSE!="" )
		{
			std::cout<<"- Write the scalar result image..."<<std::endl;
			typedef itk::ImageFileWriter<FloatImageType>	FloatImageWriterType;
			FloatImageWriterType::Pointer floatWriter = FloatImageWriterType::New();
			floatWriter->SetFileName(fileMSE );
			floatWriter->SetInput(measureImage);
			floatWriter->Update();
		}

		//------------------------------------------
		// If the CC mask is available compute stats
		//------------------------------------------
		if ( ccMask.GetPointer()!=NULL )
		{
			std::cout<<"- Compute statistics on the CC..."<<std::endl;

			typedef itk::LabelStatisticsImageFilter<FloatImageType,LabelImageType> LabelStatisticsFilter;
			LabelStatisticsFilter::Pointer statFilter = LabelStatisticsFilter::New();
			statFilter->SetLabelInput(ccMask);
			statFilter->SetInput(measureImage);
			statFilter->Update();

			//------------------------------------------
			// Threshold the MSE
			//------------------------------------------
			double threshold = statFilter->GetMean(1)+ sigmaFactor * statFilter->GetSigma(1);
			std::cout<<"  Mean = "<<statFilter->GetMean(1)<<" - Std = "<<statFilter->GetSigma(1)<<std::endl;
			std::cout<<"  Sigma factor = "<< sigmaFactor << std::endl;
			std::cout<<"  Threshold = "<<threshold<<std::endl;

			typedef itk::BinaryThresholdImageFilter<FloatImageType, LabelImageType> BinaryFilterType;
			BinaryFilterType::Pointer thresholdFilter=BinaryFilterType::New();
			thresholdFilter->SetInput(statFilter->GetOutput());
			thresholdFilter->SetLowerThreshold(threshold);
			thresholdFilter->SetInsideValue(1);
			thresholdFilter->SetOutsideValue(0);
			thresholdFilter->Update();

			
			std::cout<<"- Write the 2T mask result..."<<std::endl;
			typedef itk::ImageFileWriter<LabelImageType>	MaskImageWriterType;
			MaskImageWriterType::Pointer maskWriter = MaskImageWriterType::New();
			maskWriter->SetFileName(fileOutputMask );
			maskWriter->SetInput(thresholdFilter->GetOutput());
			maskWriter->Update();
		}
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
