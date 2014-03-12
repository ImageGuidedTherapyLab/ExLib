#include "configuration.h"
#include "tclap/CmdLine.h"
#include <iostream>

#include "crlFileName.h"
#include "crlTaskProgress.h"
#include "crlDocumentedToolOutput.h"

#include "crlConvertAndWriteScalarImage.h"
#include "crlDWIStudy.h"

#include "crlMFMTensorEstimateFilter.h"
#include "itkOneTensorEstimateFilter.h"	//tmp?
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkTimeProbe.h"

#include "GlobalDefinitions.h"

#include "clock.h"

#if defined(WIN32) && defined(_DEBUG)
  #define _CRTDBG_MAP_ALLOC
  #include <stdlib.h>
  #include <crtdbg.h>
  #define DEBUG_NEW new( _NORMAL_BLOCK, __FILE__, __LINE__ )
  #define new DEBUG_NEW
#endif


// Just for the boundary condition
crl::DWI::SInternalVoxelData* NumericTraits<crl::DWI::SInternalVoxelData* >::Zero = NULL;


using namespace std;

/**********************************************************************************************//**
 * \fn	int main(int argc, char **argv)
 *
 * \brief	Main entry-point for this application. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \param	argc	Number of command-line arguments. 
 * \param	argv	Array of command-line argument strings. 
 *
 * \return	Exit-code for the process - 0 for success, else an error code. 
*************************************************************************************************/
int main(int argc, char **argv)
{
	string nhrdfile, tensorfile, input1T, mask2Tfile, optAlgo;
	bool residuals, radialBL, radialIso, randomInit, multiRestart, useIsotropic, outputInitTensors, firstLowB, estFractions, useRobust, medinria1T, videoMode;
	unsigned int nbproc, nbTensors, maxIter, verboseLevel, overboseLevel;
	float threshold, regWeight, anisoK, robustParam, fastdiffusionfraction;
	float ricianNoiseVariance, pruneThreshold;
	string bbox;

	std::vector<std::string> constraints;

#if defined(_MSC_VER) && defined(_DEBUG)
    //_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF  );
#endif

	/*-------------------------------------------------
	 First parse the command line (with TCLAP)
	-------------------------------------------------*/
	try {
		TCLAP::CmdLine cmd("", ' ', CRKIT_VERSION_STRING );
		crl::DocumentedToolOutput crlToolOutput(
			"CRL Diffusion Weighted Imaging ToolBox.\n"
			"Estimate two-tensor models.\n"
			"For better results, the RAW data should be eddy-current corrected."
			"For faster results, the RAW data should be masked with a brain mask.",
			"Benoit Scherrer", 
			"2010" );
		cmd.setOutput( &crlToolOutput );


		// Accept only a limited number of options for the optimization algo
		std::vector<std::string> allowedAlgoStr = crl::OptimizerUtils::GetOptimizerSetStringVector(crl::OptimizerUtils::OPTIMIZER_POWELL + crl::OptimizerUtils::OPTIMIZER_NEWUOA + crl::OptimizerUtils::OPTIMIZER_BOBYQA);
		TCLAP::ValuesConstraint<std::string> allowedAlgoTCLAP( allowedAlgoStr );
		TCLAP::ValueArg<std::string> argOptimiAlgo("a","algo","Optimization algorithm",false,"newuoa",&allowedAlgoTCLAP, cmd);

		TCLAP::SwitchArg				argVideoMode	("","video",		"Output intermediate multi-tensor solutions at each iteration to further produce a video.", cmd, false);
		TCLAP::SwitchArg				argRandomInit	("","randominit",		"Use random angle for two-tensors initialization.", cmd, false);
		TCLAP::SwitchArg				argMultiRestart	("","multirestart",		"Restart several times the optimization to ensure a better minimum (has effect only in non-regularized mode).", cmd, false);
		TCLAP::SwitchArg				argFirstLowB	("","firstlowb",	"If set makes a first pass using only the low b-values.", cmd, false);
		
		TCLAP::SwitchArg				argRadialBL	("","radialBL",	"If set, estimates a radial baseline for each tensor.", cmd, false);
		TCLAP::SwitchArg				argRadialIso("","radialIso",	"If set, estimates a radial isotropic part for each tensor.", cmd, false);
		

		TCLAP::ValueArg<float>			argPrune	("","prune", "Prune tensors whose fractions are smaller than the given number (default: 0.025)", false, 0.025, "minimum fraction", cmd);

		TCLAP::SwitchArg				argResiduals	("","residuals",	"If set, export the residuals in the vector image [outputbasefilename]_residuals.nrrd.", cmd, false);


		TCLAP::ValueArg<float>			argFastDiffusionFraction	("","fastdfrac", "Sets the fast diffusion fraction", false, 0, "fraction", cmd);
		TCLAP::ValueArg<float>			argRicianNoise	("","rician",		"Use Rician noise modeling and set the noise variance", false, 0, "noise variance", cmd);
		TCLAP::ValueArg<unsigned int>	argMaxIter      ("","maxiter",		"Maximum number of iterations",	false, 3000, "MaxIter", cmd);
		TCLAP::ValueArg<string>			argBBox			("","bbox",			"Crop the input image with a bounding box before computing the tensor models",	false, "", "x,y,z,sx,sy,sz", cmd);
		TCLAP::SwitchArg				argOutputInitT	("","inittensors",	"If set output the multi-tensor models used to initialize the minimization.", cmd, false);
		TCLAP::ValueArg<unsigned int>	argV			("v","verbose",		"Set the verbose level (0:no verbose)", false, 0, "level", cmd);
		TCLAP::ValueArg<unsigned int>	argOV			("","overbose",		"Set the optimization verbose level (0:no verbose)", false, 0, "level", cmd);
		TCLAP::ValueArg<float>			argThreshold	("t","threshold",	"Set the threshold on the B0 image", false, 80, "value", cmd);
		TCLAP::ValueArg<unsigned int>	argProc			("p","proc",		"Set the number of processors",	false, 1, "Nb Proc", cmd);

		//-----------------------------------------------
		// Define the constraints
		//-----------------------------------------------
		std::vector<std::string> allowedConstraintsStr;
		allowedConstraintsStr.push_back("none");
		allowedConstraintsStr.push_back("majorEV");
		allowedConstraintsStr.push_back("minorEV");
		allowedConstraintsStr.push_back("cylinder");
		allowedConstraintsStr.push_back("sameplane");
		TCLAP::ValuesConstraint<std::string> allowedConstraintsTCLAP( allowedConstraintsStr );
		TCLAP::MultiArg<std::string> argConstraints("c","constraint","Set the constraints used for the estimation (several constraints can be defined by several -c options). Default: cylinder.",false,&allowedConstraintsTCLAP, cmd);

		TCLAP::SwitchArg				argF			("f","fractions",	"If set estimate the fractions", cmd, false);
		TCLAP::SwitchArg				argIsotropic	("e","pve",			"If set estimation an isotropic compartment", cmd, false);
		TCLAP::ValueArg<float>			argRobustParam	("","robustParam",		"Set the Geman-McLure parameter (C=robustParam * median||.||)", false, 3.0f, "value", cmd);
		TCLAP::SwitchArg				argRobust		("r","robust",		"If set use robust estimation", cmd, false);
		TCLAP::ValueArg<float>			argAnisoK	("","anisoK",	"Set the value of the parameter K for the anisotropic regularisation", false, 0.001, "value", cmd);
		TCLAP::ValueArg<float>			argReg			("","reg",			"If defined, use regularization and set its weight", false, 0, "Weight", cmd);
		TCLAP::ValueArg<unsigned int>	argNbTensors	("n","ntensors",	"Set the number of tensors (1 or 2)",	false, 1, "Nb Tensors", cmd);
		TCLAP::SwitchArg				argInput1TMI	("","1medinria",	"Set to indicate that the input 1T tensor comes from MedINRIA (divides the eigen values by 1000)", cmd, false);
		TCLAP::ValueArg<string>			arg2TMask		("","2Tmask",		"Set the mask providing the voxels where to estimate two-tensors.",	false, "", "image", cmd);
		TCLAP::ValueArg<string>			argInput1T		("1","input1T",		"If set define the one-tensor solution used to initialize the multi-tensor estimation.",	false, "", "Tensor file", cmd);
		TCLAP::ValueArg<string>			argOutput		("o","output",		"Set the output tensor base filename.",	true, "", "Tensor file", cmd);
		TCLAP::ValueArg<string>			argInput		("i","input",		"Set the input DWI raw NHRD file.",	true, "", "NHDR file", cmd);
	
		// Parse !
		cmd.parse(argc,argv);

		// Get parameters
		nhrdfile				= argInput.getValue();
		tensorfile				= argOutput.getValue();
		verboseLevel			= argV.getValue();
		nbproc					= argProc.getValue();
		nbTensors				= argNbTensors.getValue();
		useIsotropic			= argIsotropic.getValue();
		constraints				= argConstraints.getValue();
		useRobust				= argRobust.getValue();
		outputInitTensors		= argOutputInitT.getValue();
		estFractions			= argF.getValue();
		maxIter					= argMaxIter.getValue();
		optAlgo					= argOptimiAlgo.getValue();
		bbox					= argBBox.getValue();
		threshold				= argThreshold.getValue();
		overboseLevel			= argOV.getValue();
		regWeight				= argReg.getValue();
		medinria1T				= argInput1TMI.getValue();
		input1T					= argInput1T.getValue();
		firstLowB				= argFirstLowB.getValue();
		ricianNoiseVariance		= argRicianNoise.getValue();
		anisoK					= argAnisoK.getValue();
		robustParam				= argRobustParam.getValue();
		videoMode				= argVideoMode.getValue();
		multiRestart			= argMultiRestart.getValue();
		randomInit				= argRandomInit.getValue();
		fastdiffusionfraction	= argFastDiffusionFraction.getValue();
		radialBL				= argRadialBL.getValue();
		radialIso				= argRadialIso.getValue();
		pruneThreshold			= argPrune.getValue();
		residuals				= argResiduals.getValue();
		mask2Tfile				= arg2TMask.getValue();

	}
	catch ( TCLAP::ArgException& e)
	{
		std::cerr << "Error: "<<e.error()<< " for argument " <<e.argId() << std::endl;
		exit(1);
	}

	//=====================================================
	//	Beginning of the program !
	//=====================================================
	try {
		typedef Image< unsigned short, 3 >            GradientImageType;
		typedef VectorImage< unsigned short, 3 >      GradientImagesType;
		typedef Image< unsigned char, 3 >             MaskImageType;
		typedef crl::DWI::MFMTensorEstimateFilter< unsigned short, float > TensorEstimatorFilterType;

		//-----------------------------------------------
		// Some checking...
		//-----------------------------------------------
		if ( videoMode && nbproc!=1 )
		{
			cout<<"The video mode only works with one thread (-p 1). Cannot continue."<<endl;
			exit(1);
		}

		//-----------------------------------------------
		// Some default parameters
		//-----------------------------------------------
		if ( constraints.size()==0 && nbTensors>=2 )
			constraints.push_back("cylinder");

		//=====================================================
		// Load the input DWI data
		//=====================================================
		cout<<"Load..."<<endl;
		crl::DWIStudy<unsigned short> inputStudy;
		inputStudy.LoadStudy(nhrdfile);

		//-----------------------------------------------
		// Prepare to crop the DWI data if necessary
		//-----------------------------------------------
		typedef itk::RegionOfInterestImageFilter<GradientImagesType, GradientImagesType> GradientImageExtractType;
		GradientImageExtractType::InputImageRegionType cropRegion = inputStudy.DWIData->GetLargestPossibleRegion ();
		GradientImagesType::Pointer InputImages = inputStudy.DWIData;

		//-----------------------------------------------
		// If a bounding box was specified, crop the input images
		//-----------------------------------------------
		int x,y,z,sx,sy,sz;
		GradientImageExtractType::Pointer extractor1 = GradientImageExtractType::New();
		bool dataCropped=false;
		if ( bbox!="" )
		{
			if ( sscanf(bbox.c_str(),"%d,%d,%d,%d,%d,%d", &x, &y, &z, &sx, &sy, &sz)==6)
			{
				cout<<"Crop..."<<endl;
				cropRegion.SetIndex(0, x);
				cropRegion.SetIndex(1, y);
				cropRegion.SetIndex(2, z);
				cropRegion.SetSize(0, sx);
				cropRegion.SetSize(1, sy);
				cropRegion.SetSize(2, sz);

				// Crop
				extractor1->SetInput(inputStudy.DWIData);
				extractor1->SetRegionOfInterest(cropRegion);
				extractor1->Update();

				InputImages = extractor1->GetOutput();
				dataCropped = true;
			}
		}

		//=====================================================
		// One tensor estimate
		//=====================================================

		//-----------------------------------------------
		// If a single tensor solution was not specified,
		// compute it with the itk filter
		//
		// WARNING. ONLY WORKS WITH A SINGLE SHELL ACQUISITION
		//-----------------------------------------------
		TensorEstimatorFilterType::TensorImageType::Pointer OneTensorImage;
		if ( input1T=="" )
		{
			//----------------------------------------
			// Convert the gradient list to a itk container
			//----------------------------------------
			typedef vnl_vector_fixed< double,3 > GradientDirectionType;  
			typedef itk::VectorContainer< unsigned int, GradientDirectionType > GradientDirectionContainerType;
			GradientDirectionContainerType::Pointer itkGradContainer = GradientDirectionContainerType::New();
			for ( unsigned int i=0; i<inputStudy.Gradients.size(); i++ )
				itkGradContainer->InsertElement( i, inputStudy.Gradients[i] );

			//----------------------------------------
			// Computes the one tensor solution
			//----------------------------------------
			cout<<"- Compute initial one-tensor solution..."<<endl;
			typedef itk::OneTensorEstimateFilter< unsigned short, unsigned short, float > OneTensorEstimateFilterType;
			OneTensorEstimateFilterType::Pointer OneTensorFilter = OneTensorEstimateFilterType::New();
			OneTensorFilter->SetGradientImage( itkGradContainer, InputImages );

			//-----------------------------------------------------
			// 1 thread is necessary until we fix netlib/dsvdc.c
			//-----------------------------------------------------
			OneTensorFilter->SetNumberOfThreads( 1 );
			OneTensorFilter->SetBValue(inputStudy.NominalBValue);  
			OneTensorFilter->SetThreshold( static_cast<OneTensorEstimateFilterType::ReferencePixelType >(threshold));
			OneTensorFilter->Update();

			OneTensorImage = OneTensorFilter->GetOutput();
			OneTensorFilter->GetOutput()->DisconnectPipeline();
		}
		//-----------------------------------------------
		// Else open the specified one-tensor solution
		//-----------------------------------------------
		else
		{
			cout<<"- Open initial one-tensor solution..."<<endl;

			//-----------------------------------------------------
			// Open the file
			//-----------------------------------------------------
			typedef itk::ImageFileReader<TensorEstimatorFilterType::TensorImageType > TensorReaderType;
			TensorReaderType::Pointer tensorReader = TensorReaderType::New();
			tensorReader->SetFileName( input1T );
			tensorReader->Update();
			OneTensorImage = tensorReader->GetOutput();

			if ( OneTensorImage->GetLargestPossibleRegion().GetSize() != inputStudy.DWIData->GetLargestPossibleRegion().GetSize() )
			{
				std::cout<<"Error. The DWI data and the provided one tensor image sizes do not match."<<std::endl;
				exit(1);
			}

			//-----------------------------------------------------
			// Crop the opened one-tensor solution
			//-----------------------------------------------------
			if ( dataCropped )
			{
				typedef itk::RegionOfInterestImageFilter<TensorEstimatorFilterType::TensorImageType, TensorEstimatorFilterType::TensorImageType> TensorImageExtractType;
				TensorImageExtractType::Pointer tensorExtractor = TensorImageExtractType::New();
				tensorExtractor->SetInput(OneTensorImage);
				tensorExtractor->SetRegionOfInterest(cropRegion);
				tensorExtractor->Update();
				OneTensorImage = tensorExtractor->GetOutput();
				tensorExtractor->GetOutput()->DisconnectPipeline();
			}
		}

		// ---------------------------------------
		// Write out the one-tensors image 
		// ---------------------------------------
		typedef itk::ImageFileWriter<TensorEstimatorFilterType::TensorImageType > TensorWriterType;
		TensorWriterType::Pointer tensorWriter = TensorWriterType::New();
		tensorWriter->SetFileName( tensorfile );
		tensorWriter->SetInput( OneTensorImage );
		tensorWriter->Update();


		//=====================================================
		// Load the mask where to estimate two tensors
		//=====================================================
		MaskImageType::Pointer mask2T;
		if ( mask2Tfile!="" )
		{
			if ( nbTensors!=2 )
			{
				cout<<"- Warning: a mask was set with the --2Tmask but this is compatible with only two tensors. Ignore the mask." <<std::endl; 
			}
			else
			{
				cout<<"- Load the two-tensor mask..."<<endl;
				typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
				MaskReaderType::Pointer maskReader = MaskReaderType::New();
				maskReader->SetFileName( mask2Tfile );
				maskReader->Update();
				mask2T = maskReader->GetOutput();
				maskReader->GetOutput()->DisconnectPipeline();

				if ( mask2T->GetLargestPossibleRegion().GetSize() != inputStudy.DWIData->GetLargestPossibleRegion().GetSize() )
				{
					std::cout<<"Error. The DWI data and the two-tensor mask image sizes do not match."<<std::endl;
					exit(1);
				}
				
				//-----------------------------------------------------
				// Crop the mask if needed
				//-----------------------------------------------------
				if ( dataCropped )
				{
					typedef itk::RegionOfInterestImageFilter<MaskImageType, MaskImageType> MaskImageExtractType;
					MaskImageExtractType::Pointer maskExtractor = MaskImageExtractType::New();
					maskExtractor->SetInput(mask2T);
					maskExtractor->SetRegionOfInterest(cropRegion);
					maskExtractor->Update();
					mask2T = maskExtractor->GetOutput();
					maskExtractor->GetOutput()->DisconnectPipeline();
				}
			}
		}

		//=====================================================
		// Now set up the multiple-tensors estimator
		//=====================================================
		cout<<"- Compute Log-Euclidean estimate..."<<endl;
		TensorEstimatorFilterType::Pointer estimateMFM = TensorEstimatorFilterType::New();

		//-----------------------------------------------
		// Set parameters
		//-----------------------------------------------
		estimateMFM->SetVerboseLevel(verboseLevel);
		estimateMFM->SetOVerboseLevel(overboseLevel);
		estimateMFM->SetNumberOfThreads( nbproc );
		estimateMFM->SetNumberOfTensors(nbTensors);
		estimateMFM->SetRobustEstimation(useRobust);

		crl::DWI::EstimationStrategyType estimStrategy = 0;

		bool noConstraints=false;
		for ( unsigned int i=0; i<constraints.size(); i++ )
		{
			if ( constraints[i]=="majorEV" ) estimStrategy += crl::DWI::ESTIMATION_SAME_MAJOR_EIGEN_VALUE;
			if ( constraints[i]=="minorEV" ) estimStrategy += crl::DWI::ESTIMATION_SAME_MINOR_EIGEN_VALUE;
			if ( constraints[i]=="cylinder" ) estimStrategy += crl::DWI::ESTIMATION_EACH_CYLINDER;
			if ( constraints[i]=="sameplane" ) estimStrategy += crl::DWI::ESTIMATION_SAME_PLANE;
			if ( constraints[i]=="none" ) noConstraints=true;
		}
		if ( noConstraints ) estimStrategy=0;

		if ( useIsotropic ) estimStrategy += crl::DWI::ESTIMATION_F_ISO;
		if ( estFractions ) estimStrategy += crl::DWI::ESTIMATION_F_TENSORS;
		if ( radialBL ) estimStrategy += crl::DWI::ESTIMATION_RADIALBASELINE;
		if ( radialIso ) estimStrategy += crl::DWI::ESTIMATION_RADIALISO;


		estimateMFM->SetEstimationStrategy(estimStrategy);

		estimateMFM->SetOutputInitTensors(outputInitTensors);
		estimateMFM->SetMaxIter(maxIter);
		estimateMFM->SetFirstPassWithLowBValuesOnly(firstLowB);
		estimateMFM->SetInitTensorImage(OneTensorImage);
		estimateMFM->SetUseRicianNoiseModeling(ricianNoiseVariance!=0);
		estimateMFM->SetRicianNoiseVariance(ricianNoiseVariance);
		estimateMFM->SetInitialOneTensorFromMedINRIA(medinria1T);
		estimateMFM->SetRegularizationAnisotropicK(anisoK);
		estimateMFM->SetRobustParam((double)robustParam);
		estimateMFM->SetThreshold((unsigned short)threshold);
		estimateMFM->SetOptimizationMultiRestart(multiRestart);
		estimateMFM->SetVideoMode(videoMode);
		estimateMFM->SetRandomInit(randomInit);
		estimateMFM->SetFastDiffusionFraction(fastdiffusionfraction);
		estimateMFM->SetPruneThreshold(pruneThreshold);
		estimateMFM->SetExportResiduals(residuals);
		if ( mask2T.GetPointer()!=NULL ) estimateMFM->Set2TMaskImage(mask2T);

		//-----------------------------------------------
		// Set up the regularization
		//-----------------------------------------------
		if ( regWeight>0 )
		{
			estimateMFM->SetUseRegularization(true);
			estimateMFM->SetRegularizationWeight(regWeight);
		}
		else
			estimateMFM->SetUseRegularization(false);
	
		//-----------------------------------------------
		// Set up the optimization algorithm
		//-----------------------------------------------
		estimateMFM->SetOptimizationAlgo(crl::OptimizerUtils::GetOptimizerIdFromString(optAlgo));



		//=====================================================
		// MFM Estimator
		//=====================================================

		//-----------------------------------------------
		// Init elapsed time
		//-----------------------------------------------
		itk::TimeProbe elapsedTime;
		elapsedTime.Start(); 		
		Clock ti; ti.start();
	  
		//-----------------------------------------------
		// ESTIMATE THE TENSORS!
		//-----------------------------------------------
		inputStudy.DWIData = InputImages;		//if cropped...
		estimateMFM->SetInputStudy( inputStudy );
		estimateMFM->Update();

		//-----------------------------------------------
		// Show elapsed time
		//-----------------------------------------------
		elapsedTime.Stop(); ti.stop();
		string strProcessingTime = crl::TaskProgress::ConvertTimeToString(elapsedTime.GetTotal());
		cout << endl<<"Total computation time: " << strProcessingTime.c_str() << endl;

		char szBuffer[128];
		ti.formatTime(ti.laps(), szBuffer);
		string ttime = "Processor-Time : "+string(szBuffer) + ", Real Time elapsed: ";
		ti.formatTime(ti.real_laps(), szBuffer);
		ttime = ttime + string(szBuffer) +".";
		cout << ttime.c_str() << endl;
			

		//=====================================================
		// OUTPUT
		//=====================================================

		//-----------------------------------------------
		// Output infos in a log file 
		//-----------------------------------------------
		std::fstream infofile;
		crl::FileName infofilename(tensorfile);
		string info_filename = infofilename.getPath()+infofilename.getFileName()+"_info.txt";
		infofile.open ( info_filename.c_str(), std::fstream::out);
		if  ( infofile.is_open() ) 
		{
			infofile<<"--------------------------------"<<endl;
			infofile<<" crlMFMTensorEstimate"<<endl;
			infofile<<" (c) Benoit SCHERRER, 2010"<<endl;
			infofile<<"--------------------------------"<<endl;
			infofile<< endl<<"COMPUTATION TIME: " << elapsedTime.GetTotal() << "s     ("<<strProcessingTime.c_str()<<")"<<endl;
			infofile<< ttime.c_str()<<endl;

			infofile<<endl<< "BOUNDING BOX : " << bbox.c_str() << endl << endl;

			estimateMFM->Print(infofile);
			infofile.close();
		}

		//-----------------------------------------------
		// Write out the tensor files
		//-----------------------------------------------
		for ( unsigned int k=0 ; k<nbTensors ; k++ )
		{
			crl::FileName ofile(tensorfile);
			char szBuffer[256]; sprintf(szBuffer, "_t%d", k);

			TensorWriterType::Pointer tensorWriter2 = TensorWriterType::New();
			tensorWriter2->SetFileName( ofile.getCompleteFilePath_WithSufix(string(szBuffer)) );
			tensorWriter2->SetInput( estimateMFM->GetOutput(k) );
			tensorWriter2->Update();
		}

		//-----------------------------------------------
		// Write out the fraction file if needed
		//-----------------------------------------------
		if ( useIsotropic || nbTensors>1 )
		{
			crl::FileName ofile(tensorfile);
			typedef itk::ImageFileWriter<itk::Image<float, 4> > FractionsWriterType;
			FractionsWriterType::Pointer fractionsWriter = FractionsWriterType::New();
			fractionsWriter->SetFileName( ofile.getCompleteFilePath_WithSufix("_fractions" ) );
			fractionsWriter->SetInput( estimateMFM->GetFractionsImage() );
			fractionsWriter->Update();

		}

		//-----------------------------------------------
		// Write out the residuals if needed
		//-----------------------------------------------
		if ( residuals )
		{
			crl::FileName ofile(tensorfile);
			typedef itk::ImageFileWriter< TensorEstimatorFilterType::ResidualsImageType > ResidualsWriterType;
			ResidualsWriterType::Pointer residualsWriter = ResidualsWriterType::New();
			residualsWriter->SetFileName( ofile.getCompleteFilePath_WithSufix("_residuals" ) );
			residualsWriter->SetInput( estimateMFM->GetResidualsImage() );
			residualsWriter->Update();
		}


		//-----------------------------------------------
		// Write out the initial multi-tensors 
		// (before the optimization) if needed
		//-----------------------------------------------
		if ( outputInitTensors )
		{
			for ( unsigned int k=0 ; k<nbTensors ; k++ )
			{
				crl::FileName ofile(tensorfile);
				char szBuffer[256]; sprintf(szBuffer, "_init_t%d", k);

				TensorWriterType::Pointer tensorWriter2 = TensorWriterType::New();
				tensorWriter2->SetFileName( ofile.getCompleteFilePath_WithSufix(string(szBuffer)) );
				tensorWriter2->SetInput( estimateMFM->GetInitTensors()[k] );
				tensorWriter2->Update();
			}
		}


	}

	/*-------------------------------------
	If catched an exception, show a message
	-------------------------------------*/
	catch (itk::ExceptionObject& e)
	{
		cout<< "  Exception: "<<endl<< e.GetDescription()<< endl;
		return 1;
	}

  return 0;
}
