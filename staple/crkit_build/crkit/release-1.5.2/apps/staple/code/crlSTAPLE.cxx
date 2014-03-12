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

#include <iostream>
#include <string>
#include <vector>
#include "crlMSTAPLEImageFilter.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkNumericTraits.h>
#include <itkSize.h>

#include <tclap/CmdLine.h>
#include "configuration.h"

#include "MSTAPLE.h"

int crlSTAPLEImageFilterMain( int argc, char * argv[])
{
  MSTAPLEBase *staple;

  unsigned int fileDimensionality = 0;
  int maxiterations = -1;
  double relativeConvergence = 5e-07;
  
  // This parameter needs adjustment to the problem.
  double stationaryWeight = 0.5; 

  std::string *outputImage = new std::string("");
  std::vector<std::string> files;
  int underflowProtection = 0;
  bool useCompression = false;
  std::string stationaryPrior;
  bool stationaryPriorSet = false;
  std::string initialExpertPerformance;
  bool initialExpertPerformanceSet = false;
  bool assignConsensusVoxels = true;
  bool startAtEStep = false;
  bool stapleMAP = false;
  double alpha = 0.0;
  double beta = 0.0;

  // We use an image definition to read in the initial image, and identify
  // the file type and dimensionality.
  static int const ImageDimension = 3;
  typedef   unsigned char PixelType;
  typedef   itk::Image<PixelType,ImageDimension>  ImageType;
  typedef   itk::ImageFileReader< ImageType >    ImageReaderType;

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',
			 CRKIT_VERSION_STRING );
      TCLAP::ValueArg<std::string> outputImageArg("o", "outputImage",
						  "Output Image File Name",true,"",
						  "Output Image File Name",cmd);

      TCLAP::ValueArg<std::string> stationaryPriorArg("s","stationaryPrior",
						      "optional: Provide the stationary prior for each tissue class label from 0 to the largest in the input files.",false,"","optional: Provide the stationary prior for each tissue class label from 0 to the largest in the input files.",
						      cmd);

      TCLAP::ValueArg<std::string> initialExpertPerformanceArg("",
           "initialExpertPerformance",
           "optional: Provide the initial on-diagonal performance parameter for each input segmentation.",false,"","optional: Provide the initial on-diagonal performance parameter for each input segmentation.", cmd);

      TCLAP::ValueArg<unsigned int> fileDimensionalityArg("f", 
                    "fileDimensionality",
                    "File Dimensionality (2 for 2D, 3 for 3D)",false,3,
                  "optional: file Dimensionality (2 for 2D, 3 for 3D)",cmd);

      TCLAP::ValueArg<int > underflowProtectionArg("u","underflowProtection","Underflow protection : 0 none, 1 strong, 2 extreme.",false,0,"optional: Controls computation with none, some, or extremely extended precision. Useful with large numbers of input segmentations.",cmd);

      TCLAP::ValueArg<double> relativeconvergence("r","relativeConvergence","Relative convergence threshold, used to terminate the E-M algorithm when it has converged. Convergence is defined by relative changes in the mean trace of expert performance below this level.",false,5e-07,"optional relative convergence threshold",cmd);

	  TCLAP::ValueArg<double> stationaryweight("w","stationaryWeight","Stationary prior weight, used to set the weight of stationary prior for each tissue class label with respect to spatially varying prior",false,0.01,"optional stationary prior weight",cmd);

      TCLAP::ValueArg<int> maxiterationsArg("i","maxIterations",
					    "Maximum number of iterations. The E-M algorithm is terminated after this number of iterations has been computed.",false,-1,"optional maximum number of iterations",cmd);

      TCLAP::SwitchArg compressionArg("c","compressOutput",
				      "optional: Write out the reference standard using compression.",
				      cmd, false);

      TCLAP::SwitchArg startAtEStepArg("", "startAtEStep", "Start at the E Step by estimating the reference standard.", cmd, false);

      TCLAP::SwitchArg startAtMStepArg("", "startAtMStep", "Start at the M Step by estimating the parameters from the initial reference standard.", cmd, false);

      // Set to be true by default - useful for label fusion.
      TCLAP::ValueArg<bool> consensusVoxelsArg("", "assignConsensusVoxels", 
              "Assign consensus voxels directly.", false, true, 
       "Determines if voxels with label estimates that are the same amongst all inputs are assigned directly or used in the computation.", cmd);

      TCLAP::SwitchArg mapArg("m","map","optional: Use MAP (maximum a posteriori) estimation.",cmd,false);
      TCLAP::ValueArg<double> alphaArg("a","alpha","optional: Alpha parameter of the beta distribution for MAP.",false,0,"optional alpha parameter",cmd);
      TCLAP::ValueArg<double> betaArg("b","beta","optional: Beta parameter of the beta distribution for MAP.",false,0,"optional beta parameter",cmd);

      TCLAP::UnlabeledMultiArg<std::string> filesArg("files","Input File Names",true,"Input segmentation files: file1 file2 ... fileN",cmd);

      cmd.parse(argc,argv);

      if (maxiterationsArg.isSet()) maxiterations = maxiterationsArg.getValue();
      if (fileDimensionalityArg.isSet()) fileDimensionality = fileDimensionalityArg.getValue();
      if (outputImageArg.isSet()) outputImage = new std::string(outputImageArg.getValue());
      if (stationaryPriorArg.isSet()) {
        stationaryPriorSet = true;
        stationaryPrior = stationaryPriorArg.getValue();
      };
      if (initialExpertPerformanceArg.isSet()) {
        initialExpertPerformanceSet = true;
        initialExpertPerformance = initialExpertPerformanceArg.getValue();
      };

      if (filesArg.isSet()) files = filesArg.getValue();
      useCompression = compressionArg.getValue();
      if (relativeconvergence.isSet()) {
        relativeConvergence = relativeconvergence.getValue();
      }
	  if (stationaryweight.isSet()) {
        stationaryWeight = stationaryweight.getValue();
      }

      if (underflowProtectionArg.isSet()) {
	underflowProtection = underflowProtectionArg.getValue();
      }
      if (consensusVoxelsArg.isSet()) {
        assignConsensusVoxels = consensusVoxelsArg.getValue();
      }
      startAtEStep = startAtEStepArg.getValue();
      // Only one of these two can be set.
      if (startAtMStepArg.isSet()) {
        if (!startAtEStepArg.isSet()) {
          startAtEStep = !startAtMStepArg.getValue();
        }
      }

      if (mapArg.isSet()) {
	stapleMAP = mapArg.getValue();
      }
      if (alphaArg.isSet()) {
	alpha = alphaArg.getValue();
      }
      if (betaArg.isSet()) {
	beta = betaArg.getValue();
      }

      if (stapleMAP && (!alphaArg.isSet() || !betaArg.isSet()))
	{
	  std::cerr << " Alpha and beta have to be set for MAP estimation." << std::endl;
	  return EXIT_FAILURE;
	}

    }
  catch (TCLAP::ArgException &e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
    }

  // Identify the input file dimensionality
  ImageType::Pointer image1;
  image1 = ImageType::New();

  ImageReaderType::Pointer image1reader;
  image1reader = ImageReaderType::New();
  image1reader->SetFileName( files[0] ); // Read data from first input file.

  try {
    image1reader->GenerateOutputInformation();
  } catch ( itk::ExceptionObject &err )
    {
      std::cerr << "ExceptionObject caught!" << std::endl;
      std::cerr << "Failed reading from file " << files[0] << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }

  itk::ImageIOBase *imageIO = image1reader->GetImageIO();
  if (fileDimensionality != 0) { // It has been supplied as an input arg
    if (fileDimensionality != imageIO->GetNumberOfDimensions()) {
      std::cerr << "File dimensionality " << fileDimensionality <<
	" doesn't match that of input file " << files[0] << std::endl;
      return EXIT_FAILURE;
    }
  } else { // Take the dimensionality from the input files.
    fileDimensionality = imageIO->GetNumberOfDimensions();
  }

  if (imageIO->GetPixelType() != itk::ImageIOBase::SCALAR) {
    std::cerr << "Input images are expected to be of scalar type." <<
      std::endl;
    return EXIT_FAILURE;
  }

  switch(fileDimensionality) {
  case 1 : {
    staple = new MSTAPLE<1>;
    break;
  }
  case 2 : {
    staple = new MSTAPLE<2>;
    break;
  }
  case 3 : {
    staple = new MSTAPLE<3>;
    break;
  }
  case 4 : {
    staple = new MSTAPLE<4>;
    break;
  }
  default :  {
    std::cerr << "Data with dimensionality " << fileDimensionality <<
      " is not currently supported" << std::endl;
    return EXIT_FAILURE;
  }
  }

  staple->SetUnderflowProtection( underflowProtection );
  staple->SetUseWriteCompression( useCompression );
  staple->SetAssignConsensusVoxels( assignConsensusVoxels );
  staple->SetStartAtEStep( startAtEStep );

  staple->SetMAPStaple(stapleMAP);
  if (stapleMAP)
    {
      staple->SetMAPAlpha(alpha);
      staple->SetMAPBeta(beta);
    }

  if (stationaryPriorSet) {
    // Convert the string to a vector<double>
    std::vector<double> prior;
    double buffer;
    std::stringstream ss(stationaryPrior);
    while (ss >> buffer) {
      prior.push_back(buffer);
    }
    staple->SetStationaryPrior( prior );/////////////
  }

  if (initialExpertPerformanceSet) {
    // Convert the string to a vector<double>
    std::vector<double> perf;
    double buffer;
    std::stringstream ss(initialExpertPerformance);
    while (ss >> buffer) {
      perf.push_back(buffer);
    }
    staple->SetInitialExpertPerformanceParameters( perf );
  }

  std::vector<std::string>::iterator it;
  for (it = files.begin(); it != files.end(); it++)
    {
      staple->AddFileName( it->c_str() );
    }

  staple->SetOutputFileName( outputImage->c_str() );

  // Now apply the optional arguments to the object.
  if (maxiterations != -1) {
    staple->SetMaximumIterations(maxiterations);
  }

  staple->SetRelativeConvergenceThreshold( relativeConvergence );
  staple->SetStationaryPriorWeight( stationaryWeight );
  if ((stationaryWeight>1) || (stationaryWeight<0))
  {
    std::cerr << "Weight " << fileDimensionality 
                  << " is not between 0 and 1" << std::endl;
    return EXIT_FAILURE;
  }
  // Execute the stapler
  int ret = staple->Execute();
  if ( ret != 0 ) {
    std::cerr << "MSTAPLE failed!  Returned error code " << ret
	      << "." << std::endl;
    return -3;
  }

  std::cout << "Number of elapsed iterations = "
            << staple->GetElapsedIterations() << std::endl;

  //  std::cout.precision(5);
  //  std::cout.setf(ios::fixed, ios::floatfield);

  delete staple;

  return 0;
}

int main(int argc, char *argv[])
{
  return crlSTAPLEImageFilterMain( argc, argv );
}


