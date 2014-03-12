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

/* Generate a certain 3D volume with randomly assigned voxel values. */

#ifdef WIN32
#include <missingFunctionsWindows.h>
#else
#include <sys/time.h>
#endif

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

#include "configuration.h"
#include <tclap/CmdLine.h>
#include <crlDocumentedToolOutput.h>

int main(int argc, char *argv[])
{

  try {
     TCLAP::CmdLine cmd("", ' ', CRKIT_VERSION_STRING );
     crl::DocumentedToolOutput my("Generate a random image based on a predetermined true segmentation simulated a rater with the given sensitivity and specificity.]n", "Simon K. Warfield", "2010");
     my.SetKeepOriginalDataType( true );
     cmd.setOutput( &my );

     TCLAP::ValueArg<unsigned int> nrowArg("","nrow", "Number of rows", 
       false, 256, "integer", cmd);
     TCLAP::ValueArg<unsigned int> ncolArg("","ncol", "Number of cols", 
       false, 256, "integer", cmd);
     TCLAP::ValueArg<unsigned int> nsliceArg("","nslice", "Number of slices", 
       false, 1, "integer", cmd);
     TCLAP::ValueArg<double> sensitivityArg("p","sensitivity", "Sensitivity", 
       false, 0.95, "double", cmd);
     TCLAP::ValueArg<double> specificifyArg("q","specificity", "Specificity", 
       false, 0.90, "double", cmd);
     TCLAP::ValueArg<std::string> argOutput("o","output", "Output file", 
       true, "output.nrrd", "file name", cmd);
  }
  catch ( TCLAP::ArgException& e)
  {
    std::cerr << "Error: "<<e.error()<< " for argument " <<
        e.argId() << std::endl;
    exit(1);
  }


  if (argc != 7) {
    fprintf(stdout, "Usage: crlRandomSegmentationImageGenerator ncol nrow nslice p q imageout \n");
    exit(1);
  }

  static int const ImageDimension =  3;
  typedef float PixelType;

  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer randGen =
    itk::Statistics::MersenneTwisterRandomVariateGenerator::New();

  // select a seed based on the clock, to microsecond precision.
  unsigned int seed = 0;
  struct timeval preciseTimeStart;
  gettimeofday(&preciseTimeStart, NULL);
  seed = (unsigned int)(preciseTimeStart.tv_usec & 0xFFFFFFFF);

  randGen->SetSeed(seed);

  typedef   itk::Image<PixelType,ImageDimension>  ImageType;
  typedef   itk::ImageFileReader< ImageType >    ImageReaderType;
  typedef   itk::ImageFileWriter< ImageType >    ImageWriterType;

  typedef ImageType::RegionType RegionType;
  typedef ImageType::SizeType   SizeType;
  typedef ImageType::IndexType  IndexType;
  typedef itk::ImageRegionIteratorWithIndex< ImageType >
                     RegionIteratorType;

  ImageWriterType::Pointer writer = ImageWriterType::New();

  // Make an image of the desired size.
  ImageType::Pointer image = ImageType::New();

  ImageType::SpacingType spacing;
  spacing.Fill(1.0);
  image->SetSpacing( spacing );

  ImageType::PointType origin;
  origin.Fill(0.0);
  image->SetOrigin( origin );

  ImageType::SizeType size;
  size[0] = atoi(argv[1]);
  size[1] = atoi(argv[2]);
  size[2] = atoi(argv[3]);

  ImageType::IndexType start;
  start.Fill(0);
  ImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );

  image->SetRegions( region );
  image->Allocate();

  ImageType::PixelType  initialValue = 0.0;
  image->FillBuffer( initialValue );

  float p = atof(argv[4]);
  float q = atof(argv[5]);

  IndexType index;
  if (index.GetIndexDimension() == 3) {
    for (index[2] = start[2]; index[2] < (unsigned int)(start[2] + size[2]); index[2]++) {
      for (index[1] = start[1]; index[1] < (unsigned int)(start[1] + size[1]); index[1]++) {
        for (index[0] = start[0]; index[0] < (unsigned int)(start[0]+size[0]); index[0]++) {
          double ra =  randGen->GetVariate();
          if ( index[0] < static_cast<float>(size[0])/2.0 ) {
            if (ra < q) { // Correct answer will be 0.0
              (*image)[index] = 0.0;
	      // lambdai[index] = log((1-p)/q));
            } else {
              (*image)[index] = 1.0;
	      // lambdai[index] = log((p)/(1-q));
            }
          } else {          // Correct answer will be 1.0
            if (ra < p) {
              (*image)[index] = 1.0;
              // lambdai[index] = log((p)/(1-q));
            } else {
              (*image)[index] = 0.0;
              //lambdai[index] = log((1-p)/(q));
            }
          }
        }
      }
    }
  } else {
    std::cerr << "Not implemented for dimension " << index.GetIndexDimension() << std::endl;
    exit(1);
  }

  try
  {
    writer->SetFileName ( argv[6] );
    writer->SetInput( image );
    writer->Update();
  }
  catch( itk::ExceptionObject & exp )
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
  }

/*
  if (!lambdai.writeVolume(argv[7])) {
	  fprintf(stderr, "Failed to write image out to %s.\n", argv[7]);
	  exit(1);
  }
*/

  return 0;
}
