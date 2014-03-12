/*
 * Copyright 2007 Simon Warfield simon.warfield@childrens.harvard.edu
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#include "crlBinaryMorphology.h"

int main(int argc, char *argv[])
{

  typedef unsigned int PixelType;

  if (argc != 6) {
    std::cout << 
      "Usage: " << argv[0] << " inputFileName operation "
      << "foregroundValue seRadius outputFileName "
      << std::endl;
    return EXIT_FAILURE;
  }

  crlBinaryMorphologyBase::Pointer crlBinaryMorphologyPtr = 0;

  typedef itk::Image<PixelType, 3 >  ImageType;
  typedef itk::ImageFileReader< ImageType >    ImageReaderType;

  ImageType::Pointer image1;
  image1 = ImageType::New();

  ImageReaderType::Pointer image1reader;
  image1reader = ImageReaderType::New();
  image1reader->SetFileName( argv[1] );

  try {
    image1reader->GenerateOutputInformation();
  } catch (itk::ExceptionObject &err )
  {
    std::cout << "Exception caught." << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }

  itk::ImageIOBase *imageIO = image1reader->GetImageIO();
  if (imageIO->GetNumberOfDimensions() == 1) {
    crlBinaryMorphologyPtr = crlBinaryMorphology<1, PixelType>::New();
  } else if (imageIO->GetNumberOfDimensions() == 2) {
    crlBinaryMorphologyPtr = crlBinaryMorphology<2, PixelType>::New();
  } else if (imageIO->GetNumberOfDimensions() == 3) {
    crlBinaryMorphologyPtr = crlBinaryMorphology<3, PixelType>::New();
  } else if (imageIO->GetNumberOfDimensions() == 4) {
    crlBinaryMorphologyPtr = crlBinaryMorphology<4, PixelType>::New();
  };

  if (!crlBinaryMorphologyPtr) {
    std::cerr << "No template specialization for dimensionality " 
                    << imageIO->GetNumberOfDimensions() << std::endl;
    return EXIT_FAILURE;
  };

  crlBinaryMorphologyPtr->SetInputFileName(  argv[1] );
  crlBinaryMorphologyPtr->SetOutputFileName( argv[5] );

  // CHECK THE OPERATION TYPE 
  std::string *operation = new std::string( argv[2] );
  if ( (operation->compare(0,5,"dilat") == 0) )         {
    crlBinaryMorphologyPtr->SetOperationToDilation( );
  } else if ( (operation->compare(0,4,"erod") == 0)       ) {
    crlBinaryMorphologyPtr->SetOperationToErosion( );
  } else if ( (operation->compare(0,4,"open") == 0)       ) {
    crlBinaryMorphologyPtr->SetOperationToOpening( );
  } else if ( (operation->compare(0,4,"clos") == 0)       ) {
    crlBinaryMorphologyPtr->SetOperationToClosing( );
  } else {
    std::cerr << "crlBinaryMorphology doesn't know operation " <<
      argv[2] << " - instead use dilate erode opening closing " << 
      std::endl;
    return EXIT_FAILURE;
  }

  crlBinaryMorphologyPtr->SetForegroundValue( std::atof(argv[3]) );
  crlBinaryMorphologyPtr->SetRadius(          std::atoi(argv[4]) );

  int ret = crlBinaryMorphologyPtr->Execute();

  if (ret != EXIT_SUCCESS) 
  {
    std::cerr << argv[0] << "Returned error code " << ret << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

