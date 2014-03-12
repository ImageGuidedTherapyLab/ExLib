/*
 * Copyright 2010 Simon Warfield simon.warfield@childrens.harvard.edu
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
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#include "crlDensityWeightedStats.h"
#include <tclap/CmdLine.h>
#include "configuration.h"

int main( int argc, char ** argv )
{
  static int const ImageDimension = 3;
  typedef   float PixelType;

  typedef   itk::Image<PixelType,ImageDimension>  ImageType;
  typedef   itk::ImageFileReader< ImageType >    ImageReaderType;

  std::string *densityImageFileName = new std::string("");
  std::string *scalarImageFileName = new std::string("");

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ', CRKIT_VERSION_STRING );
      TCLAP::UnlabeledValueArg<std::string> densityImageFileArg("densityImageFileName","Density Image File Name",true,"","density image file name",cmd);
      TCLAP::UnlabeledValueArg<std::string> scalarImageFileArg("scalarImageFileName","Scalar Image File Name",true,"","scalar image file name",cmd);

      cmd.parse(argc,argv);

      if (densityImageFileArg.isSet()) {
        densityImageFileName = new std::string(densityImageFileArg.getValue());
      }
      if (scalarImageFileArg.isSet()) {
        scalarImageFileName = new std::string(scalarImageFileArg.getValue());
      }

      }
  catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(EXIT_FAILURE);
    }

  crl::DensityWeightedStatsBase::Pointer crlDensityWeightedStatsPtr = 0;
  ImageType::Pointer image1;
  image1 = ImageType::New();

  ImageReaderType::Pointer image1reader;
  image1reader = ImageReaderType::New();
  image1reader->SetFileName( densityImageFileName->c_str() );

  try {
    image1reader->GenerateOutputInformation();
  } catch ( itk::ExceptionObject &err )
  {
    std::cout << "ExceptionObject caught!" << std::endl;
    std::cout << err << std::endl;
    return -1;
  }

  itk::ImageIOBase *imageIO = image1reader->GetImageIO();
  if (imageIO->GetPixelType() == itk::ImageIOBase::SCALAR) {
     switch(imageIO->GetComponentType()) {
       case itk::ImageIOBase::UCHAR :
       case itk::ImageIOBase::CHAR :
       case itk::ImageIOBase::USHORT :
       case itk::ImageIOBase::SHORT :
       case itk::ImageIOBase::UINT :
       case itk::ImageIOBase::INT :
       case itk::ImageIOBase::ULONG :
       case itk::ImageIOBase::LONG :
       case itk::ImageIOBase::FLOAT :
         if (imageIO->GetNumberOfDimensions() == 1) {
crlDensityWeightedStatsPtr = crl::DensityWeightedStats< 1, float >::New();
         } else if (imageIO->GetNumberOfDimensions() == 2) {
crlDensityWeightedStatsPtr = crl::DensityWeightedStats< 2, float >::New();
         } else if (imageIO->GetNumberOfDimensions() == 3) {
crlDensityWeightedStatsPtr = crl::DensityWeightedStats< 3, float >::New();
         } else if (imageIO->GetNumberOfDimensions() == 4) {
crlDensityWeightedStatsPtr = crl::DensityWeightedStats< 4, float >::New();
         }
         break;
       case itk::ImageIOBase::DOUBLE :
         if (imageIO->GetNumberOfDimensions() == 1) {
crlDensityWeightedStatsPtr = crl::DensityWeightedStats< 1, float >::New();
         } else if (imageIO->GetNumberOfDimensions() == 2) {
crlDensityWeightedStatsPtr = crl::DensityWeightedStats< 2, float >::New();
         } else if (imageIO->GetNumberOfDimensions() == 3) {
crlDensityWeightedStatsPtr = crl::DensityWeightedStats< 3, float >::New();
         } else if (imageIO->GetNumberOfDimensions() == 4) {
crlDensityWeightedStatsPtr = crl::DensityWeightedStats< 4, float >::New();
         }
         break;
      default : break;
     }
   }

  // In case the switch statement doesn't initialize the pointer
  if (!crlDensityWeightedStatsPtr) {
    std::cout << "crlDensityWeightedStatsPtr not implemented for pixel type " <<
        imageIO->GetPixelTypeAsString(imageIO->GetPixelType()) <<
         " with dimensions " << imageIO->GetNumberOfDimensions() <<
        std::endl;
    return EXIT_FAILURE;
  }

  crlDensityWeightedStatsPtr->SetDensityImageFileName( densityImageFileName->c_str() );
  crlDensityWeightedStatsPtr->SetScalarImageFileName(  scalarImageFileName->c_str() );

  int ret = crlDensityWeightedStatsPtr->Execute();
  if (ret != 0) {
    std::cerr << "crlDensityWeightedStatsPtr returned an error code " << ret <<
    std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

}
