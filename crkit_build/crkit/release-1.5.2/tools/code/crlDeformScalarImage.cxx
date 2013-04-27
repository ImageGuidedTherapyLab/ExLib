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

#include "crlDeformScalarImage.h"

int main(int argc, char *argv[])
{
  // Read in the images to be analysed and instantiate the appropriate types
  static int const ImageDimension = 3;
  typedef   float PixelType;
  typedef   itk::Image<PixelType,ImageDimension>  ImType;
  typedef   itk::ImageFileReader< ImType >    ImReaderType;

  crlDeformScalarImageBase::Pointer crlDeformScalarImagePtr = 0;

  if ((argc != 4) && (argc != 5)) {
    std::cout << "Usage: " << argv[0] << " inputImageFile " <<
     "inputDeformationField outputImageFile ?nearest|linear|spline|sinc?" <<
     std::endl;
    return EXIT_FAILURE;
  }

  ImType::Pointer image1;
  image1 = ImType::New();
  ImReaderType::Pointer image1reader;
  image1reader = ImReaderType::New();
  image1reader->SetFileName( argv[1] );

  try {
    image1reader->GenerateOutputInformation();
  } catch ( itk::ExceptionObject &err )
  {
    std::cout << "ExceptionObject caught!" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }

  itk::ImageIOBase *imageIO = image1reader->GetImageIO();
  if (imageIO->GetPixelType() == itk::ImageIOBase::SCALAR) {
     switch(imageIO->GetComponentType()) {
       case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE :
       case itk::ImageIOBase::UCHAR :
       case itk::ImageIOBase::CHAR :
       case itk::ImageIOBase::USHORT :
       case itk::ImageIOBase::SHORT :
       case itk::ImageIOBase::UINT :
       case itk::ImageIOBase::INT :
       case itk::ImageIOBase::ULONG :
       case itk::ImageIOBase::LONG :
       case itk::ImageIOBase::FLOAT :
         if (imageIO->GetNumberOfDimensions() == 2) {
crlDeformScalarImagePtr = crlDeformScalarImage< 2, float >::New();
         } else if (imageIO->GetNumberOfDimensions() == 3) {
crlDeformScalarImagePtr = crlDeformScalarImage< 3, float >::New();
         } else if (imageIO->GetNumberOfDimensions() == 4) {
crlDeformScalarImagePtr = crlDeformScalarImage< 4, float >::New();
         } else if (imageIO->GetNumberOfDimensions() == 5) {
crlDeformScalarImagePtr = crlDeformScalarImage< 5, float >::New();
         } else if (imageIO->GetNumberOfDimensions() == 6) {
crlDeformScalarImagePtr = crlDeformScalarImage< 6, float >::New();
         }
         break;
       case itk::ImageIOBase::DOUBLE :
         if (imageIO->GetNumberOfDimensions() == 2) {
crlDeformScalarImagePtr = crlDeformScalarImage< 2, double >::New();
         } else if (imageIO->GetNumberOfDimensions() == 3) {
crlDeformScalarImagePtr = crlDeformScalarImage< 3, double >::New();
         } else if (imageIO->GetNumberOfDimensions() == 4) {
crlDeformScalarImagePtr = crlDeformScalarImage< 4, double >::New();
         } else if (imageIO->GetNumberOfDimensions() == 5) {
crlDeformScalarImagePtr = crlDeformScalarImage< 5, double >::New();
         } else if (imageIO->GetNumberOfDimensions() == 6) {
crlDeformScalarImagePtr = crlDeformScalarImage< 6, double >::New();
         }
         break;
       default : std::cerr << "No template for crlDeformScalarImage " <<
              imageIO->GetNumberOfDimensions() << " " <<
              imageIO->GetComponentTypeAsString(imageIO->GetComponentType())
              << std::endl;
     }
   }

  // In case the above statement doesn't initialize the pointer
  if (!crlDeformScalarImagePtr) {
    std::cout << "crlDeformScalarImagePtr not implemented for pixel type " <<
        imageIO->GetPixelTypeAsString(imageIO->GetPixelType()) <<
         " with dimensions " << imageIO->GetNumberOfDimensions() <<
        std::endl;
    exit(1);
  }

  crlDeformScalarImagePtr->SetInputFileName( argv[1] );
  crlDeformScalarImagePtr->SetDeformationFieldFileName( argv[2] );
  crlDeformScalarImagePtr->SetOutputFileName( argv[3] );
  if (argc == 5) {
    crlDeformScalarImagePtr->SetInterpolationType( argv[4] );
  }

  int ret = crlDeformScalarImagePtr->Execute();
  if (ret != 0) {
    std::cerr << "crlDeformScalarImage returned an error code " 
              << ret << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

}
