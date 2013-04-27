/*
 * Copyright 2008 Simon Warfield simon.warfield@childrens.harvard.edu
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 * http://crl.med.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#include "crlExtractSmallerImageFromImage.h"

// These macros require certain variables to have been initialized 
// by calls to imageIOBase
#define INSTANTIATECLASS(dim, scalar) do { \
           if (imageIODimensionality == dim) { \
             crlExtractSmallerImageFromImagePtr = \
       crl::ExtractSmallerImageFromImage< dim, scalar >::New(); \
           } \
         } while (0)

#define INSTANTIATECLASSOVERDIM(scalar) do { \
    INSTANTIATECLASS(2 , scalar); \
    INSTANTIATECLASS(3 , scalar); \
    INSTANTIATECLASS(4 , scalar); \
    INSTANTIATECLASS(5 , scalar); \
    INSTANTIATECLASS(6 , scalar); \
  } while (0)

#define INSTANTIATECLASSOVERDIMANDTYPE(componentType) do { \
  switch(componentType) { \
    case itk::ImageIOBase::UCHAR : \
         INSTANTIATECLASSOVERDIM(unsigned char); break; \
    case itk::ImageIOBase::CHAR : \
         INSTANTIATECLASSOVERDIM(char); break; \
    case itk::ImageIOBase::USHORT : \
         INSTANTIATECLASSOVERDIM(unsigned short); break; \
    case itk::ImageIOBase::SHORT : \
         INSTANTIATECLASSOVERDIM(signed short); break; \
    case itk::ImageIOBase::UINT :  \
         INSTANTIATECLASSOVERDIM(unsigned int); break; \
    case itk::ImageIOBase::INT :   \
         INSTANTIATECLASSOVERDIM(signed int); break; \
    case itk::ImageIOBase::ULONG : \
         INSTANTIATECLASSOVERDIM(unsigned long); break; \
    case itk::ImageIOBase::LONG :  \
         INSTANTIATECLASSOVERDIM(signed long); break; \
    case itk::ImageIOBase::FLOAT : \
         INSTANTIATECLASSOVERDIM(float); break; \
    case itk::ImageIOBase::DOUBLE : \
         INSTANTIATECLASSOVERDIM(double); break; \
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE : \
    default : break; \
  } \
  } while (0)
  

int main(int argc, char *argv[])
{

  itk::ImageIOBase::IOPixelType pixelType;
  itk::ImageIOBase::IOComponentType componentType;
  itk::ImageIOBase::Pointer imageIO;
  unsigned int imageIODimensionality;


  int lowIndexComponent = -1;
  int highIndexComponent = -1;
  int axisToExtractFrom = 0;
  bool axisToExtractFromSet = false;
  std::string *inputImage = new std::string("");
  std::string *outputImage = new std::string("");

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',
              CRKIT_VERSION_STRING);
      TCLAP::ValueArg<std::string> outputImageArg("o", "outputImage",
        "Output Image File Name", 
        true,"",
        "Output image file name", cmd);

      TCLAP::ValueArg<std::string> inputImageArg("i", "inputImage",
          "Input Image File Name",
          true,"",
          "Input image file name", cmd);

      TCLAP::ValueArg<int> highIndexComponentArg("m", "highIndexComponent",
          "Highest Component Index",
           true,-1,
           "integer value of highest component index",cmd);

      TCLAP::ValueArg<int> lowIndexComponentArg("l", "lowIndexComponent",
           "Lowest Component Index",
           true,-1,
           "integer value of lowest component index", cmd);
      TCLAP::ValueArg<int> axisToExtractFromArg("a","axis",
         "Axis to extract from",
          false,0,
          "integer value of axis to extract from", cmd);

      cmd.parse(argc,argv);
      if (inputImageArg.isSet()) inputImage = new std::string(inputImageArg.getValue());
      if (outputImageArg.isSet()) outputImage = new std::string(outputImageArg.getValue());
      if (lowIndexComponentArg.isSet()) lowIndexComponent = lowIndexComponentArg.getValue();
      if (highIndexComponentArg.isSet()) highIndexComponent = highIndexComponentArg.getValue();
      if (axisToExtractFromArg.isSet()) {
        axisToExtractFrom = axisToExtractFromArg.getValue();
        axisToExtractFromSet = true;
      };
    }
  catch (TCLAP::ArgException &e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(1);
    }

  if ( (lowIndexComponent == -1) && (highIndexComponent == -1) ) {
    std::cout << "Need a component index or range from zero up." << std::endl;
    exit(1);
  }
  if ( (lowIndexComponent == -1) && (highIndexComponent != -1) ||
       (lowIndexComponent != -1) && (highIndexComponent == -1) ) {
    if (lowIndexComponent == -1) lowIndexComponent = highIndexComponent;
    if (highIndexComponent == -1) highIndexComponent = lowIndexComponent;
  }

  if (lowIndexComponent > highIndexComponent) {
    int tmp = highIndexComponent;
    highIndexComponent = lowIndexComponent;
    lowIndexComponent = tmp;
  }

  crl::ExtractSmallerImageFromImageBase::Pointer 
     crlExtractSmallerImageFromImagePtr = 0;

  GetImageType( (*inputImage), pixelType, componentType,
                imageIO, imageIODimensionality);

  if (pixelType == itk::ImageIOBase::SCALAR) {
     INSTANTIATECLASSOVERDIMANDTYPE(componentType);
  }

  // In case the switch statement doesn't initialize the pointer
  if (!crlExtractSmallerImageFromImagePtr) {
    std::cout << argv[0] << " not implemented for pixel type " <<
        imageIO->GetPixelTypeAsString(imageIO->GetPixelType()) <<
         " with dimensions " << imageIO->GetNumberOfDimensions() <<
        std::endl;
    exit(1);
  }

  crlExtractSmallerImageFromImagePtr->SetInputImageFileName( (*inputImage) );
  crlExtractSmallerImageFromImagePtr->SetOutputImageFileName( (*outputImage) );
  // If set, use it, else default to last axis.
  if (axisToExtractFromSet) {
    crlExtractSmallerImageFromImagePtr->SetAxis( axisToExtractFrom );
  }
  crlExtractSmallerImageFromImagePtr->SetLowIndex( lowIndexComponent );
  crlExtractSmallerImageFromImagePtr->SetHighIndex( highIndexComponent );

  int ret = crlExtractSmallerImageFromImagePtr->Execute();
  if (ret != 0) {
    std::cerr << argv[0] << " return error code " << ret
              << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

