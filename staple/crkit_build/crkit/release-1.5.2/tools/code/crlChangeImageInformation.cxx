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

#include "crlChangeImageInformation.h"

// These macros require certain variables to have been initialized 
// by calls to imageIOBase
#define INSTANTIATECLASS(dim, scalar) do { \
           if (imageIODimensionality == dim) { \
             crlChangeImageInformationPtr = \
       crl::ChangeImageInformation< dim, scalar >::New(); \
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

  std::string *inputImage = new std::string("");
  std::string *referenceImage = new std::string("");
  bool useReferenceImage = false;
  std::string *outputImage = new std::string("");
  bool m_ChangeSpacing = false;
  bool m_ChangeOrigin = false;
  bool m_ChangeDirection = false;
  bool m_ChangeRegion = false;

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',
              CRKIT_VERSION_STRING);

      TCLAP::ValueArg<std::string> outputImageArg("o", "outputImage",
        "Output Image File Name", 
        true,"",
        "Output image file name", cmd);

      TCLAP::SwitchArg changeRegionArg("c", "changeRegion", 
        "Change region", cmd, false);
      TCLAP::SwitchArg changeDirectionArg("d", "changeDirection", 
        "Change direction", cmd, false );
      TCLAP::SwitchArg changeSpacingArg("s", "changeSpacing", 
        "Change spacing", cmd, false);
      TCLAP::SwitchArg changeOriginArg("x", "changeOrigin", 
        "Change origin", cmd, false);

      TCLAP::ValueArg<std::string> referenceImageArg("r", "referenceImage",
        "Reference Image File Name", 
        false,"", "Reference image file name", cmd);

      TCLAP::ValueArg<std::string> inputImageArg("i", "inputImage",
          "Input Image File Name",
          true,"",
          "Input image file name", cmd);

      cmd.parse(argc,argv);
      if (inputImageArg.isSet()) inputImage = new std::string(inputImageArg.getValue());
      if (referenceImageArg.isSet()) {
        referenceImage = new std::string(referenceImageArg.getValue());
        useReferenceImage = true;
        std::cout << "Using reference image " << (*referenceImage) << std::endl;
      }
      if (changeSpacingArg.isSet()) { m_ChangeSpacing = true; }
      if (changeOriginArg.isSet()) { m_ChangeOrigin = true; }
      if (changeDirectionArg.isSet()) { m_ChangeDirection = true; }
      if (changeRegionArg.isSet()) { m_ChangeRegion = true; }

      if (outputImageArg.isSet()) outputImage = new std::string(outputImageArg.getValue());
    }
  catch (TCLAP::ArgException &e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(1);
    }

  crl::ChangeImageInformationBase::Pointer crlChangeImageInformationPtr = 0;

  GetImageType( (*inputImage), pixelType, componentType,
                imageIO, imageIODimensionality);

  if (pixelType == itk::ImageIOBase::SCALAR) {
     INSTANTIATECLASSOVERDIMANDTYPE(componentType);
  }

  // In case the switch statement doesn't initialize the pointer
  if (!crlChangeImageInformationPtr) {
    std::cout << argv[0] << " not implemented for pixel type " <<
        imageIO->GetPixelTypeAsString(imageIO->GetPixelType()) <<
         " with dimensions " << imageIO->GetNumberOfDimensions() <<
        std::endl;
    exit(1);
  }

  crlChangeImageInformationPtr->SetInputImageFileName( (*inputImage) );
  std::cout << "Set input file name to be " << (*inputImage) << std::endl;

  if (useReferenceImage) {
    crlChangeImageInformationPtr->SetReferenceImageFileName((*referenceImage));
  }
  crlChangeImageInformationPtr->SetOutputImageFileName( (*outputImage) );

  crlChangeImageInformationPtr->SetChangeSpacing( m_ChangeSpacing );
  crlChangeImageInformationPtr->SetChangeOrigin( m_ChangeOrigin );
  crlChangeImageInformationPtr->SetChangeDirection( m_ChangeDirection );
  crlChangeImageInformationPtr->SetChangeRegion( m_ChangeRegion );

  int ret = crlChangeImageInformationPtr->Execute();
  if (ret != 0) {
    std::cerr << argv[0] << " return error code " << ret
              << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

