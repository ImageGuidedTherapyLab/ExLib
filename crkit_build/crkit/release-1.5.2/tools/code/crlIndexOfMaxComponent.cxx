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

#include "crlIndexOfMaxComponent.h"


// These macros require certain variables to have been initialized
// by calls to imageIOBase
#define INSTANTIATECLASS(dim, scalar) do { \
           if (imageIODimensionality == dim) { \
             crlIndexOfMaxComponentPtr = \
       crl::IndexOfMaxComponent< dim, scalar >::New(); \
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

  if (argc != 3) {
    std::cout << "Usage: " << argv[0] << " inimage outimage " << std::endl;
    std::cout << "Create a new image where each voxel has the " <<
        "value of the index of the component with the largest value." 
        << std::endl;
    return EXIT_FAILURE;
  }

  crl::IndexOfMaxComponentBase::Pointer crlIndexOfMaxComponentPtr = 0;

  GetImageType( argv[1], pixelType, componentType, 
                imageIO, imageIODimensionality);

  if (pixelType == itk::ImageIOBase::SCALAR) {
     INSTANTIATECLASSOVERDIMANDTYPE(componentType);
  }

  // In case the switch statement doesn't initialize the pointer
  if (!crlIndexOfMaxComponentPtr) {
    std::cout << "crlIndexOfMaxComponent not implemented for pixel type " <<
        imageIO->GetPixelTypeAsString(imageIO->GetPixelType()) <<
         " with dimensions " << imageIO->GetNumberOfDimensions() <<
        std::endl;
    exit(1);
  }

  crlIndexOfMaxComponentPtr->SetInputFileName( argv[1] );
  crlIndexOfMaxComponentPtr->SetOutputFileName( argv[2] );

  int ret = crlIndexOfMaxComponentPtr->Execute();
  if (ret != 0) {
    std::cerr << "crlIndexOfMaxComponent return error code " << ret 
              << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
