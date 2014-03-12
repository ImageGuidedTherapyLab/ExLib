/*
 * Copyright (c) 2008, Children's Hospital Boston.
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

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkSymmetricSecondRankTensor.h>
#include <itkDiffusionTensor3D.h>

#include <iostream>

#include "crlRootMeanSquare.h"

int GetImageType (std::string fileName,
                     itk::ImageIOBase::IOPixelType &pixelType,
                     itk::ImageIOBase::IOComponentType &componentType,
                     unsigned int &numberOfComponents,
                     itk::ImageIOBase::Pointer &imageIO,
                     unsigned int &dimensionality )
{
  typedef itk::Image<unsigned char, 2> ImageType;
  itk::ImageFileReader<ImageType>::Pointer imageReader =
                                itk::ImageFileReader<ImageType>::New();
  imageReader->SetFileName(fileName.c_str());
  try {
    imageReader->UpdateOutputInformation();
  } catch ( itk::ExceptionObject &err )
  {
    std::cout << "ExceptionObject caught!" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }

  pixelType = imageReader->GetImageIO()->GetPixelType();
  componentType = imageReader->GetImageIO()->GetComponentType();
  imageIO = imageReader->GetImageIO();
  numberOfComponents = imageIO->GetNumberOfComponents();
  dimensionality = imageIO->GetNumberOfDimensions();
  return EXIT_SUCCESS;
};

// These macros require certain variables to have been initialized
// by calls to imageIOBase
#define INSTANTIATECLASS(dim, scalar) do { \
           if (imageIODimensionality == dim) { \
             crlRootMeanSquarePtr = \
       crlRootMeanSquare< dim, scalar >::New(); \
           } \
         } while (0)

#define INSTANTIATECLASSOVERDIM(scalar) do { \
    INSTANTIATECLASS(2 , scalar); \
    INSTANTIATECLASS(3 , scalar); \
    INSTANTIATECLASS(4 , scalar); \
    INSTANTIATECLASS(5 , scalar); \
    INSTANTIATECLASS(6 , scalar); \
  } while (0)

// Root mean square calculation should involve a floating point type,
// so all types get promoted to type float or double:
#define INSTANTIATECLASSOVERDIMANDTYPE(componentType) do { \
  switch(componentType) { \
    case itk::ImageIOBase::UCHAR : \
    case itk::ImageIOBase::CHAR : \
    case itk::ImageIOBase::USHORT : \
    case itk::ImageIOBase::SHORT : \
    case itk::ImageIOBase::UINT :  \
    case itk::ImageIOBase::INT :   \
    case itk::ImageIOBase::ULONG : \
    case itk::ImageIOBase::LONG :  \
    case itk::ImageIOBase::FLOAT : \
         INSTANTIATECLASSOVERDIM(float); break; \
    case itk::ImageIOBase::DOUBLE : \
         INSTANTIATECLASSOVERDIM(double); break; \
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE : \
    default : break; \
  } \
  } while (0)

#define INSTANTIATETENSORCLASS(dim, componenttype) do { \
           if (imageIODimensionality == dim) { \
             crlRootMeanSquarePtr = \
       crlRootMeanSquare< dim, \
               itk::DiffusionTensor3D< componenttype > >::New(); \
           } \
         } while (0)

#define INSTANTIATETENSORCLASSOVERDIM(componenttype) do { \
    INSTANTIATETENSORCLASS(2 , componenttype); \
    INSTANTIATETENSORCLASS(3 , componenttype); \
    INSTANTIATETENSORCLASS(4 , componenttype); \
    INSTANTIATETENSORCLASS(5 , componenttype); \
    INSTANTIATETENSORCLASS(6 , componenttype); \
  } while (0)

// Root mean square calculation should involve a floating point type,
// so all types get promoted to type float or double:
#define INSTANTIATETENSORCLASSOVERDIMANDTYPE(componentType) do { \
  switch(componentType) { \
    case itk::ImageIOBase::UCHAR : \
    case itk::ImageIOBase::CHAR : \
    case itk::ImageIOBase::USHORT : \
    case itk::ImageIOBase::SHORT : \
    case itk::ImageIOBase::UINT :  \
    case itk::ImageIOBase::INT :   \
    case itk::ImageIOBase::ULONG : \
    case itk::ImageIOBase::LONG :  \
    case itk::ImageIOBase::FLOAT : \
         INSTANTIATETENSORCLASSOVERDIM(float); break; \
    case itk::ImageIOBase::DOUBLE : \
         INSTANTIATETENSORCLASSOVERDIM(double); break; \
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE : \
    default : break; \
  } \
  } while (0)

#define INSTANTIATEVECTORCLASS(spatialdim, componenttype, vectordim) do { \
           if (imageIODimensionality == spatialdim) { \
             if (numComponents == vectordim) { \
             crlRootMeanSquarePtr = \
       crlRootMeanSquare< spatialdim, \
               itk::Vector< componenttype, vectordim > >::New(); \
             } \
           } \
         } while (0)

#define INSTANTIATEVECTORCLASSOVERDIM(componenttype) do { \
    INSTANTIATEVECTORCLASS(2 , componenttype, 1); \
    INSTANTIATEVECTORCLASS(2 , componenttype, 2); \
    INSTANTIATEVECTORCLASS(2 , componenttype, 3); \
    INSTANTIATEVECTORCLASS(2 , componenttype, 4); \
    INSTANTIATEVECTORCLASS(2 , componenttype, 5); \
    INSTANTIATEVECTORCLASS(2 , componenttype, 6); \
    INSTANTIATEVECTORCLASS(2 , componenttype, 7); \
    INSTANTIATEVECTORCLASS(2 , componenttype, 8); \
    INSTANTIATEVECTORCLASS(2 , componenttype, 9); \
    INSTANTIATEVECTORCLASS(3 , componenttype, 1); \
    INSTANTIATEVECTORCLASS(3 , componenttype, 2); \
    INSTANTIATEVECTORCLASS(3 , componenttype, 3); \
    INSTANTIATEVECTORCLASS(3 , componenttype, 4); \
    INSTANTIATEVECTORCLASS(3 , componenttype, 5); \
    INSTANTIATEVECTORCLASS(3 , componenttype, 6); \
    INSTANTIATEVECTORCLASS(3 , componenttype, 7); \
    INSTANTIATEVECTORCLASS(3 , componenttype, 8); \
    INSTANTIATEVECTORCLASS(3 , componenttype, 9); \
    INSTANTIATEVECTORCLASS(4 , componenttype, 1); \
    INSTANTIATEVECTORCLASS(4 , componenttype, 2); \
    INSTANTIATEVECTORCLASS(4 , componenttype, 3); \
    INSTANTIATEVECTORCLASS(4 , componenttype, 4); \
    INSTANTIATEVECTORCLASS(4 , componenttype, 5); \
    INSTANTIATEVECTORCLASS(4 , componenttype, 6); \
    INSTANTIATEVECTORCLASS(4 , componenttype, 7); \
    INSTANTIATEVECTORCLASS(4 , componenttype, 8); \
    INSTANTIATEVECTORCLASS(4 , componenttype, 9); \
  } while (0)

// Root mean square calculation should involve a floating point type,
// so all types get promoted to type float or double:
#define INSTANTIATEVECTORCLASSOVERDIMANDTYPE(componentType) do { \
  switch(componentType) { \
    case itk::ImageIOBase::UCHAR : \
    case itk::ImageIOBase::CHAR : \
    case itk::ImageIOBase::USHORT : \
    case itk::ImageIOBase::SHORT : \
    case itk::ImageIOBase::UINT :  \
    case itk::ImageIOBase::INT :   \
    case itk::ImageIOBase::ULONG : \
    case itk::ImageIOBase::LONG :  \
    case itk::ImageIOBase::FLOAT : \
         INSTANTIATEVECTORCLASSOVERDIM(float); break; \
    case itk::ImageIOBase::DOUBLE : \
         INSTANTIATEVECTORCLASSOVERDIM(double); break; \
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE : \
    default : break; \
  } \
  } while (0)


int main(int argc, char *argv[])
{

  itk::ImageIOBase::IOPixelType pixelType;
  itk::ImageIOBase::IOComponentType componentType;
  unsigned int numComponents;
  itk::ImageIOBase::Pointer imageIO;
  unsigned int imageIODimensionality;

  if (argc != 3) {
    std::cout << 
      "Usage: crlRootMeanSquare image1 image2" << std::endl;
    exit(1);
  }

  crlRootMeanSquareBase::Pointer crlRootMeanSquarePtr = 0;

  GetImageType( argv[1], pixelType, componentType, numComponents,
                imageIO, imageIODimensionality);

  if (pixelType == itk::ImageIOBase::SCALAR) {
     INSTANTIATECLASSOVERDIMANDTYPE(componentType);
  } else if (pixelType == itk::ImageIOBase::VECTOR) {
     INSTANTIATEVECTORCLASSOVERDIMANDTYPE(componentType);
  } else if (
    (imageIO->GetPixelType() == itk::ImageIOBase::SYMMETRICSECONDRANKTENSOR) ||
    (imageIO->GetPixelType() == itk::ImageIOBase::DIFFUSIONTENSOR3D)
    ) {
     INSTANTIATETENSORCLASSOVERDIMANDTYPE(componentType);
  }

  // In case the switch statement doesn't initialize the pointer
  if (!crlRootMeanSquarePtr) {
    std::cout << "crlRootMeanSquare not implemented for pixel type " <<
        imageIO->GetPixelTypeAsString(imageIO->GetPixelType()) << 
         " with dimensions " << imageIO->GetNumberOfDimensions() <<
        std::endl;
    exit(1);
  }

  crlRootMeanSquarePtr->SetImage1FileName( argv[1] );
  crlRootMeanSquarePtr->SetImage2FileName( argv[2] );

  int ret = crlRootMeanSquarePtr->Execute();
  if (ret != 0) {
    std::cerr << "crlRootMeanSquare returned an error code " << ret << 
    std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Root mean square difference " << 
           crlRootMeanSquarePtr->GetRootMeanSquare() << std::endl;

  exit(0); // success
}
