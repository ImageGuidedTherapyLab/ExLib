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

#include "crlImageAlgebra.h"

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
             crlImageAlgebraPtr = \
       crlImageAlgebra< dim, scalar >::New(); \
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
             crlImageAlgebraPtr = \
       crlImageAlgebra< dim, \
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
             crlImageAlgebraPtr = \
       crlImageAlgebra< spatialdim, \
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

  if (argc != 5) {
    std::cout << 
      "Usage: " << argv[0] << " image1 operator image2 outimage" << std::endl;
    std::cout << "Computes : outimage = image1 operator image2 " << std::endl;
    std::cout << "Valid operators are : " << "add " << "subtract " 
              << "multiply " << "divide " << "difference " 
              << "squareddifference " << std::endl;
    return EXIT_FAILURE;
  }

  crlImageAlgebraBase::Pointer crlImageAlgebraPtr = 0;

  GetImageType( argv[1], pixelType, componentType, numComponents,
                imageIO, imageIODimensionality);

  if (pixelType == itk::ImageIOBase::SCALAR) {
     INSTANTIATECLASSOVERDIMANDTYPE(componentType);
  } else if ( (pixelType == itk::ImageIOBase::VECTOR) ||
              (pixelType == itk::ImageIOBase::RGB) ||
              (pixelType == itk::ImageIOBase::RGBA) ||
              (pixelType == itk::ImageIOBase::COVARIANTVECTOR) ||
              (pixelType == itk::ImageIOBase::POINT)
    ) {
     INSTANTIATEVECTORCLASSOVERDIMANDTYPE(componentType);
  } else if (
    (imageIO->GetPixelType() == itk::ImageIOBase::SYMMETRICSECONDRANKTENSOR) ||
    (imageIO->GetPixelType() == itk::ImageIOBase::DIFFUSIONTENSOR3D)
    ) {
     ; // INSTANTIATETENSORCLASSOVERDIMANDTYPE(componentType);
  }

  // In case the above statement doesn't initialize the pointer
  if (!crlImageAlgebraPtr) {
    std::cout << argv[0] << " not implemented for pixel type " <<
        imageIO->GetPixelTypeAsString(imageIO->GetPixelType()) <<
         " with dimensions " << imageIO->GetNumberOfDimensions() <<
        std::endl;
    exit(1);
  }

  crlImageAlgebraPtr->SetInput1FileName(  argv[1]);
  crlImageAlgebraPtr->SetOperationName(argv[2]);
  crlImageAlgebraPtr->SetInput2FileName(  argv[3]);
  crlImageAlgebraPtr->SetOutputFileName(  argv[4]);

  int ret = crlImageAlgebraPtr->Execute();
  if (ret != 0) {
    std::cerr << argv[0] << " returned an error code "
              << ret << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
