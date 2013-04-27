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
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#include "crlSulciDetection.h"


int GetImageType (std::string fileName,
                     itk::ImageIOBase::IOPixelType &pixelType,
                     itk::ImageIOBase::IOComponentType &componentType,
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
  dimensionality = imageIO->GetNumberOfDimensions();
  return EXIT_SUCCESS;
};

// These macros require certain variables to have been initialized
// by calls to imageIOBase
#define INSTANTIATECLASS(dim, scalar) do { \
           if (imageIODimensionality == dim) { \
             crlImageLaplacePtr = \
       crlImageLaplace< dim, scalar >::New(); \
           } \
         } while (0)

#define INSTANTIATECLASSOVERDIM(scalar) do { \
    INSTANTIATECLASS(2 , scalar); \
    INSTANTIATECLASS(3 , scalar); \
/*    INSTANTIATECLASS(4 , scalar); \
    INSTANTIATECLASS(5 , scalar); \
    INSTANTIATECLASS(6 , scalar); */ \
  } while (0)

// Laplace calculation requires a floating point type for the output image,
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

// We can't calculate a laplacian with a tensor or vector component image.

int main(int argc, char *argv[])
{

  itk::ImageIOBase::IOPixelType         pixelType;
  itk::ImageIOBase::IOComponentType componentType;
  itk::ImageIOBase::Pointer             imageIO;
  unsigned int imageIODimensionality;

  if (argc != 15) {
    std::cout << 
      "Usage: " << argv[0] << " labelimage innerlabel outerlabel"
       << " solveregionlabel numIniticalCortialLayers excludeSulci checkTopology expandSulci "
       << " maxIterations outputSegmentation outputLaplaceImage outputTangentField" 
       << " thicknessImageFileName displacementFieldFileName " << std::endl;
    return EXIT_FAILURE;
  }
  
  std::cout << " crlSulciDetection v1.0" << std::endl;

  crlImageLaplaceBase::Pointer crlImageLaplacePtr = 0;

  GetImageType( argv[1], pixelType, componentType,
                imageIO, imageIODimensionality);

  if (pixelType == itk::ImageIOBase::SCALAR) {
     INSTANTIATECLASSOVERDIMANDTYPE(componentType);
  } 
  else if (
    (imageIO->GetPixelType() == itk::ImageIOBase::SYMMETRICSECONDRANKTENSOR) ||
    (imageIO->GetPixelType() == itk::ImageIOBase::DIFFUSIONTENSOR3D)
    ) {
     // This type of input can't be used with this class.
     // INSTANTIATETENSORCLASSOVERDIMANDTYPE(componentType);
     crlImageLaplacePtr = 0;
  }

  // In case the switch statement doesn't initialize the pointer
  if (!crlImageLaplacePtr) {
    std::cout << "crlImageLaplace not implemented for pixel type " <<
        imageIO->GetPixelTypeAsString(imageIO->GetPixelType()) << 
         " with dimensions " << imageIO->GetNumberOfDimensions() <<
        std::endl;
    exit(1);
  }

  crlImageLaplacePtr->SetLabelImage( argv[1] );
  crlImageLaplacePtr->SetInnerBorderLabel( std::atoi(argv[2]) );
  crlImageLaplacePtr->SetOuterBorderLabel( std::atoi(argv[3]) );
  crlImageLaplacePtr->SetSolveRegionLabel( std::atoi(argv[4]) );
  crlImageLaplacePtr->SetInitialCorticalLayers( std::atoi( argv[5]) );
  crlImageLaplacePtr->SetExcludeSulci(std::atoi( argv[6]) );
  crlImageLaplacePtr->SetCheckTopology(std::atoi( argv[7]) );
  crlImageLaplacePtr->SetExpandSulci(std::atoi( argv[8]));
  crlImageLaplacePtr->SetMaxIterations(std::atoi (argv[9]));
  crlImageLaplacePtr->SetOutputSegFileName( argv[10] );
  crlImageLaplacePtr->SetOutputLaplaceFileName( argv[11] );
  crlImageLaplacePtr->SetOutputTangentFieldFileName( argv[12] );
  crlImageLaplacePtr->SetOutputThicknessFileName( argv[13] );
  crlImageLaplacePtr->SetOutputDisplacementFieldFileName( argv[14] );
  crlImageLaplacePtr->SetEList();

  int ret = crlImageLaplacePtr->Execute();
  if (ret != 0) {
    std::cerr << argv[0] << " returned an error code "
              << ret << " from Execute()" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
