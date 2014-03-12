/* 
 * Neil Weisenfeld neil.weisenfeld@childrens.harvard.edu
 */

#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>

#include "crlGaussianSmooth.h"

// These macros require certain variables to have been initialized
// by calls to imageIOBase
#define INSTANTIATECLASS(dim, scalar) do { \
           if (imageIODimensionality == dim) { \
             crlGaussianSmoothPtr = \
       crl::GaussianSmooth< dim, scalar >::New(); \
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


int main( int argc, char * argv[] )
{
  if( argc != 4 && argc != 6) 
    { 
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << 
      "  inputImageFile outputImageFile sigma [clamp-lo clamp-hi]" << std::endl;

    return EXIT_FAILURE;
    }

  const char* infile = argv[1];
  const char* outfile = argv[2];
  const double sigma = atof( argv[3] );

  itk::ImageIOBase::IOPixelType pixelType;
  itk::ImageIOBase::IOComponentType componentType;
  itk::ImageIOBase::Pointer imageIO;
  unsigned int imageIODimensionality;

  crl::GaussianSmoothBase::Pointer crlGaussianSmoothPtr = 0;
  GetImageType( infile, pixelType, componentType,
                imageIO, imageIODimensionality);

  if (pixelType == itk::ImageIOBase::SCALAR) {
     INSTANTIATECLASSOVERDIMANDTYPE(componentType);
  }

  // In case the switch statement doesn't initialize the pointer
  if (!crlGaussianSmoothPtr) {
    std::cout << argv[0] << " not implemented for pixel type " <<
        imageIO->GetPixelTypeAsString(imageIO->GetPixelType()) <<
         " with dimensions " << imageIO->GetNumberOfDimensions() <<
        std::endl;
    exit(1);
  }

  std::cout << "smoothing " << infile << " --> " << outfile << std::endl
    << " with sigma=" << sigma << std::endl;

  crlGaussianSmoothPtr->SetNormalizeAcrossScale( false );
  crlGaussianSmoothPtr->SetInputImageFileName( infile );
  crlGaussianSmoothPtr->SetSigma( sigma );
  crlGaussianSmoothPtr->SetOutputImageFileName( outfile );

  if (argc == 6) {
    crlGaussianSmoothPtr->SetClamp( true );
    crlGaussianSmoothPtr->SetMinValue( atof(argv[4]) );
    crlGaussianSmoothPtr->SetMaxValue( atof(argv[5]) );
  }

  int ret = crlGaussianSmoothPtr->Execute();
  if (ret != 0) {
    std::cerr << argv[0] << " return error code " << ret
              << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

