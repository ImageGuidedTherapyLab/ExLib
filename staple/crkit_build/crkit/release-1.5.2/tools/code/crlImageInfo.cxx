#include <stdlib.h>
#include "itkExceptionObject.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkOrientedImage.h"

int main( int argc, char* argv[] )
{
  if ( argc < 2 )
    {
    std::cerr << "usage: image-file [image-file...]" << std::endl;
    return 1;
    }


  for (int i = 1; i < argc; i++ )
    {
    typedef itk::OrientedImage<float,3> ImageType;
    typedef itk::ImageFileReader< ImageType > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(argv[i]);
    typedef itk::ExceptionObject ExceptionType;
    try
      {
      reader->Update();
      }
    catch (ExceptionType& err )
      {
      std::cerr << "Exception thrown: " << err << std::endl;
      return 1;
      }

    ImageType::Pointer image = reader->GetOutput();

    std::cout << "============ " << argv[i] << " ============" <<
      std::endl;

    ImageType::RegionType region = image->GetLargestPossibleRegion();
    ImageType::PointType origin = image->GetOrigin();
    ImageType::SpacingType spacing = image->GetSpacing();
    ImageType::DirectionType dir = image->GetDirection();
    ImageType::IndexType index = region.GetIndex();
    ImageType::SizeType size = region.GetSize();

    std::cout <<
      "Size: " << size << std::endl << 
      "Index: " << index << std::endl <<
      "Spacing: " << spacing << std::endl <<
      "Origin: " << origin << std::endl <<
      "Direction: " << std::endl << dir << std::endl;
    }

  return 0;
}
