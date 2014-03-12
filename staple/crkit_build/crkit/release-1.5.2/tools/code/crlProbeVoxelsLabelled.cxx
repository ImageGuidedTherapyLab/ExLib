
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIteratorWithIndex.h>

#include <iostream>
// Enable control of printed number precision
#include <iomanip>

int main(int argc, char *argv[])
{

  static int const ImageDimension = 3;
  typedef   float PixelType;

  typedef   itk::Image<PixelType,ImageDimension>  ImageType;
  typedef   itk::ImageFileReader< ImageType >    ImageReaderType;

  typedef ImageType::RegionType RegionType;
  typedef ImageType::SizeType   SizeType;
  typedef ImageType::IndexType  IndexType;

  if (argc != 4) {
    std::cout << 
      "Usage: crlProbeVoxelsLabelled inimage labelimage labelvalue" << std::endl;
    std::cout << "Report the values of all of the voxel corresponding to "
              << "the label value in the labelimage." << std::endl;
    exit(1);
  }

  PixelType labelvalue = ::atof( argv[3] );

  ImageReaderType::Pointer reader;
  ImageReaderType::Pointer lreader;

  reader = ImageReaderType::New();
  lreader = ImageReaderType::New();
  reader->SetFileName(argv[1]);
  lreader->SetFileName(argv[2]);

  try {
    reader->Update(); // read in the file, catch errors ...
  } catch ( itk::ExceptionObject &err )
  {
    std::cout << "ExceptionObject caught!" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }

  try {
    lreader->Update(); // read in the file, catch errors ...
  } catch ( itk::ExceptionObject &err )
  {
    std::cout << "ExceptionObject caught!" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }

  ImageType::Pointer image = reader->GetOutput();
  ImageType::Pointer limage = lreader->GetOutput();
  PixelType p1 = static_cast<PixelType>( 0 );
  PixelType p2 = static_cast<PixelType>( 0 );

  std::cout << setiosflags(std::ios::fixed) << std::setprecision(15);

  typedef   itk::ImageRegionConstIteratorWithIndex<ImageType> ImageIterator;
  ImageIterator iter1( lreader->GetOutput(),  
                       lreader->GetOutput()->GetLargestPossibleRegion() );
  std::cout << "Voxel values for label : " << labelvalue << std::endl;
  for (iter1.GoToBegin(); !iter1.IsAtEnd(); ++iter1) {
    ImageType::IndexType index = iter1.GetIndex();
    p1 = limage->GetPixel(index);
    if (p1 == labelvalue) {
      p2 = image->GetPixel(index);
      std::cout << p2 << std::endl;
    }
  }

  exit(0); // success
}

