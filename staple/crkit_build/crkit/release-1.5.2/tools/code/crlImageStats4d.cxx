
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkStatisticsImageFilter.h>

#include <iostream>
// Enable control of printed number precision
#include <iomanip>


int main(int argc, char *argv[])
{

  static int const ImageDimension = 4;
  typedef   float PixelType;

  typedef   itk::Image<PixelType,ImageDimension>  ImageType;
  typedef   itk::ImageFileReader< ImageType >    ImageReaderType;
  typedef   itk::StatisticsImageFilter<ImageType> FilterType;

  if (argc != 2) {
    std::cout << 
      "Usage: crlImageStats inimage" << std::endl;
    exit(1);
  }

  ImageReaderType::Pointer r;
  FilterType::Pointer f;

  r = ImageReaderType::New();
  f = FilterType::New();

  r->SetFileName(argv[1]);

  try {
    r->Update();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  f->SetInput( r->GetOutput() );
  f->Update();

  std::cout << setiosflags(std::ios::fixed) << std::setprecision(15);

  std::cout << "Intensity properties: " << std::endl;
  std::cout << "Minimum " << f->GetMinimum();
  std::cout << " Maximum " << f->GetMaximum();
  std::cout << " Mean " << f->GetMean();
  std::cout << " Variance " << f->GetVariance();
  std::cout << " Sum " << f->GetSum();
  std::cout << std::endl;

  ImageType::RegionType inputRegion =
                   r->GetOutput()->GetLargestPossibleRegion();
  ImageType::SizeType size = inputRegion.GetSize();
  ImageType::IndexType start = inputRegion.GetIndex();
  ImageType::SpacingType spacing = r->GetOutput()->GetSpacing();
  

  std::cout << "Geometric properties: " << std::endl;
  std::cout << "size: " << size << std::endl;
  std::cout << "start: " << start << std::endl;
  std::cout << "spacing: " << r->GetOutput()->GetSpacing() << std::endl;
  std::cout << "origin: " << r->GetOutput()->GetOrigin() << std::endl;

  exit(0); // success
}

