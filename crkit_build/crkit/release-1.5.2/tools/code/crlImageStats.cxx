
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkStatisticsImageFilter.h>

#include <iostream>
// Enable control of printed number precision
#include <iomanip>


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


int main(int argc, char *argv[])
{

  itk::ImageIOBase::IOPixelType pixelType;
  itk::ImageIOBase::IOComponentType componentType;
  itk::ImageIOBase::Pointer imageIO;
  unsigned int imageIODimensionality;

  if (argc < 2) {
    std::cout << 
      "Usage: crlImageStats inimage [inimage2] ..." << std::endl;
    exit(1);
  }

  GetImageType( argv[1], pixelType, componentType,
                imageIO, imageIODimensionality);

  static int const ImageDimension = 3;
  typedef   float PixelType;

  typedef   itk::Image<PixelType,ImageDimension>  ImageType;
  typedef   itk::ImageFileReader< ImageType >    ImageReaderType;
  typedef   itk::StatisticsImageFilter<ImageType> FilterType;


  ImageReaderType::Pointer r;
  FilterType::Pointer f;

  for (int i = 1; i < argc; i++) {

  r = ImageReaderType::New();
  f = FilterType::New();

  r->SetFileName(argv[i]);

  try {
    r->Update();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  f->SetInput( r->GetOutput() );
  f->Update();

  std::cout << setiosflags(std::ios::fixed) << std::setprecision(15);

  std::cout << "============ " << argv[i] << " ============" << std::endl;

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

  std::cout << "Geometric properties: " << std::endl;
  std::cout << "Size: " << size << std::endl;
  std::cout << "StartIndex: " << start << std::endl;
  std::cout << "Spacing: " << r->GetOutput()->GetSpacing() << std::endl;
  std::cout << "Origin: " << r->GetOutput()->GetOrigin() << std::endl;
  std::cout << "Direction: " << std::endl << r->GetOutput()->GetDirection() << std::endl;

  }

  exit(0); // success
}

