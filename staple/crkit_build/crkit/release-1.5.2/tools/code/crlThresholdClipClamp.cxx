#include <itkThresholdImageFilter.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

int main( int argc, char * argv[] )
{

  if ( argc != 5 )
    {
      std::cerr << "Usage: " << argv[0] << " inputImage outputImage " <<
        " thresholdBelow thresholdAbove " << std::endl;
      return EXIT_FAILURE;
    }
  
  static int const ImageDimension = 3;
  typedef  float PixelType;
  typedef itk::Image< PixelType,  ImageDimension >   ImageType;

  typedef itk::ThresholdImageFilter< ImageType >  FilterType;

  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  FilterType::Pointer filter = FilterType::New();

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  reader->SetFileName( argv[1] );

  filter->SetInput( reader->GetOutput() );
  filter->SetOutsideValue( 0 );

  filter->ThresholdOutside( atof(argv[3]), atof(argv[4]) );
  filter->Update();

  writer->SetFileName( argv[2] );
  writer->Update();

  return EXIT_SUCCESS;
}

