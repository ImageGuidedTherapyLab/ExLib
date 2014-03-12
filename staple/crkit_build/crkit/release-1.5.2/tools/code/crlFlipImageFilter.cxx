
#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkFlipImageFilter.h>

int main( int argc, char * argv[] )
{
  if( argc < 6 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImageFile outputImageFile flipAxisX flipAxisY flipAxisZ" << std::endl;
    return EXIT_FAILURE;
    }

  typedef   float PixelType;
  static int const Dimension = 3;

  typedef itk::OrientedImage< PixelType,  Dimension >   ImageType;

  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );

  typedef itk::FlipImageFilter< ImageType >  FilterType;
  FilterType::Pointer filter = FilterType::New();

  typedef FilterType::FlipAxesArrayType     FlipAxesArrayType;
  
  FlipAxesArrayType flipArray;

  flipArray[0] = atoi( argv[3] );
  flipArray[1] = atoi( argv[4] );
  flipArray[2] = atoi( argv[5] );
  filter->SetFlipAxes( flipArray );

  filter->SetInput( reader->GetOutput() );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}
