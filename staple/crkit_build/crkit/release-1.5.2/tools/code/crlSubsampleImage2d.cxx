
#include "itkOrientedImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkShrinkImageFilter.h"


int main( int argc, char * argv[] )
{
  static unsigned int const ImageDimension = 2;
  typedef   float PixelType;

  if( argc != (3+ImageDimension) ) 
    { 
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << 
      " inputImageFile outputImageFile subsamplex subsampley ..." << std::endl;
    return EXIT_FAILURE;
    }

  const char* infile = argv[1];
  const char* outfile = argv[2];
  std::vector<int> shrinkfactors(ImageDimension);
  for (unsigned int i = 0; i < ImageDimension; i++) {
    shrinkfactors[i] = atoi(argv[3 + i]);
  }

  typedef itk::OrientedImage< PixelType, ImageDimension >   ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;
  typedef itk::ShrinkImageFilter<
                        ImageType, ImageType >  FilterType;

  try 
    {
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( infile );
    reader->Update();

    FilterType::Pointer filter = FilterType::New();
    filter->SetInput( reader->GetOutput() );
    for (unsigned int i = 0; i < ImageDimension; i++) {
      filter->SetShrinkFactor(i, shrinkfactors[i]);
    }

    filter->Update();

    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outfile );
    writer->SetInput( filter->GetOutput() );
    writer->Update();
    } 
  catch ( itk::ExceptionObject& err )
    {
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

