
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>

int main( int argc, char * argv[] )
{
  if( argc != 4 ) 
    { 
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  infile radius outputImageFile " << std::endl;

    return EXIT_FAILURE;
    }

  const char* infile = argv[1];
  double radius = atof(argv[2]);
  const char* outfile = argv[3];

  const unsigned int Dimension = 2;
  typedef float PixelType;
  typedef itk::Image< PixelType, Dimension >   ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;
  typedef itk::ImageRegionIteratorWithIndex< ImageType >  IteratorType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  ImageType::Pointer image = ImageType::New();

  try 
    {
    reader->SetFileName( infile );
    reader->Update();
    } 
  catch ( itk::ExceptionObject& err )
    {
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

    image = reader->GetOutput();

    IteratorType it( image,
                     image->GetLargestPossibleRegion() );

    for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
      ImageType::IndexType ind = it.GetIndex();
      double dist = 0.0;
      for (unsigned int i = 0; i < Dimension; i++) {
        dist += 
                (ind[i] - image->GetLargestPossibleRegion().GetSize()[i]/2.0)*
                (ind[i] - image->GetLargestPossibleRegion().GetSize()[i]/2.0);
      }
      if (sqrt(dist) < radius) {
        it.Set(1.0);
      } else {
        it.Set(0.0);
      }
    }

    try {
    writer->SetFileName( outfile );
    writer->SetInput( image );
    writer->Update();
    } 
  catch ( itk::ExceptionObject& err )
    {
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

