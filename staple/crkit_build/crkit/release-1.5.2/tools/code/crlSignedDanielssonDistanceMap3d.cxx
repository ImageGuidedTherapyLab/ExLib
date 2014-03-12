#include "itkSignedDanielssonDistanceMapImageFilter.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImageFile outputDistanceMapImageFile ";
    std::cerr << " outputVoronoiMapImageFilter ";
    std::cerr << " outputVectorMapImageFilter ";
    std::cerr << std::endl;  
    return EXIT_FAILURE;
    }

  // Software Guide : BeginCodeSnippet
  typedef  float InputPixelType;
  typedef  float  OutputPixelType;
  const unsigned int Dimension = 3;

  typedef itk::Image< InputPixelType,  Dimension >   InputImageType;
  typedef itk::Image< OutputPixelType, Dimension >   OutputImageType;

  typedef itk::SignedDanielssonDistanceMapImageFilter<
                                         InputImageType, 
                                         OutputImageType >  FilterType;
  
  FilterType::Pointer filter = FilterType::New();

  // Reader and Writer types are instantiated.
  //
  typedef itk::ImageFileReader< InputImageType  >  ReaderType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );

  filter->SetInput( reader->GetOutput() );
  writer->SetInput( filter->GetOutput() );

  writer->Update();
  
  const char * voronoiMapFileName = argv[3];

  //  The Voronoi map is obtained with the \code{GetVoronoiMap()} method. 
  //
  writer->SetInput( filter->GetVoronoiMap() );
  writer->SetFileName( voronoiMapFileName );
  writer->Update();

  //  The distance filter also produces an image of \doxygen{Offset} pixels
  //  representing the vectorial distance to the closest object in the scene.
  //  The type of this output image is defined by the VectorImageType
  //  trait of the filter type.

  typedef FilterType::VectorImageType   OffsetImageType;

  //  We can use this type for instantiating an \doxygen{ImageFileWriter} type
  //  and creating an object of this class in the following lines.

  typedef itk::ImageFileWriter< OffsetImageType >  WriterOffsetType;
  WriterOffsetType::Pointer offsetWriter = WriterOffsetType::New();
  offsetWriter->SetInput(  filter->GetVectorDistanceMap()  );
  offsetWriter->SetFileName( argv[4]  );

  try
    {
    offsetWriter->Update();
    }
  catch( itk::ExceptionObject exp )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr <<     exp    << std::endl;
    }

  return EXIT_SUCCESS;
}

