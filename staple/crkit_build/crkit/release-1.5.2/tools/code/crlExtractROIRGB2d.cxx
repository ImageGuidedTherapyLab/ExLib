
/* Facilitate the selection of a region of interest from the data. */

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImage.h"
#include "itkRGBPixel.h"

int main( int argc, char ** argv )
{
  // Verify the number of parameters in the command line
  if( argc != 7 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImageFile  outputImageFile " << std::endl;
    std::cerr << " startX startY sizeX sizeY " << std::endl;
    return EXIT_FAILURE;
    }

  // typedef unsigned char InputPixelType;
  // typedef unsigned char OutputPixelType;
  typedef itk::RGBPixel< unsigned char > InputPixelType;
  typedef itk::RGBPixel< unsigned char > OutputPixelType;
  const   unsigned int        Dimension = 2;

  typedef itk::Image< InputPixelType,  Dimension >    InputImageType;
  typedef itk::Image< OutputPixelType, Dimension >    OutputImageType;


  typedef itk::ImageFileReader< InputImageType  >  ReaderType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  typedef itk::RegionOfInterestImageFilter< InputImageType, 
                                            OutputImageType > FilterType;

  FilterType::Pointer filter = FilterType::New();

  OutputImageType::IndexType start;
  start[0] = atoi( argv[3] );
  start[1] = atoi( argv[4] );

  OutputImageType::SizeType size;
  size[0] = atoi( argv[5] );
  size[1] = atoi( argv[6] );

  OutputImageType::RegionType desiredRegion;
  desiredRegion.SetSize(  size  );
  desiredRegion.SetIndex( start );

  //  Then the region is passed to the filter using SetRegionOfInterest()
  filter->SetRegionOfInterest( desiredRegion );

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  const char * inputFilename  = argv[1];
  const char * outputFilename = argv[2];


  reader->SetFileName( inputFilename  );
  writer->SetFileName( outputFilename );

  filter->SetInput( reader->GetOutput() );
  writer->SetInput( filter->GetOutput() );

  try 
    { 
    writer->Update(); 
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    } 


  return EXIT_SUCCESS;
}

