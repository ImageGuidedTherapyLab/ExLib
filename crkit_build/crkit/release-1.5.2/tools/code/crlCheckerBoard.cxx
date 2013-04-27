
#include <sstream>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkIdentityTransform.h"

#include "itkResampleImageFilter.h"
#include "itkCheckerBoardImageFilter.h"

int main( int argc, char *argv[] )
{
  if ( argc != 4 && argc != 5 )
    {
      std::cerr << "Missing Parameters " << std::endl;
      std::cerr << "Usage: " << argv[0];
      std::cerr << " image1File image2File ";
      std::cerr << " outputImagefile [numchecks (def: 4)]";
      std::cerr << std::endl;
      return EXIT_FAILURE;
    }

  int nchecks = 4;

  if ( argc == 5 )
    {
    nchecks = atoi( argv[4] );
    }
  std::cerr << "using " << nchecks << " checks per dimension " << std::endl;
  
  const    unsigned int    Dimension = 3;
  typedef  float PixelType;
  
  typedef itk::Image< PixelType, Dimension >  ImageType;

  typedef itk::ImageFileReader< ImageType  > ImageReaderType;

  ImageReaderType::Pointer  imageReader1 = ImageReaderType::New();
  ImageReaderType::Pointer imageReader2 = ImageReaderType::New();

  imageReader1->SetFileName(  argv[1] );
  imageReader2->SetFileName( argv[2] );

  typedef itk::ResampleImageFilter< 
                            ImageType, 
                            ImageType >    ResampleFilterType;

  // By default this uses an identity transform.
  ResampleFilterType::Pointer resample = ResampleFilterType::New();
  resample->SetInput( imageReader2->GetOutput() );

  ImageType::Pointer image1 = imageReader1->GetOutput();
  imageReader1->Update();
  imageReader2->Update();

  PixelType backgroundGrayLevel = 0.0;

  resample->SetSize(    image1->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputDirection(  image1->GetDirection() );
  resample->SetOutputParametersFromImage( image1 );
  resample->SetDefaultPixelValue( backgroundGrayLevel );

  // Generate checkerboards before and after registration
  //
  typedef itk::CheckerBoardImageFilter< ImageType > CheckerBoardFilterType;

  CheckerBoardFilterType::Pointer checker = CheckerBoardFilterType::New();
  typedef itk::ImageFileWriter< ImageType > ImageWriterType;
  ImageWriterType::Pointer      writer =  ImageWriterType::New();

  CheckerBoardFilterType::PatternArrayType pat;
  pat.Fill( nchecks );

  checker->SetInput1( image1 );
  checker->SetInput2( resample->GetOutput() );
  checker->SetCheckerPattern( pat );
  writer->SetInput( checker->GetOutput()   );
  writer->SetFileName( argv[3] );
  writer->Update();
  
  return EXIT_SUCCESS;
}

