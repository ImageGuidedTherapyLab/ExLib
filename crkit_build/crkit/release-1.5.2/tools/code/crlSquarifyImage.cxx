
#include <itkOrientedImage.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkImageRegionIteratorWithIndex.h>
#include "tclap/CmdLine.h"
#include "vnl/vnl_math.h"

#include "configuration.h"


int main( int argc, char *argv[] )
{
  std::string inputFile;
  std::string outFile;

  try
    {
    TCLAP::CmdLine cmd("Pad image to make first two dimensions of image square",' ', CRKIT_VERSION_STRING );
    TCLAP::UnlabeledValueArg<std::string> inputFileArg( "inputFile",
      "a valid file pathname", true, "", "input filename", cmd );
    TCLAP::UnlabeledValueArg<std::string> outFileArg( "outputFile",
      "a valid file pathname", true, "", "output filename", cmd );

    cmd.parse(argc,argv);

    inputFile = inputFileArg.getValue();
    outFile = outFileArg.getValue();
    }
  catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(EXIT_FAILURE);
    }

  typedef itk::OrientedImage<float,3> ImageType;
  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  typedef itk::ImageFileWriter< ImageType > ImageWriterType;

  ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName( inputFile );
  try 
    {
    reader->Update();
    }
  catch ( itk::ExceptionObject& err )
    {
    std::cerr << "error reading input: " << inputFile << ": " << err <<
      std::endl;
    return 1;
    }

  ImageType::SizeType oldsize;
  ImageType::SizeType newsize;

  ImageType::Pointer input = reader->GetOutput();
  ImageType::RegionType oldRegion = input->GetBufferedRegion();

  oldsize = oldRegion.GetSize();

  newsize[2] = oldsize[2];
  newsize[0] = vnl_math_max( oldsize[0], oldsize[1] );
  newsize[1] = vnl_math_max( oldsize[0], oldsize[1] );

  // allocate new image
  ImageType::IndexType start;
  start = oldRegion.GetIndex();

  ImageType::RegionType newRegion;
  newRegion.SetIndex( start );
  newRegion.SetSize( newsize );


  ImageType::Pointer output = ImageType::New();
  output->CopyInformation( input );
  output->SetRegions( newRegion );
  output->Allocate();


  // populate new image
  itk::ImageRegionIteratorWithIndex< ImageType > iter( output, output->GetBufferedRegion() );

  for ( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    if ( oldRegion.IsInside( iter.GetIndex() ) )
      {
      iter.Value() = input->GetPixel( iter.GetIndex() );
      }
    else
      {
      iter.Value() = 0;
      }
    }


  // write new image
  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput( output );
  writer->SetFileName( outFile );
  writer->UseCompressionOn();
  try 
    {
    writer->Update();
    }
  catch (itk::ExceptionObject& err )
    {
    std::cerr << "error writing output file: " << outFile << ": " << err
      << std::endl;
    return 1;
    }

  // echo the pad factors
  std::cout << newsize[0] - oldsize[0] << " " << newsize[1] - oldsize[1] << std::endl;

  return 0;
}
