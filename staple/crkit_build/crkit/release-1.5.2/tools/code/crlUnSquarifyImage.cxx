
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
  unsigned int xpad, ypad;

  try
    {
    TCLAP::CmdLine cmd("Remove padding added by crlSquarifyImage",' ',
             CRKIT_VERSION_STRING );
    TCLAP::UnlabeledValueArg<std::string> inputFileArg( "inputFile",
      "a valid file pathname", true, "", "input filename", cmd );
    TCLAP::UnlabeledValueArg<unsigned int> xpadArg( "xpadding", 
      "the first value returned by crlSquarifyImage", true, 0, "number voxels", cmd );
    TCLAP::UnlabeledValueArg<unsigned int> ypadArg( "ypadding", 
      "the first value returned by crlSquarifyImage", true, 0, "number voxels", cmd );
    TCLAP::UnlabeledValueArg<std::string> outFileArg( "outputFile",
      "a valid file pathname", true, "", "output filename", cmd );

    cmd.parse(argc,argv);

    inputFile = inputFileArg.getValue();
    xpad = xpadArg.getValue();
    ypad = ypadArg.getValue();
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
  newsize[0] = oldsize[0] - xpad;
  newsize[1] = oldsize[1] - ypad;

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
    iter.Value() = input->GetPixel( iter.GetIndex() );
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


  return 0;
}
