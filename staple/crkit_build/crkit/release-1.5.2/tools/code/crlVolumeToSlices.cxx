
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageSeriesWriter.h>
#include <itkNumericSeriesFileNames.h>
// #include <itkGE5ImageIO.h>
#include <itkRawImageIO.h>

#include <tclap/CmdLine.h>

#include "configuration.h"

int main( int argc, char *argv[] )
{
  std::string *inputFile = new std::string("");
  std::string *outputPrefix = new std::string("");
  std::string *outputExtension = new std::string("");

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',
                   CRKIT_VERSION_STRING );
      TCLAP::UnlabeledValueArg<std::string> inputFileArg("inputFile","Input File Name",true,"","input file name",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputPrefixArg("outputPrefix","Output File Prefix",true,"","output file prefix",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputExtensionArg("outputExtension","Output File Extension",true,"","output file extension",cmd);

      cmd.parse(argc,argv);

      if (inputFileArg.isSet()) inputFile = new std::string(inputFileArg.getValue());
      if (outputPrefixArg.isSet()) outputPrefix = new std::string(outputPrefixArg.getValue());
      if (outputExtensionArg.isSet()) outputExtension = new std::string(outputExtensionArg.getValue());
    }
  catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(EXIT_FAILURE);
    }


//  The type of the input image is declared here and it is used for declaring
//  the type of the reader. This will be a conventional 3D image reader.
  typedef unsigned short PixelType;
  typedef itk::Image< PixelType, 3 >      ImageType;
  typedef itk::ImageFileReader< ImageType >   ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputFile->c_str() );

//  The type of the series writer must be instantiated taking into account that
//  the input file is a 3D volume and the output files are 2D images.
//  Additionally, the output of the reader is connected as input to the writer.

  typedef itk::Image< PixelType, 2 >     Image2DType;
  typedef itk::RawImageIO< PixelType, 2> ImageIOType;

  typedef itk::ImageSeriesWriter< ImageType, Image2DType > WriterType;

  ImageIOType::Pointer specialio = ImageIOType::New();
  WriterType::Pointer writer = WriterType::New();

  writer->SetInput( reader->GetOutput() );

  typedef itk::NumericSeriesFileNames    NameGeneratorType;
  NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();

  bool useSpecialImageIO = false;
  if (*outputExtension == "raw") {
    useSpecialImageIO = true;
  }
  std::string format = *outputPrefix;
  if (useSpecialImageIO) { // Use raw IO
    format += ".%03d";
    writer->SetImageIO( specialio );
  } else {
    format += "%03d.";
    format += *outputExtension;   // filename extension
  }

  nameGenerator->SetSeriesFormat( format.c_str() );

  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while reading the image" << std::endl;
    std::cerr << excp << std::endl;
    }

  ImageType::ConstPointer inputImage = reader->GetOutput();
  ImageType::RegionType   region     = inputImage->GetLargestPossibleRegion();
  ImageType::IndexType    start      = region.GetIndex();
  ImageType::SizeType     size       = region.GetSize();

// With this information we can find the number that will identify the first
// and last slices of the 3D data set. This numerical values are then passed to
// the filenames generator object that will compose the names of the files
// where the slices are going to be stored.
// Typically, the start[2] value is 0, and we like the slice number to be one
// more than the index i.e. it is a count rather than an index value.

  // We add one so that the slice file name is a count from the starting offset
  // rather than an index.  If the initial index is 0, the first file name will
  // start with 001.
  const unsigned int firstSlice = start[2] + 1;
  const unsigned int lastSlice  = start[2] + size[2] - 1 + 1;

  nameGenerator->SetStartIndex( firstSlice );
  nameGenerator->SetEndIndex( lastSlice );
  nameGenerator->SetIncrementIndex( 1 );

  writer->SetFileNames( nameGenerator->GetFileNames() );

//  Finally we trigger the execution of the pipeline with the Update() method
//  on the writer. At this point the slices of the image will be saved in
//  individual files containing a single slice per file. The filenames used for
//  these slices are those produced by the filenames generator.

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while writing the image" << std::endl;
    std::cerr << excp << std::endl;
    }

  return EXIT_SUCCESS;
}
