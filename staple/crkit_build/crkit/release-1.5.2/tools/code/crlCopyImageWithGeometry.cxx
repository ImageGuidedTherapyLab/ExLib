

// this program reads in an image and writes it out with metadata borrowed from a 3rd party image.  
// The idea is to use an example image to get geometry data from.
//
// We'll fix three things:
// 1. orientation
// 2. center
// 3. voxel spacing/size
//


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "configuration.h"
#include <tclap/CmdLine.h>

int main( int argc, char ** argv )
{

  TCLAP::CmdLine cmd("Computational Radiology Laboratory", ' ',CRKIT_VERSION_STRING);
  TCLAP::SwitchArg directionArg("D", "direction","copy direction (orientation) information", cmd, false);
  TCLAP::SwitchArg originArg("O", "origin","copy origin information", cmd, false);
  TCLAP::SwitchArg spacingArg("S", "spacing","copy spacing information", cmd, false);
  TCLAP::SwitchArg noCompressArg("z", "nocompress", "DON'T compress output (default is to compress)",cmd, false);
  TCLAP::SwitchArg nullifyArg("N", "nullify", "clean out orientation and origin information", cmd, false);
  TCLAP::ValueArg<std::string> outputArg("o", "out", "output image file", true, "","output image", cmd);
  TCLAP::UnlabeledValueArg<std::string> inputArg("inputImage","input image file",true,"","input image",cmd);
  TCLAP::UnlabeledValueArg<std::string> geomArg("geometryImage","geometry example image file (needed for -D -O and/or -S)",false,"","geometry image",cmd);
      
  try
    {
    cmd.parse(argc,argv);
    }
  catch (TCLAP::ArgException& e)
    {
    std::cerr << "error: " << e.error() << "for argument " << e.argId() << std::endl;
    return(1);
    }
  
  bool direction = directionArg.isSet();
  bool origin = originArg.isSet();
  bool spacing = spacingArg.isSet();
  bool nocompress = noCompressArg.isSet();
  bool nullify = nullifyArg.isSet();
  const std::string& inputName = inputArg.getValue();
  const std::string& geomName = geomArg.getValue();
  const std::string& outputName = outputArg.getValue();

  if ( (direction || origin || spacing) && !geomArg.isSet() )
    {
    std::cerr << "you must specify a geometry file if you specify -D -O  or -S" << std::endl;
    return 1;
    }
  if ( nullify )
    {
    if ( direction || origin || spacing ) 
      {
      std::cerr << "you cannot combine -N with -D, -O, or -S" << std::endl;
      return 1;
      }
    if ( geomArg.isSet() )
      {
      std::cerr << "you cannot supply a geometry image when using -N." << std::endl;
      return 1;
      }
    }
    
#if 0
  if (!orientation && !center && !spacing && !compress && !nocompress)
    {
    std::cout << "Nothing to do!  You must specify at least one of -S -C -O -Z or -z." << std::endl;
    return(1);
    }
  
  // we use two to differentiate "not specified" versus specified as compress, but the 
  // default will be to compress.
  if ( !compress && !nocompress )
    {
    compress = true;
    }
#endif 
  
  std::cout << "copying " << inputName << " to " << outputName << " with: ";
  if (direction) std::cout << "direction ";
  if (origin) std::cout << "origin ";
  if (spacing) std::cout << "spacing ";
  if (nullify) std::cout << "clear direction and origin ";
  if (nocompress) std::cout << "no"; std::cout << "compression" << std::endl;
  if ( geomName != "" )
    {
    std::cout << "    from the example file " << geomName << std::endl;
    }
 
  // TODO: add support for other pixel types -- which ones?
  typedef float PixelType;
  typedef itk::Image< PixelType,3 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  ReaderType::Pointer readeri = ReaderType::New();
  readeri->SetFileName(inputName);

  
  try
    {
    readeri->Update();
    }
  catch (itk::ExceptionObject& err)
    {
    std::cerr << "error reading input image: " << inputName << ": " << err << std::endl;
    return 1;
    }
  
  ImageType::Pointer output = readeri->GetOutput();
  output->DisconnectPipeline();	// probably not necessary in this context, but why not.
  
  
  
  if ( geomName != "" )
    {
    ReaderType::Pointer readerg = ReaderType::New();
    readerg->SetFileName(geomName);
    try
      {
      readerg->Update();
      }
    catch (itk::ExceptionObject& err)
      {
      std::cerr << "error reading geometry image: " << geomName << ": " << err << std::endl;
      return 1;
      }

    ImageType::Pointer geom = readerg->GetOutput();
  
    if (spacing)
      {
      output->SetSpacing(geom->GetSpacing());
      }
    if (direction)
      {
      output->SetDirection(geom->GetDirection());
      }
    if (origin)
      {
      output->SetOrigin(geom->GetOrigin());
      }
    }
  
  if (nullify)
      {
      ImageType::DirectionType dir;
      dir.SetIdentity();
      output->SetDirection(dir);
      ImageType::PointType point;
      point.Fill(0.0);
      output->SetOrigin(point);
      }

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputName );

  try
    {
    writer->SetInput(output);
    writer->SetUseCompression( !nocompress );
    writer->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "Error writing output file: " << outputName << std::endl;
    std::cout << err << std::endl;
    return 1;
    }

  
  
  return 0;
}
