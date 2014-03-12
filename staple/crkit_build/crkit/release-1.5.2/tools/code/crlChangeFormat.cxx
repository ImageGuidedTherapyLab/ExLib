
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkCastImageFilter.h>

#include <iostream>
#include <tclap/CmdLine.h>

static int const ImageDimension = 3;
typedef itk::Image<float,ImageDimension>  ImageTypeIn;

template <class T> 
void writeOutput( ImageTypeIn::Pointer& input, 
                    std::string outFile, bool setCompressionOn)
{
    typedef itk::Image<T,ImageDimension>                ImageTypeOut;
    typedef itk::ImageFileWriter< ImageTypeOut >        ImageWriterType;

    typedef itk::CastImageFilter< ImageTypeIn, ImageTypeOut> CasterType;

    typename ImageWriterType::Pointer w = ImageWriterType::New();
    typename CasterType::Pointer c = CasterType::New();

    c->SetInput( input );

    w->SetFileName(outFile);
    w->SetInput(c->GetOutput());
    w->SetUseCompression( setCompressionOn );
    try 
      {
      w->Update();
      w->Write();
      }
    catch ( itk::ExceptionObject& err ) 
      {
      std::cerr << "error writing output file " << outFile << ": " <<
        err << std::endl;
      exit (EXIT_FAILURE);
      }
}

int main(int argc, char *argv[])
{

  std::string inpFile;
  std::string outFile;
  std::string protoFile;
  bool setCompressionOn;
  bool s16;

  try
    {
    TCLAP::CmdLine cmd("Change image format based on filename extension", 
        ' ', "1.0" );
    TCLAP::SwitchArg setCompressionOnArg( "C", "setCompressionOn", 
        "tell the ITK writer to compress the file, if possible", 
        cmd, false);
    TCLAP::SwitchArg s16Arg( "s", "s16", 
        "tell ITK to write a 16-bit signed integer file, if possible.",
        cmd, false);
    TCLAP::UnlabeledValueArg<std::string> inpFileArg( "input",
        "the pathname of the file to convert", true, "", "input filename",
        cmd );
    TCLAP::UnlabeledValueArg<std::string> outFileArg( "output",
        "the pathname of output file", true, "", "output filename",
        cmd );
    TCLAP::ValueArg<std::string> protoFileArg( "p", "proto",
        "the pathname of a file to take orientation and metadata from",
        false, "", "input image pathname", cmd);


    cmd.parse(argc, argv);
    inpFile = inpFileArg.getValue();
    outFile = outFileArg.getValue();
    protoFile = protoFileArg.getValue();
    setCompressionOn = setCompressionOnArg.getValue();
    s16 = s16Arg.getValue();
    }
  catch (TCLAP::ArgException& e )
    {
    std::cerr << "error: " << e.error() << " for argument " << e.argId()
      << std::endl;
    exit(EXIT_FAILURE);
    }

  // The convention is for silent operation when it works, error messages
  // when it doesn't.
  // std::cerr << "input file=" << inpFile << std::endl;
  // std::cerr << "output file=" << outFile << std::endl;
  // std::cerr << "proto file=" << ( ( protoFile == "" ) ? "no proto file specified" : protoFile.c_str() ) << std::endl;
  // std::cerr << "compression=" << ( setCompressionOn ? "on" : "off" ) << std::endl;
  // std::cerr << "file type=" << ( s16 ? "S16" : "FLOAT" ) << std::endl;

  typedef itk::ImageFileReader< ImageTypeIn > ImageReaderType;
  ImageReaderType::Pointer r;

  r = ImageReaderType::New();

  ImageTypeIn::Pointer input;
  ImageTypeIn::Pointer proto;

  try 
    {
    r->SetFileName(inpFile);
    r->Update();
    input = r->GetOutput();
    input->DisconnectPipeline();

    if ( protoFile != "" )
      {
      r->SetFileName( protoFile );
      r->Update();
      proto = r->GetOutput();
      proto->DisconnectPipeline();

      input->SetSpacing( proto->GetSpacing() );
      input->SetOrigin( proto->GetOrigin() );
      input->SetDirection( proto->GetDirection() );
      }
    } 
  catch (itk::ExceptionObject &e) 
    {
      std::cerr << "Caught ITK exception reading inputs: " << e << std::endl;
      exit(1);
    }

  if ( s16 )
    {
    writeOutput<short>( input, outFile, setCompressionOn );
    }
  else
    {
    writeOutput<float>( input, outFile, setCompressionOn );
    }


  exit(0);
}
