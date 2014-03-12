
#include <itkOrientedImage.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkImageSeriesReader.h>
#include <itkNumericSeriesFileNames.h>
#include <itkGE5ImageIOFactory.h>
#include <itkGE4ImageIOFactory.h>
#include <itkGE4ImageIO.h>
#include <itkGEAdwImageIOFactory.h>
//#include <itkRawImageIO.h>

#include <tclap/CmdLine.h>

//#error "the ITK GE Signa slices reader does not seem to work, so neither does this code "

int main( int argc, char *argv[] )
{
  std::string inputFirstFile;
  std::string outFile;

  try
    {
    TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',"1.0");
    TCLAP::UnlabeledValueArg<std::string> inputFirstFileArg( "inputFirstFile",
      "the pathname of the *first* slice in the Signa series", true, "", "input filename", cmd );
    TCLAP::UnlabeledValueArg<std::string> outFileArg( "outputFile",
      "a valid file pathname", true, "", "output filename", cmd );

    cmd.parse(argc,argv);

    inputFirstFile = inputFirstFileArg.getValue();
    outFile = outFileArg.getValue();
    }
  catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(EXIT_FAILURE);
    }

  itk::ObjectFactoryBase::RegisterFactory( itk::GEAdwImageIOFactory::New() );
  itk::ObjectFactoryBase::RegisterFactory( itk::GE4ImageIOFactory::New() );
  itk::ObjectFactoryBase::RegisterFactory( itk::GE5ImageIOFactory::New() );

  typedef unsigned short PixelType;
  typedef itk::OrientedImage< PixelType, 3> ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputFirstFile );

  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outFile );


  try
    {
//    reader->SetImageIO( itk::GE4ImageIO::New() );
//    reader->SetDebug(true);
//    reader->Update();
//    reader->GetImageIO()->SetDebug(true);
    writer->SetInput( reader->GetOutput() );
    writer->Update();
  reader->GetOutput()->Print(std::cerr, itk::Indent(0));
    }
  catch (itk::ExceptionObject& err )
    {
    std::cerr << err << std::endl;
    return 1;
    }


}
