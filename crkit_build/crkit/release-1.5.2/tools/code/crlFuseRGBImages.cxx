
#include <itkOrientedImage.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include "tclap/CmdLine.h"
#include <itkRGBPixel.h>
#include <itkImageRegionIterator.h>

//#include "configuration.h"


int main( int argc, char *argv[] )
{
  std::string inputFile1, inputFile2;
  std::string outFile;
  float opacity = 0.5;

  try
    {
      TCLAP::CmdLine cmd("Fuse two RGB images",' ',"");
      TCLAP::UnlabeledValueArg<std::string> inputFileArg1( "inputFile1",
							  "a valid file pathname", true, "", "input filename 1", cmd );
      TCLAP::UnlabeledValueArg<std::string> inputFileArg2( "inputFile2",
							  "a valid file pathname", true, "", "input filename 2", cmd );
      TCLAP::UnlabeledValueArg<std::string> outFileArg( "outputFile",
							"a valid file pathname", true, "", "output filename", cmd );
      TCLAP::UnlabeledValueArg<float> opacityArg("opacity","a number between 0 and 1",false,0.5,"opacity",cmd);

      cmd.parse(argc,argv);

      inputFile1 = inputFileArg1.getValue();
      inputFile2 = inputFileArg2.getValue();
      outFile = outFileArg.getValue();
      if (opacityArg.isSet())
	opacity = opacityArg.getValue();
    }
  catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(EXIT_FAILURE);
    }

  typedef itk::OrientedImage<itk::RGBPixel<unsigned char>,3> ImageType;
  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  typedef itk::ImageFileWriter< ImageType > ImageWriterType;

  ImageReaderType::Pointer reader1 = ImageReaderType::New();
  ImageReaderType::Pointer reader2 = ImageReaderType::New();
  reader1->SetFileName(inputFile1);
  reader2->SetFileName(inputFile2);
  try
    {
      reader1->Update();
      reader2->Update();
    }
  catch ( itk::ExceptionObject& err )
    {
      std::cerr << "error reading input: " << err <<
	std::endl;
      return 1;
    }

  ImageType::Pointer output = ImageType::New();
  output->SetRegions(reader1->GetOutput()->GetLargestPossibleRegion());
  output->Allocate();
  output->CopyInformation(reader1->GetOutput());

  itk::ImageRegionIterator<ImageType> it1(reader1->GetOutput(),reader1->GetOutput()->GetLargestPossibleRegion());
  itk::ImageRegionIterator<ImageType> it2(reader2->GetOutput(),reader2->GetOutput()->GetLargestPossibleRegion());
  itk::ImageRegionIterator<ImageType> oit(output,output->GetLargestPossibleRegion());

  it1.GoToBegin();
  it2.GoToBegin();
  unsigned int i;
  itk::RGBPixel<unsigned char> in1,in2,out;

  for (oit.GoToBegin(); !oit.IsAtEnd(); ++oit)
    {
      in1 = it1.Value();
      in2 = it2.Value();
      if (in2[0] || in2[1] || in2[2])
	{
	  for (i=0;i<3;i++)
	    out[i] = (unsigned char)floor(opacity*in1[i]+(1-opacity)*in2[i]);
	  oit.Set(out);
	}
      else
	{
	  oit.Set(in1);
	}

      ++it1;
      ++it2;
    }

  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput(output);

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
