
#include <itkOrientedImage.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include "tclap/CmdLine.h"
#include <itkRGBPixel.h>
#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#include <vtkImageMapToColors.h>
#include <vtkLookupTable.h>
#include <vtkImageMapToWindowLevelColors.h>
#include <itkStatisticsImageFilter.h>

//#include "configuration.h"


int main( int argc, char *argv[] )
{
  std::string inputFile;
  std::string outFile;
  std::string type;
  bool log;
  
  try
    {
      TCLAP::CmdLine cmd("Convert scalar image to RGB image",' ',"");
      TCLAP::UnlabeledValueArg<std::string> inputFileArg( "inputFile",
							  "a valid file pathname", true, "", "input filename", cmd );
      TCLAP::UnlabeledValueArg<std::string> outFileArg( "outputFile",
							"a valid file pathname", true, "", "output filename", cmd );
      TCLAP::ValueArg<std::string> typeArg("t","imageType","image type (fmri, probmap, segmentation, grayscale)",false,"grayscale","image type",cmd);
      TCLAP::SwitchArg logArg("l","log10","use logarithmic scale for probability maps",cmd,false);
      
      cmd.parse(argc,argv);
      
      inputFile = inputFileArg.getValue();
      outFile = outFileArg.getValue();
      if (typeArg.isSet()) type = typeArg.getValue();
      log = logArg.isSet();
    }
  catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(EXIT_FAILURE);
    }
  
  typedef itk::OrientedImage<float,3> ImageType;
  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  typedef itk::OrientedImage<itk::RGBPixel<unsigned char>,3> RGBImageType;
  typedef itk::ImageFileWriter< RGBImageType > ImageWriterType;
  
  itk::ImageToVTKImageFilter<ImageType>::Pointer itktovtk = itk::ImageToVTKImageFilter<ImageType>::New();
  itk::VTKImageToImageFilter<RGBImageType>::Pointer vtktoitk = itk::VTKImageToImageFilter<RGBImageType>::New();
  
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
  
  itktovtk->SetInput(reader->GetOutput());
  itktovtk->Update();
  
  if (type == "fmri")
    {
      // fMRI lookup table
      vtkLookupTable *lut = vtkLookupTable::New();
      //lut->SetTableRange(itktovtk->GetOutput()->GetScalarRange()[0],itktovtk->GetOutput()->GetScalarRange()[1]);
      lut->SetHueRange(0.0,0.18);
      lut->SetTableRange(0.001,itktovtk->GetOutput()->GetScalarRange()[1]);
      lut->SetSaturationRange(1,1);
      lut->SetValueRange(1,1);
      lut->Build();
      lut->SetTableValue(0,0,0,0);
      
      vtkImageMapToColors *mapper = vtkImageMapToColors::New();
      mapper->SetInput(itktovtk->GetOutput());
      mapper->SetOutputFormatToRGB();
      mapper->SetLookupTable(lut);
      
      mapper->Update();
      
      vtktoitk->SetInput(mapper->GetOutput());
    }
  else if (type == "probmap")
    {
      // probabibilty map lookup table
      vtkLookupTable *lut = vtkLookupTable::New();
      lut->SetTableRange(itktovtk->GetOutput()->GetScalarRange()[0],itktovtk->GetOutput()->GetScalarRange()[1]);
      lut->SetHueRange(0.7,0.0);
      lut->SetSaturationRange(1,1);
      lut->SetValueRange(1,1);
      if (log)
	lut->SetScaleToLog10();
      lut->Build();
      lut->SetTableValue(0,0,0,0);
      
      vtkImageMapToColors *mapper = vtkImageMapToColors::New();
      mapper->SetInput(itktovtk->GetOutput());
      mapper->SetOutputFormatToRGB();
      mapper->SetLookupTable(lut);
      
      mapper->Update();
      
      vtktoitk->SetInput(mapper->GetOutput());
    }
  else if (type == "segmentation")
    {
      // segmentation lookup table
      vtkLookupTable *lut = vtkLookupTable::New();
      lut->SetNumberOfColors(256);
      lut->SetTableRange(0,255);
      lut->Build();
      lut->SetTableValue(0,0.0,0.0,0.0,0.0);
      for (unsigned int n=0; n < 255; n+=17)
	{
	  lut->SetTableValue(n+1,0.7,0.0,1.0);
	  lut->SetTableValue(n+2,0.0,1.0,1.0);
	  lut->SetTableValue(n+3,1.0,0.55,0.55);
	  lut->SetTableValue(n+4,0.59,0.59,0.59);
	  lut->SetTableValue(n+5,0.0,0.0,1.0);
	  lut->SetTableValue(n+6,1.0,0.55,0.0);
	  lut->SetTableValue(n+7,1.0,0.0,0.0);
	  lut->SetTableValue(n+8,0.97,0.97,1.0);
	  lut->SetTableValue(n+9,0.64,0.32,0.0);
	  lut->SetTableValue(n+10,1.0,0.0,1.0);
	  lut->SetTableValue(n+11,0.5,0.5,0.5);
	  lut->SetTableValue(n+12,1.0,0.66,0.14);
	  lut->SetTableValue(n+13,0.5,1.0,0.0);
	  lut->SetTableValue(n+14,0.5,1.0,0.83);
	  lut->SetTableValue(n+15,0.39,0.58,0.93);
	  lut->SetTableValue(n+16,0.54,0.17,0.89);
	  lut->SetTableValue(n+17,0.99,0.99,0.99);
	}
      
      vtkImageMapToColors *mapper = vtkImageMapToColors::New();
      mapper->SetInput(itktovtk->GetOutput());
      mapper->SetOutputFormatToRGB();
      mapper->SetLookupTable(lut);      
      mapper->Update();
      
      vtktoitk->SetInput(mapper->GetOutput());
    }
  else
    {
      // grayscale
      vtkImageMapToWindowLevelColors *mapper = vtkImageMapToWindowLevelColors::New();
      mapper->SetInput(itktovtk->GetOutput());
      mapper->SetOutputFormatToRGB();
      itk::StatisticsImageFilter<ImageType>::Pointer stats = itk::StatisticsImageFilter<ImageType>::New();
      stats->SetInput(reader->GetOutput());
      stats->Update();
      mapper->SetLevel(2.0*stats->GetMean());
      mapper->SetWindow(5.0*stats->GetSigma());
      mapper->Update();
      vtktoitk->SetInput(mapper->GetOutput());
    }
  
  vtktoitk->Update();
  
  // write new image
  ImageWriterType::Pointer writer = ImageWriterType::New();
  RGBImageType::Pointer image = const_cast<RGBImageType*>(vtktoitk->GetOutput());
  image->CopyInformation(reader->GetOutput());
  
  writer->SetInput(image);
  
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
