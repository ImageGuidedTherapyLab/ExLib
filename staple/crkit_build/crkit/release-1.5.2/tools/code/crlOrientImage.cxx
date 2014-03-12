
// Other interesting programs would involve
// PermuteImageAxes Filter and
// FlipImageAxes Filter

#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkOrientImageFilter.h>

#include "configuration.h"
#include <tclap/CmdLine.h>

int main( int argc, char ** argv )
{

  std::string *inputImageFile = new std::string("");
  std::string *outputImageFile = new std::string("");
  std::string *outputOrientation = new std::string("axial");
  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',
          CRKIT_VERSION_STRING );
      TCLAP::UnlabeledValueArg<std::string> inputImageArg("inputImageFile","Input Image File",true,"","input image file",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputImageArg("outputImageFile","Output Image File",true,"","output image file",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputOrientationArg("outputOrientation","Output Orientation: axial (default), sagittal, or coronal",false,"axial","output orientation",cmd);

      cmd.parse(argc,argv);

      if (inputImageArg.isSet()) inputImageFile = new std::string(inputImageArg.getValue());
      if (outputImageArg.isSet()) outputImageFile = new std::string(outputImageArg.getValue());
      if (outputOrientationArg.isSet()) outputOrientation = new std::string(outputOrientationArg.getValue());
    }
  catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << "for argument " << e.argId() << std::endl;
      exit(-1);
    }


  typedef float PixelType;
  typedef itk::OrientedImage< PixelType, 3 >  ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFile->c_str() );
  try
    {
    reader->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }

  itk::OrientImageFilter<ImageType,ImageType>::Pointer orienter =
    itk::OrientImageFilter<ImageType,ImageType>::New();

  orienter->UseImageDirectionOn();
  if (argc == 4) {
    if (*outputOrientation == "axial") {
      orienter->SetDesiredCoordinateOrientationToAxial();
    } else if (*outputOrientation == "sagittal") {
      orienter->SetDesiredCoordinateOrientationToSagittal();
    } else if (*outputOrientation == "coronal") {
      orienter->SetDesiredCoordinateOrientationToCoronal();
    } else if (*outputOrientation == "RIP") {
      orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP);
    } else {
      std::cerr << "Can't understand orientation of " << (*outputOrientation) << std::endl;
      std::cerr << "Use one of : axial sagittal coronal RIP" << std::endl;
      exit(1);
    }
  } else {
    // select the usual DICOM LPS convention in this way:
    orienter->SetDesiredCoordinateOrientationToAxial();
  }
  orienter->SetInput(reader->GetOutput());
  orienter->Update();

  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputImageFile->c_str() );

  try
    {
    writer->SetInput(orienter->GetOutput());
    writer->UseCompressionOn( );
    writer->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }

  return 0;

}
