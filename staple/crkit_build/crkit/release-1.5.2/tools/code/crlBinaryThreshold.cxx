
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkBinaryThresholdImageFilter.h>

#include <iostream>

int main(int argc, char *argv[])
{

  static int const ImageDimensionIn = 3;
  static int const ImageDimensionOut = 3;
  typedef   float PixelType;

  typedef   itk::Image<PixelType,ImageDimensionIn>  ImageTypeIn;
  typedef   itk::Image<PixelType,ImageDimensionOut>  ImageTypeOut;
  typedef   itk::ImageFileReader< ImageTypeIn >    ImageReaderType;
  typedef   itk::ImageFileWriter< ImageTypeOut >    ImageWriterType;
  typedef   itk::BinaryThresholdImageFilter< ImageTypeIn, ImageTypeOut> 
                                                    FilterType;

  if (argc != 7) {
    std::cout << 
      "Usage: binarythreshold inimage outimage lowthreshold hithreshold invalue outvalue" << std::endl;
    std::cout << "values equal to and between the thresholds are set to invalue"
  << std::endl;
    std::cout << "values outside the thresholds are set to outvalue" 
  << std::endl;
    exit(1);
  }

  float thresholdlow = atof(argv[3]);
  float thresholdhi = atof(argv[4]);
  float invalue = atof(argv[5]);
  float outvalue = atof(argv[6]);

  ImageReaderType::Pointer r;
  ImageWriterType::Pointer w;
  FilterType::Pointer f;

  r = ImageReaderType::New();
  r->SetFileName(argv[1]);
  r->Update();
  w = ImageWriterType::New();
  w->SetFileName(argv[2]);

  f = FilterType::New();
  f->SetInput( r->GetOutput() );
  f->SetInsideValue (invalue);
  f->SetOutsideValue (outvalue);
  f->SetLowerThreshold(thresholdlow);
  f->SetUpperThreshold(thresholdhi);

  f->Update();
  w->SetInput(f->GetOutput());

  try {
    w->Update();
    w->Write();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  exit(0); // success
}

