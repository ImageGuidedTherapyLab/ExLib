
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <iostream>

int main(int argc, char *argv[])
{

  static int const ImageDimension = 2;
  typedef   float PixelType;

  typedef   itk::Image<PixelType,ImageDimension>  ImageTypeIn;
  typedef   itk::Image<PixelType,ImageDimension>  ImageTypeOut;
  typedef   itk::ImageFileReader< ImageTypeIn >    ImageReaderType;
  typedef   itk::ImageFileWriter< ImageTypeOut >    ImageWriterType;

  if ( (argc != 3) && (argc != 4) ) {
    std::cout << 
      "Usage: crlChangeFormat inimage outimage --setCompressionOn " << std::endl;
    exit(1);
  }

  ImageReaderType::Pointer r;
  ImageWriterType::Pointer w;

  r = ImageReaderType::New();
  r->SetFileName(argv[1]);
  try {
    r->Update();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  w = ImageWriterType::New();
  w->SetFileName(argv[2]);

  w->SetInput(r->GetOutput());

  if (argc == 4) {
    // Assume we have set the flag for compression to be on.
    w->UseCompressionOn( );
  } else {
    w->UseCompressionOff( );
  }

  try {
    w->Update();
    w->Write();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  exit(0); // success
}

