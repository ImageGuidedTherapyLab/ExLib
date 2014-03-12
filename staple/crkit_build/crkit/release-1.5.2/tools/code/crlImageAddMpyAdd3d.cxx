
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>

#include <iostream>

int main(int argc, char *argv[])
{

  static int const ImageDimension = 3;
  typedef   float PixelType;

  typedef   itk::Image<PixelType,ImageDimension>  ImageType;
  typedef   itk::ImageFileReader< ImageType >    ImageReaderType;
  typedef   itk::ImageRegionIterator<ImageType> ImageIterator;
  typedef   itk::ImageFileWriter< ImageType >    ImageWriterType;

  if (argc != 6) {
    std::cout << 
      "Usage: crlImageAddMpyAdd inmage add1 mult1 add2 outimage" << std::endl;
    std::cout << 
      "Forms voxelwise outimage = (inimage + add1)*mult1 + add2" << std::endl;
    exit(1);
  }

  PixelType add1 = 0.0;
  PixelType mult1 = 0.0;
  PixelType add2 = 0.0;

  add1 = atof(argv[2]);
  mult1 = atof(argv[3]);
  add2 = atof(argv[4]);

  ImageType::Pointer image1;
  image1 = ImageType::New();

  ImageWriterType::Pointer w;
  w = ImageWriterType::New();
  w->SetFileName( argv[5] );

  ImageReaderType::Pointer image1reader;
  image1reader = ImageReaderType::New();
  image1reader->SetFileName( argv[1] );
  try {
    image1reader->Update();
    image1 = image1reader->GetOutput();
  } catch ( itk::ExceptionObject &err )
  {
    std::cout << "ExceptionObject caught!" << std::endl;
    std::cout << err << std::endl;
    return -1;
  }

  image1->Update();

  // using itk iterators. Fast and beautiful.
  PixelType outval = 0.0;
  ImageIterator iter1(image1,  image1->GetLargestPossibleRegion());
  for (iter1.GoToBegin() ; !iter1.IsAtEnd()  ; ++iter1 ) {
    PixelType p1 = iter1.Get();
    outval = (p1 + add1)*mult1 + add2;
    iter1.Set(outval);
  }

 
  w->SetInput (image1);
  w->UseCompressionOn( );
  
  try
  {
    w->Update();
  }
  catch( itk::ExceptionObject & exp )
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
  }


  exit(0); // success
}
