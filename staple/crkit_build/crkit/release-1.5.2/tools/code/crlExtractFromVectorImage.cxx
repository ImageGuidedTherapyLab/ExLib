
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkImageIterator.h>
#include <itkVectorImage.h>
#include <itkImageAdaptor.h>
#include <itkVectorImageToImageAdaptor.h>
#include <itkOrientedImage.h>

#include <iostream>

int main(int argc, char *argv[])
{

  static int const ImageDimensionIn = 3;
  static int const ImageDimensionOut = 3;
  typedef   float PixelType; // Modified from short 11.22.2007

  typedef   itk::VectorImage< PixelType , ImageDimensionIn>  ImageTypeIn;
  typedef   itk::OrientedImage<PixelType,ImageDimensionOut>  ImageTypeOut;
  typedef   itk::ImageFileReader< ImageTypeIn,
               itk::DefaultConvertPixelTraits<PixelType> > ImageReaderType;

  typedef   itk::ImageFileWriter< ImageTypeOut >    ImageWriterType;

  typedef itk::ImageRegionIterator< ImageTypeIn > InIteratorType;
  typedef itk::ImageRegionIterator< ImageTypeOut > OutIteratorType;


  if (argc != 4) {
    std::cout <<
      "Usage: " << argv[0] << " inimage outimage index" << std::endl;
    std::cout << "Extract component 'index' from vector image." << std::endl;
    exit(1);
  }

  int indexofcomponenttoextract = atoi(argv[3]);

  ImageReaderType::Pointer r;
  ImageWriterType::Pointer w;

  r = ImageReaderType::New();
  r->SetFileName(argv[1]);
  r->Update();

  ImageTypeOut::Pointer outImage = ImageTypeOut::New();
  outImage->SetRegions(r->GetOutput()->GetLargestPossibleRegion().GetSize());
  outImage->SetSpacing(r->GetOutput()->GetSpacing());
  outImage->SetOrigin(r->GetOutput()->GetOrigin());
  outImage->SetDirection(r->GetOutput()->GetDirection());

  outImage->Allocate();

  InIteratorType inItr(r->GetOutput(),
                       r->GetOutput()->GetLargestPossibleRegion());
  OutIteratorType outItr(outImage, outImage->GetLargestPossibleRegion());

  for (inItr.GoToBegin(), outItr.GoToBegin(); !outItr.IsAtEnd();
                  ++inItr, ++outItr) {
    outItr.Set( inItr.Get()[indexofcomponenttoextract] );
  }

  w = ImageWriterType::New();
  w->SetFileName(argv[2]);
  w->SetInput( outImage );

  try {
    w->Update();
    w->Write();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  exit(0); // success
}
