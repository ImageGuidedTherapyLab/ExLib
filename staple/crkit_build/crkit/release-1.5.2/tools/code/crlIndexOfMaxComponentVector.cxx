
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkVariableLengthVector.h>

#include <iostream>

int main(int argc, char *argv[])
{

  static int const InImageDimension = 3;
  static int const OutImageDimension = 3;
  typedef   float PixelType;

  typedef   itk::VectorImage<PixelType,InImageDimension>  InImageType;
  typedef   itk::Image<PixelType,OutImageDimension>  OutImageType;
  typedef   itk::ImageFileReader< InImageType >    ImageReaderType;
  typedef   itk::ImageFileWriter< OutImageType >    ImageWriterType;

  typedef OutImageType::RegionType OutRegionType;
  typedef OutImageType::SizeType   OutSizeType;
  typedef OutImageType::IndexType  OutIndexType;
  typedef OutImageType::Pointer    OutImagePointer;
  typedef itk::ImageRegionIterator< OutImageType >
                     OutRegionIteratorType;

  typedef InImageType::RegionType InRegionType;
  typedef InImageType::SizeType   InSizeType;
  typedef InImageType::IndexType  InIndexType;
  typedef InImageType::Pointer    InImagePointer;
  typedef itk::ImageRegionConstIterator< InImageType >
                     InRegionIteratorType;

  if (argc != 3) {
    std::cout << 
      "Usage: " << argv[0] << " inimage outimage " << std::endl;
    std::cout << "Create a new image where each voxel has the " <<
      "value of the index of the component with the largest value." << std::endl;
    exit(1);
  }

  ImageReaderType::Pointer r;
  ImageWriterType::Pointer w;

  r = ImageReaderType::New();
  r->SetFileName(argv[1]);
  r->Update(); // read in the file , could catch errors here...
  InImagePointer inImage = r->GetOutput();
  inImage->DisconnectPipeline();	// why not...

  InSizeType inSize;
  inSize = inImage->GetLargestPossibleRegion().GetSize();

  w = ImageWriterType::New();
  w->SetFileName(argv[2]);

  
  OutImageType::RegionType outRegion = inImage->GetLargestPossibleRegion();
  OutImageType::Pointer outImage = OutImageType::New();
  outImage->CopyInformation( inImage );
  outImage->SetRegions( outRegion );
  outImage->Allocate();

  // Here in the code, we have the input image read from disk and the output
  // image of the desired size allocated.
  OutRegionIteratorType outItr(outImage, outRegion);
  InRegionIteratorType inItr(inImage, inImage->GetLargestPossibleRegion());

  itk::VariableLengthVector<PixelType> vtemp;
  PixelType max;
  unsigned int imax;
  for (outItr.Begin(), inItr.Begin(); !outItr.IsAtEnd() && !inItr.IsAtEnd(); ++outItr, ++inItr) 
    {
    vtemp = inItr.Get();
    max = vtemp[0];
    imax = 0;
    for ( unsigned int i = 1; i < vtemp.GetSize(); i++)
      {
      if ( vtemp[i] > max )
	{
	imax = i;
	max = vtemp[i];
	}
      }
    outItr.Set(imax);
    }

  w->SetInput( outImage );

  try {
    w->Update();
    w->UseCompressionOn();
    w->Write();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  exit(0); // success
}

