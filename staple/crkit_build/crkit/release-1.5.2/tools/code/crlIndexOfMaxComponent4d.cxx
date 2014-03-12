
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>

#include <iostream>

int main(int argc, char *argv[])
{

  static int const InImageDimension = 4;
  static int const OutImageDimension = 3;
  typedef   float PixelType;

  typedef   itk::Image<PixelType,InImageDimension>  InImageType;
  typedef   itk::Image<PixelType,OutImageDimension>  OutImageType;
  typedef   itk::ImageFileReader< InImageType >    ImageReaderType;
  typedef   itk::ImageFileWriter< OutImageType >    ImageWriterType;

  typedef OutImageType::RegionType OutRegionType;
  typedef OutImageType::SizeType   OutSizeType;
  typedef OutImageType::IndexType  OutIndexType;
  typedef OutImageType::Pointer    OutImagePointer;
  typedef itk::ImageRegionIteratorWithIndex< OutImageType >
                     OutRegionIteratorType;

  typedef InImageType::RegionType InRegionType;
  typedef InImageType::SizeType   InSizeType;
  typedef InImageType::IndexType  InIndexType;
  typedef InImageType::Pointer    InImagePointer;
  typedef itk::ImageRegionIteratorWithIndex< InImageType >
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

  InSizeType inSize;
  inSize = r->GetOutput()->GetLargestPossibleRegion().GetSize();

  w = ImageWriterType::New();
  w->SetFileName(argv[2]);

  OutImageType::IndexType outIndex;
  OutImageType::SizeType outSize;
  OutImageType::SpacingType outSpacing;
  InImageType::SpacingType inSpacing;
  inSpacing = r->GetOutput()->GetSpacing();

  InImageType::PointType  inputOrigin;
  OutImageType::PointType  outputOrigin;
  inputOrigin = r->GetOutput()->GetOrigin();

  InImageType::DirectionType  inputDirection;
  OutImageType::DirectionType  outputDirection;
  inputDirection = r->GetOutput()->GetDirection();

  for (unsigned int i = 0; i < OutImageType::ImageDimension; i++) {
    outSize[i] = r->GetOutput()->GetLargestPossibleRegion().GetSize()[i];
    outIndex[i] = r->GetOutput()->GetLargestPossibleRegion().GetIndex()[i];
    outputOrigin[i] = inputOrigin[i];
    outSpacing[i] = inSpacing[i];
  }

  for (unsigned int r = 0; r < OutImageType::ImageDimension; r++) {
    for (unsigned int c = 0; c < OutImageType::ImageDimension; c++) {
      outputDirection(r,c) = inputDirection(r,c);
    }
  }

  OutImageType::RegionType outRegion;
  outRegion.SetIndex( outIndex );
  outRegion.SetSize( outSize );

  OutImageType::Pointer outImage = OutImageType::New();
  // Adjust for the change in size
  outImage->SetRegions( outRegion );
  outImage->Allocate();

  // Copy geometry information from the input without change
  outImage->SetSpacing( outSpacing );
  outImage->SetDirection( outputDirection );

  // This program is designed for taking an ND array of probability values
  // and constructing an (N-1)D array with the index of the maximum 
  // probability at each voxel. As such, the Nth axis does not contain
  // useful geometry such as spacing, or direction cosines. Here we use the
  // (N-1)D geometry information to set the values for the output.
  outImage->SetOrigin(outputOrigin);
  
  // Here in the code, we have the input image read from disk and the output
  // image of the desired size allocated.
  OutRegionIteratorType outItr(outImage, outRegion);

  InImageType::IndexType inIndex;

  for (outItr.Begin(); !outItr.IsAtEnd(); ++outItr) {
    outIndex = outItr.GetIndex();
    for (unsigned int i = 0; i < outIndex.GetIndexDimension(); i++) {
      inIndex[i] = outIndex[i];
    }
    PixelType max = itk::NumericTraits<PixelType>::min();
    unsigned int maxindex = 0;
    for (unsigned int i = 0; i < inSize[InImageDimension-1]; i++) {
      inIndex[InImageDimension-1] = i;
      if (r->GetOutput()->GetPixel(inIndex) > max) {
        max = r->GetOutput()->GetPixel(inIndex);
        maxindex = i;
      }
    }
    outImage->SetPixel(outIndex, maxindex);
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

