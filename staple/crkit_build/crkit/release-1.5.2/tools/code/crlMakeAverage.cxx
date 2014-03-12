/*
 * This program reads in files, and writes out an image with the arithmetic means of all voxels
 */

#include <iostream>
#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>


int main(int argc , char * argv [])
{
  static int const ImageDimension =  3;

  typedef float PixelType;
  typedef itk::OrientedImage< PixelType, ImageDimension > ImageType;

  typedef   itk::ImageFileReader< ImageType >    ImageReaderType;
  typedef   itk::ImageFileWriter< ImageType >    ImageWriterType;

  if (argc < 3) {
    std::cout << "Usage: " << argv[0] <<
      " image1FileName [image2FileName] [...] [imageNFileName] outputImageFileName" << std::endl;
    return EXIT_FAILURE;
  }


  std::vector<ImageReaderType::Pointer> inputReaders;

  for (int n = 0; n < argc-2; n++)
    {
      ImageReaderType::Pointer tempReader = ImageReaderType::New();
      tempReader->SetFileName(argv[n+1]);
      tempReader->Update();
      inputReaders.push_back(tempReader);
    }

  ImageType::Pointer output = ImageType::New();
  output->CopyInformation(inputReaders[0]->GetOutput());
  output->SetRegions(inputReaders[0]->GetOutput()->GetLargestPossibleRegion());
  output->Allocate();

  std::vector<itk::ImageRegionIterator<ImageType> > iterators;

  for (int n = 0; n < argc-2; n++)
    {
      itk::ImageRegionIterator<ImageType> it(inputReaders[n]->GetOutput(),inputReaders[n]->GetOutput()->GetLargestPossibleRegion());
      iterators.push_back(it);
    }

  itk::ImageRegionIterator<ImageType> oit(output,output->GetLargestPossibleRegion());

  for (int n = 1; n < argc-2; n ++)
    {
      iterators[n].GoToBegin();
    }
  oit.GoToBegin();

  for (iterators[0].GoToBegin(); !iterators[0].IsAtEnd(); ++iterators[0])
    {
    double sum = 0.0;
    for (int n = 0; n < argc-2; n++)
      {
      PixelType& val=iterators[n].Value();
      sum += val;
      }
    oit.Set(sum/(argc-2));
    for (int n = 1; n < argc-2; n++)
      {
      ++iterators[n];
      }

    ++oit;
    }

  ImageWriterType::Pointer writer = ImageWriterType::New();

  writer->SetInput( output );
  writer->SetFileName( argv[argc-1] );
  writer->UseCompressionOn();

  try {
    writer->Update();
    writer->Write();
  } catch (itk::ExceptionObject &e) {
    std::cerr << "Caught ITK exception: " << e << std::endl;
    exit(1);
  }

  return EXIT_SUCCESS;

}
