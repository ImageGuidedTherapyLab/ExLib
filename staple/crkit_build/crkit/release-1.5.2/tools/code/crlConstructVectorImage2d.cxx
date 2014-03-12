/*
 * This program reads in files, consisting of components of a vector,
 * and writes out a vector image
 */

#include <iostream>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkVariableLengthVector.h>


int main(int argc , char * argv [])
{
  static int const ImageDimension =  2;

  typedef float PixelType;
  typedef itk::Image< PixelType, ImageDimension > InputImageType;

  typedef itk::VectorImage< PixelType, ImageDimension > OutputImageType;

  typedef   itk::ImageFileReader< InputImageType >    ImageReaderType;
  typedef   itk::ImageFileWriter< OutputImageType >    ImageWriterType;

  if (argc < 3) {
    std::cout << "Usage: " << argv[0] <<
      " component1FileName [component2FileName] [componentNFileName] [...] outputVectorFileName" << std::endl;
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

  OutputImageType::Pointer output = OutputImageType::New();
  output->CopyInformation(inputReaders[0]->GetOutput());
  output->SetRegions(inputReaders[0]->GetOutput()->GetLargestPossibleRegion());
  output->SetNumberOfComponentsPerPixel(argc-2);
  output->Allocate();

  std::vector<itk::ImageRegionIterator<InputImageType> > iterators;

  for (int n = 0; n < argc-2; n++)
    {
      itk::ImageRegionIterator<InputImageType> it(inputReaders[n]->GetOutput(),inputReaders[n]->GetOutput()->GetLargestPossibleRegion());
      iterators.push_back(it);
    }

  itk::ImageRegionIterator<OutputImageType> oit(output,output->GetLargestPossibleRegion());

  for (int n = 1; n < argc-2; n ++)
    {
      iterators[n].GoToBegin();
    }

  oit.GoToBegin();

  for (iterators[0].GoToBegin(); !iterators[0].IsAtEnd(); ++iterators[0])
    {
      itk::VariableLengthVector<PixelType> v;
      v.SetSize(argc-2);
      for (int n = 0; n < argc-2; n++)
	{
	  v.SetElement(n,iterators[n].Get());
	}
      oit.Set(v);
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
