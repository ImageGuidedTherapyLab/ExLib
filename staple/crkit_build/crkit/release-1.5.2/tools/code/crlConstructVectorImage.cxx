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
  static int const ImageDimension =  3;

  typedef float PixelType;
  typedef itk::Image< PixelType, ImageDimension > InputImageType;

  typedef itk::VectorImage< PixelType, ImageDimension > OutputImageType;

  typedef   itk::ImageFileReader< InputImageType >    ImageReaderType;
  typedef   itk::ImageFileWriter< OutputImageType >    ImageWriterType;

  bool normalize = false;

  if ( argc > 1 && argv[1][0] == '-' && argv[1][1] == 'n' )
    {
      std::cout << "will normalize each vector to [0,1] and sum to 1.0" << std::endl;
      normalize = true;
      argc--;
      argv++;
    }

  if (argc < 3) {
    std::cout << "Usage: " << argv[0] <<
      " [-n] component1FileName [component2FileName] [componentNFileName] [...] outputVectorFileName" << std::endl;
    std::cout << "if -n is the FIRST argument, will normalize each vector to [0,1]" << std::endl;
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

  itk::VariableLengthVector<PixelType> v;	// moved out here for speed.
  v.SetSize(argc-2);

  for (iterators[0].GoToBegin(); !iterators[0].IsAtEnd(); ++iterators[0])
    {
      double sum = 0.0;
      for (int n = 0; n < argc-2; n++)
	{
	  PixelType& val=iterators[n].Value();
	  if ( normalize && val < 0.0 )
	    {
	      std::cerr << "WARNING: negative number clamped to zero." << std::endl;
	      val = 0.0;
	    }
	  v[n] = val;
	  sum += val;
	}
      if ( normalize )
	{
	  if (sum != 0.0)
	    v /= sum;
	  else
	    v.Fill(0.0);
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
