

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkVectorImageToImageAdaptor.h"
#include "itkCastImageFilter.h"
#include "itkVectorImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include"itkNthElementImageAdaptor.h"


#include "itkImageRegionIterator.h"

int main( int argc, char * argv[] )
{
  if( argc < 6 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr <<  argv[0] << " inputImageFile   outputImageFile   start#   end#   addOnesToBegin#" << std::endl;
    std::cerr <<  "Parameters:  start# \t the first component to copy (0...N-1)"<< std::endl;
    std::cerr <<  "Parameters:  end# \t the last component to copy (0...N-1)" << std::endl;
    std::cerr <<  "Parameters: addOnesToBegin# \t number of components (all ones) to add." << std::endl;
    
    return 1;
    }

  char * infile = argv[1];
  char * outfile = argv[2];
  unsigned int start = atoi(argv[3]);
  unsigned int end = atoi(argv[4]);
  unsigned int addOnes = atoi(argv[5]);
std::cout << "resize to " << addOnes+end-start+1 << " components" << std::endl;

  /** typedefs for the kdt filter and image reader/writer */

  typedef   float  ComponentType;

  const unsigned int Dimension = 3;

  typedef itk::VectorImage< ComponentType,  Dimension >   InputImageType;
  typedef itk::VectorImage< ComponentType, Dimension >   OutputImageType;


  typedef itk::ImageFileReader< InputImageType,itk::DefaultConvertPixelTraits< ComponentType > >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( infile );
  reader->Update();
  
  if (reader->GetOutput()->GetVectorLength()>end) {
  } else {
    std::cerr << "index is too high, there are only " << reader->GetOutput()->GetVectorLength() << " in this image." << std::endl;
  }
  

  OutputImageType::Pointer output = OutputImageType::New();
  output->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  output->CopyInformation( reader->GetOutput() );
  output->SetVectorLength(addOnes+end-start+1);
  output->Allocate();
  
  typedef itk::ImageRegionIterator<InputImageType> InputImageIteratorType;
  typedef itk::ImageRegionIterator<OutputImageType> OutputImageIteratorType;
  InputImageIteratorType itIN( reader->GetOutput(), reader->GetOutput()->GetBufferedRegion() );
  OutputImageIteratorType itOUT( output, output->GetBufferedRegion() );
  
  for (itIN.GoToBegin(), itOUT.GoToBegin(); !itIN.IsAtEnd(); ++itIN, ++itOUT)  {
    InputImageType::PixelType pin = itIN.Get();
    OutputImageType::PixelType pout(addOnes+end-start+1);
    pout.Fill(1);
    for (unsigned int ind=start; ind<end+1; ind++) {
      pout[addOnes+ind-start] = pin[ind];
    }
    itOUT.Set(pout);
  }
  
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outfile );
  writer->SetInput( output );
  writer->Update();


  return 0;
}

