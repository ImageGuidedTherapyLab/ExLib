
#ifndef _CRLVECTOR3TORGB_TXX_INCLUDED
#define _CRLVECTOR3TORGB_TXX_INCLUDED 1

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkRGBPixel.h>
#include <math.h>

template <const unsigned int TImageDimension, class TComponentType>
void crl::Vector3ToRGB<TImageDimension, TComponentType>
::EstimateRGB(PixelType &in, ColorPixelType &color)
{

    // Let's convert the vector to be unit magnitude.
    double vectornorm = in[0]*in[0] + in[1]*in[1]+ in[2]*in[2];

    color[0] = static_cast<unsigned char>( 
                  floor(fabs(in[0]*255.0/vectornorm))
               );
    color[1] = static_cast<unsigned char>(
                  floor(fabs(in[1]*255.0/vectornorm))
               );
    color[2] = static_cast<unsigned char>(
                  floor(fabs(in[2]*255.0/vectornorm))
               );
}

template <const unsigned int TImageDimension, class TComponentType>
int crl::Vector3ToRGB<TImageDimension, TComponentType>
::Execute()
{
  typedef   itk::ImageRegionConstIteratorWithIndex<ImageType> ImageIterator;

  typename ImageReaderType::Pointer  reader = ImageReaderType::New();

  reader->SetFileName( m_InputFileName );

  try {
    reader->GenerateOutputInformation();
  } catch (itk::ExceptionObject &excp)
  {
    std::cerr << "Failed to read from file " << m_InputFileName << std::endl;    
    std::cerr << "Exception error message is : " << excp << std::endl;
    return EXIT_FAILURE;
  }

  reader->Update();

  // Allocate the new output image:
  typename ColorImageType::Pointer outputImage = ColorImageType::New();
  outputImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  outputImage->CopyInformation( reader->GetOutput() );
  outputImage->Allocate();

  /* Iterate over the Vector3 data , computing the color coding. */

  ImageIterator it(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());
  PixelType in;
  ColorPixelType color;

  typename ImageType::IndexType index;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    in = it.Get();
    index = it.GetIndex();
    EstimateRGB(in, color);
    outputImage->SetPixel(index, color);
  }

  typename ImageWriterType::Pointer writer = ImageWriterType::New();

  writer->SetInput ( outputImage );
  writer->SetFileName( m_OutputFileName );
  try {
    writer->Write(); // 
  } catch ( itk::ExceptionObject & excp ) {
    std::cerr << "Failed to write out the image." << std::endl;
    std::cerr << "Exception error message is : " << excp << std::endl;
  }

  exit(0); // success
}

#endif

