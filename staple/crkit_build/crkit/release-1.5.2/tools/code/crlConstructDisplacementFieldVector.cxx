/*
 * This program reads in three files, consisting of the x, y and z
 * components of a displacement field, and writes out a vector format
 * displacement field compatible with crlDeformImage.
 */

#include <iostream>
#include <itkImage.h>
#include <itkCompose3DVectorImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>


int main(int argc , char * argv [])
{
  static int const ImageDimension =  3;

  typedef float PixelType;
  typedef itk::Image< PixelType, ImageDimension > InputImageType;

  typedef itk::Compose3DVectorImageFilter< InputImageType >  FilterType;
  typedef FilterType::OutputImageType  OutputImageType;

  typedef   itk::ImageFileReader< InputImageType >    ImageReaderType;
  typedef   itk::ImageFileWriter< OutputImageType >    ImageWriterType;

  if (argc != 5) {
    std::cout << "Usage: " << argv[0] << 
      " displacementXFileName displacementYFileName displacementZFileName" <<
      " outputVectorFileName" << std::endl;
    return EXIT_FAILURE;
  }


  FilterType::Pointer filter = FilterType::New();

  ImageReaderType::Pointer zeroImage  = ImageReaderType::New();
  ImageReaderType::Pointer oneImage   = ImageReaderType::New();
  ImageReaderType::Pointer twoImage   = ImageReaderType::New();

  zeroImage->SetFileName( argv[1] );
  zeroImage->Update();
  oneImage->SetFileName(  argv[2] );
  oneImage->Update();
  twoImage->SetFileName(  argv[3] );
  twoImage->Update();
 
  filter->SetInput1( zeroImage->GetOutput() );
  filter->SetInput2( oneImage->GetOutput() );
  filter->SetInput3( twoImage->GetOutput() );

  try
    {
    filter->Update();
    }
  catch( itk::ExceptionObject & excp )
   {
   std::cerr << "Exception caught !" << std::endl;
   std::cerr << excp << std::endl;
   return EXIT_FAILURE;
   }

  ImageWriterType::Pointer writer = ImageWriterType::New();

  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[4] ) ;
  
  try {
    writer->Update();
    writer->Write();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }
  
  return EXIT_SUCCESS;

}

