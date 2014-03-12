
#include <vector>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIteratorWithIndex.h"

int main( int argc, char ** argv )
{
  if ( argc != 6 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0]
              << " inputOutputImage roiImage"
              << " startx starty startz"
              << std::endl;
      std::cerr << "Copies the entire image roiImage over a region of " 
              << "the first image." << std::endl;
      return -1;
    }

  typedef float PixelType;
  typedef itk::Image< PixelType, 3 >  ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  typedef itk::ImageRegionIteratorWithIndex< ImageType>        IteratorType;

  ReaderType::Pointer reader1 = ReaderType::New();
  ReaderType::Pointer reader2 = ReaderType::New();
  reader1->SetFileName( argv[1] );
  reader2->SetFileName( argv[2] );
  try
    {
    reader1->Update();
    reader2->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }

  // Region in the roi image
  ImageType::RegionType roiRegion;
  roiRegion = reader2->GetOutput()->GetRequestedRegion();

  // Region in the input/output image
  ImageType::RegionType overwriteRegion;
  ImageType::RegionType::IndexType inputStart;
  inputStart[0] = atoi(argv[3]);
  inputStart[1] = atoi(argv[4]);
  inputStart[2] = atoi(argv[5]);
  ImageType::RegionType::SizeType  size = roiRegion.GetSize();
  overwriteRegion.SetSize( size );
  overwriteRegion.SetIndex( inputStart );

  IteratorType it1( reader1->GetOutput(),
                               overwriteRegion );

  IteratorType it2( reader2->GetOutput(), 
                    reader2->GetOutput()->GetRequestedRegion());

  for (it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2)
    {
      it1.Set( it2.Get() );
    }


  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[1] );

  writer->SetInput( reader1->GetOutput() );
    try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

