
/* Calculate the rectangular extent of a particular label in an image. */

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImage.h"

int main( int argc, char ** argv )
{
  // Verify the number of parameters in the command line
  if( argc != 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile  label" << std::endl;
      return EXIT_FAILURE;
    }

  typedef float InputPixelType;
  const   unsigned int        Dimension = 3;

  typedef itk::Image< InputPixelType,  Dimension >    InputImageType;
  typedef itk::ImageRegionIteratorWithIndex< InputImageType> IteratorType;

  typedef itk::ImageFileReader< InputImageType  >  ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  InputPixelType labelValue = atof( argv[2] );

  try {
    reader->Update();
  }
  catch ( itk::ExceptionObject &err) {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }

  InputImageType::IndexType lowIndex;
  InputImageType::IndexType highIndex;
  InputImageType::IndexType curIndex;

  for (unsigned int i = 0; i < InputImageType::ImageDimension; i++) {
    lowIndex[i] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[i];
    highIndex[i] = reader->GetOutput()->GetLargestPossibleRegion().GetIndex()[i];
  }


  bool foundTheLabel = false;
  IteratorType it(reader->GetOutput(), 
      reader->GetOutput()->GetLargestPossibleRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    // If the image label matches the specified value, then
    // adjust the region of interest.
    if (it.Value() == labelValue) {
      foundTheLabel = true;
      curIndex = it.GetIndex();
      for (unsigned int i = 0; i < Dimension; i++) {
        if (curIndex[i] <  lowIndex[i])  lowIndex[i] = curIndex[i];
        if (curIndex[i] > highIndex[i]) highIndex[i] = curIndex[i];
      }
    }
  }

  std::cout << "Lower limit: " << lowIndex << ", Upper limit: " << highIndex
    << std::endl;

  // Indicate whether this was worthwhile or not.
  if (foundTheLabel) {
    return EXIT_SUCCESS;
  } else {
    return EXIT_FAILURE;
  }
}

