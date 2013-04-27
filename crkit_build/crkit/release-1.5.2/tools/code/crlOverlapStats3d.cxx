
#include <vector>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIteratorWithIndex.h"

int main( int argc, char ** argv )
{
  if ( argc != 4 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << " inputImageFile1 inputImageFile2 "
              << " foregroundValue "
              << std::endl;
    std::cerr << "Compute binary spatial overlap measures." << std::endl;
    return -1;
    }

  float foregroundValue = atof(argv[3]);

  typedef float PixelType;
  typedef itk::Image< PixelType, 3 >  ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;

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

  /* Now iterator over reader1 and compare the pixel values in the
   * tentative output.
   */
  IteratorType it1( reader1->GetOutput(),
                               reader1->GetOutput()->GetRequestedRegion() );

  IteratorType it2( reader2->GetOutput(), 
                    reader2->GetOutput()->GetRequestedRegion());

  unsigned long intersection = 0;
  unsigned long unioncount = 0;
  unsigned long im1foreground = 0;
  unsigned long im1background = 0;
  unsigned long im2foreground = 0;
  unsigned long im2background = 0;
  for (it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2)
    {
      float vox1 = it1.Value();
      float vox2 = it2.Value();
      if (vox1 == foregroundValue) {
        ++unioncount;
        ++im1foreground;
        if (vox2 == foregroundValue) {
          ++intersection;
          ++im2foreground;
        } else {
          ++im2background;
        }
      } else {
        ++im1background;
        if (vox2 == foregroundValue) {
          ++unioncount;
          ++im2foreground;
        } else {
          ++im2background;
        }
      }
    }

  double jaccard = static_cast<double>(intersection)/
                   static_cast<double>(unioncount);
  double dice = (2.0 * jaccard)/(jaccard + 1.0);

  std::cout 
            << " intersection " << intersection 
            << " union " << unioncount 
            << " jaccard " << jaccard 
            << " dice " << dice 
            << " im1background " << im1background
            << " im1foreground " << im1foreground
            << " im2background " << im2background
            << " im2foreground " << im2foreground
            << std::endl;

  return 0;

}

