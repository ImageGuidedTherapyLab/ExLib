#include <iostream>
#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkPointSet.h>
#include <itkEuclideanDistancePointMetric.h>
#include <itkIdentityTransform.h>
#include <itkImageRegionIterator.h>
#include "crlHausdorffDistanceImageFilter.h"
#include <itkHausdorffDistanceImageFilter.h>


int main(int argc, char *argv[])
{
  if ( argc < 3 ) {
    std::cout << "Usage: " << argv[0] << " imageFile1 imageFile2 [percentage [percentage] [...]]" << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::Image<float,3> ImageType;
  typedef itk::PointSet<float,3> PointType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader1 = ReaderType::New();
  ReaderType::Pointer reader2 = ReaderType::New();

  reader1->SetFileName(argv[1]);
  reader2->SetFileName(argv[2]);

  try
    {
      reader1->Update();
    }
  catch (...)
    {
      std::cerr << "Failed to read input image: " << argv[1] << std::endl;
      return EXIT_FAILURE;
    }

  try
    {
      reader2->Update();
    }
  catch(...)
    {
      std::cerr << "Failed to read input image: " << argv[2] << std::endl;
      return EXIT_FAILURE;
    }

  crlHausdorffDistanceImageFilter<ImageType,ImageType>::Pointer filter = crlHausdorffDistanceImageFilter<ImageType,ImageType>::New();

  filter->SetInput(0,reader1->GetOutput());
  filter->SetInput(1,reader2->GetOutput());
  filter->Update();

  std::cout << "Hausdorff distance: " << filter->GetHausdorffDistance() << std::endl;
  std::cout << "Average Hausdorff distance: " << filter->GetAverageHausdorffDistance() << std::endl;

  for (unsigned int i = 3; i < (unsigned int)argc; i++)
    {
      std::cout << argv[i] << "% Hausdorff distance: ";
      std::cout << filter->GetHausdorffDistance(atoi(argv[i])) << std::endl;
    }

  return EXIT_SUCCESS;
}
