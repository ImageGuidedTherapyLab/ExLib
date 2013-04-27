
#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIteratorWithIndex.h>

#include <itkMeanSquaresImageToImageMetric.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>

#include <itkTranslationTransform.h>
#include <itkEuler2DTransform.h>

#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include <iostream>

int main(int argc, char *argv[])
{

  static unsigned int const ImageDimension = 2;
  typedef   double PixelType;

  typedef   itk::OrientedImage<PixelType,ImageDimension>  ImageType;
  typedef   itk::ImageFileReader< ImageType >    ImageReaderType;

  typedef itk::TranslationTransform< PixelType, ImageDimension> TransformType;
  typedef itk::NearestNeighborInterpolateImageFunction< ImageType, 
                double > InterpolatorType;
  typedef itk::NormalizedCorrelationImageToImageMetric<ImageType, ImageType> 
    MetricType;

  typedef itk::ImageRegionConstIteratorWithIndex< ImageType > IteratorType;

  typedef MetricType::TransformType TransformBaseType;
  typedef TransformBaseType::ParametersType ParametersType;

  typedef ImageType::RegionType RegionType;
  typedef ImageType::SizeType   SizeType;
  typedef ImageType::IndexType  IndexType;

  if (argc != 4) {
    std::cout << 
      "Usage: crlCorrelator targetImageName movingImageName scoreImageName" 
               << std::endl;
    std::cout << "Compute SSD between the two images under all possible"
              << " translations." << std::endl;
    exit(1);
  }

  ImageReaderType::Pointer treader;
  ImageReaderType::Pointer mreader;

  treader = ImageReaderType::New();
  mreader = ImageReaderType::New();
  treader->SetFileName(argv[1]);
  mreader->SetFileName(argv[2]);

  try {
    treader->Update(); // read in the file, catch errors ...
  } catch ( itk::ExceptionObject &err )
  {
    std::cout << "ExceptionObject caught!" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }

  try {
    mreader->Update(); // read in the file, catch errors ...
  } catch ( itk::ExceptionObject &err )
  {
    std::cout << "ExceptionObject caught!" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }

  ImageType::Pointer timage = treader->GetOutput();
  ImageType::Pointer mimage = mreader->GetOutput();

  InterpolatorType::Pointer interp = InterpolatorType::New();
  TransformType::Pointer transform = TransformType::New();
  MetricType::Pointer metric = MetricType::New();

  // Wow, the type hierarchy of these classes is totally broken.
  metric->SetTransform( transform.GetPointer() );
  metric->SetInterpolator( interp.GetPointer() );
  metric->SetMovingImage( mimage.GetPointer() );
  metric->SetFixedImage( timage.GetPointer() );
  metric->SetFixedImageRegion( timage->GetLargestPossibleRegion() );
  interp->SetInputImage( mimage.GetPointer() );


  //------------------------------------------------------------
  // This call is mandatory before start querying the Metric
  // This method do all the necesary connections between the
  // internal components: Interpolator, Transform and Images
  //------------------------------------------------------------
  try {
    metric->Initialize();
  } catch( itk::ExceptionObject & e ) {
    std::cout << "Metric initialization failed" << std::endl;
    std::cout << "Reason " << e.GetDescription() << std::endl;
    return EXIT_FAILURE;
  }

  // Magical incantation only for very special metrics.
  metric->SubtractMeanOn();

  ImageType::PointType point;
  IteratorType iter1( timage,  
                       timage->GetLargestPossibleRegion() );
  ParametersType parameters( transform->GetNumberOfParameters() );
  for (unsigned int k = 0; k < ImageDimension; k++) {
    parameters[k] = 0.0f;
  }

  double metricValue = 0.0;
  for (iter1.GoToBegin(); !iter1.IsAtEnd(); ++iter1) {
    ImageType::IndexType index = iter1.GetIndex();
    timage->TransformIndexToPhysicalPoint( index, point );
    std::cout << "Working on " << point << std::endl;
    for (unsigned int k = 0; k < point.GetPointDimension(); k++) {
      parameters[k] = point[k];
    }
    // transform->SetParametersByValue( parameters );
    metricValue = metric->GetValue( parameters );
    std::cout << "At index " << index << " metric value is " << 
         metricValue << std::endl;
  }

  exit(0); // success
}

