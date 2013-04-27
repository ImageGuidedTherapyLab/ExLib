#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <vector>


int main( int argc, char * argv [] )
{
  const   unsigned int   Dimension = 2;
   
  typedef unsigned short  PixelType;
  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef itk::ImageFileReader< ImageType >   ReaderType;

  ReaderType::Pointer reader = ReaderType::New();

  ImageType::Pointer inputImage = NULL;

  if ( (argc == 1) || (strncmp("makemeone", argv[1], strlen(argv[1])) == 0)) {
    ImageType::SpacingType m_Spacing;
    ImageType::SizeType m_Size = { {100, 100} };
    ImageType::IndexType m_Start = { {0, 0} };
    ImageType::PointType m_Origin;

    m_Origin[0] = -0.1; m_Origin[1] = -0.1;
    m_Spacing[0] = 1.0; m_Spacing[1] = 1.0;

    ImageType::RegionType region;
        region.SetSize( m_Size );
        region.SetIndex( m_Start );

    inputImage = ImageType::New();
    inputImage->SetOrigin(m_Origin);
    inputImage->SetSpacing(m_Spacing);
    inputImage->SetRegions( region );
    inputImage->Allocate( );
    inputImage->FillBuffer( 100 );

  } else {

std::cout << "Reading the image." << std::endl;

  reader->SetFileName( argv[1] );

  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

    inputImage = reader->GetOutput();
  }


  typedef itk::LinearInterpolateImageFunction< ImageType > InterpolatorType;
  typedef itk::NearestNeighborInterpolateImageFunction< ImageType > NearestInterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  NearestInterpolatorType::Pointer nnInterpolator = NearestInterpolatorType::New();

  interpolator->SetInputImage( inputImage );
  nnInterpolator->SetInputImage( inputImage );

  typedef ImageType::PointType    PointType;
  PointType p;
  InterpolatorType::ContinuousIndexType index;

  p = inputImage->GetOrigin();

  // Let's compare the nearest neighbor and linear interpolator at the origin.
  if ( inputImage->TransformPhysicalPointToContinuousIndex( 
         p, index ) ) {
    std::cout << "Origin is inside the image." << std::endl;
  } else {
    std::cout << "Origin is outside the image." << std::endl;
  }
std::cout << "Origin is " << p << " and Continuous Index is " << index << std::endl;
  std::cout << "Nearest neighbor interpolator at origin: " << 
            nnInterpolator->EvaluateAtContinuousIndex( index ) << std::endl;

  std::cout << "Linear interpolator at origin: " << 
            interpolator->EvaluateAtContinuousIndex( index ) << std::endl;

  p[0] = p[0] - 0.1;
  p[1] = p[1] - 0.1;

  // Let's compare the nearest neighbor and linear interpolator near
  // the origin.
  if ( inputImage->TransformPhysicalPointToContinuousIndex( 
         p, index ) ) {
    std::cout << "Point is inside the image." << std::endl;
  } else {
    std::cout << "Point is outside the image." << std::endl;
  }
  std::cout << "Point is " << p << " and C. Index is " << index << std::endl;
  std::cout << "Nearest neighbor interpolator at point : " << 
            nnInterpolator->EvaluateAtContinuousIndex( index ) << std::endl;

  std::cout << "Linear interpolator at point : " << 
            interpolator->EvaluateAtContinuousIndex( index ) << std::endl;

  typedef std::vector< PointType >   ListOfPointsType;
  typedef std::vector< PixelType >   ListOfValuesType;

  ListOfPointsType listOfPoints;
  ListOfValuesType listOfValues;

  // Fill up somehow your list of listOfPoints
  p[0] = 0.0;  // X coordinate
  p[1] = 0.0;  // Y coordinate
  listOfPoints.push_back( p );

  p[0] = -0.6;  // X coordinate
  p[1] = -0.6;  // Y coordinate
  listOfPoints.push_back( p );

  p[0] = -0.45;  // X coordinate
  p[1] = -0.45;  // Y coordinate
  listOfPoints.push_back( p );

  p[0] = 0.45;  // X coordinate
  p[1] = 0.45;  // Y coordinate
  listOfPoints.push_back( p );

  p[0] = 99.4;  // X coordinate
  p[1] = 99.4;  // Y coordinate
  listOfPoints.push_back( p );

  p[0] = 99.5;  // X coordinate
  p[1] = 99.5;  // Y coordinate
  listOfPoints.push_back( p );

  p[0] = 99.6;  // X coordinate
  p[1] = 99.6;  // Y coordinate
  listOfPoints.push_back( p );

  // Todo: add more points, maybe reading them from a file ?

  ListOfPointsType::const_iterator pitr = listOfPoints.begin();
  ListOfPointsType::const_iterator pend = listOfPoints.end();

  while( pitr != pend )
    {
    const PointType & point = *pitr;
    bool isInside;
    isInside = inputImage->TransformPhysicalPointToContinuousIndex( point, index );
    if (!isInside) {
      std::cout << "Point " << point << " is NOT inside the image, " <<
         " so we skip this point." <<std::endl;
      ++pitr;
      continue;
    }
    PixelType pixelValue = static_cast<PixelType>( 
            interpolator->EvaluateAtContinuousIndex( index ) );
    listOfValues.push_back( pixelValue );
    std::cout << point << " = " << pixelValue << std::endl;
    ++pitr;
    }
  

  // At this point the listOfValues std::vector has the list
  // of pixel values at the locations of the list of points.

  return EXIT_SUCCESS;
}
