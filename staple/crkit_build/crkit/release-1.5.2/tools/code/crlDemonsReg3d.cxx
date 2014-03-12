
#include "itkImageFileReader.h" 
#include "itkImageFileWriter.h" 

#include "itkImageRegionIterator.h"

#include "itkDemonsRegistrationFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"


//  The following section of code implements a Command observer
//  that will monitor the evolution of the registration process.
//
  class CommandIterationUpdate : public itk::Command 
  {
  public:
    typedef  CommandIterationUpdate   Self;
    typedef  itk::Command             Superclass;
    typedef  itk::SmartPointer<CommandIterationUpdate>  Pointer;
    itkNewMacro( CommandIterationUpdate );
  protected:
    CommandIterationUpdate() {};

    const static unsigned int Dimension = 3;

    typedef itk::Image< float, Dimension> InternalImageType;
    typedef itk::Vector< float, Dimension >    VectorPixelType;
    typedef itk::Image<  VectorPixelType, Dimension > DeformationFieldType;

    typedef itk::DemonsRegistrationFilter<
                                InternalImageType,
                                InternalImageType,
                                DeformationFieldType>   RegistrationFilterType;

  public:

    void Execute(itk::Object *caller, const itk::EventObject & event)
      {
        Execute( (const itk::Object *)caller, event);
      }

    void Execute(const itk::Object * object, const itk::EventObject & event)
      {
         const RegistrationFilterType * filter = 
          dynamic_cast< const RegistrationFilterType * >( object );
        if( !(itk::IterationEvent().CheckEvent( &event )) )
          {
          return;
          }
        std::cout << filter->GetMetric() << std::endl;
      }
  };


int main( int argc, char *argv[] )
{
  if( argc < 6 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile movingImageFile ";
    std::cerr << " numberOfIterations standardDeviation ";
    std::cerr << " outputImageFile " << std::endl;
    std::cerr << " [outputDeformationFieldFile] " << std::endl;
    return 1;
    }

  const unsigned int Dimension = 3;
  typedef float PixelType;

  unsigned int numberOfIterations = atoi(argv[3]); // 50 is a good default
  float standardDeviation = atof(argv[4]);         // 1.0 is a good default

  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  typedef itk::Image< PixelType, Dimension >  MovingImageType;

  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

  FixedImageReaderType::Pointer fixedImageReader   = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName( argv[1] );
  movingImageReader->SetFileName( argv[2] );

  typedef float InternalPixelType;
  typedef itk::Image< InternalPixelType, Dimension > InternalImageType;
  typedef itk::CastImageFilter< FixedImageType, 
                                InternalImageType > FixedImageCasterType;
  typedef itk::CastImageFilter< MovingImageType, 
                                InternalImageType > MovingImageCasterType;

  FixedImageCasterType::Pointer fixedImageCaster   = FixedImageCasterType::New();
  MovingImageCasterType::Pointer movingImageCaster = MovingImageCasterType::New();

  fixedImageCaster->SetInput( fixedImageReader->GetOutput() );
  movingImageCaster->SetInput( movingImageReader->GetOutput() ); 

  // The demons algorithm relies on the assumption that pixels representing the
  // same homologous point on an object have the same intensity on both the
  // fixed and moving images to be registered. In this example, we will
  // preprocess the moving image to match the intensity between the images
  // using the \doxygen{HistogramMatchingImageFilter}. 
  //
  // \index{itk::HistogramMatchingImageFilter}
  //
  // The basic idea is to match the histograms of the two images at a user-specified number of quantile values. For robustness, the histograms are
  // matched so that the background pixels are excluded from both histograms.
  // For MR images, a simple procedure is to exclude all gray values that are
  // smaller than the mean gray value of the image.
  //

  typedef itk::HistogramMatchingImageFilter<
                                    InternalImageType,
                                    InternalImageType >   MatchingFilterType;
  MatchingFilterType::Pointer matcher = MatchingFilterType::New();

  // For this example, we set the moving image as the source or input image and
  // the fixed image as the reference image.
  //

  matcher->SetInput( movingImageCaster->GetOutput() );
  matcher->SetReferenceImage( fixedImageCaster->GetOutput() );


  // We then select the number of bins to represent the histograms and the
  // number of points or quantile values where the histogram is to be
  // matched.
  //

  matcher->SetNumberOfHistogramLevels( 1024 );
  matcher->SetNumberOfMatchPoints( 7 );


  // Simple background extraction is done by thresholding at the mean
  // intensity.
  //

  matcher->ThresholdAtMeanIntensityOn();


  // In the \doxygen{DemonsRegistrationFilter}, the deformation field is
  // represented as an image whose pixels are floating point vectors.
  //

  typedef itk::Vector< float, Dimension >    VectorPixelType;
  typedef itk::Image<  VectorPixelType, Dimension > DeformationFieldType;
  typedef itk::DemonsRegistrationFilter<
                                InternalImageType,
                                InternalImageType,
                                DeformationFieldType>   RegistrationFilterType;
  RegistrationFilterType::Pointer filter = RegistrationFilterType::New();


  // Create the Command observer and register it with the registration filter.
  //
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  filter->AddObserver( itk::IterationEvent(), observer );

  // The input fixed image is simply the output of the fixed image casting
  // filter.  The input moving image is the output of the histogram matching
  // filter.
  //

  filter->SetFixedImage( fixedImageCaster->GetOutput() );
  filter->SetMovingImage( matcher->GetOutput() );


  //
  // The demons registration filter has two parameters: the number of
  // iterations to be performed and the standard deviation of the Gaussian
  // smoothing kernel to be applied to the deformation field after each
  // iteration.

  filter->SetNumberOfIterations( numberOfIterations );
  filter->SetStandardDeviations( standardDeviation );

  // The registration algorithm is triggered by updating the filter. The
  // filter output is the computed deformation field.

  filter->Update();


  // The \doxygen{WarpImageFilter} can be used to warp the moving image with
  // the output deformation field. Like the \doxygen{ResampleImageFilter},
  // the WarpImageFilter requires the specification of the input image to be
  // resampled, an input image interpolator, and the output image spacing and
  // origin.
  //

  typedef itk::WarpImageFilter<
                          MovingImageType, 
                          MovingImageType,
                          DeformationFieldType  >     WarperType;
  typedef itk::LinearInterpolateImageFunction<
                                   MovingImageType,
                                   double          >  InterpolatorType;
  WarperType::Pointer warper = WarperType::New();
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

  warper->SetInput( movingImageReader->GetOutput() );
  warper->SetInterpolator( interpolator );
  warper->SetOutputSpacing( fixedImage->GetSpacing() );
  warper->SetOutputOrigin( fixedImage->GetOrigin() );

  // Unlike the ResampleImageFilter, the WarpImageFilter
  // warps or transform the input image with respect to the deformation field
  // represented by an image of vectors.  The resulting warped or resampled
  // image is written to file as per previous examples.
  //
  warper->SetDeformationField( filter->GetOutput() );

  typedef  float OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::CastImageFilter< 
                        MovingImageType,
                        OutputImageType > CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();

  writer->SetFileName( argv[5] );
  
  caster->SetInput( warper->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();

  // It may be also desirable to write the deformation field as an image of
  // vectors.  This can be done with the following code.
  //

  if( argc > 6 ) // if a fourth line argument has been provided...
    {

  typedef itk::ImageFileWriter< DeformationFieldType > FieldWriterType;
  FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
  fieldWriter->SetFileName( argv[6] );
  fieldWriter->SetInput( filter->GetOutput() );

  fieldWriter->Update();

  //
  // Note that the file format used for writing the deformation field must be
  // capable of representing multiple components per pixel. This is the case
  // for the MetaImage and VTK file formats for example.
  //

    }

  return 0;
}

