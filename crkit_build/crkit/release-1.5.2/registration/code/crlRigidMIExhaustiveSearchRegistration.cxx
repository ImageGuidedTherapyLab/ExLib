/*
 * Copyright 2007 Simon Warfield simon.warfield@childrens.harvard.edu
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#ifdef WIN32
#include <missingFunctionsWindows.h>
#else
#include <sys/time.h>
#endif

#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>

#include <itkImageRegistrationMethod.h>
#include <itkMultiResolutionImageRegistrationMethod.h>
#include <itkMultiResolutionPyramidImageFilter.h>
#include <itkMeanSquaresImageToImageMetric.h>
#include <itkMattesMutualInformationImageToImageMetric.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkImage.h>
#include <itkOrientedImage.h>
#include <itkOrientImageFilter.h>

#include <itkVersorRigid3DTransform.h>
#include <itkCenteredTransformInitializer.h>
#include <itkCenteredVersorTransformInitializer.h>

#include <itkVersorRigid3DTransformOptimizer.h>
#include <itkExhaustiveOptimizer.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkResampleImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkExtractImageFilter.h>

#include <tclap/CmdLine.h>
#include "configuration.h"

//  The following section of code implements a Command observer
//  that will monitor the evolution of the registration process.
//
#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:

  CommandIterationUpdate() { 
    // These will be initialized after the observer object is registered 
    // with the optimization process.
    m_SearchImage = 0; 
    m_CurrentTransform = 0;
  };

  ~CommandIterationUpdate() { 
  }

public:
  // typedef itk::VersorRigid3DTransformOptimizer     OptimizerType;
  typedef itk::ExhaustiveOptimizer     OptimizerType;
  typedef   const OptimizerType   *    OptimizerPointer;
  typedef itk::VersorRigid3DTransform< double > TransformType;

  // Duplicate definition redone here
  typedef itk::Image< float, 6> SearchImageType;
  SearchImageType::Pointer m_SearchImage;
  TransformType::Pointer m_CurrentTransform;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      OptimizerPointer m_Optimizer =
        dynamic_cast< OptimizerPointer >( object );
      if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
        return;
        }
      std::cout << m_Optimizer->GetCurrentPosition() << std::endl;
// We could store the set of exhaustively evaluated positions in an image,
// or we could keep a running tally of the best position.
      // For now, we just print out the search values.
      for (unsigned int i = 0; 
                    i < m_Optimizer->GetCurrentPosition().GetSize(); i++) {
        std::cout << m_Optimizer->GetCurrentPosition()[i] << " " ;
      }
      std::cout << m_Optimizer->GetCurrentValue() << std::endl;
// This doesn't let me know what the fixed parameters of the initial transform
// are.
      if (m_CurrentTransform) {
        m_CurrentTransform->SetParameters( m_Optimizer->GetCurrentPosition() );
        std::cout << "The current transform is :";
        m_CurrentTransform->Print( std::cout );
        std::cout << std::endl;
      }

      SearchImageType::IndexType m_SearchImageIndex;
      for (unsigned int i = 0; i < m_SearchImage->GetImageDimension(); i++) {
        m_SearchImageIndex[i] = static_cast<long int>( 
                                    m_Optimizer->GetCurrentIndex()[i] );
      }
      m_SearchImage->SetPixel( m_SearchImageIndex, 
                               m_Optimizer->GetCurrentValue() );
    }
};


int main( int argc, char *argv[] )
{
  unsigned int numberOfBins = 32; // 0 is the default == 32 bins
  unsigned int numberOfMetricSamples = 0; // 0 is the default == useAllSamples.
  unsigned int numberOfMultiresolutionLevels = 1; // (good default 3)
  //unsigned int numberIterationsOptimizer = 200; // (good default 200)
  unsigned int transformInitializationScheme = 0; // 0 - moments
                                               // 1 - geometry
                                               // 2 - identity
                                               // 3 - from file
  //float maxOptimizerStepLength = 0.2; // (good default between 0.2 - 4.0)
  //float minOptimizerStepLength = 0.0001; // (good default 0.0001)
  std::string *fixedImageFile = new std::string("");
  std::string *movingImageFile = new std::string("");
  std::string *outputImageFile = new std::string("");
  std::string *outputTransformFile = new std::string("");
  std::string *searchImageFile = new std::string("searchImage.nrrd");

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',
                         CRKIT_VERSION_STRING );
      TCLAP::UnlabeledValueArg<std::string> fixedImageArg("fixedImage","Fixed Image File Name",true,"","fixed image file name",cmd);
      TCLAP::UnlabeledValueArg<std::string> movingImageArg("movingImage","Moving Image File Name",true,"","moving image file name",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputImageArg("outputImage","Output File Name",true,"","output image file name",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputTransformArg("outputTransform","Output File for Transform",true,"","transform file name",cmd);
      TCLAP::ValueArg<unsigned int> binsArg("b","numberOfBins","Number of Bins",false,numberOfBins,"bins (uint)",cmd);
      TCLAP::ValueArg<unsigned int> samplesArg("s","numberOfMetricSamples","Number of Metric Samples",false,numberOfMetricSamples,"samples (uint)",cmd);
      TCLAP::ValueArg<unsigned int> pyramidArg("p","pyramidLevels","Number of Multiresolution Levels",false,numberOfMultiresolutionLevels,"pyramid levels (uint)",cmd);
      TCLAP::ValueArg<unsigned int> transformInitializationArg("t","transformInitializationScheme","Mechanism for initializing the transform",false,transformInitializationScheme,"transform initialization scheme (uint : 0,1,2,3)",cmd);
      TCLAP::ValueArg<std::string> searchImageFileArg( "i", "searchImageFile",
       "File to write the metric search into", false, 
        "searchImage.nrrd", "filename (string)", cmd);

      cmd.parse(argc,argv);

      if (binsArg.isSet()) numberOfBins = binsArg.getValue();
      if (samplesArg.isSet()) numberOfMetricSamples = samplesArg.getValue();
      if (pyramidArg.isSet()) numberOfMultiresolutionLevels = pyramidArg.getValue();
      if (transformInitializationArg.isSet()) transformInitializationScheme =
               transformInitializationArg.getValue();
      if (fixedImageArg.isSet()) fixedImageFile = new std::string(fixedImageArg.getValue());
      if (movingImageArg.isSet()) movingImageFile = new std::string(movingImageArg.getValue());
      if (outputImageArg.isSet()) outputImageFile = new std::string(outputImageArg.getValue());
      if (outputTransformArg.isSet()) outputTransformFile = new std::string(outputTransformArg.getValue());
      if (searchImageFileArg.isSet()) searchImageFile = new std::string(searchImageFileArg.getValue());
    }
  catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(1);
    }

  const    unsigned int    Dimension = 3;
  typedef  float           PixelType;

  typedef itk::OrientedImage< PixelType, Dimension >  FixedImageType;
  typedef itk::OrientedImage< PixelType, Dimension >  MovingImageType;

  typedef itk::VersorRigid3DTransform< double > TransformType;

  typedef itk::VersorRigid3DTransformOptimizer           VersorOptimizerType;
  typedef itk::ExhaustiveOptimizer           OptimizerType;

  // Define an image type to store the rigid search results.
  typedef itk::Image< float, 6> SearchImageType;
  typedef itk::ImageFileWriter< SearchImageType >  SearchWriterType;

  typedef itk::MeanSquaresImageToImageMetric<
                                    FixedImageType,
                                    MovingImageType >
      MeanSquaresMetricType;

  typedef itk::MattesMutualInformationImageToImageMetric<
                                    FixedImageType,
                                    MovingImageType >
      MattesMutualInformationMetricType;

  typedef MattesMutualInformationMetricType MetricType;

  typedef itk:: LinearInterpolateImageFunction<
                                    MovingImageType,
                                    double          >    LinearInterpolatorType;

  typedef itk:: NearestNeighborInterpolateImageFunction<
                                    MovingImageType,
                                    double          >
      NearestNeighborInterpolatorType;

  typedef LinearInterpolatorType InterpolatorType;
  // typedef NearestNeighborInterpolatorType InterpolatorType;

  /*
   * A MultiResolutionPyramidImageFilter is used to create a pyramid of
   * downsampled images. The size of each downsampled image is specified
   * by the user in the form of a schedule of shrink factors. A default
   * schedule can be used, where the shrink changes by a factor 2 at each
   * level.
   */
  typedef itk::MultiResolutionImageRegistrationMethod<
                                    FixedImageType,
                                    MovingImageType >    RegistrationType;

  typedef   float     InternalPixelType;
  typedef itk::OrientedImage< InternalPixelType, Dimension > InternalImageType;

  typedef itk::MultiResolutionPyramidImageFilter<
                                    InternalImageType,
                                    InternalImageType >   FixedImagePyramidType;

  typedef itk::MultiResolutionPyramidImageFilter<
                                    InternalImageType,
                                    InternalImageType >  MovingImagePyramidType;

  MetricType::Pointer         metric        = MetricType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  FixedImagePyramidType::Pointer fixedImagePyramid =
    FixedImagePyramidType::New();
  MovingImagePyramidType::Pointer movingImagePyramid =
    MovingImagePyramidType::New();

  TransformType::Pointer  transform = TransformType::New();

  registration->SetOptimizer(     optimizer     );
  registration->SetTransform( transform );
  registration->SetMetric(        metric        );
  registration->SetInterpolator(  interpolator  );
  registration->SetFixedImagePyramid(  fixedImagePyramid );
  registration->SetMovingImagePyramid(  movingImagePyramid );

  // Set up key parameters for the metric
  metric->SetNumberOfHistogramBins( numberOfBins );
  if (numberOfMetricSamples == 0) {
    // samples every voxel - not so long on modern hardware, and much more
    // robust than the sampling strategy.
    metric->UseAllPixelsOn();
  } else {
    metric->SetNumberOfSpatialSamples( numberOfMetricSamples );
  }

  //  Given that the Mattes Mutual Information metric uses a random iterator in
  //  order to collect the samples from the images, it is usually convenient to
  //  initialize the seed of the random number generator.
  // metric->ReinitializeSeed( 76926294 );
  // select a seed based on the clock, to microsecond precision.
  unsigned int seed = 0;
  struct timeval preciseTimeStart;
  gettimeofday(&preciseTimeStart, NULL);
  seed = (unsigned int)(preciseTimeStart.tv_usec & 0xFFFFFFFF);
  // std::cout << "Mattes MI random seed is " << seed << std::endl;
  metric->ReinitializeSeed( seed );

  //  The transform object is constructed below and passed to the registration
  //  method.

  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  itk::OrientImageFilter<FixedImageType,FixedImageType>::Pointer fixedImageOrienter = itk::OrientImageFilter<FixedImageType,FixedImageType>::New();

  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
  itk::OrientImageFilter<MovingImageType,MovingImageType>::Pointer movingImageOrienter = itk::OrientImageFilter<MovingImageType,MovingImageType>::New();

  fixedImageReader->SetFileName(  fixedImageFile->c_str() );
  fixedImageOrienter->UseImageDirectionOn();
  fixedImageOrienter->SetDesiredCoordinateOrientationToAxial();
  fixedImageOrienter->SetInput( fixedImageReader->GetOutput() );
  fixedImageOrienter->Update();

  movingImageReader->SetFileName( movingImageFile->c_str() );
  movingImageOrienter->UseImageDirectionOn();
  movingImageOrienter->SetDesiredCoordinateOrientationToAxial();
  movingImageOrienter->SetInput( movingImageReader->GetOutput() );
  movingImageOrienter->Update();

  // Put the data through an orient image filter so that they are appropriately
  // corresponding when we run the registration.

  registration->SetFixedImage(    fixedImageOrienter->GetOutput()    );
  registration->SetMovingImage(   movingImageOrienter->GetOutput()   );

  registration->SetFixedImageRegion(
     fixedImageOrienter->GetOutput()->GetBufferedRegion() );

// Using a versor transform initializer
  typedef itk::CenteredVersorTransformInitializer<
                                             FixedImageType,
                                             MovingImageType
                                                 >  TransformInitializerType;

  TransformInitializerType::Pointer initializer =
                                          TransformInitializerType::New();

  //  The initializer is now connected to the transform and to the fixed and
  //  moving images.

  initializer->SetTransform(   transform );
  initializer->SetFixedImage(  fixedImageOrienter->GetOutput() );
  initializer->SetMovingImage( movingImageOrienter->GetOutput() );

  // This initializer is designed to utilize either the center of mass
  // computed from the images, (there has long been a bug that causes 
  // initialization of rotations to fail in MomentsOn() mode.)

    if (transformInitializationScheme == 0) {  // moments, default
      initializer->MomentsOn();
      initializer->InitializeTransform();
    } else if (transformInitializationScheme == 1) { // geometry centers
      initializer->GeometryOn();
      initializer->InitializeTransform();
    } else if (transformInitializationScheme == 2) { // identity transform
      ; // Here we simply don't call the initializer, so the transform
        // is used with its default initialization.
    } else if (transformInitializationScheme == 3) { // read from a file
      std::cout << "Warning: Transform initialization from a file "
         << "not yet implemented." << std::endl;
    }

  /* For debugging uncomment this: */
     std::cout << "Initial transform from initializer is " 
               << transform << std::endl;

  //  We now pass the parameters of the current transform as the initial
  //  parameters to be used when the registration process starts.

  registration->SetInitialTransformParameters( transform->GetParameters() );

  //const double translationScale = 1.0 / 1000.0;

  // This allows us to specify search steps for the optimizer
  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );

/* Use of optimizerScales to modify the search grip:
 * If you wish to use different step lengths along each parametric axis,
 * you can use the SetScales() method. This accepts an array, each element
 * represents the number of subdivisions per step length. For instance scales
 * of [0 0 0 0.5 1 4 ] along with a step length of 2 will cause the optimizer
 * to search the metric space on a grid with x,y,z spacing of [1 2 8].
 */
  // These get reset below based on the voxel size and unit degree rotation.
  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = 1.0 ;
  optimizerScales[4] = 1.0 ;
  optimizerScales[5] = 1.0 ;

  optimizer->SetScales( optimizerScales );
  OptimizerType::StepsType steps( transform->GetNumberOfParameters() );
// For a unit magnitude quaternion (i.e. a versor ), ITK uses the
// parameters x,y,z and ensures unit magnitude by setting w.
//   The parameters x,y,z are parallel to the axis of rotation,
// and their magnitude = sin( theta / 2 ), where
// theta is the angle of rotation.
// Therefore, a 1 degree rotation around a basis vector can be achieved with
// a single parameter value set to be 0.0934159274.
// Recall 1 degree = 0.0174532925 radians .
//
// Steps + and - this amount, as well as computing step at 0, 
// so this is a total of 21 x 21 x 21 evaluations, with the parameter changing
// by a factor of 2 at each time.
//
  steps[0] = 0; // corresponds to +/- 0 steps of versor search at given size.
  steps[1] = 0;
  steps[2] = 0;
  steps[3] = 10; // 10 translation x steps.
  steps[4] = 10; // 10 translation y steps.
  steps[5] = 10; // 10 translation z steps.
  optimizer->SetNumberOfSteps( steps );

  // Make the minimum step length match the smallest voxel size of 
  // the fixed image.
  FixedImageType::SpacingType spacing = fixedImageOrienter->GetOutput()->GetSpacing();
  optimizer->SetStepLength( spacing[0] );
  for (unsigned int i = 1; i < spacing.GetNumberOfComponents(); i++) {
    if (spacing[i] < optimizer->GetStepLength()) {
      optimizer->SetStepLength( spacing[i] );
    }
  }

  // This is the step size we want for small rotations to be considered.
  double unitdegreerotation = 0.0934159274;

  // We set the translation scale to be a multiple of this, based on the
  // ratio of the smallest voxel size to the rotation increment.
  //   This way, we always step at least one voxel or one degree of rotation.
  if (0) { // 0 means we ignore rotations.
    optimizerScales[0] = spacing[0] / unitdegreerotation;
    optimizerScales[1] = spacing[1] / unitdegreerotation;
    optimizerScales[2] = spacing[2] / unitdegreerotation;
    if (optimizer->GetStepLength() > unitdegreerotation) {
      optimizer->SetStepLength( unitdegreerotation );
    }
  } else {
    std::cout << "Turning off rotation searches..." << std::endl;
    optimizerScales[0] = optimizerScales[1] = optimizerScales[2] = 0.0;
  }

std::cout << "Step Length is " << optimizer->GetStepLength() << std::endl;

  // Now alloc the image which stores these search results:
  SearchImageType::Pointer searchImage = SearchImageType::New(); 
  SearchImageType::IndexType searchImageIndex;
  SearchImageType::SizeType searchImageSize;
  SearchImageType::SpacingType searchImageSpacing;
  SearchImageType::RegionType region;
  for (unsigned int i = 0; i < searchImage->GetImageDimension(); i++) {
    searchImageSize[i] = steps[i]*2 + 1;
    searchImageIndex[i] = 0;
    searchImageSpacing[i] = optimizer->GetStepLength();
  }
  region.SetSize(  searchImageSize  );
  region.SetIndex( searchImageIndex );
  searchImage->SetRegions( region );
  searchImage->Allocate( );

  SearchImageType::PixelType initialValue = 0.0;
  searchImage->FillBuffer( initialValue );

  // optimizer->SetMaximumStepLength( maxOptimizerStepLength );
  // optimizer->SetMinimumStepLength( minOptimizerStepLength );
  // optimizer->SetNumberOfIterations( numberIterationsOptimizer );

  // Initialze the observer with the initial transform, so that we can
  // know how the fixed parameters (the center for rotations) has been
  // initialized.
  //
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );
  observer->m_SearchImage = searchImage;
  observer->m_CurrentTransform = TransformType::New();
  observer->m_CurrentTransform->SetCenter( transform->GetCenter() );
  observer->m_CurrentTransform->SetParameters( transform->GetParameters() );

  // Set the number of levels of the multiresolution pyramid.
  // FUTURE IMPROVEMENT: Allow generation of the multiresolution pyramid
  // on the basis of the physical extents of the data.
  registration->SetNumberOfLevels( numberOfMultiresolutionLevels );

  try
    {
    registration->StartRegistration();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return -1;
    }

  OptimizerType::ParametersType finalParameters =
                    optimizer->GetMinimumMetricValuePosition();

  const double versorX              = finalParameters[0];
  const double versorY              = finalParameters[1];
  const double versorZ              = finalParameters[2];
  const double finalTranslationX    = finalParameters[3];
  const double finalTranslationY    = finalParameters[4];
  const double finalTranslationZ    = finalParameters[5];

//  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  const unsigned int numberOfIterations = optimizer->GetMaximumNumberOfIterations();

  // To find the best value, we need to loop over the estimated metric function.
  const double bestValue = optimizer->GetMinimumMetricValue();

  /* For debugging: */
  std::cout << std::endl << std::endl;
  std::cout << "Result = " << std::endl;
  std::cout << " versor X      = " << versorX  << std::endl;
  std::cout << " versor Y      = " << versorY  << std::endl;
  std::cout << " versor Z      = " << versorZ  << std::endl;
  std::cout << " Translation X = " << finalTranslationX  << std::endl;
  std::cout << " Translation Y = " << finalTranslationY  << std::endl;
  std::cout << " Translation Z = " << finalTranslationZ  << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;
  

  transform->SetParameters( finalParameters );

  TransformType::MatrixType matrix = transform->GetRotationMatrix();
  TransformType::OffsetType offset = transform->GetOffset();

  typedef itk::ResampleImageFilter<
                            MovingImageType,
                            FixedImageType >    ResampleFilterType;

  TransformType::Pointer finalTransform = TransformType::New();
  finalTransform->SetCenter( transform->GetCenter() );
  finalTransform->SetParameters( finalParameters );

  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  InterpolatorType::Pointer   resampleInterpolator  = InterpolatorType::New();
  resampler->SetInterpolator( resampleInterpolator );

  resampler->SetTransform( finalTransform );
  resampler->SetInput( movingImageOrienter->GetOutput() );

  FixedImageType::Pointer fixedImage = fixedImageOrienter->GetOutput();

  /* Since we may be working with oriented images, ensure that the output
   * direction is set properly.
   */
  resampler->SetOutputDirection( fixedImage->GetDirection() );
  resampler->SetOutputParametersFromImage( fixedImage );
  resampler->SetDefaultPixelValue( 0.0 );

  typedef  float OutputPixelType;

  typedef itk::OrientedImage< OutputPixelType, Dimension > OutputImageType;

  typedef itk::CastImageFilter<
                        FixedImageType,
                        OutputImageType > CastFilterType;

  typedef itk::ImageFileWriter< OutputImageType >  WriterType;


  SearchWriterType::Pointer      searchImageWriter =  SearchWriterType::New();
  searchImageWriter->SetFileName( *searchImageFile );
  searchImageWriter->SetInput( searchImage );
  searchImageWriter->Update();

  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();

  writer->SetFileName( outputImageFile->c_str() );
  writer->UseCompressionOn();

  caster->SetInput( resampler->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();

  // Now write the transform to a file.
  itk::TransformFileWriter::Pointer trsfWriter;
  trsfWriter = itk::TransformFileWriter::New();
  trsfWriter->SetInput( transform );
  trsfWriter->SetFileName( outputTransformFile->c_str() );

  try
  {
    trsfWriter->Update();
  }
  catch( itk::ExceptionObject & excp)
  {
    std::cerr << "Error while saving the transform to file " 
              << (*outputTransformFile) << std::endl;
    std::cerr << excp << std::endl;
    return 0;
  }

  return 0;
}
