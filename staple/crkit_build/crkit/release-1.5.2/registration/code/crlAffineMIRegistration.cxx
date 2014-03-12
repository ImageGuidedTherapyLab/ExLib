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


/*
 * This would be improved if there was simply one registration program
 * that could handle different types of transform.
 *
 * Missing:
 *   add a capacity to handle different types of transforms:
 *      translation only
 *      versor 3D     (        versor + translation)
 *      similarity 3D (scale + versor + translation)
 *      scale skew versor 3d transform -
 *        parameters: 0-2 versor, 3-5 translation, 6-8 Scale, 9-14 Skew
 *
 * We need to :
 *   add test data for a regression test
 *
 */

#ifdef WIN32
#include <missingFunctionsWindows.h>
#else
#include <sys/time.h>
#endif

#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegistrationMethod.h>
#include <itkMultiResolutionImageRegistrationMethod.h>
#include <itkMultiResolutionPyramidImageFilter.h>
#include <itkOrientedImage.h>
#include <itkOrientImageFilter.h>
#include <itkCastImageFilter.h>

#include <itkMeanSquaresImageToImageMetric.h>
#include <itkMattesMutualInformationImageToImageMetric.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>
#include <itkNormalizedMutualInformationHistogramImageToImageMetric.h>
#include <itkMatchCardinalityImageToImageMetric.h>

#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include <itkVersorRigid3DTransform.h>
#include <itkCenteredTransformInitializer.h>
#include <itkCenteredVersorTransformInitializer.h>
#include <itkScaleSkewVersor3DTransform.h>
#include <itkAffineTransform.h>
#include <itkResampleImageFilter.h>

#include <itkVersorRigid3DTransformOptimizer.h>
#include <itkPowellOptimizer.h>

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
  CommandIterationUpdate() {};
public:
  typedef itk::VersorRigid3DTransformOptimizer     VersorRigid3DOptimizerType;
  typedef itk::PowellOptimizer     PowellOptimizerType;
  typedef PowellOptimizerType OptimizerType;

  typedef   const OptimizerType   *    OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      OptimizerPointer optimizer =
        dynamic_cast< OptimizerPointer >( object );
      if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
        return;
        }
      std::cout << optimizer->GetCurrentIteration() << "   ";
      std::cout << optimizer->GetValue() << "   ";
      std::cout << optimizer->GetCurrentPosition() << std::endl;
    }
};


int main( int argc, char *argv[] )
{

  unsigned int numberOfBins = 32; // 0 is the default == 32 bins
  unsigned int numberOfMetricSamples = 0; // 0 is the default == useAllSamples.
  unsigned int numberOfMultiresolutionLevels = 3; // (good default 3)
  unsigned int numberOptimizerIterations = 200; // (good default 200)
  float PowellStepLength = 0.01; // (good default is 0.01)
  float PowellStepTolerance = 0.001; // (good default 0.001)
  std::string *fixedImageFile = new std::string("");
  std::string *movingImageFile = new std::string("");
  std::string *outputImageFile = new std::string("");
  std::string *outputTransformFile = new std::string("");
  std::string *initialTransformFile = new std::string("");

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',
                          CRKIT_VERSION_STRING );
      TCLAP::UnlabeledValueArg<std::string> fixedImageArg("fixedImage","Fixed Image File Name",true,"","fixed image file name",cmd);
      TCLAP::UnlabeledValueArg<std::string> movingImageArg("movingImage","Moving Image File Name",true,"","moving image file name",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputImageArg("outputImage","Output File Name",true,"","output image file name",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputTransformArg("outputTransform","Output File for Transform",true,"","transform file name",cmd);
      TCLAP::UnlabeledValueArg<std::string> initialTransformArg("initialTransform","Optional Initial Rigid Body Transform File",false,"","optional transform file name",cmd);
      TCLAP::ValueArg<unsigned int> binsArg("b","numberOfBins","Number of Bins",false,numberOfBins,"unsigned integer",cmd);
      TCLAP::ValueArg<unsigned int> samplesArg("s","numberOfMetricSamples","Number of Metric Samples",false,numberOfMetricSamples,"unsigned integer",cmd);
      TCLAP::ValueArg<unsigned int> pyramidArg("p","pyramidLevels","Number of Multiresolution Levels",false,numberOfMultiresolutionLevels,"unsigned integer",cmd);
      TCLAP::ValueArg<unsigned int> iterationsArg("i","numberOptimizerIterations","Number of Optimizer Iterations",false,numberOptimizerIterations,"unsigned integer",cmd);
      TCLAP::ValueArg<float> PowellStepLengthArg("l","PowellStepLength","Powell Optimizer Step Length",false,PowellStepLength,"floating point number",cmd);
      TCLAP::ValueArg<float> PowellStepToleranceArg("t","PowellStepTolerance","Powell Optimizer Step Tolerance",false,PowellStepTolerance,"floating point number",cmd);


      cmd.parse(argc,argv);

      if (binsArg.isSet()) numberOfBins = binsArg.getValue();
      if (samplesArg.isSet()) numberOfMetricSamples = samplesArg.getValue();
      if (pyramidArg.isSet()) numberOfMultiresolutionLevels = pyramidArg.getValue();
      if (iterationsArg.isSet()) numberOptimizerIterations = iterationsArg.getValue();
      if (PowellStepLengthArg.isSet()) PowellStepLength = PowellStepLengthArg.getValue();
      if (PowellStepToleranceArg.isSet()) PowellStepTolerance = PowellStepToleranceArg.getValue();
      if (fixedImageArg.isSet()) fixedImageFile = new std::string(fixedImageArg.getValue());
      if (movingImageArg.isSet()) movingImageFile = new std::string(movingImageArg.getValue());
      if (outputImageArg.isSet()) outputImageFile = new std::string(outputImageArg.getValue());
      if (outputTransformArg.isSet()) outputTransformFile = new std::string(outputTransformArg.getValue());
      if (initialTransformArg.isSet()) initialTransformFile = new std::string(initialTransformArg.getValue());
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


  // I want to be able to select at run time to do one or the other of these.
  typedef itk::VersorRigid3DTransform< double > VersorRigid3DTransformType;
  typedef itk::ScaleSkewVersor3DTransform< double > ScaleSkewVersor3DTransformType;
  typedef ScaleSkewVersor3DTransformType TransformType;

  typedef itk::VersorRigid3DTransformOptimizer     VersorRigid3DOptimizerType;
  typedef itk::PowellOptimizer     PowellOptimizerType;
  typedef PowellOptimizerType OptimizerType;

  typedef itk::ImageToImageMetric<FixedImageType, MovingImageType >
                 MetricBaseType;

  typedef itk::MeanSquaresImageToImageMetric<
                                    FixedImageType,
                                    MovingImageType >
      MeanSquaresMetricType;

  typedef itk::MatchCardinalityImageToImageMetric<
                                    FixedImageType,
                                    MovingImageType >
      MatchCardinalityMetricType;

  typedef itk::MattesMutualInformationImageToImageMetric<
                                    FixedImageType,
                                    MovingImageType >
      MattesMutualInformationMetricType;

  typedef itk::NormalizedCorrelationImageToImageMetric<
                                    FixedImageType,
                                    MovingImageType >
      NormalizedCorrelationMetricType;

  typedef itk::NormalizedMutualInformationHistogramImageToImageMetric<
                                    FixedImageType,
                                    MovingImageType >
      NormalizedMutualInformationMetricType;

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

  MetricBaseType::Pointer         metricBase;
  MattesMutualInformationMetricType::Pointer metricMattesMI =
                                    MattesMutualInformationMetricType::New();

  NormalizedCorrelationMetricType::Pointer metricNormCorr =
                                    NormalizedCorrelationMetricType::New();

  MeanSquaresMetricType::Pointer metricSSD =
                                    MeanSquaresMetricType::New();

  NormalizedMutualInformationMetricType::Pointer metricNMI =
                               NormalizedMutualInformationMetricType::New();

  MatchCardinalityMetricType::Pointer metricMatchCardinality =
          MatchCardinalityMetricType::New();

  metricBase = static_cast< MattesMutualInformationMetricType *>(
                      (metricMattesMI).GetPointer() );

  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  FixedImagePyramidType::Pointer fixedImagePyramid =
    FixedImagePyramidType::New();
  MovingImagePyramidType::Pointer movingImagePyramid =
    MovingImagePyramidType::New();

  TransformType::Pointer  transform = TransformType::New();

  std::cout << "Optimizer will identify the best matching " <<
        transform->GetNameOfClass() << " transform." << std::endl;

  registration->SetOptimizer(     optimizer     );
  registration->SetTransform( transform );
  registration->SetMetric(        metricBase        );
  registration->SetInterpolator(  interpolator  );
  registration->SetFixedImagePyramid(  fixedImagePyramid );
  registration->SetMovingImagePyramid(  movingImagePyramid );

  // Set up key parameters for the metric
  // Determine how we initialize the metric based on what the metric
  // actually is.
  if (!strcmp(
        metricBase->GetNameOfClass(),"NormalizedCorrelationImageToImageMetric")
     ) {
     // The Normalized Correlation metric does not use random sampling.
     // Correlation is high at good overlap.
     // The class must be returning minus one times normalized correlation.
     // We should remove the mean signal level.
    metricNormCorr->SetSubtractMean( true );
  } else if (!strcmp(
    metricBase->GetNameOfClass(),"MattesMutualInformationImageToImageMetric")
     ) {
    metricMattesMI->SetNumberOfHistogramBins( numberOfBins );
    if (numberOfMetricSamples == 0) {
      // samples every voxel - not so long on modern hardware, and much more
      // robust than the sampling strategy.
      metricMattesMI->UseAllPixelsOn( );
    } else {
      // FUTURE : this should depend on the size of the image
      metricMattesMI->SetNumberOfSpatialSamples( numberOfMetricSamples );
      std::cout << "Using " << numberOfMetricSamples << " samples." << std::endl;
    }

    // select a seed based on the clock, to microsecond precision.
    unsigned int seed = 0;
    struct timeval preciseTimeStart;
    gettimeofday(&preciseTimeStart, NULL);
    seed = (unsigned int)(preciseTimeStart.tv_usec & 0xFFFFFFFF);

    // std::cout << "Mattes MI random seed is " << seed << std::endl;
    metricMattesMI->ReinitializeSeed( seed );
  } else if (!strcmp(
        metricBase->GetNameOfClass(),"MeanSquaresImageToImageMetric")
     ) {
      ; // No special initialization is needed for mean squares.
  } else if (!strcmp(
        metricBase->GetNameOfClass(),"MatchCardinalityImageToImageMetric")
     ) {
      ; // No special initialization is needed for match cardinality.
  } else if (!strcmp(
        metricBase->GetNameOfClass(),
             "NormalizedMutualInformationHistogramImageToImageMetric")
     ) {
     // Lot's of special initialization needed for this NMI implementation.
   NormalizedMutualInformationMetricType::HistogramType::SizeType histogramSize;
     histogramSize[0] = numberOfBins;  // 2 elements, since two scalar images
     histogramSize[1] = numberOfBins;
     metricNMI->SetHistogramSize( histogramSize );
     // The NMI metric does not use random selection of seeds.
     // Ensure the optimizer knows which way we want to go:
     optimizer->SetMaximize( true ); // Push joint entropy H(A,B) to be small,
                                     // so the metric gets larger.
  } else {
    std::cerr << "Error: Don't know how to initialize the metric of type " <<
      metricBase->GetNameOfClass() << "." << std::endl;
    return EXIT_FAILURE;
  }

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
  fixedImageOrienter->Update();
  movingImageOrienter->Update();

  registration->SetFixedImageRegion(
     fixedImageOrienter->GetOutput()->GetBufferedRegion() );

  // This code here deals with establishing the initial transform parameters.
  // We want to support either estimating the initial transform from the
  // data, or else starting with a specified initial transform.
  // Using a versor transform initializer
  typedef itk::CenteredVersorTransformInitializer<
                                             FixedImageType,
                                             MovingImageType
                                                 >  TransformInitializerType;

  TransformInitializerType::Pointer initializer =
                                          TransformInitializerType::New();

  if ((*initialTransformFile) != "") {
    itk::TransformFileReader::Pointer trsfReader;
    trsfReader = itk::TransformFileReader::New();
    trsfReader->SetFileName( initialTransformFile->c_str() );
    try {
      trsfReader->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
      std::cerr << "Error while reading transform file" << (*initialTransformFile) << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }

    typedef itk::TransformFileReader::TransformListType * TransformListType;
    TransformListType transforms = trsfReader->GetTransformList();
    std::cout << "Number of transforms = " << transforms->size() << std::endl;

    typedef itk::AffineTransform< double, 3> AffineTransformType;
    typedef itk::VersorRigid3DTransform< double > VersorRigid3DTransformType;
    VersorRigid3DTransformType::Pointer trsf_read =
                 VersorRigid3DTransformType::New();

    if (transforms->size() == 1) {
      itk::TransformFileReader::TransformListType::const_iterator trsfit =
           transforms->begin();
      if (!strcmp((*trsfit)->GetNameOfClass(),"VersorRigid3DTransform"))
      {
        trsf_read = static_cast<VersorRigid3DTransformType*>((*trsfit).GetPointer());
        std::cout << "Just loaded a versor rigid 3d transform." << std::endl;
        trsf_read->Print(std::cout);
        // This is the recommended procedure for initializing when the
        // transforms have different mappings of actual parameter to index
        // in the array of parameters.
        transform->SetCenter( trsf_read->GetCenter() );
        transform->SetTranslation( trsf_read->GetTranslation() );
        // API BUG: itk::VersorRigid3DTransform<double>' has no member
        //   named 'GetRotation'
        // API BUG: ScaleSkewVersor3DTransform has no SetVersor()
        transform->SetRotation( trsf_read->GetVersor() );
      } else if ( (!strcmp((*trsfit)->GetNameOfClass(),"ScaleSkewVersor3DTransform")) ) {
        trsf_read = static_cast<ScaleSkewVersor3DTransformType*>(
                                                  (*trsfit).GetPointer());
        std::cout << "Just loaded a ScaleSkewVersor3DTransform." << std::endl;
        trsf_read->Print(std::cout);
        // API BUG:
        // Amazingly, the ScaleSkewVersor3DTransform does not support the
        // critical SetMatrix API. Instead we do this:
        transform->SetFixedParameters( trsf_read->GetFixedParameters() );
        transform->SetParameters(      trsf_read->GetParameters()      );
      } else {
        std::cerr << "Can't initialize from transform of type " <<
                ( (*trsfit)->GetNameOfClass() ) << " ." << std::endl;
      return EXIT_FAILURE;
      }
    } else {
      std::cerr << "Can't initialize from " << transforms->size() <<
             " transforms." << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "Initializing using the transform read in. " << std::endl;

    // The transform reading code initializers the input transform.
    registration->SetInitialTransformParameters( transform->GetParameters() );

  } else {

      std::cout << "Initializing using moments." << std::endl;

      //  The initializer is now connected to the transform and to the fixed and
      //  moving images.
      initializer->SetTransform(   transform );
      initializer->SetFixedImage(  fixedImageOrienter->GetOutput() );
      initializer->SetMovingImage( movingImageOrienter->GetOutput() );

      initializer->MomentsOn();

      // This appears to initialize just the center and the translation.
      // Although principal axes are computed, they are not used to initialize
      // the versor.
      initializer->InitializeTransform();

      std::cout << "Initial transform from initializer is " <<
         transform << std::endl;

    //  We now pass the parameters of the current transform as the initial
    //  parameters to be used when the registration process starts.
    //
    registration->SetInitialTransformParameters( transform->GetParameters() );

  }

  std::cout << "transform as initialized : " <<
         transform << std::endl;

  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
  const double translationScale = 1.0 / 1000.0;

  for (unsigned int i = 0; i < transform->GetNumberOfParameters(); i++) {
    optimizerScales[i] = 1.0;
  }
  if (transform->GetNumberOfParameters() == 6 ) {
    // Versor3DRigidTransform
  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = translationScale;
  optimizerScales[4] = translationScale;
  optimizerScales[5] = translationScale;
  } else if (transform->GetNumberOfParameters() == 15) {
     // The translation parameters for a ScaleSkewVersor3DTransform
     // are in positions 3,4,5.
     // We could check transform->GetNameOfClass() to be sure.
    for (unsigned int i = 0; i < transform->GetNumberOfParameters(); i++) {
      optimizerScales[i] = 1.0;
    }
    optimizerScales[3] = translationScale;
    optimizerScales[4] = translationScale;
    optimizerScales[5] = translationScale;
  } else {
    std::cerr <<"Don't know how to scale optimizer for transform."<<std::endl;
  }

  optimizer->SetScales( optimizerScales );

  /* For Versor3DRigidTransformOptimizer:
  optimizer->SetMaximumStepLength( 0.2000  );
  optimizer->SetMinimumStepLength( 0.0001 );
  optimizer->SetNumberOfIterations( 200 );
  */

  /* For PowellOptimizer: */
  optimizer->SetMaximumIteration( numberOptimizerIterations ) ;
  optimizer->SetStepLength( PowellStepLength ) ;
  optimizer->SetStepTolerance( PowellStepTolerance ) ;

  // Create the Command observer and register it with the optimizer.
  //
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  // Set the number of levels of the multiresolution pyramid.
  // FUTURE: Expand options for generating the multiresolution pyramid
  // on the basis of the physical extents of the data.
  registration->SetNumberOfLevels(numberOfMultiresolutionLevels);

  try
    {
    registration->StartRegistration();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    std::cerr << "No result being generated." << std::endl;
    return -1;
    }

  OptimizerType::ParametersType finalParameters =
                    registration->GetLastTransformParameters();

  const double versorX              = finalParameters[0];
  const double versorY              = finalParameters[1];
  const double versorZ              = finalParameters[2];
  const double finalTranslationX    = finalParameters[3];
  const double finalTranslationY    = finalParameters[4];
  const double finalTranslationZ    = finalParameters[5];

  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();

  const double bestValue = optimizer->GetValue();

  // Print out results
  //
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

  std::cout << "Matrix = " << std::endl << matrix << std::endl;
  std::cout << "Offset = " << std::endl << offset << std::endl;

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

  /* Since we are working with oriented images, ensure that the output
   * direction is set properly.
   */
  resampler->SetOutputDirection(  fixedImage->GetDirection() );
  resampler->SetOutputParametersFromImage( fixedImage );
  resampler->SetDefaultPixelValue( 0.0 );

  typedef  float OutputPixelType;

  typedef itk::OrientedImage< OutputPixelType, Dimension > OutputImageType;

  typedef itk::CastImageFilter<
                        FixedImageType,
                        OutputImageType > CastFilterType;

  typedef itk::ImageFileWriter< OutputImageType >  WriterType;


  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();

  writer->SetFileName( outputImageFile->c_str() );

  caster->SetInput( resampler->GetOutput() );
  writer->SetInput( caster->GetOutput()   );

  try
  {
    writer->UseCompressionOn( );
    writer->Update();
  }
  catch( itk::ExceptionObject & excp)
  {
    std::cerr << "Error while saving the output image " <<
             outputImageFile << std::endl;
    std::cerr << excp << std::endl;
    //return 0;
  }

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
    std::cerr << "Error while saving the transform to file " << (*outputTransformFile) <<
    std::endl;
    std::cerr << excp << std::endl;
    return 0;
  }

  return 0;
}
