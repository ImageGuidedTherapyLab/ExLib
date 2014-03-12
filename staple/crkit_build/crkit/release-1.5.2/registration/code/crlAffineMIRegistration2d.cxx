/*
 * Copyright 2008 Children's Hospital Boston
 * Computational Radiology Laboratory
 * Simon Keith Warfield, Ph.D.
 * simon.warfield@childrens.harvard.edu
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
#include <itkImage.h>
#include <itkCastImageFilter.h>
// The 2D image classes don't actually support direction cosines properly.

#include <itkMeanSquaresImageToImageMetric.h>
#include <itkMattesMutualInformationImageToImageMetric.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>
#include <itkNormalizedMutualInformationHistogramImageToImageMetric.h>
#include <itkMatchCardinalityImageToImageMetric.h>

#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>

#include <itkCenteredTransformInitializer.h>

#include <itkAffineTransform.h>
#include <itkSimilarity2DTransform.h>
#include <itkResampleImageFilter.h>

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
  typedef itk::PowellOptimizer     PowellOptimizerType;
  typedef PowellOptimizerType      OptimizerType;

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
  unsigned int numberOfMultiresolutionLevels = 5; // (good default 5)
  unsigned int numberOptimizerIterations = 200; // (good default 200)
  float PowellStepLength = 0.01; // (good default is 0.01)
  float PowellStepTolerance = 0.001; // (good default 0.001)
  std::string *fixedImageFile = new std::string("");
  std::string *movingImageFile = new std::string("");
  std::string *outputImageFile = new std::string("");
  std::string *outputTransformFile = new std::string("");
  std::string *initialTransformFile = new std::string("");
  unsigned int transformInitializationScheme = 0; // 0 - moments
                                               // 1 - geometry
                                               // 2 - identity
                                               // 3 - from file


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
      // Only allow this to be supplied through a long argument.
      TCLAP::ValueArg<float> PowellStepLengthArg("","PowellStepLength","Powell Optimizer Step Length",false,PowellStepLength,"floating point number",cmd);
      // Only allow this to be supplied through a long argument.
      TCLAP::ValueArg<float> PowellStepToleranceArg("","PowellStepTolerance","Powell Optimizer Step Tolerance",false,PowellStepTolerance,"floating point number",cmd);
      TCLAP::ValueArg<unsigned int> transformInitializationArg("t","transformInitializationScheme","Mechanism for initializing the transform",false,transformInitializationScheme,"transform initialization scheme (uint : 0 moments, 1 geometrical center,2 identity transform)",cmd);


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
      if (transformInitializationArg.isSet()) {
         transformInitializationScheme =
               transformInitializationArg.getValue();
      } else if (initialTransformArg.isSet()) {
        initialTransformFile = new std::string(initialTransformArg.getValue());
      }
    }
  catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(1);
    }

  const    unsigned int    Dimension = 2;
  typedef  float           PixelType;

  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  typedef itk::Image< PixelType, Dimension >  MovingImageType;

  // I want to be able to select at run time to do one or the other of these.
  typedef itk::AffineTransform< double, Dimension > AffineTransformType;
  typedef itk::Similarity2DTransform< double > SimilarityTransformType;

  typedef AffineTransformType TransformType;

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
  typedef itk::Image< InternalPixelType, Dimension > InternalImageType;

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
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName(  fixedImageFile->c_str() );

  movingImageReader->SetFileName( movingImageFile->c_str() );

  // Put the data through an orient image filter so that they are appropriately
  // corresponding when we run the registration.

  registration->SetFixedImage(    fixedImageReader->GetOutput()    );
  registration->SetMovingImage(   movingImageReader->GetOutput()   );
  fixedImageReader->Update();
  movingImageReader->Update();

  registration->SetFixedImageRegion(
     fixedImageReader->GetOutput()->GetBufferedRegion() );

  typedef itk::CenteredTransformInitializer<
			AffineTransformType,
                                             FixedImageType,
                                             MovingImageType
                                                 >  TransformInitializerType;

  TransformInitializerType::Pointer initializer =
                                          TransformInitializerType::New();

  // We want to support either estimating the initial transform from the
  // data, or else starting with a specified initial transform.
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

    typedef itk::Rigid2DTransform< double > Rigid2DTransformType;
    Rigid2DTransformType::Pointer trsf_readRigid2D = Rigid2DTransformType::New();
    AffineTransformType::Pointer trsf_readAffine = AffineTransformType::New();
    SimilarityTransformType::Pointer trsf_readSim = SimilarityTransformType::New();

    if (transforms->size() == 1) {
      itk::TransformFileReader::TransformListType::const_iterator trsfit =
           transforms->begin();
      if (!strcmp((*trsfit)->GetNameOfClass(),"Rigid2DTransform"))
      {
        trsf_readRigid2D = static_cast<Rigid2DTransformType*>((*trsfit).GetPointer());
        std::cout << "Loaded a rigid 2D transform." << std::endl;
        trsf_readRigid2D->Print(std::cout);
        // This is the recommended procedure for initializing when the
        // transforms have different mappings of actual parameter to index
        // in the array of parameters.
        transform->SetFixedParameters( trsf_readRigid2D->GetFixedParameters() );
        transform->SetMatrix(      trsf_readRigid2D->GetMatrix()      );
      } else if ( (!strcmp((*trsfit)->GetNameOfClass(),"AffineTransform")) ) {
        trsf_readAffine = static_cast<AffineTransformType *>(
                                                  (*trsfit).GetPointer());
        std::cout << "Loaded an AffineTransform." << std::endl;
        trsf_readAffine->Print(std::cout);
        transform->SetFixedParameters( trsf_readAffine->GetFixedParameters() );
        transform->SetParameters(      trsf_readAffine->GetParameters()      );
      } else if ( (!strcmp((*trsfit)->GetNameOfClass(),"Similarity2DTransform")) ) {
        trsf_readSim = static_cast<SimilarityTransformType *>(
                                                  (*trsfit).GetPointer());
        std::cout << "Loaded a Similarity2DTransform." << std::endl;
        trsf_readSim->Print(std::cout);
        transform->SetFixedParameters( trsf_readSim->GetFixedParameters() );
        transform->SetMatrix(      trsf_readSim->GetMatrix()      );
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
    if (transformInitializationScheme == 0) {  // moments, default

      std::cout << "Initializing using moments." << std::endl;

      //  The initializer is now connected to the transform and to the fixed and
      //  moving images.
      initializer->SetTransform(   transform );
      initializer->SetFixedImage(  fixedImageReader->GetOutput() );
      initializer->SetMovingImage( movingImageReader->GetOutput() );

      initializer->MomentsOn();

      // This appears to initialize just the center and the translation.
      // Although principal axes are computed, they are not used to initialize
      // the versor.
      initializer->InitializeTransform();

    } else if (transformInitializationScheme == 1) { // geometry centers
      // Implementation of geometrical initialization
      std::cout << "Initializing using geometry." << std::endl;
      const FixedImageType::SpacingType&
      fixedSpacing = fixedImageReader->GetOutput()->GetSpacing();

      const FixedImageType::PointType&
      fixedOrigin  = fixedImageReader->GetOutput()->GetOrigin();

      FixedImageType::SizeType fixedSize  =
          fixedImageReader->GetOutput()->GetLargestPossibleRegion().GetSize();

      TransformType::InputPointType centerFixed;
  centerFixed[0] = fixedOrigin[0] + fixedSpacing[0] * (fixedSize[0] - 1) / 2.0;
  centerFixed[1] = fixedOrigin[1] + fixedSpacing[1] * (fixedSize[1] - 1) / 2.0;

      const MovingImageType::SpacingType&
      movingSpacing = movingImageReader->GetOutput()->GetSpacing();

      const MovingImageType::PointType&
      movingOrigin  = movingImageReader->GetOutput()->GetOrigin();

      MovingImageType::SizeType movingSize  =
          movingImageReader->GetOutput()->GetLargestPossibleRegion().GetSize();

      TransformType::InputPointType centerMoving;
  centerMoving[0] = movingOrigin[0] + movingSpacing[0] * (movingSize[0]-1)/ 2.0;
  centerMoving[1] = movingOrigin[1] + movingSpacing[1] * (movingSize[1]-1)/ 2.0;

      transform->SetCenter( centerFixed );
      transform->SetTranslation( centerMoving - centerFixed );
    } else if  (transformInitializationScheme == 2) { // identity transform
        ; // Here we simply don't call the initializer, so the transform
          // is used with its default initialization.
      std::cout << "Initializing using identity transform." << std::endl;
    }

    //  We now pass the parameters of the current transform as the initial
    //  parameters to be used when the registration process starts.
    //
    registration->SetInitialTransformParameters( transform->GetParameters() );

  } 
  std::cout << "Initial transform from initializer is " <<
         transform << std::endl;

  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
  const double translationScale = 1.0 / 1000.0;

  for (unsigned int i = 0; i < transform->GetNumberOfParameters(); i++) {
    optimizerScales[i] = 1.0;
  }
  // full affine transform
  // could check transform->GetNameOfClass() to be sure.
  // The last Dimension parameters are the translation parameters.
  for (unsigned int i = 0; i < Dimension; i++) {
    optimizerScales[transform->GetNumberOfParameters()-1-i] = translationScale;
  }

  optimizer->SetScales( optimizerScales );

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

  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();

  const double bestValue = optimizer->GetValue();

  // Print out results
  //
  std::cout << std::endl;
  std::cout << "Result = " << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;

  std::cout << "Final transform computed : " << transform << std::endl;

  typedef itk::ResampleImageFilter<
                            MovingImageType,
                            FixedImageType >    ResampleFilterType;

  TransformType::Pointer finalTransform = TransformType::New();

  finalTransform->SetCenter( transform->GetCenter() );

  finalTransform->SetParameters( transform->GetParameters() );

  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  InterpolatorType::Pointer   resampleInterpolator  = InterpolatorType::New();
  resampler->SetInterpolator( resampleInterpolator );

  resampler->SetTransform( finalTransform );
  resampler->SetInput( movingImageReader->GetOutput() );

  FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

  /* Since we are working with oriented images, ensure that the output
   * direction is set properly.
   */
  resampler->SetOutputDirection(  fixedImage->GetDirection() );
  resampler->SetOutputStartIndex( fixedImage->GetLargestPossibleRegion().GetIndex());
  resampler->SetOutputParametersFromImage( fixedImage );
  resampler->SetDefaultPixelValue( 0.0 );

  typedef  float OutputPixelType;

  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

  typedef itk::CastImageFilter<
                        FixedImageType,
                        OutputImageType > CastFilterType;

  typedef itk::ImageFileWriter< OutputImageType >  WriterType;


  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();

  writer->SetFileName( outputImageFile->c_str() );
  writer->UseCompressionOn( );

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
    std::cerr << "Error while saving the transform to file " << (*outputTransformFile) <<
    std::endl;
    std::cerr << excp << std::endl;
    return 0;
  }

  return 0;
}
