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
#include <itkNormalizedCorrelationImageToImageMetric.h>

#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkImage.h>

#include <itkRigid2DTransform.h>
#include <itkCenteredRigid2DTransform.h>
#include <itkCenteredTransformInitializer.h>
#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkPowellOptimizer.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkResampleImageFilter.h>

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
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
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
  unsigned int numberOfMultiresolutionLevels = 6; // (good default 3-6)
  unsigned int numberIterationsOptimizer = 200; // (good default 200)
  unsigned int transformInitializationScheme = 5; // hard coded for now.
                                               // 0 - moments
                                               // 1 - geometry
                                               // 2 - identity
                                               // 3 - from file
  float maxOptimizerStepLength = 0.1; // (good default between 0.1 - 4.0)
  float minOptimizerStepLength = 0.001; // (good default 0.001 - 0.0001)
  std::string *fixedImageFile = new std::string("");
  std::string *movingImageFile = new std::string("");
  std::string *outputImageFile = new std::string("");
  std::string *outputTransformFile = new std::string("");

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
      TCLAP::ValueArg<unsigned int> iterationsArg("i","numberOptimizerIterations","Number of Optimizer Iterations",false,numberIterationsOptimizer,"iterations (uint)",cmd);
      TCLAP::ValueArg<float> maxStepArg("m","maxStepOptimizer","Maximum Optimizer Step Length",false,maxOptimizerStepLength,"max step (float)",cmd);
      TCLAP::ValueArg<float> minStepArg("n","minStepOptimizer","Minimum Optimizer Step Length",false,minOptimizerStepLength,"min step (float)",cmd);
      TCLAP::ValueArg<unsigned int> transformInitializationArg("t","transformInitializationScheme","Mechanism for initializing the transform",false,transformInitializationScheme,"transform initialzation scheme (uint : 0,1,2,3)",cmd);

      cmd.parse(argc,argv);

      if (binsArg.isSet()) numberOfBins = binsArg.getValue();
      if (samplesArg.isSet()) numberOfMetricSamples = samplesArg.getValue();
      if (pyramidArg.isSet()) numberOfMultiresolutionLevels = pyramidArg.getValue();
      if (iterationsArg.isSet()) numberIterationsOptimizer = iterationsArg.getValue();
      if (maxStepArg.isSet()) maxOptimizerStepLength = maxStepArg.getValue();
      if (minStepArg.isSet()) minOptimizerStepLength = minStepArg.getValue();
      if (transformInitializationArg.isSet()) transformInitializationScheme =
               transformInitializationArg.getValue();
      if (fixedImageArg.isSet()) fixedImageFile = new std::string(fixedImageArg.getValue());
      if (movingImageArg.isSet()) movingImageFile = new std::string(movingImageArg.getValue());
      if (outputImageArg.isSet()) outputImageFile = new std::string(outputImageArg.getValue());
      if (outputTransformArg.isSet()) outputTransformFile = new std::string(outputTransformArg.getValue());
    }
  catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(1);
    }

  const    unsigned int    Dimension = 2;
  typedef  float           PixelType;

  // Not an oriented image in 2D
  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  // Not an oriented image in 2D
  typedef itk::Image< PixelType, Dimension >  MovingImageType;

  // typedef itk::CenteredRigid2DTransform< double > TransformType;
  typedef itk::Rigid2DTransform< double > TransformType;

  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;

  typedef itk::MeanSquaresImageToImageMetric<
                                    FixedImageType,
                                    MovingImageType >
      MeanSquaresMetricType;

  typedef itk::MattesMutualInformationImageToImageMetric<
                                    FixedImageType,
                                    MovingImageType >
      MattesMutualInformationMetricType;

  typedef itk::NormalizedCorrelationImageToImageMetric<
                                    FixedImageType,
                                    MovingImageType >
      NormalizedCorrelationMetricType;

  typedef NormalizedCorrelationMetricType MetricType;

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

  // The Normalized Correlation metric does not use random sampling.
  // Correlation is high at good overlap.
  metric->SetSubtractMean( true ); // We should remove the mean signal level.
  // Not clear ? FIX FIX FIX  optimizer->SetMaximize( true );

  //  The transform object is constructed below and passed to the registration
  //  method.

  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();

  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName(  fixedImageFile->c_str() );
  try
    {
    fixedImageReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Failed to read from fixed image file : " <<
         (*fixedImageFile) << std::endl;
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return -1;
    }

  movingImageReader->SetFileName( movingImageFile->c_str() );
  try
    {
    movingImageReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Failed to read from fixed image file : " <<
         (*movingImageFile) << std::endl;
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return -1;
    }

  registration->SetFixedImage(    fixedImageReader->GetOutput()    );
  registration->SetMovingImage(   movingImageReader->GetOutput()   );

  registration->SetFixedImageRegion(
     fixedImageReader->GetOutput()->GetBufferedRegion() );

  /* we initialize the images based on their centers */
  FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

  const FixedImageType::SpacingType&
    fixedSpacing = fixedImage->GetSpacing();
  const FixedImageType::PointType&
    fixedOrigin  = fixedImage->GetOrigin();

  FixedImageType::SizeType fixedSize  =
          fixedImage->GetLargestPossibleRegion().GetSize();

  TransformType::InputPointType centerFixed;
  centerFixed[0] = fixedOrigin[0] + fixedSpacing[0] * fixedSize[0] / 2.0;
  centerFixed[1] = fixedOrigin[1] + fixedSpacing[1] * fixedSize[1] / 2.0;

  std::cout << "center fixed is :" << centerFixed << std::endl;
  std::cout << "fixed origin ("<< fixedOrigin[0]<<", "<<fixedOrigin[1]<< 
                              ") " << std::endl;
  std::cout << "fixed spacing ("<< fixedSpacing[0]<<", "<<fixedSpacing[1]<< 
                              ") " << std::endl;
  std::cout << "fixed size ("<< fixedSize[0]<<", "<<fixedSize[1]<< 
                              ") " << std::endl;

  MovingImageType::Pointer movingImage = movingImageReader->GetOutput();

  const MovingImageType::SpacingType&
    movingSpacing = movingImage->GetSpacing();
  const MovingImageType::PointType&
    movingOrigin  = movingImage->GetOrigin();

  MovingImageType::SizeType movingSize  =
          movingImage->GetLargestPossibleRegion().GetSize();

  TransformType::InputPointType centerMoving;
  centerMoving[0] = movingOrigin[0] + movingSpacing[0] * movingSize[0] / 2.0;
  centerMoving[1] = movingOrigin[1] + movingSpacing[1] * movingSize[1] / 2.0;

  std::cout << "center moving is :" << centerMoving << std::endl;
  std::cout << "moving origin ("<< movingOrigin[0]<<", "<<movingOrigin[1]<< 
                              ") " << std::endl;
  std::cout << "moving spacing ("<< movingSpacing[0]<<", "<<movingSpacing[1]<< 
                              ") " << std::endl;
  std::cout << "moving size ("<< movingSize[0]<<", "<<movingSize[1]<< 
                              ") " << std::endl;

  transform->SetCenter( centerFixed );
  transform->SetTranslation( centerMoving - centerFixed );
  transform->SetAngle( 0.0 );

  /* For debugging uncomment this: */
     std::cout << "Initial transform from initializer is " 
               << transform << std::endl;
  /* */

  registration->SetInitialTransformParameters( transform->GetParameters() );

  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
  const double translationScale = 1.0 / 1000.0;

  optimizerScales[0] = 1.0;
  optimizerScales[1] = translationScale;
  optimizerScales[2] = translationScale;
  // NOT FOR RIGID2D, ONLY FOR CENTERED optimizerScales[3] = translationScale;
  // ONLY FOR CENTERED optimizerScales[4] = translationScale;

  optimizer->SetScales( optimizerScales );

  optimizer->SetMaximumStepLength( maxOptimizerStepLength );
  optimizer->SetMinimumStepLength( minOptimizerStepLength );
  optimizer->SetNumberOfIterations( numberIterationsOptimizer );

  // Create the Command observer and register it with the optimizer.
  //
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

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
                    registration->GetLastTransformParameters();

  const double finalAngle 		= finalParameters[0];
  const double finalTranslationX    	= finalParameters[1];
  const double finalTranslationY    	= finalParameters[2];

  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();

  const double bestValue = optimizer->GetValue();

  const double finalAngleInDegrees = finalAngle * 45.0 / atan(1.0);

  std::cout << "Result = " << std::endl;
  std::cout << " Angle (radians)   = " << finalAngle  << std::endl;
  std::cout << " Angle (degrees)   = " << finalAngleInDegrees  << std::endl;
  std::cout << " Translation X = " << finalTranslationX  << std::endl;
  std::cout << " Translation Y = " << finalTranslationY  << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Best metric value  = " << bestValue          << std::endl;

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
  resampler->SetInput( movingImageReader->GetOutput() );

  /* Since we may be working with oriented images, ensure that the output
   * direction is set properly.
   */
  resampler->SetOutputDirection( fixedImage->GetDirection() );
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
