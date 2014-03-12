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


/* Missing:
 *   add a capacity to handle different types of transforms:
 *      translation only
 *      versor 3D     (        versor + translation)
 *      similarity 3D (scale + versor + translation)
 *
 * MERGE WITH THE AffineMIRegistration code.
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
#include <itkRecursiveMultiResolutionPyramidImageFilter.h>
#include <itkMultiResolutionPyramidImageFilter.h>
#include <itkMeanSquaresImageToImageMetric.h>
#include <itkMattesMutualInformationImageToImageMetric.h>
#include <itkNormalizedMutualInformationHistogramImageToImageMetric.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>

#include "crlLinearInterpolateImageFunction.h"
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkImage.h>
#include <itkOrientedImage.h>
#include <itkOrientImageFilter.h>


#include <itkVersor.h>
#include <itkVersorRigid3DTransform.h>
#include <itkCenteredTransformInitializer.h>
#include <itkCenteredVersorTransformInitializer.h>
#include <vnl/vnl_inverse.h>
#include <vnl/algo/vnl_determinant.h>
#include <itkImageMomentsCalculator.h>

/* In this code, we estimate an optimal VersorRigid3DTransform.
 *   We use a specialized optimizer designed for this type of transform.
 */
#include <itkVersorRigid3DTransformOptimizer.h>
#include <itkPowellOptimizer.h>

#include <itkResampleImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <tclap/CmdLine.h>
#include "configuration.h"

//  The following section of code implements a Command observer
//  that will monitor the evolution of the registration process.
//
#include "itkCommand.h"

class RegistrationIterationUpdate : public itk::Command
{
public:
  typedef  RegistrationIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:
  RegistrationIterationUpdate() { writePyramidImagesToDisk = false; };
public:

  // Turn this to true when debugging - it will save out images.
  // Turn it to false when using it in production.
  bool writePyramidImagesToDisk;
  typedef itk::OrientedImage< float, 3>  FixedImageType;
  typedef itk::OrientedImage< float, 3>  MovingImageType;
  typedef itk::ImageFileWriter< FixedImageType >  FixedWriterType;
  typedef itk::ImageFileWriter< MovingImageType >  MovingWriterType;


  typedef itk::VersorRigid3DTransformOptimizer     VersorOptimizerType;
  typedef itk::PowellOptimizer     PowellOptimizerType;
  // typedef VersorOptimizerType OptimizerType;
  typedef PowellOptimizerType OptimizerType;

  typedef   OptimizerType   *    OptimizerPointer;
  typedef itk::MultiResolutionImageRegistrationMethod<
                                    FixedImageType,
                                    MovingImageType >    RegistrationType;
  typedef RegistrationType *RegistrationPointer;

  void SetWritePyramidToDisk(bool write) {
    writePyramidImagesToDisk = write;
  }

  void Execute(const itk::Object *caller, const itk::EventObject & event)
    {
      return;
    }

  void Execute(itk::Object * object, const itk::EventObject & event)
    {
      if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
        return;
        }
      RegistrationPointer registration = dynamic_cast<RegistrationPointer>(
                            object );
      std::cout << "Multiresolution Level : " <<
           registration->GetCurrentLevel()  << std::endl;

      std::cout << "Fixed image at this level is:" << std::endl;
    FixedImageType::Pointer image =
           registration->GetFixedImagePyramid()->GetOutput(
                 registration->GetCurrentLevel() );

    FixedImageType::RegionType region = image->GetLargestPossibleRegion();
    FixedImageType::PointType origin = image->GetOrigin();
    FixedImageType::SpacingType spacing = image->GetSpacing();
    FixedImageType::DirectionType dir = image->GetDirection();
    FixedImageType::IndexType index = region.GetIndex();
    FixedImageType::SizeType size = region.GetSize();

 std::cout <<
      "Size: " << size << std::endl <<
      "Index: " << index << std::endl <<
      "Spacing: " << spacing << std::endl <<
      "Origin: " << origin << std::endl <<
      "Direction: " << std::endl << dir << std::endl;

      std::cout << std::endl;

      std::cout << "Moving image at this level is:" << std::endl;
    MovingImageType::Pointer movingImage =
           registration->GetMovingImagePyramid()->GetOutput(
                 registration->GetCurrentLevel() );

    MovingImageType::RegionType mregion = movingImage->GetLargestPossibleRegion();
    MovingImageType::PointType morigin = movingImage->GetOrigin();
    MovingImageType::SpacingType mspacing = movingImage->GetSpacing();
    MovingImageType::DirectionType mdir = movingImage->GetDirection();
    MovingImageType::IndexType mindex = mregion.GetIndex();
    MovingImageType::SizeType msize = mregion.GetSize();

    std::cout <<
      "Size: " << msize << std::endl <<
      "Index: " << mindex << std::endl <<
      "Spacing: " << mspacing << std::endl <<
      "Origin: " << morigin << std::endl <<
      "Direction: " << std::endl << mdir << std::endl;

      std::cout << std::endl;


      if (writePyramidImagesToDisk) {
        FixedWriterType::Pointer fwriter = FixedWriterType::New();
        std::stringstream ss;
        std::string fname("fixedimage-pyramid-level-");
        ss << "fixedimage-pyramid-level-" << registration->GetCurrentLevel();
        ss >> fname;
        fname += ".nrrd";
        std::cout << "Writing fixed image " << fname << std::endl;
        fwriter->SetFileName(fname);
        fwriter->SetInput(image);
        fwriter->Update();

        MovingWriterType::Pointer mwriter = MovingWriterType::New();
        ss.str(std::string()); ss.clear();
        fname = "";
        ss << "movingimage-pyramid-level-" << registration->GetCurrentLevel();
        ss >> fname;
        fname += ".nrrd";
        std::cout << "Writing moving image " << fname << std::endl;
        mwriter->SetFileName(fname);
        mwriter->SetInput( movingImage );
        mwriter->Update();
      }

      OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >(
        registration->GetOptimizer() );

      std::cout << "Optimizer is at: " << std::endl;
      std::cout << "Value " << optimizer->GetValue() << "   ";
      std::cout << "Parameters" << optimizer->GetCurrentPosition() << std::endl;
      std::cout << std::endl;
    }
};

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
  typedef itk::VersorRigid3DTransformOptimizer     VersorOptimizerType;
  typedef itk::PowellOptimizer     PowellOptimizerType;
  typedef PowellOptimizerType OptimizerType;
  typedef   const OptimizerType   *    OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
        return;
        }
      OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( object );

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
  unsigned int numberIterationsOptimizer = 200; // (good default 200)
  unsigned int transformInitializationScheme = 0; // 0 - moments
                                               // 1 - geometry
                                               // 2 - identity
                                               // 3 - from file
  float maxOptimizerStepLength = 0.2; // (good default between 0.2 - 4.0)
  float minOptimizerStepLength = 0.0001; // (good default 0.0001)
  float PowellStepLength = 0.01; // (good default is 0.01)
  float PowellStepTolerance = 0.001; // (good default 0.001)
  bool writePyramidToDisk = false;

  std::string *fixedImageFile = new std::string("");
  std::string *movingImageFile = new std::string("");
  std::string *outputImageFile = new std::string("");
  std::string *outputTransformFile = new std::string("");
  std::string *initialTransformFile = new std::string("");

  std::string *metricName = new std::string("");

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',
                         CRKIT_VERSION_STRING );
      TCLAP::UnlabeledValueArg<std::string> fixedImageArg("fixedImage","Fixed Image File Name",true,"","fixed image file name",cmd);
      TCLAP::UnlabeledValueArg<std::string> movingImageArg("movingImage","Moving Image File Name",true,"","moving image file name",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputImageArg("outputImage","Output File Name",true,"","output image file name",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputTransformArg("outputTransform","Output File for Transform",true,"","transform file name",cmd);
      TCLAP::ValueArg<std::string> metricNameArg("","metricName","Name of metric",false,"","name of metric (string: one of mi, normmi, normcorr)",cmd);
      TCLAP::ValueArg<unsigned int> binsArg("b","numberOfBins","Number of Bins",false,numberOfBins,"bins (uint)",cmd);
      TCLAP::ValueArg<unsigned int> samplesArg("s","numberOfMetricSamples","Number of Metric Samples",false,numberOfMetricSamples,"samples (uint)",cmd);
      TCLAP::ValueArg<unsigned int> pyramidArg("p","pyramidLevels","Number of Multiresolution Levels",false,numberOfMultiresolutionLevels,"pyramid levels (uint)",cmd);
      TCLAP::ValueArg<unsigned int> iterationsArg("i","numberOptimizerIterations","Number of Optimizer Iterations",false,numberIterationsOptimizer,"iterations (uint)",cmd);
      TCLAP::ValueArg<bool> writePyramidToDiskArg("","writePyramidToDisk","Write the Image Pyramid to Disk",false,writePyramidToDisk,"write pyramid to disk (bool)",cmd);
      TCLAP::ValueArg<float> maxStepArg("m","maxStepOptimizer","Maximum Optimizer Step Length",false,maxOptimizerStepLength,"max step (float)",cmd);
      TCLAP::ValueArg<float> minStepArg("n","minStepOptimizer","Minimum Optimizer Step Length",false,minOptimizerStepLength,"min step (float)",cmd);

      // Either specify an initial rigid body transform, or provide a
      // suggestion as to how to do the initialization of the transform.
      TCLAP::ValueArg<std::string> initialTransformArg("l", "initialTransform","Initial Transform File",false,"","initial transform file name",cmd);
      TCLAP::ValueArg<unsigned int> transformInitializationArg("t","transformInitializationScheme","Mechanism for initializing the transform",false,transformInitializationScheme,"transform initialization scheme (uint : 0 moments translation only, 1 geometrical center,2 identity transform, 3 moments translation and rotation)",cmd);

      cmd.parse(argc,argv);

      if (binsArg.isSet()) numberOfBins = binsArg.getValue();
      if (samplesArg.isSet()) numberOfMetricSamples = samplesArg.getValue();
      if (pyramidArg.isSet()) numberOfMultiresolutionLevels = pyramidArg.getValue();
      if (writePyramidToDiskArg.isSet()) {
        writePyramidToDisk = writePyramidToDiskArg.getValue();
      }
      if (iterationsArg.isSet()) numberIterationsOptimizer = iterationsArg.getValue();
      if (maxStepArg.isSet()) maxOptimizerStepLength = maxStepArg.getValue();
      if (minStepArg.isSet()) minOptimizerStepLength = minStepArg.getValue();
      /* We allow initialization by specifying an initial transform or one of
       * moments, geometry or identity transform initialization.
       */

      if (transformInitializationArg.isSet()) {
         transformInitializationScheme =
               transformInitializationArg.getValue();
      } else if (initialTransformArg.isSet()) {
        initialTransformFile = new std::string(initialTransformArg.getValue());
      }
      if (fixedImageArg.isSet()) fixedImageFile = new std::string(fixedImageArg.getValue());
      if (movingImageArg.isSet()) movingImageFile = new std::string(movingImageArg.getValue());
      if (outputImageArg.isSet()) outputImageFile = new std::string(outputImageArg.getValue());
      if (outputTransformArg.isSet()) outputTransformFile = new std::string(outputTransformArg.getValue());
      if (metricNameArg.isSet()) {
        metricName = new std::string(metricNameArg.getValue());
      }
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

  typedef itk::VersorRigid3DTransformOptimizer VersorOptimizerType;
  typedef itk::PowellOptimizer PowellOptimizerType;
  //typedef VersorOptimizerType OptimizerType;
  typedef PowellOptimizerType OptimizerType;

  typedef itk::ImageToImageMetric<FixedImageType, MovingImageType >
                 MetricBaseType;

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

  typedef itk::NormalizedMutualInformationHistogramImageToImageMetric<
                                    FixedImageType,
                                    MovingImageType >
      NormalizedMutualInformationMetricType;

  typedef itk:: LinearInterpolateImageFunction<
                                    MovingImageType,
                                    double          >    itkLinearInterpolatorType;

  typedef crl:: LinearInterpolateImageFunction<
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

  typedef itk::RecursiveMultiResolutionPyramidImageFilter<
                                    InternalImageType,
                                    InternalImageType >   FixedImagePyramidType;

  typedef itk::RecursiveMultiResolutionPyramidImageFilter<
                                    InternalImageType,
                                    InternalImageType >  MovingImagePyramidType;

  MetricBaseType::Pointer         metricBase;
  MattesMutualInformationMetricType::Pointer metricMattesMI =
                                  MattesMutualInformationMetricType::New();
  NormalizedCorrelationMetricType::Pointer metricNormCorr =
                                    NormalizedCorrelationMetricType::New();
  NormalizedMutualInformationMetricType::Pointer metricNMI =
                              NormalizedMutualInformationMetricType::New();
  // Default metric is Mattes Mutual Information
  if (metricName->empty()) {
    metricBase = static_cast< NormalizedCorrelationMetricType *>(
                      (metricNormCorr).GetPointer() );
    std::cout << "Metric is: NormalizedCorrelation" << std::endl;
  } else if (metricName->compare("mi") == 0) {
    metricBase = static_cast< MattesMutualInformationMetricType *>(
                      (metricMattesMI).GetPointer() );
    std::cout << "Metric is: MattesMutualInformation" << std::endl;
  } else if (metricName->compare("normmi") == 0) {
    metricBase = static_cast< NormalizedMutualInformationMetricType *>(
                      (metricNMI).GetPointer() );
    std::cout << "Metric is: NormalizedMutualInformation" << std::endl;
  } else if (metricName->compare("normcorr") == 0) {
    metricBase = static_cast< NormalizedCorrelationMetricType *>(
                      (metricNormCorr).GetPointer() );
    std::cout << "Metric is: NormalizedCorrelation" << std::endl;
  } else {
    std::cerr << "Don't know about metric called " << (*metricName) <<std::endl;
    std::cerr << "Try: mi, normmi or normcorr" << std::endl;
    return EXIT_FAILURE;
  }

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
  registration->SetMetric(        metricBase        );
  registration->SetInterpolator(  interpolator  );
  registration->SetFixedImagePyramid(  fixedImagePyramid );
  registration->SetMovingImagePyramid(  movingImagePyramid );

  // Select how we initialize the metric based on what the metric
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
    // Set up key parameters for the metric
    metricMattesMI->SetNumberOfHistogramBins( numberOfBins );
    if (numberOfMetricSamples == 0) {
      // samples every voxel - not so long on modern hardware, and much more
      // robust than the sampling strategy.
      metricMattesMI->UseAllPixelsOn();
    } else {
      metricMattesMI->SetNumberOfSpatialSamples( numberOfMetricSamples );
    }
    // Initialize the random number generator for the sampling process.
    // select a seed based on the clock, to microsecond precision.
    unsigned int seed = 0;
    struct timeval preciseTimeStart;
    gettimeofday(&preciseTimeStart, NULL);
    seed = (unsigned int)(preciseTimeStart.tv_usec & 0xFFFFFFFF);
    // std::cout << "Mattes MI random seed is " << seed << std::endl;
    metricMattesMI->ReinitializeSeed( seed );
  } else if (!strcmp(
    metricBase->GetNameOfClass(),"NormalizedMutualInformationHistogramImageToImageMetric")
     ) {
     // Set up key parameters for the normalized mutual information metric
     NormalizedMutualInformationMetricType::HistogramType::SizeType histogramSize;
     histogramSize[0] = numberOfBins;  // 2 elements, since two scalar images
     histogramSize[1] = numberOfBins;

     metricNMI->SetHistogramSize( histogramSize );
    // The Normalized Mutual Information metric does not use random sampling.

    // Push joint entropy H(A,B) to be small, so the metric gets larger.
    optimizer->SetMaximize( true );

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
  // Read in the image so that the orienter doesn't report wrong origin info.
  try {
    fixedImageReader->Update();
  } catch ( itk::ExceptionObject &err ) {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }

  fixedImageOrienter->UseImageDirectionOn();
  fixedImageOrienter->SetDesiredCoordinateOrientationToAxial();
  fixedImageOrienter->SetInput( fixedImageReader->GetOutput() );
  fixedImageOrienter->Update();

  movingImageReader->SetFileName( movingImageFile->c_str() );
  // Read in the image so that the orienter has data fast.
  try {
    movingImageReader->Update();
  } catch ( itk::ExceptionObject &err ) {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }
  movingImageOrienter->UseImageDirectionOn();
  movingImageOrienter->SetDesiredCoordinateOrientationToAxial();
  movingImageOrienter->SetInput( movingImageReader->GetOutput() );
  movingImageOrienter->Update();

  // Put the data through an orient image filter so that they are appropriately
  // corresponding when we run the registration.

  registration->SetFixedImage(    fixedImageOrienter->GetOutput()    );
  registration->SetMovingImage(   movingImageOrienter->GetOutput()   );

  registration->SetFixedImageRegion(
     fixedImageOrienter->GetOutput()->GetLargestPossibleRegion()
  );
  //     fixedImageOrienter->GetOutput()->GetBufferedRegion() );

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

  // Check if an initial transform file has been supplied, or if we are
  // estimating this on our own.
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
    AffineTransformType::Pointer affine_read = AffineTransformType::New();
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
      } else if (!strcmp((*trsfit)->GetNameOfClass(),"AffineTransform"))
      {
        affine_read = static_cast<AffineTransformType*>((*trsfit).GetPointer());
        std::cout << "Just loaded an affine transform." << std::endl;
        affine_read->Print(std::cout);
       // This is the recommended procedure for initializing when the
        // transforms have different mappings of actual parameter to index
        // in the array of parameters, as is the case here.
        transform->SetCenter( affine_read->GetCenter() );
        transform->SetTranslation( affine_read->GetTranslation() );
        // There is no built in means to infer the rotation implied by the
        // transformation, so we have to compute it.
        // Shoemake in the paper "Matrix Animation and Polar Decomposition".
        // http://citeseer.ist.psu.edu/shoemake92matrix.html
        typedef vnl_matrix<double> VnlMatrixType;
        VnlMatrixType M = affine_read->GetMatrix().GetVnlMatrix();
        VnlMatrixType PQ = M;
        VnlMatrixType NQ = M;
        VnlMatrixType PQNQDiff;

        typedef itk::Versor<double> VersorType;

        const unsigned int maximumIterations = 100;

        for(unsigned int ni = 0; ni < maximumIterations; ni++ )
        {
          // Average current Qi with its inverse transpose
          NQ = ( PQ + vnl_inverse_transpose( PQ ) ) / 2.0;
          PQNQDiff = NQ - PQ;
          if( PQNQDiff.frobenius_norm() < 1e-7 )
          {
            std::cout << "Polar decomposition used "
                 << ni << " iterations " << std::endl;
            break;
          } else {
            PQ = NQ;
          }
        }
        VnlMatrixType QMatrix;
        QMatrix = NQ;
        std::cout << "Initial Matrix = " << std::endl << M << std::endl;
        std::cout << "Q Matrix = " << std::endl << QMatrix << std::endl;
        VersorType versor;
        versor.Set( QMatrix );

        transform->SetRotation( versor );

      } else {
        std::cerr << "Can't initialize from transform of type " <<
                ( (*trsfit)->GetNameOfClass() ) << " ." << std::endl;
      return EXIT_FAILURE;
      }
    } else {
      // Future work could allow composition of all of the transforms.
      std::cerr << "Can't initialize from " << transforms->size() <<
             " transforms." << std::endl;
      return EXIT_FAILURE;
    }
  } else {
      // moments, default
      if ( (transformInitializationScheme == 0) ||
           (transformInitializationScheme == 3) ) {  // moments, default
        // The calculations the initializer does, are redone, correctly below.
        initializer->MomentsOn();
        initializer->SetComputeRotation( false );
        if (transformInitializationScheme == 3) {
          initializer->SetComputeRotation( true ); 
        }
        initializer->InitializeTransform();

      } else if (transformInitializationScheme == 1) { // geometry centers
        /* The printing of the initializer fails unless the moments have
         * been calculated.
         */
        initializer->MomentsOn();
        initializer->InitializeTransform();
        initializer->GeometryOn();
        initializer->InitializeTransform();

      } else if (transformInitializationScheme == 2) { // identity transform
        ; // Here we simply don't call the initializer, so the transform
          // is used with its default initialization.
      }
  }

  /* For debugging uncomment this: */
  /* It turns out, ITK 3.6 and below have a bug such that an initializer
   * in geometry only mode cannot print itself, as it tries to print the
   * moments information. */
     std::cout << "Initial transform from initializer is " << std::endl
               << transform << std::endl;
     std::cout << "The initial transform parameters are " << std::endl
               << transform->GetParameters() << std::endl;
     std::cout << "The fixed transform parameters are " << std::endl
               << transform->GetFixedParameters() << std::endl;
  /* */

  //  We now pass the parameters of the current transform as the initial
  //  parameters to be used when the registration process starts.

  registration->SetInitialTransformParameters( transform->GetParameters() );

  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
  const double translationScale = 1.0 / 1000.0;

  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = translationScale;
  optimizerScales[4] = translationScale;
  optimizerScales[5] = translationScale;

  optimizer->SetScales( optimizerScales );

  // This initialization of optimizer parameters depends on what the
  // underlying type of the optimizer is.
  // For the Versor Optimizer
/*
  optimizer->SetMaximumStepLength( maxOptimizerStepLength );
  optimizer->SetMinimumStepLength( minOptimizerStepLength );
  optimizer->SetNumberOfIterations( numberIterationsOptimizer );
*/

  // For the Powell Optimizer
  // The command line optimizer parameters currently do not set these.
  optimizer->SetMaximumIteration( numberIterationsOptimizer );
  optimizer->SetStepLength( PowellStepLength );
  optimizer->SetStepTolerance( PowellStepTolerance );

  optimizer->SetInitialPosition( transform->GetParameters() );

  // Create the Command observer and register it with the optimizer.
  //
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  RegistrationIterationUpdate::Pointer regCommand =
                              RegistrationIterationUpdate::New();
  registration->AddObserver( itk::IterationEvent(), regCommand );

  regCommand->SetWritePyramidToDisk( writePyramidToDisk );

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

  /* For debugging:

  const double versorX              = finalParameters[0];
  const double versorY              = finalParameters[1];
  const double versorZ              = finalParameters[2];
  const double finalTranslationX    = finalParameters[3];
  const double finalTranslationY    = finalParameters[4];
  const double finalTranslationZ    = finalParameters[5];

  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();

  const double bestValue = optimizer->GetValue();

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
  */

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

  resampler->SetOutputParametersFromImage( fixedImageOrienter->GetOutput() );
  resampler->SetDefaultPixelValue( 0.0 );
  resampler->Update();

  typedef  float OutputPixelType;

  typedef itk::OrientedImage< OutputPixelType, Dimension > OutputImageType;

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

  int returnCode = EXIT_SUCCESS;

  try
    {
      writer->Update();
    }
  catch (itk::ExceptionObject &excp)
    {
      std::cerr << "Error while saveing output image." << std::endl;
      std::cerr << excp << std::endl;
      returnCode = EXIT_FAILURE;
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
    std::cerr << "Error while saving the transform to file "
              << (*outputTransformFile) << std::endl;
    std::cerr << excp << std::endl;
    returnCode = EXIT_FAILURE;
  }

  return returnCode;
}
