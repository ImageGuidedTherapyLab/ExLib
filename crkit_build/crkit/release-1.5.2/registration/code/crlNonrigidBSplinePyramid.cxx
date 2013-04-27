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

#include <itkMultiResolutionImageRegistrationMethod.h>
#include <itkMattesMutualInformationImageToImageMetric.h>
#include <itkMeanSquaresImageToImageMetric.h>
#include <itkMattesMutualInformationImageToImageMetric.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>
#include <itkNormalizedMutualInformationHistogramImageToImageMetric.h>
#include <itkMatchCardinalityImageToImageMetric.h>

#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkImage.h>

#include <itkTimeProbesCollectorBase.h>

#include <itkBSplineDeformableTransform.h>
#include <itkRegularStepGradientDescentOptimizer.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkOrientedImage.h>
#include <itkOrientImageFilter.h>

#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSquaredDifferenceImageFilter.h"

// Reading and Writing transforms
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>

#include "itkTransformFactory.h"

#include "itkMultiResolutionPyramidImageFilter.h"

#include "itkBSplineResampleImageFunction.h"
#include "itkIdentityTransform.h"
#include "itkBSplineDecompositionImageFilter.h"

#include "itkAffineTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkScaleSkewVersor3DTransform.h"


#include "configuration.h"

static const    unsigned int    ImageDimension = 3;
typedef  float           PixelType;

typedef itk::OrientedImage< PixelType, ImageDimension >  FixedImageType;
typedef itk::OrientedImage< PixelType, ImageDimension >  MovingImageType;
typedef itk::OrientedImage< PixelType, ImageDimension > OutputImageType;
typedef itk::ImageFileWriter< OutputImageType >  WriterType;

typedef itk::Vector< float, ImageDimension >  VectorType;
typedef itk::OrientedImage< VectorType, ImageDimension >  
                                               DeformationFieldType;

static const unsigned int SpaceDimension = ImageDimension;
static const unsigned int SplineOrder = 3;
typedef double CoordinateRepType;

typedef itk::BSplineDeformableTransform<
                            CoordinateRepType,
                            SpaceDimension,
                            SplineOrder >     TransformType;

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

  typedef itk:: LinearInterpolateImageFunction<
                                    MovingImageType,
                                    double          >    LinearInterpolatorType;

  typedef itk:: NearestNeighborInterpolateImageFunction<
                                    MovingImageType,
                                    double          >
      NearestNeighborInterpolatorType;

  typedef LinearInterpolatorType InterpolatorType;
  // typedef NearestNeighborInterpolatorType InterpolatorType;

void initializeBulkTransformFromFile( const char*,  TransformType::Pointer  );

int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile  movingImageFile outputImagefile ";
    std::cerr << "NumberControlPoints MaxStepLength MinStepLength ";
    std::cerr << "OutputFileTransform NbPyramidLevel numberOfBins numberOfMetricSamples ";
    std::cerr << "nbIterations OutputFileDefField [ InitialScaleSkewVersor3DTransformFile ]" << std::endl;
    return 1;
    }

  //number of CPs in each direction
  unsigned int nctrlptsinit = 5;
  if(argc>3) nctrlptsinit = atoi(argv[4]);

  //Number of Pyramidal levels for multiresolution optimization
  signed int nbpyrlevel = 2;
  if(argc>7) nbpyrlevel = atoi(argv[8]);

  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
 itk::OrientImageFilter<FixedImageType,FixedImageType>::Pointer fixedImageOrienter = itk::OrientImageFilter<FixedImageType,FixedImageType>::New();

  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
 itk::OrientImageFilter<MovingImageType,MovingImageType>::Pointer movingImageOrienter = itk::OrientImageFilter<MovingImageType,MovingImageType>::New();

  fixedImageReader->SetFileName(  argv[1] );
  movingImageReader->SetFileName( argv[2] );

  fixedImageOrienter->SetInput( fixedImageReader->GetOutput() );
  fixedImageOrienter->UseImageDirectionOn();
  fixedImageOrienter->SetDesiredCoordinateOrientationToAxial();
  fixedImageOrienter->Update();

  movingImageOrienter->UseImageDirectionOn();
  movingImageOrienter->SetDesiredCoordinateOrientationToAxial();
  movingImageOrienter->SetInput( movingImageReader->GetOutput() );
  movingImageOrienter->Update();

  FixedImageType::ConstPointer fixedImage = fixedImageOrienter->GetOutput();
  MovingImageType::ConstPointer movingImage = movingImageOrienter->GetOutput();

  std::cout<<"Image voxel size: "<<fixedImage->GetSpacing()<<std::endl;

  std::cout<<"Image dimensions: "<<fixedImage->GetLargestPossibleRegion().GetSize()<<std::endl;

  //get image size
  FixedImageType::SizeType size = fixedImage->GetLargestPossibleRegion().GetSize();
  //get voxel size
  FixedImageType::SpacingType spacing = fixedImage->GetSpacing();

  //find the smallest dimension
  float mindim;
  float dimX = spacing[0]*(float)size[0];
  float dimY = spacing[1]*(float)size[1];
  float dimZ = spacing[2]*(float)size[2];

  std::cout<<"dimensions of the image: "<<dimX<<" "<<dimY<<" "<<dimZ<<std::endl;


  if(dimX<=dimY)
    mindim = dimX;
  else
    mindim = dimY;
  if(mindim >dimZ)
    mindim = dimZ;

  std::cout<<"minimum dimension: "<<mindim<<std::endl;

  //compute the number of CPs in each direction
  float distBetweenCPs = mindim/((float)nctrlptsinit-1.);

  std::cout<<"distance between CPs: "<<distBetweenCPs<<std::endl;

  unsigned int nctrlpts[ImageDimension];
  nctrlpts[0] = (unsigned int)((dimX/distBetweenCPs)+0.5)+1;
  nctrlpts[1] = (unsigned int)((dimY/distBetweenCPs)+0.5)+1;
  nctrlpts[2] = (unsigned int)((dimZ/distBetweenCPs)+0.5)+1;

  std::cout<<"nctrlpts in each direction: "<<nctrlpts[0]<<" "<<nctrlpts[1]<<" "<<nctrlpts[2]<<std::endl;


  typedef TransformType::ParametersType     ParametersType;
  ParametersType parametersScale( (unsigned int) ((nctrlpts[0]+3)*(nctrlpts[1]+3)*(nctrlpts[2]+3))* ImageDimension );
  parametersScale.Fill( 0.0 );


//    try
//      {


//	int sc=nbpyrlevel;

  // Beginning of scale loop
  for(int sc=nbpyrlevel;sc>0;sc--)
    {

      unsigned int nctrlptsScale[ImageDimension];
      nctrlptsScale[0] = (unsigned int) ceil(nctrlpts[0]/(pow((float)2,sc-1)));
      nctrlptsScale[1] = (unsigned int) ceil(nctrlpts[1]/(pow((float)2,sc-1)));
      nctrlptsScale[2] = (unsigned int) ceil(nctrlpts[2]/(pow((float)2,sc-1)));


      std::cout<<std::endl<<std::endl<<std::endl<<std::endl;
      std::cout<<"REGISTRATION AT SCALE "<<sc<<" WITH "<<nctrlptsScale[0]<<" "<<nctrlptsScale[1]<<" "<<nctrlptsScale[2]<<" CONTROL POINTS"<<std::endl;
      std::cout<<std::endl<<std::endl<<std::endl<<std::endl;

  //Max step length for the optimizer
  double tol = 4;
  if(argc>4) tol = (double) atof(argv[5]);

  //Min step length for the optimizer
  double acc = 0.0001;
  if(argc>5) acc = (double) atof(argv[6]);

  //Number of bins for computing the histogram
  unsigned int numberOfBins = 32;
  if(argc>8) numberOfBins = atoi(argv[9]);

  //Number of samples used to compute the metric
  unsigned int numberOfMetricSamples = 100000;
  if(argc>9) numberOfMetricSamples = atoi(argv[10]);

  //Number of iterations (optimizer)
  unsigned int nbIterations = 200;
  if(argc>10) nbIterations = atoi(argv[11]);


  //std::cout<<"Number of control points: "<<nctrlptsinit<<std::endl;

  typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;

  typedef itk::MultiResolutionImageRegistrationMethod<
                                    FixedImageType,
                                    MovingImageType >    RegistrationType;

  // Pyramide
  typedef itk::MultiResolutionPyramidImageFilter<
                                    FixedImageType,
                                    FixedImageType >   FixedImagePyramidType;

  typedef itk::MultiResolutionPyramidImageFilter<
                                    MovingImageType,
                                    MovingImageType >   MovingImagePyramidType;


  FixedImagePyramidType::Pointer fixedImagePyramid = FixedImagePyramidType::New();
  MovingImagePyramidType::Pointer movingImagePyramid = MovingImagePyramidType::New();

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


  registration->SetMetric(        metricBase        );
  registration->SetOptimizer(     optimizer     );
  registration->SetInterpolator(  interpolator  );

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
    metricMattesMI->SetNumberOfHistogramBins( numberOfBins);
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
   NormalizedMutualInformationMetricType::HistogramType::SizeType histogramSize;     histogramSize[0] = numberOfBins;  // 2 elements, since two scalar images
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

  TransformType::Pointer  transform = TransformType::New();

  if ( argc > 13 )      // the only correct argc check in the program
    {
    initializeBulkTransformFromFile( argv[13], transform);
    transform->Print(std::cout);
    }

  registration->SetTransform( transform );
  registration->SetFixedImagePyramid( fixedImagePyramid );
  registration->SetMovingImagePyramid( movingImagePyramid );


  fixedImagePyramid->SetInput( fixedImage );
  movingImagePyramid->SetInput( movingImage );


  fixedImagePyramid->SetNumberOfLevels( 1 );
  movingImagePyramid->SetNumberOfLevels( 1 );


  FixedImagePyramidType::ScheduleType scheduleFixed=fixedImagePyramid->GetSchedule();

  scheduleFixed[0][0]=(int) pow(2.,sc-1);
  scheduleFixed[0][1]=(int) pow(2.,sc-1);
  scheduleFixed[0][2]=(int) pow(2.,sc-1);

  MovingImagePyramidType::ScheduleType scheduleMoving=movingImagePyramid->GetSchedule();

  scheduleMoving[0][0]=(int) pow(2.,sc-1);
  scheduleMoving[0][1]=(int) pow(2.,sc-1);
  scheduleMoving[0][2]=(int) pow(2.,sc-1);

  fixedImagePyramid->SetSchedule(scheduleFixed);
  movingImagePyramid->SetSchedule(scheduleMoving);

  std::cout << std::endl << "Schedule fix: " << fixedImagePyramid->GetSchedule()<< std::endl;;
  std::cout << std::endl << "Schedule mov: " << movingImagePyramid->GetSchedule()<< std::endl;;


  fixedImagePyramid->GetOutput( 0 )->Update();
  movingImagePyramid->GetOutput( 0 )->Update();


  FixedImageType::ConstPointer pyrfixedImage = fixedImagePyramid->GetOutput( 2 );


  registration->SetFixedImage(  fixedImage   );
  registration->SetMovingImage(   movingImage   );

  FixedImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();

  registration->SetFixedImageRegion( fixedRegion );

  typedef TransformType::RegionType RegionType;
  RegionType bsplineRegion;
  RegionType::SizeType   gridSizeOnImage;
  RegionType::SizeType   gridBorderSize;
  RegionType::SizeType   totalGridSize;

  //  gridSizeOnImage.Fill(nctrlptsScale );
  gridSizeOnImage[0]=nctrlptsScale[0];
  gridSizeOnImage[1]=nctrlptsScale[1];
  gridSizeOnImage[2]=nctrlptsScale[2];
  std::cout<<"gridSizeOnImage: "<<gridSizeOnImage<<std::endl;

  gridBorderSize.Fill( 3 );    // Border for spline order = 3 ( 1 lower, 2 upper )
  totalGridSize = gridSizeOnImage + gridBorderSize;

  bsplineRegion.SetSize( totalGridSize );

  typedef TransformType::SpacingType SpacingType;
  SpacingType spacing = fixedImage->GetSpacing();

  typedef TransformType::OriginType OriginType;
  OriginType origin = fixedImage->GetOrigin();;

  FixedImageType::SizeType fixedImageSize = fixedRegion.GetSize();

  for(unsigned int r=0; r<ImageDimension; r++)
    {
    spacing[r] *= floor( static_cast<double>(fixedImageSize[r] - 1)  /
                  static_cast<double>(gridSizeOnImage[r] - 1) );
    origin[r]  -=  spacing[r];
    }

  transform->SetGridSpacing( spacing );
  transform->SetGridOrigin( origin );
  transform->SetGridRegion( bsplineRegion );


  const unsigned int numberOfParameters =
               transform->GetNumberOfParameters();

  ParametersType parameters( numberOfParameters );

  if(sc == nbpyrlevel)
    {
      parameters.Fill( 0.0 );
    }
  else
    {

      for(unsigned int i=0;i<numberOfParameters;i++)
	{
	  parameters[i]=parametersScale[i];
	}
    }


  transform->SetParameters( parameters );
  registration->SetInitialTransformParameters( transform->GetParameters() );

  std::cout << "Initial Parameters at scale "<< sc <<" = " << std::endl;
  std::cout << transform->GetParameters() << std::endl;


  std::cout << "Initial optimizer parameters = " << tol << "   " << acc << "   "<< "1" << std::endl;

    optimizer->SetMaximumStepLength( tol );
    optimizer->SetMinimumStepLength( acc );
    optimizer->SetNumberOfIterations( nbIterations );


  // Add a time probe
  itk::TimeProbesCollectorBase collector;

  registration->SetNumberOfLevels( 1 );



  std::cout << std::endl << "Starting Registration" << std::endl;

  try
    {
    collector.Start( "Registration" );
    registration->StartRegistration();
    collector.Stop( "Registration" );
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return -1;
    }

  OptimizerType::ParametersType finalParameters =
                    registration->GetLastTransformParameters();

  std::cout << "Last Transform Parameters" << std::endl;
  std::cout << finalParameters << std::endl;

  // Report the time taken by the registration
  collector.Report();

  transform->SetParameters( finalParameters );


  // Writing tranfo
  itk::TransformFileWriter::Pointer transwriter;
  transwriter = itk::TransformFileWriter::New();
//  transwriter->SetInput( transform );
  transwriter->AddTransform( transform );
  if (argc > 6) {
    transwriter->SetFileName( argv[7] );
  } else {
    transwriter->SetFileName( "DefaultBSplineTransformOutput.trsf" );
  }
  transwriter->Update();

  try
   {
   transwriter->Update();
   }
  catch( itk::ExceptionObject & excp )
   {
   std::cerr << "Error while saving the transforms" << std::endl;
   std::cerr << excp << std::endl;
  return 0;
  }


  typedef itk::ResampleImageFilter< MovingImageType, FixedImageType >    ResampleFilterType;
  ResampleFilterType::Pointer resample = ResampleFilterType::New();
  resample->SetTransform( transform );
  resample->SetInput( movingImageReader->GetOutput() );
  resample->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputOrigin(  fixedImage->GetOrigin() );
  resample->SetOutputSpacing( fixedImage->GetSpacing() );
  resample->SetOutputDirection( fixedImage->GetDirection() );
  resample->SetDefaultPixelValue( 0 );

  WriterType::Pointer      writer =  WriterType::New();
  writer->SetFileName( argv[3] );
  writer->UseCompressionOn( );

  writer->SetInput( resample->GetOutput()   );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return -1;
    }


  // Generate the explicit deformation field resulting from
  // the registration.


  DeformationFieldType::Pointer field = DeformationFieldType::New();
  field->SetRegions( fixedRegion );
  field->SetOrigin( fixedImage->GetOrigin() );
  field->SetSpacing( fixedImage->GetSpacing() );
  field->SetDirection( fixedImage->GetDirection() );
  field->Allocate();

  typedef itk::ImageRegionIterator< DeformationFieldType > FieldIterator;
  FieldIterator fi( field, fixedRegion );

  fi.GoToBegin();

  TransformType::InputPointType  fixedPoint;
  TransformType::OutputPointType movingPoint;
  DeformationFieldType::IndexType index;

  VectorType displacement;

  while( ! fi.IsAtEnd() )
    {
    index = fi.GetIndex();
    field->TransformIndexToPhysicalPoint( index, fixedPoint );
    movingPoint = transform->TransformPoint( fixedPoint );

    displacement = movingPoint - fixedPoint;
    fi.Set( displacement );
    ++fi;
    }

  typedef itk::ImageFileWriter< DeformationFieldType >  FieldWriterType;
  FieldWriterType::Pointer fieldWriter = FieldWriterType::New();

  fieldWriter->UseCompressionOn();
  fieldWriter->SetInput( field );

  if( argc >= 12 )
    {
    fieldWriter->SetFileName( argv[12] );
    try
      {
      fieldWriter->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Exception thrown " << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
      }
    }


  // Next scale preparation
  TransformType::Pointer  transformHigh = TransformType::New();

  for(unsigned int i=0; i<ImageDimension; i++)
    {
      if( (nctrlptsScale[i]<nctrlpts[i]) && (sc>1) ) nctrlptsScale[i] = (unsigned int) ceil(nctrlpts[i]/(pow((float)2,sc-2)));
    }

  RegionType::SizeType   gridHighSizeOnImage;

  gridHighSizeOnImage[0]=nctrlptsScale[0];
  gridHighSizeOnImage[1]=nctrlptsScale[1];
  gridHighSizeOnImage[2]=nctrlptsScale[2];
//   std::cout<<"gridSizeOnImage: "<<gridSizeOnImage<<std::endl;

  totalGridSize = gridHighSizeOnImage + gridBorderSize;

  bsplineRegion.SetSize( totalGridSize );

  SpacingType spacingHigh = fixedImage->GetSpacing();
  OriginType  originHigh  = fixedImage->GetOrigin();;

  for(unsigned int rh=0; rh<ImageDimension; rh++)
    {
    spacingHigh[rh] *= floor( static_cast<double>(fixedImageSize[rh] - 1)  /
                            static_cast<double>(gridHighSizeOnImage[rh] - 1) );
    originHigh[rh]  -=  spacingHigh[rh];
    }

  transformHigh->SetGridSpacing( spacingHigh );
  transformHigh->SetGridOrigin( originHigh );
  transformHigh->SetGridRegion( bsplineRegion );

  ParametersType parametersHigh( transformHigh->GetNumberOfParameters() );
  parametersHigh.Fill( 0.0 );

  unsigned int counter = 0;

  for ( unsigned int k = 0; k < SpaceDimension; k++ )
    {
    typedef TransformType::ImageType ParametersImageType;
    typedef itk::ResampleImageFilter<ParametersImageType,ParametersImageType> ResamplerType;
    ResamplerType::Pointer upsampler = ResamplerType::New();

    typedef itk::BSplineResampleImageFunction<ParametersImageType,double> FunctionType;
    FunctionType::Pointer function = FunctionType::New();

    typedef itk::IdentityTransform<double,SpaceDimension> IdentityTransformType;
    IdentityTransformType::Pointer identity = IdentityTransformType::New();

    upsampler->SetInput( transform->GetCoefficientImage()[k] );
    upsampler->SetInterpolator( function );
    upsampler->SetTransform( identity );
    upsampler->SetSize( transformHigh->GetGridRegion().GetSize() );
    upsampler->SetOutputSpacing( transformHigh->GetGridSpacing() );
    upsampler->SetOutputOrigin( transformHigh->GetGridOrigin() );

    typedef itk::BSplineDecompositionImageFilter<ParametersImageType,ParametersImageType>
      DecompositionType;
    DecompositionType::Pointer decomposition = DecompositionType::New();

    decomposition->SetSplineOrder( SplineOrder );
    decomposition->SetInput( upsampler->GetOutput() );
    decomposition->Update();

    ParametersImageType::Pointer newCoefficients = decomposition->GetOutput();

    // copy the coefficients into the parameter array
    typedef itk::ImageRegionIterator<ParametersImageType> Iterator;
    Iterator it( newCoefficients, transformHigh->GetGridRegion() );
    while ( !it.IsAtEnd() )
      {

	std::cout<<counter<<" ";

	parametersScale[ counter ] = it.Get();
	parametersHigh[ counter ] = it.Get();
	counter++;
	++it;
      }

    }


  transformHigh->SetParameters( parametersHigh );

  std::cout << "Estimated New Parameters at scale "<< sc <<" = " << std::endl;
    std::cout << transformHigh->GetParameters() << std::endl;


} //end of scales loop


//      }
//    catch( itk::ExceptionObject & excp )
//      {
//      std::cerr << "Exception thrown " << std::endl;
//      std::cerr << excp << std::endl;
//      return EXIT_FAILURE;
//      }
  std::cout << std::endl;
  exit(0);
}


void initializeBulkTransformFromFile( const char* file, TransformType::Pointer splines )
{
    std::cout << "reading bulk transform from file " << file << std::endl;

    itk::TransformFileReader::Pointer trsfReader;
    trsfReader = itk::TransformFileReader::New();
    trsfReader->SetFileName( file );
    try {
      trsfReader->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
      std::cerr << "Error while reading transform file" << file << std::endl;
      std::cerr << excp << std::endl;
      exit(1);
    }

    typedef itk::TransformFileReader::TransformListType * TransformListType;
    TransformListType transforms = trsfReader->GetTransformList();
    std::cout << "Number of transforms = " << transforms->size() << std::endl;

    typedef itk::AffineTransform< double, 3> AffineTransformType;
    typedef itk::ScaleSkewVersor3DTransform< double > ScaleSkewVersor3DTransformType;
    ScaleSkewVersor3DTransformType::Pointer trsf_read =
                 ScaleSkewVersor3DTransformType::New();

    if (transforms->size() == 1) {
      itk::TransformFileReader::TransformListType::const_iterator trsfit =
           transforms->begin();

      if (!strcmp((*trsfit)->GetNameOfClass(),"ScaleSkewVersor3DTransform")) {
        trsf_read = static_cast<ScaleSkewVersor3DTransformType*>((*trsfit).GetPointer());
        std::cout << "Just loaded a scale skew versor 3d transform." << std::endl;
        trsf_read->Print(std::cout);
      } else {
        std::cerr << "Can't initialize from " << (*trsfit)->GetNameOfClass() << std::endl;
        exit(1);
      }

    } else {
      std::cerr << "Can't initialize from " << transforms->size() <<
             " transforms." << std::endl;
      exit(1);
    }
    std::cout << "Initializing using the transform read in. " << std::endl;

    // Will this be the right number of parameters in the correct order ?
    // The critical representation of every transform involves:
    // 1. the center
    // 2. the translation
    // 3. the matrix
//    transform->SetCenter( trsf_read->GetCenter() );
//    transform->SetTranslation( trsf_read->GetTranslation() );
//    transform->SetRotation( trsf_read->GetVersor() );
    splines->SetBulkTransform( trsf_read );

}
