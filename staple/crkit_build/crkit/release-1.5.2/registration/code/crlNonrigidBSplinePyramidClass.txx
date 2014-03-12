#include "crlNonrigidBSplinePyramidClass.h"
#include "itkTransformFileReader.h"
#include "itkMacro.h"
#include "itkProcessObject.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkMultiResolutionPyramidImageFilter.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkResampleImageFilter.h"
#include "itkBSplineResampleImageFunction.h"
#include "itkImageRegionIterator.h"
#include "itkResampleImageFilter.h"
#include <algorithm>


namespace crl
{

template < class TImage, unsigned int TDimension > 
NonrigidBSplinePyramid<TImage,TDimension>
::NonrigidBSplinePyramid()
{
  m_ControlPoints = 20;
  m_MaxStepLength = 10.;
  m_MinStepLength = 0.01;
  m_RelaxationFactor = 0.7;
  m_NumberPyramidLevels = 2;
  m_NumberBins = 20;
  m_NumberIterations = 200;
  m_NumberSamples = 0;
  m_InitialBulkTransform = 0;
  m_FinalTransform = 0;
  m_MaskObject = 0;
  
  // filter inputs 
  this->SetNumberOfRequiredInputs(2);	// 0==fixed, 1==moving
  
  // filter outputs ( 0 == resampled image, 1 == deformation field )
  this->itk::ProcessObject::SetNumberOfRequiredOutputs(2);
  typename OutputImageType::Pointer output = OutputImageType::New();
  this->SetNthOutput(0, output.GetPointer() );
  typename DeformationFieldType::Pointer def = DeformationFieldType::New();
  this->SetNthOutput(1, def.GetPointer() );
}

template < class TImage, unsigned int TDimension >
void NonrigidBSplinePyramid<TImage,TDimension>
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
  Superclass::PrintSelf(os, indent);
}


template < class TImage, unsigned int TDimension > 
void NonrigidBSplinePyramid<TImage,TDimension>
::GenerateData()
{
  std::cerr << "we're in GenerateData()!!!!!" << std::endl;

  if ( ImageDimension != 3 )
    {
    itkExceptionMacro("we're not set up for image dimensionalities other than 3 yet");
    }

  this->AllocateOutputs();
  
  typename FixedImageType::Pointer fixedImage = this->GetFixedImage();
  typename MovingImageType::Pointer movingImage = this->GetMovingImage();  

  std::cout<<"Image voxel size: "<<fixedImage->GetSpacing()<<std::endl;
  std::cout<<"Image dimensions: "<<fixedImage->GetLargestPossibleRegion().GetSize()<<std::endl;

  //get image size
  typename FixedImageType::SizeType size = fixedImage->GetLargestPossibleRegion().GetSize();
  typename FixedImageType::SpacingType spacing = fixedImage->GetSpacing();

  //find the smallest dimension
  float dimX = spacing[0]*(float)size[0];
  float dimY = spacing[1]*(float)size[1];
  float dimZ = spacing[2]*(float)size[2];

  std::cout<<"dimensions of the image: "<<dimX<<" "<<dimY<<" "<<dimZ<<std::endl;

  float mindim = std::min(std::min(dimX,dimY),dimZ);

  std::cout<<"minimum dimension: "<<mindim<<std::endl;
  
  unsigned int& nctrlptsinit = m_ControlPoints;
  unsigned int& nbpyrlevel = m_NumberPyramidLevels;

  //compute the number of CPs in each direction
  float distBetweenCPs = mindim/((float)nctrlptsinit-1.);

  std::cout<<"distance between CPs: "<<distBetweenCPs<<std::endl;

  unsigned int nctrlpts[ImageDimension];
  nctrlpts[0] = (unsigned int)((dimX/distBetweenCPs)+0.5)+1;
  nctrlpts[1] = (unsigned int)((dimY/distBetweenCPs)+0.5)+1;
  nctrlpts[2] = (unsigned int)((dimZ/distBetweenCPs)+0.5)+1;

  std::cout<<"nctrlpts in each direction: "<<nctrlpts[0]<<" "<<nctrlpts[1]<<" "<<nctrlpts[2]<<std::endl;

  ParametersType  parametersScale( (unsigned int) ((nctrlpts[0]+3)*(nctrlpts[1]+3)*(nctrlpts[2]+3))* ImageDimension );
  parametersScale.Fill( 0.0 );

  typename TransformType::Pointer  transform;
  
  // Beginning of scale loop
  for(unsigned int sc=nbpyrlevel;sc>0;sc--)
    {

    unsigned int nctrlptsScale[ImageDimension];
    nctrlptsScale[0] = (unsigned int) ceil(nctrlpts[0]/(pow((float)2, (float)sc-1)));
    nctrlptsScale[1] = (unsigned int) ceil(nctrlpts[1]/(pow((float)2, (float)sc-1)));
    nctrlptsScale[2] = (unsigned int) ceil(nctrlpts[2]/(pow((float)2, (float)sc-1)));

    std::cout<<std::endl<<std::endl<<std::endl<<std::endl;
    std::cout<<"REGISTRATION AT SCALE "<<sc<<" WITH "<<nctrlptsScale[0]<<" "
	  <<nctrlptsScale[1]<<" "<<nctrlptsScale[2]<<" CONTROL POINTS"
	  <<std::endl;
    std::cout<<std::endl<<std::endl<<std::endl<<std::endl;

    //Max step length for the optimizer
    float& tol = m_MaxStepLength;
    //Min step length for the optimizer
    float& acc = m_MinStepLength;
    //Number of bins for computing the histogram
    unsigned int& nbBins = m_NumberBins;
    //Number of samples used to compute the metric
    unsigned int& nbSamples = m_NumberSamples;
    //Number of iterations (optimizer)
    unsigned int& nbIterations = m_NumberIterations;

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


    typename FixedImagePyramidType::Pointer fixedImagePyramid = FixedImagePyramidType::New();
    typename MovingImagePyramidType::Pointer movingImagePyramid = MovingImagePyramidType::New();


    typename MetricBaseType::Pointer         metric;
    
    if ( this->m_NormalizedMI ) 
      {
      typename NMMetricType::Pointer  tmetric = NMMetricType::New();
//      if ( nbSamples != 0 ) tmetric->SetNumberOfSpatialSamples(nbSamples); else tmetric->UseAllPixelsOn();
      metric = tmetric;
      }
    else
      {
      typename MetricType::Pointer tmetric = MetricType::New();
      tmetric->SetNumberOfHistogramBins(nbBins);
      if ( nbSamples != 0 ) tmetric->SetNumberOfSpatialSamples(nbSamples); else tmetric->UseAllPixelsOn();
      tmetric->SetUseCachingOfBSplineWeights(1);
      metric = tmetric;
      }
    typename OptimizerType::Pointer      optimizer     = OptimizerType::New();
    typename InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
    typename RegistrationType::Pointer   registration  = RegistrationType::New();

   
    registration->SetMetric(        metric        );
    registration->SetOptimizer(     optimizer     );
    registration->SetInterpolator(  interpolator  );


//    typename TransformType::Pointer  transform = TransformType::New();
    transform = TransformType::New();

    if ( m_InitialBulkTransform.IsNotNull() )      
      {
      transform->SetBulkTransform(m_InitialBulkTransform);
      }
    
    if ( m_MaskObject.IsNotNull() )
      {
      metric->SetFixedImageMask(m_MaskObject);
      }

    registration->SetTransform( transform );
    registration->SetFixedImagePyramid( fixedImagePyramid );
    registration->SetMovingImagePyramid( movingImagePyramid );

    fixedImagePyramid->SetInput( fixedImage );
    movingImagePyramid->SetInput( movingImage );

    fixedImagePyramid->SetNumberOfLevels( 1 );
    movingImagePyramid->SetNumberOfLevels( 1 );

    typename FixedImagePyramidType::ScheduleType scheduleFixed=fixedImagePyramid->GetSchedule();

    scheduleFixed[0][0]=(int) pow(2.,(double)sc-1);
    scheduleFixed[0][1]=(int) pow(2.,(double)sc-1);
    scheduleFixed[0][2]=(int) pow(2.,(double)sc-1);

    typename MovingImagePyramidType::ScheduleType scheduleMoving=movingImagePyramid->GetSchedule();

    scheduleMoving[0][0]=(int) pow(2.,(double)sc-1);
    scheduleMoving[0][1]=(int) pow(2.,(double)sc-1);
    scheduleMoving[0][2]=(int) pow(2.,(double)sc-1);

    fixedImagePyramid->SetSchedule(scheduleFixed);
    movingImagePyramid->SetSchedule(scheduleMoving);

    std::cout << std::endl << "Schedule fix: " << fixedImagePyramid->GetSchedule()<< std::endl;;
    std::cout << std::endl << "Schedule mov: " << movingImagePyramid->GetSchedule()<< std::endl;;

    fixedImagePyramid->GetOutput( 0 )->Update();
    movingImagePyramid->GetOutput( 0 )->Update();

//   FixedImageType::ConstPointer pyrfixedImage = fixedImagePyramid->GetOutput( 2 );

    registration->SetFixedImage(  fixedImage   );
    registration->SetMovingImage(   movingImage   );

    typename FixedImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();

    registration->SetFixedImageRegion( fixedRegion );

    typedef typename TransformType::RegionType RegionType;
    RegionType bsplineRegion;
    typename RegionType::SizeType   gridSizeOnImage;
    typename RegionType::SizeType   gridBorderSize;
    typename RegionType::SizeType   totalGridSize;

    //  gridSizeOnImage.Fill(nctrlptsScale );
    gridSizeOnImage[0]=nctrlptsScale[0];
    gridSizeOnImage[1]=nctrlptsScale[1];
    gridSizeOnImage[2]=nctrlptsScale[2];
    std::cout<<"gridSizeOnImage: "<<gridSizeOnImage<<std::endl;

    gridBorderSize.Fill( 3 );    // Border for spline order = 3 ( 1 lower, 2 upper )
    totalGridSize = gridSizeOnImage + gridBorderSize;

    bsplineRegion.SetSize( totalGridSize );

    typedef typename TransformType::SpacingType SpacingType;
    SpacingType spacing = fixedImage->GetSpacing();

    typedef typename TransformType::OriginType OriginType;
    OriginType origin = fixedImage->GetOrigin();;

    typename FixedImageType::SizeType fixedImageSize = fixedRegion.GetSize();

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

//    std::cout << "Initial Parameters at scale "<< sc <<" = " << std::endl;
//    std::cout << transform->GetParameters() << std::endl;

    std::cout << "Initial optimizer parameters = " << tol << "   " << acc << "   "<< "1" << std::endl;

    optimizer->SetMaximumStepLength( tol );
    optimizer->SetMinimumStepLength( acc );
    optimizer->SetRelaxationFactor( m_RelaxationFactor );
    optimizer->SetNumberOfIterations( nbIterations );

#if 0
    if ( m_NormalizedMI )
      optimizer->MaximizeOn();
#endif

    // Create the Command observer and register it with the optimizer.
    //
//    crl::CommandIterationUpdate<OptimizerType>::Pointer observer = crl::CommandIterationUpdate<OptimizerType>::New();
//    optimizer->AddObserver( itk::IterationEvent(), observer );
    
    crl::CommandIterationUpdate<OptimizerType>::CreateAndRegister(optimizer);

    // Add a time probe
    itk::TimeProbesCollectorBase collector;

    registration->SetNumberOfLevels( 1 );

    std::cout << std::endl << "Starting Registration" << std::endl;


    collector.Start( "Registration" );
    registration->StartRegistration();
    collector.Stop( "Registration" );
 
    OptimizerType::ParametersType finalParameters =
                    registration->GetLastTransformParameters();

//    std::cout << "Last Transform Parameters" << std::endl;
//    std::cout << finalParameters << std::endl;

    // Report the time taken by the registration
    collector.Report();

    transform->SetParameters( finalParameters );
    
    // okay, resample us some output
    typename OutputImageType::Pointer output = this->GetOutput(0);

    typedef itk::ResampleImageFilter< MovingImageType, FixedImageType >    ResampleFilterType;
    typename ResampleFilterType::Pointer resample = ResampleFilterType::New();
    resample->SetTransform( transform );
    resample->SetInput( movingImage );
    resample->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
    resample->SetOutputOrigin(  fixedImage->GetOrigin() );
    resample->SetOutputSpacing( fixedImage->GetSpacing() );
    resample->SetOutputDirection( fixedImage->GetDirection() );
    resample->SetDefaultPixelValue( 0 );
    resample->GraftOutput( output );
    resample->Update();

    // Generate the explicit deformation field resulting from
    // the registration.


    typename DeformationFieldType::Pointer field = this->GetDeformationField();
    
    std::cerr << "==== DEFORMATION FIELD ====" << std::endl;
    field->Print(std::cerr);
    std::cerr << std::endl;

    typedef itk::ImageRegionIterator< DeformationFieldType > FieldIterator;
//    typename FixedImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();
    FieldIterator fi( field, fixedRegion );

    fi.GoToBegin();

    typename TransformType::InputPointType  fixedPoint;
    typename TransformType::OutputPointType movingPoint;
    typename DeformationFieldType::IndexType index;

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

  

#if 0
    m_FinalTransform->SetBulkTransform( transform->GetBulkTransform() );
    m_FinalTransform->SetFixedParameters( transform->GetFixedParameters() );
    m_FinalTransform->SetParameters( transform->GetParameters() );

    std::cerr << "============== TRANSFORM ============" << std::endl;
    transform->PrintSelf(std::cerr);
    std::cerr << "============== FINAL TRANSFORM ============" << std::endl;
    m_FinalTransform->PrintSelf(std::cerr);
    exit(0);
#endif
    

#if 0 
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
#endif
    if ( sc == 1 )
      {
      std::cerr << "final pass -- setting m_FinalTransform" << std::endl;
      m_FinalTransform = TransformType::New();
      m_FinalTransform->SetBulkTransform(transform->GetBulkTransform());
      m_FinalTransform->SetGridRegion(transform->GetGridRegion());
      m_FinalTransform->SetGridSpacing(transform->GetGridSpacing());
      m_FinalTransform->SetGridOrigin(transform->GetGridOrigin());
      m_FinalTransform->SetParametersByValue(transform->GetParameters() );
      }
    else
      {
    // Next scale preparation
    typename TransformType::Pointer  transformHigh = TransformType::New();

    for(unsigned int i=0; i<ImageDimension; i++)
      {
      if( (nctrlptsScale[i]<nctrlpts[i]) && (sc>1) ) nctrlptsScale[i] = (unsigned int) ceil(nctrlpts[i]/(pow((float)2,(float)sc-2)));
      }

    typename RegionType::SizeType   gridHighSizeOnImage;

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
      typedef typename TransformType::ImageType ParametersImageType;
      typedef itk::ResampleImageFilter<ParametersImageType,ParametersImageType> ResamplerType;
      typename ResamplerType::Pointer upsampler = ResamplerType::New();

      typedef itk::BSplineResampleImageFunction<ParametersImageType,double> FunctionType;
      typename FunctionType::Pointer function = FunctionType::New();

      typedef itk::IdentityTransform<double,SpaceDimension> IdentityTransformType;
      typename IdentityTransformType::Pointer identity = IdentityTransformType::New();

      upsampler->SetInput( transform->GetCoefficientImage()[k] );
      upsampler->SetInterpolator( function );
      upsampler->SetTransform( identity );
      upsampler->SetSize( transformHigh->GetGridRegion().GetSize() );
      upsampler->SetOutputSpacing( transformHigh->GetGridSpacing() );
      upsampler->SetOutputOrigin( transformHigh->GetGridOrigin() );

      typedef itk::BSplineDecompositionImageFilter<ParametersImageType,ParametersImageType>
      	DecompositionType;
      typename DecompositionType::Pointer decomposition = DecompositionType::New();

      decomposition->SetSplineOrder( SplineOrder );
      decomposition->SetInput( upsampler->GetOutput() );
      decomposition->Update();

      typename ParametersImageType::Pointer newCoefficients = decomposition->GetOutput();

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

//    std::cout << "Estimated New Parameters at scale "<< sc <<" = " << std::endl;
//    std::cout << transformHigh->GetParameters() << std::endl;
      }

    } //end of scales loop
  

#if 0
    // okay, resample us some output
    typename OutputImageType::Pointer output = this->GetOutput(0);

    typedef itk::ResampleImageFilter< MovingImageType, FixedImageType >    ResampleFilterType;
    typename ResampleFilterType::Pointer resample = ResampleFilterType::New();
    resample->SetTransform( transform );
    resample->SetInput( movingImage );
    resample->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
    resample->SetOutputOrigin(  fixedImage->GetOrigin() );
    resample->SetOutputSpacing( fixedImage->GetSpacing() );
    resample->SetOutputDirection( fixedImage->GetDirection() );
    resample->SetDefaultPixelValue( 0 );
    resample->GraftOutput( output );
    resample->Update();
#endif




  
  
  std::cout << std::endl;
}


} // namespace crl


