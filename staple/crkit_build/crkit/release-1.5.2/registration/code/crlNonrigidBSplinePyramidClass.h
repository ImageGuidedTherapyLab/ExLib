#ifndef CRLNONRIGIDBSPLINEPYRAMIDCLASS_H_
#define CRLNONRIGIDBSPLINEPYRAMIDCLASS_H_

#include "itkImageSource.h"
#include "itkMacro.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageFileReader.h"	
#include "itkImageFileWriter.h"
#include "itkOrientedImage.h"
#include "itkBSplineDeformableTransform.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkNormalizedMutualInformationHistogramImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkProcessObject.h"
#include "itkScaleSkewVersor3DTransform.h"
#include "itkIndent.h"
#include "crlCommandIterationUpdate.h"

namespace crl 
{


template < class TImage, unsigned int TDimension=3 >
class NonrigidBSplinePyramid : public itk::ImageSource<TImage>
{
public:
  typedef NonrigidBSplinePyramid 	Self;
  typedef itk::ImageSource<TImage> 	Superclass;
  typedef itk::SmartPointer<Self> 	Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;
  
  itkNewMacro(Self);
  itkTypeMacro( NonrigidBSplinePyramid, ImageSource );
  
  typedef TImage  FixedImageType;
  typedef TImage  MovingImageType;
  typedef TImage  OutputImageType;
  static const unsigned int ImageDimension = TDimension;
  
  typedef itk::ImageMaskSpatialObject< ImageDimension > 	ImageMaskObjectType;
  typedef itk::Vector< float, ImageDimension >  		VectorType;
  
  typedef itk::OrientedImage< VectorType, ImageDimension >  
                                                 		DeformationFieldType;

  typedef typename itk::ScaleSkewVersor3DTransform<double > BulkTransformType;
  
  static const unsigned int SpaceDimension = ImageDimension;
  static const unsigned int SplineOrder = 3;
  typedef double CoordinateRepType;

  typedef itk::BSplineDeformableTransform<
                              CoordinateRepType,
                              SpaceDimension,
                              SplineOrder >     TransformType;
  
  typedef typename TransformType::ParametersType     ParametersType;


  typedef itk::MattesMutualInformationImageToImageMetric<
                                      FixedImageType,
                                      MovingImageType >    MetricType;
  
  typedef itk::NormalizedMutualInformationHistogramImageToImageMetric<
  				     FixedImageType,
  				     MovingImageType > NMMetricType;
  
  typedef itk::ImageToImageMetric< FixedImageType, MovingImageType > MetricBaseType;


  typedef itk:: LinearInterpolateImageFunction<
                                      MovingImageType,
                                      double          >    InterpolatorType;

  itkSetMacro(NormalizedMI, bool);
  itkGetMacro(NormalizedMI, bool);
  itkBooleanMacro(NormalizedMI);  
  itkSetMacro(ControlPoints, unsigned int);
  itkGetMacro(ControlPoints, unsigned int);
  itkSetMacro(RelaxationFactor, float);
  itkGetMacro(RelaxationFactor, float);
  itkSetMacro(MaxStepLength, float);
  itkGetMacro(MaxStepLength, float);
  itkSetMacro(MinStepLength, float);
  itkGetMacro(MinStepLength, float);
  itkSetMacro(NumberPyramidLevels, unsigned int );
  itkGetMacro(NumberPyramidLevels, unsigned int );
  itkSetMacro(NumberBins, unsigned int );
  itkGetMacro(NumberBins, unsigned int );
  itkSetMacro(NumberSamples, unsigned int );
  itkGetMacro(NumberSamples, unsigned int );
  itkSetMacro(NumberIterations, unsigned int );
  itkGetMacro(NumberIterations, unsigned int );
  itkSetObjectMacro(InitialBulkTransform, BulkTransformType );
  itkGetObjectMacro(InitialBulkTransform, BulkTransformType );
  itkSetObjectMacro(MaskObject,ImageMaskObjectType);
  itkGetConstObjectMacro(MaskObject,ImageMaskObjectType);
  
  itkGetConstObjectMacro(FinalTransform, TransformType);

  void SetFixedImage(const FixedImageType* f) { this->SetNthInput(0,const_cast<FixedImageType*>(f));};
  FixedImageType* GetFixedImage() { return static_cast<FixedImageType*>(this->GetInput(0));};
  void SetMovingImage(const MovingImageType* f) { this->SetNthInput(1,const_cast<MovingImageType*>(f));};
  MovingImageType* GetMovingImage() { return static_cast<MovingImageType*>(this->GetInput(1));};
  
  DeformationFieldType* GetDeformationField() { return reinterpret_cast<DeformationFieldType*>(this->GetOutput(1)); };
  
  void PrintSelf( std::ostream&, itk::Indent  ) const;
  void GenerateData();
  
protected:
  NonrigidBSplinePyramid();
  ~NonrigidBSplinePyramid() {};
  
  
private:
  NonrigidBSplinePyramid(const Self&);		// purposefully not implemented
  
//  void initializeBulkTransformFromFile( const std::string& file, typename TransformType::Pointer splines );

  // parameters
  unsigned int m_ControlPoints;
  float m_MaxStepLength;
  float m_MinStepLength;
  float m_RelaxationFactor;
  unsigned int m_NumberPyramidLevels;
  unsigned int m_NumberBins;
  unsigned int m_NumberSamples;
  unsigned int m_NumberIterations;
  bool m_NormalizedMI;
  typename BulkTransformType::Pointer m_InitialBulkTransform;
  typename TransformType::Pointer m_FinalTransform;
  typename ImageMaskObjectType::Pointer m_MaskObject;
};



} // namespace CRL

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlNonrigidBSplinePyramidClass.txx"
#endif

#endif /*CRLNONRIGIDBSPLINEPYRAMIDCLASS_H_*/

