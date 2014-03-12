#ifndef __crlMFMTensorResampleImageFilter_h
#define __crlMFMTensorResampleImageFilter_h

#include "itkTransform.h"
#include "itkImageFunction.h"
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkDiffusionTensor3D.h"
#include "itkImageRegionIterator.h"
#include "itkFixedArray.h"
#include "itkSize.h"
#include <itkSimpleFastMutexLock.h>

using namespace itk;

namespace crl
{

/** \class MFMTensorResampleImageFilter
 * \brief Resample a tensor image via a coordinate transform. Supports multi-tensors via multiple inputs definition.
 *
 * Since this filter produces an image which is a different size than
 * its input, it needs to override several of the methods defined
 * in ProcessObject in order to properly manage the pipeline execution model.
 * In particular, this filter overrides
 * ProcessObject::GenerateInputRequestedRegion() and
 * ProcessObject::GenerateOutputInformation().
 *
 * This filter is implemented as a multithreaded filter.  It provides a 
 * ThreadedGenerateData() method for its implementation.
 */
template <typename TPixel, unsigned int VImageDimension >
class ITK_EXPORT MFMTensorResampleImageFilter:
	public itk::ImageToImageFilter< itk::Image< itk::DiffusionTensor3D<TPixel>, VImageDimension>, itk::Image< itk::DiffusionTensor3D<TPixel>, VImageDimension> >
{
public:
  /** Standard class typedefs. */
  typedef MFMTensorResampleImageFilter                Self;
  typedef typename itk::Image< itk::DiffusionTensor3D<TPixel>, VImageDimension> TInputImage;
  typedef typename itk::Image< itk::DiffusionTensor3D<TPixel>, VImageDimension> TOutputImage;
  typedef typename itk::VectorImage< float, 3 >				FractionsImageType;
  typedef itk::ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef itk::SmartPointer<Self>                            Pointer;
  typedef itk::SmartPointer<const Self>                      ConstPointer;

  typedef TInputImage                           InputImageType;
  typedef TOutputImage                          OutputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename OutputImageType::Pointer     OutputImagePointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;

  typedef itk::ContinuousIndex<double, VImageDimension> ContinuousIndexType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(MFMTensorResampleImageFilter, ImageToImageFilter);

  /** Number of dimensions. */
  itkStaticConstMacro(ImageDimension, unsigned int, TOutputImage::ImageDimension);

  /** Transform typedef.
   *
   * \todo Check that input and output images have the same number of 
   * dimensions; this is required by the current implementation of 
   * AffineTransform. */
  typedef itk::Transform<double, itkGetStaticConstMacro(ImageDimension), itkGetStaticConstMacro(ImageDimension)> TransformType;
  typedef typename TransformType::ConstPointer TransformPointerType;

  /** Image size typedef. */
  typedef itk::Size<itkGetStaticConstMacro(ImageDimension)> SizeType;

  /** Image index typedef. */
  typedef typename TOutputImage::IndexType IndexType;

  /** Image pixel value typedef. */
  typedef typename TOutputImage::PixelType   PixelType;
  typedef typename PixelType::ValueType      PixelComponentType;

  /** Typedef to describe the output image region type. */
  typedef typename TOutputImage::RegionType OutputImageRegionType;

  /** Image spacing,origin and direction typedef */
  typedef typename TOutputImage::SpacingType   SpacingType;
  typedef typename TOutputImage::PointType     OriginPointType;
  typedef typename TOutputImage::DirectionType DirectionType;
  typedef typename TOutputImage::PointType PointType;
  

	typedef enum {
		  INTERP_NONE = 0,
		  INTERP_EUCLIDEAN = 1,
		  INTERP_LOGEUCLIDEAN = 2
	} InterpolationModeValue;
	typedef unsigned int InterpolationModeType;


  /** Set the coordinate transformation.
   * Set the coordinate transform to use for resampling.  Note that this
   * must be in index coordinates and is the output-to-input transform,
   * NOT the input-to-output transform that you might naively expect.
   * The default is itk::AffineTransform<TInterpolatorPrecisionType, ImageDimension>. */
  itkSetConstObjectMacro( Transform, TransformType ); 

  /** Get a pointer to the coordinate transform. */
  itkGetConstObjectMacro( Transform, TransformType );

  /** Set the size of the output image. */
  itkSetMacro( Size, SizeType );

  /** Get the size of the output image. */
  itkGetConstReferenceMacro( Size, SizeType );
     
  /** Set the pixel value when a transformed pixel is outside of the
   * image.  The default default pixel value is 0. */
  itkSetMacro(DefaultPixelValue,PixelType);

  /** Get the pixel value when a transformed pixel is outside of the image */
  itkGetConstMacro(DefaultPixelValue,PixelType);

  /** Set the output image spacing. */
  itkSetMacro(OutputSpacing, SpacingType);
  virtual void SetOutputSpacing( const double* values);

  /** Get the output image spacing. */
  itkGetConstReferenceMacro( OutputSpacing, SpacingType );

  /** Set the output image origin. */
  itkSetMacro(OutputOrigin, PointType);
  virtual void SetOutputOrigin( const double* values);

  /** Get the output image origin. */
  itkGetConstReferenceMacro( OutputOrigin, PointType );

  /** Set the output direciton cosine matrix. */
  itkSetMacro(OutputDirection, DirectionType);
  itkGetConstReferenceMacro(OutputDirection, DirectionType);

  /** Set the start index of the output largest possible region. 
   * The default is an index of all zeros. */
  itkSetMacro( OutputStartIndex, IndexType );

  /** The interpolation mode
  */
  itkSetMacro( TensorInterpolationMode, InterpolationModeType );
  itkGetConstMacro( TensorInterpolationMode, InterpolationModeType );

  /** Output the progress using crlTaskProgress
  */
  itkSetMacro( ShowProgress, bool );
  itkGetConstMacro( ShowProgress, bool );

  /** Get the start index of the output largest possible region. */
  itkGetConstReferenceMacro( OutputStartIndex, IndexType );

  itkSetMacro( FractionsImage, FractionsImageType::Pointer );
  itkGetMacro( ResampledFractionsImage, FractionsImageType::Pointer);


  /** VectorResampleImageFilter produces an image which is a different size
   * than its input.  As such, it needs to provide an implementation
   * for GenerateOutputInformation() in order to inform the pipeline
   * execution model.  The original documentation of this method is
   * below. \sa ProcessObject::GenerateOutputInformaton() */
  virtual void GenerateOutputInformation();

  /** VectorResampleImageFilter needs a different input requested region than
   * the output requested region.  As such, VectorResampleImageFilter needs
   * to provide an implementation for GenerateInputRequestedRegion()
   * in order to inform the pipeline execution model.
   * \sa ProcessObject::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion();

  /** This method is used to set the state of the filter before 
   * multi-threading. */
  virtual void BeforeThreadedGenerateData();

  /** This method is used to set the state of the filter after 
   * multi-threading. */
  virtual void AfterThreadedGenerateData();

  /** Method Compute the Modified Time based on changed to the components. */
  unsigned long GetMTime( void ) const;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<typename TInputImage::PixelType::ValueType>));
  itkConceptMacro(OutputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<PixelComponentType>));
  /** End concept checking */
#endif

protected:
  MFMTensorResampleImageFilter();
  ~MFMTensorResampleImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** VectorResampleImageFilter can be implemented as a multithreaded filter.  Therefore,
   * this implementation provides a ThreadedGenerateData() routine which
   * is called for each processing thread. The output image data is allocated
   * automatically by the superclass prior to calling ThreadedGenerateData().
   * ThreadedGenerateData can only write to the portion of the output image
   * specified by the parameter "outputRegionForThread"
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData() */
  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            int threadId );

  int SelectMatchingTensor(  const itk::Vector< double, 3> &refVector, const IndexType index, PixelType &out ) ;
  void EstimatePrimaryEigenvector( const PixelType &in, itk::Vector< double, 3> &evec) const;
  void InterpolateMatchingTensor( int idTensorIn, const PixelType &in, const ContinuousIndexType &cindex, PixelType &out, typename FractionsImageType::PixelType &outFractions ) ;

  void RotateOneTensor(const PixelType& intensor, const vnl_matrix<double>& R, PixelType& outTensor );

private:
  MFMTensorResampleImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  SizeType                m_Size;       // Size of the output image
  TransformPointerType    m_Transform;  // Coordinate transform to use
  PixelType               m_DefaultPixelValue; // default pixel value if the point is outside the image
  SpacingType             m_OutputSpacing; // output image spacing
  PointType               m_OutputOrigin;  // output image origin
  DirectionType           m_OutputDirection; // output image direction cosines
  IndexType               m_OutputStartIndex; // output start index
  InterpolationModeType	  m_TensorInterpolationMode;

  FractionsImageType::Pointer m_FractionsImage;
  FractionsImageType::Pointer m_ResampledFractionsImage;

  itk::SimpleFastMutexLock mutexAccessSVD;
  bool				m_ShowProgress;
};

  
} // end namespace crl
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "crlMFMTensorResampleImageFilter.txx"
#endif
  
#endif
