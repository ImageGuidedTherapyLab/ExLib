
#ifndef __itkDisplacementFieldJacobianRecursiveGaussianImageFilter_h
#define __itkDisplacementFieldJacobianRecursiveGaussianImageFilter_h

#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImage.h"
#include "itkMatrix.h"
#include "itkPixelTraits.h"
#include "itkCommand.h"
#include "itkProgressAccumulator.h"


namespace itk
{

/** \class DisplacementFieldJacobianRecursiveGaussianImageFilter
 * \brief Computes the Jacobian matrix of a dispalcement field image by convolution
 *        with the derivative of a Gaussian
 * 
 * This filter is implemented using the recursive gaussian
 * filters
 *
 * 
 * \ingroup GradientFilters   
 * \ingroup Singlethreaded
 */
// NOTE that the ITK_TYPENAME macro has to be used here in lieu 
// of "typename" because VC++ doesn't like the typename keyword 
// on the defaults of template parameters
template <typename TInputImage, 
          typename TOutputImage= Image< Matrix< 
  ITK_TYPENAME TInputImage::PixelType::ComponentType,
  ::itk::GetImageDimension<TInputImage>::ImageDimension, ::itk::GetImageDimension<TInputImage>::ImageDimension >,
                                        ::itk::GetImageDimension<TInputImage>::ImageDimension > >
class ITK_EXPORT DisplacementFieldJacobianRecursiveGaussianImageFilter:
    public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef DisplacementFieldJacobianRecursiveGaussianImageFilter  Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                   Pointer;
  typedef SmartPointer<const Self>        ConstPointer;
  
  
  /** Pixel Type of the input image */
  typedef TInputImage                                    InputImageType;
  typedef typename TInputImage::PixelType                PixelType;
  typedef typename PixelType::ComponentType ComponentType;


  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      ::itk::GetImageDimension<TInputImage>::ImageDimension);

  /** Scalar image typedef, for the output of the RecursiveGaussianImageFilter */
  typedef Image<ComponentType, ImageDimension>                                    ScalarImageType;


  /**  Output Image Nth Element Adaptor
   *  This adaptor allows to use conventional scalar 
   *  smoothing filters to compute each one of the 
   *  components of the gradient image pixels. */
  typedef VectorIndexSelectionCastImageFilter< TInputImage,
                                  ScalarImageType >  VectorScalarSelectionFilterType;
  typedef typename VectorScalarSelectionFilterType::Pointer VectorScalarSelectionFilterTypePointer;

  
  /**  Derivative filter type */
  typedef GradientRecursiveGaussianImageFilter<
                                ScalarImageType,
                                InputImageType
                                        >    GradientFilterType;


  /**  Pointer to a derivative filter.  */
  typedef typename GradientFilterType::Pointer  GradientFilterPointer;

  
  /**  Pointer to the Output Image */
  typedef typename TOutputImage::Pointer          OutputImagePointer;                                  


  /** Type of the output Image */
  typedef TOutputImage      OutputImageType;
  typedef typename          OutputImageType::PixelType      OutputPixelType;
  typedef typename OutputPixelType::ComponentType  OutputComponentType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set Sigma value. Sigma is measured in the units of image spacing.  */
  void SetSigma( ComponentType sigma );

  /** Define which normalization factor will be used for the Gaussian */
  void SetNormalizeAcrossScale( bool normalizeInScaleSpace );
  itkGetMacro( NormalizeAcrossScale, bool );

  /** Define which normalization factor will be used for the Gaussian */
  void SetComputeRotationMatrix( bool computeR );
  itkGetMacro( ComputeRotationMatrix, bool );
  
  /** DisplacementFieldJacobianRecursiveGaussianImageFilter needs all of the input to produce an
   * output. Therefore, DisplacementFieldJacobianRecursiveGaussianImageFilter needs to provide
   * an implementation for GenerateInputRequestedRegion in order to inform
   * the pipeline execution model.
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion() throw(InvalidRequestedRegionError);

protected:
  
  DisplacementFieldJacobianRecursiveGaussianImageFilter();
  virtual ~DisplacementFieldJacobianRecursiveGaussianImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Generate Data */
  void GenerateData( void );

  // Override since the filter produces the entire dataset
  void EnlargeOutputRequestedRegion(DataObject *output);

  
private:

  DisplacementFieldJacobianRecursiveGaussianImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  GradientFilterPointer      m_GradientFilterU;
  GradientFilterPointer      m_GradientFilterV;
  GradientFilterPointer      m_GradientFilterW;
  VectorScalarSelectionFilterTypePointer     m_ImageAdaptorU;
  VectorScalarSelectionFilterTypePointer     m_ImageAdaptorV;
  VectorScalarSelectionFilterTypePointer     m_ImageAdaptorW;

  /** Normalize the image across scale space */
  bool m_NormalizeAcrossScale; 
  bool m_ComputeRotationMatrix;

  ProgressAccumulator::Pointer  m_Progress;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDisplacementFieldJacobianRecursiveGaussianImageFilter.txx"
#endif

#endif




