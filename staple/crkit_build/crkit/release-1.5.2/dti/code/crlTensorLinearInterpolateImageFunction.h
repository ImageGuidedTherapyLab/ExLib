#ifndef __CRL_TENSORLINEARINTERPOLATEIMAGEFUNCTION_H
#define __CRL_TENSORLINEARINTERPOLATEIMAGEFUNCTION_H

#include "crlTensorInterpolateImageFunction.h"

namespace crl
{

/** 
 * \class TensorLinearInterpolateImageFunction
 * \brief Linearly interpolate a tensor or a log tensor image.
 *
 * TensorLinearInterpolateImageFunction linearly interpolates an
 * image intensity non-integer pixel position. This class is templated
 * over the input image type and the coordinate representation type.
 *
 * This function works for N-dimensional images.
 *
 * \warning This function work only for tensor images.
 *
 * \ingroup ImageFunctions ImageInterpolators
 * 
 */
template <class TInputImage, class TCoordRep = double>
class ITK_EXPORT TensorLinearInterpolateImageFunction : 
  public TensorInterpolateImageFunction<TInputImage,TCoordRep> 
{
public:
  /** Standard class typedefs. */
  typedef TensorLinearInterpolateImageFunction Self;
  typedef TensorInterpolateImageFunction<TInputImage,TCoordRep> Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(TensorLinearInterpolateImageFunction, 
    TensorInterpolateImageFunction);

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;
  typedef typename Superclass::PixelType      PixelType;
  typedef typename Superclass::ValueType      ValueType;
  typedef typename Superclass::RealType       RealType;
  typedef typename Superclass::PointType       PointType;
    
  /** Grab the vector dimension from the superclass. */
  itkStaticConstMacro(Dimension, unsigned int,
                       Superclass::Dimension);

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Output type is Tensor<double,Dimension> */
  typedef typename Superclass::OutputType OutputType;

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the linearly interpolated image intensity at a 
   * specified point position. No bounds checking is done.
   * The point and neighborhood is assumed to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex( 
    const ContinuousIndexType & index ) const;

protected:
  TensorLinearInterpolateImageFunction();
  ~TensorLinearInterpolateImageFunction(){};
  void PrintSelf(std::ostream& os, itk::Indent indent) const;

private:
  TensorLinearInterpolateImageFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Number of neighbors used in the interpolation */
  static const unsigned long  m_Neighbors;  

};

} // end namespace crl

// Define instantiation macro for this template.
#define ITK_TEMPLATE_TensorLinearInterpolateImageFunction(_, EXPORT, x, y) namespace crl { \
  _(2(class EXPORT TensorLinearInterpolateImageFunction< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef TensorLinearInterpolateImageFunction< ITK_TEMPLATE_2 x > \
                                                  TensorLinearInterpolateImageFunction##y; } \
  }

#if ITK_TEMPLATE_TXX
# include "crlTensorLinearInterpolateImageFunction.txx"
#endif

#endif
