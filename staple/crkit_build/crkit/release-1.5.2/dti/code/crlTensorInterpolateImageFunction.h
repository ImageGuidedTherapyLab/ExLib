#ifndef __CRL_TENSORINTERPOLATEIMAGEFUNCTION_H
#define __CRL_TENSORINTERPOLATEIMAGEFUNCTION_H

#include <itkImageFunction.h>
#include <itkFixedArray.h>
#include <itkSymmetricSecondRankTensor.h>

namespace crl
{

/**
 * Due to a bug in MSVC, an enum value cannot be accessed out of a template
 * parameter until the template class opens.  In order for templated classes
 * to access the dimension of a template parameter in defining their
 * own dimension, this class is needed as a work-around.
 */
template <typename T>
struct GetDimension
{
  itkStaticConstMacro(Dimension, int, T::Dimension);
}; 

  
/** \class TensorInterpolateImageFunction
 * \brief Base class for all tensor image interpolaters.
 *
 * TensorInterpolateImageFunction is the base for all ImageFunctions that
 * interpolates image with tensor pixel types. This function outputs
 * a return value of type SymmetricSecondRankTensor<double,Dimension>.
 *
 * This class is templated input image type and the coordinate
 * representation type.
 *
 * \warning This hierarchy of functions work only for images 
 * with tensor-based pixel types.
 * 
 * \sa InterpolateImageFunction
 * \ingroup ImageFunctions ImageInterpolators
 */
template <class TInputImage, class TCoordRep = double>
class ITK_EXPORT TensorInterpolateImageFunction :
  public itk::ImageFunction<
    TInputImage,
    typename TInputImage::PixelType,
    TCoordRep >
{
public:
  /** Extract the tensor dimension from the pixel template parameter. */
  itkStaticConstMacro(Dimension, unsigned int,
                      TInputImage::PixelType::InternalDimension);
  
  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard class typedefs. */
  typedef TensorInterpolateImageFunction Self;
  typedef itk::ImageFunction<TInputImage,
    typename TInputImage::PixelType,
    TCoordRep > Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(TensorInterpolateImageFunction, ImageFunction);

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;
  typedef typename InputImageType::PixelType  PixelType;
  typedef typename PixelType::ValueType       ValueType;
  typedef typename itk::NumericTraits<ValueType>::RealType  RealType;

  /** Point typedef support. */
  typedef typename Superclass::PointType PointType;

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Output type is RealType of TInputImage::PixelType. */
  typedef typename Superclass::OutputType OutputType;

  /** CoordRep typedef support. */
  typedef TCoordRep CoordRepType;

  /** Returns the interpolated image pixel value at a 
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType Evaluate( const PointType& point ) const
    {
    ContinuousIndexType index;
    this->GetInputImage()->TransformPhysicalPointToContinuousIndex( point, index );
    return ( this->EvaluateAtContinuousIndex( index ) );
    }

  /** Interpolate the image at a continuous index position
   *
   * Returns the interpolated image pixel value at a 
   * specified index position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * Subclasses must override this method.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex( 
    const ContinuousIndexType & index ) const = 0;

  /** Interpolate the image at an index position.
   * Simply returns the image pixel value at the
   * specified index position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtIndex( const IndexType & index ) const
    {
    OutputType output;
    PixelType input = this->GetInputImage()->GetPixel( index );
    for( unsigned int k = 0; k < Dimension; k++ )
      {
      output[k] = static_cast<double>( input[k] );
      }
    return ( output );
    }

protected:
  TensorInterpolateImageFunction() {}
  ~TensorInterpolateImageFunction() {}
  void PrintSelf(std::ostream& os, itk::Indent indent) const
    { Superclass::PrintSelf( os, indent ); }

private:
  TensorInterpolateImageFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace crl

#endif


