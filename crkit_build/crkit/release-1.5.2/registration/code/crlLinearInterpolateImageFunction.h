#ifndef __crlLinearInterpolateImageFunction_h
#define __crlLinearInterpolateImageFunction_h

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"

#if defined( ITK_USE_OPTIMIZED_REGISTRATION_METHODS ) && defined( LINEAR_INTERPOLATOR_FIXED )
#include "itkOptLinearInterpolateImageFunction.h"
#else

#include "itkInterpolateImageFunction.h"

namespace crl
{

/** \class LinearInterpolateImageFunction
 * \brief Linearly interpolate an image at specified positions.
 *
 * LinearInterpolateImageFunction linearly interpolates image intensity at
 * a non-integer pixel position. This class is templated
 * over the input image type and the coordinate representation type 
 * (e.g. float or double).
 *
 * This function works for N-dimensional images.
 * This function does bounds checking.
 *
 * \warning This function work only for images with scalar pixel
 * types. For vector images use VectorLinearInterpolateImageFunction.
 *
 * \sa VectorLinearInterpolateImageFunction
 *
 * \ingroup ImageFunctions ImageInterpolators 
 */
template <class TInputImage, class TCoordRep = double>
class ITK_EXPORT LinearInterpolateImageFunction : 
  public itk::InterpolateImageFunction<TInputImage,TCoordRep> 
{
public:
  /** Standard class typedefs. */
  typedef LinearInterpolateImageFunction                  Self;
  typedef itk::InterpolateImageFunction<TInputImage,TCoordRep> Superclass;
  typedef itk::SmartPointer<Self>                              Pointer;
  typedef itk::SmartPointer<const Self>                        ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(crlLinearInterpolateImageFunction, itk::InterpolateImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** OutputType typedef support. */
  typedef typename Superclass::OutputType OutputType;

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;

  /** RealType typedef support. */
  typedef typename Superclass::RealType RealType;

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the linearly interpolated image intensity at a 
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex( 
    const ContinuousIndexType & index ) const;

protected:
  LinearInterpolateImageFunction();
  ~LinearInterpolateImageFunction(){};
  void PrintSelf(std::ostream& os, itk::Indent indent) const;

private:
  LinearInterpolateImageFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  /** Number of neighbors used in the interpolation */
  static const unsigned long  m_Neighbors;  

};

} // end namespace crl

// Define instantiation macro for this template.
#define ITK_TEMPLATE_crlLinearInterpolateImageFunction(_, EXPORT, x, y) namespace crl { \
  _(2(class EXPORT LinearInterpolateImageFunction< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef LinearInterpolateImageFunction< ITK_TEMPLATE_2 x > \
                                                  LinearInterpolateImageFunction##y; } \
  }

#if ITK_TEMPLATE_TXX
# include "crlLinearInterpolateImageFunction.txx"
#endif

#endif

#endif
