
#ifndef _itkTensorScalarMultiplyFunctor_h
#define _itkTensorScalarMultiplyFunctor_h 1

#include <itkUnaryFunctorImageFilter.h>
#include <itkSymmetricSecondRankTensor.h>
#include <math.h>

namespace itk
{

namespace Functor {

template< class TInputTensor, class TOutputTensor >
class TensorScalarMultiply
{
  public:
  typedef typename TInputTensor::ValueType ValueType;
  typedef typename TInputTensor::EigenVectorsMatrixType MatrixType;
  typedef typename TInputTensor::EigenValuesArrayType EigenValuesType;

  TensorScalarMultiply() { scalar = 1.0;};
  ~TensorScalarMultiply() {};

  inline TOutputTensor operator()( const TInputTensor & x )
  {
    TOutputTensor result;

    result = x * scalar;

    return result;
  }

  void SetScalar( ValueType s) { scalar = s; }

  private:
  ValueType scalar;

};

} // namespace Functor

template <class TInputTensorImage, class TOutputTensorImage >
class ITK_EXPORT TensorScalarMultiplyImageFilter :
public 
UnaryFunctorImageFilter<TInputTensorImage,TOutputTensorImage,
    Functor::TensorScalarMultiply< 
        typename TInputTensorImage::PixelType, 
        typename TOutputTensorImage::PixelType> >
{
  public:
  /** Standard class typedefs. */
  typedef TensorScalarMultiplyImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputTensorImage, TOutputTensorImage,
    Functor::TensorScalarMultiply<
        typename TInputTensorImage::PixelType,
        typename TOutputTensorImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  void SetScalar(typename TInputTensorImage::PixelType::ValueType s) {
    scalar = s;
    this->GetFunctor().SetScalar( this->scalar );
  }

  void BeforeThreadedGenerateData(void);

  protected:
    TensorScalarMultiplyImageFilter() { scalar = 1.0;}
    virtual ~TensorScalarMultiplyImageFilter() {}

  private:
    TensorScalarMultiplyImageFilter(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

    typename TInputTensorImage::PixelType::ValueType scalar;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "TensorScalarMultiplyFunctor.txx"
#endif

#endif

