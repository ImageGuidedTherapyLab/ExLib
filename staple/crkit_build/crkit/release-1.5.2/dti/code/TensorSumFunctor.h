
#ifndef _itkTensorSumFunctor_h
#define _itkTensorSumFunctor_h 1

#include <itkBinaryFunctorImageFilter.h>
#include <itkSymmetricSecondRankTensor.h>
#include <math.h>

namespace itk
{

namespace Functor {

template< class TInputTensor1, class TInputTensor2, class TOutputTensor >
class TensorSum
{
  public:
  typedef typename TInputTensor1::ValueType ValueType;
  typedef typename TInputTensor1::EigenVectorsMatrixType MatrixType;
  typedef typename TInputTensor1::EigenValuesArrayType EigenValuesType;

  TensorSum() {};
  ~TensorSum() {};

  inline TOutputTensor operator()( const TInputTensor1 & x, 
                                   const TInputTensor2 & y )
  {
    TOutputTensor result;

    result = x + y;

    return result;
  }

};

} // namespace Functor

template <class TInputTensorImage1, class TInputTensorImage2, class TOutputTensorImage >
class ITK_EXPORT TensorSumImageFilter :
public 
BinaryFunctorImageFilter<TInputTensorImage1,TInputTensorImage2,TOutputTensorImage,
    Functor::TensorSum< 
        typename TInputTensorImage1::PixelType, 
        typename TInputTensorImage2::PixelType, 
        typename TOutputTensorImage::PixelType> >
{
  public:
  /** Standard class typedefs. */
  typedef TensorSumImageFilter Self;
  typedef BinaryFunctorImageFilter<TInputTensorImage1, TInputTensorImage2,
                                      TOutputTensorImage,
    Functor::TensorSum<
        typename TInputTensorImage1::PixelType,
        typename TInputTensorImage2::PixelType,
        typename TOutputTensorImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  protected:
    TensorSumImageFilter() {}
    virtual ~TensorSumImageFilter() {}

  private:
    TensorSumImageFilter(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#endif

