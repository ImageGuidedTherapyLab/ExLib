///itkRotateTensorsPixelwiseImageFilter.h

#ifndef __itkRotateTensorsPixelwiseImageFilter_h
#define __itkRotateTensorsPixelwiseImageFilter_h

#include "itkBinaryFunctorImageFilter.h"
#include "itkVector.h"
#include "itkMatrix.h"
#include "itkSymmetricSecondRankTensor.h"


namespace itk
{
  
/** \class RotateTensorsPixelwiseImageFilter
 * \brief Implements an operator for pixel-wise rotation of a tensor image.
 *
 *
 * \ingroup IntensityImageFilters  Multithreaded
 */
namespace Functor {  
  
template< class TInputTensor, class TRotationMatrix, class TOutputTensor >
class ReorientTensor
{
public:
  
  typedef typename TInputTensor::ValueType ValueType;
  typedef TRotationMatrix  MatrixType;
//  typedef Vector<ValueType,6>  VectorType;
//  typedef SymmetricSecondRankTensor<ValueType,3>  TensorType;

  ReorientTensor() {};
  ~ReorientTensor() {};
  inline TOutputTensor operator()( const TInputTensor & x, const TRotationMatrix & r)
  {
    MatrixType xin;
    xin(0,0) = x[0];  xin(0,1) = x[1];  xin(0,2) = x[2];
    xin(1,0) = x[1];  xin(1,1) = x[3];  xin(1,2) = x[4];
    xin(2,0) = x[2];  xin(2,1) = x[4];  xin(2,2) = x[5];
    
    //Ttmp = T*R
    MatrixType x_tmp;
    x_tmp = xin*r;

    
    //Trotated = Rt*T*R = Rt*Ttmp
    MatrixType xout;
    xout = r.GetTranspose()*x_tmp.GetVnlMatrix();
    
    TOutputTensor  result;
    result[0] = xout(0,0);  result[1] = xout(0,1);  result[2] = xout(0,2);
    result[3] = xout(1,1); result[4] = xout(1,2);
    result[5] = xout(2,2);
//std::cout << r << std::endl;    
    return result;
  }
    
  
}; 

}
template <class TInputTensorImage1, class TInputImage2, class TOutputTensorImage>
class ITK_EXPORT RotateTensorsPixelwiseImageFilter :
    public
BinaryFunctorImageFilter<TInputTensorImage1,TInputImage2,TOutputTensorImage, 
                         Functor::ReorientTensor< 
  typename TInputTensorImage1::PixelType, 
  typename TInputImage2::PixelType,
  typename TOutputTensorImage::PixelType>   >


{
public:
  /** Standard class typedefs. */
  typedef RotateTensorsPixelwiseImageFilter  Self;
  typedef BinaryFunctorImageFilter<TInputTensorImage1,TInputImage2,TOutputTensorImage, 
                                   Functor::ReorientTensor< 
    typename TInputTensorImage1::PixelType, 
    typename TInputImage2::PixelType,
    typename TOutputTensorImage::PixelType>   
  >  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
protected:
  RotateTensorsPixelwiseImageFilter() {}
  virtual ~RotateTensorsPixelwiseImageFilter() {}

private:
  RotateTensorsPixelwiseImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk


#endif
