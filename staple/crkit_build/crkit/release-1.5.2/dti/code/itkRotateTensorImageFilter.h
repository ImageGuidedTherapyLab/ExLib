#ifndef __itkRotateTensorImageFilter_h
#define __itkRotateTensorImageFilter_h

#include "itkUnaryFunctorImageFilter.h"

#include "itkVector.h"
#include "itkMatrix.h"
#include "itkSymmetricSecondRankTensor.h"

namespace itk
{

/** This functor class rotates the basis of the tensor-pixels.
 *   Tensors can be represented as:
 *   -  itk::Vector<ValueType, 6>, 
 *   -  itk::Matrix<ValueType, 3,3> and
 *   -  itk::SymmetricSecondRankTensor<ValueType,3> objects.  
 *
 *    The rotation is done using the following formula: 
 *    T_new = Rt * T * R, 
 *     where R is the rotation matrix, and Rt is the transposed Rotation Matrix. */

namespace Functor {

template< typename TInput, typename  TOutput>
class Rotate
{
public:
  //typedef typename NumericTraits< typename TInput::ValueType >::RealType RealType;
  typedef typename TInput::ValueType ValueType;
  Rotate() {}
  ~Rotate() {}
  void SetRotationMatrix( Matrix<ValueType,3,3>  v ) { m_RotationMatrix = v; }
  //itkStaticConstMacro(VectorDimension,unsigned int,TInput::Dimension);
  
  typedef Vector<ValueType,6>  VectorType;
  typedef Matrix<ValueType,3,3>  MatrixType;
  typedef SymmetricSecondRankTensor<ValueType,3>  TensorType;
  
  inline TOutput operator()( const VectorType & x )
  {

    MatrixType xin;
    xin(0,0) = x[0];  xin(0,1) = x[1];  xin(0,2) = x[2];
    xin(1,0) = x[1];  xin(1,1) = x[3];  xin(1,2) = x[4];
    xin(2,0) = x[2];  xin(2,1) = x[4];  xin(2,2) = x[5];
    
    //Ttmp = T*R
    MatrixType x_tmp;
    x_tmp = xin*m_RotationMatrix;

    
    //Trotated = Rt*T*R = Rt*Ttmp
    MatrixType xout;
    xout = m_RotationMatrix.GetTranspose()*x_tmp.GetVnlMatrix();
    
    TOutput  result;
    result[0] = xout(0,0);  result[1] = xout(0,1);  result[2] = xout(0,2);
    result[3] = xout(1,1); result[4] = xout(1,2);
    result[5] = xout(2,2);
 //if (x[1]>0)
// std::cout << x << "\t" << result << std::endl;
    return result;
  }
  
  inline TensorType operator()( const TensorType & x )
  {

    MatrixType xin;
    for (unsigned int i=0; i<3; i++)
      for(unsigned int j=0; j<3; j++)
        xin(i,j) = x(i,j);
    
    //Ttmp = T*R
    MatrixType x_tmp;
    x_tmp = xin*m_RotationMatrix;

    
    //Trotated = Rt*T*R = Rt*Ttmp
    MatrixType xout;
    xout = m_RotationMatrix.GetTranspose()*x_tmp.GetVnlMatrix();
    
    TensorType  result;
    for (unsigned int i=0; i<3; i++)
      for(unsigned int j=i; j<3; j++)
        result(i,j) = x(i,j);
 if (x[1]>0)
 //std::cout << x << "\t" << result << std::endl;

    return result;
  }
  
  inline MatrixType operator()( const MatrixType & x )
  {
    
    MatrixType  x_tmp;
    //Ttmp = T*R
    x_tmp = x*m_RotationMatrix;
    
    //Trotated = Rt*T*R = Rt*Ttmp
    MatrixType result;
    result = m_RotationMatrix.GetTranspose()*x_tmp.GetVnlMatrix();

    return result;
  }

private:
  MatrixType m_RotationMatrix;
};

}  // end namespace functor


/** \class RotateTensorImageFilter
 * \brief Applies a linear transformation to the magnitude of pixel vectors in a
 * vector Image.
 *
 * RotateTensorImageFilter applies pixel-wise a linear transformation
 * to the intensity values of input image pixels.
 *
 * All computations are performed in the precison of the input pixel's
 * ValueType.
 *
 * \sa Rotate
 *
 * \ingroup IntensityImageFilters  Multithreaded
 *
 */
template <typename  TInputImage, typename  TOutputImage=TInputImage>
class ITK_EXPORT RotateTensorImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage,
                        Functor::Rotate<
  typename TInputImage::PixelType,
  typename TOutputImage::PixelType>   >
{
public:
  /** Standard class typedefs. */
  typedef RotateTensorImageFilter  Self;
  typedef UnaryFunctorImageFilter<TInputImage,
                                  TOutputImage,
                                  Functor::Rotate< typename TInputImage::PixelType, typename TOutputImage::PixelType> >  Superclass;
  
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  typedef typename TOutputImage::PixelType  OutputPixelType;
  typedef typename TInputImage::PixelType     InputPixelType;
  typedef typename InputPixelType::ValueType        InputValueType;
  typedef typename OutputPixelType::ValueType     OutputValueType;


  typedef Matrix<InputValueType,3,3>  RotationMatrixType;
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  itkSetMacro(RotationMatrix, RotationMatrixType);
//  void SetRotationMatrix( RotationMatrixType v);
  itkGetConstReferenceMacro( RotationMatrix, RotationMatrixType );

  void PrintParameters();

  /** Process to execute before entering the multithreaded section */
  void BeforeThreadedGenerateData(void);

  /** Print internal ivars */
  void PrintSelf(std::ostream& os, Indent indent) const;

protected:
  RotateTensorImageFilter();
  virtual ~RotateTensorImageFilter() {};

private:
  RotateTensorImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  RotationMatrixType m_RotationMatrix;

};



} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRotateTensorImageFilter.txx"
#endif

#endif
