
#ifndef _itkTensorExpFunctor_h
#define _itkTensorExpFunctor_h 1

#include <itkUnaryFunctorImageFilter.h>
#include <itkSymmetricSecondRankTensor.h>
#include <math.h>

#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_det.h>

#include "crlFastOps.h"

#define TENSOR_EXP_FAST_VERSION 1

namespace itk
{

namespace Functor {

template< class TInputTensor, class TOutputTensor >
class TensorExp
{
  public:
  typedef typename TInputTensor::ValueType ValueType;
  typedef typename TInputTensor::EigenVectorsMatrixType MatrixType;
  typedef typename TInputTensor::EigenValuesArrayType EigenValuesType;

  TensorExp() {};
  ~TensorExp() {};


#ifdef TENSOR_EXP_FAST_VERSION

  typedef vnl_matrix<double>					VnlMatrixType;
  typedef vnl_matrix<double>					VnlTensorType;
  typedef vnl_diag_matrix<double>				VnlEigValueType;
  typedef vnl_matrix<double>					VnlEigVectorType;
  typedef itk::SymmetricEigenAnalysis< VnlTensorType , VnlEigValueType, VnlEigVectorType > SymmetricEigenAnalysisType;

  inline void Tensor2VnlMatrix(VnlTensorType &out, const TInputTensor &tensor ) const
  {
	  //vnl_matrix: m[row][col]

	  double *o = out.data_block();
	  const typename TInputTensor::ValueType *p = tensor.GetDataPointer ();

	  // Fast copy
	  *(o) = *(p);
	  *(o+1) = *(p+1);
	  *(o+2) = *(p+2);
	  *(o+3) = *(p+1);
	  *(o+4) = *(p+3);
	  *(o+5) = *(p+4);
	  *(o+6) = *(p+2);
	  *(o+7) = *(p+4);
	  *(o+8) = *(p+5);
  }

  inline void VnlMatrix2Tensor( TOutputTensor &out , const VnlMatrixType &m ) 
  {
	  typename TOutputTensor::ValueType *p = out.GetDataPointer ();
	  *p = m(0,0);
	  *(p+1) =  m(0,1);
	  *(p+2) = m(0,2);
	  *(p+3) = m(1,1);
	  *(p+4) = m(1,2);
	  *(p+5) = m(2,2);
  }

  inline TOutputTensor operator()( const TInputTensor & x )
  {
	VnlTensorType D(3,3);
	Tensor2VnlMatrix(D,x); 

	// Compute the eigen values/vectors
	VnlEigVectorType		eVec(3,3);
	VnlEigValueType			eVal(3) ;

	SymmetricEigenAnalysisType EigenAnalysis(3);
	EigenAnalysis.ComputeEigenValuesAndVectors(D, eVal, eVec);

	for ( int i=0 ; i<3 ; i++ )	eVal[i] =  exp ( eVal[i] );

	crl::FastOps::AtBA_Bdiag_3x3(D, eVec, eVal); 


    TOutputTensor result;
	VnlMatrix2Tensor(result, D);

    return result;
  }

#else

  inline TOutputTensor operator()( const TInputTensor & x )
  {
    MatrixType V;
    EigenValuesType e;
    TOutputTensor result;

    // Compute the matrix exponential by, 
    //  1. diagonalizing the tensor
    //  2. computing the natural exponential of the diagonal
    //  3. reversing the diagonalization

    // Each ROW of the matrix V represents an eigenvector, which is different
    // to matlab. In matlab each column is an eigenvector.
    x.ComputeEigenAnalysis(e, V);

    V = V.GetTranspose();  // Make the columns be the eigenvectors

    MatrixType xin; // Convert from symmetric tensor class to matrix
    xin(0,0) = x[0];  xin(0,1) = x[1];  xin(0,2) = x[2];
    xin(1,0) = x[1];  xin(1,1) = x[3];  xin(1,2) = x[4];
    xin(2,0) = x[2];  xin(2,1) = x[4];  xin(2,2) = x[5];

    MatrixType diag;
    diag = V.GetInverse(); // diag <- V^{-1}
    diag = diag*xin;       // diag <- V^{-1}* x
    diag = diag * V;       // diag <- V^{-1}* x * V

    diag(0,0) = exp(diag(0,0));
    diag(1,1) = exp(diag(1,1));
    diag(2,2) = exp(diag(2,2));

    MatrixType xout;
    xout = V * diag;
    xout = xout * V.GetInverse();

    // Convert from matrix type to symmetric tensor output
    result[0] = xout(0,0); result[1] = xout(0,1);  result[2] = xout(0,2);
    result[3] = xout(1,1); result[4] = xout(1,2);
    result[5] = xout(2,2);

    return result;
  }
#endif

};

} // namespace Functor

template <class TInputTensorImage, class TOutputTensorImage >
class ITK_EXPORT TensorExpImageFilter :
public 
UnaryFunctorImageFilter<TInputTensorImage,TOutputTensorImage,
    Functor::TensorExp< 
        typename TInputTensorImage::PixelType, 
        typename TOutputTensorImage::PixelType> >
{
  public:
  /** Standard class typedefs. */
  typedef TensorExpImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputTensorImage, TOutputTensorImage,
    Functor::TensorExp<
        typename TInputTensorImage::PixelType,
        typename TOutputTensorImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  protected:
    TensorExpImageFilter() {}
    virtual ~TensorExpImageFilter() {}

  private:
    TensorExpImageFilter(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#endif

