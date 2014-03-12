
#ifndef _itkTensorLogFunctor_h
#define _itkTensorLogFunctor_h 1

#include <itkUnaryFunctorImageFilter.h>
#include <itkSymmetricSecondRankTensor.h>
#include <math.h>
#include <float.h>

#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_det.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

#include "crlFastOps.h"

#define TENSOR_LOG_FAST_VERSION 1

namespace itk
{

namespace Functor {

template< class TInputTensor, class TOutputTensor >
class TensorLog
{
  public:
  typedef typename TInputTensor::ValueType ValueType;
  typedef typename TInputTensor::EigenVectorsMatrixType MatrixType;
  typedef typename TInputTensor::EigenValuesArrayType EigenValuesType;

  TensorLog() {};
  ~TensorLog() {};

#ifdef TENSOR_LOG_FAST_VERSION
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

    float minnonzero = FLT_EPSILON;   // Smallest positive number representable.
	for ( int i=0 ; i<3 ; i++ )	eVal[i] =   (    eVal[i] < minnonzero   ?    log(minnonzero) : log(eVal[i])   );

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
    // double minnonzero = 1e-14;
    float minnonzero = FLT_EPSILON;   // Smallest positive number representable.
    // float minflt = FLT_MIN;   // Smallest float representable.

    // Compute the matrix logarithm by, 
    //  1. diagonalizing the tensor
    //  2. computing the natural logarithm of the diagonal
    //  3. reversing the diagonalization

    // Invalid tensors, such as non-positive definite matrices, such be 
    // projected to become positive definite. They are projected to be the
    // positive definite matrix closest in Frobenius norm to the original
    // invalid tensor supplied.

    // Each ROW of the matrix V represents an eigenvector, which is different
    // to matlab. In matlab each column is an eigenvector.
    x.ComputeEigenAnalysis(e, V);
    // std::cout << "x is " << x << std::endl;

    V = V.GetTranspose();  // Make the columns be the eigenvectors
    // std::cout << "eigenvalues " << e << std::endl;
    // std::cout << "eigenvectors " << V << std::endl;
     
    MatrixType xin; // Convert from symmetric tensor class to matrix
    xin(0,0) = x[0];  xin(0,1) = x[1];  xin(0,2) = x[2];
    xin(1,0) = x[1];  xin(1,1) = x[3];  xin(1,2) = x[4];
    xin(2,0) = x[2];  xin(2,1) = x[4];  xin(2,2) = x[5];

    MatrixType diag;
    diag = V.GetInverse(); // diag <- V^{-1} Recall transpose = inverse
                             // for the orthogonal set of eigenvectors.
    diag = diag*xin;       // diag <- V^{-1}* x
    diag = diag * V;       // diag <- V^{-1}* x * V

    if ((diag(0,0) < FLT_EPSILON) || 
        (diag(1,1) < FLT_EPSILON) || 
        (diag(2,2) < FLT_EPSILON)) {
      // std::cerr << "Encountered a negative eigenvalue." << std::endl;
      // std::cerr << "diag matrix is " << std::endl << diag << std::endl;
      if (diag(0,0) < FLT_EPSILON) {
        diag(0,0) = minnonzero;
      }
      if (diag(1,1) < FLT_EPSILON) {
        diag(1,1) = minnonzero;
      }
      if (diag(2,2) < FLT_EPSILON) {
        diag(2,2) = minnonzero;
      }
    }
    diag(0,0) = log(diag(0,0));
    diag(1,1) = log(diag(1,1));
    diag(2,2) = log(diag(2,2));

    /* Enable this for numerical checking:
    if (isnan(diag(0,0)) || isnan(diag(1,1)) || isnan(diag(2,2)) ) {
      std::cout << "Created a NaN by taking log." << std::endl;
    }
    if (isinf(diag(0,0)) || isinf(diag(1,1)) || isinf(diag(2,2)) ) {
      std::cout << "Created an inf by taking log." << std::endl;
    }
    */
    
    // std::cout << "log diagonal matrix " << diag << std::endl;

    MatrixType xout;
    xout = V * diag;
    xout = xout * V.GetInverse();

    // std::cout << "log tensor matrix " << xout << std::endl;

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
class ITK_EXPORT TensorLogImageFilter :
public 
UnaryFunctorImageFilter<TInputTensorImage,TOutputTensorImage,
    Functor::TensorLog< 
        typename TInputTensorImage::PixelType, 
        typename TOutputTensorImage::PixelType> >
{
  public:
  /** Standard class typedefs. */
  typedef TensorLogImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputTensorImage, TOutputTensorImage,
    Functor::TensorLog<
        typename TInputTensorImage::PixelType,
        typename TOutputTensorImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  protected:
    TensorLogImageFilter() {}
    virtual ~TensorLogImageFilter() {}

  private:
    TensorLogImageFilter(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#endif

