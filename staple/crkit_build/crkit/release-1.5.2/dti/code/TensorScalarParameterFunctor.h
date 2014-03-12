
#ifndef _itkTensorScalarParamaterFunctor_h
#define _itkTensorScalarParamaterFunctor_h 1

#include <itkUnaryFunctorImageFilter.h>
#include <itkSymmetricSecondRankTensor.h>
#include <math.h>

#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_det.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

namespace itk
{

namespace Functor {

template< class TInputTensor, class TOutputScalar >
class TensorFractionalAnisotropy
{
  public:
  typedef typename TInputTensor::ValueType ValueType;
  typedef typename TInputTensor::EigenVectorsMatrixType MatrixType;
  typedef typename TInputTensor::EigenValuesArrayType EigenValuesType;

  TensorFractionalAnisotropy() {};
  ~TensorFractionalAnisotropy() {};

  inline TOutputScalar operator()( const TInputTensor & x )
  {
    MatrixType V;
    EigenValuesType e;
    TOutputScalar result;

    // Compute the tensor fractional anisotropy by, 
    //  1. computing the eigenanalysis of the tensor
    //  2. compute the parameters and write them out

// basser pierpaoli 1996
// Any tensor can be written as a sum of its isotropic part and 
//    its anisotropic part:
// D = <D>I + (D - <D>I) where I is the identity tensor, and
//                       <D> is the scalar called mean diffusivity.
// bulk mean diffusivity : trace(D)/3 = (e1+e2+e3)/3 = ebar = <D>
// The deviation tensor is :
//  {\em D} = (D - <D>I)
// 
// We can use the tensor scalar product which is
// \sqrt(D:D) = \sqrt{  \sum_i=1^3 \sum_j=1^3 D_{ij}^2 = e1^2 + e2^2 + e3^2 }
//  The tensor scalar product of the isotropic part of the tensor D is:
//  \sqrt( <D>I:<D>I ) = <D> sqrt(I:I) = <D> sqrt(3)
//  The tensor scalar product of the anisotropic part of the tensor D is:
//  \sqrt( {\em D}:{\em D} ) = \sqrt((e1 - ebar)^2 + (e2-ebar)^2 + (e3-ebar)^2)
//          = \sqrt( 3 * Var(\lamda) )

// RA = length of anisotropic part / length of isotropic part
//    = coefficient of variation of eigenvalues 
//    = sqrt( Variance(\lambda) ) / bulk mean diffusivity

//   FA has better noise properties = 
//      and is the proportion of the magnitude of D that can be ascribed to
//       anisotropic diffusion =
//      \sqrt(3 / 2) \sqrt( {\em D}:{\em D} ) /  \sqrt(D:D)
//      i.e. normalized length of anisotropic part / length of entire tensor

    x.ComputeEigenAnalysis(e, V);

    // bulk mean diffusivity : result = ( e[0] + e[1] + e[2] ) * 1.0 / 3.0;
    // bulk mean diffusivity : result = ( e[0] + e[1] + e[2] ) * 1.0 / 3.0;
    double ebar = (e[0] + e[1] + e[2])/3.0;
    double maganisotropicsquared =  (e[0] - ebar)*(e[0] - ebar) +
                             (e[1] - ebar)*(e[1] - ebar) +
                             (e[2] - ebar)*(e[2] - ebar);
    double magtensorsquared = e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
    // double magisotropic = ( e[0] + e[1] + e[2] )/ sqrt(3.0);

    // These do not depend on the ordering of the eigenvalues.
    // RA = sqrt( maganisotropicsquared ) / magisotropic 
    // FA = sqrt(3.0/2.0)*sqrt(maganisotropicsquared)/sqrt(magtensorsquared)

    // Here we return the fractional anisotropy
    result = sqrt(3.0/2.0)*sqrt(maganisotropicsquared)/sqrt(magtensorsquared);

    return result;
  }

};

template< class TInputTensor, class TOutputScalar >
class TensorAxialDiffusivity
{
  public:
  typedef typename TInputTensor::ValueType ValueType;
  typedef typename TInputTensor::EigenVectorsMatrixType MatrixType;
  typedef typename TInputTensor::EigenValuesArrayType EigenValuesType;

  TensorAxialDiffusivity() {};
  ~TensorAxialDiffusivity() {};

  inline TOutputScalar operator()( const TInputTensor & x )
  {
    MatrixType V;
    EigenValuesType e;
    TOutputScalar result;

    // Compute the tensor axial diffusivity by, 
    //  1. computing the eigenanalysis of the tensor
    //  2. compute the parameters and write them out

    // Axial diffusivity = e1 = \lambda_1
    // Radial diffusivity = 0.5*(e2+e3) = (\lambda_2 + \lambda_3)/2.0

    // The ordering of the eigenvalues generated is from smallest to largest.
    // Index 0 has the smallest eigenvalue, index 2 has the largest.
    x.ComputeEigenAnalysis(e, V);

    // Here we return the axial diffusivity - return the largest eigenvalue.
    result = e[2];

    return result;
  }

};

template< class TInputTensor, class TOutputScalar >
class TensorRadialDiffusivity
{
  public:
  typedef typename TInputTensor::ValueType ValueType;
  typedef typename TInputTensor::EigenVectorsMatrixType MatrixType;
  typedef typename TInputTensor::EigenValuesArrayType EigenValuesType;

  TensorRadialDiffusivity() {};
  ~TensorRadialDiffusivity() {};

  inline TOutputScalar operator()( const TInputTensor & x )
  {
    MatrixType V;
    EigenValuesType e;
    TOutputScalar result;

    // Compute the tensor fractional anisotropy by, 
    //  1. computing the eigenanalysis of the tensor
    //  2. compute the parameters and write them out

    // Axial diffusivity = e1 = \lambda_1
    // Radial diffusivity = 0.5*(e2+e3) = (\lambda_2 + \lambda_3)/2.0

    x.ComputeEigenAnalysis(e, V);

    // Here we return the radial diffusivity, the mean of the two smallest
    // eigenvalues.
    result = 0.5*(e[0]+e[1]);

    return result;
  }

};

template< class TInputTensor, class TOutputScalar >
class TensorMeanDiffusivity
{
  public:
  typedef typename TInputTensor::ValueType ValueType;
  typedef typename TInputTensor::EigenVectorsMatrixType MatrixType;
  typedef typename TInputTensor::EigenValuesArrayType EigenValuesType;

  TensorMeanDiffusivity() {};
  ~TensorMeanDiffusivity() {};

  inline TOutputScalar operator()( const TInputTensor & x )
  {
    MatrixType V;
    EigenValuesType e;
    TOutputScalar result;

    // Compute the tensor fractional anisotropy by, 
    //  1. computing the eigenanalysis of the tensor
    //  2. compute the parameters and write them out

    // Axial diffusivity = e1 = \lambda_1
    // Radial diffusivity = 0.5*(e2+e3) = (\lambda_2 + \lambda_3)/2.0
    // Bulk mean diffusivity = Trace(D)/3 

    x.ComputeEigenAnalysis(e, V);

    // Here we return the mean diffusivity
    result = (e[0]+e[1]+e[2])/3.0;

    return result;
  }

};

template< class TInputTensor, class TOutputScalar >
class TensorFrobeniusNorm
{
  public:
  typedef typename TInputTensor::ValueType ValueType;
  typedef typename TInputTensor::EigenVectorsMatrixType MatrixType;
  typedef typename TInputTensor::EigenValuesArrayType EigenValuesType;

  TensorFrobeniusNorm() {};
  ~TensorFrobeniusNorm() {};

  inline TOutputScalar operator()( const TInputTensor & x )
  {
    MatrixType V;
    EigenValuesType e;
    TOutputScalar result;

    // Compute the tensor fractional anisotropy by, 
    //  1. computing the eigenanalysis of the tensor
    //  2. compute the parameters and write them out

// 
// We can use the tensor scalar product which is
// \sqrt(D:D) = \sqrt{  \sum_i=1^3 \sum_j=1^3 D_{ij}^2 = e1^2 + e2^2 + e3^2 }
//  The tensor scalar product of the isotropic part of the tensor D is:
//  \sqrt( <D>I:<D>I ) = <D> sqrt(I:I) = <D> sqrt(3)
//  The tensor scalar product of the anisotropic part of the tensor D is:
//  \sqrt( {\em D}:{\em D} ) = \sqrt((e1 - ebar)^2 + (e2-ebar)^2 + (e3-ebar)^2)
//          = \sqrt( 3 * Var(\lamda) )

    x.ComputeEigenAnalysis(e, V);
    // std::cout << "x is " << x << std::endl;

    result = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);

    return result;
  }

};

template< class TInputTensor, class TOutputScalar >
class TensorLinearShapeMeasure
{
  public:
  typedef typename TInputTensor::ValueType ValueType;
  typedef typename TInputTensor::EigenVectorsMatrixType MatrixType;
  typedef typename TInputTensor::EigenValuesArrayType EigenValuesType;

  TensorLinearShapeMeasure() {};
  ~TensorLinearShapeMeasure() {};

  inline TOutputScalar operator()( const TInputTensor & x )
  {
    MatrixType V;
    EigenValuesType e;
    TOutputScalar result;

    x.ComputeEigenAnalysis(e, V);
    // std::cout << "x is " << x << std::endl;

    // We prefer this norm in order to maintain Cs + Cp + Cl = 1
    double norm = e[0] + e[1] + e[2];
    result = (e[2] - e[1])/norm; // Cl

    return result;
  }

};

template< class TInputTensor, class TOutputScalar >
class TensorPlanarShapeMeasure
{
  public:
  typedef typename TInputTensor::ValueType ValueType;
  typedef typename TInputTensor::EigenVectorsMatrixType MatrixType;
  typedef typename TInputTensor::EigenValuesArrayType EigenValuesType;

  TensorPlanarShapeMeasure() {};
  ~TensorPlanarShapeMeasure() {};

  inline TOutputScalar operator()( const TInputTensor & x )
  {
    MatrixType V;
    EigenValuesType e;
    TOutputScalar result;

    x.ComputeEigenAnalysis(e, V);
    // std::cout << "x is " << x << std::endl;

    // Note that \lambda_1 is the largest eigenvalue, but the ordering
    // returned by x.ComputeEigenAnalysis is smallest to largest.
    // Cp = 2*(\lambda_2 - \lambda_3) / sqrt(sum of squares of eigenvalues)

    // double norm = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
    // We prefer this norm in order to maintain Cs + Cp + Cl = 1
    double norm = e[0] + e[1] + e[2];
    result = 2.0*(e[1] - e[0])/norm; // Cp

    return result;
  }

};


template< class TInputTensor, class TOutputScalar >
class TensorSphericalShapeMeasure
{
  public:
  typedef typename TInputTensor::ValueType ValueType;
  typedef typename TInputTensor::EigenVectorsMatrixType MatrixType;
  typedef typename TInputTensor::EigenValuesArrayType EigenValuesType;

  TensorSphericalShapeMeasure() {};
  ~TensorSphericalShapeMeasure() {};

  inline TOutputScalar operator()( const TInputTensor & x )
  {
    MatrixType V;
    EigenValuesType e;
    TOutputScalar result;

    x.ComputeEigenAnalysis(e, V);
    // std::cout << "x is " << x << std::endl;

    // Note that \lambda_1 is the largest eigenvalue, but the ordering
    // returned by x.ComputeEigenAnalysis is smallest to largest.
    // Cs = 3*(\lambda_3) / sqrt(sum of squares of eigenvalues)

    // double norm = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
    // We prefer this norm in order to maintain Cs + Cp + Cl = 1
    double norm = e[0] + e[1] + e[2];
    result = 3.0*(e[0])/norm; // Cs

    return result;
  }

};


template< class TInputTensor, class TOutputScalar >
class TensorLargestEigenValue
{
  public:
  typedef typename TInputTensor::ValueType ValueType;
  typedef typename TInputTensor::EigenVectorsMatrixType MatrixType;
  typedef typename TInputTensor::EigenValuesArrayType EigenValuesType;

  TensorLargestEigenValue() {};
  ~TensorLargestEigenValue() {};

  inline TOutputScalar operator()( const TInputTensor & x )
  {
    MatrixType V;
    EigenValuesType e;
    TOutputScalar result;

    // Compute the tensor eigenvalues

    x.ComputeEigenAnalysis(e, V);
    // std::cout << "x is " << x << std::endl;

    result = e[2];

    return result;
  }

};

template< class TInputTensor, class TOutputScalar >
class TensorMiddleEigenValue
{
  public:
  typedef typename TInputTensor::ValueType ValueType;
  typedef typename TInputTensor::EigenVectorsMatrixType MatrixType;
  typedef typename TInputTensor::EigenValuesArrayType EigenValuesType;

  TensorMiddleEigenValue() {};
  ~TensorMiddleEigenValue() {};

  inline TOutputScalar operator()( const TInputTensor & x )
  {
    MatrixType V;
    EigenValuesType e;
    TOutputScalar result;

    // Compute the tensor eigenvalues

    x.ComputeEigenAnalysis(e, V);
    // std::cout << "x is " << x << std::endl;

    // The conventional ordering of eigenvalues is for \lambda_1 to be 
    // the largest and \lambda_3 the smallest.  The above EigenAnalysis returns
    // the results in the opposite order (smallest to largest).

    result = e[1];

    return result;
  }

};

template< class TInputTensor, class TOutputScalar >
class TensorSmallestEigenValue
{
  public:
  typedef typename TInputTensor::ValueType ValueType;
  typedef typename TInputTensor::EigenVectorsMatrixType MatrixType;
  typedef typename TInputTensor::EigenValuesArrayType EigenValuesType;

  TensorSmallestEigenValue() {};
  ~TensorSmallestEigenValue() {};

  inline TOutputScalar operator()( const TInputTensor & x )
  {
    MatrixType V;
    EigenValuesType e;
    TOutputScalar result;

    // Compute the tensor eigenvalues

    // This routine calls vnl_symmetric_eigensystem, which returns the 
    // eigenvectors as columns.
    x.ComputeEigenAnalysis(e, V);
    // std::cout << "x is " << x << std::endl;

    // The conventional ordering of eigenvalues is for \lambda_1 to be 
    // the largest and \lambda_3 the smallest.  The above EigenAnalysis returns
    // the results in the opposite order (smallest to largest).

    result = e[0];

    return result;
  }

};

template< class TInputTensor, class TOutputScalar >
class TensorComponent
{
  public:
  typedef typename TInputTensor::ValueType ValueType;
  typedef typename TInputTensor::EigenVectorsMatrixType MatrixType;
  typedef typename TInputTensor::EigenValuesArrayType EigenValuesType;

  unsigned int componentIndex;

  TensorComponent() { componentIndex = 0; };

  virtual ~TensorComponent() {}; 

   virtual void SetComponent( unsigned int newComponentIndex ) {
    componentIndex = newComponentIndex;
  }

  inline TOutputScalar operator()( const TInputTensor & x )
  {
    TOutputScalar result;
   
    result = x[componentIndex];

    return result;
  }

};

} // namespace Functor

template <class TInputTensorImage, class TOutputScalarImage >
class ITK_EXPORT TensorFractionalAnisotropyImageFilter :
public 
UnaryFunctorImageFilter<TInputTensorImage,TOutputScalarImage,
    Functor::TensorFractionalAnisotropy< 
        typename TInputTensorImage::PixelType, 
        typename TOutputScalarImage::PixelType> >
{
  public:
  /** Standard class typedefs. */
  typedef TensorFractionalAnisotropyImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputTensorImage, TOutputScalarImage,
    Functor::TensorFractionalAnisotropy<
        typename TInputTensorImage::PixelType,
        typename TOutputScalarImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  protected:
    TensorFractionalAnisotropyImageFilter() {}
    virtual ~TensorFractionalAnisotropyImageFilter() {}

  private:
    TensorFractionalAnisotropyImageFilter(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

};

template <class TInputTensorImage, class TOutputScalarImage >
class ITK_EXPORT TensorAxialDiffusivityImageFilter :
public 
UnaryFunctorImageFilter<TInputTensorImage,TOutputScalarImage,
    Functor::TensorAxialDiffusivity< 
        typename TInputTensorImage::PixelType, 
        typename TOutputScalarImage::PixelType> >
{
  public:
  /** Standard class typedefs. */
  typedef TensorAxialDiffusivityImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputTensorImage, TOutputScalarImage,
    Functor::TensorAxialDiffusivity<
        typename TInputTensorImage::PixelType,
        typename TOutputScalarImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  protected:
    TensorAxialDiffusivityImageFilter() {}
    virtual ~TensorAxialDiffusivityImageFilter() {}

  private:
    TensorAxialDiffusivityImageFilter(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

};


template <class TInputTensorImage, class TOutputScalarImage >
class ITK_EXPORT TensorRadialDiffusivityImageFilter :
public 
UnaryFunctorImageFilter<TInputTensorImage,TOutputScalarImage,
    Functor::TensorRadialDiffusivity< 
        typename TInputTensorImage::PixelType, 
        typename TOutputScalarImage::PixelType> >
{
  public:
  /** Standard class typedefs. */
  typedef TensorRadialDiffusivityImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputTensorImage, TOutputScalarImage,
    Functor::TensorRadialDiffusivity<
        typename TInputTensorImage::PixelType,
        typename TOutputScalarImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  protected:
    TensorRadialDiffusivityImageFilter() {}
    virtual ~TensorRadialDiffusivityImageFilter() {}

  private:
    TensorRadialDiffusivityImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

};

template <class TInputTensorImage, class TOutputScalarImage >
class ITK_EXPORT TensorMeanDiffusivityImageFilter :
public 
UnaryFunctorImageFilter<TInputTensorImage,TOutputScalarImage,
    Functor::TensorMeanDiffusivity< 
        typename TInputTensorImage::PixelType, 
        typename TOutputScalarImage::PixelType> >
{
  public:
  /** Standard class typedefs. */
  typedef TensorMeanDiffusivityImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputTensorImage, TOutputScalarImage,
    Functor::TensorMeanDiffusivity<
        typename TInputTensorImage::PixelType,
        typename TOutputScalarImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  protected:
    TensorMeanDiffusivityImageFilter() {}
    virtual ~TensorMeanDiffusivityImageFilter() {}

  private:
    TensorMeanDiffusivityImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

};

// Also called the Tensor Scalar Product
template <class TInputTensorImage, class TOutputScalarImage >
class ITK_EXPORT TensorFrobeniusNormImageFilter :
public 
UnaryFunctorImageFilter<TInputTensorImage,TOutputScalarImage,
    Functor::TensorFrobeniusNorm< 
        typename TInputTensorImage::PixelType, 
        typename TOutputScalarImage::PixelType> >
{
  public:
  /** Standard class typedefs. */
  typedef TensorFrobeniusNormImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputTensorImage, TOutputScalarImage,
    Functor::TensorFrobeniusNorm<
        typename TInputTensorImage::PixelType,
        typename TOutputScalarImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  protected:
    TensorFrobeniusNormImageFilter() {}
    virtual ~TensorFrobeniusNormImageFilter() {}

  private:
    TensorFrobeniusNormImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

};

template <class TInputTensorImage, class TOutputScalarImage >
class ITK_EXPORT TensorLargestEigenValueImageFilter :
public 
UnaryFunctorImageFilter<TInputTensorImage,TOutputScalarImage,
    Functor::TensorLargestEigenValue< 
        typename TInputTensorImage::PixelType, 
        typename TOutputScalarImage::PixelType> >
{
  public:
  /** Standard class typedefs. */
  typedef TensorLargestEigenValueImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputTensorImage, TOutputScalarImage,
    Functor::TensorLargestEigenValue<
        typename TInputTensorImage::PixelType,
        typename TOutputScalarImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  protected:
    TensorLargestEigenValueImageFilter() {}
    virtual ~TensorLargestEigenValueImageFilter() {}

  private:
    TensorLargestEigenValueImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

};

template <class TInputTensorImage, class TOutputScalarImage >
class ITK_EXPORT TensorMiddleEigenValueImageFilter :
public 
UnaryFunctorImageFilter<TInputTensorImage,TOutputScalarImage,
    Functor::TensorMiddleEigenValue< 
        typename TInputTensorImage::PixelType, 
        typename TOutputScalarImage::PixelType> >
{
  public:
  /** Standard class typedefs. */
  typedef TensorMiddleEigenValueImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputTensorImage, TOutputScalarImage,
    Functor::TensorMiddleEigenValue<
        typename TInputTensorImage::PixelType,
        typename TOutputScalarImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  protected:
    TensorMiddleEigenValueImageFilter() {}
    virtual ~TensorMiddleEigenValueImageFilter() {}

  private:
    TensorMiddleEigenValueImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

};


template <class TInputTensorImage, class TOutputScalarImage >
class ITK_EXPORT TensorSmallestEigenValueImageFilter :
public 
UnaryFunctorImageFilter<TInputTensorImage,TOutputScalarImage,
    Functor::TensorSmallestEigenValue< 
        typename TInputTensorImage::PixelType, 
        typename TOutputScalarImage::PixelType> >
{
  public:
  /** Standard class typedefs. */
  typedef TensorSmallestEigenValueImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputTensorImage, TOutputScalarImage,
    Functor::TensorSmallestEigenValue<
        typename TInputTensorImage::PixelType,
        typename TOutputScalarImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  protected:
    TensorSmallestEigenValueImageFilter() {}
    virtual ~TensorSmallestEigenValueImageFilter() {}

  private:
    TensorSmallestEigenValueImageFilter(const Self&);//purposely not implemented
    void operator=(const Self&); //purposely not implemented

};


template <class TInputTensorImage, class TOutputScalarImage >
class ITK_EXPORT TensorComponentImageFilter :
public 
UnaryFunctorImageFilter<TInputTensorImage,TOutputScalarImage,
    Functor::TensorComponent< 
        typename TInputTensorImage::PixelType, 
        typename TOutputScalarImage::PixelType> >
{
  public:
  /** Standard class typedefs. */
  typedef TensorComponentImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputTensorImage, TOutputScalarImage,
    Functor::TensorComponent<
        typename TInputTensorImage::PixelType,
        typename TOutputScalarImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  void SetComponent(unsigned int comp) 
  {
    this->GetFunctor().SetComponent( comp );
  }
  
  protected:
    TensorComponentImageFilter() {}
    virtual ~TensorComponentImageFilter() {}

  private:
    TensorComponentImageFilter(const Self&);//purposely not implemented
    void operator=(const Self&); //purposely not implemented

};

template <class TInputTensorImage, class TOutputScalarImage >
class ITK_EXPORT TensorLinearShapeImageFilter :
public 
UnaryFunctorImageFilter<TInputTensorImage,TOutputScalarImage,
    Functor::TensorLinearShapeMeasure< 
        typename TInputTensorImage::PixelType, 
        typename TOutputScalarImage::PixelType> >
{
  public:
  /** Standard class typedefs. */
  typedef TensorLinearShapeImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputTensorImage, TOutputScalarImage,
    Functor::TensorLinearShapeMeasure<
        typename TInputTensorImage::PixelType,
        typename TOutputScalarImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  protected:
    TensorLinearShapeImageFilter() {}
    virtual ~TensorLinearShapeImageFilter() {}

  private:
    TensorLinearShapeImageFilter(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

};


template <class TInputTensorImage, class TOutputScalarImage >
class ITK_EXPORT TensorPlanarShapeImageFilter :
public 
UnaryFunctorImageFilter<TInputTensorImage,TOutputScalarImage,
    Functor::TensorPlanarShapeMeasure< 
        typename TInputTensorImage::PixelType, 
        typename TOutputScalarImage::PixelType> >
{
  public:
  /** Standard class typedefs. */
  typedef TensorPlanarShapeImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputTensorImage, TOutputScalarImage,
    Functor::TensorPlanarShapeMeasure<
        typename TInputTensorImage::PixelType,
        typename TOutputScalarImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  protected:
    TensorPlanarShapeImageFilter() {}
    virtual ~TensorPlanarShapeImageFilter() {}

  private:
    TensorPlanarShapeImageFilter(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

};


template <class TInputTensorImage, class TOutputScalarImage >
class ITK_EXPORT TensorSphericalShapeImageFilter :
public 
UnaryFunctorImageFilter<TInputTensorImage,TOutputScalarImage,
    Functor::TensorSphericalShapeMeasure< 
        typename TInputTensorImage::PixelType, 
        typename TOutputScalarImage::PixelType> >
{
  public:
  /** Standard class typedefs. */
  typedef TensorSphericalShapeImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputTensorImage, TOutputScalarImage,
    Functor::TensorSphericalShapeMeasure<
        typename TInputTensorImage::PixelType,
        typename TOutputScalarImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  protected:
    TensorSphericalShapeImageFilter() {}
    virtual ~TensorSphericalShapeImageFilter() {}

  private:
    TensorSphericalShapeImageFilter(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#endif

