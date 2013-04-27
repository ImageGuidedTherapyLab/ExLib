
#ifndef _CRL_IMAGEALGEBRA_INCLUDED
#define _CRL_IMAGEALGEBRA_INCLUDED 1

#include <itkImage.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkCovariantVector.h>
#include <itkPoint.h>

// Used to interpolate the tangent field at continuous locations
#include <itkVectorLinearInterpolateImageFunction.h>

#ifdef WIN32
#include <missingFunctionsWindows.h>
#endif

class ITK_EXPORT crlImageLaplaceBase : public itk::LightObject
{
public:
  typedef crlImageLaplaceBase Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(crlImageLaplaceBase, LightObject);

  crlImageLaplaceBase() : 
   m_MaxNumberOfIterations(500),
   m_RelativeErrorTolerance(1e-03),
   m_InnerBorderVoltage(0.0),
   m_OuterBorderVoltage(10000.0)
   { };

  ~crlImageLaplaceBase() {};

  void SetLabelImage(std::string name) {
    m_LabelImageFileName = name;
  };
  void SetInnerBorderLabel(signed int label) {
    m_InnerBorderLabel = label;
  };
  void SetOuterBorderLabel(signed int label) {
    m_OuterBorderLabel = label;
  };
  void SetSolveRegionLabel(signed int label) {
    m_SolveRegionLabel = label;
  };
  void SetOutputFileName(std::string name) {
    m_OutputFileName = name;
  };

  void SetOutputTangentFieldFileName(std::string name) {
    m_OutputTangentFieldFileName = name;
  };

  void SetOutputDisplacementFieldFileName(std::string name) {
    m_OutputDisplacementFieldFileName = name;
  };

  void SetOutputThicknessFileName(std::string name) {
    m_OutputThicknessFileName = name;
  };

  virtual int Execute() = 0;

protected:
    double m_SolveRegionLabel;
    double m_OuterBorderLabel;
    double m_InnerBorderLabel;
    std::string m_LabelImageFileName;
    std::string m_OutputFileName;
    std::string m_OutputTangentFieldFileName;
    std::string m_OutputThicknessFileName;
    std::string m_OutputDisplacementFieldFileName;
    unsigned int m_MaxNumberOfIterations;
    double m_RelativeErrorTolerance;

    double m_InnerBorderVoltage;
    double m_OuterBorderVoltage;

private:
  crlImageLaplaceBase(const crlImageLaplaceBase&); // purposely not implemented
  void operator=(const crlImageLaplaceBase &); // purposely not implemented
};

template <unsigned int Dimension, class PixelType>
class ITK_EXPORT crlImageLaplace : public crlImageLaplaceBase 
{
public:
  /** Standard class typedefs. */
  typedef crlImageLaplace Self;
  typedef crlImageLaplaceBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

    typedef   itk::Image< PixelType, Dimension> ImageType;
    typedef   itk::ImageFileReader< ImageType >    ImageReaderType;
    typedef   itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
    typedef   itk::ImageFileWriter< ImageType >    ImageWriterType;
    typedef itk::CovariantVector< float, Dimension  >   VectorPixelType;  
    typedef itk::Image< VectorPixelType, Dimension  >   VectorImageType;
    typedef   itk::ImageRegionIteratorWithIndex<VectorImageType> VectorIteratorType;
    typedef   itk::ImageFileWriter< VectorImageType >    VectorImageWriterType;
    typedef   itk::VectorLinearInterpolateImageFunction<
                       VectorImageType, double >  VectorInterpolatorType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  crlImageLaplace()
  {
    m_Reader = ImageReaderType::New();
    m_Writer = ImageWriterType::New();
    m_ThicknessWriter = ImageWriterType::New();
    m_TangentFieldWriter = VectorImageWriterType::New();
    m_DisplacementFieldWriter = VectorImageWriterType::New();
    m_InnerBorderLabel = m_OuterBorderLabel = m_SolveRegionLabel = 0.0;
    m_Image = ImageType::New();
    m_ThicknessImage = ImageType::New();
    m_L0Image = ImageType::New();
    m_L1Image = ImageType::New();
    m_OutputTangentImage = VectorImageType::New();
    m_OutputDisplacementFieldImage = VectorImageType::New();
    m_VecInterpolator = VectorInterpolatorType::New();
    m_LabelImage = 0;
  }

  ~crlImageLaplace()
  {
    m_Image = 0;
    m_ThicknessImage = 0;
    m_L0Image = 0;
    m_L1Image = 0;
    m_OutputTangentImage = 0;
    m_OutputDisplacementFieldImage = 0;

    m_Writer = 0;
    m_ThicknessWriter = 0;
    m_Reader = 0;
    m_TangentFieldWriter = 0;
    m_DisplacementFieldWriter = 0;
    m_VecInterpolator = 0;
    m_LabelImage = 0;
  }

  int Execute();

  bool Derivative(double s,
           itk::FixedArray<double, Dimension> &y,
           itk::FixedArray<double, Dimension> &dydx);

  bool OneVoxelStep(itk::Point<double, Dimension> point,
                    itk::Point<double, Dimension> &newpoint)
  {
    typename VectorImageType::IndexType vindex;
    bool isinside = m_OutputTangentImage->TransformPhysicalPointToIndex(point,
                              vindex);
    if (!isinside) {
      return false; // No step can be taken.
    }

    typename VectorInterpolatorType::PixelType actualVector;

    // Get the tangent vector for the voxel closest to the input point.
    actualVector = m_OutputTangentImage->GetPixel( vindex );
    // Take a step
    for (unsigned int i = 0; i < Dimension; i++) {
      newpoint[i] = point[i] +  
                         actualVector[i]*m_OutputTangentImage->GetSpacing()[i];
    }
    return true;
  }

  // Initiate a stream line calculation at the specified start point,
  // follow it by numerical integration along the tangent field until
  // the path reaches out of the solution region, project the end point
  // back to the boundary of the solution region, and return the set of
  // steps from the start to the end representing the stream line.
  int ComputeStreamLine(itk::Point<double, Dimension> start,
    std::vector< itk::Point<double, Dimension> > &ysave);

  protected :

// In order to access the image of gradients, we define a pixel accessor.
class VectorPixelAccessor
{
public:

  void operator=( const VectorPixelAccessor & vpa )
    {
      m_Index = vpa.m_Index;
    }

  PixelType Get( const VectorPixelType& input ) const
    {
      return static_cast<PixelType>( input[ m_Index ] );
    }

  void SetIndex( unsigned int index )
    {
      m_Index = index;
    }
private:
  unsigned int m_Index;
};

    typename ImageReaderType::Pointer m_Reader;
    typename ImageWriterType::Pointer m_Writer;
    typename ImageWriterType::Pointer m_ThicknessWriter;
    typename VectorImageWriterType::Pointer m_TangentFieldWriter;
    typename VectorImageWriterType::Pointer m_DisplacementFieldWriter;

    typename ImageType::Pointer m_LabelImage;
    typename ImageType::Pointer m_OutImage;
    typename VectorImageType::Pointer m_OutputTangentImage;
    typename VectorImageType::Pointer m_OutputDisplacementFieldImage;
    typename VectorInterpolatorType::Pointer m_VecInterpolator;

private:
  crlImageLaplace(const Self&); // purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typename ImageType::Pointer m_Image;
  typename ImageType::Pointer m_ThicknessImage;
  typename ImageType::Pointer m_L0Image;
  typename ImageType::Pointer m_L1Image;

};

/* In order to compute stream lines across the Laplace solution, we
 * will use a numerical integration scheme.
 *  We will use a Carp-Cash Runge-Kutta numerical integration scheme with
 * adaptive step size control.
 *
 *  Given values of n variables y[0,...,n-1] and their derivatives 
 * dydx[0,...,n-1] known at x, use the fifth order Cash-Karp Runge-Kutta
 * method to advance the solution over an interval h and return the 
 * incremented variables as yout[0,...,n-1]. Also return an estimate of
 * the local truncation error in yout using the embedded fourth-order method.
 * The routine derivates(x,y,dydx) is a user-supplied routine that
 *   returns derivatives dydx at x.
 */
template <unsigned int Dimension, class PixelType>
class ITK_EXPORT crlRungeKutta
{

public:
  static const double a2 , a3 , a4 , a5 , a6 ;
  static const double b21 , b31 , b32  ,
         b41 , b42 ,     b43 ,
         b51 , b52 , b53 , b54 ,
         b61 , b62 , b63 ,
         b64 , b65 , 
         c1 , 
         c3 , 
         c4 , 
         c6 ;
  static const double dc5 ;
  static const double dc1 , // c1 - 
                      dc3 , // c3 - 
                      dc4 , // c4 - 
                      dc6 ; // c6 - 

  static const double SAFETY;
  static const double PGROW;
  static const double PSHRINK;
  static const double ERRCON; // (5/SAFETY) to power (1/PGROW).

protected:
  itk::FixedArray<double, Dimension> ak2;
  itk::FixedArray<double, Dimension> ak3;
  itk::FixedArray<double, Dimension> ak4;
  itk::FixedArray<double, Dimension> ak5;
  itk::FixedArray<double, Dimension> ak6;
  itk::FixedArray<double, Dimension> ytemp;

public:
  bool rkck(itk::FixedArray<double, Dimension> &y, 
            itk::FixedArray<double, Dimension> &dydx, 
            double x, double h, 
            itk::FixedArray<double, Dimension> &yout, 
            itk::FixedArray<double, Dimension> &yerr) 
  {

    for (unsigned int i = 0; i < Dimension; i++) {
      ytemp[i] = y[i] + b21 * h * dydx[i];
    }
    if (!m_crlImageLaplacePtr->Derivative(x+a2*h, ytemp, ak2)) return false;
    for (unsigned int i = 0; i < Dimension; i++) {
      ytemp[i] = y[i] + h*(b31*dydx[i] + b32*ak2[i]);
    }
    if (!m_crlImageLaplacePtr->Derivative(x+a3*h, ytemp, ak3)) return false;
    for (unsigned int i = 0; i < Dimension; i++) {
      ytemp[i] = y[i] + h*(b41*dydx[i] + b42*ak2[i] + b43*ak3[i]);
    }
    if (!m_crlImageLaplacePtr->Derivative(x+a4*h, ytemp, ak4)) return false;
    for (unsigned int i = 0; i < Dimension; i++) {
      ytemp[i] = y[i] + h*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] + b54*ak4[i]);
    }
    if (!m_crlImageLaplacePtr->Derivative(x+a5*h, ytemp, ak5)) return false;
    for (unsigned int i = 0; i < Dimension; i++) {
      ytemp[i] = y[i] + h*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] + b64*ak4[i]
                            + b65*ak5[i]);
    }
    if (!m_crlImageLaplacePtr->Derivative(x+a6*h, ytemp, ak6)) return false;
    // Accumulate the temporary results with the correct weights.
    for (unsigned int i = 0; i < Dimension; i++) {
      yout[i] = y[i] + h*(c1*dydx[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]);
    }
    // Estimate error as difference between fourth and fifth order methods.
    for (unsigned int i = 0; i < Dimension; i++) {
      yerr[i] = h*(dc1*dydx[i] + dc3*ak3[i] + dc4*ak4[i] 
                   + dc5*ak5[i] + dc6*ak6[i] );
    }
    
    return true; // Successful calculation.
  }

  // This routine uses a fifth-order Runge Kutta integration scheme,
  // with a Cash-Karp Runge Kutta step which monitors the step size
  // to ensure accuracy and to allow adaptation of the step size.
  //   On output, y and x are replaced with their newly estimated values,
  // hdid holds the step size actually taken, and hnext holds the estimate
  // of the next step size to take.
  int rkqs(itk::FixedArray<double, Dimension> &y, 
            itk::FixedArray<double, Dimension> &dydx,
       double & x, double htry, double eps, 
       itk::FixedArray<double, Dimension> &yscal, 
       double &hdid, double &hnext) {
    double errmax, h, htemp, xnew;
    itk::FixedArray<double, Dimension> ytemp;
    itk::FixedArray<double, Dimension> yerr;

    h = htry; // Set stepsize to the initial trial value.
    errmax = 2.0;
    while (errmax > 1.0) {
      // Take a step using rkck.
      if (!rkck(y, dydx, x, h, ytemp, yerr)) {
        errmax = 2.0; // Failure return indicates calculation outside of the
                      // image, or the solve region.
                      // Shrink the step size and try again.
      } else {
        errmax = 0.0; // Evaluate the accuracy
        for (unsigned int i = 0; i < Dimension; i++) {
          errmax = fmax(errmax, fabs(yerr[i]/yscal[i]));
        }
        errmax /= eps;  // Scale relative to required tolerance.
      }
      if (errmax > 1.0) { 
        // Truncation error too large - reduce stepsize.
        htemp = SAFETY * h * pow(errmax, PSHRINK);
        // Change by no more than a factor of ten.
        h = (h >= 0.0 ? fmax(htemp, 0.1*h) : fmin(htemp, 0.1*h) );
        xnew = x + h;
        if (xnew == x) {
          std::cerr << "Underflow in rkqs." << std::endl;
          return EXIT_FAILURE;
        }
      }
    }
    // Step succeeded when errmax <= 1.0. Compute size of next step.
    if (errmax > ERRCON) {
      hnext = SAFETY*h*pow(errmax, PGROW);
    } else {
      hnext = 5.0*h;  // No more than a factor of 5 increase.
    }
    hdid = h;
    x += h;
    for (unsigned int i = 0; i < Dimension; i++) {
      y[i] = ytemp[i];
    }
    return EXIT_SUCCESS;
  }

  public:
    crlImageLaplace<Dimension, PixelType> *m_crlImageLaplacePtr;
    crlRungeKutta() {
      m_crlImageLaplacePtr = 0;
    };
    crlRungeKutta(crlImageLaplace<Dimension, PixelType> *crlImageLaplacePtr) {
      m_crlImageLaplacePtr = crlImageLaplacePtr;
    };
    ~crlRungeKutta() { m_crlImageLaplacePtr = 0; };

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlImageLaplace.txx"
#endif

#endif

