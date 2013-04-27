
#ifndef _CRL_TRACTGENERATOR_H_INCLUDED 
#define _CRL_TRACTGENERATOR_H_INCLUDED 1

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkDiffusionTensor3D.h>
#include <itkRGBPixel.h>
#include <itkMatrix.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkSmartPointer.h>
#include <itkLightObject.h>

#include <itkPolyLineParametricPath.h>
#include <itkVectorContainer.h>
#include <itkVector.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

#ifdef WIN32
#include <missingFunctionsWindows.h>
#else
#include <sys/time.h>
#endif

#include "crlTensorLinearInterpolateImageFunction.h"

namespace crl {

class ITK_EXPORT TractGeneratorBase : public itk::Object
{
  // Defines the public API for this class.
  // but has no templated arguments.  These are provided by classes that
  // inherit this API.
  //  This allows us to create a pointer to a class of this type, 
  // and then change it at run time to match the types we actually need to 
  // deal with.

  public:
    typedef TractGeneratorBase Self;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(TractGeneratorBase, LightObject);

    TractGeneratorBase() {
      m_InputFileName = "";
      m_ROIFileName = "";
      m_OutputFileName = "";
      m_MaxTractLength = 600;
      m_MinimumFA = 0.15;
      m_FractionalAnisotropyMomentum = 0.0;
      pi = acos(-1.0);
      m_MaximumAngle = 80.0; // In degrees
      m_AngleThresholdMomentum = 0.5;
      m_DirectionMomentum = 0.1;
      m_TendFraction = 0.9;
      m_TractsPerVoxel = 1;
      m_StepsPerVoxel = 3;
      m_RandGen = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
      // select a seed based on the clock, to microsecond precision.
      unsigned int seed = 0;
      struct timeval preciseTimeStart;
      gettimeofday(&preciseTimeStart, NULL);
      seed = (unsigned int)(preciseTimeStart.tv_usec & 0xFFFFFFFF);
      m_RandGen->SetSeed(seed);
    }

    virtual ~TractGeneratorBase() {
    }

    void SetInputFileName( const char *s ) { m_InputFileName = s; }
    void SetROIFileName( const char *s ) { m_ROIFileName = s; }
    void SetOutputFileName( const char *s ) { m_OutputFileName = s; }

    void SetMaximumTractLength( unsigned int maxTractLength ) { 
      m_MaxTractLength = maxTractLength; 
    }

    itkSetMacro( ROILabel, unsigned int);
    itkGetMacro( ROILabel, unsigned int);

    itkSetMacro( MinimumFA, double);
    itkGetMacro( MinimumFA, double);

    void SetMaximumAngle( double maxAngle ) { 
      // We store the cosine of the maximum angle,
      // and we expect the input is in degrees, so we convert here to radians.
      m_MaximumAngle = cos(maxAngle * pi / 180.0);
    }
    double GetMaximumAngle( ) { 
      return ( acos(m_MaximumAngle)*180.0/pi );
    }

    itkSetMacro( DirectionMomentum, double);
    itkGetMacro( DirectionMomentum, double);

    itkSetMacro( AngleThresholdMomentum, double);
    itkGetMacro( AngleThresholdMomentum, double);

    itkSetMacro( FractionalAnisotropyMomentum, double);
    itkGetMacro( FractionalAnisotropyMomentum, double);

    itkSetMacro( TendFraction, double);
    itkGetMacro( TendFraction, double);

    itkSetMacro( TractsPerVoxel, unsigned int);
    itkGetMacro( TractsPerVoxel, unsigned int);

    itkSetMacro( StepsPerVoxel, unsigned int);
    itkGetMacro( StepsPerVoxel, unsigned int);

    virtual int Execute() = 0;

  private:
    TractGeneratorBase(const TractGeneratorBase &); // purposely not implemented
    void operator=(const TractGeneratorBase &); //purposely not implemented

  protected:
    std::string m_InputFileName;
    std::string m_ROIFileName;
    std::string m_OutputFileName;

    unsigned int m_MaxTractLength;
    double m_MinimumFA;
    double m_FractionalAnisotropyMomentum;
    double m_MaximumAngle;
    double m_AngleThresholdMomentum;
    unsigned int m_ROILabel;
    double m_DirectionMomentum;
    double m_TendFraction;
    double pi;
    unsigned int m_TractsPerVoxel;
    unsigned int m_StepsPerVoxel;
    itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer m_RandGen;


};

template <const unsigned int TImageDimension, class TComponentType >
class ITK_EXPORT TractGenerator : 
                    public crl::TractGeneratorBase
{
  public: 
  /** Standard class typedefs. */
  typedef TractGenerator Self;
  typedef TractGeneratorBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::DiffusionTensor3D<TComponentType> PixelType;
  typedef   itk::OrientedImage< PixelType, TImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> TensorImageReaderType;

  typedef   itk::OrientedImage< unsigned short, TImageDimension> ROIImageType;
  typedef   itk::ImageFileReader<ROIImageType> ROIImageReaderType;

  typedef itk::PolyLineParametricPath< TImageDimension > TractType;
  typedef itk::VectorContainer< unsigned int, typename TractType::Pointer >
    			TractContainerType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  TractGenerator() {
    m_TensorInterpolator = crl::TensorLinearInterpolateImageFunction< 
                                         ImageType, double >::New();
  };

  ~TractGenerator() {};

  void EstimateTract(typename ImageType::PointType &point, 
                      typename TractType::Pointer &tract);

  void EstimatePrimaryEigenvector(typename ImageType::PixelType &in, 
                               itk::Vector< TComponentType, 3> &evec);

  int Execute();

  private:
    TractGenerator(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

    typename ROIImageType::Pointer m_ROIImage;
    typename ImageType::Pointer m_TensorImage;
    typename ImageType::Pointer m_LogTensorImage;
    typename crl::TensorLinearInterpolateImageFunction< 
               ImageType, double >::Pointer m_TensorInterpolator;
};

} // end namespace crl

#endif

