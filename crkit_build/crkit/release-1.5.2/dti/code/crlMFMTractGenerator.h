
#ifndef _CRL_MFMTRACTGENERATOR_H_INCLUDED 
#define _CRL_MFMTRACTGENERATOR_H_INCLUDED 1

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkDiffusionTensor3D.h>
#include <itkRGBPixel.h>
#include <itkMatrix.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkSmartPointer.h>
#include <itkLightObject.h>
#include <itkImageRegionConstIterator.h>

#include <itkPolyLineParametricPath.h>
#include <itkVectorContainer.h>
#include <itkVector.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"

#ifdef WIN32
#include <missingFunctionsWindows.h>
#else
#include <sys/time.h>
#endif

#include "crlTensorLinearInterpolateImageFunction.h"

namespace crl {

class ITK_EXPORT MFMTractGeneratorBase : public itk::Object
{
  // Defines the public API for this class.
  // but has no templated arguments.  These are provided by classes that
  // inherit this API.
  //  This allows us to create a pointer to a class of this type, 
  // and then change it at run time to match the types we actually need to 
  // deal with.

  public:
	typedef enum {
		  INTERP_NONE = 0,
		  INTERP_EUCLIDEAN = 1,
		  INTERP_LOGEUCLIDEAN = 2
	} InterpolationModeValue;
	typedef unsigned int InterpolationModeType;

    typedef MFMTractGeneratorBase Self;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(MFMTractGeneratorBase, LightObject);

    MFMTractGeneratorBase() {
      m_ROIFileName = "";
      m_OutputFileName = "";
	  m_FractionFileName = "";
      m_MaxTractLength = 2000;
      m_MinimumFA = 0.15;
      m_FractionalAnisotropyMomentum = 0.0;
      pi = acos(-1.0);
      m_MaximumAngle = 80.0; // In degrees
      m_AngleThresholdMomentum = 0.5;
      m_DirectionMomentum = 0.1;
      m_TendFraction = 0.9;
      m_TractsPerVoxel = 1;
      m_StepsPerVoxel = 3;
      m_TractSimplificationMinAngle = 0;
	  m_ShowProgress = true;
	  m_ROILabels.clear();
      m_RandGen = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
      // select a seed based on the clock, to microsecond precision.
      unsigned int seed = 0;
      struct timeval preciseTimeStart;
      gettimeofday(&preciseTimeStart, NULL);
      seed = (unsigned int)(preciseTimeStart.tv_usec & 0xFFFFFFFF);
      m_RandGen->SetSeed(seed);
    }

    virtual ~MFMTractGeneratorBase() {
    }

    void SetInputFileNames( std::vector<std::string> inames) { 
      m_InputFileNames = inames; 
      // m_TensorImageVector.resize(    inames.size() );
      // m_LogTensorImageVector.resize( inames.size() );
    }

	void SetROIFileName( const std::string& s ) { m_ROIFileName = s; }
    void SetOutputFileName( const std::string& s ) { m_OutputFileName = s; }
    void SetFractionFileName( const std::string& s ) { m_FractionFileName = s; }

    void SetMaximumTractLength( unsigned int maxTractLength ) { 
      m_MaxTractLength = maxTractLength; 
    }

    //itkSetMacro( ROILabel, unsigned int);
    //itkGetMacro( ROILabel, unsigned int);
    
	void AddROILabel(unsigned int l)
	{
		m_ROILabels.push_back(l);
	}
	const std::vector<unsigned int>& GetROILabels() const { return m_ROILabels; }

    itkSetMacro( MinimumFA, double);
    itkGetMacro( MinimumFA, double);

    itkSetMacro( ShowProgress, bool);
    itkGetMacro( ShowProgress, bool);


    void SetMaximumAngle( double maxAngle ) { 
      // We store the cosine of the maximum angle,
      // and we expect the input is in degrees, so we convert here to radians.
      m_MaximumAngle = cos(maxAngle * pi / 180.0);
    }
    double GetMaximumAngle( ) const { 
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

    itkSetMacro( TensorInterpolationMode, InterpolationModeType );
    itkGetMacro( TensorInterpolationMode, InterpolationModeType );

    itkSetMacro( TractSimplificationMinAngle, double );
    itkGetMacro( TractSimplificationMinAngle, double );

    virtual int Execute() = 0;

  private:
    MFMTractGeneratorBase(const MFMTractGeneratorBase &); // purposely not implemented
    void operator=(const MFMTractGeneratorBase &); //purposely not implemented

  protected:
    std::vector<std::string> m_InputFileNames;
    std::string m_ROIFileName;
    std::string m_OutputFileName;
	std::string m_FractionFileName;

    unsigned int m_MaxTractLength;
    double m_MinimumFA;
    double m_FractionalAnisotropyMomentum;
    double m_MaximumAngle;
    double m_AngleThresholdMomentum;
	std::vector<unsigned int> m_ROILabels;	/// The list of labels to use for the seeding ROI
    double m_DirectionMomentum;
    double m_TendFraction;
    double pi;
    unsigned int m_TractsPerVoxel;
    unsigned int m_StepsPerVoxel;
    unsigned int m_TensorIndex; /// Record the tensor volume that is active.
	double m_TractSimplificationMinAngle;
	bool	m_ShowProgress;

	InterpolationModeType m_TensorInterpolationMode; /// Tensor interpolation mode
    itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer m_RandGen;


};

template <const unsigned int TImageDimension, class TComponentType >
class ITK_EXPORT MFMTractGenerator : 
                    public crl::MFMTractGeneratorBase
{
  public: 
  /** Standard class typedefs. */
  typedef MFMTractGenerator Self;
  typedef MFMTractGeneratorBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::DiffusionTensor3D<TComponentType> PixelType;
  typedef   itk::OrientedImage< PixelType, TImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> TensorImageReaderType;
  typedef   itk::ImageRegionConstIterator< ImageType> TensorImageIteratorType;

  typedef   itk::OrientedImage< unsigned short, TImageDimension> ROIImageType;
  typedef   itk::ImageFileReader<ROIImageType> ROIImageReaderType;

  typedef	itk::VectorImage<float, 3> FractionImageType;	/** Defines the 4D image of fractions */

  typedef itk::PolyLineParametricPath< TImageDimension > TractType;
  typedef itk::VectorContainer< unsigned int, typename TractType::Pointer >
    			TractContainerType;
  typedef itk::VectorContainer< unsigned int, vnl_vector_fixed<float, 3> > TractColorContainerType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  MFMTractGenerator() {
  };

  ~MFMTractGenerator() {};

  int Execute();

protected:
  int LoadMFM();
  int LoadROI();
  void InitializeAlgo( int nbThreads=-1 );

  void EstimateTract( const typename ImageType::PointType &point, 
                      typename TractType::Pointer &tract);

  typename TractType::Pointer SimplifyTract ( double angleMin, typename TractType::Pointer &inTract );

  void EstimatePrimaryEigenvector(const typename ImageType::PixelType &in, 
                               itk::Vector< TComponentType, 3> &evec) const;

  int SelectMatchingTensor(
        const typename ImageType::PixelType &in, 
        const typename ImageType::IndexType &index, 
        itk::Vector< TComponentType, 3> &evec,
        typename ImageType::PixelType &out) const;

  void InterpolateMatchingTensor(
        const typename ImageType::PixelType &in, 
        const typename TractType::ContinuousIndexType &cindex, 
        itk::Vector< TComponentType, 3> &evec,
        typename ImageType::PixelType &out) const;

  void GetFractionValues(const typename ImageType::IndexType &index, std::vector<float> &fractions ) const;

  unsigned int SelectTensorDependingOnFractions(const typename TractType::ContinuousIndexType& cindex) const;

  vtkPolyData* createPolyDataFromITK( typename TractContainerType::Pointer itkTracts, typename TractColorContainerType::Pointer itkTractColors) const;
  int WriteTracts( typename TractContainerType::Pointer tracts, typename TractColorContainerType::Pointer tractColors, const std::string& fileName) const;


  private:
    MFMTractGenerator(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

  protected:
    typename ROIImageType::Pointer m_ROIImage;
	typename FractionImageType::Pointer m_FractionImage;
	std::vector<typename ImageType::Pointer> m_TensorImageVector;		/// This is a std::vector to support multi-tensors
    std::vector<typename ImageType::Pointer> m_LogTensorImageVector;	/// Precomputed image of log-tensor for each tensors
    std::vector<typename crl::TensorLinearInterpolateImageFunction<
               ImageType, double >::Pointer> m_TensorInterpolatorVector;
	unsigned long int m_NbVoxelsInROI;
};

} // end namespace crl

#endif

