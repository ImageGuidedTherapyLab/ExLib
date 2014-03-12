#ifndef _ITK_CONTINUOUS_STAPLE_IMAGE_FILTER_H_
#define _ITK_CONTINUOUS_STAPLE_IMAGE_FILTER_H_

#include <iostream>
#include <fstream>
#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <vnl/vnl_matrix.h>
#include <vector>

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

namespace itk
{
  class ITK_EXPORT ContinuousSTAPLEImageFilter : 
    public ImageToImageFilter< itk::VectorImage <double, 3>, itk::VectorImage <double, 3> >
  {
    public:
      /** Standard class typedefs. */
      typedef ContinuousSTAPLEImageFilter Self;
      typedef itk::VectorImage <double, 3>  TInputImage;
      typedef itk::VectorImage <double, 3>  TOutputImage;
      typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
      typedef SmartPointer<Self> Pointer;
      typedef SmartPointer<const Self>  ConstPointer;
      
      /** Method for creation through the object factory. */
      itkNewMacro(Self);
      
      /** Run-time type information (and related methods) */
      itkTypeMacro(ContinuousSTAPLEImageFilter, ImageToImageFilter);

      typedef TOutputImage::PixelType OutputPixelType;
      typedef TInputImage::PixelType InputPixelType;
     
      /** Image typedef support */
      typedef TInputImage  InputImageType;
      typedef TOutputImage OutputImageType;
      typedef InputImageType::Pointer InputImagePointer;
      typedef OutputImageType::Pointer OutputImagePointer;

      /** Superclass typedefs. */
      typedef Superclass::InputImageRegionType InputImageRegionType;
      typedef Superclass::OutputImageRegionType OutputImageRegionType;

      typedef std::vector <double> BiasType;
      typedef vnl_matrix <double> CovMatrixType;

      VariableLengthVector<double> GetBias(unsigned int i) const
      {
	if (i > this->GetNumberOfInputs())
	  itkExceptionMacro(<< "Array reference out of bounds.");

	VariableLengthVector<double> resVal(m_NDim);
	for (unsigned int j = 0;j < m_NDim;++j)
	  resVal[j] = m_Bias[i][j];

	return resVal;
      }

      const CovMatrixType GetCovariance(unsigned int i) const
      {
	if (i > this->GetNumberOfInputs())
	  itkExceptionMacro(<< "Array reference out of bounds.");

	return m_Covariance[i];
      }

      const CovMatrixType GetReferenceCovariance() const
      {
	return m_ReferenceCovariance;
      }

      void InitializeExpertBias(double biasValue);
      void InitializeExpertCovariance(double covValue);

      itkSetMacro(Verbose, bool);

      /** Set/Get the maximum number of iterations after which the STAPLE algorithm
       *  will be considered to have converged. */
      itkSetMacro(MaximumIterations, unsigned int);
      itkGetMacro(MaximumIterations, unsigned int);

      itkSetMacro(NDim, unsigned int);
      itkGetMacro(NDim, unsigned int);
      
      /** Set/Get the threshold for which a change in the maximum of the difference between parameters
       * shall be so small as to trigger termination of the estimation procedure.
       */
      itkSetMacro(RelativeConvergenceThreshold, double);
      itkGetMacro(RelativeConvergenceThreshold, double);

      typedef itk::Image <unsigned char, 3> MaskImageType;
      typedef MaskImageType::RegionType MaskRegionType;
      typedef MaskImageType::Pointer MaskImagePointer;

      /** Set/Get the mask on which to compute STAPLE estimate. */
      itkSetMacro(ComputationMask, MaskImagePointer);
      itkGetMacro(ComputationMask, MaskImagePointer);
      
      itkSetMacro(OverrideComputationRegion, bool);

      /** Get the number of elapsed iterations of the iterative E-M algorithm. */
      itkGetMacro(ElapsedIterations, unsigned int);

      /** Compute the E-step of the algorithm 
       *  (i.e. the reference standard from the parameters and the data).
       */
      void EstimateReferenceStandard();
      void EstimateReferenceStandard(TOutputImage::RegionType &region);

      /** Compute the M-step of the algorithm 
       *  (i.e. the parameters of each expert from the reference standard and the data).
       */
      void EstimatePerformanceParameters();
      void EstimatePerformanceParameters(unsigned int minExp, unsigned int maxExp);
      
      bool endConditionReached();    
      void PrintPerformanceParameters(std::string &matlabFileName);

    protected:
      ContinuousSTAPLEImageFilter()
      {
	m_ElapsedIterations = 0;
	m_NDim = 6;
	m_MaximumIterations = NumericTraits<unsigned int>::max();
	m_RelativeConvergenceThreshold = 0;
	m_Bias.clear();
	m_Covariance.clear();
        m_InvCovariances.clear();
	m_OldBias.clear();
	m_OldCovariance.clear();
	m_ComputationMask = NULL;
	m_Verbose = true;
        
        m_lowerLimits.clear();
        m_upperLimits.clear();
        
        m_OverrideComputationRegion = false;
      }

      virtual ~ContinuousSTAPLEImageFilter() {}

      void GenerateData();
      int SplitRequestedRegion(int i, int num, OutputImageRegionType& splitRegion);

      void CreateComputationMaskFromInputs();
      void InitializeComputationRegionFromMask();
      void InitializeSplitRegionsFromMask();

      void PrintSelf(std::ostream&, itk::Indent) const;

      struct EMStepThreadStruct
      {
        Pointer Filter;
      };

      // Does the splitting and calls EstimateReferenceStandard on a sub region
      static ITK_THREAD_RETURN_TYPE ThreadEstimateRefStd( void *arg );

      // Does the splitting and calls EstimatePerformanceParameters on a sub sample of experts
      static ITK_THREAD_RETURN_TYPE ThreadEstimatePerfParams( void *arg );

    private:
      ContinuousSTAPLEImageFilter(const Self&); //purposely not implemented
      void operator=(const Self&); //purposely not implemented

      unsigned int m_ElapsedIterations;
      unsigned int m_MaximumIterations;
      unsigned int m_NDim;
      double       m_RelativeConvergenceThreshold;
      bool         m_Verbose, m_OverrideComputationRegion;

      //For overriding the split image region with something more intelligent taking into account the mask
      std::vector < unsigned int > m_lowerLimits, m_upperLimits;
      
      std::vector <BiasType> m_Bias;
      std::vector <CovMatrixType> m_Covariance;
      std::vector <CovMatrixType> m_InvCovariances;

      std::vector <BiasType> m_OldBias;
      std::vector <CovMatrixType> m_OldCovariance;

      MaskImagePointer m_ComputationMask;
      // Optimization of multithread code, compute only on region defined from mask... Uninitialized in constructor.
      MaskRegionType m_ComputationRegion;
      
      CovMatrixType m_ReferenceCovariance;
  };

} // end namespace itk

#endif
