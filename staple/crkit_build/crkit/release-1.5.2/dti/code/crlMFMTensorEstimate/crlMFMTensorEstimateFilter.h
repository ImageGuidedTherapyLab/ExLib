/**********************************************************************************************//**
 * \file	crlMFMTensorEstimate\crlMFMTensorEstimateFilter.h
 *
 * \brief	Declares the crlMFMTensorEstimateFilter class. 
*************************************************************************************************/

#ifndef CRL_TENSOR_ESTIMATEREG_FILTER_H
#define CRL_TENSOR_ESTIMATEREG_FILTER_H

#include "GlobalDefinitions.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/algo/vnl_svd.h"

#include "itkVectorContainer.h"
#include "itkVectorImage.h"
#include "itkBarrier.h"
#include "itkImageToImageFilter.h"
#include "itkDiffusionTensor3D.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"

#include "itkSingleValuedNonLinearOptimizer.h"


using namespace itk;

//------------------------------------------------------------
// Define some useful macros for the filter
// (modify the MFMTensorEstimateFilter::Parameters structure)
//------------------------------------------------------------
/** Set built-in type.  Creates member Set"name"() (e.g., SetVisibility()); */
#define mySetMacro(name,type) \
  virtual void Set##name (const type _arg) \
  { \
	itkDebugMacro("setting " #name " to " << _arg); \
	if (this->m_AlgoParameters.m_##name != _arg) \
    { \
	  this->m_AlgoParameters.m_##name = _arg; \
	  this->Modified(); \
    } \
  } 

/** Get built-in type.  Creates member Get"name"() (e.g., GetVisibility()); */
#define myGetMacro(name,type) \
	virtual type Get##name () \
{ \
	itkDebugMacro("returning " << #name " of " << this->m_AlgoParameters.m_##name ); \
	return this->m_AlgoParameters.m_##name; \
}

#define myGetSetMacros(name,type) \
  virtual void Set##name (const type _arg) \
  { \
	itkDebugMacro("setting " #name " to " << _arg); \
	if (this->m_AlgoParameters.m_##name != _arg) \
    { \
	  this->m_AlgoParameters.m_##name = _arg; \
	  this->Modified(); \
    } \
  } \
  virtual type Get##name () \
  { \
	itkDebugMacro("returning " << #name " of " << this->m_AlgoParameters.m_##name ); \
	return this->m_AlgoParameters.m_##name; \
  }


namespace crl {
	namespace DWI {

		class MFMTensorEstimateCostFunction;

		/**********************************************************************************************//**
		 * \class	MFMTensorEstimateFilter
		 *
		 * \brief	Estimate (multi-)tensor models using a variational approach in the log-euclidean
		 * 			framework.
		 * 			
		 * 			This filter takes two inputs:.
		 * 			-  Input 0 : The set of gradient images (GradientImagesType)  SET BY SetGradientImage
		 * 			-  Input 1 : The initial one-tensor solution. 
		 *
		 * \author	Benoit Scherrer
		 * \date	March 2010
		*************************************************************************************************/
		template< class TGradientImagePixelType, class TTensorPixelType=double >
		class ITK_EXPORT MFMTensorEstimateFilter :
			public ImageToImageFilter< Image< TGradientImagePixelType, 3 >,  Image< DiffusionTensor3D< TTensorPixelType >, 3 > >
		{

		public:

			//------------------------------------------------------------
			// Definite customized types
			//------------------------------------------------------------
			typedef MFMTensorEstimateFilter						Self;
			typedef itk::SmartPointer<Self>                          Pointer;
			typedef itk::SmartPointer<const Self>                    ConstPointer;
			typedef itk::ImageToImageFilter< Image< TGradientImagePixelType, 3>, Image< DiffusionTensor3D< TTensorPixelType >, 3 > > Superclass;
			typedef typename Superclass::OutputImageRegionType	OutputImageRegionType;

			// Pixel types
			typedef TGradientImagePixelType						GradientPixelType;
			typedef itk::DiffusionTensor3D< TTensorPixelType >	TensorPixelType;

			// Images
			typedef itk::Image< TensorPixelType, 3 >			TensorImageType;
			typedef TensorImageType								OutputImageType;
			typedef itk::Image<float, 4>						OutputFractionsImage;	/** Typedef defining the 4D image of fractions */
			typedef itk::Image< GradientPixelType, 3 >			GradientImageType;		/** Typedef defining one (of the many) gradient images.  */
			typedef itk::VectorImage< GradientPixelType, 3 >	GradientImagesType;		/** Typedef defining a set of gradient images */
			typedef itk::Image<InternalVoxelData *, 3>			TensorDataImageType;	/// Image of internal data
			typedef itk::VectorImage<float>						ResidualsImageType;
			typedef itk::Image< unsigned char, 3 >				MaskImageType;

			typedef vnl_vector_fixed< double,3 >				GradientDirectionType;  /** Holds each magnetic field gradient used to acquire one DWImage */

			// Iterators
			typedef typename	itk::ConstantBoundaryCondition<TensorDataImageType> ConstantBoundaryConditionType;
			typedef typename	itk::ShapedNeighborhoodIterator< TensorDataImageType, ConstantBoundaryConditionType > TensorDataImageIteratorType;

			// General optimizer
			typedef itk::SingleValuedNonLinearOptimizer	OptimizerType;

			//------------------------------------------------------------
			// ITK stuff
			//------------------------------------------------------------
			itkNewMacro(Self);						/** Method for creation through the object factory. */
			itkTypeMacro(MFMTensorEstimateFilter, ImageToImageFilter);	/** Runtime information support. */

public:
			void SetInputStudy( crl::DWIStudy<TGradientImagePixelType>& inputStudy );
			void SetInitTensorImage(  TensorImageType *image );
			void Set2TMaskImage( MaskImageType::Pointer maskImage );

			OutputFractionsImage*									GetFractionsImage() { return this->m_FractionsImage.GetPointer(); }
			const std::vector<typename OutputImageType::Pointer > &	GetInitTensors() const { return this->m_InitTensorImages; }
			typename ResidualsImageType::Pointer					GetResidualsImage() const{ return m_ResidualsImage; }

			myGetSetMacros( Threshold, double );							/** Threshold on the B=0 image data. The output tensor will be a null tensor for pixels in the reference image that have a value less than this threshold.*/
			myGetSetMacros( VerboseLevel, int );							/** Set/Get the verbose level */
			myGetSetMacros( OVerboseLevel, int );							/** Set/Get the optimization verbose level */
			myGetSetMacros( UseRicianNoiseModeling, bool );					/** Set/Get the use of Rician noise modeling */
			myGetSetMacros( RicianNoiseVariance, double );					/** Set/Get the estimtated Rician noise variance */
			myGetSetMacros( RobustEstimation, bool );						/** Set/Get the use of robust estimation */	
			myGetSetMacros( RobustParam, double );							/** Set/Get the robust parameter */
			myGetSetMacros( UseRegularization, bool );						/** Set/Get the use of regularization */
			myGetSetMacros( RegularizationWeight, double );					/** Set/Get the weight for the regularization */
			myGetSetMacros( RegularizationAnisotropicK, double );			/** Set/Get the weight for the regularization */
			myGetSetMacros( NumberOfTensors, int );							/** Set/Get the number of tensors to estimate */
			myGetSetMacros( MaxIter, unsigned int );						/** Set/Get the maximum number of iterations */
			myGetSetMacros( OutputInitTensors, bool );						/** Set/Get the output of init tensors */
			myGetSetMacros( OptimizationAlgo, crl::OptimizerUtils::OptimizerIdType );		/** Set/Get the estimation algorithm */
			myGetSetMacros( InitialOneTensorFromMedINRIA, bool );			/** Set/Get if the initial one tensor sol comes from MedINRIA (divides the eigen values by 1000) */
			myGetSetMacros( FirstPassWithLowBValuesOnly, bool );			/** Set/Get if we do a first pass with low b-values only */
			myGetSetMacros( OptimizationMultiRestart, bool );				/** Set/Get the multi restart mode */	
			myGetSetMacros( VideoMode, bool );								/** Set/Get the video mode */	
			myGetSetMacros( RandomInit, bool );								/** Set/Get the random initialization mode */	
			myGetSetMacros( EstimationStrategy, EstimationStrategyType );	/** Set/Get the random initialization mode */	
			myGetSetMacros( FastDiffusionFraction, double );				/** Set/Get The fraction of fast diffusion (slow diffusion fraction = 1-m_FastDiffusionFraction).*/
			myGetSetMacros( PruneThreshold, double );						/** Set/Get The minimum fraction (else prune the tensor) */
			myGetSetMacros( ExportResiduals, bool );						/** Set/Get if the filter should export the residuals */

		protected:
		//------------------------------------------------------------
		// Protected functions
		//------------------------------------------------------------
			MFMTensorEstimateFilter();
			~MFMTensorEstimateFilter() {};
			void PrintSelf(std::ostream& os, Indent indent) const;
	
			// Use the multithreaded filters
			virtual void GenerateOutputInformation();
			virtual void BeforeThreadedGenerateData();
			virtual void ThreadedGenerateData( const OutputImageRegionType &outputRegionForThread, int);
			virtual void AfterThreadedGenerateData();

			void OnePassWithSmallBValues( TensorDataImageIteratorType& tdit, OptimizerType *optimizer, MFMTensorEstimateCostFunction *costFunction );
			void SaveForVideo( int passID, TensorDataImageIteratorType& tdit, MFMTensorEstimateCostFunction *costFunction );

		private:
		//------------------------------------------------------------
		// Private data
		//------------------------------------------------------------
			GlobalAlgoParameters								m_AlgoParameters;				/** All the parameters for the algorithm */

			typename TensorImageType::Pointer					m_InitTensorImage;				/** Image of the initial one-tensor model */
			typename TensorDataImageType::Pointer				m_TensorDataImage;				/** Image of internal data for each voxel */
			typename OutputFractionsImage::Pointer				m_FractionsImage;				/** Images of the fractions */
			typename ResidualsImageType::Pointer				m_ResidualsImage;				/** Image of the residuals */

			std::vector<typename OutputImageType::Pointer >		m_InitTensorImages;				/** Set of images of the initial multi-tensor models */
			MaskImageType::Pointer								m_2TMaskImage;

			std::vector<unsigned int> 							m_BaselineInd; 					// contains the indicies of the baseline images
			std::vector<unsigned int> 							m_GradientInd; 					// contains the indicies of the gradient images

			itk::Barrier::Pointer								m_BarrierSync;					/** Barrier synchronization for threads */
		};

	}
}
#ifndef ITK_MANUAL_INSTANTIATION
#include "crlMFMTensorEstimateFilter.txx"
#endif

#endif
