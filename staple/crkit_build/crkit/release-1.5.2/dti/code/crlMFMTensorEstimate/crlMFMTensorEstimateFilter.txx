/*
 * Copyright (c) 2008-2009 Children's Hospital Boston.
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 * http://www.crl.med.harvard.edu
 * For more information contact: simon.warfield@childrens.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#ifndef CRL_TENSOR_ESTIMATEREG_FILTER_TXX
#define CRL_TENSOR_ESTIMATEREG_FILTER_TXX

#define DEBUG_SAVE_EVOLUTION_VIDEO

#include <math.h>
#include "crlMFMTensorEstimateFilter.h"
#include "crlMFMTensorEstimateCostFunction.h"
#include "crlTaskProgress.h"
#include "crlCommonMath.h"
#include "crlOptimizerUtils.h"

#include "itkArray.h"
#include "itkMacro.h"
#include "itkConjugateGradientOptimizer.h"
#include "itkLBFGSOptimizer.h"
#include "itkFRPROptimizer.h"
#include "itkPowellOptimizer.h"
#include "itkNLOPTOptimizers.h"
#include "itkImageFileWriter.h"

#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"


namespace itk {

/**********************************************************************************************//**
 * \class	NumericTraits<struct crl::DWI::SInternalVoxelData *>
 *
 * \brief	NumericTraits to be able to use itk::ConstantBoundaryCondition on  
 * 			crl::DWI::SInternalVoxelData data. 
 *
 * \author	Benoit Scherrer
 * \date	April 2010
*************************************************************************************************/
template <>
class NumericTraits<struct crl::DWI::SInternalVoxelData *> : public vcl_numeric_limits<int> {
public:
  typedef crl::DWI::SInternalVoxelData*           ValueType;
  typedef crl::DWI::SInternalVoxelData*           PrintType;
  typedef unsigned int  AbsType;
  typedef long          AccumulateType;
  typedef double        RealType;
  typedef RealType      ScalarRealType;
  typedef float         FloatType;

  //static crl::DWI::SInternalVoxelData* ITKCommon_EXPORT Zero;
  static crl::DWI::SInternalVoxelData* ITK_EXPORT Zero;

};

}

using namespace itk;

namespace crl {
	namespace DWI {

			/**********************************************************************************************//**
			 * \brief	Constructor 
			 *
			 * \author	Benoit Scherrer
			 * \date	March 2010
			 *
			 * \return	. 
			*************************************************************************************************/
		template< class TGradientImagePixelType, class TTensorPixelType >
			MFMTensorEstimateFilter< TGradientImagePixelType, TTensorPixelType >
			::MFMTensorEstimateFilter()
		{
			m_AlgoParameters.m_MaxIter = 2000;
			m_AlgoParameters.m_VerboseLevel = 0;
			m_AlgoParameters.m_OVerboseLevel = 0;
			m_AlgoParameters.m_NumberOfTensors = 0;
			m_AlgoParameters.m_RobustEstimation = true;
			m_AlgoParameters.m_RobustParam = 3.0;

			m_AlgoParameters.m_FastDiffusionFraction = 0;

			m_AlgoParameters.m_UseRegularization = false;
			m_AlgoParameters.m_RegularizationAnisotropicK = 0.001;
			m_AlgoParameters.m_RegularizationWeight = 0.0;

			m_AlgoParameters.m_OutputInitTensors = false;
			m_AlgoParameters.m_Threshold = NumericTraits< GradientPixelType >::min();

			m_AlgoParameters.m_OptimizationAlgo = crl::OptimizerUtils::OPTIMIZER_BOBYQA;
			m_AlgoParameters.m_OptimizationMultiRestart = true;

			m_AlgoParameters.m_InitialOneTensorFromMedINRIA = false;

			m_AlgoParameters.m_UseRicianNoiseModeling = false;
			m_AlgoParameters.m_RicianNoiseVariance = 1;

			m_AlgoParameters.m_FirstPassWithLowBValuesOnly = true;
			m_AlgoParameters.m_LowBValueMaxIndex = 0;

			m_AlgoParameters.m_EstimationStrategy = ESTIMATION_F_TENSORS;

			m_AlgoParameters.m_VideoMode = false;
			m_AlgoParameters.m_RandomInit = false;

			m_AlgoParameters.m_FastDiffusionFraction = 0;

			// At least 2 inputs is necessary: a vector image, and the init tensor image
			this->SetNumberOfRequiredInputs( 2 ); 

			// At least two output, and creates the second one
			this->SetNumberOfRequiredOutputs( 2 );
			this->SetNthOutput( 1, OutputImageType::New());
			}

			template< class TGradientImagePixelType, class TTensorPixelType >
			void MFMTensorEstimateFilter< TGradientImagePixelType, TTensorPixelType >
				::GenerateOutputInformation()
			{
				// Creates the additional outputs if required
				unsigned int prevNum = this->GetNumberOfOutputs();
				unsigned int nbOutputs = m_AlgoParameters.m_NumberOfTensors;
				this->SetNumberOfOutputs( nbOutputs );
				for ( unsigned int i = prevNum; i < nbOutputs; ++i )
				{
					this->SetNthOutput( i, this->MakeOutput( i ).GetPointer() );
				}

				// call the superclass' implementation of this method
				Superclass::GenerateOutputInformation();

				//// get pointers to the input and output
				//typename TOutputImage::RegionType outputLargestPossibleRegion;
				//outputLargestPossibleRegion.SetSize( m_Size );
				//outputLargestPossibleRegion.SetIndex( m_OutputStartIndex );

				//for ( unsigned int i=0; i<this->GetNumberOfOutputs() ; i++ )
				//{
				//	OutputImagePointer outputPtr = this->GetOutput(i);
				//	if ( outputPtr )
				//	{
				//		// Set the size of the output region
				//		outputPtr->SetLargestPossibleRegion( outputLargestPossibleRegion );

				//		// Set spacing and origin
				//		outputPtr->SetSpacing( m_OutputSpacing );
				//		outputPtr->SetOrigin( m_OutputOrigin );
				//		outputPtr->SetDirection( m_OutputDirection );

				//		//outputPtr->Allocate();
				//	}
				//}
			}

		template< class TGradientImagePixelType, class TTensorPixelType >
			void MFMTensorEstimateFilter< TGradientImagePixelType, TTensorPixelType >
			::BeforeThreadedGenerateData()
		{
			/*------------------------------------------
			Check if the gradients have been set, and check some
			other stuff
			-------------------------------------------*/
			if ( m_AlgoParameters.m_EstimationStrategy&ESTIMATION_RADIALBASELINE &&  m_AlgoParameters.m_EstimationStrategy&ESTIMATION_RADIALISO )
				itkExceptionMacro( << "ESTIMATION_RADIALBASELINE and ESTIMATION_RADIALISO cannot be set at the same time" );

			if ( (m_AlgoParameters.m_EstimationStrategy&ESTIMATION_RADIALISO) &&  (!(m_AlgoParameters.m_EstimationStrategy&ESTIMATION_F_TENSORS)) )
				itkExceptionMacro( << "ESTIMATION_F_TENSORS must be set to use ESTIMATION_RADIALISO" );

			if ( this->m_AlgoParameters.m_NumberOfTensors<1 )
				itkExceptionMacro( << "Invalid m_NumberOfTensors<1 specified!" );

			if ( this->GetInput(0)->GetLargestPossibleRegion().GetSize() != this->GetInput(1)->GetLargestPossibleRegion().GetSize() )
				itkExceptionMacro( << "Error. The DWI study image size does not match the one-tensor image size." );

			if ( this->m_2TMaskImage.GetPointer()!=NULL &&
				m_2TMaskImage->GetLargestPossibleRegion().GetSize() != this->GetInput(0)->GetLargestPossibleRegion().GetSize() )
				itkExceptionMacro( << "Error. The DWI study image size does not match the two-tensor mask image size." );

			/*------------------------------------------
			Allocate the internal data structure for each voxels
			-------------------------------------------*/
			m_TensorDataImage = TensorDataImageType::New();
			m_TensorDataImage->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
			m_TensorDataImage->Allocate();

			//----------------------------------------------
			// Allocate the residuals vector image if necessary
			//----------------------------------------------
			if ( m_AlgoParameters.m_ExportResiduals )
			{
				m_ResidualsImage = ResidualsImageType::New();
				m_ResidualsImage->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
				m_ResidualsImage->SetNumberOfComponentsPerPixel ( m_AlgoParameters.m_BValues.size() );
				m_ResidualsImage->Allocate();
				m_ResidualsImage->SetOrigin(this->GetInput(0)->GetOrigin());
				m_ResidualsImage->SetSpacing(this->GetInput(0)->GetSpacing());
				m_ResidualsImage->SetDirection(this->GetInput(0)->GetDirection());
			}

			/*------------------------------------------
			Get the size of the fractions image
			-------------------------------------------*/
			OutputFractionsImage::IndexType ofIndex;
			OutputFractionsImage::SizeType ofSize;
			OutputFractionsImage::PointType ofOrigin;
			OutputFractionsImage::SpacingType ofSpacing;
			OutputFractionsImage::DirectionType ofDir;
			ofIndex.Fill(0); ofSize.Fill(0); ofOrigin.Fill(0); ofSpacing.Fill(1); ofDir.Fill(0); 

			for ( int i=0 ; i<3 ; i++ )
			{
				ofIndex[i]=this->GetInput(0)->GetLargestPossibleRegion().GetIndex()[i];
				ofSize[i]=this->GetInput(0)->GetLargestPossibleRegion().GetSize()[i];
				ofSpacing[i]=this->GetInput(0)->GetSpacing()[i];
				ofOrigin[i] = this->GetInput(0)->GetOrigin()[i];
			}

			for ( int i=0; i<3 ; i++ ) {
				for ( int j=0; j<3 ; j++ )
				{
					ofDir(i,j) = this->GetInput(0)->GetDirection()(i,j);
				}
			}
			ofDir(3,3)=1.0;

			//----------------------------------------------
			// Compute the nb of fractions
			// (the cost function is not created yet)
			//----------------------------------------------
			int nbF;
			bool useIsotropic = ((m_AlgoParameters.m_EstimationStrategy&ESTIMATION_F_ISO)!=0);
			if ( m_AlgoParameters.m_NumberOfTensors>1 )
			{
				if ( useIsotropic ) nbF = m_AlgoParameters.m_NumberOfTensors+1;
				else nbF = (m_AlgoParameters.m_NumberOfTensors);
			}
			else if ( useIsotropic )
				nbF = 2;
			else
				nbF = 1;

			if ( m_AlgoParameters.m_EstimationStrategy&ESTIMATION_RADIALBASELINE )
				nbF+=m_AlgoParameters.m_NumberOfTensors;
			if ( m_AlgoParameters.m_EstimationStrategy&ESTIMATION_RADIALISO )
				nbF+=m_AlgoParameters.m_NumberOfTensors;

			ofSize[3] = nbF;

			//----------------------------------------------
			// Allocate the fractions image
			//----------------------------------------------
			OutputFractionsImage::RegionType region;
			region.SetIndex(ofIndex);
			region.SetSize(ofSize);
			m_FractionsImage = OutputFractionsImage::New();
			m_FractionsImage->SetRegions(region);
			m_FractionsImage->SetOrigin(ofOrigin);
			m_FractionsImage->SetSpacing(ofSpacing);
			m_FractionsImage->SetDirection(ofDir);
			m_FractionsImage->Allocate();


			/*------------------------------------------
			If we want to output the init tensors,
			allocate the required tensor images
			-------------------------------------------*/
			if ( m_AlgoParameters.m_OutputInitTensors )
			{
				m_InitTensorImages.resize(m_AlgoParameters.m_NumberOfTensors);
				for ( int i=0 ; i<m_AlgoParameters.m_NumberOfTensors ; i++ )
				{
					m_InitTensorImages[i] = OutputImageType::New();
					m_InitTensorImages[i]->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
					m_InitTensorImages[i]->SetOrigin(this->GetInput(0)->GetOrigin());
					m_InitTensorImages[i]->SetSpacing(this->GetInput(0)->GetSpacing());
					m_InitTensorImages[i]->SetDirection(this->GetInput(0)->GetDirection());
					m_InitTensorImages[i]->Allocate();
				}
			}

			//----------------------------------------------
			// Create the barrier to synchronize the 
			// initialization of all threads
			//----------------------------------------------
			m_BarrierSync = itk::Barrier::New();
			OutputImageRegionType  	fakeSplitRegion;
			int realNbThreads = SplitRequestedRegion (0, this->GetNumberOfThreads(), fakeSplitRegion );
			m_BarrierSync->Initialize(realNbThreads);

		}


		// 
		//
		// 
		template< class TGradientImagePixelType, class TTensorPixelType >
			void MFMTensorEstimateFilter< TGradientImagePixelType, TTensorPixelType >
			::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
			int procID ) 
		{

			//--------------------------------------------------------------
			//--------------------------------------------------------------
			// INITIALIZATION FOR THAT THREAD
			//--------------------------------------------------------------
			//--------------------------------------------------------------

			/*------------------------------------------
			Gets the output images pointers
			-------------------------------------------*/
			std::vector<typename OutputImageType::Pointer > outputImages;
			outputImages.resize(m_AlgoParameters.m_NumberOfTensors);
			for ( int i=0 ; i<m_AlgoParameters.m_NumberOfTensors ; i++ )
				outputImages[i] = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(i));

			/*------------------------------------------
			Creates the output images iterators for that thread
			-------------------------------------------*/
			std::vector< itk::ImageRegionIterator< OutputImageType > > oits;
			for ( int i=0 ; i<m_AlgoParameters.m_NumberOfTensors ; i++ )
			{
				oits.push_back( itk::ImageRegionIterator< OutputImageType >(outputImages[i], outputRegionForThread));
				oits[oits.size()-1].GoToBegin();
			}

			/*------------------------------------------
			Creates the iterator type for the tensor interal data
			(with the good boundary condition type)
			-------------------------------------------*/
			typename			TensorDataImageIteratorType::RadiusType radius;
			radius[0]=1;
			radius[1]=1;
			radius[2]=1;
			TensorDataImageIteratorType tdit(radius, m_TensorDataImage, outputRegionForThread);

			/*------------------------------------------
			Sets the shape of the shaped iterator
			-------------------------------------------*/
			TensorDataImageIteratorType::OffsetType o1 = {{-1,0,0}};
			TensorDataImageIteratorType::OffsetType o2 = {{1,0,0}};
			TensorDataImageIteratorType::OffsetType o3 = {{0,1,0}};
			TensorDataImageIteratorType::OffsetType o4 = {{0,-1,0}};
			TensorDataImageIteratorType::OffsetType o5 = {{0,0,1}};
			TensorDataImageIteratorType::OffsetType o6 = {{0,0,-1}};
			tdit.ActivateOffset( o1);
			tdit.ActivateOffset( o2);
			tdit.ActivateOffset( o3);
			tdit.ActivateOffset( o4);
			tdit.ActivateOffset( o5);
			tdit.ActivateOffset( o6);
			tdit.GoToBegin();

			/*------------------------------------------
			Gets the input images
			-------------------------------------------*/
			typename GradientImagesType::Pointer gradientImagePointer = NULL;
			typename TensorImageType::Pointer initTensorImagePointer = NULL;

			gradientImagePointer = static_cast< GradientImagesType * >( this->ProcessObject::GetInput(0) );
			initTensorImagePointer = static_cast< TensorImageType * >( this->ProcessObject::GetInput(1) );

			/*------------------------------------------
			Creates the input iterators
			-------------------------------------------*/
			typedef ImageRegionConstIterator< TensorImageType >		InitTensorIteratorType;
			typedef ImageRegionConstIterator< GradientImagesType >	GradientIteratorType;
			typedef typename GradientImagesType::PixelType			GradientVectorType;

			GradientIteratorType git(gradientImagePointer, outputRegionForThread );
			InitTensorIteratorType itit(initTensorImagePointer, outputRegionForThread);
			git.GoToBegin();
			itit.GoToBegin();

			//--------------------------------------------------------
			// Creates the iterator for the two-tensor mask
			//--------------------------------------------------------
			typedef ImageRegionConstIterator< MaskImageType >		MaskIteratorType;
			MaskIteratorType itMask;
			bool use2TMask = (m_2TMaskImage.GetPointer()!=NULL);
			if ( use2TMask ) itMask = MaskIteratorType(m_2TMaskImage, outputRegionForThread);

			//--------------------------------------------------------------
			//--------------------------------------------------------------
			// STEP 1. Creates all precompute data for each voxel
			//--------------------------------------------------------------
			//--------------------------------------------------------------
			if ( procID==0 ) cout<<"Initialize the algorithm..."<<endl;
			while ( !tdit.IsAtEnd() )
			{
				InternalVoxelData *tensorData = new InternalVoxelData;

				//----------------------------------------
				// Compute S0 only one time
				//----------------------------------------
				GradientVectorType Si = git.Get();
				double S0 = 0;
				for(unsigned int i = 0; i < m_BaselineInd.size(); ++i) S0 += Si[m_BaselineInd[i]];
				S0 /= ((double)m_BaselineInd.size());
				tensorData->S0 = S0;

				tensorData->RobustEstimation_C = 0;
				tensorData->HasConverged = false;

				//----------------------------------------
				// Check that the init tensor is not null
				//----------------------------------------
				TensorPixelType initTensor = itit.Value();
				TTensorPixelType *p = initTensor.GetDataPointer();
				bool ok=false;
				for ( int i=0 ; i<6 && (!ok); i++ )
				{
					if ( *p!=0 ) ok=true;
					p++;
				}

				//----------------------------------------
				// Check if valid tensor
				//----------------------------------------
				if ( ok && S0 >=m_AlgoParameters.m_Threshold ) 
				{
					// Copy the Si in the structure
					tensorData->Si.clear();
					for ( unsigned int i=0 ; i<m_GradientInd.size() ; i++ )
						tensorData->Si.push_back(Si[m_GradientInd[i]]);
				}
				//----------------------------------------
				// Else if invalid tensor don t copy
				//----------------------------------------
				else 
				{
					delete tensorData;
					tensorData = NULL;
				}

				// Set the tensor data
				tdit.SetCenterPixel(tensorData);
				++tdit;
				++git;
				++itit;
			}
			tdit.GoToBegin();
			git.GoToBegin();
			itit.GoToBegin();


			//--------------------------------------------------------------
			//--------------------------------------------------------------
			// STEP 2. Creates the costfunction and the optimizer
			//--------------------------------------------------------------
			//--------------------------------------------------------------

			//--------------------------------------------------------------
			// Create the cost function and set parameters
			//--------------------------------------------------------------
			MFMTensorEstimateCostFunction::Pointer costFunction = MFMTensorEstimateCostFunction::New();
			costFunction->SetProcID(procID);
			costFunction->SetAlgoParameters(&m_AlgoParameters);
			costFunction->SetUseRobustEstimation(false);

			//--------------------------------------------------------------
			// Create the optimizer
			//--------------------------------------------------------------
			unsigned int MAX_ITER_PER_PASS = m_AlgoParameters.m_MaxIter;			//m_MaxIter/4; //m_MaxIter;
			if ( m_AlgoParameters.m_UseRegularization || m_AlgoParameters.m_RobustEstimation ) 
				MAX_ITER_PER_PASS=100;
			else if ( m_AlgoParameters.m_OptimizationMultiRestart )
				MAX_ITER_PER_PASS = 200;
			if ( m_AlgoParameters.m_VideoMode ) MAX_ITER_PER_PASS=20;

			OptimizerUtils::GenericOptimizerType::Pointer optimizer = OptimizerUtils::CreateOptimizer(m_AlgoParameters.m_OptimizationAlgo, MAX_ITER_PER_PASS);
			OptimizerUtils::SetOptimizerParameters(optimizer, 
					1e-4,	// double stepLength 
					1e-7,	// double stepTol 
					1e-7,	// double valTol
					1e-5,	// double gradTol 
					MAX_ITER_PER_PASS,
					false );

			//--------------------------------------------------------------
			//--------------------------------------------------------------
			// STEP 3. FIRST INIT ALL THE TENSORS
			//--------------------------------------------------------------
			//--------------------------------------------------------------
			if ( procID==0 ) {
				if ( m_AlgoParameters.m_NumberOfTensors==2 ) cout<<"Initialize all two-tensors..."<<endl;
				else cout<<"Initialize all tensors..."<<endl;
			}


			//--------------------------------------------------------
			// SYNCHRONIZATION OF ALL THE THREADS:
			// WAIT THAT EACH OF THEM TO HAVE CREATED THE InternalVoxelData STRUCTURE TO DEAL
			// WITH NEIGHBORS
			//--------------------------------------------------------
			m_BarrierSync->Wait();

			//--------------------------------------------------------
			// 1. Initialize each tensor 
			// 2. Output the init tensors in a map if necessary
			// 3. Initialize the InternalVoxelData::Neighbors if regularization
			//--------------------------------------------------------

			// For the initialization be sure to use all the b-values
			costFunction->SetCurrentMaxBValueIndexToUse(m_AlgoParameters.m_BValues.size());
			std::vector<TensorPixelType > initTensors;

			while( !itit.IsAtEnd() )
			{				
				InternalVoxelData* tensorData = tdit.GetCenterPixel();
				TensorPixelType initTensor = itit.Value();

				// If tensor OK
				if ( tensorData )
				{
					bool forceSingleTensor=false;
					if ( use2TMask ) forceSingleTensor = (itMask.Value()==0);

					if ( forceSingleTensor ) tensorData->nbTensorAtVoxel = 1;
					else tensorData->nbTensorAtVoxel = m_AlgoParameters.m_NumberOfTensors;
					
					// Compute the initial tensors data (in tensorData)
					costFunction->SetInternalVoxelData(tensorData);
					costFunction->InitAlgo(initTensor, forceSingleTensor );
					costFunction->ComputeFinalTensors(initTensors);
				}
				// Else just output null tensor
				else
				{
					initTensors.resize( m_AlgoParameters.m_NumberOfTensors);
					for ( int i=0 ; i< m_AlgoParameters.m_NumberOfTensors ; i++ )
						initTensors[i].Fill(0);

					// Even if the tensor is invalid, just allocate fractions with 0
					/*int nbF = costFunction->GetNumberOfFractionsToEstimate();
					if ( nbF!=0 )
					{
						if ( m_AlgoParameters.m_OptimizationAlgo==crl::OptimizerUtils::OPTIMIZER_BOBYQA ) nbF++;
						tensorData->Fractions.set_size( nbF );
					}
					else
						tensorData->Fractions.set_size(1);
					tensorData->Fractions.fill(0);

					tensorData->normalizedFractions.set_size(tensorData->Fractions.size());
					tensorData->normalizedFractions.fill(0);*/
				}

				// If needed, output init tensors
				if (  m_AlgoParameters.m_OutputInitTensors )
				{
					typename InitTensorIteratorType::IndexType idx = itit.GetIndex();
					for ( int i=0 ; i< m_AlgoParameters.m_NumberOfTensors ; i++ )
						m_InitTensorImages[i]->SetPixel(idx, initTensors[i]);
				}

				// If regularization, init the neighbors structure
				if (  m_AlgoParameters.m_UseRegularization && tensorData )
				{
					const TensorDataImageIteratorType::IndexListType &list = tdit.GetActiveIndexList ();
					TensorDataImageIteratorType::IndexListType::const_iterator lit(list.begin());
					tensorData->Neighbors.resize(6);
					for ( int i=0; i<6 ; i++ ) tensorData->Neighbors[i]=NULL;

					 for ( unsigned int i=0; i<tdit.GetActiveIndexListSize (); i++ )
					 {
						 TensorDataImageIteratorType::OffsetType o = tdit.GetOffset(*lit);

						 if ( o[0] == -1 )  tensorData->Neighbors[0]=tdit.GetPixel(*lit);
						 else if ( o[0] == 1 )  tensorData->Neighbors[1]=tdit.GetPixel(*lit);
						 else if ( o[1] == -1 )  tensorData->Neighbors[2]=tdit.GetPixel(*lit);
						 else if ( o[1] == 1 )  tensorData->Neighbors[3]=tdit.GetPixel(*lit);
						 else if ( o[2] == -1 )  tensorData->Neighbors[4]=tdit.GetPixel(*lit);
						 else if ( o[2] == 1 )  tensorData->Neighbors[5]=tdit.GetPixel(*lit);
						 else throw itk::ExceptionObject(__FILE__, __LINE__, "Invalid offset","crlMFMTensorEstimateFilter");

						 // tensorData->Neighbors[1+o[2]][1+o[1]][1+o[0]] = tdit.GetPixel(*lit);
						 lit++;
					 }
				}

				++itit;
				++tdit;
				if ( use2TMask ) ++itMask;
			}
			tdit.GoToBegin();
			itit.GoToBegin();

			//--------------------------------------------------------------
			// SYNCHRONIZATION OF ALL THE THREADS:
			// WAIT THAT EACH OF THEM HAS FINISHED ITS INITIALIZATION
			//--------------------------------------------------------------
			m_BarrierSync->Wait();
			if ( procID==0 ) cout<<">>RUN<<"<<endl;
			m_BarrierSync->Wait();		// Rewait (that the thread 0 has performed its cout

			//--------------------------------------------------------------
			// ONE INITIALIZATION PASS WITH ONLY SMALL B-VALUES ?
			//--------------------------------------------------------------
			if (  m_AlgoParameters.m_FirstPassWithLowBValuesOnly )
				OnePassWithSmallBValues( tdit, optimizer, costFunction );

			// NOW be sure to use all the b-values
			costFunction->SetCurrentMaxBValueIndexToUse( m_AlgoParameters.m_BValues.size() );

//#define OUTPUT_EVAL_IN_FILE 1

#ifdef OUTPUT_EVAL_IN_FILE
			std::fstream debugfile;
			char szBuffer[256];
			sprintf(szBuffer,"EVAL_p%d.txt",procID);
			debugfile.open ( szBuffer, std::fstream::out );
			if  ( !debugfile.is_open() ) 
			{
				cout<<"Cannot open the file"<<szBuffer<<endl;
				return;
			}
#endif



			//--------------------------------------------------------------
			//--------------------------------------------------------------
			// STEP 4. MAIN LOOP
			//--------------------------------------------------------------
			//--------------------------------------------------------------

			int NbMaxPasses;
			//if ( m_UseRegularization || m_RobustEstimation ) 
				NbMaxPasses = (int)m_AlgoParameters.m_MaxIter/MAX_ITER_PER_PASS;
			//else if (m_RobustEstimation) 
			//	NbMaxPasses=2;
			//else
			//	NbMaxPasses=1;
				
			if ( NbMaxPasses<1 ) NbMaxPasses=1;	//to be sure...

			// For the first pass, don't use reg or robust
			costFunction->SetUseRegularization(false);
			costFunction->SetUseRobustEstimation(false);

			//-----------------------------------------------------
			// Prepare the NLOPT pointer in case we use bounds
			//-----------------------------------------------------
			itk::NLOPTOptimizers *nlopt = NULL;
			if ( this->m_AlgoParameters.m_OptimizationAlgo==crl::OptimizerUtils::OPTIMIZER_BOBYQA )
			{
				nlopt = dynamic_cast<itk::NLOPTOptimizers *>(optimizer.GetPointer());
			}


			TaskProgress taskProgress(procID, this->GetNumberOfThreads(), NbMaxPasses * outputRegionForThread.GetSize()[0]*outputRegionForThread.GetSize()[1]*outputRegionForThread.GetSize()[2]);

			//--------------------------------------------------------------
			// MAIN LOOP : perform all the passess
			//--------------------------------------------------------------
			for ( int icmPass = 0; icmPass<NbMaxPasses ; icmPass++ )
			{
				//------------------------------------------
				// tmp tmp...
				//------------------------------------------
#ifdef DEBUG_SAVE_EVOLUTION_VIDEO
				if ( m_AlgoParameters.m_VideoMode ) 
				{
					SaveForVideo(icmPass, tdit, costFunction );

					OptimizerUtils::SetOptimizerMaxNbIter(optimizer, 20+2*(icmPass-NbMaxPasses/2) );
				}
#endif

				//------------------------------------------
				// If use regularization, only actually use it after
				// the first pass, and gradually (linearly) 
				// increase the value of the regularization weight
				// until its max value m_AlgoParameters.m_RegularizationWeight
				//------------------------------------------
				if (  m_AlgoParameters.m_UseRegularization )
				{
					// Don t use regularization for the first pass
					if ( icmPass<=0 ) costFunction->SetUseRegularization(false);
					else 
					{
						costFunction->SetUseRegularization(true);
						costFunction->SetCurrentRegularizationWeight(  m_AlgoParameters.m_RegularizationWeight * my_MIN(1.0, icmPass/((double)NbMaxPasses)) ); //*((double)icmPass+1)/((double)NbMaxPasses));
					}
				}

				//--------------------------------------------------
				// Iterate on all pixels of the region
				//---------------------------------------------------
				tdit.GoToBegin();
				while( !tdit.IsAtEnd() )
				{
					taskProgress.Update();

					//----------------------------------------
					// Get a pointer to the tensor data
					//----------------------------------------
					InternalVoxelData* tensorData = tdit.GetCenterPixel();
					costFunction->SetInternalVoxelData(tensorData);

					//----------------------------------------
					// Launch the optimization for that tensor
					//----------------------------------------
					if ( tensorData && ((!tensorData->HasConverged) || icmPass==NbMaxPasses-1) ) //just allow the last pass
					{
						//---------------------------------------------------------------
						// Setup the scales for the optimizer
						// (the number of tensors can be variable from one voxel to an other)
						//---------------------------------------------------------------
						optimizer->SetCostFunction(costFunction);
						optimizer->Modified();

						OptimizerType::ScalesType					scales;
						scales.set_size(costFunction->GetNumberOfParameters());
						for ( unsigned int i=0; i<scales.size() ; i++ ) 
							scales[i]=1; //   1/0.01;

						if ( costFunction->GetNumberOfFractionsToEstimate()!=0 ) 
							for ( int i=0; i<costFunction->GetNumberOfFractionsToEstimate() ; i++ )
								scales[scales.size()-i-1]=0.1;   //  1/0.1

						optimizer->SetScales(scales);

						try {
							//-----------------------------------------------------
							// If robust estimation, we actually use the robust estimation
							// only after the third pass
							//-----------------------------------------------------
							if (  m_AlgoParameters.m_RobustEstimation )
							{
								if ( icmPass<2 ) costFunction->SetUseRobustEstimation(false);
								else  if ( icmPass==2) 
								{
									costFunction->ComputeC_ForRobustEstimation();
									costFunction->SetUseRobustEstimation(true);
								}
								else if (icmPass>2)
								{
									//costFunction->ComputeC_ForRobustEstimation();
									costFunction->SetUseRobustEstimation(true);
								}
							}

							//-----------------------------------------------------
							// Get initial parameters from the value of the L's
							//-----------------------------------------------------
							MFMTensorEstimateCostFunction::ParametersType initP = costFunction->GetOptimizerParametersFromPrecomputeData();
							optimizer->SetInitialPosition( initP );

							//-----------------------------------------------------
							// Setup the bounds if necessary
							//-----------------------------------------------------
							if ( nlopt!=NULL )
							{
								nlopt->SetLowerBoundParameters(tensorData->lowerBound);
								nlopt->SetUpperBoundParameters(tensorData->upperBound);
							}

							//-----------------------------------------------------
							// Start the optimization !
							//-----------------------------------------------------
							optimizer->StartOptimization();

							//-----------------------------------------------------
							// Get the output optimizer position 
							//-----------------------------------------------------
							MFMTensorEstimateCostFunction::ParametersType endP = optimizer->GetCurrentPosition();
							costFunction->EndAlgoIteration(endP);
						}
						//-----------------------------------------------------
						// If error....
						//-----------------------------------------------------
						catch (...)
						{
							cout<<">>> Exception f="<<costFunction->GetValue(optimizer->GetCurrentPosition())<<endl;
						}

						//-----------------------------------------------------
						// VERBOSE: show f/dF
						//-----------------------------------------------------
						if ( m_AlgoParameters.m_VerboseLevel>=3 )
						{
							MFMTensorEstimateCostFunction::ParametersType df;
							costFunction->GetDerivative(optimizer->GetCurrentPosition(), df);

							double norm = 0 ;
							for ( unsigned int i=0; i<df.size() ; i++ ) norm += df[i]*df[i];
							norm = sqrt(norm);

							cout<<"  End pos f="<<costFunction->GetValue(optimizer->GetCurrentPosition());
							cout << " |dF|=" << norm<<"  ( dF=" ;
							for ( unsigned int i=0; i<df.size() ; i++ ) cout<< df[i]<< " ";
							cout<<")"<<endl;
						}

						//-----------------------------------------------------
						// TEST IF NAN
						//-----------------------------------------------------
						if ( isnan(costFunction->GetValue(optimizer->GetCurrentPosition()) ) )
						{
							cout<<"NAN REACHED"<<endl<<"Initial one tensor solution: "<<endl;
							//double *p = initTensor.GetDataPointer();
							//cout<<"  "<<*(p)<<"\t"<<*(p+1)<<"\t"<<*(p+2)<<endl;
							//cout<<"  "<<*(p+1)<<"\t"<<*(p+3)<<"\t"<<*(p+4)<<endl;
							//cout<<"  "<<*(p+2)<<"\t"<<*(p+4)<<"\t"<<*(p+5)<<endl;

							//typedef itk::SymmetricEigenAnalysis< vnl_matrix<double> , vnl_diag_matrix<double>, vnl_matrix<double> > SymEigenAnalysisType;
							//SymEigenAnalysisType eigenD(3);

							//vnl_matrix<double> D(3,3);
							//vnl_matrix<double>		eVec(3,3);
							//vnl_diag_matrix<double> eVal(3) ;

							//costFunction->OneTensor2VnlMatrix(D, initTensor );
							//eigenD.ComputeEigenValuesAndVectors(D, eVal, eVec);
							//cout<<"Eigen values: "<<endl<<"  "<<eVal[0]<<" \t "<<eVal[1]<<" \t "<<eVal[2]<<endl<<endl;
						}

#ifdef OUTPUT_EVAL_IN_FILE
						//debugfile<<"  End pos f="<<costFunction->GetValue(optimizer->GetCurrentPosition())<< " IT:" << optimizer->GetStopConditionDescription().c_str() << endl;
						optimizer->Print(debugfile);
#endif
					}
					++tdit;

				} //while( !tdit.IsAtEnd() )
			
#ifdef OUTPUT_EVAL_IN_FILE
				debugfile<<endl<<endl<<endl<<endl<<endl<<endl<<endl;
#endif
				//--------------------------------------------
				// If regularization, wait for all the threads to
				// have finish this pass
				// (and also wait before going to STEP4 where
				//  neighbouring voxel data could be deleted)
				//--------------------------------------------
				if (  m_AlgoParameters.m_UseRegularization )
					m_BarrierSync->Wait();

			} //for ( int icmPass = 0; icmPass<NbMaxPasses ; icmPass++ )



			//--------------------------------------------------------------
			//--------------------------------------------------------------
			// STEP 4. AFTER THE MAIN LOOP
			//--------------------------------------------------------------
			//--------------------------------------------------------------

#ifdef OUTPUT_EVAL_IN_FILE
			debugfile.close();
#endif

			//--------------------------------------------------------------
			// Create the iterator for the residuals
			//--------------------------------------------------------------
			typedef itk::ImageRegionIterator< ResidualsImageType > ResidualsIteratorType;
			ResidualsIteratorType oitResiduals;
			
			if ( m_AlgoParameters.m_ExportResiduals )
				oitResiduals = ResidualsIteratorType( m_ResidualsImage, outputRegionForThread );

			//--------------------------------------------------------------
			// STEP 3. OUTPUT IN MAPS and free each InternalVoxelData
			//--------------------------------------------------------------
			std::vector<TensorPixelType > finalTensors;
			tdit.GoToBegin();
			for ( int i=0 ; i< m_AlgoParameters.m_NumberOfTensors ; i++ )	oits[i].GoToBegin(); 

			while ( !tdit.IsAtEnd() )
			{
				//--------------------------------------------
				// Some initialization for that voxel
				//--------------------------------------------
				typename TensorDataImageType::IndexType idx = tdit.GetIndex();
				typename OutputFractionsImage::IndexType idxF;
				idxF[0]=idx[0]; idxF[1]=idx[1]; idxF[2]=idx[2]; 

				InternalVoxelData* tensorData = tdit.GetCenterPixel();
				costFunction->SetInternalVoxelData(tensorData);

				//--------------------------------------------
				// Compute or not final tensor result depending if the tensor is valid
				//--------------------------------------------
				if ( tensorData )
				{
					// Compute the final values of the tensors, and output them in the image
					costFunction->ComputeFinalTensors(finalTensors);
				}
				else
				{ 
					finalTensors.resize( m_AlgoParameters.m_NumberOfTensors);
					for ( int i=0 ; i< m_AlgoParameters.m_NumberOfTensors ; i++ )
						finalTensors[i].Fill(0);
				}

				//--------------------------------------------
				// If more than 1 tensor, order them by fraction!
				//--------------------------------------------
				if ( tensorData && (tensorData->nbTensorAtVoxel==2) && (tensorData->normalizedFractions.size()>=2) )
				{
					// RQ: the possible isotropic part is in the last element of normalizedFractions
					// So there is no problem to order the first nbTensorAtVoxel tensors
					if ( tensorData->normalizedFractions[1] > tensorData->normalizedFractions[0] )
					{
						TensorPixelType tmp = finalTensors[0];
						finalTensors[0] = finalTensors[1];
						finalTensors[1] = tmp;

						float tmp2 = tensorData->normalizedFractions[0];
						tensorData->normalizedFractions[0]=tensorData->normalizedFractions[1];
						tensorData->normalizedFractions[1]=tmp2;
					}
				}

				//--------------------------------------------
				// Output the isotropic/anisotropic parts of the tensors
				//--------------------------------------------
				for ( int i=0 ; i< m_AlgoParameters.m_NumberOfTensors ; i++ )
				{
					oits[i].Set(finalTensors[i]);

					//// Put isotropic part/anisotropic part. even if tensor null, OK !
					//TensorPixelType I, A;
					//A = finalTensors[i]; I.Fill(0); 
					//double trace = A(0,0)+A(1,1)+A(2,2);
					//for ( int j=0 ; j<3 ; j++ )
					//{
					//	I(j,j) = 1.0/3.0 * trace;
					//	A(j,j) = A(j,j) - I(j,j);
					//}
					//oitsA[i].Set(A);
					//oitsI[i].Set(I);
				}

				//--------------------------------------------
				// Output fractions
				//--------------------------------------------
				if ( tensorData )
				{
					VnlFractionsType finalFractions;
					costFunction->ComputeFinalFractions(finalFractions);
					for ( unsigned int i=0; i<finalFractions.size() ; i++ )
					{
						idxF[3]=i;
						m_FractionsImage->SetPixel(idxF, finalFractions[i] );
					}

					if ( m_AlgoParameters.m_EstimationStrategy&ESTIMATION_RADIALBASELINE ) 
					{
						for ( int i=0 ; i< m_AlgoParameters.m_NumberOfTensors ; i++ )
						{
							idxF[3]=i + finalFractions.size();
							if ( i<tensorData->nbTensorAtVoxel )
								m_FractionsImage->SetPixel(idxF, tensorData->current_RadialCst[i]);
							else
								m_FractionsImage->SetPixel(idxF, 0);
						}
					}
				}
				else
				{
					for ( unsigned int i=0; i<m_FractionsImage->GetLargestPossibleRegion().GetSize()[3] ; i++ )
					{
						idxF[3]=0;
						m_FractionsImage->SetPixel(idxF, 0 );
					}
				}

				//--------------------------------------------
				// Output the residuals if necessary
				//--------------------------------------------
				if ( m_AlgoParameters.m_ExportResiduals )
				{
					typename ResidualsImageType::PixelType p(m_AlgoParameters.m_g.size());
					if ( tensorData )
					{
						std::vector<double> residuals;
						costFunction->ComputeResiduals(residuals);
						for ( unsigned int i=0; i<residuals.size() ; i++ )
							p[i] = residuals[i];
					}
					else
						p.Fill(0);


					oitResiduals.Set(p);
					++oitResiduals;
				}


				if ( tensorData )
				{
					//costFunction->TestCostFunction(tensorData->lastParameters); //costFunction->GetOptimizerParametersFromPrecomputeData());
				}

				//--------------------------------------------
				// Increment the iterators
				//--------------------------------------------
				++tdit;
				for ( int i=0 ; i< m_AlgoParameters.m_NumberOfTensors ; i++ )
					++oits[i]; 

				//------------------------------
				// Free all the memory for that voxel
				//------------------------------
				if ( tensorData )
				{
					tensorData->Neighbors.clear();
					tensorData->SiLs.clear();
					tensorData->current_GradientRadialProjection.clear();
					tensorData->current_GradientTensorProjection.clear();
					tensorData->current_RadialMatrix.clear();
					tensorData->current_RadialCst.clear();
					tensorData->Fractions.clear();
					tensorData->stabilized_L.clear();
					tensorData->current_L.clear();
					tensorData->current_EigenVectorsLE.clear();
					tensorData->current_EigenValuesLE.clear();
					//tensorData->D.clear();
					tensorData->normalizedFractions.clear();
					delete tensorData;
				}
			}


			m_BarrierSync->Wait();

		}	

		/**********************************************************************************************//**
		 * \brief	After threaded generate data 
		 *
		 * \author	Benoit Scherrer
		 * \date	July 2010
		 *
		 * \return	. 
		*************************************************************************************************/
		template< class TGradientImagePixelType, class TTensorPixelType >
		void MFMTensorEstimateFilter< TGradientImagePixelType, TTensorPixelType >
			::AfterThreadedGenerateData()
		{

		}


	

		/**********************************************************************************************//**
		 * \brief	Do one iteration with only small b values (defined by
		 * 			m_AlgoParameters.m_LowBValueMaxIndex). 
		 *
		 * \author	Benoit Scherrer
		 * \date	March 2010
		 *
		 * \param [in,out]	tdit			The tdit. 
		 * \param [in,out]	optimizer		If non-null, the optimizer. 
		 * \param [in,out]	costFunction	If non-null, the cost function. 
		 *
		 * \return	. 
		*************************************************************************************************/
		template< class TGradientImagePixelType, class TTensorPixelType >
		void MFMTensorEstimateFilter< TGradientImagePixelType, TTensorPixelType >
			::OnePassWithSmallBValues( TensorDataImageIteratorType& tdit, OptimizerType *optimizer, MFMTensorEstimateCostFunction *costFunction )
		{
			// Init the iterators for that pass
			tdit.GoToBegin();

			costFunction->SetUseRegularization(false);
			costFunction->SetUseRobustEstimation(false);
			costFunction->SetCurrentMaxBValueIndexToUse( m_AlgoParameters.m_LowBValueMaxIndex );

			// Pass on all pixels in the region
			while( !tdit.IsAtEnd() )
			{
				//taskProgress.Update();

				// Get a pointer to the tensor data
				InternalVoxelData* tensorData = tdit.GetCenterPixel();
				costFunction->SetInternalVoxelData(tensorData);

				//----------------------------------------
				// Launch the optimization for that tensor
				//----------------------------------------
				if ( tensorData )
				{
					try {
						// Get initial parameters from the value of the L's
						MFMTensorEstimateCostFunction::ParametersType initP = costFunction->GetOptimizerParametersFromPrecomputeData();
						optimizer->SetInitialPosition( initP );

						// Start the optimization !
						optimizer->StartOptimization();

						// Get last position
						MFMTensorEstimateCostFunction::ParametersType endP = optimizer->GetCurrentPosition();
						costFunction->EndAlgoIteration(endP);
					}
					catch (...)
					{
						cout<<">>> Exception f="<<costFunction->GetValue(optimizer->GetCurrentPosition())<<endl;
					}

					// TEST IF NAN
					if ( isnan(costFunction->GetValue(optimizer->GetCurrentPosition()) ) )
					{
						cout<<"NAN REACHED"<<endl<<"Initial one tensor solution: "<<endl;
					}

				}
				++tdit;

			} //while( !tdit.IsAtEnd() )

			tdit.GoToBegin();

		}

		/**********************************************************************************************//**
		 * \fn	template< class TGradientImagePixelType,
		 * 		class TTensorPixelType > void MFMTensorEstimateFilter< TGradientImagePixelType,
		 * 		TTensorPixelType > ::SaveForVideo( int passID, TensorDataImageIteratorType& tdit,
		 * 		MFMTensorEstimateCostFunction *costFunction )
		 *
		 * \brief	Save the current state of the two-tensors in a nrrd file.
		 * 			
		 * 			\warning Do work only with one thread. 
		 *
		 * \author	Benoit Scherrer
		 * \date	March 2010
		 *
		 * \param	passID					Identifier for the pass. 
		 * \param [in,out]	tdit			The tdit. 
		 * \param [in,out]	costFunction	If non-null, the cost function. 
		 *
		 * \return	. 
		*************************************************************************************************/
		template< class TGradientImagePixelType, class TTensorPixelType >
		void MFMTensorEstimateFilter< TGradientImagePixelType, TTensorPixelType >
			::SaveForVideo( int passID, TensorDataImageIteratorType& tdit, MFMTensorEstimateCostFunction *costFunction )
		{
			//!! WARNING. MUST USE ONLY ONE THREAD TO OUTPUT VIDEO
			for ( int k=0; k<m_AlgoParameters.m_NumberOfTensors ; k++ )
			{
				typename TensorImageType::Pointer imgVideo = TensorImageType::New(); 
				imgVideo->SetRegions(tdit.GetRegion());
				imgVideo->Allocate();

				typedef ImageRegionIterator< TensorImageType >		VideoIteratorType;
				VideoIteratorType vit(imgVideo, tdit.GetRegion() );
				vit.GoToBegin();
				tdit.GoToBegin();

				std::vector<TensorPixelType > finalTensors;

				// Pass on all pixels in the region
				while( !tdit.IsAtEnd() )
				{
					// Get a pointer to the tensor data
					InternalVoxelData* tensorData = tdit.GetCenterPixel();
					costFunction->SetInternalVoxelData(tensorData);

					if ( tensorData )
					{
						costFunction->ComputeFinalTensors(finalTensors);

					}
					else
					{
						finalTensors.resize( m_AlgoParameters.m_NumberOfTensors);
						for ( int i=0 ; i< m_AlgoParameters.m_NumberOfTensors ; i++ )
							finalTensors[i].Fill(0);
					}

					vit.Set(finalTensors[k]);

					++vit;
					++tdit;
				}

				typedef itk::ImageFileWriter<TensorImageType > VideoImgWriterType;
				typename VideoImgWriterType::Pointer imgWriter = VideoImgWriterType::New();
				char szImgFileName[256];
				sprintf(szImgFileName,"video_%03d_t%d.nrrd",passID,k);
				imgWriter->SetFileName(szImgFileName);
				imgWriter->SetInput(imgVideo);
				imgWriter->Update();


			}

			tdit.GoToBegin();

		}

//#define DEBUG_SAVE_GRAD
#ifdef DEBUG_SAVE_GRAD
				typedef itk::Image<float,3> DebugImageType;

				DebugImageType::Pointer debugImg = DebugImageType::New();
				debugImg->SetRegions(outputRegionForThread);
				debugImg->Allocate();
				typedef ImageRegionIterator< DebugImageType >		DebugIteratorType;

				DebugIteratorType dit(debugImg, outputRegionForThread );
				dit.GoToBegin();
				tdit.GoToBegin();

				// Pass on all pixels in the region
				while( !tdit.IsAtEnd() )
				{
					taskProgress.Update();

					// Get a pointer to the tensor data
					InternalVoxelData* tensorData = tdit.GetCenterPixel();



					if ( tensorData )
					{
						double regterm1=0;
						double regterm2=0;
						for ( int i=0; i<3 ; i++ )
						{
							InternalVoxelData *LA = tensorData->Neighbors[2*i];
							InternalVoxelData *LB = tensorData->Neighbors[2*i+1];

							if ( LA==NULL ) LA = tensorData;
							if ( LB==NULL ) LB = tensorData;

							//if ( !LA->validTensor ) LA = tensorData;
							//if ( !LB->validTensor ) LB = tensorData;


							vnl_matrix<double> m0A0 = (tensorData->L[0] - LA->L[0] );
							vnl_matrix<double> m0A1 = (tensorData->L[0] - LA->L[1] );
							vnl_matrix<double> m0B0 = (tensorData->L[0] - LB->L[0] );
							vnl_matrix<double> m0B1 = (tensorData->L[0] - LB->L[1] );

							vnl_matrix<double> m1A0 = (tensorData->L[1] - LA->L[0] );
							vnl_matrix<double> m1A1 = (tensorData->L[1] - LA->L[1] );
							vnl_matrix<double> m1B0 = (tensorData->L[1] - LB->L[0] );
							vnl_matrix<double> m1B1 = (tensorData->L[1] - LB->L[1] );

							double d0A0 = m0A0.frobenius_norm();
							double d0A1 = m0A1.frobenius_norm();
							double d0B0 = m0B0.frobenius_norm();
							double d0B1 = m0B1.frobenius_norm();				

							double d1A0 = m1A0.frobenius_norm();
							double d1A1 = m1A1.frobenius_norm();
							double d1B0 = m1B0.frobenius_norm();
							double d1B1 = m1B1.frobenius_norm();				

							// L[0] is associed with LA[0]->L[0]
							vnl_matrix<double> Reg0(3,3), Reg1(3,3);

							if ( d0A0<d0A1 )	Reg0 = LA->L[0];
							else				Reg0 = LA->L[1];
							if ( d0B0<d0B1 )	Reg0 = 0.5*(Reg0 - LB->L[0]);
							else				Reg0 = 0.5*(Reg0 - LB->L[1]);

							if ( d1A0<d1A1 )	Reg1 = LA->L[0];
							else				Reg1 = LA->L[1];
							if ( d1B0<d1B1 )	Reg1 = 0.5*(Reg1 - LB->L[0]);
							else				Reg1 = 0.5*(Reg1 - LB->L[1]);

							double d = (Reg0.frobenius_norm());
							regterm1 += d*d;
							d = (Reg1.frobenius_norm());
							regterm2 += d*d;

						}

						regterm1 = sqrt(regterm1);
						regterm2 = sqrt(regterm2);

						int nbN=0;
						for ( int i=0; i<6 ; i++ ) 
							if ( tensorData->Neighbors[i]!=NULL && (tensorData->Neighbors[i])) 
								nbN++;
						
						


						dit.Set(nbN);


						//dit.Set(sqrt(1+regterm1*regterm1/(0.25*0.25) ));
					}
					else
						dit.Set(0);

					++tdit;
					++dit;
				}

				tdit.GoToBegin();

				typedef itk::ImageFileWriter<DebugImageType > DebugWriterType;
				DebugWriterType::Pointer debugWriter = DebugWriterType::New();
				char szDebugFileName[256];
				sprintf(szDebugFileName,"Debug_p%d_iter%d.nrrd",procID,icmPass);
				debugWriter->SetFileName(szDebugFileName);
				debugWriter->SetInput(debugImg);
				debugWriter->Update();

				//debugWriter->Delete();
				//debugImg->Delete();
#endif

		template< class TGradientImagePixelType, class TTensorPixelType >
		void MFMTensorEstimateFilter< TGradientImagePixelType, TTensorPixelType >
			::SetInitTensorImage(  TensorImageType *image )
		{
			this->m_InitTensorImage = image;
			this->ProcessObject::SetNthInput( 1, 
				const_cast< TensorImageType* >(image) );

		}
		
		template< class TGradientImagePixelType, class TTensorPixelType >
		void MFMTensorEstimateFilter< TGradientImagePixelType, TTensorPixelType >
			::Set2TMaskImage( MaskImageType::Pointer maskImage )
		{
			m_2TMaskImage = maskImage;
		}

		/**********************************************************************************************//**
		 * \fn	template< class TGradientImagePixelType,
		 * 		class TTensorPixelType > void MFMTensorEstimateFilter<TGradientImagePixelType,
		 * 		TTensorPixelType > ::SetInputStudy( crl::DWIStudy<TGradientImagePixelType>& inputStudy )
		 *
		 * \brief	Sets the input study. 
		 *
		 * \author	Benoit Scherrer
		 * \date	December 2010
		 *
		 * \param [in,out]	inputStudy	The input study. 
		 *
		 * \return	. 
		*************************************************************************************************/
		template< class TGradientImagePixelType, class TTensorPixelType >
		void MFMTensorEstimateFilter<TGradientImagePixelType, TTensorPixelType >
			::SetInputStudy( crl::DWIStudy<TGradientImagePixelType>& inputStudy )
		{
			cout<<"Nominal B-Value = "<< inputStudy.NominalBValue<<endl;
			cout<<"Gradients:"<<endl;

			//--------------------------------------------------------------
			// Initialize
			//--------------------------------------------------------------
			m_AlgoParameters.m_BValues.clear();
			m_AlgoParameters.m_g.clear();
			m_AlgoParameters.m_G.clear();
			m_GradientInd.clear();
			m_BaselineInd.clear();

			//--------------------------------------------------------------
			// Prepare a table of indexes to short the b-values
			//--------------------------------------------------------------
			std::vector<int> indexes(inputStudy.BValues.size());
			for( unsigned int i=0; i<inputStudy.BValues.size(); i++ )
				indexes[i]=i;

			//--------------------------------------------------------------
			// First order by increasing b-values
			//--------------------------------------------------------------
			std::vector<double> tmpBVals = inputStudy.BValues;
			for ( unsigned int i=0; i<tmpBVals.size() ; i++ )
			{
				for ( unsigned int j=i+1; j<tmpBVals.size() ; j++ )
				{
					if ( tmpBVals[j] < tmpBVals[i] )
					{
						double v = tmpBVals[i];
						tmpBVals[i]=tmpBVals[j];
						tmpBVals[j] = v;

						int d=indexes[i];
						indexes[i]=indexes[j];
						indexes[j]=d;
					}
				}
			}

			double minB = 4000000;
			double maxB = 0.0;

			//--------------------------------------------------------------
			// Extract the indices and values of the baseline/gradient images
			// (a baseline image is a one with zero-norm gradients)
			//--------------------------------------------------------------
			for ( unsigned int i=0; i<tmpBVals.size() ; i++ )
			{
				//------------------------------------------
				// Base line image (B=0)
				//------------------------------------------
				if(tmpBVals[i] <= 0.0) 
				{
					m_BaselineInd.push_back(indexes[i]);
					cout<<inputStudy.Gradients[indexes[i]]<<"\tB-Value=0"<<endl;
				}
				//------------------------------------------
				// Gradient sensitized image
				//------------------------------------------
				else
				{
					m_GradientInd.push_back(indexes[i]);
					m_AlgoParameters.m_g.push_back( inputStudy.Gradients[indexes[i]] );
					m_AlgoParameters.m_BValues.push_back( inputStudy.BValues[indexes[i]] );
		
					// precomputed G = g g^t for each g
					VnlMatrixType G(3,3);
					crl::FastOps::ggt(G, inputStudy.Gradients[indexes[i]] );
					m_AlgoParameters.m_G.push_back(G);

					double b = inputStudy.BValues[indexes[i]];
					if ( b>maxB ) maxB = b;
					if ( b<minB ) minB = b;

					cout << inputStudy.Gradients[indexes[i]] << "\tB-Value=" << b << endl;
				}
			}

			// Set the threshold for what is considered a low b-value
			if ( m_AlgoParameters.m_FirstPassWithLowBValuesOnly )
			{
				m_AlgoParameters.m_LowBValueMaxIndex=0;
				double lowBValueThreshold = minB+minB*0.1;
				for ( unsigned int i=0; i<m_AlgoParameters.m_BValues.size() ; i++ )
				{
					if ( m_AlgoParameters.m_BValues[i]>lowBValueThreshold )
					{
						m_AlgoParameters.m_LowBValueMaxIndex = i;
						break;
					}
				}
			
				if ( m_AlgoParameters.m_LowBValueMaxIndex==0 ) 
				{
					m_AlgoParameters.m_LowBValueMaxIndex=m_AlgoParameters.m_BValues.size();
					m_AlgoParameters.m_FirstPassWithLowBValuesOnly=false;
					cout << " - Did not find maximum b-value for the first pass."<<endl;
				}
				else
					cout<<" - Maximum B-Value for first pass: "<<m_AlgoParameters.m_BValues[m_AlgoParameters.m_LowBValueMaxIndex-1]<<endl;

			}
			else
				m_AlgoParameters.m_LowBValueMaxIndex = m_AlgoParameters.m_BValues.size();



			//--------------------------------------------------------------
			// Set the gradient images to be the input 0
			//--------------------------------------------------------------
			this->ProcessObject::SetNthInput( 0, const_cast< GradientImagesType* >(inputStudy.DWIData.GetPointer()) );
		}


		template< class TGradientImagePixelType, class TTensorPixelType >
			void MFMTensorEstimateFilter< TGradientImagePixelType, TTensorPixelType >
			::PrintSelf(std::ostream& os, Indent indent) const
		{
			//Superclass::PrintSelf(os,indent);

			os << indent << "Threshold for reference B0 image: " << m_AlgoParameters.m_Threshold << std::endl << std::endl;

			os << indent << "OptimizationAlgo: ";
			os << crl::OptimizerUtils::GetOptimizerStringFromId(m_AlgoParameters.m_OptimizationAlgo)<<std::endl;

			os << indent << "ROBUST   : "<< m_AlgoParameters.m_RobustEstimation<<std::endl;
			os << indent << "MAX ITER: "<< m_AlgoParameters.m_MaxIter<<std::endl<<std::endl;

			os << indent << "FIRST PASS WITH ONLY LOWER B-VALUE: "<< m_AlgoParameters.m_FirstPassWithLowBValuesOnly<<std::endl<<std::endl;

			os << indent << "USE RICIAN MODELING:    "<< m_AlgoParameters.m_UseRicianNoiseModeling<<std::endl;
			os << indent << "USE RICIAN VARIANCE:    "<< m_AlgoParameters.m_RicianNoiseVariance<<std::endl<<std::endl;

			os << indent << "USE REGULARIZATION:    "<< m_AlgoParameters.m_UseRegularization<<std::endl;
			os << indent << "REGULARIZATION WEIGHT: "<< m_AlgoParameters.m_RegularizationWeight<<std::endl<<std::endl;

			os << indent << "ESTIMATION STRATEGY: "<< (int)m_AlgoParameters.m_EstimationStrategy<<std::endl;
			os << indent << "   SAME_MAJOR_EIGEN_VALUE = " << (int)(m_AlgoParameters.m_EstimationStrategy&ESTIMATION_SAME_MAJOR_EIGEN_VALUE)<<std::endl;
			os << indent << "   SAME_MINOR_EIGEN_VALUE = " << (int)(m_AlgoParameters.m_EstimationStrategy&ESTIMATION_SAME_MINOR_EIGEN_VALUE)<<std::endl;
			os << indent << "   ESTIMATION_EACH_CYLINDER = " << (int)(m_AlgoParameters.m_EstimationStrategy&ESTIMATION_EACH_CYLINDER)<<std::endl;
			os << indent << "   ESTIMATION_SAME_PLANE = " << (int)(m_AlgoParameters.m_EstimationStrategy&ESTIMATION_SAME_PLANE)<<std::endl;
			os << indent << "   F_ISO = " << (int)(m_AlgoParameters.m_EstimationStrategy&ESTIMATION_F_ISO)<<std::endl;
			os << indent << "   ESTIMATION_F_TENSORS = " << (int)(m_AlgoParameters.m_EstimationStrategy&ESTIMATION_F_TENSORS)<<std::endl;

		}

	}
}
#endif
