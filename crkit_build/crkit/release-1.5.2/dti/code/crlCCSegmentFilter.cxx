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

#include "itkSymmetricEigenAnalysis.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

#include "crlDWICommon.h"
#include "crlCCSegmentFilter.h"

#define safe_log(x) ((x>1e-15?log(x):log(1e-15)))


using namespace itk;

namespace crl { 

	/**********************************************************************************************//**
	 * \fn	CCSegmentFilter::CCSegmentFilter()
	 *
	 * \brief	Default constructor. 
	 *
	 * \author	Benoit Scherrer
	 * \date	April 2011
	*************************************************************************************************/
	CCSegmentFilter::CCSegmentFilter()
	{
		this->SetNumberOfRequiredInputs( 1 ); 
	}

	/**********************************************************************************************//**
	 * \fn	CCSegmentFilter::~CCSegmentFilter()
	 *
	 * \brief	Destructor. 
	 *
	 * \author	Benoit Scherrer
	 * \date	April 2011
	*************************************************************************************************/
	CCSegmentFilter::~CCSegmentFilter()
	{
	}	

	/**********************************************************************************************//**
	 * \fn	void CCSegmentFilter::BeforeThreadedGenerateData()
	 *
	 * \brief	Before threaded generate data. 
	 *
	 * \author	Benoit Scherrer
	 * \date	April 2011
	*************************************************************************************************/
	void CCSegmentFilter::BeforeThreadedGenerateData()
	{
		m_InputTensor = (TensorImageType *) this->GetInput();

		m_OutputSegm = this->GetOutput();
		m_OutputSegm = OutputImageType::New();
		m_OutputSegm->SetOrigin(m_InputTensor->GetOrigin());
		m_OutputSegm->SetSpacing(m_InputTensor->GetSpacing());
		m_OutputSegm->SetDirection(m_InputTensor->GetDirection());
		m_OutputSegm->SetRegions(m_InputTensor->GetLargestPossibleRegion());
		m_OutputSegm->Allocate();
	}

	/**********************************************************************************************//**
	 * \fn	void CCSegmentFilter::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
	 * 		int procID)
	 *
	 * \brief	Threaded generate data. 
	 *
	 * \author	Benoit Scherrer
	 * \date	April 2011
	 *
	 * \param	outputRegionForThread	The output region for thread. 
	 * \param	procID					Identifier for the proc. 
	*************************************************************************************************/
	void CCSegmentFilter::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int procID) 
	{
		SymmetricEigenAnalysisType	EigenAnalysis(3);

		//---------------------------------------------------------------
		// Create the iterators
		//---------------------------------------------------------------
		typedef itk::ImageRegionIterator<TensorImageType> TensorIteratorType;
		TensorIteratorType itIn(m_InputTensor, outputRegionForThread);

		typedef itk::ImageRegionIterator<OutputImageType> OutputIteratorType;
		OutputIteratorType itOut(m_OutputSegm, outputRegionForThread);

		//---------------------------------------------------------------
		// Prepare some data
		//---------------------------------------------------------------
		VnlTensorType			D(3,3);
		TensorPixelType			tensor;
		VnlEigValueType			eVal(3);
		VnlEigVectorType		eVec(3,3);

		TaskProgress taskProgress(procID, this->GetNumberOfThreads(), outputRegionForThread.GetSize()[0]*outputRegionForThread.GetSize()[1]*outputRegionForThread.GetSize()[2]);

		//---------------------------------------------------------------
		// For all the tensor of the input region
		//---------------------------------------------------------------
		while ( !itOut.IsAtEnd() )
		{
			taskProgress.Update();

			//--------------------------------------
			// Take the tensor 
			//--------------------------------------
			tensor = itIn.Value();
			crl::DWI::Tensor2VnlMatrix(D,tensor);	
			EigenAnalysis.ComputeEigenValuesAndVectors(D, eVal, eVec ); 

			double FA = tensor.GetFractionalAnisotropy();
			//double lambda1 = eVal[2];

			//double red = eVec[2][0]/(eVec[2][0]+eVec[2][1]+eVec[2][2]);
			double norm = sqrt(eVec[2][0]*eVec[2][0]+eVec[2][1]*eVec[2][1]+eVec[2][2]*eVec[2][2]);
			double red;

 			if ( norm > 1e-2 ) red = eVec[2][0]/norm;
			else red=0;

			if ( red<0 ) red=-red;

			//---------------------------------------------------
			// 
			//---------------------------------------------------
			if ( red*FA>0.5 )  itOut.Set(1);
			else itOut.Set(0);

			//---------------------------------------------------
			// Increment the iterators
			//---------------------------------------------------
			++itIn;
			++itOut;
		}
	}

	/**********************************************************************************************//**
	 * \fn	void CCSegmentFilter::AfterThreadedGenerateData ()
	 *
	 * \brief	After threaded generate data. 
	 *
	 * \author	Benoit Scherrer
	 * \date	April 2011
	*************************************************************************************************/
	void CCSegmentFilter::AfterThreadedGenerateData ()
	{
		typedef itk::ConnectedComponentImageFilter< OutputImageType, OutputImageType, OutputImageType > ConnectedComponentsFilterType;
		ConnectedComponentsFilterType::Pointer filter = ConnectedComponentsFilterType::New();
		filter->SetInput(m_OutputSegm);
		filter->SetBackgroundValue(0);
		filter->Update();

		typedef itk::RelabelComponentImageFilter<OutputImageType, OutputImageType> RelabelFilterType;
		RelabelFilterType::Pointer filter2 = RelabelFilterType::New();
		filter2->SetInput(filter->GetOutput());
		//filter2->SetMinimumObjectSize(40);
		filter2->Update();

		unsigned long int maxId=1;
		for ( unsigned int i=1; i<filter2->GetNumberOfObjects(); i++ )
		{
			if ( filter2->GetSizeOfObjectsInPixels ()[i] > filter2->GetSizeOfObjectsInPixels ()[maxId] )
			{
				maxId = i;
			}
		}
		std::cout<<"- Max component = "<<maxId<<"/"<< filter2->GetNumberOfObjects() <<std::endl;

		typedef itk::BinaryThresholdImageFilter<OutputImageType, OutputImageType> BinaryFilterType;
		BinaryFilterType::Pointer filter3=BinaryFilterType::New();
		filter3->SetInput(filter2->GetOutput());
		filter3->SetUpperThreshold(maxId);
		filter3->SetLowerThreshold(maxId);
		filter3->SetInsideValue(1);
		filter3->SetOutsideValue(0);
		filter3->Update();

		this->SetNthOutput(0, filter3->GetOutput());
	}

	/**********************************************************************************************//**
	 * \fn	void CCSegmentFilter::PrintSelf(std::ostream& os, Indent indent) const
	 *
	 * \brief	Print self. 
	 *
	 * \author	Benoit Scherrer
	 * \date	April 2011
	 *
	 * \param [in,out]	os	The operating system. 
	 * \param	indent		The indent. 
	*************************************************************************************************/
	void CCSegmentFilter::PrintSelf(std::ostream& os, Indent indent) const
	{
		//Superclass::PrintSelf(os,indent);
		os << indent << "TO DO  " << std::endl;
	}

}   // end namespace crl

