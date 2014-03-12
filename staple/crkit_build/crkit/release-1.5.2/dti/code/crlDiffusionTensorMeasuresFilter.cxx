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
#include "crlDiffusionTensorMeasuresFilter.h"

#define safe_log(x) ((x>1e-15?log(x):log(1e-15)))


using namespace itk;

namespace crl { 

	/**********************************************************************************************//**
	 * \fn	DiffusionTensorMeasuresFilter::DiffusionTensorMeasuresFilter()
	 *
	 * \brief	Default constructor. 
	 *
	 * \author	Benoit Scherrer
	 * \date	April 2011
	*************************************************************************************************/
	DiffusionTensorMeasuresFilter::DiffusionTensorMeasuresFilter()
	{
		this->SetNumberOfRequiredInputs( 1 ); 
	}

	/**********************************************************************************************//**
	 * \fn	DiffusionTensorMeasuresFilter::~DiffusionTensorMeasuresFilter()
	 *
	 * \brief	Destructor. 
	 *
	 * \author	Benoit Scherrer
	 * \date	April 2011
	*************************************************************************************************/
	DiffusionTensorMeasuresFilter::~DiffusionTensorMeasuresFilter()
	{
	}	

	/**********************************************************************************************//**
	 * \fn	void DiffusionTensorMeasuresFilter::BeforeThreadedGenerateData()
	 *
	 * \brief	Before threaded generate data. 
	 *
	 * \author	Benoit Scherrer
	 * \date	April 2011
	*************************************************************************************************/
	void DiffusionTensorMeasuresFilter::BeforeThreadedGenerateData()
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
	 * \fn	void DiffusionTensorMeasuresFilter::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
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
	void DiffusionTensorMeasuresFilter::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int procID) 
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

			//--------------------------------------
			// Compute the linear component
			//--------------------------------------
			double linearComponent;
			double sum = eVal[0]+eVal[1]+eVal[2];
			if ( sum<=1e-5 ) linearComponent=0;
			else linearComponent = 1-( eVal[2]-eVal[1] ) / sum;
		
			itOut.Set(linearComponent);

			//---------------------------------------------------
			// Increment the iterators
			//---------------------------------------------------
			++itIn;
			++itOut;
		}
		
	}

	void DiffusionTensorMeasuresFilter::AfterThreadedGenerateData ()
	{
		this->SetNthOutput(0, m_OutputSegm );
	}

	/**********************************************************************************************//**
	 * \fn	void DiffusionTensorMeasuresFilter::PrintSelf(std::ostream& os, Indent indent) const
	 *
	 * \brief	Print self. 
	 *
	 * \author	Benoit Scherrer
	 * \date	April 2011
	 *
	 * \param [in,out]	os	The operating system. 
	 * \param	indent		The indent. 
	*************************************************************************************************/
	void DiffusionTensorMeasuresFilter::PrintSelf(std::ostream& os, Indent indent) const
	{
		//Superclass::PrintSelf(os,indent);
		os << indent << "TO DO  " << std::endl;
	}

}   // end namespace crl

