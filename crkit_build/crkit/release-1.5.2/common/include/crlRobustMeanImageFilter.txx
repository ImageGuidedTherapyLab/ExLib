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


#ifndef __crlRobustMeanImageFilter_txx
#define __crlRobustMeanImageFilter_txx

#include "crlRobustMeanImageFilter.h"
#include "crlRobustMeanCostFunction.h"
#include "crlTaskProgress.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkProgressReporter.h"
#include "itkConjugateGradientOptimizer.h"
#include "itkFRPROptimizer.h"

namespace crl{

/**********************************************************************************************//**
*
 * \brief	Robust mean image filter constructor. 
 *
 * \author	Benoit Scherrer
 * \date	April 2010
 *
 * \return	. 
*************************************************************************************************/
template <class TInputImage, class TOutputImage>
RobustMeanImageFilter<TInputImage,TOutputImage>
::RobustMeanImageFilter()
{
	this->SetNumberOfOutputs ( 1 );
	this->SetNumberOfRequiredInputs(1);
	m_OutputMeanImage = false;
	m_OutputDiffImage = false;
	m_ShowTaskProgress = true;
}

/**********************************************************************************************//**
 * \fn	template <class TInputImage, class TOutputImage> RobustMeanImageFilter<TInputImage,
 * 		TOutputImage> ::~RobustMeanImageFilter()
 *
 * \brief	Destructor. 
 *
 * \author	Benoit Scherrer
 * \date	April 2010
*************************************************************************************************/
template <class TInputImage, class TOutputImage>
RobustMeanImageFilter<TInputImage,TOutputImage>
::~RobustMeanImageFilter()
{
}

/**********************************************************************************************//**
 * \brief	Function called after the buffer has been allocated but before threads are spawned
 *
 * \author	Benoit Scherrer
 * \date	April 2010
 *
 * \exception	itk::ExceptionObject	Thrown when exception. 

*************************************************************************************************/
template <class TInputImage, class TOutputImage>
void 
RobustMeanImageFilter<TInputImage,TOutputImage>	
::BeforeThreadedGenerateData()
{
	if ( this->GetNumberOfInputs() < 1 ) 
		throw itk::ExceptionObject(__FILE__,__LINE__,"No images were added to compute the robust mean.", "RobustMeanImageFilter");

	// Take the first image as a reference 
	InputImageConstPointer oneImage = this->GetInput(0);

	/*-----------------------------------------
	Allocate the robust mean image
	-----------------------------------------*/
	m_RobustMeanEstimate = this->GetOutput();
	m_RobustMeanEstimate->SetRegions(oneImage->GetLargestPossibleRegion());
	m_RobustMeanEstimate->SetOrigin(oneImage->GetOrigin());
	m_RobustMeanEstimate->SetSpacing(oneImage->GetSpacing());
	m_RobustMeanEstimate->SetDirection(oneImage->GetDirection());
	m_RobustMeanEstimate->Allocate();

	/*-----------------------------------------
	Allocate mean image
	-----------------------------------------*/
	m_MeanImage = OutputImageType::New();
	if ( m_OutputMeanImage )
	{
		m_MeanImage->SetRegions(oneImage->GetLargestPossibleRegion());
		m_MeanImage->SetOrigin(oneImage->GetOrigin());
		m_MeanImage->SetSpacing(oneImage->GetSpacing());
		m_MeanImage->SetDirection(oneImage->GetDirection());
		m_MeanImage->Allocate();
	}

	/*-----------------------------------------
	Allocate diff image
	-----------------------------------------*/
	m_DiffImage = OutputImageType::New();
	if ( m_OutputDiffImage )
	{
		m_DiffImage->SetRegions(oneImage->GetLargestPossibleRegion());
		m_DiffImage->SetOrigin(oneImage->GetOrigin());
		m_DiffImage->SetSpacing(oneImage->GetSpacing());
		m_DiffImage->SetDirection(oneImage->GetDirection());
		m_DiffImage->Allocate();
	}

	if ( m_ShowTaskProgress ) TaskProgress::InitTaskProgress();
}

#define CONJGRAD 1

/**********************************************************************************************//**
 *
 * \brief	Implement the threaded robust mean image filer
 *
 * \author	Benoit Scherrer
 * \date	April 2010
 *
 * \param	outputRegionForThread	The output region for that thread. 
 * \param	procID					Identifier for the thread. 
*************************************************************************************************/
template <class TInputImage, class TOutputImage>
void 
RobustMeanImageFilter<TInputImage,TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int procID ) 
{
	TaskProgress taskProgress(procID, this->GetNumberOfThreads(), outputRegionForThread.GetSize()[0]*outputRegionForThread.GetSize()[1]*outputRegionForThread.GetSize()[2]);
	taskProgress.SetLinePrefix("     ");

	/*-----------------------------------------
	Creates the cost function for that thread
	-----------------------------------------*/
	RobustMeanCostFunction::Pointer costFunction = RobustMeanCostFunction::New();

#ifndef CONJGRAD
	/*-----------------------------------------
	Use the FRPR optimization method of ITK
	-----------------------------------------*/

	typedef itk::FRPROptimizer					OptimizerType;
	OptimizerType::Pointer optimizer = OptimizerType::New();
	optimizer->SetCostFunction(costFunction);
	optimizer->SetStepLength(1e-5);
	optimizer->SetStepTolerance (1e-7);
	optimizer->SetValueTolerance(1e-7); 
	optimizer->SetMaximumIteration(1000);
	optimizer->SetMaximize(false);
#else
	/*-----------------------------------------
	Use the conjuguate gradient method of ITK
	-----------------------------------------*/
	typedef itk::ConjugateGradientOptimizer			OptimizerType;
	typedef  OptimizerType::InternalOptimizerType	vnlOptimizerType;

	OptimizerType::Pointer optimizer = OptimizerType::New();
	optimizer->SetCostFunction(costFunction);
	optimizer->MinimizeOn();

	vnlOptimizerType * vnlOptimizer = optimizer->GetOptimizer();
	const double F_Tolerance      = 1e-3;  // Function value tolerance
	const double G_Tolerance      = 1e-4;  // Gradient magnitude tolerance
	const double X_Tolerance      = 1e-8;  // Search space tolerance
	const double Epsilon_Function = 1e-10; // Step
	const int    Max_Iterations   =   1000; // Maximum number of iterations
	vnlOptimizer->set_f_tolerance( F_Tolerance );
	vnlOptimizer->set_g_tolerance( G_Tolerance );
	vnlOptimizer->set_x_tolerance( X_Tolerance );
	vnlOptimizer->set_epsilon_function( Epsilon_Function );
	vnlOptimizer->set_max_function_evals( Max_Iterations );
	vnlOptimizer->set_check_derivatives(1 );
#endif

	/*-----------------------------------------
	Prepare iterators
	-----------------------------------------*/
	// Output iterators
	ImageRegionIterator<TOutputImage> itOutput(m_RobustMeanEstimate, outputRegionForThread );

	OutputImageRegionType nullRegion;
	ImageRegionIterator<TOutputImage> itMean(m_MeanImage, m_OutputMeanImage ? outputRegionForThread:nullRegion );
	ImageRegionIterator<TOutputImage> itDiff(m_DiffImage, m_OutputDiffImage ? outputRegionForThread:nullRegion );

	// As many input iterators as inputs
	std::vector< ImageRegionConstIterator<TInputImage> > itsInput;
	for ( unsigned int i=0; i<this->GetNumberOfInputs() ; i++ )
		itsInput.push_back(ImageRegionConstIterator<TInputImage>(this->GetInput(i), outputRegionForThread ));


	/*-----------------------------------------
	 Compute the robust mean for each voxels!
	-----------------------------------------*/
	itk::Array< double > optimizerInitPos(1);
	while( !itOutput.IsAtEnd() )
	{
		if ( m_ShowTaskProgress ) taskProgress.Update();

		double sumValue = 0.0;
		double nbValues = 0;
		double result = 0;

		//----------------------------------------------------
		// If only one input....
		//----------------------------------------------------
		if ( this->GetNumberOfInputs()==1 )
		{
			nbValues = 1;
			result = itsInput[0].Value();
			sumValue = result;

			++itsInput[0];
		}
		//----------------------------------------------------
		// If more than one input....
		//----------------------------------------------------
		else
		{
			//----------------------------------------------------
			// First compute the mean to initialize the non-linear optimization scheme
			// (and fill the observations for the metric)
			//-----------------------------------------------------
			costFunction->ClearObservations();
			for ( unsigned int i=0; i<this->GetNumberOfInputs() ; i++ )
			{
				double obs = itsInput[i].Value();
				costFunction->AddObservation(obs);	

				sumValue += obs;
				nbValues += 1.0;
				++itsInput[i];
			}
		
			//----------------------------------------------------
			// Estimate the robust parameter from the observations
			// ( =1.48 median||errs||)
			//----------------------------------------------------
			costFunction->EstimateT();

			//----------------------------------------------------
			// Initialize the gradient descent with the mean
			// and start the optimization
			//----------------------------------------------------
			optimizerInitPos[0] = sumValue/nbValues;
			optimizer->SetInitialPosition( optimizerInitPos );
			optimizer->StartOptimization();

			//----------------------------------------------------
			// Get the result
			//----------------------------------------------------
			result = optimizer->GetCurrentPosition()[0];
		}

		//----------------------------------------------------
		// Output the value
		//----------------------------------------------------
		itOutput.Set( (typename TOutputImage::PixelType)result );
		++itOutput;

		//----------------------------------------------------
		// If needed, output the (non robust) mean value
		//----------------------------------------------------
		if ( m_OutputMeanImage )
		{
			itMean.Set( (typename TOutputImage::PixelType) (sumValue/nbValues));
			++itMean;
		}

		//----------------------------------------------------
		// If needed, output the diff between the mean and the robust mean
		//----------------------------------------------------
		if ( m_OutputDiffImage )
		{
			itDiff.Set( (typename TOutputImage::PixelType) std::abs(result - sumValue/nbValues));
			++itDiff;
		}

	}
}

/**********************************************************************************************//**
 * \brief	Print information about the object, including superclasses
 *
 * \author	Benoit Scherrer
 * \date	April 2010
 *
 * \param [in,out]	os	The output stream. 
 * \param	indent		The indent. 
*************************************************************************************************/
template <class TInputImage, class TOutputImage>
void 
RobustMeanImageFilter<TInputImage,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "TO DO " << std::endl;
}


} /* end namespace crl. */

#endif





  




