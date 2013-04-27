/**********************************************************************************************//**
 * \file	crlTaskProgress.h
 *
 * \brief	Declares the crl::TaskProgress class. 
*************************************************************************************************/

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


#ifndef h_CRL_TASK_PROGRESS
#define h_CRL_TASK_PROGRESS

#ifdef WIN32
#else
#include <pthread.h>
#endif

#include <string>
#include <itkTimeProbe.h>
#include <itkCommand.h>
#include <itkProcessObject.h>

namespace crl {

	/**********************************************************************************************//**
	 * \class	TaskProgress
	 *
	 * \brief	Class used to show the progress of a (multi-threaded) task. 
	 * It's really easy to use. For example, for an multi-threaded itk::ImageToImageFilter:
	 *
	 * \code
	 * template <class TInputImage, class TOutputImage>
	 * void 
	 * MyImageToImageFilter<TInputImage,TOutputImage>
	 * ::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int procID ) 
	 * {
	 *     // Create a crl::TaskProgress object for each thread
	 *     unsigned long int maxCnt = outputRegionForThread.GetSize()[0]*outputRegionForThread.GetSize()[1]*outputRegionForThread.GetSize()[2]);
	 *     crl::TaskProgress taskProgress(procID, this->GetNumberOfThreads(), maxCnt);
	 *     taskProgress.SetLinePrefix("     ");
	 * 
	 *     (...)
	 *  
	 *     // if iterator is an iterator on outputRegionForThread, then update the crl::TaskProgress object at each iteration
	 *     while( !iterator.IsAtEnd() )
	 *     {
	 *          taskProgress.Update();
	 *          (...)
	 *     }
	 * \endcode
	 *
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	*************************************************************************************************/
	class TaskProgress
	{
	public:
		TaskProgress(unsigned int threadID, unsigned int nbThreads, unsigned long int counterMaxVal, bool showRemainingTime=true, bool showInitialProgress=true );
		~TaskProgress();
		
		void		InitShowProgress();
		void		Update();
		void		SetCurrentCounter(unsigned int counter);
	
		void		SetShowEstimatedRemainingTime(bool b);
		bool		GetShowEstimatedRemainingTime() const;

		static std::string ConvertTimeToString(itk::RealTimeClock::TimeStampType t) ;
		static		void InitTaskProgress();

		void		SetLinePrefix(const std::string& prefix);
		std::string GetLinePrefix() const;


	protected:
		void		ShowProgress();
		void		SetCursorPosition();

	private:
		class		PrivateData;
		PrivateData *d;

	};

	/**********************************************************************************************//**
	 * \class	ItkFilterProgressObserver
	 *
	 * \brief	Defines a itk::Command that can be plugged to a ITK filter to show the progress of
	 * 			the filter.
	 * 			
	 * 			EXAMPLE:  
	 * 			ItkFilterProgressObserver::Pointer progressObs = ItkFilterProgressObserver::New();  
	 * 			progressObs->SetProgressLinePrefix("Processing... ");  
	 *			filter->AddObserver(ProgressEvent(), progressObs );  
	 * 			filter->AddObserver(StartEvent(), progressObs );  
	 *	
	 * \author	Benoit Scherrer
	 * \date	June 2010
	*************************************************************************************************/
	class ItkFilterProgressObserver : public itk::Command 
	{
	public:
		typedef  ItkFilterProgressObserver  Self;
		typedef  itk::Command				Superclass;
		typedef  itk::SmartPointer<Self>	Pointer;
		typedef   itk::ProcessObject		ProcessType;
		typedef   const ProcessType   *		ProcessPointer;

		itkNewMacro( Self );

	public:
		void	Execute(itk::Object *caller, const itk::EventObject & event);
		void	Execute(const itk::Object * object, const itk::EventObject & event);
		void	SetProgressLinePrefix(const std::string& str);

	protected:
		ItkFilterProgressObserver();
		TaskProgress m_TaskProgress;
	};

}

#endif
