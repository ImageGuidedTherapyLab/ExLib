/**********************************************************************************************//**
 * \file	crlTaskProgress.cxx
 *
 * \brief	The crl::TaskProgress class. 
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


#include "crlTaskProgress.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <cstdio>

#ifdef WIN32
  #include <windows.h>
  #pragma warning( disable: 4996 )
#endif

using namespace std;

namespace crl {

/*--------------------------------------------------------------
 Some VT100 macros (don t work for win32)
--------------------------------------------------------------*/
#define reset "\33c"
#define clear_screen "\33[2J"
#define cursor_up(    count ) "\33[" << count << "A"
#define cursor_down(  count ) "\33[" << count << "B"
#define cursor_right( count ) "\33[" << count << "C"
#define cursor_left(  count ) "\33[" << count << "D"


/*--------------------------------------------------------------
 Some global variables for all threads
--------------------------------------------------------------*/
int g_LastProcToCout = -1;		
int g_NbFinished = 0;			

// Declare a mutex (plateform dependent)
#ifdef WIN32
	int		g_InitYLine = 0;
	HANDLE  g_Mutex = CreateMutex(NULL, FALSE, NULL);
#else
	pthread_mutex_t g_Mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

/**********************************************************************************************//**
 * \class	TaskProgress::PrivateData
 *
 * \brief	The private data for FileName. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
*************************************************************************************************/
class TaskProgress::PrivateData
	{
	public:
		PrivateData(unsigned int threadID, unsigned int nbThreads, unsigned long int counterMaxVal, bool showRemainingTime):
			m_ThreadID(threadID),
			m_NbThreads(nbThreads),
			m_CounterMaxVal(counterMaxVal),
			m_CurrentCounter(0),
			m_CounterBeforeUpdate((long int)(m_CounterMaxVal/100 + 0.5)),
			m_ShowRemainingTime(showRemainingTime),
			m_LinePrefix("")
			{}

		unsigned int		m_ThreadID;
		unsigned int		m_NbThreads;
		unsigned long int	m_CounterMaxVal;
		unsigned long int	m_CurrentCounter;
		long int			m_CounterBeforeUpdate;

		itk::TimeProbe		m_TimeProbe;

		bool				m_ShowRemainingTime;
		std::string			m_LinePrefix;
	};

	/**********************************************************************************************//**
	 * \fn	TaskProgress::TaskProgress(unsigned int threadID, unsigned int nbThreads,
	 * 		unsigned long int counterMaxVal, bool showRemainingTime)
	 *
	 * \brief	Constructor. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \param	threadID			Thread identifier. 
	 * \param	nbThreads			The number of threads. 
	 * \param	counterMaxVal		The maximum value of the counter (corresponding to 100%). 
	 * \param	showRemainingTime	If true, will show the estimated remaining time. 
	*************************************************************************************************/
	TaskProgress::TaskProgress(unsigned int threadID, unsigned int nbThreads, unsigned long int counterMaxVal, bool showRemainingTime, bool showInitialProgress)
	{
		d = new TaskProgress::PrivateData(threadID, nbThreads, counterMaxVal, showRemainingTime);

		if ( showInitialProgress )
			InitShowProgress();

	}

	/**********************************************************************************************//**
	 * \fn	TaskProgress::~TaskProgress()
	 *
	 * \brief	Destructor. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	*************************************************************************************************/
	TaskProgress::~TaskProgress()
	{
		delete d;
	}

	/**********************************************************************************************//**
	 * \fn	void TaskProgress::InitTaskProgress()
	 *
	 * \brief	Initialises the task progress class (initialize global variables shared among
	 * 			each instance of crl::TaskProgress in each threads) 
	 *
	 * \author	Benoit Scherrer
	 * \date	April 2010
	*************************************************************************************************/
	void TaskProgress::InitTaskProgress()
	{
		g_LastProcToCout = -1;		
		g_NbFinished = 0;			

#ifdef WIN32
		g_InitYLine = 0;
#endif
	}

	//void TaskProgress::ClearScreen()
	//{
	//	cout<<clear_screen;
	//}
	void TaskProgress::InitShowProgress()
	{
		d->m_TimeProbe.Start();
		ShowProgress();

	}

	/**********************************************************************************************//**
	 * \fn	void TaskProgress::Update()
	 *
	 * \brief	Increment the internal counter and print out the result if needed
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	*************************************************************************************************/
	void TaskProgress::Update()
	{
		if ( d->m_CurrentCounter==0 ) d->m_TimeProbe.Start();

		d->m_CounterBeforeUpdate--;
		d->m_CurrentCounter++;

		if ( d->m_CounterBeforeUpdate<=0 || d->m_CurrentCounter==d->m_CounterMaxVal ) 
		{
			ShowProgress();
			d->m_CounterBeforeUpdate = (long int)(d->m_CounterMaxVal/100 + 0.5);
		}
	}

	/**********************************************************************************************//**
	 * \fn	void TaskProgress::SetCurrentCounter(unsigned int counter)
	 *
	 * \brief	Sets the value of the counter and print out the result. 
	 *
	 * \author	Benoit Scherrer
	 * \date	June 2010
	 *
	 * \param	counter	The counter. 
	*************************************************************************************************/
	void TaskProgress::SetCurrentCounter(unsigned int counter)
	{
		if ( d->m_CurrentCounter==0 ) d->m_TimeProbe.Start();

		d->m_CurrentCounter = counter;
		ShowProgress();
	}


	/**********************************************************************************************//**
	 * \fn	void TaskProgress::ShowProgress()
	 *
	 * \brief	Shows the progress for the current thread. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	*************************************************************************************************/
	void TaskProgress::ShowProgress()
	{
#ifdef WIN32
		WaitForSingleObject(g_Mutex, INFINITE);
#else
		pthread_mutex_lock( &g_Mutex );
#endif

		// Set the position of the cursor
		SetCursorPosition();

		// Update the value of the %
		int p ;
		if ( d->m_CounterMaxVal==0 ) p=0;
		else p = (int) ( d->m_CurrentCounter * 100.0 / ((double) d->m_CounterMaxVal));

		char szBuffer[128];
		if ( d->m_NbThreads>1 )	sprintf(szBuffer,"#%02d: %03d", d->m_ThreadID, p);
		else          sprintf(szBuffer,"%03d", p);
		cout<< d->m_LinePrefix.c_str() << szBuffer <<"%  "; 

		// If the task is finished for that thread, show computation time
		if ( d->m_CurrentCounter==d->m_CounterMaxVal ) 
		{
			d->m_TimeProbe.Stop();
			int t = (int)d->m_TimeProbe.GetTotal();
			string strTime = ConvertTimeToString(t);
			cout<< "("<<strTime <<")                  ";
			g_NbFinished++;
		}
		// Show the estimated left time if needed and if possible (p>=1)
		else if ( d->m_ShowRemainingTime && p>=1 )
		{
			d->m_TimeProbe.Stop();
			int t = (int)(100.0*d->m_TimeProbe.GetTotal()/((double)p)) - (int)d->m_TimeProbe.GetTotal();
			d->m_TimeProbe.Start();

			string strTime = ConvertTimeToString(t);
			cout<< "(Estimated left: "<<strTime.c_str()<<")";

		}


		flush(cout);

#ifdef WIN32
		cout<<endl;
#else
		// If it was the last thread, put cursor at the end
		if ( g_NbFinished == (int)d->m_NbThreads )
			cout<<cursor_down(d->m_NbThreads-g_LastProcToCout+1)<<endl;
#endif

#ifdef WIN32
		ReleaseMutex(g_Mutex);
#else
		pthread_mutex_unlock( &g_Mutex );
#endif
	}

	/**********************************************************************************************//**
	 * \fn	std::string GetTotalTime() const
	 *
	 * \brief	Gets the total time. 
	 *
	 * \author	Benoit Scherrer
	 * \date	April 2010
	 *
	 * \return	The total time. 
	*************************************************************************************************/
	std::string TaskProgress::ConvertTimeToString(itk::RealTimeClock::TimeStampType t) 
	{
		char szBuffer[512];
	
		int ss = ((int)t)%60; 
		int mm = (int)((t-ss)/60)%60;
		int hh = (int)(t / 3600);
		sprintf(szBuffer, "%02dh %02dm %02ds", hh, mm, ss );

		return string(szBuffer);
	}

	/**********************************************************************************************//**
	 * \fn	void TaskProgress::SetLinePrefix(const std::string& prefix)
	 *
	 * \brief	Sets a line prefix (used for cout). 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \param	prefix	The prefix. 
	*************************************************************************************************/
	void TaskProgress::SetLinePrefix(const std::string& prefix)
	{
		d->m_LinePrefix = prefix;
	}

	/**********************************************************************************************//**
	 * \fn	std::string TaskProgress::GetLinePrefix()
	 *
	 * \brief	Gets the line prefix. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \return	The line prefix. 
	*************************************************************************************************/
	std::string TaskProgress::GetLinePrefix() const
	{
		return d->m_LinePrefix;
	}

	/**********************************************************************************************//**
	 * \fn	void TaskProgress::SetCursorPosition()
	 *
	 * \brief	Sets the cursor position. 
	 * \warning This function must be called inside the mutex protection
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	*************************************************************************************************/
	void TaskProgress::SetCursorPosition()
	{
		//----------------------------------------------------------
		// If it is the first time this function is called (for all threads)
		// Init all the lines
		//----------------------------------------------------------
		if ( g_LastProcToCout==-1 )
		{
#ifdef WIN32
			//-------------------------------------
			// If we are under win32, we can use 
			//functions to get/set the cursor pos
			//-------------------------------------
			CONSOLE_SCREEN_BUFFER_INFO csbi;
			if ( GetConsoleScreenBufferInfo( GetStdHandle(STD_OUTPUT_HANDLE), &csbi))
				g_InitYLine = csbi.dwCursorPosition.Y;
			else
				g_InitYLine = 0;
#endif
			//-------------------------------------
			// Init all the lines
			//-------------------------------------
			for ( unsigned int i=0; i<(d->m_NbThreads-1) ; i++ ) cout<<" "<<endl;
			cout<<" ";
			g_LastProcToCout = d->m_NbThreads-1;
		}

		//-------------------------------------
		// Under win32: set the cursor position
		//-------------------------------------
#ifdef WIN32
		COORD pos;
		pos.X=0;
		pos.Y=g_InitYLine + d->m_ThreadID;
		SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), pos);
#else
		//-------------------------------------
		// Under linux: move the cursor (relative move)
		//-------------------------------------

		// Move the vertical position of the cursor
		int voffset = (d->m_ThreadID-g_LastProcToCout);
		if ( voffset<0 ) cout<<cursor_up(-voffset);
		else if ( voffset>0 ) cout<<cursor_down(voffset);

		// Hz position: doesn t really matter the exact number of chars...
		cout<<cursor_left(d->m_LinePrefix.size()+45);
#endif

		g_LastProcToCout = d->m_ThreadID;
	}

	/**********************************************************************************************//**
	 * \fn	void TaskProgress::SetShowEstimatedRemainingTime(bool b)
	 *
	 * \brief	Sets if it should show the estimated remaining time. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \param	b	true to. 
	*************************************************************************************************/
	void TaskProgress::SetShowEstimatedRemainingTime(bool b)
	{
		d->m_ShowRemainingTime = b;
	}

	/**********************************************************************************************//**
	 * \fn	bool TaskProgress::GetShowEstimatedRemainingTime()
	 *
	 * \brief	Gets if it shows the estimated remaining time. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \return	true if it succeeds, false if it fails. 
	*************************************************************************************************/
	bool TaskProgress::GetShowEstimatedRemainingTime() const
	{
		return d->m_ShowRemainingTime;
	}







	ItkFilterProgressObserver::ItkFilterProgressObserver():
	m_TaskProgress(0,1,100,true, false) 
	{ }

	void ItkFilterProgressObserver::Execute(itk::Object *caller, const itk::EventObject & event)
	{
		Execute( (const itk::Object *)caller, event);
	}

	void ItkFilterProgressObserver::Execute(const itk::Object * object, const itk::EventObject & event)
	{
		ProcessPointer filter = dynamic_cast< ProcessPointer >( object );
		if( typeid( event ) == typeid( itk::ProgressEvent ) )
			m_TaskProgress.SetCurrentCounter((unsigned int)(filter->GetProgress()*100));

		if( typeid( event ) == typeid( itk::StartEvent ) )
		{
			TaskProgress::InitTaskProgress();
			m_TaskProgress.SetCurrentCounter(0);
		}
		//			if( typeid( event ) == typeid( itk::EndEvent ) )
	}

	void ItkFilterProgressObserver::SetProgressLinePrefix(const std::string& str)
	{
		m_TaskProgress.SetLinePrefix(str);
	}
}
