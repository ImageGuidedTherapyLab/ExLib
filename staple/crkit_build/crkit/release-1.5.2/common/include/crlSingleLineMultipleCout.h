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

#ifndef CRL_SINGLE_LINE_MULTIPLE_COUT_H
#define CRL_SINGLE_LINE_MULTIPLE_COUT_H

#include <iostream>
#include <cstdio>

#ifdef WIN32
  #include <windows.h>
  #pragma warning( disable: 4996 )
#endif

namespace crl{

	class SingleLineMultipleCout
	{
	public:
		SingleLineMultipleCout()
		{
			ResetToCurrentLine();
		}

		void ResetToCurrentLine()
		{
#ifdef WIN32
			//-------------------------------------
			// If we are under win32, we can use 
			//functions to get/set the cursor pos
			//-------------------------------------
			CONSOLE_SCREEN_BUFFER_INFO csbi;
			if ( GetConsoleScreenBufferInfo( GetStdHandle(STD_OUTPUT_HANDLE), &csbi))
				m_LastLineIndex = csbi.dwCursorPosition.Y;
			else
				m_LastLineIndex = 0;
#else
			std::cout<<std::endl;
#endif
		}

		void cout(const std::string& str)
		{
			cout((const char*)str.c_str());
		}

		void cout(const char* str)
		{
			//-------------------------------------
			// Under window: set the position
			//-------------------------------------
#ifdef WIN32
			COORD pos;
			pos.X=0;
			pos.Y=m_LastLineIndex;
			SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), pos);

			//-------------------------------------
			// Under linux: move the cursor (relative move)
			//-------------------------------------
#else
			// Move the vertical position of the cursor -1
			std::cout<<"\33[" << 1 << "A";		// up 1
			std::cout<<"\33[" << 45 << "D";	//left 45
#endif
			std::cout<< str << std::endl;
		}

	protected:
#ifdef WIN32
		int		m_LastLineIndex;
#endif
	};


} // end namespace crl.


#endif
