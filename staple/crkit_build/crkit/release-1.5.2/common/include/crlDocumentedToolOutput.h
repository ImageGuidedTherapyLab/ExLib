/**********************************************************************************************//**
 * \file	crlDocumentedToolOutput.h
 *
 * \brief	Declares the crl::DocumentedToolOutput class. 
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

#include <tclap/StdOutput.h>
#include <iostream>
#include <string>
#include "crlFileName.h"

namespace crl {

/**********************************************************************************************//**
 * \class	DocumentedToolOutput
 *
 * \brief	The default DocumentedToolOutput class for CRL tools.
 * 			
 * 			Example of use (from crlMaskImages):  
 * 	\code
 * try {
 * 	  TCLAP::CmdLine cmd("", ' ', CRKIT_VERSION_STRING );  
 *  
 * 	  crl::DocumentedToolOutput my("Mask an input image with an other and output the result.\n"
 * 						"Works with scalar images, vector images and diffusion vector images.",  
 * 	 					"Benoit	Scherrer", "2010" );  
 * 	  my.SetKeepOriginalDataType(true);
 * 	  cmd.setOutput( &my );
 *  
 * 	  TCLAP::ValueArg<string> argOutput  ("o","output",	"Output file (default: [inputfile]_masked)", false, "", "file", cmd);  
 * 	  TCLAP::ValueArg<string> argMask    ("m","mask",		"Mask file", true, "", "file", cmd);  
 *    TCLAP::ValueArg<string> argInput   ("i","input",	"Input file to mask", true, "", "file", cmd);  
 * }
 * // If parsing error
 * catch ( TCLAP::ArgException& e)
 * {
 *    std::cerr << "Error: "<<e.error()<< " for argument " <<e.argId() << std::endl;
 *    exit(1);
 * }
 * (etc...)
 * \endcode 
 *
 * \author	Benoit Scherrer.
 * \date	April 2010.
*************************************************************************************************/
class DocumentedToolOutput : public TCLAP::StdOutput
{
public:

	/**********************************************************************************************//**
	 *
	 * \brief	Constructor. 
	 *
	 * \author	Benoit Scherrer
	 * \date	April 2010
	 *
	 * \param	shortDescription	Information describing the tool (what's the purpose, etc...). 
	 * \param	authorName			Name of the author. 
	 * \param	year				The year (for the copyright). 
	*************************************************************************************************/
	DocumentedToolOutput(const std::string& shortDescription="", const std::string& authorName="", const std::string& year=""):
	  m_ShortDescription(shortDescription),
	  m_AuthorName(authorName),
	  m_Year(year),
	  m_KeepOriginalDataType(false)
	{
	}

    /**********************************************************************************************//**
     * \fn	void SetKeepOriginalDataType(bool b)
     *
     * \brief	Sets if that tool keeps the original data type or not. If yes, a special comment will
     * 			be added in the detailed usage of the tool. 
     *
     * \author	Benoit Scherrer
     * \date	April 2010
     *
     * \param	b	true to show the message. 
    *************************************************************************************************/
    void SetKeepOriginalDataType(bool b)
	{
		m_KeepOriginalDataType = b;
	}

	void showProgNameAndVersion(TCLAP::CmdLineInterface& _cmd, std::ostream& c)
	{
		crl::FileName prgm(_cmd.getProgramName());
		std::string version = _cmd.getVersion();
		if ( version!="" ) version=" (CRKIT v"+version+")";

		c << "=========================================================================="<<std::endl;
		c << prgm.getFileName() << version <<std::endl;
		if ( m_ShortDescription!="" ) c << std::endl << m_ShortDescription << std::endl;

		c << std::endl;
		if ( m_Year!="" ) c << "Copyright (c) "<<m_Year<<", ";
		else c << "Copyright (c) ";
		c << "Children's Hospital Boston."<<std::endl;
		c << "    Computational Radiology Laboratory (CRL)."<<std::endl;
		if ( m_AuthorName!="") c<<"    ("<<m_AuthorName<<")"<<std::endl;
		c << "=========================================================================="<<std::endl;

	}

	virtual void failure(TCLAP::CmdLineInterface& _cmd, TCLAP::ArgException& e)
	{ 
		showProgNameAndVersion(_cmd, std::cerr);

		std::string progName = _cmd.getProgramName();

		std::cerr << "PARSE ERROR: " << e.argId() << std::endl
			<< "             " << e.error() << std::endl << std::endl;

		if ( _cmd.hasHelpAndVersion() )
		{
			std::cerr << "Brief USAGE: " << std::endl;

			_shortUsage( _cmd, std::cerr );	

			std::cerr << std::endl << "For complete USAGE and HELP type: " 
				<< std::endl << "   " << progName << " --help" 
				<< std::endl << std::endl;
		}
		else
			usage(_cmd);
	}

	virtual void usage(TCLAP::CmdLineInterface& _cmd)
	{
		showProgNameAndVersion(_cmd, std::cerr);

		std::cout << std::endl << "USAGE: " << std::endl << std::endl; 

		_shortUsage( _cmd, std::cout );

		std::cout << std::endl << std::endl << "Where: " << std::endl << std::endl;

		_longUsage( _cmd, std::cout );

		if ( m_KeepOriginalDataType )
		{
			std::cout << "REMARK: This tool keeps the original data type."<<std::endl; 
		}

		//cout << "my usage message:" << endl;
		//list<TCLAP::Arg*> args = c.getArgList();
		//for (TCLAP::ArgListIterator it = args.begin(); it != args.end(); it++)
		//	cout << (*it)->longID() 
		//		 << "  (" << (*it)->getDescription() << ")" << endl;
	}

	virtual void version(TCLAP::CmdLineInterface& _cmd)
	{
		std::string progName = _cmd.getProgramName();
		std::string version = _cmd.getVersion();

		std::cout << std::endl << progName << "  version: " 
			<< version << std::endl << std::endl;
	}

protected:
	std::string m_ShortDescription;
	std::string m_AuthorName;
	std::string m_Year;

	bool	m_KeepOriginalDataType;
};

}


