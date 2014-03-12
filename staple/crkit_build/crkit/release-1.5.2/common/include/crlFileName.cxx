/**********************************************************************************************//**
 * \file	crlFileName.cxx
 *
 * \brief	The crl::FileName class. 
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


#include "crlFileName.h"

using namespace std;

namespace crl {

/**********************************************************************************************//**
 * \class	FileName::PrivateData
 *
 * \brief	The private data for FileName. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
*************************************************************************************************/
class FileName::PrivateData
{
public:
	PrivateData():
		strDirectory(""),
		strExtension(""),
		strFileName("")
		{}

	std::string	strDirectory;		// directory, including last '/' char
	std::string	strExtension;
	std::string	strFileName;
};

/**********************************************************************************************//**
 * \fn	FileName::FileName()
 *
 * \brief	Default constructor. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
*************************************************************************************************/
FileName::FileName()
{
	d = new crl::FileName::PrivateData();
}

/**********************************************************************************************//**
 * \fn	FileName::FileName(const std::string& strPathAndFileName)
 *
 * \brief	Constructor, with a full path and file name. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \param	strPathAndFileName	Filename. 
*************************************************************************************************/
FileName::FileName(const std::string& strPathAndFileName)
{
	d = new FileName::PrivateData();
	setCompleteFilePath(strPathAndFileName);
}

/**********************************************************************************************//**
 * \fn	FileName::FileName( char *strPathAndFileName)
 *
 * \brief	Constructor, with a full path and file name. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \param [in,out]	strPathAndFileName	If non-null, filename. 
*************************************************************************************************/
FileName::FileName( char *strPathAndFileName)
{
	d = new FileName::PrivateData();
	setCompleteFilePath(string(strPathAndFileName));
}

/**********************************************************************************************//**
 * \fn	FileName::~FileName()
 *
 * \brief	Destructor. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
*************************************************************************************************/
FileName::~FileName()
{
	delete d;
}

/**********************************************************************************************//**
 * \fn	void FileName::setCompleteFilePath(const string& filename)
 *
 * \brief	Sets the full filename (path+filename) 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \param	filename	Filename of the file. 
*************************************************************************************************/
void FileName::setCompleteFilePath(const string& filename)
{
	// File directory
	string::size_type p = filename.find_last_of('/');
	if ( p==string::npos ) p = filename.find_last_of('\\');		// for win32
	if ( p==string::npos )
	{
		d->strDirectory="";
		d->strFileName = filename;
	}
	else 
	{	
		d->strDirectory = filename.substr(0, p+1) ;
		d->strFileName = filename.substr(p+1);
	}

	// File extension
	d->strExtension = "";
	p = d->strFileName.find_last_of('.');
	if ( p!=string::npos ) 
	{
		d->strExtension = d->strFileName.substr(p+1);
		if ( d->strExtension.length()>4 ) d->strExtension="";
		else {
			d->strFileName = d->strFileName.substr(0, p );

			// If gz, check if there is an other extension before
			// (ex : nii.gz, tar.gz, ...)
			if ( d->strExtension=="gz" )
			{
				p = d->strFileName.find_last_of('.');
				if ( p!=string::npos )
				{
					if ( d->strFileName.substr(p+1).length()<=3 )
					{
						d->strExtension = d->strFileName.substr(p+1) + "." + d->strExtension;
						d->strFileName = d->strFileName.substr(0, p );
					}
				}
			}
		}
	}
}

/**********************************************************************************************//**
 * \fn	string FileName::getCompleteFilePath(void) const
 *
 * \brief	Gets the full filename (path+filename) 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \return	The path and file name. 
*************************************************************************************************/
string FileName::getCompleteFilePath(void) const
{
	if ( d->strExtension=="" )
		return d->strDirectory+d->strFileName; 
	else
		return d->strDirectory+d->strFileName+"."+d->strExtension; 
}

/**********************************************************************************************//**
 * \fn	string FileName::getCompleteFilePath_WithPrefix(const string& prefix) const
 *
 * \brief	Returns a new full filename string with an additional prefix to the filename. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \param	prefix	The prefix. 
 *
 * \return	The path and file name with prefix. 
*************************************************************************************************/
string FileName::getCompleteFilePath_WithPrefix(const string& prefix) const
{
	if ( d->strExtension=="" )
		return d->strDirectory+prefix+d->strFileName;
	else
		return d->strDirectory+prefix+d->strFileName+"."+d->strExtension;
}

/**********************************************************************************************//**
 * \fn	string FileName::getCompleteFilePath_WithSufix(const string& sufix) const
 *
 * \brief	Gets a path and file name with sufix. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \param	sufix	The sufix. 
 *
 * \return	The path and file name with sufix. 
*************************************************************************************************/
string FileName::getCompleteFilePath_WithSufix(const string& sufix) const
{
	if ( d->strExtension=="" )
		return d->strDirectory+d->strFileName+sufix;
	else
		return d->strDirectory+d->strFileName+sufix+"."+d->strExtension;
}

/**********************************************************************************************//**
 * \fn	std::string FileName::getWithNewFilename(const std::string& newFileName) const
 *
 * \brief	Gets the current complete file name with a new filename (and keeps the original
 * 			extension) 
 *
 * \author	Benoit Scherrer
 * \date	April 2010
 *
 * \param	newFileName	Filename of the new file. 
 *
 * \return	The with new filename. 
*************************************************************************************************/
std::string FileName::getWithNewFilename(const std::string& newFileName) const
{
	if ( d->strExtension=="" )
		return d->strDirectory + newFileName;
	else
		return d->strDirectory + newFileName + "." + d->strExtension;
}

/**********************************************************************************************//**
 *
 * \brief	Gets the current complete file name with a new filename (containing the extension)
 *
 * \author	Benoit Scherrer
 * \date	April 2010
 *
 * \param	newFileName	Filename of the new file. 
 *
 * \return	. 
*************************************************************************************************/
std::string FileName::getWithNewFilenameAndExt(const std::string& newFileName) const
{
	return d->strDirectory + newFileName;
}

/**********************************************************************************************//**
 * \fn	std::string FileName::getWithNewPath(const std::string& newPath) const
 *
 * \brief	Gets the current complete file name with a new path. 
 *
 * \author	Benoit Scherrer
 * \date	April 2010
 *
 * \param	newPath	New pathname. 
 *
 * \return	The with new path. 
*************************************************************************************************/
std::string FileName::getWithNewPath(const std::string& newPath) const
{
	string path = newPath;

	if ( path.at(path.length()-1)!='/' && path.at(path.length()-1)!='\\' )
		path = path+"/";

	if ( d->strExtension=="" )
		return path + d->strFileName;
	else
		return path + d->strFileName + "." + d->strExtension;

}

/**********************************************************************************************//**
 * \fn	std::string FileName::getWithNewExtension(const std::string& newExt) const
 *
 * \brief	Gets the current complete filename with a new extension.  
 *
 * \author	Benoit Scherrer
 * \date	April 2010
 *
 * \param	newExt	The new file extension, without the first '.'. For example 'nrrd', 'nii.gz', ...
 *
 * \return	The with new extension. 
*************************************************************************************************/
std::string FileName::getWithNewExtension(const std::string& newExt) const
{
	if ( newExt=="" )
		return d->strDirectory+d->strFileName;
	else
		return d->strDirectory+d->strFileName+"."+newExt;

}


/**********************************************************************************************//**
 * \fn	string FileName::getPath(void) const
 *
 * \brief	Gets the path, including the last '/'. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \return	The path. 
*************************************************************************************************/
string	FileName::getPath(void) const 
{ 
	return d->strDirectory; 
}

/**********************************************************************************************//**
 * \fn	string FileName::getExtension(void) const
 *
 * \brief	Gets the file extension (without the first '.'. for example: "nrrd"). 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \return	The extension. 
*************************************************************************************************/
string	FileName::getExtension(void) const 
{ 
	return d->strExtension; 
}

/**********************************************************************************************//**
 * \fn	string FileName::getFileName(void) const
 *
 * \brief	Gets the file name. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \return	The file name. 
*************************************************************************************************/
string	FileName::getFileName(void) const 
{ 
	return d->strFileName; 
}

/**********************************************************************************************//**
 * \fn	string FileName::getFileNameWithExt(void) const
 *
 * \brief	Gets the file name with the extension. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \return	The file name with extent. 
*************************************************************************************************/
string	FileName::getFileNameWithExt(void) const 
{ 
	if ( d->strExtension=="" )
		return d->strFileName;
	else
		return d->strFileName+"."+d->strExtension; 
}

/**********************************************************************************************//**
 * \fn	void FileName::setExtension(const string& ext)
 *
 * \brief	Sets the extension. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \param	ext	The extent. 
*************************************************************************************************/
void FileName::setExtension(const string& ext)
{ 
	if ( ext=="" ) d->strExtension="";
	else if ( ext.data()[0]=='.' ) d->strExtension = ext.substr(1);
	else d->strExtension = ext; 
}

/**********************************************************************************************//**
 * \fn	void FileName::setPath(const string& path)
 *
 * \brief	Sets the path (Add an additional '/' at the end if necessary). 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \param	path	Full pathname. 
*************************************************************************************************/
void FileName::setPath(const string& path)
{ 
	if ( path=="" ) d->strDirectory = "";
	else if ( path.data()[path.size()-1]!='/' && path.data()[path.size()-1]!='\\' )  d->strDirectory = path+"/";
	else d->strDirectory = path;
}

/**********************************************************************************************//**
 * \fn	void FileName::setFileName(const string& filename)
 *
 * \brief	Sets the file name. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \param	filename	Filename of the file. 
*************************************************************************************************/
void FileName::setFileName(const string& filename) 
{ 
	d->strFileName = filename; 
}

void FileName::addFileNamePrefix( const std::string& prefix )
{
	d->strFileName = prefix + d->strFileName;
}

void FileName::addFileNameSuffix( const std::string& suffix )
{
	d->strFileName = d->strFileName + suffix;
}


}
