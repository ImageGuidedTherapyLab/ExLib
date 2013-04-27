/**********************************************************************************************//**
 * \file	crlFileName.h
 *
 * \brief	Declares the crl::FileName class. 
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


#ifndef h_FILE_NAME_UTIL
#define h_FILE_NAME_UTIL

#include <string>

namespace crl {

/**********************************************************************************************//**
 * \class	FileName
 *
 * \brief	Class to manipulate filenames. You can set or
 * 			extract the file extension, the filename or the path. You can also construct a new
 * 			file name with a sufix or a prefix. 
 *
 * <b>Examples:</b>
 * \code
 * crl::FileName fn("/usr/local/test.txt");		
 *
 * fn.getExtension();                         // returns "txt"
 * fn.getPath();                              // returns "/usr/local/"
 * fn.getFileName();                          // returns "test"
 * fn.getFileNameWithExt();                   // returns "test.txt"
 * fn.getCompleteFilePath_WithPrefix("a_");   // returns "/usr/local/a_test.txt"
 * fn.getCompleteFilePath_WithSufix("_01");   // returns "/usr/local/test_01.txt"
 * \endcode
 *
 * \remark The file extension is defined by the string following the last '.' 
 * and is of maximum 3 characters. However, for 'gz' files, getExtension()
 * looks if there is an other valid extension before 'gz'. 
 * \code
 * crl::FileName fn("/usr/local/test.nii.gz");	
 * fn.getExtension();                         // returns "nii.gz"
 * \endcode
 *
 * \author	Benoit Scherrer
 * \date	December 2009
*************************************************************************************************/
class FileName
{
public:
	FileName();
	FileName(const std::string& strPathAndFileName );
	FileName(char *strPathAndFileName);
	~FileName();

	void	setCompleteFilePath(const std::string& filename);
	void	setExtension(const std::string& ext); 
	void	setPath(const std::string& path) ;
	void	setFileName( const std::string& filename);

	void	addFileNamePrefix( const std::string& prefix );
	void	addFileNameSuffix( const std::string& suffix );


	std::string	getPath(void) const ;
	std::string	getFileName(void) const ;
	std::string	getExtension(void) const ;
	std::string	getFileNameWithExt(void) const ;

	std::string	getCompleteFilePath(void) const; 
	std::string	getCompleteFilePath_WithPrefix(const std::string& prefix) const;
	std::string	getCompleteFilePath_WithSufix(const std::string& sufix) const;

	std::string getWithNewFilename(const std::string& newFileName) const;
	std::string getWithNewFilenameAndExt(const std::string& newFileName) const;
	std::string getWithNewPath(const std::string& newPath) const;
	std::string getWithNewExtension(const std::string& newExt) const;

private:
	class		PrivateData;
	PrivateData *d;
};

}

#endif
