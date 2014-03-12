/*
 * Copyright (c) 2008-2011 Children's Hospital Boston.
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

// Remove annoying warning with Visual Studio
#ifdef _MSC_VER
  #pragma warning( disable: 4996 )
#endif

#include "crlDWIRaw.h"


#include "itkImageFileReader.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileReader.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>


namespace crl {
namespace DWI {

template <class TPixel> 
DWIRaw<TPixel>::DWIRaw( void ):
m_Verbose(true),
m_NumberOfImages(0),
m_NumberOfGradientImages(0),
m_ReadB0(false),
m_B0(0)
{

}

template <class TPixel> 
DWIRaw<TPixel>::~DWIRaw()
{
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> void DWIRaw<TPixel>::LoadData(const string& filename)
 *
 * \brief	Loads a data file. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \exception	Thrown when error. 
 *
 * \typeparam	TPixel	. 
 * \param	filename	Filename of the file. 
*************************************************************************************************/
template <class TPixel> 
void DWIRaw<TPixel>::LoadData(const std::string& filename)
{
	/*---------------------------------
	 Open the nhrd file with ITK
	 ---------------------------------*/
	typename itk::ImageFileReader<RawImageType>::Pointer reader = itk::ImageFileReader<RawImageType>::New();

	// Set the properties for NrrdReader
	reader->SetFileName(filename);

	// Read in the nrrd data. The file contains the reference image and the gradient
	// images. (could raise an exception)
	reader->Update();
	m_Data = reader->GetOutput();

	//-------------------------------------
	// Open input header file to read the list of file names
	//-------------------------------------
	m_FileNames.clear();
	std::ifstream ifile;
	ifile.open ( filename.c_str(), std::ofstream::in);
	if  ( !ifile.is_open() ) 
		throw itk::ExceptionObject(__FILE__,__LINE__,"Cannot open the input header file for reading", "");

	// Look for the string "data file: LIST"
	char szBuffer[4096];
	bool next=true;
	while ( (!ifile.eof()) && next)
	{
		ifile.getline(szBuffer, 4095);
		std::string line = std::string(szBuffer);
		if ( line.find("data file: LIST")!=std::string::npos )
			next=false;
	}

	// If there are the filenames in the NHDR
	if  ( !next ) 
	{
		// Now read the filenames
		while ( (!ifile.eof()) && m_FileNames.size()< m_Data->GetNumberOfComponentsPerPixel() )
		{
			ifile.getline(szBuffer, 4095);
			std::string line = std::string(szBuffer);
			if ( line!="" )
				m_FileNames.push_back(line);
		}

		ifile.close();

		if ( m_FileNames.size()!=m_Data->GetNumberOfComponentsPerPixel() )
			throw itk::ExceptionObject(__FILE__,__LINE__,"Error. Invalid number of files ", "");
	}
	else
		m_FileNames.clear();

	// -------------------------------------------------------------------------
	// Parse the Nrrd headers to get the B value and the gradient directions used
	// for diffusion weighting. 
	// 
	// The Nrrd headers should look like :
	// The tags specify the B value and the gradient directions. If gradient 
	// directions are (0,0,0), it indicates that it is a reference image. 
	//
	// DWMRI_b-value:=800
	// DWMRI_gradient_0000:= 0 0 0
	// DWMRI_gradient_0001:=-1.000000       0.000000        0.000000
	// DWMRI_gradient_0002:=-0.166000       0.986000        0.000000
	// DWMRI_gradient_0003:=0.110000        0.664000        0.740000
	// ...
	// 
	itk::MetaDataDictionary imgMetaDictionary = m_Data->GetMetaDataDictionary();    
	std::vector<std::string> imgMetaKeys = imgMetaDictionary.GetKeys();
	std::vector<std::string>::const_iterator itKey = imgMetaKeys.begin();
	std::string metaString;

	m_NumberOfImages = 0;
	m_NumberOfGradientImages = 0;
	m_ReadB0 = false;
	m_B0 = 0;

	GradientDirectionType vect3d;
	m_DiffusionVectors = GradientDirectionContainerType::New();


	for (; itKey != imgMetaKeys.end(); itKey ++)
	{
		double x,y,z;

		itk::ExposeMetaData<std::string> (imgMetaDictionary, *itKey, metaString);
		if (itKey->find("DWMRI_gradient") != std::string::npos)
		{ 
			if ( m_Verbose ) std::cout << *itKey << " ---> " << metaString << std::endl;      
			sscanf(metaString.c_str(), "%lf %lf %lf\n", &x, &y, &z);
			vect3d[0] = x; vect3d[1] = y; vect3d[2] = z;
			m_DiffusionVectors->InsertElement( m_NumberOfImages, vect3d );
			++m_NumberOfImages;

			// If the direction is 0.0, this is a reference image
			if (vect3d[0] == 0.0 &&
				vect3d[1] == 0.0 &&
				vect3d[2] == 0.0)
			{
				continue;
			}
			++m_NumberOfGradientImages;;
		}
		else if (itKey->find("DWMRI_b-value") != std::string::npos)
		{
			if ( m_Verbose ) std::cout << *itKey << " ---> " << metaString << std::endl;      
			m_ReadB0 = true;
			m_B0 = atof(metaString.c_str());
		}
	}
	if ( m_Verbose ) 
	{
		std::cout << "Number of gradient images: "
		<< m_NumberOfGradientImages
		<< " and Number of reference images: "
		<< m_NumberOfImages - m_NumberOfGradientImages
		<< std::endl;
	}

	if(!m_ReadB0)
		throw itk::ExceptionObject(__FILE__,__LINE__,"BValue not specified in header file", "");

}

/**********************************************************************************************//**
 * \fn	template <class TPixel> typename DWIRaw<TPixel>::GradientDirectionContainerType* DWIRaw<TPixel>::GetDiffusionVectors() const
 *
 * \brief	Gets the diffusion vectors. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \typeparam	TPixel The pixel type for the gradient images. 
 *
 * \return	null if it fails, else the diffusion vectors. 
*************************************************************************************************/
template <class TPixel> 
typename DWIRaw<TPixel>::GradientDirectionContainerType* DWIRaw<TPixel>::GetDiffusionVectors() const 
{ 
	return m_DiffusionVectors.GetPointer(); 
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> typename DWIRaw<TPixel>::RawImageType* DWIRaw<TPixel>::GetData() const
 *
 * \brief	Gets the raw data as a VectorImage. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \typeparam	TPixel	. 
 *
 * \return	null if it fails, else the data. 
*************************************************************************************************/
template <class TPixel> 
typename DWIRaw<TPixel>::RawImageType* DWIRaw<TPixel>::GetData() const 
{ 
	return m_Data.GetPointer(); 
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> double DWIRaw<TPixel>::GetB0() const
 *
 * \brief	Gets the maximum B-value. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \typeparam	TPixel	. 
 *
 * \return	The b 0. 
*************************************************************************************************/
template <class TPixel> 
double DWIRaw<TPixel>::GetB0() const  
{ 
	return m_B0; 
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> const vector<string>& DWIRaw<TPixel>::GetFileNames() const
 *
 * \brief	Gets the file names associed with a NHDR (the return vector can be empty). 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \typeparam	TPixel	. 
 *
 * \return	The file names. 
*************************************************************************************************/
template <class TPixel> 
const std::vector<std::string>& DWIRaw<TPixel>::GetFileNames() const 
{ 
	return m_FileNames; 
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> int DWIRaw<TPixel>::GetNumberOfB0Images() const
 *
 * \brief	Gets the number of B=0 images. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \typeparam	TPixel	. 
 *
 * \return	The number of b 0 images. 
*************************************************************************************************/
template <class TPixel> 
int DWIRaw<TPixel>::GetNumberOfB0Images()  const
{ 
	return (m_NumberOfImages -m_NumberOfGradientImages); 
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> int DWIRaw<TPixel>::GetNumberOfGradientImages() const
 *
 * \brief	Gets the total number of gradient images. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \typeparam	TPixel	. 
 *
 * \return	The number of gradient images. 
*************************************************************************************************/
template <class TPixel> 
int DWIRaw<TPixel>::GetNumberOfGradientImages() const
{ 
	return m_NumberOfGradientImages; 
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> void DWIRaw<TPixel>::SetVerbose( bool v )
 *
 * \brief	Sets the verbose mode. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \typeparam	TPixel	. 
 * \param	v	true to v. 
*************************************************************************************************/
template <class TPixel> 
void DWIRaw<TPixel>::SetVerbose( bool v ) 
{ 
	m_Verbose = v; 
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> bool DWIRaw<TPixel>::GetVerbose() const
 *
 * \brief	Gets the verbose mode. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \typeparam	TPixel	. 
 *
 * \return	true if it succeeds, false if it fails. 
*************************************************************************************************/
template <class TPixel> 
bool DWIRaw<TPixel>::GetVerbose() const
{ 
	return m_Verbose; 
}


} // end namespace DWI.
} // end namespace crl.

//***********************************************************************
