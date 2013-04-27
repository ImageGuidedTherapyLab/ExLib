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

#ifndef CRL_DWI_CONVERT_FSL_TO_NHDR
#define CRL_DWI_CONVERT_FSL_TO_NHDR

#include <iostream>
#include <stdlib.h>

#include <itkLightObject.h>
#include <itkImageIOBase.h>

#include "crlGeneralizedInstantiatorMacros.h"
#include "crlVectorImageUtils.h"
#include "crlNHDRFileUtils.h"
#include "crlDWICommon.h"

using namespace crl::DWI;

namespace crl {

/**********************************************************************************************//**
 * \class	crlDWIConvertFSLToNHDR_Base
 *
 * \brief	Base class for the instanciator. 
 *
 * \author	Benoit Scherrer
 * \date	January 2011
*************************************************************************************************/
class ITK_EXPORT crlDWIConvertFSLToNHDR_Base : public itk::LightObject
{
	CRL_INSTANCIATOR_BASE(crlDWIConvertFSLToNHDR_Base);

public:
	// Constructor/Destructor
	crlDWIConvertFSLToNHDR_Base():
	m_B0max(5)
	{}

	~crlDWIConvertFSLToNHDR_Base() {}

	// Algo parameters
	void SetInputDataFileName (const std::string& fileName ) { m_InputDataFileName=fileName; }
	void SetInputBVecsFileName (const std::string& fileName ) { m_InputBVecsFileName=fileName; }
	void SetInputBValsFileName (const std::string& fileName ) { m_InputBValsFileName=fileName; }
	void SetOutputFileName (const std::string& fileName ) { m_OutputFileName=fileName; }
	void SetB0max(float b) {m_B0max=b; }

	// Function to implement 
	virtual int Execute() = 0;

	/**********************************************************************************************//**
	 * \fn	void string2vector(const std::string& str, vector<double> &outputVals )
	 *
	 * \brief	Converts a space-separated list of numbers in a string to a vector<double> 
	 *
	 * \param	str					The string. 
	 * \param [in,out]	outputVals	The output vals. 
	*************************************************************************************************/
	void string2vector(const std::string& str, vector<double> &outputVals )
	{
		stringstream ss(str);
		string val;
		outputVals.clear();

		getline(ss,val,' ');
		while ( ss.rdstate() == ios::goodbit )
		{
			outputVals.push_back( atof(val.c_str()) );
			getline(ss,val,' ');
		}
		if ( val!="" && val!=" " )
			outputVals.push_back( atof(val.c_str()) );		//don't forget the last one
	}


protected:
	std::string		m_InputDataFileName;
	std::string		m_InputBVecsFileName;
	std::string		m_InputBValsFileName;
	std::string		m_OutputFileName;
	float			m_B0max;
};

/**********************************************************************************************//**
 * \class	crlDWIConvertFSLToNHDR_Algo
 *
 * \brief	Main class 
 *
 * \author	Benoit Scherrer
 * \date	January 2011
*************************************************************************************************/
template <typename TPixel, unsigned int VImageDimension  >
class ITK_EXPORT crlDWIConvertFSLToNHDR_Algo : public crlDWIConvertFSLToNHDR_Base
{
	CRL_INSTANCIATOR_CLASS(crlDWIConvertFSLToNHDR_Base, crlDWIConvertFSLToNHDR_Algo);

public:
	crlDWIConvertFSLToNHDR_Algo() {}
	~crlDWIConvertFSLToNHDR_Algo() {}

	virtual int Execute()
	{
		//-------------------------------------
		//  Read the gradient vectors
		//-------------------------------------
		cout<<"- Read gradient vectors from <"<<m_InputBVecsFileName.c_str()<<">..."<<endl;
		std::ifstream ifile_bvecs;
		ifile_bvecs.open ( m_InputBVecsFileName.c_str(), ifstream::in);
		if  ( !ifile_bvecs.is_open() ) 
		{
			crl::FileName fname(this->m_InputDataFileName);
			fname.setExtension("bvec");
			ifile_bvecs.open (fname.getCompleteFilePath().c_str(), ifstream::in);
			if  ( !ifile_bvecs.is_open() ) 
				throw itk::ExceptionObject(__FILE__,__LINE__,"Cannot open the gradient vector file for reading", "");
		}
		std::string sgX, sgY, sgZ;
		std::getline(ifile_bvecs, sgX);
		std::getline(ifile_bvecs, sgY);
		std::getline(ifile_bvecs, sgZ);

		std::vector<double> vgX, vgY, vgZ;
		string2vector(sgX, vgX );
		string2vector(sgY, vgY );
		string2vector(sgZ, vgZ );
		ifile_bvecs.close();

		// Check the different sizes
		if ( vgX.size()!=vgY.size() ||  vgX.size()!=vgZ.size() || vgY.size()!=vgZ.size())
			throw itk::ExceptionObject(__FILE__,__LINE__,"Error while reading the gradient vectors. The X,Y and Z components doesn't have the same size.", "");

		// Construct a std::vector<GradientDirectionType>
		std::vector<GradientDirectionType> bvecs;
		for ( unsigned int i=0; i<vgX.size() ; i++ )
		{	
			GradientDirectionType v;
			v[0] = vgX[i];
			v[1] = vgY[i];
			v[2] = vgZ[i];
			bvecs.push_back(v);
		}

		//-------------------------------------
		//  Read the B-Values
		//-------------------------------------
		cout<<"- Read B-Values from <"<<m_InputBValsFileName.c_str()<<">..."<<endl;
		std::ifstream ifile_bvals;
		ifile_bvals.open ( m_InputBValsFileName.c_str(), ifstream::in);
		if  ( !ifile_bvals.is_open() ) 
		{
			crl::FileName fname(this->m_InputDataFileName);
			fname.setExtension("bval");
			ifile_bvals.open (fname.getCompleteFilePath().c_str(), ifstream::in);
			if  ( !ifile_bvals.is_open() ) 
				throw itk::ExceptionObject(__FILE__,__LINE__,"Cannot open the b-values file for reading", "");
		}

		std::string sgB;
		std::getline(ifile_bvals, sgB);
		std::vector<double> vgB;
		string2vector(sgB, vgB );

		ifile_bvals.close();

		// Check the size
		if ( vgB.size()!=vgX.size()  )
			throw itk::ExceptionObject(__FILE__,__LINE__,"Error while reading the B-values. The number of b-values is not equal to the number of gradient directions.", "");

		//-------------------------------------------------
		// Threshold all small b-values
		//-------------------------------------------------
		if ( m_B0max!=0 ) std::cout<<"  Be Aware: B-values lower than "<<m_B0max<<"s/mm2 are considered to be B=0s/mm2"<<std::endl;
		GradientDirectionType vNull;
		vNull.fill(0);
		for ( unsigned int i=0; i<vgB.size(); i++ )
		{
			if ( vgB[i]<=m_B0max )
				bvecs[i] = vNull;
		}

		//-------------------------------------------------
		// Open the FSL 4D data file
		//-------------------------------------------------
		typedef typename itk::Image<TPixel, VImageDimension>	FSLImageType;
		typedef typename itk::ImageFileReader<FSLImageType> ReaderType;
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( m_InputDataFileName );
		reader->Update();

		// Check the size
		if ( vgB.size()!=reader->GetOutput()->GetLargestPossibleRegion().GetSize()[3]  )
			throw itk::ExceptionObject(__FILE__,__LINE__,"Error. The number of B-values/gradient vectors is not equal to the number of images in the input file.", "");


		/*-------------------------------------------------
		 Create the VectorImage from the FSL 4D file
		/-------------------------------------------------*/
		cout<<"- Convert N-D FSL data to VectorImage data..."<<endl;
		typedef typename itk::VectorImage< TPixel, (VImageDimension-1) > OutputImageType;
		typename OutputImageType::Pointer img = ConvertNDToVectorImage<TPixel, VImageDimension>( reader->GetOutput()  );
	
		/*-------------------------------------------------
		 Export the VectorImage with the correct b-value/gradient vectors
		/-------------------------------------------------*/
		cout<<"- Export to NHDR in <"<<m_OutputFileName.c_str()<<">..."<<endl;
		NHDR_WriteFromVectorImage< TPixel, (VImageDimension-1) >( img, bvecs, vgB, m_OutputFileName );
		return 0;
	}
};



}

#endif
