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

#ifndef CRL_DWI_CONVERT_NHDR_FOR_FSL_H 
#define CRL_DWI_CONVERT_NHDR_FOR_FSL_H

#ifdef _MSC_VER
  #pragma warning( disable: 4996 )
#endif

#include <iostream>
#include <stdlib.h>

#include <itkLightObject.h>
#include <itkImageIOBase.h>
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include "crlDWICommon.h"
#include "crlGeneralizedInstantiatorMacros.h"
#include "crlVectorImageUtils.h"
#include "crlFileName.h"
#include "crlDWIRaw.h"


namespace crl {

/**********************************************************************************************//**
 * \class	crlDWIConvertNHDRForFSL_Base
 *
 * \brief	Base class for the generalized instantiator 
 *
 * \author	Benoit Scherrer
 * \date	January 2011
*************************************************************************************************/
class ITK_EXPORT crlDWIConvertNHDRForFSL_Base : public itk::LightObject
{
	CRL_INSTANCIATOR_BASE(crlDWIConvertNHDRForFSL_Base);

public:
	// Constructor/Destructor
	crlDWIConvertNHDRForFSL_Base():
	m_B0max(5)
	{}
	~crlDWIConvertNHDRForFSL_Base() {}

	// Algo parameters
	void SetInputFileName (const std::string& fileName ) { m_InputFileName=fileName; }
	void SetOutputFileName (const std::string& fileName ) { m_OutputFileName=fileName; }
	void SetOutputBVecsFileName (const std::string& fileName ) { m_OutputBVecsFileName=fileName; }
	void SetOutputBValsFileName (const std::string& fileName ) { m_OutputBValsFileName=fileName; }
	void SetOutputB0FileName (const std::string& fileName ) { m_OutputB0FileName=fileName; }
	void SetB0max(float b) {m_B0max=b; }

	// Function to implement 
	virtual int Execute() = 0;

protected:
	std::string		m_InputFileName;
	std::string		m_OutputFileName;
	std::string		m_OutputBVecsFileName;
	std::string		m_OutputBValsFileName;
	std::string		m_OutputB0FileName;
	float			m_B0max;
};

/**********************************************************************************************//**
 * \class	crlDWIConvertNHDRForFSL_Algo
 *
 * \brief	Main class
 *
 * \author	Benoit Scherrer
 * \date	January 2011
*************************************************************************************************/
template <typename TPixel, unsigned int VImageDimension  >
class ITK_EXPORT crlDWIConvertNHDRForFSL_Algo : public crlDWIConvertNHDRForFSL_Base
{
	CRL_INSTANCIATOR_CLASS(crlDWIConvertNHDRForFSL_Base, crlDWIConvertNHDRForFSL_Algo);

public:
	crlDWIConvertNHDRForFSL_Algo() {}
	~crlDWIConvertNHDRForFSL_Algo() {}

	virtual int Execute()
	{
		//-------------------------------------
		//  Read the diffusion data
		//-------------------------------------
		typedef typename crl::DWI::DWIRaw<TPixel>  DWIRawType;
		DWIRawType rawData;
		rawData.SetVerbose(false);
		cout<<"- Read <"<<m_InputFileName.c_str()<<">..."<<endl;
		rawData.LoadData(m_InputFileName);
	
		//-------------------------------------
		//  Open bvecs file for writing
		//-------------------------------------
		std::ofstream ofile_bvecs;
		ofile_bvecs.open ( m_OutputBVecsFileName.c_str(), ofstream::out);
		if  ( !ofile_bvecs.is_open() ) 
			throw itk::ExceptionObject(__FILE__,__LINE__,"Cannot open the bvecs file for writing", "");

		//-------------------------------------
		//  Open bvals file for writing
		//-------------------------------------
		std::ofstream ofile_bvals;
		ofile_bvals.open ( m_OutputBValsFileName.c_str(), ofstream::out);
		if  ( !ofile_bvals.is_open() ) 
			throw itk::ExceptionObject(__FILE__,__LINE__,"Cannot open the bvals file for writing", "");


		//-------------------------------------
		// Loop on all the gradients and write data in bvecs and bvals
		//-------------------------------------
		cout<<"- Export gradient vectors in <"<<m_OutputBVecsFileName.c_str()<<">..."<<endl;
		cout<<"- Export B-Values in <"<<m_OutputBValsFileName.c_str()<<">..."<<endl;
		if ( m_B0max!=0 ) std::cout<<"  Be Aware: B-values lower than "<<m_B0max<<"s/mm2 are considered to be B=0s/mm2"<<std::endl;

		typename DWIRawType::GradientDirectionContainerType* gradContainer = rawData.GetDiffusionVectors();
		string list_gX = "";
		string list_gY = "";
		string list_gZ = "";
		string list_B = "";
	
		bool b0already = false;
		int id=0;
		char szBuffer[128];
		for( typename DWIRawType::GradientDirectionContainerType::Iterator it = gradContainer->Begin(); it != gradContainer->End(); it++)
		{
			vnl_vector_fixed<double, 3> v = it.Value();

			if (list_gX != "" ) list_gX = list_gX + " ";
			if (list_gY != "" ) list_gY = list_gY + " ";
			if (list_gZ != "" ) list_gZ = list_gZ + " ";
			if (list_B != "" ) list_B = list_B + " ";

			double nn = v.two_norm();
			
			// Normalize the gradient vector
			vnl_vector_fixed<double, 3> vnn = v;
			if (nn!=0) vnn=v/nn;

			// Compute the B-value from the norm of the vector
			double B=rawData.GetB0() * nn * nn;

			if ( B<=m_B0max ) {
				B=0;
				vnn[0]=0;
				vnn[1]=0;
				vnn[2]=0;
			}

			sprintf(szBuffer, "%f", vnn[0]); list_gX += string(szBuffer);
			sprintf(szBuffer, "%f", vnn[1]); list_gY += string(szBuffer);
			sprintf(szBuffer, "%f", vnn[2]); list_gZ += string(szBuffer);

			sprintf(szBuffer, "%d", (int)B); list_B += string(szBuffer);

			// Write one B=0 image in a file to use brain extraction after
			if ( (m_OutputB0FileName!="") && (!b0already) && (nn==0) )
			{
				crl::WriteComponentOfVectorImage<TPixel, VImageDimension>( rawData.GetData(), id, m_OutputB0FileName );
				b0already = true;
			} 
			id++;
		}

		ofile_bvecs << list_gX.c_str() << endl;
		ofile_bvecs << list_gY.c_str() << endl;
		ofile_bvecs << list_gZ.c_str() << endl;

		ofile_bvals << list_B.c_str() << endl;

		ofile_bvecs.close();
		ofile_bvals.close(); 

		//-------------------------------------
		// Now convert all nrrd files into one 4D image for FSL
		//-------------------------------------
		cout<<"- Convert the VectorImage to a N-D image in <"<<m_OutputFileName.c_str()<<">..."<<endl;

		ConvertVectorImageToNDImage<TPixel, VImageDimension>( rawData.GetData(), m_OutputFileName  );
		return 0;
	}

};



}

#endif
