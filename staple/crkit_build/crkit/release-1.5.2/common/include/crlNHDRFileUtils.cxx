/**********************************************************************************************//**
 * \file	crlNHDRFileUtils.cxx
 *
 * \brief	Useful function to manipulate files 
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



#include "crlNHDRFileUtils.h"
#include "crlFileName.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <itkExceptionObject.h>
#include <itkMetaDataObject.h>
#include <itkAffineTransform.h>

using namespace std;

namespace crl {

	/**********************************************************************************************//**
	 * \fn	void NHDR_IsValid( const string& fileHdr, int &nbImages )
	 *
	 * \brief	Check the validity of a NHDR header and returns the number of images. 
	 *
	 * \author	Benoit Scherrer
	 * \date	December 2009
	 *
	 * \exception itk::ExceptionObject	Thrown when error. 
	 *
	 * \param	fileHdr				The file header. 
	 * \param [in,out]	nbImages	The number of images. 
	*************************************************************************************************/
	void NHDR_IsValid( const string& fileHdr, int &nbImages )
	{
		nbImages = -1;

		//-------------------------------------
		// Open input header file to read the list of file names
		//-------------------------------------
		std::ifstream ifile;
		ifile.open ( fileHdr.c_str(), ofstream::in);
		if  ( !ifile.is_open() ) 
			throw itk::ExceptionObject(__FILE__,__LINE__,"Cannot open the input header file for reading", "NHDR_IsValid");

		//-------------------------------------
		// Read the first line = NRRD signature
		//-------------------------------------
		char szBuffer[4096];

		ifile.getline(szBuffer, 4095);
		if ( szBuffer[0]!='N' || szBuffer[1]!='R' || szBuffer[2]!='R' || szBuffer[3]!='D' )
			throw itk::ExceptionObject(__FILE__,__LINE__,"Error. Invalid NRRD signature in the input reference header.", "NHDR_IsValid");

		//-------------------------------------
		// Read until string "data file: LIST"
		//-------------------------------------
		bool next=true;
		while ( (!ifile.eof()) && next)
		{
			ifile.getline(szBuffer, 4095);
			string line = string(szBuffer);
			if ( line.find("data file:")!=string::npos )
				next=false;

			int x, y, z, t;
			// Get the number of images
			if ( sscanf(szBuffer, "sizes: %d %d %d %d", &x, &y, &z, &t) ) 
				nbImages = t;
		}

		// If there are the filenames in the NHDR
		if  ( next ) 
			throw itk::ExceptionObject(__FILE__,__LINE__,"Error. Cannot find 'data file:' in the input header file", "NHDR_IsValid");

		// If there are the filenames in the NHDR
		if  ( nbImages==-1 ) 
			throw itk::ExceptionObject(__FILE__,__LINE__,"Error. Cannot find 'sizes:' in the input header file", "NHDR_IsValid");

	}

	/**********************************************************************************************//**
	 * \fn	void NHDR_CopyAndModifyHeader(const string& refHdr, const string& outputHdr,
	 * 		const std::vector<string>& filenames, const std::vector<vnl_vector_fixed< double,
	 * 		3 > >& gvectors )
	 *
	 * \brief	Copy an input NHDR header in a new NHDR header with different data files and gradient
	 * 			vectors. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \exception	itk::ExceptionObject	Thrown when exception. 
	 *
	 * \param	refHdr		The reference header. 
	 * \param	outputHdr	The output header. 
	 * \param	filenames	The filenames. 
	 * \param	gvectors	The new gradient vectors (already scaled to contain the b-value). If the
	 * 						size is null, then keep the original gradient vectors. 
	 *
	 * ### exception	itk	:ExceptionObject	Thrown when error. 
	*************************************************************************************************/
	void NHDR_CopyAndModifyHeader(const string& refHdr, const string& outputHdr, const std::vector<string>& filenames, const std::vector<vnl_vector_fixed< double,3 > >& gvectors )
	{
		bool modifyGradients = (gvectors.size()!=0) ;

		//-------------------------------------
		// Some checkings
		//-------------------------------------
		if ( modifyGradients && filenames.size() != gvectors.size() )
			throw itk::ExceptionObject(__FILE__,__LINE__,"The list of filenames and the list of gradient vectors don't have the same size", "NHDR_CopyAndModifyHeader");

		//-------------------------------------
		// Open input header file to read the list of file names
		//-------------------------------------
		std::ifstream ifile;
		ifile.open ( refHdr.c_str(), ofstream::in);
		if  ( !ifile.is_open() ) 
			throw itk::ExceptionObject(__FILE__,__LINE__,"Cannot open the input header file for reading", "NHDR_CopyAndModifyHeader");

		//-------------------------------------
		// Open output header file for writing
		//-------------------------------------
		std::ofstream ofile;
		ofile.open ( outputHdr.c_str(), ofstream::out );
		if  ( !ofile.is_open() ) 
			throw itk::ExceptionObject(__FILE__,__LINE__,"Cannot open the output header file for writing", "NHDR_CopyAndModifyHeader");

		//-------------------------------------
		// Check that the input NHDR is valid,
		// and get the number of gradient dirs
		//-------------------------------------
		int nbRefImgs;
		NHDR_IsValid(refHdr, nbRefImgs);

		// Check the number of files only if no new gradient vectors are defined
		if ( (!modifyGradients) && nbRefImgs!=(int)filenames.size() )
			throw itk::ExceptionObject(__FILE__,__LINE__,"Error. The number of images in the reference is not the same. Cannot continue (the gradients will be false)", "NHDR_CopyAndModifyHeader");

		//-------------------------------------
		// Copy until the string "data file:" or "DWMRI_gradient_" )
		//-------------------------------------
		char szBuffer[4096];
		bool next=true;
		while ( (!ifile.eof()) && next)
		{
			ifile.getline(szBuffer, 4095);
			string line = string(szBuffer);

			// If we arrived to the files definition, stop copying the header 
			if ( line.find("data file:")!=string::npos )
				next=false;

			// If we currently are on a gradient vector definition
			else if ( modifyGradients && (line.find("DWMRI_gradient_") !=string::npos)  )
				next=false;			// Stop it, we are going to modify the gradients

			// If we modify both the filenames and the gradients, modify the sizes
			else if ( modifyGradients & (line.find("sizes:")!=string::npos) )
			{
				int x, y, z, t;
				// Get and modify the sizes
				if ( sscanf(szBuffer, "sizes: %d %d %d %d", &x, &y, &z, &t)==4 )
					ofile << "sizes: "<<x<<" "<<y<<" "<<z<<" "<<gvectors.size()<<endl;
				else
				{
					string errMsg = "Error while reading the image sizes (LINE="+line+").";
					throw itk::ExceptionObject(__FILE__,__LINE__,errMsg, "NHDR_CopyAndModifyHeader");
				}
			}
			else
				ofile << szBuffer << endl;
		}

		ifile.close();

		// If we arrived at the end of the NHDR, and next==true, this is an error!
		if  ( next ) 
		{
			ofile.close();
			throw itk::ExceptionObject(__FILE__,__LINE__,"Error. Cannot find 'data file:' in the input header file", "");
		}
		//-------------------------------------
		// If we have to modify the gradients
		//-------------------------------------
		if ( modifyGradients )
		{
			for ( unsigned int g=0 ; g<gvectors.size() ; g++ )
			{
				vnl_vector<double> v = gvectors[g];
				ofile << "DWMRI_gradient_" << std::right << std::setw(4) << std::setfill('0') << g << ":= " << v[0] <<" " << v[1] << " " << v[2] << " " << std::endl;
			}
		}

		//-------------------------------------
		// Now modify the files
		//-------------------------------------
		ofile << "data file: LIST"<<endl;		//ensure it is a 'list' mode
		for ( unsigned int i=0; i<filenames.size() ; i++ )
		{
			crl::FileName fn(filenames[i]);
			ofile << fn.getFileNameWithExt() << std::endl;
		}
		ofile.close();

	}


	/**********************************************************************************************//**
	 * \fn	void NHDR_CopyAndModifyHeader(const string& refHdr, const string& outputHdr,
	 * 		const std::vector<string>& filenames)
	 *
	 * \brief	Copy a NHDR header in a new NHDR header with different data files. 
	 * \warning The number of new data files must be the same as in the original header.
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \exception	Thrown if error . 
	 *
	 * \param	refHdr		The reference header filename. 
	 * \param	outputHdr	The output header filename. 
	 * \param	filenames	A vector containing the new filenames.  
	*************************************************************************************************/
	void NHDR_CopyAndModifyHeader(const string& refHdr, const string& outputHdr, const std::vector<string>& filenames)
	{
		std::vector< vnl_vector_fixed< double,3 > > emptyVecs;
		NHDR_CopyAndModifyHeader(refHdr, outputHdr, filenames, emptyVecs);
	}




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
	void NHDR_GetGradientVectors( const itk::MetaDataDictionary& imgMetaDictionary,
		std::vector< vnl_vector_fixed< double,3 >  >& gVectors,
		std::vector< double >& bValues,
		bool includeNullGradients )
	{
		//--------------------------------------------
		// Some initialization
		//--------------------------------------------
		typedef vnl_vector_fixed< double,3 > GradientVectorType;
		GradientVectorType vect3d;

		bool m_ReadB0 = false;
		double m_B0 = 0;
		int nbNullGrad = 0;

		std::vector< vnl_vector_fixed< double,3 >  > tmpVectors;
		tmpVectors.clear();
		gVectors.clear();
		bValues.clear();

		//--------------------------------------------
		// Read the keys in the meta dictionary
		//--------------------------------------------
		std::vector<std::string> imgMetaKeys = imgMetaDictionary.GetKeys();
		std::vector<std::string>::const_iterator itKey = imgMetaKeys.begin();
		std::string metaString;

		for (; itKey != imgMetaKeys.end(); itKey ++)
		{
			double x,y,z;

			itk::ExposeMetaData<std::string> (imgMetaDictionary, *itKey, metaString);

			//--------------------------------------------
			// Is it a gradient vector?
			//--------------------------------------------
			if (itKey->find("DWMRI_gradient") != std::string::npos)
			{ 
				//if ( m_Verbose ) std::cout << *itKey << " ---> " << metaString << std::endl;      
				sscanf(metaString.c_str(), "%lf %lf %lf\n", &x, &y, &z);
				vect3d[0] = x; vect3d[1] = y; vect3d[2] = z;

				if ( vect3d.one_norm()<1e-8 ) nbNullGrad++;
				tmpVectors.push_back(vect3d);
			}

			//--------------------------------------------
			// Is it the nominal B-value?
			//--------------------------------------------
			else if (itKey->find("DWMRI_b-value") != std::string::npos)
			{
				//if ( m_Verbose ) std::cout << *itKey << " ---> " << metaString << std::endl;      
				m_ReadB0 = true;
				m_B0 = atof(metaString.c_str());
			}
		}

		//--------------------------------------------
		// Check that the nominal b-value has been read
		//--------------------------------------------
		if(!m_ReadB0)
			throw itk::ExceptionObject(__FILE__,__LINE__,"Nominal B-Value not found in the NHDR file", "");
	
		std::cout << "- Found " << tmpVectors.size()-nbNullGrad << " gradient vectors, "<<nbNullGrad<<" B=0 and a nominal B-value of "<<m_B0<<"."<<std::endl;

		//--------------------------------------------
		// Now normalize the gradients and compute the b-values
		//--------------------------------------------
		for ( unsigned int i=0; i<tmpVectors.size(); i++ )
		{
			vect3d = tmpVectors[i];
			
			if ( vect3d.one_norm() <= 1e-8 )
			{
				vect3d[0]=vect3d[1]=vect3d[2]=0;
				if ( includeNullGradients ) 
				{
					gVectors.push_back(vect3d);
					bValues.push_back(0);
				}
			}
			else
			{
				double nn = vect3d.two_norm();
				double b = m_B0 * nn * nn;

				gVectors.push_back( vect3d / nn );
				bValues.push_back(b);
			}
		}
	}

	/**********************************************************************************************//**
	 * \fn	itk::Matrix<double,3,
	 * 		3> NHDR_GetMeasurementFrame( const itk::MetaDataDictionary& imgMetaDictionary,
	 * 		bool verbose )
	 *
	 * \brief	Get the measurement frame of the NHDR/NRRD file as an itk::Matrix. If not there print
	 * 			a warning and return the identity matrix. 
	 *
	 * \author	Benoit Scherrer
	 * \date	July 2010
	 *
	 * \param	imgMetaDictionary	Dictionary of image metas. 
	 * \param	verbose				true to print out the extracted measurement frame. 
	 *
	 * \return	The measurement frame. 
	*************************************************************************************************/
	itk::Matrix<double,3,3> NHDR_GetMeasurementFrame( const itk::MetaDataDictionary& imgMetaDictionary, bool verbose )
	{
		itk::Matrix<double,3,3> measurementFrame;
		measurementFrame.SetIdentity();

		if (imgMetaDictionary.HasKey("NRRD_measurement frame"))
		{
			std::vector<std::vector<double> > frame;
			itk::ExposeMetaData<std::vector<std::vector<double> > >(imgMetaDictionary,"NRRD_measurement frame",frame);
			if ( verbose )
			{
				std::cout << "Measurement frame found in NRRD file: (";
				std::cout << frame[0][0] << " " << frame[0][1] << " " << frame[0][2] << ") (";
				std::cout << frame[1][0] << " " << frame[1][1] << " " << frame[1][2] << ") (";
				std::cout << frame[2][0] << " " << frame[2][1] << " " << frame[2][2] << ")" << std::endl;
			}
			measurementFrame[0][0] = frame[0][0];
			measurementFrame[0][1] = frame[0][1];
			measurementFrame[0][2] = frame[0][2];
			measurementFrame[1][0] = frame[1][0];
			measurementFrame[1][1] = frame[1][1];
			measurementFrame[1][2] = frame[1][2];
			measurementFrame[2][0] = frame[2][0];
			measurementFrame[2][1] = frame[2][1];
			measurementFrame[2][2] = frame[2][2];
		}
		else
		{
			std::cout << "No measurement frame found in NRRD file, assuming identity." << std::endl;
		}

		return measurementFrame;
	}

	/**********************************************************************************************//**
	 * \fn	itk::MetaDataDictionary& NHDR_SetMeasurementFrame( itk::MetaDataDictionary& imgMetaDictionary,
	 * 		const itk::Matrix<double,3,3>& measurementFrame )
	 *
	 * \brief	Set the measurement frame of the NHDR/NRRD file from an itk::Matrix in a meta fata
	 * 			dictionary. Modify the input imgMetaDictionary (and return its reference for
	 * 			convenience). 
	 *
	 * \author	Benoit Scherrer
	 * \date	July 2010
	 *
	 * \param [in,out]	imgMetaDictionary	Dictionary of image metas. 
	 * \param	measurementFrame			The measurement frame. 
	 *
	 * \return	. 
	*************************************************************************************************/
	itk::MetaDataDictionary& NHDR_SetMeasurementFrame( itk::MetaDataDictionary& imgMetaDictionary, const itk::Matrix<double,3,3>& measurementFrame )
	{
		std::vector<std::vector<double> > frame(3);

		for ( int i=0; i<3; i++ )
		{
			frame[i].resize(3);
			for( int j=0; j<3; j++ )
				frame[i][j] = measurementFrame[i][j];
		}

		itk::EncapsulateMetaData<std::vector<std::vector<double> > >(imgMetaDictionary,"NRRD_measurement frame",frame);
		return imgMetaDictionary;
	}
	
	itk::Matrix<double,3,3> NHDR_ApplyTransformToMeasurementFrame(const itk::Matrix<double,3,3>& measurementFrame, const itk::AffineTransform<double, 3>::Pointer & transform )
	{
		return NHDR_ApplyTransformToMeasurementFrame<double, 3, 3>( measurementFrame, transform.GetPointer() );

	}

	/**********************************************************************************************//**
	 * \fn	vnl_vector_fixed<double, 3> NHDR_ApplyMeasurementFrameToOneVector(const itk::Matrix<double,3,
	 * 		3>& measurementFrame, const vnl_vector_fixed<double, 3>& gvector )
	 *
	 * \brief	Apply the measurement frame to a vector and returns a new vector. 
	 *
	 * \author	Benoit Scherrer
	 * \date	July 2010
	 *
	 * \param	measurementFrame	The measurement frame. 
	 * \param	gvector				The input gradient vector. 
	 *
	 * \return	a new gradient vector. 
	*************************************************************************************************/
	vnl_vector_fixed<double, 3> NHDR_ApplyMeasurementFrameToOneVector(const itk::Matrix<double,3,3>& measurementFrame, const vnl_vector_fixed<double, 3>& gvector )
	{
		//------------------------------------------
		// Compute the inverse of the measurement frame
		// matrix to apply to the gradient vectors
		//------------------------------------------
		typedef itk::AffineTransform< double, 3 >  AffineTransformType;
		AffineTransformType::Pointer mfTransform = AffineTransformType::New();
		AffineTransformType::Pointer inverse = AffineTransformType::New();
		mfTransform->SetMatrix(measurementFrame);
		mfTransform->GetInverse(inverse);

		//------------------------------------------
		// Transform the vector
		//------------------------------------------
		vnl_vector_fixed<double, 3> newV;
		AffineTransformType::InputVnlVectorType g = inverse->TransformVector(gvector);
		newV[0] = g[0];
		newV[1] = g[1];
		newV[2] = g[2];

		return newV;

	}

	/**********************************************************************************************//**
	 * \fn	void NHDR_ApplyMeasurementFrameToVectorSet(const itk::Matrix<double,3,
	 * 		3>& measurementFrame, const std::vector<vnl_vector_fixed<double, 3> > gvectors )
	 *
	 * \brief	Apply the measurement frame to a vector set and returns a new vector set. 
	 *
	 * \author	Benoit Scherrer
	 * \date	July 2010
	 *
	 * \param	measurementFrame	The measurement frame. 
	 * \param	gvectors			The gvectors. 
	*************************************************************************************************/
	std::vector<vnl_vector_fixed<double, 3> > NHDR_ApplyMeasurementFrameToVectorSet(const itk::Matrix<double,3,3>& measurementFrame, const std::vector<vnl_vector_fixed<double, 3> >& gvectors )
	{
		//------------------------------------------
		// Compute the inverse of the measurement frame
		// matrix to apply to the gradient vectors
		//------------------------------------------
		typedef itk::AffineTransform< double, 3 >  AffineTransformType;
		AffineTransformType::Pointer mfTransform = AffineTransformType::New();
		AffineTransformType::Pointer inverse = AffineTransformType::New();
		mfTransform->SetMatrix(measurementFrame);
		mfTransform->GetInverse(inverse);

		//------------------------------------------
		// Transform all vectors
		//------------------------------------------
		std::vector<vnl_vector_fixed<double, 3> > newVectors;
		for ( unsigned int i=0; i<gvectors.size() ; i++ )
		{
			AffineTransformType::InputVnlVectorType g = inverse->TransformVector(gvectors[i]);
			newVectors.push_back(g);
		}		

		return newVectors;
	}
}


