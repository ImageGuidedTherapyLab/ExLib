/**********************************************************************************************//**
 * \file	crlNHDRFileUtils.h
 *
 * \brief	Declare useful functions for NHDR file manipulations
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


#ifndef CRL_NHDR_FILE_UTILS
#define CRL_NHDR_FILE_UTILS


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <vnl/vnl_vector.h>

#include <itkMetaDataDictionary.h>
#include <itkTransform.h>
#include <itkAffineTransform.h>

#include "crlVectorImageUtils.h"
#include "crlRobustMeanImageFilter.h"

//#ifndef DEPRECATED
//#if defined(_MSC_VER)
//#define DEPRECATED __declspec(deprecated("Warning. This function is deprecated"))
//#else
//#define DEPRECATED
//#endif
//#endif

namespace crl {
	extern void NHDR_CopyAndModifyHeader(const std::string& refHdr, const std::string& outputHdr, const std::vector<std::string>& filenames, const std::vector<vnl_vector_fixed< double,3 > >& gvectors );
	extern void NHDR_CopyAndModifyHeader(const std::string& refHdr, const std::string& outputHdr, const std::vector<std::string>& filenames);
	extern void NHDR_IsValid( const std::string& fileHdr, int &nbImages );

	extern void NHDR_GetGradientVectors( const itk::MetaDataDictionary& imgMetaDictionary,
		std::vector< vnl_vector_fixed< double,3 >  >& gVectors,
		std::vector< double >& bValues,
		bool includeNullGradients=true);

	extern itk::Matrix<double,3,3> NHDR_GetMeasurementFrame( const itk::MetaDataDictionary& imgMetaDictionary, bool verbose=true );
	extern itk::MetaDataDictionary& NHDR_SetMeasurementFrame( itk::MetaDataDictionary& imgMetaDictionary, const itk::Matrix<double,3,3>& measurementFrame );
	vnl_vector_fixed<double, 3> NHDR_ApplyMeasurementFrameToOneVector(const itk::Matrix<double,3,3>& measurementFrame, const vnl_vector_fixed<double, 3>& gvector );
	std::vector<vnl_vector_fixed<double, 3> > NHDR_ApplyMeasurementFrameToVectorSet(const itk::Matrix<double,3,3>& measurementFrame, const std::vector<vnl_vector_fixed<double, 3> >& gvectors );

	/**********************************************************************************************//**
	 * \fn	template < typename TScalarType, unsigned int NInputDimensions,
	 * 		unsigned int NOutputDimensions> itk::Matrix<double,3,
	 * 		3> NHDR_ApplyTransformToMeasurementFrame(const itk::Matrix<double,3,3>& mesurementFrame,
	 * 		const itk::TransformType<TScalarType, NInputDimensions, NOutputDimensions>& transform )
	 *
	 * \brief	Apply a transform to a measurement frame. Warning: if the transform comes from a
	 * 			registration (for example), the INVERSE of the registration transform should be
	 * 			provided to that function.
	 * 			
	 * 			\remark This function is a template to be able to be compatible with any itk
	 * 			transform. 
	 *
	 * \author	Benoit Scherrer
	 * \date	July 2010
	 *
	 * \param	mesurementFrame	The input mesurement frame. 
	 * \param	transform		The input transform. 
	 *
	 * \return	The new mesurement frame. 
	*************************************************************************************************/
	template < typename TScalarType, unsigned int NInputDimensions, unsigned int NOutputDimensions>
	itk::Matrix<double,3,3> NHDR_ApplyTransformToMeasurementFrame(const itk::Matrix<double,3,3>& measurementFrame, const itk::Transform<TScalarType, NInputDimensions, NOutputDimensions> *transform )
	{
		itk::Matrix<double,3,3> outMeasurementFrame;

		//---------------------------------------------
		// Apply the transform
		//---------------------------------------------
		itk::Vector<double,3> vec;
		vec = transform->TransformVector(measurementFrame[0]);
		outMeasurementFrame[0][0] = vec[0];
		outMeasurementFrame[0][1] = vec[1];
		outMeasurementFrame[0][2] = vec[2];

		vec = transform->TransformVector(measurementFrame[1]);
		outMeasurementFrame[1][0] = vec[0];
		outMeasurementFrame[1][1] = vec[1];
		outMeasurementFrame[1][2] = vec[2];

		vec = transform->TransformVector(measurementFrame[2]);
		outMeasurementFrame[2][0] = vec[0];
		outMeasurementFrame[2][1] = vec[1];
		outMeasurementFrame[2][2] = vec[2];
		
		return outMeasurementFrame;
	}

	/**********************************************************************************************//**
	 * \fn	itk::Matrix<double,3,3> NHDR_ApplyTransformToMeasurementFrame(const itk::Matrix<double,3,
	 * 		3>& measurementFrame, const itk::AffineTransform<double, 3>::Pointer & transform )
	 *
	 * \brief	Apply an affine transform to a measurement frame (provided for convenience)
	 *
	 * \author	Benoit Scherrer
	 * \date	July 2010
	 *
	 * \param	measurementFrame	The measurement frame. 
	 * \param	transform			The affine transform. 
	 *
	 * \return	The new measurement frame. 
	*************************************************************************************************/
	itk::Matrix<double,3,3> NHDR_ApplyTransformToMeasurementFrame(const itk::Matrix<double,3,3>& measurementFrame, const itk::AffineTransform<double, 3>::Pointer & transform );

	/**********************************************************************************************//**
	 * \enum	NHDR_LoadOptions
	 *
	 * \brief	Options for  NHDR_Load
	*************************************************************************************************/
	typedef enum {
		NHDR_KEEP_ALL_NULL_B0,
		NHDR_COMPUTE_MEAN_B0,
		NHDR_COMPUTE_ROBUSTMEAN_B0
	} NHDR_LoadOptions;
	typedef unsigned int NHDR_LoadOptionsType;

	/**********************************************************************************************//**
	 * \fn	template <typename TPixel,
	 * 		unsigned int VImageDimension > typename itk::VectorImage< TPixel,
	 * 		VImageDimension >::Pointer NHDR_Load ( const std::string& filename,
	 * 		std::vector<vnl_vector_fixed< double,3 > > &gvecs, std::vector<double>& bvalues,
	 * 		NHDR_LoadOptionsType options)
	 *
	 * \brief	Load the DWI RAW data from a NHDR/NRRD file and extract the gradient vectors and b-
	 * 			values. The parameter 'options' defines the behaviour when several B=0 images are
	 * 			present in the NHDR. If set to NHDR_COMPUTE_ROBUSTMEAN_B0, the robust mean is
	 * 			actually computed only if there are at least 3 B=0 images. 
	 *
	 * \author	Benoit Scherrer
	 * \date	July 2010
	 *
	 * \param	filename		Filename of the file. 
	 * \param [in,out]	gvecs	The gvecs. 
	 * \param [in,out]	bvalues	The bvalues. 
	 * \param	options			Options for controlling the loading. 
	 *
	 * \return	. 
	*************************************************************************************************/
	template  <typename TPixel, unsigned int VImageDimension >
	typename itk::VectorImage< TPixel, VImageDimension >::Pointer
		NHDR_Load ( const std::string& filename, 
								std::vector<vnl_vector_fixed< double,3 > > &gvecs,
								std::vector<double>& bvalues,
								NHDR_LoadOptionsType options)
	{
		//---------------------------------------
		// First load the RAW vector image
		//---------------------------------------
		std::cout<<"- Loading the DWI file <" << filename << ">..." << std::endl;

		typedef itk::VectorImage<TPixel, VImageDimension> DWIStudyType;
		typedef itk::ImageFileReader<DWIStudyType> DWIStudyReaderType;
		typename DWIStudyReaderType::Pointer reader = DWIStudyReaderType::New();
		reader->SetFileName(filename);
		reader->Update();
		typename DWIStudyType::Pointer dwiStudy = reader->GetOutput();

		//---------------------------------------
		// Then load the gradients/bvalues
		//---------------------------------------
		NHDR_GetGradientVectors( dwiStudy->GetMetaDataDictionary(), gvecs, bvalues, true );

		//---------------------------------------
		// If we want to keep all null B0, just returns the result
		//---------------------------------------
		if ( options==NHDR_KEEP_ALL_NULL_B0 )
			return dwiStudy;

		//-----------------------------------------
		// Else extracts all the B=0 images and connect them to the mean filter
		//-----------------------------------------
		std::cout<<"- Extract B=0 images..."<<std::endl;
		typedef itk::Image<TPixel, VImageDimension> DWIOneImageType;
		typedef crl::RobustMeanImageFilter< DWIOneImageType, DWIOneImageType >  RobustMeanImageFilterType;
		typename RobustMeanImageFilterType::Pointer MeanB0Filter = RobustMeanImageFilterType::New();

		int nbB0 = 0; 
		int nbNonNullB0 = 0;
		for( unsigned int b=0; b<bvalues.size(); b++ )
		{
			// extract zero-norm gradients
			if( bvalues[b]==0 ) 
			{
				typename DWIOneImageType::Pointer oneImage = crl::ExtractComponentOfVectorImage< TPixel, VImageDimension>( dwiStudy, b );
				MeanB0Filter->SetInput(nbB0, oneImage);
				nbB0++;
			}
			else
				nbNonNullB0++;
		}

		//-----------------------------------------
		// If no null b-value was found, just return the original study
		//-----------------------------------------
		if ( nbB0==0 )
		{
			std::cout<<"  Warning. No B=0 image were found."<<std::endl;
			return dwiStudy;
		}

		//-----------------------------------------
		// if only one was found, and in the first position,
		// just returns
		//-----------------------------------------
		if ( nbB0==1 && bvalues[0]==0 )
			return dwiStudy;

		//-----------------------------------------
		// Else compute the mean/robust mean
		//-----------------------------------------
		if ( options==NHDR_COMPUTE_ROBUSTMEAN_B0 && nbB0>=3 ) 
			std::cout<<"- Compute the robust mean B=0 image..."<<std::endl;
		else
			std::cout<<"- Compute the mean B=0 image..."<<std::endl;

		MeanB0Filter->SetOutputMeanImage(true);
		MeanB0Filter->Update();
		
		//-----------------------------------------
		// Depending on the options, use the mean or the robust mean
		//-----------------------------------------
		typename DWIOneImageType::Pointer meanB0;
		if ( options==NHDR_COMPUTE_ROBUSTMEAN_B0 && nbB0>=3 ) 
			meanB0 = MeanB0Filter->GetOutput();
		else
			meanB0 = MeanB0Filter->GetMeanImage();

		//-------------------------------------------------------
		// Creates the new vector file
		//-------------------------------------------------------
		typename DWIStudyType::Pointer newStudy = DWIStudyType::New();
		newStudy->SetRegions(dwiStudy->GetLargestPossibleRegion());
		newStudy->SetOrigin(dwiStudy->GetOrigin());
		newStudy->SetSpacing(dwiStudy->GetSpacing());
		newStudy->SetDirection(dwiStudy->GetDirection());
		newStudy->SetNumberOfComponentsPerPixel(nbNonNullB0+1);
		newStudy->Allocate();

		//-------------------------------------------------------
		// Prepare the copy of the data
		//-------------------------------------------------------
		typedef itk::ImageRegionConstIterator<DWIStudyType>		InIteratorType;
		typedef itk::ImageRegionIterator<DWIStudyType>			OutIteratorType;
		typedef itk::ImageRegionConstIterator<DWIOneImageType>	InB0IteratorType;

		InIteratorType inItr(dwiStudy, dwiStudy->GetLargestPossibleRegion() );
		OutIteratorType outItr(newStudy, newStudy->GetLargestPossibleRegion() );
		InB0IteratorType inB0Itr ( meanB0, meanB0->GetLargestPossibleRegion()); 
		inItr.GoToBegin(); inB0Itr.GoToBegin();

		typename DWIStudyType::PixelType inP(dwiStudy->GetNumberOfComponentsPerPixel() ) ;
		typename DWIStudyType::PixelType outP(newStudy->GetNumberOfComponentsPerPixel() ) ;

		//-------------------------------------------------------
		// Copy the data!
		//-------------------------------------------------------
		for ( outItr.GoToBegin() ; !outItr.IsAtEnd() ; ++outItr )
		{
			//-------------------------------
			// First set the B=0 value
			//-------------------------------
			outP[0] = inB0Itr.Value();

			//-------------------------------
			// Then copy only a part of the vector
			//-------------------------------
			inP = inItr.Get();
			int j = 1;
			for ( unsigned int i=0; i<bvalues.size(); i++ )
			{
				if ( bvalues[i]!=0 )
				{
					outP[j] = inP[i];
					j++;
				}
			}
			outItr.Set(outP);

			++inItr;
			++inB0Itr;
		}

		//-----------------------------------------
		// Put only one null b-value (and one null grad) in 
		// the first component of the output bvalues/gvecs
		//-----------------------------------------
		std::vector<vnl_vector_fixed< double,3 > > new_gvecs;
		std::vector<double> new_bvalues;

		vnl_vector_fixed< double,3 > vnull;
		vnull.fill(0);
		new_bvalues.push_back(0);
		new_gvecs.push_back(vnull);

		for ( unsigned int i=0; i<bvalues.size(); i++ )
		{
			if ( bvalues[i]!=0 )
			{
				new_bvalues.push_back(bvalues[i]);
				new_gvecs.push_back(gvecs[i]);
			}
		}

		bvalues.clear(); gvecs.clear();
		for ( unsigned int i=0; i<new_bvalues.size(); i++ )
		{
			bvalues.push_back(new_bvalues[i]);
			gvecs.push_back(new_gvecs[i]);
		}

		//-----------------------------------------
		// Returns the new study (with the dictionnary!)
		//-----------------------------------------
		newStudy->SetMetaDataDictionary(dwiStudy->GetMetaDataDictionary());
		return (newStudy);

	}

	/**********************************************************************************************//**
	 *
	 * \brief	Exports a itk::VectorImage in a NHDR file by using a reference NHDR file for the NHDR
	 * 			options.  
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \exception	itk::ExceptionObject	Thrown when exception. 
	 *
	 * \param	vectorImage			The vector image. 
	 * \param	gvecs				The gradient vectors. 
	 * \param	outputBaseFilename	Filename of the output base file. 
	 * \param	refNHDRHeader		The reference nhdr header. 
	 *
	 * ### exception	Thrown	when error. 
	*************************************************************************************************/
	template  <typename TPixel, unsigned int VImageDimension >
	void NHDR_WriteFromVectorImage ( const typename itk::VectorImage< TPixel, VImageDimension >::Pointer& vectorImage, 
								const std::vector<vnl_vector_fixed< double,3 > > &gvecs,
								const std::string& outputBaseFilename,
								const std::string& refNHDRHeader,
                                                                bool useCompression=false )
	{
		typedef itk::VectorImage<TPixel, VImageDimension >	InputImageType;
		typedef itk::ImageFileWriter<InputImageType> WriterType;

		// if gvecs.size==0, we use the same gradients as before
		if ( (gvecs.size()!=0) && (vectorImage->GetNumberOfComponentsPerPixel()!=gvecs.size()) )
			throw itk::ExceptionObject(__FILE__,__LINE__,"The number of gradient vectors must be the same as the number of direction images.", "NHDR_WriteFromVectorImage");


		crl::FileName fnOutput(outputBaseFilename);
		fnOutput.setExtension("nhdr");			//ensure correct extension

		//-------------------------------------------------------
		// First extract all images of the vector image
		//-------------------------------------------------------
		std::vector<std::string> createdFileNames;
		WriteAllComponentsOfVectorImage<TPixel, VImageDimension> ( vectorImage,  fnOutput.getCompleteFilePath(), createdFileNames, useCompression );

		//-------------------------------------------------------
		// Then create the NHDR by considering the reference NHDR
		//-------------------------------------------------------
		NHDR_CopyAndModifyHeader(refNHDRHeader, fnOutput.getCompleteFilePath(), createdFileNames, gvecs );


	}

	/**********************************************************************************************//**
	 *
	 * \brief	Exports a itk::VectorImage in a complete new NHDR file. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \exception	itk::ExceptionObject	Thrown when error. 
	 *
	 * \param	vectorImage			The vector image. 
	 * \param	bvecs				The gradient vectors. 
	 * \param	bvals				The b-values for each gradient vector. 
	 * \param	outputBaseFilename	Filename of the output base file. 
	 * \param	measurementFrame	If defined, the measurement frame (ex: "(1,0,0) (0,1,0) (0,0,1)"). 
	 *
	*************************************************************************************************/
	template  <typename TPixel, unsigned int VImageDimension >
	void NHDR_WriteFromVectorImage ( const typename itk::VectorImage< TPixel, VImageDimension >::Pointer& vectorImage, 
								const std::vector<vnl_vector_fixed< double,3 > > &bvecs,
								const std::vector<double> &bvals,
								const std::string& outputBaseFilename,
                                                                bool useCompression=false,
								const std::string& measurementFrame="" )
	{
		typedef itk::VectorImage<TPixel, VImageDimension >	InputImageType;
		typedef itk::ImageFileWriter<InputImageType> WriterType;

		if ( bvecs.size()!=bvals.size() )
			throw itk::ExceptionObject(__FILE__,__LINE__,"The gradient vectors and the B-value list should have the same size.", "NHDR_WriteFromVectorImage");

		crl::FileName fnOutput(outputBaseFilename);
		fnOutput.setExtension("nhdr");			//ensure correct extension

		//-------------------------------------------------------
		// First extract all images of the vector image
		//-------------------------------------------------------
		std::vector<std::string> createdFileNames;
		WriteAllComponentsOfVectorImage<TPixel, VImageDimension> ( vectorImage,  fnOutput.getCompleteFilePath(), createdFileNames, useCompression );

		//-------------------------------------------------------
		// Prepare some strings for writing the header...
		//-------------------------------------------------------
		std::stringstream ssDim("");
		ssDim << "dimension: " <<(VImageDimension+1) ;

		// Data Size
		std::stringstream ssSize("");
		ssSize<<"sizes:";
		for ( unsigned int i=0; i<VImageDimension ; i++ ) ssSize << " " << vectorImage->GetLargestPossibleRegion().GetSize()[i];
		ssSize << " " << vectorImage->GetNumberOfComponentsPerPixel();

		// Origin
		std::stringstream ssOrigin("");
		ssOrigin<<"space origin: (";
		for ( unsigned int i=0; i<VImageDimension ; i++ ) 
		{
			if ( i!=0 ) ssOrigin<<",";
			ssOrigin << " " << vectorImage->GetOrigin()[i];
		}
		ssOrigin<<")";

		// Data type
		std::stringstream ssType("");
		ssType<<"type: ";
		TPixel trick=0;
		try { throw trick; }
		catch ( char  ) { ssType<<"char"; }
		catch ( unsigned char  ) { ssType<<"unsigned char"; }
		catch ( short  ) { ssType<<"short"; }
		catch ( unsigned short  ) { ssType<<"unsigned short"; }
		catch ( int  ) { ssType<<"int"; }
		catch ( unsigned int  ) { ssType<<"unsigned int"; }
		catch ( long int  ) { ssType<<"long"; }
		catch ( unsigned long int  ) { ssType<<"unsigned long"; }
		catch ( float  ) { ssType<<"float"; }
		catch ( double  ) { ssType<<"double"; }
		catch (...)
		{
			throw itk::ExceptionObject(__FILE__,__LINE__,"Unknown file type", "");
		}

		// Direction
		std::stringstream ssDir("");
		ssDir << "space directions:";
		typename InputImageType::DirectionType imgDir = vectorImage->GetDirection();
		typename InputImageType::SpacingType imgSpacing = vectorImage->GetSpacing();

		for ( unsigned int i=0; i<VImageDimension ; i++ )
		{
			ssDir << " (" ;
			for ( unsigned int j=0; j<VImageDimension ; j++ )
			{
				if ( j!=0 ) ssDir<<",";
				ssDir << (imgDir[j][i] * imgSpacing[i]);
			}
			ssDir << ")";
		}
		//cout << ssDir.str().c_str()<<endl;

		// B-value
		std::stringstream ssBVal("");
		std::vector<double>::const_iterator posMaxB = std::max_element (bvals.begin(), bvals.end());
		double maxBValue = *posMaxB ;
		ssBVal<< "DWMRI_b-value:=" << maxBValue ;

		//-------------------------------------------------------
		// Open the header for writing
		//-------------------------------------------------------
		std::ofstream nhdrFile;
		nhdrFile.open ( fnOutput.getCompleteFilePath().c_str(), std::ofstream::out);
		if  ( !nhdrFile.is_open() ) 
			throw itk::ExceptionObject(__FILE__,__LINE__,"Cannot open the NHDR file for writing", "");

		//-------------------------------------------------------
		// Write the header
		//-------------------------------------------------------
	    // Header infos
		nhdrFile << "NRRD0005" << std::endl;
	    nhdrFile << "# Complete NRRD file format specification at:" << std::endl;
	    nhdrFile << "# http://teem.sourceforge.net/nrrd/format.html" << std::endl;
	    nhdrFile << "content: epi-dti-asset" << std::endl;
		nhdrFile << ssDim.str().c_str() << std::endl;
		nhdrFile << ssSize.str().c_str() << std::endl;
	    nhdrFile << ssType.str().c_str() << std::endl;
	    nhdrFile << "encoding: raw" << std::endl;
	    nhdrFile << "endian: little" << std::endl;
	    nhdrFile << "byte skip: -1" << std::endl;
	    nhdrFile << "space: left-posterior-superior" << std::endl;
	    nhdrFile << "centers: cell cell cell none" << std::endl;
	    nhdrFile << ssDir.str().c_str() << " none" << std::endl;
	    nhdrFile << "space units: \" mm\" \"mm\" \"mm\"" << std::endl;
		nhdrFile << ssOrigin.str().c_str() << std::endl;
		nhdrFile << "kinds: space space space list" << std::endl;
		nhdrFile << "measurement frame: " << ((measurementFrame=="")?"(1,0,0) (0,1,0) (0,0,1)":measurementFrame) << std::endl;
	    nhdrFile << "modality:=DWMRI" << std::endl;
	    nhdrFile << ssBVal.str().c_str() << std::endl;

		// gradients vectors
	    for ( unsigned int g=0; g<bvecs.size() ; g++ )
		{
			vnl_vector_fixed< double,3 > v;
			
			// Get the normalized gradient vector
			v = bvecs[g];
			double nn = v.two_norm();
			if ( nn!=0 ) v = v/nn;

			// Put the good norm (depending on the max B-value)
			v  = v * (sqrt(bvals[g]/maxBValue));

			nhdrFile << "DWMRI_gradient_" << std::right << std::setw(4) << std::setfill('0') << g << ":= " << v[0] <<" " << v[1] << " " << v[2] << " " << std::endl;
		}		
		// files
	    nhdrFile << "data file: LIST" << std::endl;
		for ( unsigned int i=0; i<createdFileNames.size() ; i++ )
		{
			crl::FileName fn(createdFileNames[i]);
			nhdrFile << fn.getFileNameWithExt() << std::endl;
		}
	    nhdrFile.close();
	}
}

#endif
