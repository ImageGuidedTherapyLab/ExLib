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

/**********************************************************************************************//**
 * \file	crlMFMUtils.h
 *
 * \brief	Declare useful functions for MFM/Multi-tensors
*************************************************************************************************/

#ifndef CRL_MFM_UTILS
#define CRL_MFM_UTILS

#ifdef _MSC_VER
  #pragma warning( disable: 4996 )
#endif

#include <iostream>
#include <iomanip>
#include <string>
#include <itkImageIOBase.h>
#include <itkImageFileReader.h>
#include <itkImage.h>
#include <itkVectorImage.h>

#include "crlVectorImageUtils.h"

namespace crl {

	class crlIOUtils
	{
	public:
		static itk::ImageIOBase::Pointer GetImageIO (const std::string& fileName )
		{
			itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO( fileName.c_str(), itk::ImageIOFactory::ReadMode );
			if ( imageIO.IsNull() )
			{
				std::string err = "Can't make image IO for image <"+fileName+">. Possibly an unsupported type.";
				throw itk::ExceptionObject(__FILE__, __LINE__, err);
			}
			imageIO->SetFileName( fileName );		// necessary?  ITK does it.
			imageIO->ReadImageInformation();

			return imageIO;
		}
	};

	template  <typename TPixel >
	typename itk::Image< TPixel, 4 >::Pointer
	MFM_OpenFractionsImageAs4DImage ( const std::string& filename )
	{
		//------------------------------------------
		// Gets image informations
		//------------------------------------------
		itk::ImageIOBase::Pointer imageIO = crlIOUtils::GetImageIO(filename);
		typename itk::ImageIOBase::IOPixelType pixelType = imageIO->GetPixelType();
		unsigned int dimensionality = imageIO->GetNumberOfDimensions();
		int nbComponents = imageIO->GetNumberOfComponents ();

		typedef typename itk::Image< TPixel, 4 > Image4DType;
		typedef typename itk::VectorImage< TPixel, 3 > VectorImageType;

		//------------------------------------------
		// Open directly the 4D image
		//------------------------------------------
		if ( (pixelType==itk::ImageIOBase::SCALAR) && (dimensionality==4) && (nbComponents==1) )
		{
			typedef typename itk::ImageFileReader<Image4DType> Image4DReaderType;
			typename Image4DReaderType::Pointer reader = Image4DReaderType::New();
			reader->SetFileName( filename );
			reader->Update();
			return reader->GetOutput();
		}

		//------------------------------------------
		// Open as a VectorImage and converts to 4D
		//------------------------------------------
		else if ( (pixelType==itk::ImageIOBase::VECTOR) && (dimensionality==3) && (nbComponents>1) )
		{
			typedef typename itk::ImageFileReader<VectorImageType> VectorImageReaderType;
			typename VectorImageReaderType::Pointer reader = VectorImageReaderType::New();
			reader->SetFileName( filename );
			reader->Update();

			return crl::ConvertVectorImageToNDImage<TPixel, 3>(reader->GetOutput());
		}
		else
		{
			std::string err = "Can't open the fractions file <"+filename+">: Invalid pixel type";
			throw itk::ExceptionObject(__FILE__, __LINE__, err);
		}

	}



	template  <typename TPixel >
	typename itk::VectorImage< TPixel, 3 >::Pointer
	MFM_OpenFractionsImageAsVectorImage ( const std::string& filename )
	{

		//------------------------------------------
		// Gets image informations
		//------------------------------------------
		itk::ImageIOBase::Pointer imageIO = crlIOUtils::GetImageIO(filename);
		typename itk::ImageIOBase::IOPixelType pixelType = imageIO->GetPixelType();
		unsigned int dimensionality = imageIO->GetNumberOfDimensions();
		int nbComponents = imageIO->GetNumberOfComponents ();

		typedef typename itk::Image< TPixel, 4 > Image4DType;
		typedef typename itk::VectorImage< TPixel, 3 > VectorImageType;

		//------------------------------------------
		// Open directly the vector image
		//------------------------------------------
		if ( (pixelType==itk::ImageIOBase::VECTOR) && dimensionality==3 && nbComponents>1  )
		{
			typedef typename itk::ImageFileReader<VectorImageType> VectorImageReaderType;
			typename VectorImageReaderType::Pointer reader = VectorImageReaderType::New();
			reader->SetFileName( filename );
			reader->Update();
			return reader->GetOutput();
		}

		//------------------------------------------
		// Open as a VectorImage and converts to 4D
		//------------------------------------------
		else if ( (pixelType==itk::ImageIOBase::SCALAR) && dimensionality==4 && nbComponents==1 )
		{
			typedef typename itk::ImageFileReader<Image4DType> Image4DReaderType;
			typename Image4DReaderType::Pointer reader = Image4DReaderType::New();
			reader->SetFileName( filename );
			reader->Update();

			return crl::ConvertNDToVectorImage<TPixel, 4>(reader->GetOutput());
		}
		else
		{
			std::string err = "Can't open the fractions file <"+filename+">: Invalid pixel type";
			throw itk::ExceptionObject(__FILE__, __LINE__, err);
		}
	}

	/**********************************************************************************************//**
	 * \fn	template <typename TensorImageType > void MFM_LoadTensors ( const std::vector<std::string> & filenames,
	 * 		std::vector<typename TensorImageType::Pointer>& tensorImages )
	 *
	 * \brief	Load multiple tensors in a std::vector  
	 *
	 * \author	Benoit Scherrer
	 * \date	July 2010
	 *
	 * \param	filenames				The filenames. 
	 * \param [in,out]	tensorImages	The loaded tensor images. 
	*************************************************************************************************/
	template  <typename TensorImageType >
	void MFM_LoadTensors ( const std::vector<std::string> & filenames, std::vector<typename TensorImageType::Pointer>& tensorImages  )
	{
		tensorImages.clear();

		typedef itk::ImageFileReader<TensorImageType> TensorImageReaderType;
		typename TensorImageReaderType::Pointer reader = TensorImageReaderType::New();
		for ( unsigned int i=0; i<filenames.size() ; i++ )
		{
			reader->SetFileName( filenames[i] );
			reader->Update();
			tensorImages.push_back(reader->GetOutput());
			reader->GetOutput()->DisconnectPipeline();
		}
	}

}

#endif
