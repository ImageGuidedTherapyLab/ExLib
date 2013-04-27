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

#ifndef CRL_TRANSFORM_UTILS
#define CRL_TRANSFORM_UTILS

#ifdef _MSC_VER
  #pragma warning( disable: 4996 )
#endif

#include <iostream>
#include <iomanip>
#include <string>

#include <itkTransformFileReader.h>
#include <itkTranslationTransform.h>
#include <itkVersorRigid3DTransform.h>
#include <itkSimilarity3DTransform.h>
#include <itkScaleSkewVersor3DTransform.h>
#include <itkAffineTransform.h>

namespace crl {

	template  <typename TScalarType, unsigned int NDimensions >
	typename itk::AffineTransform< TScalarType, NDimensions >::Pointer
	ConvertAnyTransformToAffine ( typename itk::TransformBase::Pointer inputTransform )
	{
		typedef itk::AffineTransform< TScalarType, NDimensions >		AffineTransformType;
		typedef itk::ScaleSkewVersor3DTransform< TScalarType >			ScaleSkewVersor3DTransformType;
		typedef itk::VersorTransform< TScalarType >						GenericVersorTransformType;
		typedef itk::Similarity3DTransform< TScalarType >  				Similarity3DTransformType;
		typedef itk::VersorRigid3DTransform< TScalarType >  			VersorRigid3DTransformType;
		typedef itk::TranslationTransform< TScalarType, NDimensions >	TranslationTransformType;

		typename AffineTransformType::Pointer outputTransform = AffineTransformType::New();
		outputTransform->SetIdentity();

		if (!strcmp(inputTransform->GetNameOfClass(), "AffineTransform"))
		{
			typename AffineTransformType::Pointer affine_read = static_cast<AffineTransformType*>(inputTransform.GetPointer());

			// Setting the matrix and offset is sufficient to duplicate the transform 
			// performance operating on points, 
			// but does not set the center appropriately for follow on registration.
			outputTransform->SetTranslation(affine_read->GetTranslation());
			outputTransform->SetCenter(affine_read->GetCenter());
			outputTransform->SetMatrix(affine_read->GetMatrix());

		} 
		else if (!strcmp(inputTransform->GetNameOfClass(), "TranslationTransform")) 
		{
			typename TranslationTransformType::Pointer  translation_read = static_cast<TranslationTransformType*>(inputTransform.GetPointer());
			outputTransform->SetOffset(translation_read->GetOffset());
		} 
		else if (!strcmp(inputTransform->GetNameOfClass(), "VersorRigid3DTransform")) 
		{
			typename VersorRigid3DTransformType::Pointer  versorrigid_read = static_cast<VersorRigid3DTransformType*>(inputTransform.GetPointer());
			outputTransform->SetCenter(versorrigid_read->GetCenter());
			outputTransform->SetTranslation(versorrigid_read->GetTranslation());
			outputTransform->SetMatrix(versorrigid_read->GetMatrix());
			// The above is sufficient to duplicate the transform performance
			// but does not set the center appropriately for follow on registration.

		} 
		else if (!strcmp(inputTransform->GetNameOfClass(), "Similarity3DTransform"))
		{
			typename Similarity3DTransformType::Pointer  similarity_read = static_cast<Similarity3DTransformType*>(inputTransform.GetPointer());
			outputTransform->SetCenter(similarity_read->GetCenter());
			outputTransform->SetTranslation(similarity_read->GetTranslation());
			outputTransform->SetMatrix(similarity_read->GetMatrix());
		} 
		else if (!strcmp(inputTransform->GetNameOfClass(), "ScaleSkewVersor3DTransform")) 
		{
			typename ScaleSkewVersor3DTransformType::Pointer  scaleskewversor_read = static_cast<ScaleSkewVersor3DTransformType*>(inputTransform.GetPointer());
			outputTransform->SetCenter(scaleskewversor_read->GetCenter());
			outputTransform->SetTranslation(scaleskewversor_read->GetTranslation());
			outputTransform->SetMatrix(scaleskewversor_read->GetMatrix());
		} 
		else if (!strcmp(inputTransform->GetNameOfClass(), "VersorTransform"))
		{
			typename GenericVersorTransformType::Pointer  genversor_read = static_cast<GenericVersorTransformType*>(inputTransform.GetPointer());
			outputTransform->SetCenter(genversor_read->GetCenter());
			outputTransform->SetTranslation(genversor_read->GetTranslation());
			outputTransform->SetMatrix(genversor_read->GetMatrix());
		} 
		else 
		{
			std::string errStr = "Don't know how to convert a " +
				std::string(inputTransform->GetNameOfClass()) + " transform to an affine transform." ;
			throw itk::ExceptionObject(__FILE__,__LINE__,errStr, "");
		}
		return outputTransform;
	}


	template  <typename TScalarType, unsigned int NDimensions >
	typename itk::AffineTransform< TScalarType, NDimensions >::Pointer
		LoadAnyTransformAsAffine ( const std::string& filename )
	{
		//---------------------------------------
		// Open the file
		//---------------------------------------
		typename itk::TransformFileReader::Pointer trsfreader = itk::TransformFileReader::New();
		trsfreader->SetFileName( filename.c_str() );
		try {
			trsfreader->Update();
		} catch ( itk::ExceptionObject & excp )
		{
			// Display error from reading the transform file.
			throw itk::ExceptionObject(__FILE__,__LINE__,
				"Error while reading the transform file "+filename+"\n" 
				+ std::string(excp.GetDescription()));
		}

		//---------------------------------------
		// Now try to work out how many and what type of transforms were read.
		// We only want to get one transform.
		//---------------------------------------
		typedef itk::TransformFileReader::TransformListType * TransformListType;
		TransformListType transforms = trsfreader->GetTransformList();
		//std::cout << "  Number of transforms = " << transforms->size() << std::endl;
		itk::TransformFileReader::TransformListType::const_iterator it = transforms->begin();

		// Check that only one transform
		if (transforms->size() <= 0 || transforms->size() > 1) 
		{
			throw itk::ExceptionObject(__FILE__,__LINE__,"Error. Only support file with a single-transform in it.");
		}

		//---------------------------------------
		// Converts and returns the transform
		//---------------------------------------
		return ConvertAnyTransformToAffine<TScalarType, NDimensions> (*it);
	}

}

#endif
