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

#ifndef CRL_DWI_STUDY_H
#define CRL_DWI_STUDY_H

#include <itkVectorImage.h>
#include <itkImage.h>

namespace crl{

	/**********************************************************************************************//**
	 * \class	DWIStudy
	 *
	 * \brief	DWI study. 
	 *
	 * \author	Benoit Scherrer
	 * \date	September 2010
	*************************************************************************************************/
	template <class TPixel>
	class DWIStudy
	{
	public:
		typedef itk::VectorImage<TPixel,3>			DWIGradientImageSetType;
		typedef itk::Image<TPixel,3>				DWIGradientImageType;
		typedef vnl_vector_fixed< double,3 >		GradientVectorType;
		typedef std::vector< GradientVectorType >	GradientVectorSetType;
		typedef std::vector<double>					BValueSetType;

		typedef enum {
			GRADIENTS_SCANNER_COORDINATES = 0x1,
			GRADIENTS_MFRAME_COORDINATES = 0x2,
		} GradientCoordinatesType;

	public:
		DWIStudy();
		~DWIStudy();
		
		DWIStudy(const DWIStudy<TPixel>&);		
		void operator=(const DWIStudy<TPixel>&);	

		void									LoadStudy ( const std::string& fileName );
		void									WriteStudy( const std::string& fileName, bool useCompression=false );
		void									WriteStudyWithNewData( const std::string& fileName, typename DWIGradientImageSetType::Pointer data, bool useCompression=false );

		void									OrientStudyInAxial ();
		typename DWIGradientImageType::Pointer	ComputeMeanB0();
		typename DWIGradientImageType::Pointer	ComputeMeanBNonNull();

		void									SetGradientCoordinates(GradientCoordinatesType gc);
		GradientCoordinatesType					GetGradientCoordinates();

		unsigned int							GetNumberOfB0() const;
		unsigned int							GetNumberOfNonNullVectors() const;
		unsigned int							GetTotalNumberOfImages() const;
		double									GetMaximumBValue() const;

		void									ComputeScaledGradients( GradientVectorSetType& scaledGradients ) const;

	protected:
		void									LoadFSLStudy ( const std::string& fileName );
		void									LoadNHDRStudy ( const std::string& fileName );

	public:
		typename DWIGradientImageSetType::Pointer	DWIData;
		GradientVectorSetType						Gradients;
		BValueSetType								BValues;
		double										NominalBValue;

		GradientVectorSetType						OriginalGradients;	/// Contains the original non-flipped/non-permuted gradients (same as Gradients if OrientStudyInAxial() has not been called)

	protected:
		GradientCoordinatesType						m_GradientCoordinatesMode;
	};


} // end namespace crl.

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlDWIStudy.txx"
#endif

#endif
