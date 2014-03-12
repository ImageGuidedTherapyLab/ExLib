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

#ifndef CRL_DWI_RAW_H
#define CRL_DWI_RAW_H

#include <iostream>
#include <vector>
#include <string>

#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkVectorContainer.h"

namespace crl {
namespace DWI {

/**********************************************************************************************//**
 * \class	DWIRaw
 *
 * \brief	Template class to read Diffusion-Weighted Image RAW data 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
*************************************************************************************************/
template <class TPixel>
class DWIRaw 
{	
public:
	/** Define some typedefs. */
	typedef TPixel												GradientImagePixelType ;
	typedef typename itk::VectorImage< GradientImagePixelType, 3 > 		RawImageType;
	typedef typename itk::Image< GradientImagePixelType, 3 > 			GradientOneImageType;

	typedef vnl_vector_fixed< double,3 > 								GradientDirectionType;
	typedef typename itk::VectorContainer< unsigned int, GradientDirectionType > GradientDirectionContainerType;

public:
	DWIRaw( void );
	~DWIRaw();

	void LoadData(const std::string& filename);
	
	GradientDirectionContainerType*	GetDiffusionVectors() const ;
	RawImageType*					GetData() const ;
	const std::vector<std::string>&	GetFileNames() const ;

	double							GetB0() const ;

	int 							GetNumberOfB0Images() const;
	int 							GetNumberOfGradientImages() const;

	void							SetVerbose( bool v );
	bool							GetVerbose() const;

protected:
	typename RawImageType::Pointer						m_Data;
	typename GradientDirectionContainerType::Pointer	m_DiffusionVectors;

	bool						m_Verbose;
	std::vector<std::string>	m_FileNames;
	unsigned int				m_NumberOfImages ;
	unsigned int				m_NumberOfGradientImages ;
	bool						m_ReadB0;
	double						m_B0 ;
};

} // end namespace DWI.
} // end namespace crl.

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlDWIRaw.txx"
#endif

#endif


