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

#ifndef CRL_DWI_GRADIENTS_H
#define CRL_DWI_GRADIENTS_H

#include <iostream>
using namespace std;

#include "crlDWICommon.h"


namespace crl {
namespace DWI {

/**********************************************************************************************//**
 * \class	DWIGradients
 *
 * \brief	Class to create DWI gradient schemes.
 *
 * \author	Benoit Scherrer
 * \date	August 2010
*************************************************************************************************/
class DWIGradients 
{	
public:
	typedef vnl_vector_fixed< double,3 > GradientDirectionType;
	typedef itk::VectorContainer< unsigned int, GradientDirectionType > GradientSchemeITKType;
	typedef std::vector< GradientDirectionType > GradientSchemeStdType;
	
	DWIGradients( void );
	~DWIGradients();

	void						SetVerbose(bool b) { m_Verbose = b; }

	static void					GetStringDescriptionHelp(std::string& shortHelp, std::string& longHelp);
	void						SetFromString(const std::string& gradientScheme);

	void						SetNominalB0(double _B0);
	double						GetNominalB0() { return m_NominalB0;}

	void						AddOneShell(int nbDirections, double _B0);
	void						AddSpiralShell(int nbDirections, int NbBValues, double MinB, double MaxB);
	void						AddCubeCorners(int nbRepeatCubeN2, int nbRepeatCubeN3 );
	void						AddRandomGradientsOnCubeEdges(int nbGradients);

	void						AddNullGradients(int nb);

	void						ClearAllGradientDirections();
	GradientSchemeITKType*		GetSchemeAsITKContainer() const;
	GradientSchemeStdType		GetSchemeAsStdVector() const;

	void						ExtractNormalizedGradientsAndBValues(GradientSchemeStdType& grads, std::vector<double>& bvalues, bool includeB0=true );

protected:
	const double*				GetGradientList(int nbDirections);

protected:
	GradientSchemeITKType::Pointer		m_DiffusionVectors;

	double			m_NominalB0 ;
	bool			m_Verbose;
};

} // end namespace DWI.
} // end namespace crl.

#endif


