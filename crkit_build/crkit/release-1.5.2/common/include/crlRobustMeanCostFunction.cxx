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

#include <algorithm>
#include "crlRobustMeanCostFunction.h"

namespace crl{ 

/**********************************************************************************************//**
 * \fn	RobustMeanCostFunction::RobustMeanCostFunction()
 *
 * \brief	Default constructor. 
 *
 * \author	Benoit Scherrer
 * \date	April 2010
*************************************************************************************************/
RobustMeanCostFunction::RobustMeanCostFunction()
{
	m_ParamT = 20.0;
	m_NormType = LORENTZIAN;
}

/**********************************************************************************************//**
 * \fn	RobustMeanCostFunction::~RobustMeanCostFunction()
 *
 * \brief	Destructor. 
 *
 * \author	Benoit Scherrer
 * \date	April 2010
*************************************************************************************************/
RobustMeanCostFunction::~RobustMeanCostFunction()
{
}

/**********************************************************************************************//**
 * \fn	void RobustMeanCostFunction::PrintSelf(std::ostream& os, Indent indent) const
 *
 * \brief	Print information about the object, including superclasses. 
 *
 * \author	Benoit Scherrer
 * \date	April 2010
 *
 * \param [in,out]	os	The output stream. 
 * \param	indent		The indent. 
*************************************************************************************************/
void RobustMeanCostFunction::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "TO DO  " << std::endl;
}

/**********************************************************************************************//**
 * \fn	unsigned int RobustMeanCostFunction::GetNumberOfParameters() const
 *
 * \brief	Gets the number of parameters of this cost function (equal to 1). 
 *
 * \author	Benoit Scherrer
 * \date	April 2010
 *
 * \return	The number of parameters. 
*************************************************************************************************/
unsigned int RobustMeanCostFunction::GetNumberOfParameters() const 
{
	return 1 ;
}

/**********************************************************************************************//**
 * \fn	void RobustMeanCostFunction::EstimateT()
 *
 * \brief	Estimate the constant for the robust M-Estimator (equal to 1.48*median||errs|| )
 *
 * \author	Benoit Scherrer
 * \date	April 2010
*************************************************************************************************/
void RobustMeanCostFunction::EstimateT()
{
	/*-------------------------------------
	 Compute the mean value
	-------------------------------------*/
	double mean_value=0;
	for ( unsigned int i=0; i<m_Observations.size() ; i++ )
		mean_value += m_Observations[i];
	mean_value = mean_value / ((double)m_Observations.size());

	/*-------------------------------------
	 Compute the list of abs errors to the mean value
	-------------------------------------*/
	ObservationType	listErrors;
	listErrors.resize(m_Observations.size());

	for( unsigned int i = 0; i<m_Observations.size() ; i++ )
	{
		double error = mean_value -m_Observations[i];
		listErrors[i] = (error>0) ? error:(-error);
	}

	/*-------------------------------------
	 Sort these errors and take the median
	-------------------------------------*/
	std::sort(listErrors.begin(), listErrors.end());
	double median = listErrors[(int)(listErrors.size()/2)];

	m_ParamT = 1.48*median;
	
#ifdef OUTPUT_IN_DEBUG_FILES
	std::fstream ofile;
	char szBuffer[256];
 	sprintf(szBuffer,"ROBUST_p%d.txt",m_ProcID);
	ofile.open ( szBuffer, std::fstream::out | std::fstream::app);
	if  ( !ofile.is_open() ) 
	{
		cout<<"Cannot open the file"<<szBuffer<<endl;
		return;
	}
	for ( unsigned int i=0 ; i<listErrors.size() ; i++ )
		ofile<<listErrors[i]<<" - ";
	ofile<<endl<<"MEDIAN = "<<median<<endl<<endl;
	ofile.close();
#endif

	listErrors.clear();
}

/**********************************************************************************************//**
 * \fn	RobustMeanCostFunction::MeasureType RobustMeanCostFunction::GetValue(const ParametersType &parameters) const
 *
 * \brief	Gets the value of the cost function. 
 *
 * \author	Benoit Scherrer
 * \date	April 2010
 *
 * \param	parameters	Options for controlling the operation. 
 *
 * \return	The value. 
*************************************************************************************************/
RobustMeanCostFunction::MeasureType
RobustMeanCostFunction::GetValue(const ParametersType &parameters) const
{
	double value = 0.0 ;
	double z, x, T2;

	z = parameters[0];
	T2 = m_ParamT*m_ParamT;				// precompute the square

	for ( unsigned int k=0; k<m_Observations.size() ; k++ )
	{
		x = z - m_Observations[k];

		switch ( m_NormType )
		{
		case NORM_L1:
			value += ((x>0)?x:-x);
			break;

		case NORM_L2:
			value += x*x;
			break;

		case LORENTZIAN:
			value += log(1 + 0.5 * x*x/T2 );
			break;
		}
	}

	return value ;
}

/**********************************************************************************************//**
 * \fn	void RobustMeanCostFunction::GetDerivative( const ParametersType &parameters,
 * 		DerivativeType &derivative ) const
 *
 * \brief	Gets the derivative of the cost function. 
 *
 * \author	Benoit Scherrer
 * \date	April 2010
 *
 * \param	parameters			Input position where to compute the derivate. 
 * \param [in,out]	derivative	The derivative. 
*************************************************************************************************/
void RobustMeanCostFunction::GetDerivative( const ParametersType &parameters, DerivativeType &derivative ) const 
{
	double value = 0.0 ;
	double z, x, T2;

	z = parameters[0];
	T2 = m_ParamT*m_ParamT;				// precompute the square

	for ( unsigned int k=0; k<m_Observations.size() ; k++ )
	{
		x = z - m_Observations[k];

		switch ( m_NormType )
		{
		case NORM_L1:
			value += ((x>0)?1:-1);
			break;

		case NORM_L2:
			value += 2*x;
			break;

		case LORENTZIAN:
			value += (2*x)/(2*T2 + x*x ) ;
			break;
		}
	}

	derivative.SetSize(1);
	derivative[0]=value;
}


} // end of namespace crl

