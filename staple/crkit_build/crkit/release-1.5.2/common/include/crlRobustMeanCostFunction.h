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

#ifndef CRL_ROBUSTMEAN_COSTFUNCTION_H
#define CRL_ROBUSTMEAN_COSTFUNCTION_H

#include "itkSingleValuedCostFunction.h"

using namespace itk;

namespace crl{ 

/**********************************************************************************************//**
 * \class	RobustMeanCostFunction
 *
 * \brief	Cost function for the crl::RobustMeanImageFilter. 
 *
 * \author	Benoit Scherrer
 * \date	April 2010
*************************************************************************************************/
class ITK_EXPORT RobustMeanCostFunction 
  : public SingleValuedCostFunction 
{
public:
	typedef enum
	{
		NORM_L1,
		NORM_L2,
		LORENTZIAN
	} RobustNormType;

  /** Standard class typedefs */
  typedef RobustMeanCostFunction			Self;
  typedef SingleValuedCostFunction			Superclass;
  typedef SmartPointer< Self >				Pointer;
  typedef SmartPointer< const Self >		ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(RobustMeanCostFunction, SingleValuedCostFunction) ;

  /** Method for creation through the object factory. */
  itkNewMacro(Self) ;

  /**  ParametersType typedef.
   *  It defines a position in the optimization search space. */
  typedef Superclass::ParametersType ParametersType ;
  typedef Superclass::DerivativeType DerivativeType ;

  /**  MeasureType typedef.
   *  It defines a type used to return the cost function value. */
  typedef Superclass::MeasureType MeasureType ;

  typedef std::vector<double> ObservationType;

  /** Set/Get the value of the norm type  */
  itkSetMacro(NormType, RobustNormType);
  itkGetMacro(NormType, RobustNormType);

  /** Set/Get the value of parameter T  */
  itkSetMacro(ParamT, double);
  itkGetMacro(ParamT, double);

  void ClearObservations() {m_Observations.clear(); }
  void AddObservation( double value ) {m_Observations.push_back(value); }
  
  virtual unsigned int GetNumberOfParameters() const ;

  /** This method returns the value of the cost function corresponding
    * to the specified parameters. */ 
  virtual MeasureType GetValue( const ParametersType & parameters ) const ;

  /** This method returns the derivative of the cost function corresponding
    * to the specified parameters.   */ 
  virtual void GetDerivative( const ParametersType &,
                                    DerivativeType & ) const ;

  void EstimateT();

protected:
  RobustMeanCostFunction() ;
  virtual ~RobustMeanCostFunction() ;
  virtual void PrintSelf(std::ostream& os, Indent indent) const ;


private:

	ObservationType		m_Observations;
	RobustNormType		m_NormType;
	double				m_ParamT;

} ; // end of class


} // end of namespace crl


#endif

