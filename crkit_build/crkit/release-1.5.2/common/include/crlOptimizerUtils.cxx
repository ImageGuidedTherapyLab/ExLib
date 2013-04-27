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
 * \file	crlOptimizerUtils.cxx
 *
 * \brief	Useful function to manipulate optimizers 
*************************************************************************************************/


#include "crlOptimizerUtils.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <itkExceptionObject.h>
#include "itkConjugateGradientOptimizer.h"
#include "itkLBFGSOptimizer.h"
#include "itkFRPROptimizer.h"
#include "itkPowellOptimizer.h"
#include "itkNLOPTOptimizers.h"

namespace crl {

	const unsigned int OptimizerUtils::OPTIMIZER_NB = 6;
	const std::string OptimizerUtils::OptimizerStringId[] =
		{ 
			"powell",
			"newuoa",
			"bobyqa",
			"frpr",
			"lbfgs",
			"conjgrad"
		};

	/**********************************************************************************************//**
	 * \fn	OptimizerUtils::OptimizerIdType OptimizerUtils::GetOptimizerIdFromString(const std::string& id)
	 *
	 * \brief	Gets the optimizer identifier OptimizerIdType from a string. 
	 *
	 * \author	Benoit Scherrer
	 * \date	July 2010
	 *
	 * \exception	itk::ExceptionObject	Thrown if unknown string identifier. 
	 *
	 * \param	strId	The identifier. 
	 *
	 * \return	The optimizer identifier from string. 
	 *
	 * \sa	OptimizerIdType
	*************************************************************************************************/
	OptimizerUtils::OptimizerIdType OptimizerUtils::GetOptimizerIdFromString(const std::string& strId)
	{
		int id=1;

		for ( unsigned int i=0; i<OPTIMIZER_NB ; i++ )
		{
			if ( OptimizerStringId[i] == strId )
				return (OptimizerIdType)id;
			id=id<<1;
		}

		std::stringstream s;
		s << "Invalid optimizer name: <"<<id<<">"<<std::endl;
		throw itk::ExceptionObject(__FILE__,__LINE__, s.str().c_str(), "GetOptimizerIdFromString" );
	}

	/**********************************************************************************************//**
	 * \fn	std::string OptimizerUtils::GetOptimizerStringFromId(OptimizerIdType optimizerId)
	 *
	 * \brief	Gets the optimizer string description from an internal identifier. 
	 *
	 * \author	Benoit Scherrer
	 * \date	July 2010
	 *
	 * \exception	itk::ExceptionObject	Thrown when exception. 
	 *
	 * \param	optimizerId	Identifier for the optimizer. 
	 *
	 * \return	The optimizer string from identifier. 
	*************************************************************************************************/
	std::string OptimizerUtils::GetOptimizerStringFromId(OptimizerIdType optimizerId)
	{
		int id=1;
		for ( unsigned int i=0; i<OPTIMIZER_NB ; i++ )
		{	
			if ( id == (int)optimizerId )
				return OptimizerStringId[i];
			id =id<<1;
		}

		std::stringstream s;
		s << "Invalid optimizer internal identifier: <"<<(int)optimizerId<<">"<<std::endl;
		throw itk::ExceptionObject(__FILE__,__LINE__, s.str().c_str(), "GetOptimizerStringFromId" );
	}

	/**********************************************************************************************//**
	 * \fn	std::vector<std::string> OptimizerUtils::GetOptimizerSetStringVector(OptimizerIdSetType optimizerIds )
	 *
	 * \brief	Creates a std::vector of strings containing the names of the specified optimizers
	 * 			
	 * 			It is useful to generate a user-friendly TCLAP argument. For example:  
	 * 			
	 * 			\code  
	 *  std::vector<std::string> allowedAlgoStr 
	 *        = crl::OptimizerUtils::GetOptimizerSetStringVector(crl::OptimizerUtils::OPTIMIZER_NEWUOA + crl::OptimizerUtils::OPTIMIZER_BOBYQA);
	 * 	TCLAP::ValuesConstraint<std::string> allowedAlgoTCLAP( allowedAlgoStr );
	 * 	TCLAP::ValueArg<std::string> argOptimiAlgo("a","algo","Optimization algorithm",false,"newuoa",&allowedAlgoTCLAP, cmd);
	 * 			\endcode.  
	 *
	 * \author	Benoit Scherrer
	 * \date	July 2010
	 *
	 * \param	optimizerIds	A combination of the OPTIMIZER_XXX identifier types. 
	 *
	 * \return	the vector of strings. 
	*************************************************************************************************/
	std::vector<std::string> OptimizerUtils::GetOptimizerSetStringVector(OptimizerIdSetType optimizerIds )
	{
		unsigned int id = 1;
		std::vector<std::string> res;

		for ( unsigned int i=0; i<OPTIMIZER_NB ; i++ )
		{
			if ( optimizerIds&id ) 
				res.push_back(OptimizerStringId[i]);

			id = id<<1;
		}

		return res;
	}

	/**********************************************************************************************//**
	 * \fn	OptimizerUtils::GenericOptimizerType::Pointer OptimizerUtils::CreateOptimizer(OptimizerIdType optimizerId,
	 * 		int nbMaxIter, bool maximize )
	 *
	 * \brief	Creates an optimizer given its internal \a OptimizerIdType identifier. 
	 *
	 * \author	Benoit Scherrer
	 * \date	July 2010
	 *
	 * \exception	itk::ExceptionObject	Thrown when the identifier is invalid. 
	 *
	 * \param	optimizerId	Identifier for the optimizer. 
	 * \param	nbMaxIter	The maximum number of iterations. 
	 * \param	maximize	true to maximize (else minimize). 
	 *
	 * \return	The optimizer. 
	 *
	 * \sa	GetOptimizerIdFromString
	 * \sa  SetOptimizerParameters
	*************************************************************************************************/
	OptimizerUtils::GenericOptimizerType::Pointer OptimizerUtils::CreateOptimizer(OptimizerIdType optimizerId, int nbMaxIter, bool maximize )
	{
		GenericOptimizerType::Pointer optimizer;

		//---------------------------------------------------------------
		// LBFGS
		//---------------------------------------------------------------
		if ( optimizerId==OPTIMIZER_LBFGS )
		{
			typedef itk::LBFGSOptimizer							LBFGSOptimizerType;
			typedef LBFGSOptimizerType::InternalOptimizerType	vnlLBFGSOptimizerType;
			LBFGSOptimizerType::Pointer LBFGSoptimizer = LBFGSOptimizerType::New();
			optimizer = LBFGSoptimizer;
		}

		//---------------------------------------------------------------
		// FRPR
		//---------------------------------------------------------------
		else if ( optimizerId==OPTIMIZER_FRPR )
		{
			typedef itk::FRPROptimizer	FRPROptimizerType;
			FRPROptimizerType::Pointer FRPRoptimizer = FRPROptimizerType::New();
			optimizer=FRPRoptimizer;
		}

		//---------------------------------------------------------------
		// CONJUGATE GRADIENT DESCENT 
		//---------------------------------------------------------------
		else if ( optimizerId==OPTIMIZER_CONJGRAD )
		{
			typedef itk::ConjugateGradientOptimizer		CONJGRADOptimizerType;
			typedef CONJGRADOptimizerType::InternalOptimizerType	vnlCONJGRADOptimizerType;
			CONJGRADOptimizerType::Pointer CONJGRADoptimizer = CONJGRADOptimizerType::New();
			optimizer = CONJGRADoptimizer;
		}

		//---------------------------------------------------------------
		// POWELL
		//---------------------------------------------------------------
		else if ( optimizerId==OPTIMIZER_POWELL)
		{
			typedef itk::PowellOptimizer				POWELLOptimizerType;
			POWELLOptimizerType::Pointer POWELLoptimizer = POWELLOptimizerType::New();
			optimizer = POWELLoptimizer;
		}

		//---------------------------------------------------------------
		// NEWUOA - from NLOPT
		//---------------------------------------------------------------
		else if ( optimizerId==OPTIMIZER_NEWUOA )
		{
			typedef itk::NLOPTOptimizers				NLOPTOptimizerType;
			NLOPTOptimizerType::Pointer NLOPTOptimiz = NLOPTOptimizerType::New();
			optimizer = NLOPTOptimiz;
			NLOPTOptimiz->SetAlgorithm(NLOPTOptimizerType::NLOPT_LN_NEWUOA);
		}

		//---------------------------------------------------------------
		// BOBYQA - from NLOPT
		//---------------------------------------------------------------
		else if ( optimizerId==OPTIMIZER_BOBYQA )
		{
			typedef itk::NLOPTOptimizers NLOPTOptimizerType;
			NLOPTOptimizerType::Pointer NLOPTOptimiz = NLOPTOptimizerType::New();
			optimizer = NLOPTOptimiz;
			NLOPTOptimiz->SetAlgorithm(NLOPTOptimizerType::NLOPT_LN_BOBYQA);
		}
		else
		{
			std::stringstream s;
			s << "Error. Unkown optimizer internal identifier: <"<<optimizerId<<">"<<std::endl;
			throw itk::ExceptionObject(__FILE__,__LINE__, s.str().c_str(), "CreateOptimizer" );
		}

		SetOptimizerParameters(optimizer, 1e-3, 1e-5, 1e-5, 1e-5, nbMaxIter, maximize );
		return optimizer;
	}

	/**********************************************************************************************//**
	 * \fn	void OptimizerUtils::SetOptimizerMaxNbIter(GenericOptimizerType* optimizer, int nbMaxIter)
	 *
	 * \brief	Sets the number of max iteration for an optimizer. 
	 *
	 * \author	Benoit Scherrer
	 * \date	July 2010
	 *
	 * \exception	itk::ExceptionObject	Thrown when the type of the optimizer is not recognized. 
	 *
	 * \param [in,out]	optimizer	If non-null, the optimizer. 
	 * \param	nbMaxIter			The nb maximum iterator. 
	 *
	*************************************************************************************************/
	void OptimizerUtils::SetOptimizerMaxNbIter(GenericOptimizerType* optimizer, int nbMaxIter)
	{
		//---------------------------------------------------------------
		// POWELL
		//---------------------------------------------------------------
		typedef itk::PowellOptimizer				POWELLOptimizerType;
		POWELLOptimizerType *powell = dynamic_cast<POWELLOptimizerType *>(optimizer);
		if ( powell!=NULL )
		{
			powell->SetMaximumIteration(nbMaxIter);
			return;
		}

		//---------------------------------------------------------------
		// NLOPT based optimizers
		//---------------------------------------------------------------
		typedef itk::NLOPTOptimizers				NLOPTOptimizerType;
		NLOPTOptimizerType *nlopt = dynamic_cast<NLOPTOptimizerType *>(optimizer);
		if ( nlopt!=NULL )
		{
			nlopt->SetMaxEval(nbMaxIter);
			return;
		}

		//---------------------------------------------------------------
		// LBFGS
		//---------------------------------------------------------------
		typedef itk::LBFGSOptimizer					LBFGSOptimizerType;
		LBFGSOptimizerType *lbfgs = dynamic_cast<LBFGSOptimizerType *>(optimizer);
		if ( lbfgs!=NULL )
		{
			lbfgs->GetOptimizer()->set_max_function_evals( nbMaxIter );
			return;
		}

		//---------------------------------------------------------------
		// FRPR
		//---------------------------------------------------------------
		typedef itk::FRPROptimizer					FRPROptimizerType;
		FRPROptimizerType *frpr = dynamic_cast<FRPROptimizerType *>(optimizer);
		if ( frpr!=NULL )
		{
			frpr->SetMaximumIteration(nbMaxIter);
			return;
		}

		//---------------------------------------------------------------
		// CONJUGATE GRADIENT DESCENT 
		//---------------------------------------------------------------
		typedef itk::ConjugateGradientOptimizer		CONJGRADOptimizerType;
		CONJGRADOptimizerType *conjgrad = dynamic_cast<CONJGRADOptimizerType *>(optimizer);
		if ( conjgrad!=NULL )
		{
			conjgrad->GetOptimizer()->set_max_function_evals( nbMaxIter );
			return;
		}

		//---------------------------------------------------------------
		// If we are there there is an error
		//---------------------------------------------------------------
		std::stringstream s;
		s << "Error. Unkown optimizer class: <"<<optimizer->GetNameOfClass ()<<">"<<std::endl;
		throw itk::ExceptionObject(__FILE__,__LINE__, s.str().c_str(), "SetOptimizerMaxNbIter" );
	}

	/**********************************************************************************************//**
	 * \fn	void OptimizerUtils::SetOptimizerParameters(GenericOptimizerType* optimizer,
	 * 		double stepLength, double stepTol, double valTol, double gradTol, int maxNbIter,
	 * 		bool maximize )
	 *
	 * \brief	Sets the most common parameters for an optimizer.
	 * 			
	 * 			Not all optimizer use all the most common parameters:
	 * 			-  OPTIMIZER_POWELL : stepLength, stepTol, valTol, maxNbIter,  maximize
	 * 			-  NLOPT based (OPTIMIZER_NEWUOA, OPTIMIZER_BOBYQA) : stepTol, valTol, maxNbIter,
	 * 			maximize
	 * 			-  OPTIMIZER_FRPR : stepTol, valTol, maxNbIter,  maximize
	 * 			-  OPTIMIZER_LBFGS : stepLength, stepTol, valTol, gradTol, maxNbIter,  maximize
	 * 			-  OPTIMIZER_CONJGRAD : stepLength, stepTol, valTol, gradTol, maxNbIter,  maximize. 
	 *
	 * \author	Benoit Scherrer
	 * \date	July 2010
	 *
	 * \exception	itk::ExceptionObject	Thrown when the type of the optimizer is not recognized. 
	 *
	 * \param [in,out]	optimizer	If non-null, the optimizer. 
	 * \param	stepLength			The length of each step. 
	 * \param	stepTol				The search space (step) tolerance. 
	 * \param	valTol				The function value tolerance. 
	 * \param	gradTol				The gradient magnitude tolerance (only for gradient-based
	 * 								methods) 
	 * \param	maxNbIter			The maximum number of iterations. 
	 * \param	maximize			true to maximize. 
	 *
	*************************************************************************************************/
	void OptimizerUtils::SetOptimizerParameters(GenericOptimizerType* optimizer, 
		 double stepLength, double stepTol, double valTol, double gradTol, int maxNbIter, bool maximize )
	{
		typedef itk::PowellOptimizer				POWELLOptimizerType;
		POWELLOptimizerType *powell = dynamic_cast<POWELLOptimizerType *>(optimizer);
		if ( powell!=NULL )
		{
			powell->SetStepLength(stepLength);
			powell->SetStepTolerance(stepTol);
			powell->SetValueTolerance(valTol);
			powell->SetMaximumIteration(maxNbIter);
			powell->SetMaximize(maximize);
			return;
		}

		typedef itk::NLOPTOptimizers NLOPTOptimizerType;
		NLOPTOptimizerType *nlopt = dynamic_cast<NLOPTOptimizerType *>(optimizer);
		if ( nlopt!=NULL )
		{
			nlopt->SetXTolAbs(stepTol);
			nlopt->SetFTolAbs(valTol);
			nlopt->SetMaxEval(maxNbIter);
			nlopt->SetMaximize(maximize);
			return;
		}

		typedef itk::LBFGSOptimizer					LBFGSOptimizerType;
		LBFGSOptimizerType *lbfgs = dynamic_cast<LBFGSOptimizerType *>(optimizer);
		if ( lbfgs!=NULL )
		{
			typedef LBFGSOptimizerType::InternalOptimizerType	vnlLBFGSOptimizerType;
			vnlLBFGSOptimizerType	* vnlOptimizer = lbfgs->GetOptimizer();
			vnlOptimizer->set_f_tolerance( valTol );
			vnlOptimizer->set_g_tolerance( gradTol );
			vnlOptimizer->set_x_tolerance( stepTol );
			vnlOptimizer->set_epsilon_function( stepLength );
			vnlOptimizer->set_max_function_evals( maxNbIter );
			lbfgs->SetDefaultStepLength(stepLength);
			lbfgs->SetGradientConvergenceTolerance(gradTol);

			if ( !maximize ) lbfgs->MinimizeOn();
			else lbfgs->MaximizeOn();
			return;
		}

		typedef itk::FRPROptimizer FRPROptimizerType;
		FRPROptimizerType *frpr = dynamic_cast<FRPROptimizerType *>(optimizer);
		if ( frpr!=NULL )
		{
			frpr->SetStepLength(stepLength);
			frpr->SetStepTolerance (stepTol);
			frpr->SetValueTolerance(valTol); 
			frpr->SetMaximumIteration(maxNbIter);
			frpr->SetMaximize(maximize);
			return;
		}

		typedef itk::ConjugateGradientOptimizer		CONJGRADOptimizerType;
		CONJGRADOptimizerType *conjgrad = dynamic_cast<CONJGRADOptimizerType *>(optimizer);
		if ( conjgrad!=NULL )
		{
			typedef CONJGRADOptimizerType::InternalOptimizerType	vnlCONJGRADOptimizerType;
			vnlCONJGRADOptimizerType	* vnlOptimizer = conjgrad->GetOptimizer();
			
			vnlOptimizer->set_f_tolerance( valTol );
			vnlOptimizer->set_g_tolerance( gradTol );
			vnlOptimizer->set_x_tolerance( stepTol );
			vnlOptimizer->set_epsilon_function( stepLength );
			vnlOptimizer->set_max_function_evals( maxNbIter );

			if ( !maximize ) conjgrad->MinimizeOn();
			else conjgrad->MaximizeOn();

			return;
		}

		//---------------------------------------------------------------
		// If we are there there is an error
		//---------------------------------------------------------------
		std::stringstream s;
		s << "Error. Unkown optimizer class: <"<<optimizer->GetNameOfClass ()<<">"<<std::endl;
		throw itk::ExceptionObject(__FILE__,__LINE__, s.str().c_str(), "SetOptimizerParameters" );
	}

	//==============================================================
	// The following functions may be defined in the future to
	// set additional specific parameters for each optimizer
	//==============================================================
	//void OptimizerUtils::SetPowellAdditionalParameters(GenericOptimizerType* optimizer, 
	//	)
	//{
	//	typedef itk::PowellOptimizer				POWELLOptimizerType;
	//	POWELLOptimizerType *powell = dynamic_cast<POWELLOptimizerType *>(optimizer);
	//	if ( powell!=NULL )
	//	{
	//	}
	//}

	//void OptimizerUtils::SetNLOPTAdditionalParameters(GenericOptimizerType* optimizer, 
	//	)
	//{
	//	typedef itk::NLOPTOptimizers NLOPTOptimizerType;
	//	NLOPTOptimizerType *nlopt = dynamic_cast<NLOPTOptimizerType *>(optimizer);
	//	if ( nlopt!=NULL )
	//	{

	//	}
	//}

	//void OptimizerUtils::SetLBFGSAdditionalParameters(GenericOptimizerType* optimizer, 
	// )
	//{
	//	typedef itk::LBFGSOptimizer					LBFGSOptimizerType;
	//	LBFGSOptimizerType *lbfgs = dynamic_cast<LBFGSOptimizerType *>(optimizer);
	//	if ( lbfgs!=NULL )
	//	{
	//		//LBFGSoptimizer->SetLineSearchAccuracy(0.8);
	//		//vnlOptimizer->set_check_derivatives(1 );
	//	}
	//}

	//void OptimizerUtils::SetFRPRAdditionalParameters(GenericOptimizerType* optimizer,
	//	)
	//{
	//	typedef itk::FRPROptimizer					FRPROptimizerType;
	//	FRPROptimizerType *frpr = dynamic_cast<FRPROptimizerType *>(optimizer);
	//	if ( frpr!=NULL )
	//	{
	//		//FRPRoptimizer->SetMaximumLineIteration();		
	//		//FRPRoptimizer->SetUseUnitLengthGradient(true);
	//	}
	//}

	//void OptimizerUtils::SetCONJGRADAdditionalParameters(GenericOptimizerType* optimizer, 
	// )
	//{
	//	typedef itk::ConjugateGradientOptimizer		CONJGRADOptimizerType;
	//	CONJGRADOptimizerType *conjgrad = dynamic_cast<CONJGRADOptimizerType *>(optimizer);
	//	if ( conjgrad!=NULL )
	//	{
	//		typedef CONJGRADOptimizerType::InternalOptimizerType	vnlCONJGRADOptimizerType;
	//		vnlCONJGRADOptimizerType	* vnlOptimizer = conjgrad->GetOptimizer();
	//	}
	//}



}
