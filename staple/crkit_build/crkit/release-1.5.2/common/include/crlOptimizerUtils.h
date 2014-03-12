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
 * \file	crlOptimizerUtils.h
 *
 * \brief	Declare useful functions for itk optimizers
*************************************************************************************************/

#ifndef CRL_OPTIMIZERS_UTILS
#define CRL_OPTIMIZERS_UTILS

#ifdef _MSC_VER
  #pragma warning( disable: 4996 )
#endif

#include <string>
#include <itkSingleValuedNonLinearOptimizer.h>

namespace crl {

	/**********************************************************************************************//**
	 * \class	OptimizerUtils
	 *
	 * \brief	Class containing useful functions to create and use different optimizers. 
	 *
	 * \section Example
	 * 			
	 * 			\code  
	 * //---------------------------------------------
	 * // Initialize a TCLAP option -a newuoa/bobyqa
	 * //---------------------------------------------
	 * std::vector<std::string> allowedAlgoStr  
	 *    = crl::OptimizerUtils::GetOptimizerSetStringVector(crl::OptimizerUtils::OPTIMIZER_NEWUOA + crl::OptimizerUtils::OPTIMIZER_BOBYQA);  
	 * TCLAP::ValuesConstraint<std::string> allowedAlgoTCLAP( allowedAlgoStr );  
	 * TCLAP::ValueArg<std::string> argOptimiAlgo("a","algo","Optimization algorithm",false,"newuoa",&allowedAlgoTCLAP, cmd);  
	 *
	 * (...)
	 *
	 * //---------------------------------------------
	 * // Creates the optimizer
	 * //---------------------------------------------
	 * crl::OptimizerUtils::OptimizerIdType idAlgo = crl::OptimizerUtils::GetOptimizerIdFromString(argOptimiAlgo.getValue())
	 * crl::OptimizerUtils::GenericOptimizerType::Pointer optimizer = crl::OptimizerUtils::CreateOptimizer(idAlgo);
	 * crl::OptimizerUtils::SetOptimizerParameters(optimizer, 1e-3, 1e-5, 1e-5, 1e-5, 2000, false );
	 * optimizer->SetCostFunction(costFunction);
	 * 			\endcode.  
     *
	 * \author	Benoit Scherrer
	 * \date	July 2010
	*************************************************************************************************/
	class OptimizerUtils
	{
	public:
		typedef itk::SingleValuedNonLinearOptimizer GenericOptimizerType;	/**< Defines the type of a generic optimizer */
		static const unsigned int OPTIMIZER_NB; //= 6;								/**< Current number of supported optimizers*/

		/**********************************************************************************************//**
		 * \enum	OptimizerIdType
		 *
		 * \brief	Internal identifier for each optimizer. 
		*************************************************************************************************/
		typedef enum {
			OPTIMIZER_POWELL = 0x1,		/**< & */
			OPTIMIZER_NEWUOA = 0x2,		/**< & */
			OPTIMIZER_BOBYQA = 0x4,		/**< & */
			OPTIMIZER_FRPR = 0x8,		/**< & */
			OPTIMIZER_LBFGS = 0x10,		/**< & */
			OPTIMIZER_CONJGRAD = 0x20,	/**< & */
		} OptimizerIdType;
		typedef unsigned int OptimizerIdSetType;	/**< Can contain several OPTIMIZER_XX definitions (ex: OPTIMIZER_POWELL+OPTIMIZER_NEWUOA) */

		/**********************************************************************************************//**
		 * \brief	 String name of each supported optimizer
		*************************************************************************************************/
		static const std::string OptimizerStringId[];

		static OptimizerIdType GetOptimizerIdFromString(const std::string& id);
		static std::string GetOptimizerStringFromId(OptimizerIdType optimizerId);
		static std::vector<std::string> GetOptimizerSetStringVector(OptimizerIdSetType optimizerIds);

		static GenericOptimizerType::Pointer CreateOptimizer(OptimizerIdType optimizerId, int nbMaxIter=5000, bool maximize=false );
		static void SetOptimizerParameters(GenericOptimizerType* optimizer, double stepLength, double stepTol, double valTol, double gradTol, int maxNbIter, bool maximize );
		static void SetOptimizerMaxNbIter(GenericOptimizerType* optimizer, int nbMaxIter);
	};


}

#endif
