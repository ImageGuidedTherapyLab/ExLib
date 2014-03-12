/**********************************************************************************************//**
 * \file	crlGeneralizedInstantiatorMacros.h
 *
 * \brief	Declare some useful macros to simply the writing of classes/base classes for the 
 *			crlGeneralizedInstantiator. 
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


#ifndef CRL_GENERALIZED_INSTANTIATOR_MACROS_H
#define CRL_GENERALIZED_INSTANTIATOR_MACROS_H

/**********************************************************************************************//**
 * \brief	Macro used to simplify the code of a base class for the instanciator
 *
 * \remarks	Benoit Scherrer, June 2010. 
 *
 * \param	className	Name of the class. 
*************************************************************************************************/
#define CRL_INSTANCIATOR_BASE(className) \
public:\
	/** Standard class typedefs. */\
	typedef className						Self;\
	typedef itk::SmartPointer<Self>			Pointer;\
	typedef itk::SmartPointer<const Self>	ConstPointer;\
\
	/** Run-time type information (and related methods). */ \
	itkTypeMacro(className, LightObject); \
\
private:\
	className(const Self &);		/** purposely not implemented */  \
	void operator=(const Self &);			/** purposely not implemented */ 

/**********************************************************************************************//**
 * \def	CRL_INSTANCIATOR_CLASS(superclassName, className) public: typedef className Self
 *
 * \brief	Macro used to simplify the code of a class for the instanciator
 *
 * \remarks	Benoit Scherrer, June 2010. 
 *
 * \param	superclassName	Name of the superclass. 
 * \param	className		Name of the class. 
*************************************************************************************************/
#define CRL_INSTANCIATOR_CLASS(superclassName, className) \
public:\
	/** Standard class typedefs. */\
	typedef className						Self;\
	typedef	superclassName					Superclass;\
	typedef itk::SmartPointer<Self>			Pointer;\
	typedef itk::SmartPointer<const Self>	ConstPointer;\
\
	/** Method for creation through the object factory. */ \
	itkNewMacro(Self); \
\
private:\
	className(const Self &);		/** purposely not implemented */  \
	void operator=(const Self &);	/** purposely not implemented */


#endif
