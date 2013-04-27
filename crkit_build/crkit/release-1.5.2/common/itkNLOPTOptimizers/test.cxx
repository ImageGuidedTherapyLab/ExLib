/*
 * ITK Wrapper for the NLOPT library - test program
 *
 * Copyright (c) 2010-2011 Children's Hospital Boston.
 * Benoit Scherrer, CRL (Computational Radiology Laboratory), Harvard Medical School
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

#include <iostream>
#include "itkNLOPTOptimizers.h"	
#include "itkSingleValuedCostFunction.h"

using namespace std;

/*------------------------------------------------------------
 The test cost-function
 f(x,y) = cos(x)+(y-1)*(y-1)
------------------------------------------------------------*/
namespace itk {
class ITK_EXPORT MyCostFunction 
	: public itk::SingleValuedCostFunction 
{
public:
  /** Standard class typedefs */
  typedef MyCostFunction					Self;
  typedef SingleValuedCostFunction			Superclass;
  typedef SmartPointer< Self >				Pointer;
  typedef SmartPointer< const Self >		ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(MyCostFunction, SingleValuedCostFunction) ;

  /** Method for creation through the object factory. */
  itkNewMacro(Self) ;

  virtual unsigned int GetNumberOfParameters() const { return 2; }

  /** This method returns the value of the cost function corresponding
    * to the specified parameters. */ 
  virtual MeasureType GetValue( const ParametersType & parameters ) const 
  {
		return (parameters[0]-1)*(parameters[0]-1)+(parameters[1]-2)*(parameters[1]-2);
//		return cos(parameters[0]) + (parameters[1]-1)*(parameters[1]-1);
  }

  /** This method returns the derivative of the cost function corresponding
    * to the specified parameters. 
	NOT USED IN NEWUOA */ 
  virtual void GetDerivative( const ParametersType &,
                                    DerivativeType & ) const 
  {
  }


protected:
	MyCostFunction() {}
	virtual ~MyCostFunction() {}
  virtual void PrintSelf(std::ostream& os, Indent indent) const 
  {}

private:

} ; // end of class
}

/*------------------------------------------------------------
Main entry point of the program
------------------------------------------------------------*/
int main(int argc, char **argv)
{
	try {
		// Create the cost function
		itk::MyCostFunction::Pointer costFun = itk::MyCostFunction::New();

		// Create the optimizer
		itk::NLOPTOptimizers::Pointer optimizer = itk::NLOPTOptimizers::New();
		optimizer->SetCostFunction(costFun);

		// Set the parameters
		optimizer->SetMaxEval(500);
		optimizer->SetAlgorithm(itk::NLOPTOptimizers::NLOPT_LN_NEWUOA);

		// Set the initial position
		itk::NLOPTOptimizers::ParametersType pinit(2);
		pinit[0]=pinit[1] = 4;
		optimizer->SetInitialPosition(pinit);

		// Verbose level: 0,1,2,3
		optimizer->SetVerboseLevel(1);

		optimizer->StartOptimization();

		std::cout<<"Return code "<<(int)optimizer->GetErrorCode();
		std::cout<<" = "<<optimizer->GetErrorCodeDescription()<<endl;

		if ( !optimizer->isSuccessful() ) 
		{
			std::cout<<"Fatal error"<<std::endl<<std::endl;
		}
		else
		{
			std::cout<<"Position = "<<optimizer->GetCurrentPosition()<<endl;
			std::cout<<"Value = "<<optimizer->GetCurrentCost()<<endl;
			std::cout<<"Time = "<<optimizer->GetTimeForOptimization()<<" seconds"<<endl<<std::endl;
			// no number of iterations in NLOPT?
			//std::cout<<"Number of iterations = "<<optimizer->GetNumberOfIterations()<<endl;
		}

	}

	/*-------------------------------------
	If catched an exception, show a message
	-------------------------------------*/
	catch (itk::ExceptionObject& e)
	{
		cout<< "  Exception: "<<endl<< e.GetDescription();
	}

  return 0;
}
