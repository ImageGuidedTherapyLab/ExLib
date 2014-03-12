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

#ifndef CRL_MFM_TENSOR_ESTIMATE_COSTFUNCTION_H
#define CRL_MFM_TENSOR_ESTIMATE_COSTFUNCTION_H

#include "GlobalDefinitions.h"

#include "itkSingleValuedCostFunction.h"
#include "itkMatrix.h"
#include "itkDiffusionTensor3D.h"
#include "itkSymmetricEigenAnalysis.h"

using namespace itk;

namespace crl{ 
	namespace DWI{ 

		/**********************************************************************************************//**
		 * \class	MFMTensorEstimateCostFunction
		 *
		 * \brief	Tensor estimation cost function. 
		 *
		 * \author	Benoit Scherrer
		 * \date	March 2010
		*************************************************************************************************/
		class ITK_EXPORT MFMTensorEstimateCostFunction : public SingleValuedCostFunction 
		{
		public:
			/** Standard class typedefs */
			typedef MFMTensorEstimateCostFunction			Self;
			typedef SingleValuedCostFunction			Superclass;
			typedef SmartPointer< Self >				Pointer;
			typedef SmartPointer< const Self >			ConstPointer;

			typedef itk::SymmetricEigenAnalysis< VnlTensorType , VnlEigValueType, VnlEigVectorType > SymmetricEigenAnalysisType;

			/** Run-time type information (and related methods). */
			itkTypeMacro(MFMTensorEstimateCostFunction, SingleValuedCostFunction) ;

			/** Method for creation through the object factory. */
			itkNewMacro(Self) ;

			/**  ParametersType typedef.
			*  It defines a position in the optimization search space. */
			typedef float TensorPixelType;
			typedef itk::DiffusionTensor3D< TensorPixelType >    TensorType;
			typedef Superclass::ParametersType ParametersType ;
			typedef Superclass::DerivativeType DerivativeType ;


			/**  MeasureType typedef.
			*  It defines a type used to return the cost function value. */
			typedef Superclass::MeasureType MeasureType ;

			virtual unsigned int GetNumberOfParameters() const ;

			/** This method returns the value of the cost function corresponding
			* to the specified parameters. */ 
			virtual MeasureType GetValue( const ParametersType & parameters ) const ;

			/** This method returns the derivative of the cost function corresponding
			* to the specified parameters.   */ 
			virtual void GetDerivative( const ParametersType &,
				DerivativeType & ) const ;

			/** Set/Get the number of tensors to estimate */
			itkSetMacro( ProcID, int );
			itkGetMacro( ProcID, int );


			/** Set/Get the current use of robust estimation */
			itkSetMacro( UseRobustEstimation, bool );
			itkGetMacro( UseRobustEstimation, bool );

			/** Set/Get the current use of regularization */
			itkSetMacro( UseRegularization, bool );
			itkGetMacro( UseRegularization, bool );

			/** Set/Get the weight for the regularization */
			itkSetMacro( CurrentRegularizationWeight, double );
			itkGetMacro( CurrentRegularizationWeight, double );

			/** Set/Get the maximum b-value index to use (used to do a first pass with only low b-values)*/
			itkSetMacro( CurrentMaxBValueIndexToUse, int );
			itkGetMacro( CurrentMaxBValueIndexToUse, int );

			void			SetInternalVoxelData ( InternalVoxelData *d );
			void			SetAlgoParameters ( crl::DWI::GlobalAlgoParameters *p ) ;

			void			InitAlgo(const	TensorType &tensor, bool forceSingleTensor=false);
			ParametersType	GetOptimizerParametersFromPrecomputeData();
			void			EndAlgoIteration(ParametersType &parameters);

			void			ComputeC_ForRobustEstimation();
			void			ComputeResiduals(std::vector<double> & residuals);

			void			OneTensor2VnlMatrix(VnlTensorType &out, const TensorType &tensor ) const;

			void			ComputeFinalTensors(std::vector<TensorType > &tensor ) const;
			void			ComputeFinalFractions(VnlFractionsType &fractions ) const;

			VnlMatrixType	RotationMatrix(double angle, int r) const;

			int GetNumberOfFractionsToEstimate() const;


			void TestCostFunction(const ParametersType &parameters) const; //tmp!

		protected:
			MFMTensorEstimateCostFunction() ;
			virtual ~MFMTensorEstimateCostFunction() ;
			virtual void PrintSelf(std::ostream& os, Indent indent) const ;

			void ComputeSiLs( void ) const;
			void Compute_ggt( VnlMatrixType &out, const VnlGradientType& g );

			void ApplySoftMaxToFractions() const;

			void InitComputeNumberOfFractionsToEstimate() ;
			void InitComputeNumberOfParameters()  ;

			void Parameters2TensorAndFractions( VnlTensorSetType &tensors,  VnlFractionsType &fractions, const ParametersType &parameters ) const;
			void GetTensorFromSphericalParameterization( VnlTensorType &out, VnlEigValueType &outEigenValues, VnlEigVectorType &outEigenVectors, double alpha, double beta, double gamma, double l1, double l2, double l3  ) const;

		private:
			GlobalAlgoParameters*		m_AlgoParameters;
			InternalVoxelData*			m_InternalVoxelData;

			SymmetricEigenAnalysisType	EigenAnalysis;
			int							m_NbFractionsToEstimate;
			int							m_NbOfParameters;

			int							m_ProcID;
			bool						m_UseRobustEstimation;
			bool						m_UseRegularization;
			double						m_CurrentRegularizationWeight;
			int							m_CurrentMaxBValueIndexToUse;

			// Copy the values of the m_AlgoParameters:
			bool						m_UseIsotropicPart;
			EstimationStrategyType		m_EstimationStrategy;

		} ; // end of class

	} // end of namespace DWI
} // end of namespace crl

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlMFMTensorEstimateCostFunction.txx"
#endif

#endif

