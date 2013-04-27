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

#ifndef GLOBAL_DEFINITIONS_H
#define GLOBAL_DEFINITIONS_H

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_diag_matrix.h"
#include "vnl/vnl_matrix_fixed.h"

#include <itkSimpleFastMutexLock.h>
#include <itkArray.h>

#include "crlOptimizerUtils.h"

namespace crl{
	namespace DWI {

//------------------------------------------------------------
// Some enum for the algorithm parameters
//------------------------------------------------------------

// Type for the SetConstraintsVector method
typedef enum {
  ESTIMATION_SAME_MAJOR_EIGEN_VALUE = 1,
  ESTIMATION_SAME_MINOR_EIGEN_VALUE = 2,
  ESTIMATION_EACH_CYLINDER = 4,
  ESTIMATION_SAME_PLANE = 8,
  ESTIMATION_F_ISO = 16,
  ESTIMATION_F_TENSORS = 32,
  ESTIMATION_RADIALBASELINE = 64,
  ESTIMATION_RADIALISO = 128
} EstimationStrategyValue;
typedef unsigned int EstimationStrategyType;

typedef vnl_matrix<double>					VnlMatrixType;
typedef vnl_matrix<double>					VnlTensorType;
typedef vnl_diag_matrix<double>				VnlEigValueType;
typedef vnl_matrix<double>					VnlEigVectorType;
typedef std::vector<VnlTensorType >			VnlTensorSetType;
typedef std::vector<VnlEigValueType >		VnlEigValueSetType;
typedef std::vector<VnlEigVectorType >		VnlEigVectorSetType;

typedef vnl_vector<double>					VnlFractionsType;

typedef vnl_vector<double>					VnlGradientType;
typedef std::vector<VnlGradientType >		VnlGradientSetType;

//------------------------------------------------------------
// Global parameters to share with the cost function
//------------------------------------------------------------
typedef struct {
	double												m_Threshold; /** Threshold on the reference image data */

	bool												m_RobustEstimation;
	double												m_RobustParam;

	bool												m_UseRicianNoiseModeling;
	double												m_RicianNoiseVariance;

	bool												m_UseRegularization;
	double												m_RegularizationWeight;
	double												m_RegularizationAnisotropicK;

	bool												m_FirstPassWithLowBValuesOnly;
	int													m_LowBValueMaxIndex;	// set in SetGradientImage;

	EstimationStrategyType								m_EstimationStrategy;

	bool												m_InitialOneTensorFromMedINRIA;
	int 												m_VerboseLevel, m_OVerboseLevel;
	unsigned int										m_MaxIter;

	int													m_NumberOfTensors;

	bool												m_OutputInitTensors;

	crl::OptimizerUtils::OptimizerIdType				m_OptimizationAlgo;
	bool												m_OptimizationMultiRestart;

	std::vector<double>									m_BValues;		/// The list of b-values
	VnlGradientSetType									m_g;			/// The list of unit-length gradients
	std::vector< vnl_matrix<double> > 					m_G;			/// The precomputed values of G = g g^t for each g (for the computation of the matrix derivative)

	bool												m_VideoMode;
	bool												m_RandomInit;

	double												m_FastDiffusionFraction;

	double												m_PruneThreshold;

	bool												m_ExportResiduals;
} GlobalAlgoParameters;


//------------------------------------------------------------
// Internal data for each voxel
//------------------------------------------------------------
typedef struct SInternalVoxelData
{
	double 	S0;														/// The S0 value for that voxel
	std::vector< double > 					Si;						/// The list of the intensity values for each gradients

	VnlEigVectorType						OneTensorBasis;			/// The one-tensor basis (for the constrained algo)
	VnlEigValueType							OneTensorEigenValues;	/// The list of the one-tensor eigen values (for the constrained algo)

	int										nbTensorAtVoxel;

	vnl_matrix< double >					SiLs;					/// The current values of the SiLs (without the fraction) for each dir, each tensor
	std::vector< double >					SiL_Isotropic;				/// The precomputed constant value for isotropic part for each bvalue
	VnlFractionsType						Fractions;				/// The current values of the fractions before the softmax transform
	VnlFractionsType						normalizedFractions;	/// The current values of the fractions after the softmax transform

	double							fastDiffusionFraction;


	itk::Array<double>						lowerBound;				/// The lower bound of the parameters (for BOBYQA)
	itk::Array<double>						upperBound;				/// The upper bound of the parameters (for BOBYQA)

	itk::SimpleFastMutexLock				mutexAccessL;
	itk::Array<double>						lastParameters;			/// The last list of parameters, used to relaunch the minimization (mainly for the regularization)
	double									RobustEstimation_C;		/// The current value of the C parameter for robust estimation for that voxel

	VnlTensorSetType						stabilized_L;			/// The value of the L's after the 

	VnlTensorSetType						current_L;				/// The current value of the L's of each tensors for these the iterations
	VnlEigVectorSetType				 		current_EigenVectorsLE;
	VnlEigValueSetType						current_EigenValuesLE;
	VnlTensorSetType						current_RadialMatrix;	/// The current value of the radial matrix
	std::vector<double>						current_RadialCst;		/// One radial cst per tensor
	vnl_matrix< double >					current_GradientRadialProjection;	/// The projection of the gradient on the radial matrix g'Cg
	vnl_matrix< double >					current_GradientTensorProjection;	/// The projection of the gradient on the radial matrix g'Cg

	bool									HasConverged;

	std::vector<SInternalVoxelData *>		Neighbors;				/// The list of pointers to the internal structure of the neighbors

} InternalVoxelData;

}
}

#endif
