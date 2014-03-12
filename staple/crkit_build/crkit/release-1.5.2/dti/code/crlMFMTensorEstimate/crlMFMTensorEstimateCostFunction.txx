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

#ifndef CRL_MFM_TENSOR_ESTIMATE_COSTFUNCTION_TXX
#define CRL_MFM_TENSOR_ESTIMATE_COSTFUNCTION_TXX

#include <fstream>
#include <limits>
#include <vnl/vnl_bessel.h>

#include <itkSimpleFastMutexLock.h>
#include "itkGradientMagnitudeImageFilter.h" //tmp

#include "crlMFMTensorEstimateCostFunction.h"
#include "crlCommonMath.h"
#include "crlDWICommon.h"
#include "crlFastOps.h"

using namespace std;

namespace crl{ 
namespace DWI{ 

#ifndef my_MIN
#define my_MIN(a,b) (a<b?a:b)
#endif
#ifndef my_MAX
#define my_MAX(a,b) (a>b?a:b)
#endif

#define safe_log(x) ((x>1e-15?log(x):log(1e-15)))

#ifndef _ASSERT
#define _ASSERT(a)
#endif

//#define OUTPUT_IN_DEBUG_FILES 1
#define NEW_FASTSLOWMODEL 0

template <class T>
 T vnl_trace(vnl_matrix<T> const& M)
 {
   T sum(0);
   const unsigned int N = M.rows()<M.cols() ? M.rows() : M.cols();
   for (unsigned int i=0; i<N; ++i)
     sum += M(i, i);
   return sum;
 }

/**********************************************************************************************//**
 * \fn	MFMTensorEstimateCostFunction::MFMTensorEstimateCostFunction()
 *
 * \brief	Default constructor. 
 *
 * \author	Benoit Scherrer
 * \date	April 2010
*************************************************************************************************/
MFMTensorEstimateCostFunction::MFMTensorEstimateCostFunction():
EigenAnalysis(3)
{
	m_InternalVoxelData = NULL;
	m_AlgoParameters = NULL;

	m_NbFractionsToEstimate = -1;
	m_NbOfParameters = -1;

	m_ProcID = 0;
	m_UseRobustEstimation = false;
	m_UseRegularization = false;
	m_CurrentRegularizationWeight = 0.0;
	m_CurrentMaxBValueIndexToUse = 0;

	m_UseIsotropicPart = false;
	m_EstimationStrategy = 0;

}

/**********************************************************************************************//**
 * \fn	MFMTensorEstimateCostFunction::~MFMTensorEstimateCostFunction()
 *
 * \brief	Destructor. 
 *
 * \author	Benoit Scherrer
 * \date	April 2010
*************************************************************************************************/
MFMTensorEstimateCostFunction::~MFMTensorEstimateCostFunction()
{
}

void MFMTensorEstimateCostFunction::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "TO DO  " << std::endl;
}

/**********************************************************************************************//**
 * \fn	unsigned int MFMTensorEstimateCostFunction::GetNumberOfParameters() const
 *
 * \brief	Gets the number of parameters.  
 * 			
 * \author	Benoit Scherrer
 * \date	February 2010
 *
 * \exception	itk::ExceptionObject	Thrown when exception. 
 *
 * \return	The number of parameters. 
 *
 * \sa	GetNumberOfFractionsToEstimate
*************************************************************************************************/
unsigned int MFMTensorEstimateCostFunction::GetNumberOfParameters() const 
{
	if ( m_NbOfParameters<0 ) 
		throw itk::ExceptionObject(__FILE__, __LINE__, "Error. The number of parameters to estimate has not been computed yet.","MFMTensorEstimateCostFunction::GetNumberOfFractionsToEstimate()");

	else return (unsigned int)m_NbOfParameters;
}

/**********************************************************************************************//**
 * \fn	void MFMTensorEstimateCostFunction::InitComputeNumberOfParameters()
 *
 * \brief	Initialises the value of  m_NbOfParameters depending on the tensor estimation
 * 			strategy and the fraction estimation strategy.  Must be called once before being able
 * 			to use GetNumberOfParameters() 
 *
 * \author	Benoit Scherrer
 * \date	April 2010
 *
 * \exception	itk::ExceptionObject	Thrown when error. 
 *
 * \sa	GetNumberOfParameters()
*************************************************************************************************/
void MFMTensorEstimateCostFunction::InitComputeNumberOfParameters() 
{
	if ( m_InternalVoxelData==NULL ) {
		m_NbOfParameters=-1;
		return;
	}

	int nb=0;

	nb = 6;
	if ( m_EstimationStrategy & ESTIMATION_EACH_CYLINDER ) nb -= 1;
	if ( m_EstimationStrategy & ESTIMATION_SAME_PLANE ) nb -= 2;

	if (m_InternalVoxelData->nbTensorAtVoxel>=2)
	{
		int nbPerTensor = 6;
		if ( m_EstimationStrategy & ESTIMATION_SAME_MAJOR_EIGEN_VALUE ) nbPerTensor-=1;
		if ( m_EstimationStrategy & ESTIMATION_SAME_MINOR_EIGEN_VALUE ) nbPerTensor-=1;
		if ( m_EstimationStrategy & ESTIMATION_EACH_CYLINDER ) nbPerTensor-=1;
		if ( m_EstimationStrategy & ESTIMATION_SAME_PLANE ) nbPerTensor-=2;

		nb += nbPerTensor*(m_InternalVoxelData->nbTensorAtVoxel-1);
	}

	m_NbOfParameters = ( nb+GetNumberOfFractionsToEstimate());

	if ( m_EstimationStrategy & ESTIMATION_RADIALBASELINE )
		m_NbOfParameters += m_InternalVoxelData->nbTensorAtVoxel;
}

/**********************************************************************************************//**
 * \fn	int MFMTensorEstimateCostFunction::GetNumberOfFractionsToEstimate() const
 *
 * \brief	Gets the number of fractions to estimate. 
 *
 * \author	Benoit Scherrer
 * \date	February 2010
 *
 * \exception	itk::ExceptionObject	Thrown when exception. 
 *
 * \return	The number of fractions to estimate. 
*************************************************************************************************/
int MFMTensorEstimateCostFunction::GetNumberOfFractionsToEstimate() const
{
	if ( m_NbFractionsToEstimate<0 ) 
		throw itk::ExceptionObject(__FILE__, __LINE__, "Error. The number of fractions to estimate has not been computed yet.","MFMTensorEstimateCostFunction::GetNumberOfFractionsToEstimate()");

	else return m_NbFractionsToEstimate;
}

/**********************************************************************************************//**
 * \fn	void MFMTensorEstimateCostFunction::InitComputeNumberOfFractionsToEstimate()
 *
 * \brief	Initialises the value of  m_NbFractionsToEstimate depending on the fraction
 * 			estimation strategy.  Must be called once before being able to use
 * 			GetNumberOfFractionsToEstimate() 
 *
 * \author	Benoit Scherrer
 * \date	April 2010
 *
 * \exception	itk::ExceptionObject	Thrown when error. 
 *
 * \sa	GetNumberOfFractionsToEstimate()
*************************************************************************************************/
void MFMTensorEstimateCostFunction::InitComputeNumberOfFractionsToEstimate()
{
	m_NbFractionsToEstimate=0;

	if ( m_InternalVoxelData==NULL ) {
		m_NbFractionsToEstimate=-1;
		return;
	}

	if ( (m_EstimationStrategy & (ESTIMATION_F_ISO+ESTIMATION_F_TENSORS)) == (ESTIMATION_F_ISO+ESTIMATION_F_TENSORS) )
		m_NbFractionsToEstimate = 1+m_InternalVoxelData->nbTensorAtVoxel;
	else if ( m_EstimationStrategy & ESTIMATION_F_TENSORS )
	{
		// If only one tensor and no iso, STRANGE ERROR ! (should not happen)
		// we just can't estimate any fraction! 
		if ( m_InternalVoxelData->nbTensorAtVoxel>1 ) 
			m_NbFractionsToEstimate = m_InternalVoxelData->nbTensorAtVoxel;
		else 
		{
			std::cout<<"WARNING. Strange error. Try to estimate tensor fractions, but only one tensor"<<std::endl;
			m_NbFractionsToEstimate = 0;
		}
	}
	else if ( m_EstimationStrategy & ESTIMATION_F_ISO )
		m_NbFractionsToEstimate = 2;	// 1 for the fractions, 1 for all other tensors

	// One more fraction for each tensor
	if ( m_EstimationStrategy&ESTIMATION_RADIALISO )
		m_NbFractionsToEstimate += m_InternalVoxelData->nbTensorAtVoxel;

	//!!bobyqa fractions
	//if ( m_NbFractionsToEstimate>1 && m_AlgoParameters->m_OptimizationAlgo==crl::OptimizerUtils::OPTIMIZER_BOBYQA )
	//	m_NbFractionsToEstimate-=1;
}

/**********************************************************************************************//**
 * \fn	void MFMTensorEstimateCostFunction::OneTensor2VnlMatrix(vnl_matrix<double> &out,
 * 		const TensorType &tensor ) const
 *
 * \brief	Converts a TensorType to a vnl_matrix<double>. 
 *
 * \author	Benoit Scherrer
 * \date	February 2010
 *
 * \param [in,out]	out	The out. 
 * \param	tensor		The tensor. 
*************************************************************************************************/
void MFMTensorEstimateCostFunction::OneTensor2VnlMatrix(VnlTensorType &out, const TensorType &tensor ) const
{
	//vnl_matrix: m[row][col]

	double *o = out.data_block();
	const TensorPixelType *p = tensor.GetDataPointer ();

	// Fast copy
	*(o) = *(p);
	*(o+1) = *(p+1);
	*(o+2) = *(p+2);
	*(o+3) = *(p+1);
	*(o+4) = *(p+3);
	*(o+5) = *(p+4);
	*(o+6) = *(p+2);
	*(o+7) = *(p+4);
	*(o+8) = *(p+5);
}

/**********************************************************************************************//**
 * \fn	void MFMTensorEstimateCostFunction::GetTensorFromSphericalParameterization( vnl_matrix<double> &out,
 * 		double alpha, double beta, double gamma, double l1, double l2, double l3 ) const
 *
 * \brief	Gets a 3x3 tensor from its spherical parameterization parameters. 
 *
 * \author	Benoit Scherrer
 * \date	May 2010
 *
 * \param [in,out]	out	The resulting tensor. 
 * \param	alpha		The alpha angle. 
 * \param	beta		The beta angle. 
 * \param	gamma		The gamma angle. 
 * \param	l1			The first (largest) eigen value. 
 * \param	l2			The second eigen value. 
 * \param	l3			The third (smallest) eigen value. 
*************************************************************************************************/
void MFMTensorEstimateCostFunction::GetTensorFromSphericalParameterization( VnlTensorType &out, VnlEigValueType &outEigenValues, VnlEigVectorType &outEigenVectors, double alpha, double beta, double gamma, double l1, double l2, double l3  ) const
{
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// first angle alpha: MUST BE around the minor eigen value (to have valid initialization & plane constraint)
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	VnlEigVectorType	R(3,3);

	//--------------------------------------
	// Compute the rotation matrix from the Euler angles
	// x-y-z convention
	//--------------------------------------
	double phi = beta;
	double theta=gamma;
	double psi=alpha;

	R(0,0) = cos(theta) * cos(phi);
	R(0,1) = cos(theta) * sin(phi);
	R(0,2) = -sin(theta);

	R(1,0) = sin(psi)*sin(theta)*cos(phi) - cos(psi)*sin(phi);
	R(1,1) = sin(psi)*sin(theta)*sin(phi) + cos(psi)*cos(phi);
	R(1,2) = cos(theta)*sin(psi);

	R(2,0) = cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi);
	R(2,1) = cos(psi)*sin(theta)*sin(phi) - sin(psi)*cos(phi);
	R(2,2) = cos(theta)*cos(psi);

	//--------------------------------------
	// Set the eigen values
	//--------------------------------------
	outEigenValues[0]=l3;
	outEigenValues[1]=l2;
	outEigenValues[2]=l1;

	//--------------------------------------
	// Compute the eigen vectors and the
	// reconstruced tensor
	//--------------------------------------
	out.set_size(3,3);
	vnl_fastops::AB(outEigenVectors, R, m_InternalVoxelData->OneTensorBasis); //, R);
	crl::FastOps::AtBA_Bsym_3x3(out, outEigenVectors, outEigenValues);
}

/**********************************************************************************************//**
 * \fn	void MFMTensorEstimateCostFunction::Parameters2Tensor( int idTensor,
 * 		vnl_matrix<double> &out, const ParametersType &parameters ) const
 *
 * \brief	Extract from the parameters a diffusion tensor  (vnl_matrix) 
 *
 * 			Organization of the parameters in memory: [Tensor1] [Tensor2] [Frac1] [Frac2]. The
 * 			number of values for each element depends on the estimation settings.
 *          
 * \author	Benoit Scherrer
 * \date	February 2010
 *
 * \exception	itk::ExceptionObject	Thrown when exception. 
 *
 * \param	idTensor	The tensor identifier . 
 * \param [in,out]	out	out	The output diffusion tensor. 
 * \param	parameters	The input optimizer parameters. 
*************************************************************************************************/
void MFMTensorEstimateCostFunction::Parameters2TensorAndFractions( VnlTensorSetType &tensors,  VnlFractionsType &fractions, const ParametersType &parameters ) const
{
	const double *p = parameters.data_block ();

	//-------------------------------------------
	// THE TENSORS
	//-------------------------------------------
	tensors.resize(m_InternalVoxelData->nbTensorAtVoxel);

	//-----------------
	// The first tensor
	//-----------------
	double a, b, c, T1_l1, T1_l2, T1_l3;
	
	// -- The angles --
	a=*(p++);
	if (m_EstimationStrategy & ESTIMATION_SAME_PLANE) { b=0; c=0; }
	else { b = *(p++); c=*(p++); }

	// -- The eigen values --
	T1_l1 = *(p++);
	T1_l2 = *(p++);
	if ( m_EstimationStrategy & ESTIMATION_EACH_CYLINDER ) T1_l3 = T1_l2;
	else T1_l3 = *(p++);

	GetTensorFromSphericalParameterization(tensors[0], m_InternalVoxelData->current_EigenValuesLE[0], m_InternalVoxelData->current_EigenVectorsLE[0], a, b, c, T1_l1, T1_l2, T1_l3 );

	//-----------------
	// The other tensors
	//-----------------
	if (m_InternalVoxelData->nbTensorAtVoxel>1)
	{
		for (int i=1; i<m_InternalVoxelData->nbTensorAtVoxel ; i++)
		{
			a=*(p++);
			if (m_EstimationStrategy & ESTIMATION_SAME_PLANE) { b=0; c=0; }
			else { b = *(p++); c=*(p++); }

			double l1, l2, l3;
			if (m_EstimationStrategy & ESTIMATION_SAME_MAJOR_EIGEN_VALUE) l1=T1_l1;
			else l1=*(p++);
			
			if ( m_EstimationStrategy & ESTIMATION_EACH_CYLINDER )
			{
				if (m_EstimationStrategy & ESTIMATION_SAME_MINOR_EIGEN_VALUE) { l2=T1_l3; l3=T1_l3; }
				else {
					l2 = *(p++);
					l3 = l2;
				}
			}
			else
			{
				l2=*(p++);
				if (m_EstimationStrategy & ESTIMATION_SAME_MINOR_EIGEN_VALUE) l3=T1_l3;
				else l3=*(p++);
			}

			GetTensorFromSphericalParameterization(tensors[i], m_InternalVoxelData->current_EigenValuesLE[i], m_InternalVoxelData->current_EigenVectorsLE[i],a, b, c, l1, l2, l3);
		}
	}

	//-------------------------------------------
	// The radial part
	//-------------------------------------------
	if ( m_EstimationStrategy & ESTIMATION_RADIALBASELINE )
	{
		VnlEigValueType		eValRadial(3); 

		for (int i=0; i<m_InternalVoxelData->nbTensorAtVoxel ; i++)
		{
			m_InternalVoxelData->current_RadialCst[i] = exp(*(p++));
			
			eValRadial[0] = m_InternalVoxelData->current_RadialCst[i];
			eValRadial[1] = m_InternalVoxelData->current_RadialCst[i];
			eValRadial[2] = 0;
			crl::FastOps::AtBA_Bdiag_3x3(m_InternalVoxelData->current_RadialMatrix[i], m_InternalVoxelData->current_EigenVectorsLE[i], eValRadial);
		}
	}
	//-------------------------------------------
	// Or the radial iso?
	//-------------------------------------------
	if ( m_EstimationStrategy & ESTIMATION_RADIALISO )
	{
		VnlEigValueType		eValRadial(3); 

		for (int i=0; i<m_InternalVoxelData->nbTensorAtVoxel ; i++)
		{	
			eValRadial[0] = g_MediumDiffCoeff;
			eValRadial[1] = g_MediumDiffCoeff;
			eValRadial[2] = 0;
			crl::FastOps::AtBA_Bdiag_3x3(m_InternalVoxelData->current_RadialMatrix[i], m_InternalVoxelData->current_EigenVectorsLE[i], eValRadial);
		}
	}

	//-------------------------------------------
	// NOW THE FRACTIONS
	//-------------------------------------------

	EstimationStrategyType fracEstimation = m_EstimationStrategy & (ESTIMATION_F_ISO+ESTIMATION_F_TENSORS);

	//-------------------------------------------
	// If only F_iso, just read 2 fractions from
	// the parameters: the isotropic part, and the other for all tensors
	//-------------------------------------------
	if ( fracEstimation == ESTIMATION_F_ISO )
	{
		// 2 fractions in parameters, but allocate the complete m_InternalVoxelData->nbTensorAtVoxel+1 (for easy access after)
		fractions.set_size(m_InternalVoxelData->nbTensorAtVoxel+1);

		//fill all the tensors fractions to the same value
		int k;
		for ( k=0 ; k<(m_InternalVoxelData->nbTensorAtVoxel) ; k++ )
			fractions[k] = *p;
		
		//!!bobyqa fractions
		//if ( m_AlgoParameters->m_OptimizationAlgo==crl::OptimizerUtils::OPTIMIZER_BOBYQA )
		//{
		//	fractions[m_NumberOfTensors]=1-fractions[0];
		//	for ( int k2=0 ; k2<(m_NumberOfTensors) ; k2++ )
		//		fractions[k2]=fractions[k2]/m_NumberOfTensors;
		//}
		//else
		{
			p++;

			// and the last fraction: for the isotropic part 
			fractions[k++] = *p;
		}
	}
	/*-----------------------------------------------
	 Estimate both isotropic part and tensor fractions, or only fractions
	-----------------------------------------------*/
	else if ( fracEstimation !=0 )
	{	
		// NB: ESTIMATION_RADIALISO is accepted only if ESTIMATION_F_TENSORS is also set.
		// ( else exception in MFMTensorEstimateFilter::BeforeThreadedGenerateData()  )

		int nb = m_InternalVoxelData->nbTensorAtVoxel;							// tensor fractions
		if ( fracEstimation&ESTIMATION_F_ISO ) nb+=1;							// + iso part?
		if ( m_EstimationStrategy&ESTIMATION_RADIALISO ) nb+=m_InternalVoxelData->nbTensorAtVoxel;	// + radial iso?
		fractions.set_size(nb);
		
		//------------------------------------------
		// BOBYQA uses bounded optimization so read n-1 fractions
		//------------------------------------------
		//!!bobyqa fractions
		/*if ( m_AlgoParameters->m_OptimizationAlgo==crl::OptimizerUtils::OPTIMIZER_BOBYQA )
		{
			double sum=0;
			for ( unsigned int k=0 ; k<fractions.size()-1 ; k++ )
			{
				fractions[k] = *p;
				sum+= *p;
				p++;
			}
			if ( sum<1 )
				fractions[fractions.size()-1]=1-sum;
			else
			{
			}
		}
		else*/
		{
			for ( unsigned int k=0 ; k<fractions.size() ; k++ )
			{
				fractions[k] = *p;
				p++;
			}		
		}
	}
	//-------------------------------------------
	// Else, if no fraction estimation, still allocate
	// the fractions vector with the same value
	// for each fraction
	//-------------------------------------------
	else
	{
		fractions.set_size(m_InternalVoxelData->nbTensorAtVoxel);
		fractions.fill(1.0);
	}

	/*-------------------------------------------
	 Apply the softmax transform to the fractions
	-------------------------------------------*/
	ApplySoftMaxToFractions();

}


/**********************************************************************************************//**
 * \fn	void MFMTensorEstimateCostFunction::InitAlgo(const TensorType &tensor)
 *
 * \brief	Initialises the algorithm. 
 *
 * \author	Benoit Scherrer
 * \date	February 2010
 *
 * \exception	itk::ExceptionObject	Thrown when exception. 
 *
 * \param	tensor	The tensor. 
*************************************************************************************************/
void MFMTensorEstimateCostFunction::InitAlgo(const TensorType &tensor, bool forceSingleTensor)
{

	/*-----------------------------------
	 First check the input tensor
	-----------------------------------*/
	// Check if negative or null eigen value for the input tensor
	VnlTensorType D(3,3), L(3,3);
	OneTensor2VnlMatrix(D,tensor); 

	// Compute the eigen values/vectors
	VnlEigVectorType		eVec(3,3);
	VnlEigValueType			eVal(3) ;
	VnlEigValueType			eValLE(3);
	EigenAnalysis.ComputeEigenValuesAndVectors(D, eVal, eVec);

	//--------------------------------------
	// If the initial one-tensor solution comes from MedINRIA, divides the eigen values by 1000
	//--------------------------------------
	if ( m_AlgoParameters->m_InitialOneTensorFromMedINRIA )
	{
		for ( int i=0; i<3 ; i++ ) eVal[i]=eVal[i]/1000;
	}

	//--------------------------------------
	// Look for the minimum positive eigen value
	//--------------------------------------
	double minPosEigenValue=1e10;
	for ( int i=0 ; i<3 ; i++ )	
		if ( eVal[i]<minPosEigenValue && eVal[i]>0 ) minPosEigenValue = eVal[i];
	if ( minPosEigenValue==1e10 ) minPosEigenValue=1e-4;

	//--------------------------------------
	// Then replace all negative or null eigen values  
	//--------------------------------------
	bool modifiedD=false;
	for ( int i=0 ; i<3 ; i++ )
	{
		if ( eVal[i]<=0 ) 
		{	
			if ( modifiedD ) eVal[i] = minPosEigenValue/3;
			else eVal[i] = minPosEigenValue/2;

			modifiedD=true;
		}
	}
	if ( m_AlgoParameters->m_InitialOneTensorFromMedINRIA ) modifiedD=true;

	//--------------------------------------
	// If modified, recompute !
	//--------------------------------------
	if ( modifiedD )
	{
		crl::FastOps::AtBA_Bsym_3x3(D, eVec, eVal); 
		// and redecompose
		EigenAnalysis.ComputeEigenValuesAndVectors(D, eVal, eVec);	
	}

	//--------------------------------------
	// Save the one-tensor solution basis
	// and eigen values
	//--------------------------------------
	m_InternalVoxelData->OneTensorBasis.set_size(3,3);
	m_InternalVoxelData->OneTensorBasis = eVec;
	m_InternalVoxelData->OneTensorEigenValues.set_size(3);
	m_InternalVoxelData->OneTensorEigenValues = eVal;
	
	/*-----------------------------------
	 Precompute the values for the isotropic part
	 (constant medium diff coeff)
	-----------------------------------*/
	m_InternalVoxelData->SiL_Isotropic.clear();
	if ( m_EstimationStrategy & ESTIMATION_F_ISO )
	{
		for ( unsigned int i=0; i<m_AlgoParameters->m_BValues.size() ; i++ )
			m_InternalVoxelData->SiL_Isotropic.push_back( exp(- m_AlgoParameters->m_BValues[i] * g_MediumDiffCoeff ));
	}

	m_InternalVoxelData->fastDiffusionFraction = m_AlgoParameters->m_FastDiffusionFraction;

	/*----------------------------------
	 Initializes the tensors
	----------------------------------*/
	m_InternalVoxelData->stabilized_L.clear();
	m_InternalVoxelData->lastParameters.set_size(this->GetNumberOfParameters());
	m_InternalVoxelData->lowerBound.set_size(this->GetNumberOfParameters());
	m_InternalVoxelData->upperBound.set_size(this->GetNumberOfParameters());

	double lmin, lmax;
	lmin = 1e-5; lmax = 2*eVal[2];
	if ( lmax<lmin ) lmin = lmax/5;
	lmin = safe_log(lmin);
	lmax = safe_log(lmax);

	/*-----------------------------------
	If one tensor,just push L and set init fractions
	-----------------------------------*/
	if ( m_AlgoParameters->m_NumberOfTensors==1 || forceSingleTensor)
	{
		//--------------------------------------
		// To convert in log-Euclidean
		//--------------------------------------
		for ( int i=0 ; i<3 ; i++ )	eValLE[i] = safe_log(eVal[i]);

		//--------------------------------------
		// Set the init parameters
		//--------------------------------------
		int id=0;
		int idBounds=0;

		//--------------------------------------
		// Set the tensor
		//--------------------------------------
		m_InternalVoxelData->nbTensorAtVoxel = 1;
		m_InternalVoxelData->lastParameters[id++]=0;
		m_InternalVoxelData->lowerBound[idBounds] = -M_PI; m_InternalVoxelData->upperBound[idBounds++] = M_PI;

		if ((m_EstimationStrategy & ESTIMATION_SAME_PLANE)==0) 
		{
			m_InternalVoxelData->lastParameters[id++]=0;
			m_InternalVoxelData->lastParameters[id++]=0;
			m_InternalVoxelData->lowerBound[idBounds] = -M_PI/2; m_InternalVoxelData->upperBound[idBounds++] = M_PI/2;
			m_InternalVoxelData->lowerBound[idBounds] = -M_PI/2; m_InternalVoxelData->upperBound[idBounds++] = M_PI/2;
		}

		m_InternalVoxelData->lastParameters[id++]=eValLE[2];
		m_InternalVoxelData->lastParameters[id++]=eValLE[1];
		m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;
		m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;

		if ((m_EstimationStrategy & ESTIMATION_EACH_CYLINDER)==0)
		{
			m_InternalVoxelData->lastParameters[id++]=eValLE[0];
			m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;
		}

		//--------------------------------------------
		// Set the radial baseline
		//--------------------------------------------
		if ( m_AlgoParameters->m_EstimationStrategy & ESTIMATION_RADIALBASELINE ) 
		{
			m_InternalVoxelData->lastParameters[id++]=m_InternalVoxelData->OneTensorEigenValues[0]; // the baseline init value for the matrix //0.05;	// the radial baseline fraction
			m_InternalVoxelData->lowerBound[idBounds] = -100; m_InternalVoxelData->upperBound[idBounds++] = 5;
		}

		//--------------------------------------
		// For one tensor: fraction only if isotropic part
		//--------------------------------------
		if ( m_AlgoParameters->m_EstimationStrategy & ESTIMATION_F_ISO ) 
		{
			//!!bobyqa fractions
			//if ( m_AlgoParameters->m_OptimizationAlgo==crl::OptimizerUtils::OPTIMIZER_BOBYQA )
			//{
			//	m_InternalVoxelData->lastParameters[id++]=0.8;	// the tensor fractions
			//	m_InternalVoxelData->lowerBound[idBounds] = 0; m_InternalVoxelData->upperBound[idBounds++] = 1;
			//}
			//else
			{
				m_InternalVoxelData->lastParameters[id++]=0.8;	// the tensor fraction
				m_InternalVoxelData->lastParameters[id++]=0.2;	// the isotropic part
				m_InternalVoxelData->lowerBound[idBounds] = 0; m_InternalVoxelData->upperBound[idBounds++] = 1;
				m_InternalVoxelData->lowerBound[idBounds] = 0; m_InternalVoxelData->upperBound[idBounds++] = 1;
			}

			if ( m_AlgoParameters->m_EstimationStrategy & ESTIMATION_RADIALISO )
			{
				m_InternalVoxelData->lastParameters[id-1] -= 0.01;
				m_InternalVoxelData->lastParameters[id++]= 0.01;	// the radial iso fractions
				m_InternalVoxelData->lowerBound[idBounds] = 0; m_InternalVoxelData->upperBound[idBounds++] = 1;
			}
		}

		//--------------------------------------
		// a small checking
		//--------------------------------------
		if (id!=(int)this->GetNumberOfParameters())
			throw itk::ExceptionObject(__FILE__, __LINE__, "Internal Error. Not all parameters have been initialized","");
		if (idBounds!=(int)this->GetNumberOfParameters())
			throw itk::ExceptionObject(__FILE__, __LINE__, "Internal Error. Not all bound parameters have been initialized","");

	}

	/*-----------------------------------
	If two tensors, compute two	perpendicular tensors in the plane
	of the largest eigen values
	-----------------------------------*/
	else if ( m_AlgoParameters->m_NumberOfTensors>=2 )
	{	
		m_InternalVoxelData->nbTensorAtVoxel = m_AlgoParameters->m_NumberOfTensors;

		//---------------------------------------------
		// Define the angle between the two tensors
		//---------------------------------------------
		double angle;
		if ( m_AlgoParameters->m_RandomInit ) angle = (double)(rand()%4000)/4000.0;
		else angle = eVal[1]/eVal[2] * 3.14156/4.0;

		//---------------------------------------------
		// Setup initial parameters for the optimizer
		//---------------------------------------------

		int id=0;
		int idBounds=0;

		for ( int tensor=0; tensor<m_AlgoParameters->m_NumberOfTensors; tensor++ )
		{
			//-------------------------------------------
			// First two tensors: 
			// First angle: + and - angle /other: random
			//-------------------------------------------
			if ( tensor<2 )
				m_InternalVoxelData->lastParameters[id++]=(1-tensor)*tensor - tensor*angle;
			else
				m_InternalVoxelData->lastParameters[id++]=(double)(rand()%4000)/4000.0;
			m_InternalVoxelData->lowerBound[idBounds] = -M_PI; m_InternalVoxelData->upperBound[idBounds++] = M_PI;

			//-------------------------------------------
			// Other euler angles
			//-------------------------------------------
			if ((m_EstimationStrategy & ESTIMATION_SAME_PLANE)==0) 
			{
				m_InternalVoxelData->lastParameters[id++]=0;
				m_InternalVoxelData->lastParameters[id++]=0;
				m_InternalVoxelData->lowerBound[idBounds] = -M_PI/2; m_InternalVoxelData->upperBound[idBounds++] = M_PI/2;
				m_InternalVoxelData->lowerBound[idBounds] = -M_PI/2; m_InternalVoxelData->upperBound[idBounds++] = M_PI/2;
			}

			if ( tensor==0 )
			{
				m_InternalVoxelData->lastParameters[id++]=safe_log(eVal[2]);
				m_InternalVoxelData->lastParameters[id++]=safe_log(eVal[0]+1e-5);
				m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;
				m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;
				if ((m_EstimationStrategy & ESTIMATION_EACH_CYLINDER)==0)
				{
					m_InternalVoxelData->lastParameters[id++]=safe_log(eVal[0]);
					m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;
				}
			}
			else
			{
				if ((m_EstimationStrategy & ESTIMATION_SAME_MAJOR_EIGEN_VALUE)==0) 
				{
					m_InternalVoxelData->lastParameters[id++]=safe_log(eVal[2]);
					m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;
				}

				if ((m_EstimationStrategy & ESTIMATION_EACH_CYLINDER))
				{
					if ((m_EstimationStrategy & ESTIMATION_SAME_MINOR_EIGEN_VALUE)==0) 
					{
						m_InternalVoxelData->lastParameters[id++]=safe_log(eVal[0]+1e-5);
						m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;
					}
				}
				else
				{
					m_InternalVoxelData->lastParameters[id++]=safe_log(eVal[0]+1e-5);
					m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;
					if ((m_EstimationStrategy & ESTIMATION_SAME_MINOR_EIGEN_VALUE)==0) 
					{
						m_InternalVoxelData->lastParameters[id++]=safe_log(eVal[0]+1e-5);	
						m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;
					}
				}
			}
		}

		////--------------------------------------
		//// First tensor
		////--------------------------------------
		//m_InternalVoxelData->lastParameters[id++]=angle;
		//m_InternalVoxelData->lowerBound[idBounds] = -M_PI; m_InternalVoxelData->upperBound[idBounds++] = M_PI;

		//if ((m_EstimationStrategy & ESTIMATION_SAME_PLANE)==0) 
		//{
		//	m_InternalVoxelData->lastParameters[id++]=0;
		//	m_InternalVoxelData->lastParameters[id++]=0;
		//	m_InternalVoxelData->lowerBound[idBounds] = -M_PI/2; m_InternalVoxelData->upperBound[idBounds++] = M_PI/2;
		//	m_InternalVoxelData->lowerBound[idBounds] = -M_PI/2; m_InternalVoxelData->upperBound[idBounds++] = M_PI/2;
		//}

		//m_InternalVoxelData->lastParameters[id++]=safe_log(eVal[2]);
		//m_InternalVoxelData->lastParameters[id++]=safe_log(eVal[0]+1e-5);
		//m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;
		//m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;
		//if ((m_EstimationStrategy & ESTIMATION_EACH_CYLINDER)==0)
		//{
		//	m_InternalVoxelData->lastParameters[id++]=safe_log(eVal[0]);
		//	m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;
		//}

		////--------------------------------------
		//// Second tensor
		////--------------------------------------
		//m_InternalVoxelData->lastParameters[id++]=-angle;
		//m_InternalVoxelData->lowerBound[idBounds] = -M_PI; m_InternalVoxelData->upperBound[idBounds++] = M_PI;
		//if ((m_EstimationStrategy & ESTIMATION_SAME_PLANE)==0) 
		//{
		//	m_InternalVoxelData->lastParameters[id++]=0;
		//	m_InternalVoxelData->lastParameters[id++]=0;
		//	m_InternalVoxelData->lowerBound[idBounds] = -M_PI/2; m_InternalVoxelData->upperBound[idBounds++] = M_PI/2;
		//	m_InternalVoxelData->lowerBound[idBounds] = -M_PI/2; m_InternalVoxelData->upperBound[idBounds++] = M_PI/2;
		//}

		//if ((m_EstimationStrategy & ESTIMATION_SAME_MAJOR_EIGEN_VALUE)==0) 
		//{
		//	m_InternalVoxelData->lastParameters[id++]=safe_log(eVal[2]);
		//	m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;
		//}

		//if ((m_EstimationStrategy & ESTIMATION_EACH_CYLINDER))
		//{
		//	if ((m_EstimationStrategy & ESTIMATION_SAME_MINOR_EIGEN_VALUE)==0) 
		//	{
		//		m_InternalVoxelData->lastParameters[id++]=safe_log(eVal[0]+1e-5);
		//		m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;
		//	}
		//}
		//else
		//{
		//	m_InternalVoxelData->lastParameters[id++]=safe_log(eVal[0]+1e-5);
		//	m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;
		//	if ((m_EstimationStrategy & ESTIMATION_SAME_MINOR_EIGEN_VALUE)==0) 
		//	{
		//		m_InternalVoxelData->lastParameters[id++]=safe_log(eVal[0]+1e-5);	
		//		m_InternalVoxelData->lowerBound[idBounds] = lmin; m_InternalVoxelData->upperBound[idBounds++] = lmax;
		//	}
		//}

		//--------------------------------------
		// The baseline parameters
		//--------------------------------------
		if ( m_AlgoParameters->m_EstimationStrategy & ESTIMATION_RADIALBASELINE ) 
		{
			for ( int tensor=0; tensor<m_AlgoParameters->m_NumberOfTensors; tensor++ )
				m_InternalVoxelData->lastParameters[id++]=m_InternalVoxelData->OneTensorEigenValues[0]; // the baseline init value for the matrix // the radial baseline fraction
		
			for ( int tensor=0; tensor<m_AlgoParameters->m_NumberOfTensors; tensor++ )
				m_InternalVoxelData->lowerBound[idBounds] = -100; m_InternalVoxelData->upperBound[idBounds++] = 5;
		}

		////--------------------------------------
		//// The baseline parameters
		////--------------------------------------
		//if ( m_AlgoParameters->m_EstimationStrategy & ESTIMATION_RADIALBASELINE ) 
		//{
		//	m_InternalVoxelData->lastParameters[id++]=m_InternalVoxelData->OneTensorEigenValues[0]; // the baseline init value for the matrix // the radial baseline fraction
		//	m_InternalVoxelData->lastParameters[id++]=m_InternalVoxelData->OneTensorEigenValues[0]; // the baseline init value for the matrix // the radial baseline fraction

		//	m_InternalVoxelData->lowerBound[idBounds] = -100; m_InternalVoxelData->upperBound[idBounds++] = 5;
		//	m_InternalVoxelData->lowerBound[idBounds] = -100; m_InternalVoxelData->upperBound[idBounds++] = 5;
		//}

		//--------------------------------------
		// The fractions
		//--------------------------------------
		EstimationStrategyType fracEstimation = m_EstimationStrategy & (ESTIMATION_F_ISO+ESTIMATION_F_TENSORS);

		if ( fracEstimation == ESTIMATION_F_ISO )
		{
			m_InternalVoxelData->lastParameters[id++]=0.8;	//all tensors
			m_InternalVoxelData->lowerBound[idBounds] = 0; m_InternalVoxelData->upperBound[idBounds++] = 1;

			//!!bobyqa fractions
			//if ( m_AlgoParameters->m_OptimizationAlgo!=crl::OptimizerUtils::OPTIMIZER_BOBYQA )
			{
				m_InternalVoxelData->lastParameters[id++]=0.2;	//iso
				m_InternalVoxelData->lowerBound[idBounds] = 0; m_InternalVoxelData->upperBound[idBounds++] = 1;
			}
		}
		else if (fracEstimation==ESTIMATION_F_TENSORS)
		{
			for ( int tensor=0; tensor<m_AlgoParameters->m_NumberOfTensors; tensor++ )
			{
				m_InternalVoxelData->lastParameters[id++]=1.0/((double)m_AlgoParameters->m_NumberOfTensors);	//tensor 1
				m_InternalVoxelData->lowerBound[idBounds] = 0; m_InternalVoxelData->upperBound[idBounds++] = 1;
			}

			//m_InternalVoxelData->lastParameters[id++]=0.5;	//tensor 1
			//m_InternalVoxelData->lowerBound[idBounds] = 0; m_InternalVoxelData->upperBound[idBounds++] = 1;

			////!!bobyqa fractions
			////if ( m_AlgoParameters->m_OptimizationAlgo!=crl::OptimizerUtils::OPTIMIZER_BOBYQA )
			//{
			//	m_InternalVoxelData->lastParameters[id++]=0.5;	//tensor 2
			//	m_InternalVoxelData->lowerBound[idBounds] = 0; m_InternalVoxelData->upperBound[idBounds++] = 1;
			//}
		}
		else if (fracEstimation==ESTIMATION_F_TENSORS+ESTIMATION_F_ISO)
		{
			for ( int tensor=0; tensor<m_AlgoParameters->m_NumberOfTensors; tensor++ )
			{
				m_InternalVoxelData->lastParameters[id++]=0.9/((double)m_AlgoParameters->m_NumberOfTensors);	//tensor 1
				m_InternalVoxelData->lowerBound[idBounds] = 0; m_InternalVoxelData->upperBound[idBounds++] = 1;
			}

			//m_InternalVoxelData->lastParameters[id++]=0.45;	//tensor1
			//m_InternalVoxelData->lastParameters[id++]=0.45;	//tensor2
			//m_InternalVoxelData->lowerBound[idBounds] = 0; m_InternalVoxelData->upperBound[idBounds++] = 1;
			//m_InternalVoxelData->lowerBound[idBounds] = 0; m_InternalVoxelData->upperBound[idBounds++] = 1;

			//!!bobyqa fractions
			//if ( m_AlgoParameters->m_OptimizationAlgo!=crl::OptimizerUtils::OPTIMIZER_BOBYQA )
			{
				m_InternalVoxelData->lastParameters[id++]=0.1;	//iso
				m_InternalVoxelData->lowerBound[idBounds] = 0; m_InternalVoxelData->upperBound[idBounds++] = 1;
			}
		}

		if ( m_AlgoParameters->m_EstimationStrategy & ESTIMATION_RADIALISO )
		{
			if ( m_AlgoParameters->m_NumberOfTensors!=2 )
				throw itk::ExceptionObject(__FILE__, __LINE__, "Error. Initialization of radial part and more than 2 tensors to do","");

			m_InternalVoxelData->lastParameters[id-1] -= 0.01;
			m_InternalVoxelData->lastParameters[id-2] -= 0.01;
			m_InternalVoxelData->lastParameters[id++]= 0.01;	// the radial iso fractions
			m_InternalVoxelData->lastParameters[id++]= 0.01;	// the radial iso fractions
			m_InternalVoxelData->lowerBound[idBounds] = 0; m_InternalVoxelData->upperBound[idBounds++] = 1;
			m_InternalVoxelData->lowerBound[idBounds] = 0; m_InternalVoxelData->upperBound[idBounds++] = 1;
		}

		//--------------------------------------
		// a small checking
		//--------------------------------------
		if (id!=(int)this->GetNumberOfParameters())
			throw itk::ExceptionObject(__FILE__, __LINE__, "Error. Not all bound parameters have been initialized","");
		if (idBounds!=(int)this->GetNumberOfParameters())
			throw itk::ExceptionObject(__FILE__, __LINE__, "Error. Not all bound parameters have been initialized","");

	}
	else
		throw itk::ExceptionObject(__FILE__, __LINE__, "Initialization currently works for only one or two tensors","");


	//---------------------------------------------
	// Initialize some shared structures
	//---------------------------------------------
	int nbT = m_InternalVoxelData->nbTensorAtVoxel;
	m_InternalVoxelData->stabilized_L.resize(nbT);
	m_InternalVoxelData->current_L.resize(nbT);
	m_InternalVoxelData->current_EigenValuesLE.resize(nbT);
	m_InternalVoxelData->current_EigenVectorsLE.resize(nbT);
	m_InternalVoxelData->current_RadialMatrix.resize(nbT);
	m_InternalVoxelData->current_RadialCst.resize(nbT);
	for (int tensor=0; tensor<nbT ; tensor++ )
	{
		m_InternalVoxelData->current_EigenValuesLE[tensor].set_size(3);
		m_InternalVoxelData->current_EigenVectorsLE[tensor].set_size(3,3);
		m_InternalVoxelData->current_L[tensor].set_size(3,3);
		m_InternalVoxelData->stabilized_L[tensor].set_size(3,3);
		m_InternalVoxelData->current_RadialMatrix[tensor].set_size(3,3);
	}
	m_InternalVoxelData->current_GradientRadialProjection.set_size(nbT, m_AlgoParameters->m_g.size());
	m_InternalVoxelData->current_GradientTensorProjection.set_size(nbT, m_AlgoParameters->m_g.size());
	m_InternalVoxelData->SiLs.set_size(nbT, m_AlgoParameters->m_g.size());

	//---------------------------------------------
	// Convert that initial parameters to the initial values of L and Fractions
	// (because neighboring Ls must be computed)
	//---------------------------------------------
	Parameters2TensorAndFractions( m_InternalVoxelData->stabilized_L, m_InternalVoxelData->Fractions, m_InternalVoxelData->lastParameters );
	for ( unsigned int i=0; i<m_InternalVoxelData->stabilized_L.size(); i++ )
		m_InternalVoxelData->current_L[i] = m_InternalVoxelData->stabilized_L[i];

	//-------------------------------------------------------
	// Initializes the Sils for the values of the Ls
	//-------------------------------------------------------

	// Compute the SiLs from m_InternalVoxelData->current_L, m_InternalVoxelData->current_EigenValuesLE, m_InternalVoxelData->current_EigenVectorsLE
	ComputeSiLs();
}

/**********************************************************************************************//**
 * \fn	void MFMTensorEstimateCostFunction::ApplySoftMaxToFractions() const
 *
 * \brief	Applies the softmax transform to the fractions, to ensure they belong to [0 1] and
 * 			all sum to 1.
 * 			
 * 			The softmax transform is defined by:  
 * 			\f$ \eta_i = \frac{\exp f_i }{\sum_j \exp f_j} \f$. 
 *
 * \author	Benoit Scherrer
 * \date	February 2010
*************************************************************************************************/
void MFMTensorEstimateCostFunction::ApplySoftMaxToFractions() const
{
	double sum = 0;
	m_InternalVoxelData->normalizedFractions.set_size(m_InternalVoxelData->Fractions.size());
	
	unsigned int nbTensorFractions = m_InternalVoxelData->Fractions.size();
	
	//---------------------------------------------
	// IF BOBYQA don't use softmax transform, just normalize
	//---------------------------------------------
	//!!bobyqa fractions
	if ( m_AlgoParameters->m_OptimizationAlgo==crl::OptimizerUtils::OPTIMIZER_BOBYQA )
	{
		for ( unsigned int i=0; i<nbTensorFractions ; i++ )
			sum += m_InternalVoxelData->Fractions[i];

		for ( unsigned int i=0; i<nbTensorFractions ; i++ )
			m_InternalVoxelData->normalizedFractions[i]=m_InternalVoxelData->Fractions[i]/sum;
	}
	//---------------------------------------------
	// Else apply softmax
	//---------------------------------------------
	else
	{
		//---------------------------------------------
		// Compute the sum of the exp of the fractions
		//---------------------------------------------

		for ( unsigned int i=0; i<nbTensorFractions ; i++ )
		{
			m_InternalVoxelData->normalizedFractions[i]=exp(m_InternalVoxelData->Fractions[i]);
			sum += m_InternalVoxelData->normalizedFractions[i];
		}

		//---------------------------------------------
		// Compute the soft-max fractions
		//---------------------------------------------
		for ( unsigned int i=0; i<nbTensorFractions ; i++ )
			m_InternalVoxelData->normalizedFractions[i]=m_InternalVoxelData->normalizedFractions[i]/sum;
	}
}

/**********************************************************************************************//**
 * \fn	MFMTensorEstimateCostFunction::ParametersType MFMTensorEstimateCostFunction::GetOptimizerParametersFromPrecomputeData()
 *
 * \brief	Returns the ParametersType structure (optimizer parameters) from the current value of
 * 			m_InternalVoxelData. 
 *
 * \author	Benoit Scherrer
 * \date	February 2010
 *
 * \return	The optimizer position from precompute data. 
*************************************************************************************************/
MFMTensorEstimateCostFunction::ParametersType MFMTensorEstimateCostFunction::GetOptimizerParametersFromPrecomputeData()
{
	ParametersType p(this->GetNumberOfParameters());
	m_InternalVoxelData->lastParameters.set_size(p.size());

	for (unsigned int i=0; i<p.size() ; i++ )
		p[i]=m_InternalVoxelData->lastParameters[i];

	return p;
}

/**********************************************************************************************//**
 * \fn	void MFMTensorEstimateCostFunction::EndAlgoIteration(ParametersType &parameters)
 *
 * \brief	Finishes an iteration of the algorithm from the parameters: update the
 * 			m_InternalVoxelData structure:
 * 			- Recomputes the final SiLs for each tensors (from the parameters)
 * 			- Update the log-euclidean eigen vectors/eigen values to make them available to
 * 			  neighbors (for the regularization)
 * 			- Recompute the fractions from the parameters. 
 *
 * \author	Benoit Scherrer
 * \date	February 2010
 *
 * \param [in,out]	parameters	The input optimizer parameters. 
*************************************************************************************************/
void MFMTensorEstimateCostFunction::EndAlgoIteration(ParametersType &parameters)
{	
	//---------------------------------------------
	// Use to see if we stop the estimation for that tensor (convergence?)
	//---------------------------------------------
	double									diffTensorMetric = 0;

	//---------------------------------------------
	// Copy the list of parameters  (the size of lastParameters is already OK)
	// to save them
	//---------------------------------------------
	for (unsigned int i=0; i<parameters.size() ; i++ )
		m_InternalVoxelData->lastParameters[i] = parameters[i];

	//---------------------------------------------
	// Transform the parameters to tensors/fractions
	//---------------------------------------------
	Parameters2TensorAndFractions( m_InternalVoxelData->current_L, m_InternalVoxelData->Fractions, parameters ) ;

	//---------------------------------------------
	// Compute the SiLs from m_InternalVoxelData->current_L, m_InternalVoxelData->current_EigenValuesLE, m_InternalVoxelData->current_EigenVectorsLE
	// Needed in case of robust estimation (for ComputeC_ForRobustEstimation)
	//---------------------------------------------
	ComputeSiLs();

	//---------------------------------------------
	// Measure the difference between the new tensor and the old one
	//---------------------------------------------
	for (int i=0 ; i<m_InternalVoxelData->nbTensorAtVoxel ; i++ )
		diffTensorMetric += (m_InternalVoxelData->current_L[i]-m_InternalVoxelData->stabilized_L[i]).frobenius_norm();

	//---------------------------------------------
	// Save in m_InternalVoxelData the new value of L (result of the optimization)
	// We only write here the value of L, so protect with the mutex
	//---------------------------------------------
	if ( m_UseRegularization ) m_InternalVoxelData->mutexAccessL.Lock();
	for (int i=0 ; i<m_InternalVoxelData->nbTensorAtVoxel ; i++ ) m_InternalVoxelData->stabilized_L[i]  = m_InternalVoxelData->current_L[i] ;
	if ( m_UseRegularization ) m_InternalVoxelData->mutexAccessL.Unlock();

	//---------------------------------------------
	// Determine if the convergence criterion has been reached
	//---------------------------------------------
	if ( diffTensorMetric < 0.01 ) m_InternalVoxelData->HasConverged = true; 
}

/**********************************************************************************************//**
 * \fn	void MFMTensorEstimateCostFunction::ComputeFinalTensors(std::vector<TensorType > &tensors ) const
 *
 * \brief	Calculates the final tensors estimate from the m_InternalVoxelData->L,  
 * 			m_InternalVoxelData->normalizedFractions and m_InternalVoxelData->OneTensorBasis. 
 *			Output the results in tensors.
 *    RETURNS m_AlgoParameters->m_NumberOfTensors tensors (will null tensors if only one was estimated)
 *
 * \author	Benoit Scherrer
 * \date	February 2010
 *
 * \param [in,out]	tensors	The tensors. 
*************************************************************************************************/
void MFMTensorEstimateCostFunction::ComputeFinalTensors(std::vector<TensorType > &tensors ) const
{
	tensors.clear();
	VnlFractionsType	fractions;
	VnlTensorType		D(3,3);
	VnlEigValueType		eVal(3) ;

	//---------------------------------------------
	// Be sure to have the last version of the estimated tensor
	//---------------------------------------------
	Parameters2TensorAndFractions( m_InternalVoxelData->stabilized_L, fractions, m_InternalVoxelData->lastParameters ) ;
	
	//---------------------------------------------
	// FOR EACH OUTPUT tensor
	//---------------------------------------------
	for ( int k=0 ; k<m_AlgoParameters->m_NumberOfTensors ; k++ )
	{
		//---------------------------------------------
		// Prune the tensor?
		//---------------------------------------------
		if ( k>=m_InternalVoxelData->nbTensorAtVoxel
			|| m_InternalVoxelData->normalizedFractions[k]<=m_AlgoParameters->m_PruneThreshold )
			D.fill(0);

		//---------------------------------------------
		// If not...
		//---------------------------------------------
		else
		{
			// Compute the eigen values/vectors
			//vnl_matrix<double>		eVec(3,3);
			// Use the 3D version of L
			//EigenAnalysis.ComputeEigenValuesAndVectors(m_InternalVoxelData->stabilized_L[k], eVal, eVec);
			// Put it in eucl
			//for ( int i=0 ; i<3 ; i++ )	eVal[i] = exp(eVal[i]);
			//crl::FastOps::AtBA_Bdiag_3x3(D, eVec, eVal);

			//---------------------------------------------
			// Converts log-Eucl -> Eucl
			//---------------------------------------------
			for ( int i=0 ; i<3 ; i++ )	eVal[i] = exp(m_InternalVoxelData->current_EigenValuesLE[k][i]);
			crl::FastOps::AtBA_Bdiag_3x3(D, m_InternalVoxelData->current_EigenVectorsLE[k], eVal);
		}

		//---------------------------------------------
		//  and add to the list of tensors
		//---------------------------------------------
		TensorType tensor;
		tensor[0] = D(0,0);
		tensor[1] = D(0,1);
		tensor[2] = D(0,2);
		tensor[3] = D(1,1);
		tensor[4] = D(1,2);
		tensor[5] = D(2,2);
		tensors.push_back(tensor);
	}
}
		
void MFMTensorEstimateCostFunction::ComputeFinalFractions(VnlFractionsType &fractions ) const
{
	int nbFractions = m_AlgoParameters->m_NumberOfTensors;
	if ( m_UseIsotropicPart ) nbFractions += 1;

	fractions.set_size(nbFractions);

	for ( int i=0; i< m_AlgoParameters->m_NumberOfTensors; i++ )
	{
		if ( i<m_InternalVoxelData->nbTensorAtVoxel ) fractions[i] = m_InternalVoxelData->normalizedFractions[i];
		else fractions[i] = 0;
	}
	if ( m_UseIsotropicPart ) 
		fractions[m_AlgoParameters->m_NumberOfTensors] = m_InternalVoxelData->normalizedFractions[m_InternalVoxelData->normalizedFractions.size()-1];
}

/**********************************************************************************************//**
 * \fn	void MFMTensorEstimateCostFunction::ComputeSiLs() const
 *
 * \brief	Compute in m_InternalVoxelData->SiLs the SiL (without fraction) for ONE tensor L for
 * 			EACH direction (and also return eigen values/vectors)
 * 			            Sil() = exp (-bk * gk^t * Di * gk)
 * 			
 * 			WARNING : Compute the SiLs values from m_InternalVoxelData->current_L,
 * 			m_InternalVoxelData-
 * 			>current_EigenValuesLE, m_InternalVoxelData->current_EigenVectorsLE.  So these values
 * 			must be valid. 
 *
 * \author	Benoit Scherrer
 * \date	February 2010
*************************************************************************************************/
void MFMTensorEstimateCostFunction::ComputeSiLs() const
{
	_ASSERT( m_CurrentMaxBValueIndexToUse!=0 );
	
	VnlEigValueType		eVal(3); 
	VnlTensorType		expL(3,3);

	//---------------------------------------------
	// For each tensor
	//---------------------------------------------
	for ( int tensor=0 ; tensor<m_InternalVoxelData->nbTensorAtVoxel ; tensor++ )
	{

		// Compute D = exp(L) for this tensor 
		//EigenAnalysis.ComputeEigenValuesAndVectors(m_InternalVoxelData->current_L[tensor], 
		//	m_InternalVoxelData->current_EigenValuesLE[tensor], 
		//	m_InternalVoxelData->current_EigenVectorsLE[tensor] ); 

		//---------------------------------------------
		// Convert log-Eucl -> Eucl: D = exp(L)
		//---------------------------------------------
		for ( int i=0 ; i<3 ; i++ ) eVal[i] = exp(m_InternalVoxelData->current_EigenValuesLE[tensor][i]);
		crl::FastOps::AtBA_Bdiag_3x3(expL, m_InternalVoxelData->current_EigenVectorsLE[tensor], eVal);		// really more fast than	expL = eVecLE.transpose() * eVal * eVecLE;

		//---------------------------------------------
		// Update SiLs for each direction: new model: put log(Sil)
		//---------------------------------------------
#if NEW_FASTSLOWMODEL
		for( int i = 0; i< m_CurrentMaxBValueIndexToUse ; i++ )
			m_InternalVoxelData->SiLs(tensor,i) =  -m_AlgoParameters->m_BValues[i] * crl::FastOps::btAb_Asym_3x3(expL, m_AlgoParameters->m_g[i] ) ;
#else
		for( int i = 0; i< m_CurrentMaxBValueIndexToUse ; i++ )
		{
			m_InternalVoxelData->current_GradientTensorProjection = crl::FastOps::btAb_Asym_3x3(expL, m_AlgoParameters->m_g[i] );
			m_InternalVoxelData->SiLs(tensor,i) = exp( -m_AlgoParameters->m_BValues[i] * m_InternalVoxelData->current_GradientTensorProjection(tensor, i)  )  ;
		}
#endif

		//---------------------------------------------
		// Update the radial parts SiLs ?
		//---------------------------------------------
		if ( m_EstimationStrategy&ESTIMATION_RADIALBASELINE )
		{	
			for( int i = 0; i< m_CurrentMaxBValueIndexToUse ; i++ )
				m_InternalVoxelData->current_GradientRadialProjection(tensor,i) = crl::FastOps::btAb_Asym_3x3(m_InternalVoxelData->current_RadialMatrix[tensor], m_AlgoParameters->m_g[i] );
		}
	}
}

/**********************************************************************************************//**
 * \fn	void MFMTensorEstimateCostFunction::ComputeResiduals(std::vector& residuals)
 *
 * \brief	Calculates the residuals. 
 * \warning ComputeSiLs must have been called before, for all tensors
 *
 * \author	Benoit Scherrer
 * \date	July 2010
 *
 * \exception	itk::ExceptionObject	Thrown when exception. 
 *
 * \param [in,out]	residuals	The residuals. 
*************************************************************************************************/
void MFMTensorEstimateCostFunction::ComputeResiduals(std::vector<double> & residuals)
{
	residuals.resize(m_CurrentMaxBValueIndexToUse);

	//---------------------------------------------------------
	// For EACH gradient direction
	//---------------------------------------------------------
	_ASSERT( m_CurrentMaxBValueIndexToUse!=0 );

	for( int i = 0; i< m_CurrentMaxBValueIndexToUse ; i++ )
	{
		double estSi = 0;

#if NEW_FASTSLOWMODEL
		double fast=1;
		double slow=0;
		for ( int k=0 ; k<m_InternalVoxelData->nbTensorAtVoxel ; k++ )
		{
			fast *= exp( m_InternalVoxelData->normalizedFractions[k] * m_InternalVoxelData->SiLs(k,i)  )  ;
			slow += m_InternalVoxelData->normalizedFractions[k] *exp( m_InternalVoxelData->SiLs(k,i)  )  ;
		
			throw itk::ExceptionObject(__FILE__,__LINE__,"ERROR. Code source fast/slow and radial base line not updated", ""); 
			//!!if ( m_EstimationStrategy&ESTIMATION_RADIALBASELINE )
			//	slow += m_InternalVoxelData->normalizedFractions[k+m_InternalVoxelData->nbTensorAtVoxel] *exp( m_InternalVoxelData->SiLs_Radial(k,i)  )  ;
		}
		if ( m_UseIsotropicPart )
		{
			double fiso = m_InternalVoxelData->normalizedFractions[m_InternalVoxelData->normalizedFractions.size()-1];

			fast *= exp(-fiso*m_AlgoParameters->m_BValues[i] * g_MediumDiffCoeff);
			slow += fiso*exp(-m_AlgoParameters->m_BValues[i] * g_MediumDiffCoeff);
		}

		estSi = (1-m_InternalVoxelData->fastDiffusionFraction)*slow   +   m_InternalVoxelData->fastDiffusionFraction*fast;
#else
		
		for ( int k=0 ; k<m_InternalVoxelData->nbTensorAtVoxel ; k++ )	
		{
			if ( m_EstimationStrategy&ESTIMATION_RADIALBASELINE )
			{
				double projC = m_InternalVoxelData->current_GradientRadialProjection(k, i);
				double projD = m_InternalVoxelData->current_GradientTensorProjection(k, i);

				//double v = exp(-projC);
				//double v2=m_InternalVoxelData->SiLs(k,i) ;
				//double w_parallel = 0.1 + 0.9*exp(projD)/(exp(projD)+exp(projC));

				//estSi = m_InternalVoxelData->normalizedFractions[k] *
				//	(  w_parallel * v2    +   (1-w_parallel) * v );

				estSi = m_InternalVoxelData->SiLs(k,i) + projC * exp(-m_AlgoParameters->m_BValues[0]*projC );
			}
			else
			{
				estSi += m_InternalVoxelData->normalizedFractions[k] * m_InternalVoxelData->SiLs(k,i);
			}

		}
		if ( m_UseIsotropicPart )
			estSi += m_InternalVoxelData->normalizedFractions[m_InternalVoxelData->normalizedFractions.size()-1] *  m_InternalVoxelData->SiL_Isotropic[i];

#endif

		residuals[i] = (m_InternalVoxelData->S0*estSi) - m_InternalVoxelData->Si[i];
	}
}

/**********************************************************************************************//**
 * \fn	void MFMTensorEstimateCostFunction::ComputeC_ForRobustEstimation()
 *
 * \brief	Calculates the m_RobustEstimation_C for robust estimation and update the current InternalVoxelData. 
 *
 * \author	Benoit Scherrer
 * \date	February 2010
*************************************************************************************************/
void MFMTensorEstimateCostFunction::ComputeC_ForRobustEstimation()
{
	std::vector<double>						listErrors;
	

	//---------------------------------------------
	// Compute all the residuals
	//---------------------------------------------
	ComputeResiduals(listErrors);
	
	//---------------------------------------------
	// Compute their absolute value
	//---------------------------------------------
	for ( unsigned int i=0; i<listErrors.size(); i++ )
	{
		if ( listErrors[i]<0 ) listErrors[i] = -listErrors[i];
	}

	//---------------------------------------------
	// Sort these residuals and take the median
	//---------------------------------------------
	std::sort(listErrors.begin(), listErrors.end());
	double median = listErrors[(int)(listErrors.size()/2)];
	double C = m_AlgoParameters->m_RobustParam*median; //1.48*median; //(1.48*median);
	
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
	this->m_InternalVoxelData->RobustEstimation_C = C;

}




/**********************************************************************************************//**
 * \fn	MFMTensorEstimateCostFunction::MeasureType MFMTensorEstimateCostFunction::GetValue(const ParametersType &parameters) const
 *
 * \brief	Gets the value of the cost function for the given parameters. 
 *
 * \author	Benoit Scherrer
 * \date	February 2010
 *
 * \param	parameters	The input optimizer parameters. 
 *
 * \return	The value. 
*************************************************************************************************/
MFMTensorEstimateCostFunction::MeasureType
MFMTensorEstimateCostFunction::GetValue(const ParametersType &parameters) const
{
	//_ASSERT(m_InternalVoxelData!=NULL);
	

	//---------------------------------------------------------
	// Compute the SiLs for each tensor (from parameters)
	//---------------------------------------------------------
	// Extract the tensors and fractions from the parameters
	Parameters2TensorAndFractions( m_InternalVoxelData->current_L, m_InternalVoxelData->Fractions, parameters ) ;

	// Actually compute the SiLs from m_InternalVoxelData->current_L, m_InternalVoxelData->current_EigenValuesLE, m_InternalVoxelData->current_EigenVectorsLE
	ComputeSiLs();

	// Now return the value of the function
	double value = 0;
	double C2=m_InternalVoxelData->RobustEstimation_C * m_InternalVoxelData->RobustEstimation_C;

	//---------------------------------------------------------
	// For EACH gradient direction
	//---------------------------------------------------------
	_ASSERT( m_CurrentMaxBValueIndexToUse!=0 );
	for( int i = 0; i< m_CurrentMaxBValueIndexToUse ; i++ )
	{
		double estSi = 0;

#if NEW_FASTSLOWMODEL
		double fast=1;
		double slow=0;
		for ( int k=0 ; k<m_InternalVoxelData->nbTensorAtVoxel ; k++ )
		{
			fast *= exp( m_InternalVoxelData->normalizedFractions[k] * m_InternalVoxelData->SiLs(k,i)  )  ;
			slow += m_InternalVoxelData->normalizedFractions[k] *exp( m_InternalVoxelData->SiLs(k,i)  )  ;
		
			throw itk::ExceptionObject(__FILE__,__LINE__,"ERROR. Code source fast/slow and radial base line not updated", ""); 
			if ( m_EstimationStrategy&ESTIMATION_RADIALBASELINE )
				slow += m_InternalVoxelData->normalizedFractions[k+m_InternalVoxelData->nbTensorAtVoxel] *exp( m_InternalVoxelData->SiLs_Radial(k,i)  )  ;
		}
		if ( m_UseIsotropicPart )
		{
			double fiso = m_InternalVoxelData->normalizedFractions[m_InternalVoxelData->normalizedFractions.size()-1];

			fast *= exp(-fiso*m_AlgoParameters->m_BValues[i] * g_MediumDiffCoeff);
			slow += fiso*exp(-m_AlgoParameters->m_BValues[i] * g_MediumDiffCoeff);
		}

		estSi = (1-m_InternalVoxelData->fastDiffusionFraction)*slow   +   m_InternalVoxelData->fastDiffusionFraction*fast;
#else
		
		for ( int k=0 ; k<m_InternalVoxelData->nbTensorAtVoxel ; k++ )	
		{

			if ( m_EstimationStrategy&ESTIMATION_RADIALBASELINE )
			{		
				// Get the projections g'Cg  and  g'Dg = g' exp(L) g
				double projC = m_InternalVoxelData->current_GradientRadialProjection(k, i);
				double projD = m_InternalVoxelData->current_GradientTensorProjection(k, i);

				//estSi = (1-projC)*m_InternalVoxelData->SiLs(k,i) + projC;

				//double w_parallel = 0.5 + 0.5*exp(projD)/(exp(projD)+exp(projC));
				//estSi = m_InternalVoxelData->normalizedFractions[k] *
				//	(  w_parallel * m_InternalVoxelData->SiLs(k,i)   +   (1-w_parallel) * exp(-m_AlgoParameters->m_BValues[0]*projC) );
				
				// GOOD ONE
				//estSi = m_InternalVoxelData->SiLs(k,i) + exp(projC)/(exp(projD)+exp(projC)) * exp(-m_AlgoParameters->m_BValues[0]*projC );
				
				estSi = m_InternalVoxelData->SiLs(k,i) + projC * exp(-m_AlgoParameters->m_BValues[0]*projC );

			}
			else
			{
				estSi += m_InternalVoxelData->normalizedFractions[k] * m_InternalVoxelData->SiLs(k,i);
			}
			

			if ( m_EstimationStrategy&ESTIMATION_RADIALISO )
			{
				double proj = crl::FastOps::btAb_Asym_3x3(m_InternalVoxelData->current_RadialMatrix[k], m_AlgoParameters->m_g[i] );
				estSi += m_InternalVoxelData->normalizedFractions[k+m_InternalVoxelData->nbTensorAtVoxel] * exp(-m_AlgoParameters->m_BValues[i]*proj);
			}
		}
		if ( m_UseIsotropicPart )
		{
			int id=m_InternalVoxelData->normalizedFractions.size()-1;
			if ( m_EstimationStrategy&ESTIMATION_RADIALISO ) id-=2;
			estSi += m_InternalVoxelData->normalizedFractions[id] *  m_InternalVoxelData->SiL_Isotropic[i];
		}

#endif


		//------------------------
		// Rician Noise Modeling
		//------------------------
		if ( m_AlgoParameters->m_UseRicianNoiseModeling )
		{
			if ( m_UseRobustEstimation )
				throw itk::ExceptionObject(__FILE__,__LINE__,"Rician noise modeling and robust estimation are currently not compatible.", "");

			double Smeasured = m_InternalVoxelData->Si[i];
			double Smodel = m_InternalVoxelData->S0*estSi;

			double sigma2 = m_AlgoParameters->m_RicianNoiseVariance;
			#define BESSEL vnl_bessel0 //BesselIzero //
			double pdf = Smeasured/sigma2 * exp( - (Smeasured*Smeasured+Smodel*Smodel)/(2*sigma2) ) * BESSEL(Smeasured*Smodel/sigma2);
			value = value - safe_log(pdf);
		}
		//------------------------
		// Gaussian Noise Modeling
		//------------------------
		else
		{
			double error = (m_InternalVoxelData->S0*estSi) - m_InternalVoxelData->Si[i];
			
			//double error = log(m_InternalVoxelData->S0) + estSi - log(m_InternalVoxelData->Si[i]);
			if ( m_AlgoParameters->m_RobustEstimation )
			{
				double e2=error*error;

				// Geman-McLure
				//value += e2/(e2+C2);

				// Lorentzian
				//value += log(1+0.5*e2/C2);

				value += e2/(1+e2/C2);
			}
			else
				value+=error*error;
		}
	}

	//---------------------------------------------------------
	// Test: add the l1 norm?
	//---------------------------------------------------------
	//VnlEigValueType eVal(3);
	//VnlTensorType expL;
	//for ( int k=0 ; k<m_InternalVoxelData->nbTensorAtVoxel ; k++ )
	//{
	//	for ( int i=0 ; i<3 ; i++ ) 
	//		eVal[i] = exp(m_InternalVoxelData->current_EigenValuesLE[k][i]);
	//	
	//	crl::FastOps::AtBA_Bdiag_3x3(expL, m_InternalVoxelData->current_EigenVectorsLE[k], eVal);		// really more fast than	expL = eVecLE.transpose() * eVal * eVecLE;

	//	for ( int i=0; i<3; i++ )
 //		for (int j=0; j<3; j++ )
	//		value += std::fabs(expL[i][j]);
	//}

	//---------------------------------------------------------
	// Add the regularization term
	//---------------------------------------------------------
	if ( m_UseRegularization )
	{
		vnl_vector<double> regtermj(m_InternalVoxelData->nbTensorAtVoxel);
		regtermj.fill(0);

		//long double regterm2=0;
		//long double l0 = sqrtl(eValLE[0][0]*eValLE[0][0]+eValLE[0][1]*eValLE[0][1]+eValLE[0][2]*eValLE[0][2]);
		//long double l1 = sqrtl(eValLE[1][0]*eValLE[1][0]+eValLE[1][1]*eValLE[1][1]+eValLE[1][2]*eValLE[1][2]);
		//long double meanL = 0;
		//long double L1 = (l0+l1)/2;
		double K2 = m_AlgoParameters->m_RegularizationAnisotropicK * m_AlgoParameters->m_RegularizationAnisotropicK;
		
		vnl_matrix_fixed<double,3,3> Regj;

		//vnl_matrix_fixed<double,3,3> mjA0, mjA1, mjB0, mjB1, Regj;
		//double  djA0, djA1, djB0, djB1;


		//---------------------------------------------
		// For the three direction X/Y/Z
		//---------------------------------------------
		for ( int i=0; i<3 ; i++ )
		{
			//---------------------------------------------
			// Gets the two neighbors
			//---------------------------------------------	
			InternalVoxelData *nA = m_InternalVoxelData->Neighbors[2*i];
			InternalVoxelData *nB = m_InternalVoxelData->Neighbors[2*i+1];

			//---------------------------------------------
			// Protect the access to the neighbors
			// because we are going to read nA->stabilized_L
			// and nB->stabilized_L
			//---------------------------------------------
			if ( nA==NULL ) nA = m_InternalVoxelData; //continue;
			else nA->mutexAccessL.Lock();

			if ( nB==NULL ) nB = m_InternalVoxelData; //continue;
			else nB->mutexAccessL.Lock();

			//---------------------------------------------
			// Compute the regularization term for that direction for each tensor
			//---------------------------------------------
			for ( int j=0;j<m_InternalVoxelData->nbTensorAtVoxel;j++ )
			{
//!!!tmp todo
				//---------------------------------------- 
				// Compute argmin_Lq || Lq(x-xm) - Lj(x) ||_LE
				// and argmin_Lq || Lq(x+xm) - Lj(x) ||_LE
				//---------------------------------------- 
				crl::DWI::VnlTensorType *LA = &(nA->stabilized_L[0]);
				if ( nA->nbTensorAtVoxel > 1 )
				{
					double minDA = (m_InternalVoxelData->current_L[j] - nA->stabilized_L[0] ).frobenius_norm();
					for ( int jj=1; jj<nA->nbTensorAtVoxel; jj++ )
					{
						double newDA=(m_InternalVoxelData->current_L[j] - nA->stabilized_L[jj] ).frobenius_norm();
						if ( newDA<minDA ) {
							minDA = newDA;
							LA = &(nA->stabilized_L[jj]);
						}
					}
				}

				crl::DWI::VnlTensorType *LB = &(nB->stabilized_L[0]);
				if ( nB->nbTensorAtVoxel > 1 )
				{
					double minDB = (m_InternalVoxelData->current_L[j] - nB->stabilized_L[0] ).frobenius_norm();
					for ( int jj=1; jj<nB->nbTensorAtVoxel; jj++ )
					{
						double newDB=(m_InternalVoxelData->current_L[j] - nB->stabilized_L[jj] ).frobenius_norm();
						if ( newDB<minDB ) {
							minDB = newDB;
							LB = &(nB->stabilized_L[jj]);
						}
					}
				}


				//Regj = (LA - 2.0*m_InternalVoxelData->current_L[j]);
				Regj = 0.5 * ( (*LA) + (-2.0)*m_InternalVoxelData->current_L[j] + (*LB));

		
				double d = (Regj.frobenius_norm());

				regtermj[j] += d*d;
			}

			//---------------------------------------------
			// Unlock the mutex
			//---------------------------------------------
			if ( nA != m_InternalVoxelData ) nA->mutexAccessL.Unlock();
			if ( nB != m_InternalVoxelData ) nB->mutexAccessL.Unlock();
		}
	
		//---------------------------------------------
		// Actually compute the regularization term
		//---------------------------------------------
		for ( unsigned int j=0;j<regtermj.size();j++ )
			value += m_CurrentRegularizationWeight*( sqrt(1+regtermj[j]/K2) );
	} // end if ( m_UseRegularization )



	//---------------------------------------------------------
	// Add the constraint on the largest eigen value
	// (only if not bobyqa, because bobyqa performs 
	//  bounded optimization)
	//---------------------------------------------------------
	if ( m_AlgoParameters->m_OptimizationAlgo!=crl::OptimizerUtils::OPTIMIZER_BOBYQA )
	{
		for ( int j=0;j<m_InternalVoxelData->nbTensorAtVoxel;j++ )
		{
			double x = (m_InternalVoxelData->OneTensorEigenValues[2])/(m_InternalVoxelData->current_EigenValuesLE[j][2]);
			value += (x-1)*(x-1)/(0.5*0.5);
			//if ( x>5 ) x=5;
			//value += exp( 10*(x-1) );
		}
	}
	
	//---------------------------------------------------------
	// Add constraint on the angle depending on the FA
	//---------------------------------------------------------
	//if ( m_NumberOfTensors==2 )
	//{
	//	vnl_vector<double> v1(3), v2(3);

	//	// Take the primary eigen vectors of the two tensors
	//	v1[0] = m_InternalVoxelData->current_EigenVectorsLE[0][2][0];
	//	v1[1] = m_InternalVoxelData->current_EigenVectorsLE[0][2][1];
	//	v1[2] = m_InternalVoxelData->current_EigenVectorsLE[0][2][2];

	//	v2[0] = m_InternalVoxelData->current_EigenVectorsLE[1][2][0];
	//	v2[1] = m_InternalVoxelData->current_EigenVectorsLE[1][2][1];
	//	v2[2] = m_InternalVoxelData->current_EigenVectorsLE[1][2][2];

	//	// Get the angle between 0 pi/2  (symmetry)
	//	double alpha = angle<double>(v1, v2);
	//	if ( alpha>M_PI/2 ) alpha = M_PI-alpha;

	//	// Get the FA of the one tensor solution
	//	double l1 = m_InternalVoxelData->OneTensorEigenValues[0];
	//	double l2 = m_InternalVoxelData->OneTensorEigenValues[1];
	//	double l3 = m_InternalVoxelData->OneTensorEigenValues[2];
	//	double fa = sqrt ((l1-l2)*(l1-l2)+(l1-l3)*(l1-l3)+(l2-l3)*(l2-l3)) / sqrt(2*(l1*l1+l2*l2+l3*l3));

	//	static const double stdFA = 1;
	//	static const double stdAngle = 1;

	//	double constraint = 1 - exp(-(alpha-M_PI/2)*(alpha-M_PI/2)/(stdAngle*stdAngle)) * exp(-(fa-1)*(fa-1)/(stdFA*stdFA));
	//	value -= 100*safe_log(constraint);

	//}


	//---------------------------------------------------------
	// Test Constraint on the FA ?
	//---------------------------------------------------------
	/*for ( int j=0;j<m_NumberOfTensors;j++ )
	{
		double l1 = exp(m_InternalVoxelData->current_EigenValuesLE[j][0]);
		double l2 = exp(m_InternalVoxelData->current_EigenValuesLE[j][1]);
		double l3 = exp(m_InternalVoxelData->current_EigenValuesLE[j][2]);

		double fa = sqrt ((l1-l2)*(l1-l2)+(l1-l3)*(l1-l3)+(l2-l3)*(l2-l3)) / sqrt(2*(l1*l1+l2*l2+l3*l3));
		
		// fa should be at least 0.3
		value += 10*m_InternalVoxelData->Si.size()*exp(-(fa*fa)/(0.3*0.3));
	}*/

	//---------------------------------------------------------
	// Some verbose?
	//---------------------------------------------------------
	if ( m_AlgoParameters->m_VerboseLevel>=3 ) 
	{
		std::cout<<"f="<<value<<"  ";
		for( unsigned int i = 0; i< parameters.size() ; i++ )
		{
			std::cout << parameters[i]<<"  ";
		}
		std::cout<<std::endl;
	}

	//---------------------------------------------------------
	// Return the value
	//---------------------------------------------------------
	if (isnan(value) ) value = std::numeric_limits<double>::infinity();
	return value ;
}

/**********************************************************************************************//**
 * \fn	void MFMTensorEstimateCostFunction::GetDerivative( const ParametersType &parameters,
 * 		DerivativeType &derivative ) const
 *
 * \brief	Gets the derivative of the cost function. 
 *
 * \author	Benoit Scherrer
 * \date	July 2010
 *
 * \exception	itk::ExceptionObject	Thrown when exception. 
 *
 * \param	parameters			Options for controlling the operation. 
 * \param [in,out]	derivative	The derivative. 
*************************************************************************************************/
void MFMTensorEstimateCostFunction::GetDerivative( const ParametersType &parameters, DerivativeType &derivative ) const 
{
	// Transform the parameters to a vnl matrix
	VnlTensorType	 			D;
	VnlTensorType		  		DeltaSims[2];

	throw itk::ExceptionObject(__FILE__,__LINE__,"Derivative computation to be checked. Probably doesn t work anymore.", "");

	if ( this->m_UseRegularization )
		throw itk::ExceptionObject(__FILE__,__LINE__,"Gradient descent and regularization are currently not compatible.", "");
	if ( m_AlgoParameters->m_UseRicianNoiseModeling )
		throw itk::ExceptionObject(__FILE__,__LINE__,"Rician noise modeling and gradient descent are currently not compatible.", "");


	// Extract the tensors and fractions from the parameters
	Parameters2TensorAndFractions( m_InternalVoxelData->current_L, m_InternalVoxelData->Fractions, parameters ) ;

	// Actually compute the SiLs from m_InternalVoxelData->current_L, m_InternalVoxelData->current_EigenValuesLE, m_InternalVoxelData->current_EigenVectorsLE
	ComputeSiLs();

	// Some init for the estimation of a tensor
	VnlTensorType RGRt(3,3), deriv(3,3), M(3,3);
	for ( int i=0; i<m_InternalVoxelData->nbTensorAtVoxel ; i++ )
	{
		DeltaSims[i].set_size(3,3);
		DeltaSims[i].fill(0);
	}	


	// Some init for the estimation of the fractions
	VnlFractionsType dF;	
	int nbFractionsToEstimate=GetNumberOfFractionsToEstimate();
	if ( nbFractionsToEstimate>0 )
	{
		dF.set_size(nbFractionsToEstimate);
		dF.fill(0);
	}

	_ASSERT( m_CurrentMaxBValueIndexToUse!=0 );
	for( int i = 0; i< m_CurrentMaxBValueIndexToUse ; i++ )
	{
		double Si=0;
		for ( int k=0 ; k<m_InternalVoxelData->nbTensorAtVoxel ; k++ )
			Si += m_InternalVoxelData->normalizedFractions[k]*m_InternalVoxelData->SiLs(k,i);

		if ( m_UseIsotropicPart )
			Si += m_InternalVoxelData->normalizedFractions[m_InternalVoxelData->normalizedFractions.size()-1] *  m_InternalVoxelData->SiL_Isotropic[i];

		double error = m_InternalVoxelData->S0*Si - m_InternalVoxelData->Si[i];
		double weight;

		if ( m_AlgoParameters->m_RobustEstimation )
		{
			double C2=m_InternalVoxelData->RobustEstimation_C * m_InternalVoxelData->RobustEstimation_C;
			//Geman-McLure
			weight = (2*error*C2)/((error*error+C2)*(error*error+C2));

			// Lorentzian
			//weight = 2*weight/(2*C2+weight*weight);


			//weight = 2*weight/((1+weight*weight/C2)*(1+weight*weight/C2));
		}
		else
			weight = 2 * error;

		//---------------------------------------------------------
		//  I. Compute the derivative for all tensors
		//---------------------------------------------------------
		for ( int tensor=0; tensor<m_InternalVoxelData->nbTensorAtVoxel ; tensor++ )
		{
			//---------------------------------------------------------
			// Compute the exp derivative
			// See Pierre Fillard TMI 2007	
			// Compute RGRt
			//---------------------------------------------------------
			vnl_fastops::ABAt(RGRt, m_InternalVoxelData->current_EigenVectorsLE[tensor], m_AlgoParameters->m_G[i] );
			const double *p = m_InternalVoxelData->current_EigenValuesLE[tensor].data_block();
			double sm, sl;

			for ( int l=0 ; l<3 ; l++ )
				for ( int m=0 ; m<3 ; m++ )
				{
					sm = *(p+m);
					sl = *(p+l);

					if ( sm!=sl )
						M(l,m) = RGRt(l,m) * (exp(sm) - exp(sl))/(sm-sl);
					else
						M(l,m) = RGRt(l,m) * exp(sm);
				}

			crl::FastOps::AtBA_3x3(deriv, m_InternalVoxelData->current_EigenVectorsLE[tensor], M );

			DeltaSims[tensor] = DeltaSims[tensor] + weight * m_AlgoParameters->m_BValues[i] * m_InternalVoxelData->SiLs(tensor,i)*deriv;
		}

		//---------------------------------------------------------
		// II. Compute the derivative for the fractions
		//---------------------------------------------------------
	
		// If isotropic part, add the last term
		if ( nbFractionsToEstimate>0 )
		{
			EstimationStrategyType estimF = m_EstimationStrategy&(ESTIMATION_F_ISO+ESTIMATION_F_TENSORS);
			if ( (estimF==ESTIMATION_F_ISO) && (m_InternalVoxelData->nbTensorAtVoxel==2) )
			{
				double e1 = exp(m_InternalVoxelData->Fractions[0]);
				double e3 = exp(m_InternalVoxelData->Fractions[2]);
				double sum_e = e1+e1+e3;
				double sum_e2 = sum_e*sum_e;
				
				dF[0] += weight * m_InternalVoxelData->S0  * ( (e1*sum_e - 2*e1*e1)/sum_e2 * ( m_InternalVoxelData->SiLs(0,i) + m_InternalVoxelData->SiLs(1,i) ) - 2*e3*e1/sum_e2 * m_InternalVoxelData->SiL_Isotropic[i] );
				dF[1] += weight * m_InternalVoxelData->S0  * ( -e1*e3/sum_e2 * ( m_InternalVoxelData->SiLs(0,i)  + m_InternalVoxelData->SiLs(1,i)) + (e3*sum_e-e3*e3)/sum_e2 * m_InternalVoxelData->SiL_Isotropic[i]   );
			}

			else if ( (estimF==(ESTIMATION_F_ISO+ESTIMATION_F_TENSORS)) && (m_InternalVoxelData->nbTensorAtVoxel==2) )
			{
				double e1 = exp(m_InternalVoxelData->Fractions[0]);
				double e2 = exp(m_InternalVoxelData->Fractions[1]);
				double e3 = exp(m_InternalVoxelData->Fractions[2]);
				double sum_e = e1+e2+e3;
				double sum_e2 = sum_e*sum_e;

				dF[0] += weight * m_InternalVoxelData->S0  * ( (e1*sum_e - e1*e1)/sum_e2  * m_InternalVoxelData->SiLs(0,i) - e2*e1/sum_e2 * m_InternalVoxelData->SiLs(1,i) - e3*e1/sum_e2 * m_InternalVoxelData->SiL_Isotropic[i] );
				dF[1] += weight * m_InternalVoxelData->S0  * ( -e1*e2/sum_e2 * m_InternalVoxelData->SiLs(0,i)  +  (e2*sum_e-e2*e2)/sum_e2 * m_InternalVoxelData->SiLs(1,i) - e3*e2/sum_e2 * m_InternalVoxelData->SiL_Isotropic[i]   );
				dF[2] += weight * m_InternalVoxelData->S0  * ( -e1*e3/sum_e2 * m_InternalVoxelData->SiLs(0,i)  -  e2*e3/sum_e2 * m_InternalVoxelData->SiLs(1,i) - (e3*sum_e-e3*e3)/sum_e2 * m_InternalVoxelData->SiL_Isotropic[i]   );

			}
	
			// NO isotropic part, only the fractions
			else if ( (estimF==ESTIMATION_F_TENSORS) && (m_InternalVoxelData->nbTensorAtVoxel==2) )
			{
				double e1 = exp(m_InternalVoxelData->Fractions[0]);
				double e2 = exp(m_InternalVoxelData->Fractions[1]);
				double sum_e = e1+e2;
				double sum_e2 = sum_e*sum_e;

				dF[0] += weight * m_InternalVoxelData->S0  * ( (e1*sum_e - e1*e1)/sum_e2  * m_InternalVoxelData->SiLs(0,i) - e1*e2/sum_e2 * m_InternalVoxelData->SiLs(1,i));
				dF[1] += weight * m_InternalVoxelData->S0  * ( -e1*e2/sum_e2 * m_InternalVoxelData->SiLs(0,i)  +  (e2*sum_e-e2*e2)/sum_e2 * m_InternalVoxelData->SiLs(1,i)    );
			}

			else
				throw itk::ExceptionObject(__FILE__,__LINE__,"Gradient descent is currently only compatible with exactlty 2 tensors.", "");

		}
	}
	
	// 
	for ( int tensor=0; tensor<m_InternalVoxelData->nbTensorAtVoxel ; tensor++ )
		DeltaSims[tensor] = DeltaSims[tensor] * (- m_InternalVoxelData->normalizedFractions[tensor]*m_InternalVoxelData->S0 );
	

	// Now convert result to parameters
	derivative.SetSize(GetNumberOfParameters());
	int offset=0;
	for ( int i=0; i<m_InternalVoxelData->nbTensorAtVoxel ; i++ )
	{

		derivative[offset]   = DeltaSims[i](0,0);
		derivative[offset+1] = DeltaSims[i](0,1);
		derivative[offset+2] = DeltaSims[i](0,2);
		derivative[offset+3] = DeltaSims[i](1,1);
		derivative[offset+4] = DeltaSims[i](1,2);
		derivative[offset+5] = DeltaSims[i](2,2);
		offset+=6;

	}
	if ( nbFractionsToEstimate>0 )
	{
		for ( int i=0 ; i<nbFractionsToEstimate ; i++ )
			derivative[offset+i] = dF[i];
	}

	if ( m_AlgoParameters->m_VerboseLevel>=3 )
	{
		std::cout<<"dF = ";
		for ( unsigned int i=0; i<derivative.size() ; i++ ) std::cout<<derivative[i]<<"  ";
		std::cout<<std::endl;
	}

	if ( m_AlgoParameters->m_VerboseLevel==4 )
	{
		double r=1.1;
		for ( long int i=0; i<2000000;i++) { r += (float)rand()/458998; }
	}

}

/**********************************************************************************************//**
 * \fn	vnl_matrix<double> MFMTensorEstimateCostFunction::RotationMatrix(double angle,
 * 		int r) const
 *
 * \brief	Creates a rotation matrix of a given angle around the given axe. 
 *
 * \author	Benoit Scherrer
 * \date	July 2010
 *
 * \exception	itk::ExceptionObject	Thrown when exception. 
 *
 * \param	angle	The angle. 
 * \param	r		The. 
 *
 * \return	. 
*************************************************************************************************/
VnlMatrixType MFMTensorEstimateCostFunction::RotationMatrix(double angle, int r) const
{
	VnlMatrixType m(3,3);

	if ( r==0 )
	{
		m(0,0)=1;
		m(0,1)=0;
		m(0,2)=0;
		m(1,0)=0;
		m(1,1)=cos(angle);
		m(1,2)=sin(angle);
		m(2,0)=0;
		m(2,1)=-sin(angle);
		m(2,2)=cos(angle);
	}
	else if ( r==1)
	{
		m(0,0)=cos(angle);
		m(0,1)=0;
		m(0,2)=sin(angle);
		m(1,0)=0;
		m(1,1)=1;
		m(1,2)=0;
		m(2,0)=-sin(angle);
		m(2,1)=0;
		m(2,2)=cos(angle);
	}
	else if ( r==2 )
	{
		m(0,0)=cos(angle);
		m(0,1)=sin(angle);
		m(0,2)=0;
		m(1,0)=-sin(angle);
		m(1,1)=cos(angle);
		m(1,2)=0;
		m(2,0)=0;
		m(2,1)=0;
		m(2,2)=1;
	}
	else
		throw itk::ExceptionObject(__FILE__, __LINE__, "Invalid value of r","MFMTensorEstimateCostFunction::RotationMatrix");
	
	return m;
}

/**********************************************************************************************//**
 * \fn	void MFMTensorEstimateCostFunction::SetAlgoParameters ( crl::DWI::GlobalAlgoParameters *p )
 *
 * \brief	Sets the algo parameters. 
 *
 * \author	Benoit Scherrer
 * \date	July 2010
 *
 * \param [in,out]	p	If non-null, the. 
*************************************************************************************************/
void MFMTensorEstimateCostFunction::SetAlgoParameters ( crl::DWI::GlobalAlgoParameters *p ) 
{ 
	m_AlgoParameters = p; 

	// Copy the contents of some recurrently used parameters
	// (faster that m_AlgoParameters->XX ?? not sure...)
	// (but at least more easy to read)

	m_EstimationStrategy	= m_AlgoParameters->m_EstimationStrategy;

	// to directly and more simply do some tests
	m_UseIsotropicPart			= ((m_EstimationStrategy&(ESTIMATION_F_ISO))!=0); 
}

void MFMTensorEstimateCostFunction::SetInternalVoxelData ( InternalVoxelData *d ) 
{ 
	//----------------------------------------
	// If only one tensor for that voxel, and we don't estimate
	// the isotropic diffusion, just force to NO fractions to estimate!
	//----------------------------------------
	m_EstimationStrategy	= m_AlgoParameters->m_EstimationStrategy;
	if ( (d!=NULL) && (d->nbTensorAtVoxel==1) && (!m_UseIsotropicPart) && (m_EstimationStrategy&ESTIMATION_F_TENSORS) )
		m_EstimationStrategy = m_EstimationStrategy - ESTIMATION_F_TENSORS ;

	m_InternalVoxelData = d;
	InitComputeNumberOfFractionsToEstimate();
	InitComputeNumberOfParameters();
} 

void MFMTensorEstimateCostFunction::TestCostFunction(const ParametersType &parameters) const
{
	int NB_X = 50;
	int NB_Y =50;
	int NB_F = 50;
	typedef itk::Image<float,3> TestImageType;

	static int done=0;
	done++;
	if (done!=15) 
	{
		//cout<<"no "<<d<<endl;
		return;
	}

	//-------------------------------------------
	// Creates the output image
	//-------------------------------------------
	TestImageType::IndexType ofIndex;
	TestImageType::SizeType ofSize;
	ofSize[0]=(NB_X);
	ofSize[1]=NB_Y;
	ofSize[2]=NB_F;
	ofIndex[0]=0;
	ofIndex[1]=0;
	ofIndex[2]=0;

	TestImageType::RegionType region;
	region.SetIndex(ofIndex);
	region.SetSize(ofSize);

	TestImageType::Pointer imgTest = TestImageType::New();
	imgTest->SetRegions(region);
	imgTest->Allocate();

	//-------------------------------------------
	// Creates the output iterator
	//-------------------------------------------
	typedef itk::ImageRegionIteratorWithIndex< TestImageType > TestIteratorType;
	TestIteratorType it(imgTest, region);

	//-------------------------------------------
	// Extract the tensors and fractions from the parameters
	//-------------------------------------------
	Parameters2TensorAndFractions( m_InternalVoxelData->current_L, m_InternalVoxelData->Fractions, parameters ) ;

	// Compute the SiLs from m_InternalVoxelData->current_L, m_InternalVoxelData->current_EigenValuesLE, m_InternalVoxelData->current_EigenVectorsLE
	ComputeSiLs();




	if ( m_InternalVoxelData->Fractions[0]<m_InternalVoxelData->Fractions[1] )
	{
		double tmp=m_InternalVoxelData->Fractions[0];
		m_InternalVoxelData->Fractions[0]=m_InternalVoxelData->Fractions[1];
		m_InternalVoxelData->Fractions[1] =tmp;
		
		VnlEigVectorType tmp2 = m_InternalVoxelData->current_EigenVectorsLE[0];
		m_InternalVoxelData->current_EigenVectorsLE[0]=m_InternalVoxelData->current_EigenVectorsLE[1];
		m_InternalVoxelData->current_EigenVectorsLE[1] =tmp2;

		VnlEigValueType tmp3 = m_InternalVoxelData->current_EigenValuesLE[0];
		m_InternalVoxelData->current_EigenValuesLE[0]=m_InternalVoxelData->current_EigenValuesLE[1];
		m_InternalVoxelData->current_EigenValuesLE[1] =tmp3;
	}

	double d=(m_InternalVoxelData->current_L[0]-m_InternalVoxelData->current_L[1]).frobenius_norm();
	cout<<d<<"  f="<<m_InternalVoxelData->Fractions[0]<<endl;
	cout<< exp(m_InternalVoxelData->current_EigenValuesLE[0][2])<<","<< exp(m_InternalVoxelData->current_EigenValuesLE[0][1])<<","<< exp(m_InternalVoxelData->current_EigenValuesLE[0][0])<<endl;
	cout<<m_InternalVoxelData->current_EigenVectorsLE[0][0][2]<<","<<m_InternalVoxelData->current_EigenVectorsLE[0][1][2]<<","<<m_InternalVoxelData->current_EigenVectorsLE[0][2][2]<<"   ";
	cout<<m_InternalVoxelData->current_EigenVectorsLE[1][0][2]<<","<<m_InternalVoxelData->current_EigenVectorsLE[1][1][2]<<","<<m_InternalVoxelData->current_EigenVectorsLE[1][2][2]<<"   ";
	
	cout<<"go"<<endl;
	


	//-------------------------------------------
	// Defines the positions where to evaluate the cost function
	//-------------------------------------------
	float MIN_X = -3.14157/2;
	float MAX_X=3.14157/2;
	float MIN_Y=0.1;
	float MAX_Y=5;
	float MIN_F=0.01;
	float MAX_F=0.99;

	int mode=2;

	if (mode==1)
	{
		MIN_X = -3.14157/2;
		MAX_X=3.14157/2;
		MIN_Y=0.5;
		MAX_Y=2;
		MIN_F=0.5;
		MAX_F=2;
	}
	else if ( mode==2)
	{
		MIN_X = safe_log(1.55e-3 - 5e-4);
		MAX_X = safe_log(1.55e-3 + 5e-4);
		MIN_Y = safe_log(1.39e-3 - 5e-4);
		MAX_Y = safe_log(1.39e-3 + 5e-4);
		MIN_F= 0.65;
		MAX_F= 0.75;

		/*MIN_X = safe_log(exp(m_InternalVoxelData->current_EigenValuesLE[0][2]) - 1e-4);
		MAX_X = safe_log(exp(m_InternalVoxelData->current_EigenValuesLE[0][2]) + 1e-4);
		MIN_Y = safe_log(exp(m_InternalVoxelData->current_EigenValuesLE[1][2]) - 1e-4);
		MAX_Y = safe_log(exp(m_InternalVoxelData->current_EigenValuesLE[1][2]) + 1e-4);
		MIN_F= m_InternalVoxelData->Fractions[0]-0.02;
		MAX_F= m_InternalVoxelData->Fractions[0]+0.02;*/
	}
	//-------------------------------------------
	// Initialize some variables..
	//-------------------------------------------

	VnlTensorType			L2_0(3,3),L2_1(3,3);
	VnlEigVectorType  		eVec2_0,eVec2_1, r1(3,3),r2(3,3);
	VnlEigValueType			eVal2_0, eVal2_1;

	eVec2_0.set_size(3,3);
	eVal2_0.set_size(3);
	eVec2_1.set_size(3,3);
	eVal2_1.set_size(3);	
	// Precompute in m_InternalVoxelData->SiLs the exp(-b gt expL g) for the current tensor for each direction

	eVal2_0 = m_InternalVoxelData->current_EigenValuesLE[0];
	eVal2_1 = m_InternalVoxelData->current_EigenValuesLE[1];

	// Compute D = exp(L) for this tensor 
	//EigenAnalysis.ComputeEigenValuesAndVectors(m_InternalVoxelData->current_L[0], eVal2_0, eVec2_0 ); 
	//EigenAnalysis.ComputeEigenValuesAndVectors(m_InternalVoxelData->current_L[1], eVal2_1, eVec2_1 ); 

	//-------------------------------------------
	// Computes for each output voxel!
	//-------------------------------------------
	it.GoToBegin();	
	while ( !it.IsAtEnd() )
	{
		TestImageType::IndexType id = it.GetIndex();

		if (mode==0)
		{
			r1 = RotationMatrix(id[0]*(MAX_X-MIN_X)/NB_X + MIN_X, 0);
			r2 = RotationMatrix(id[1]*(MAX_X-MIN_X)/NB_X + MIN_X, 0);
			//float zz=id[1]*(MAX_Y-MIN_Y)/NB_Y + MIN_Y;

			crl::FastOps::AtBA_Bsym_3x3(L2_0, r1, eVal2_0) ;
			crl::FastOps::AtBA_Bsym_3x3(L2_1, r2, eVal2_1) ;

			throw itk::ExceptionObject(__FILE__,__LINE__,"WARNING. CODE HAS TO BE MODIFIED SINCE THE MODIF OF ComputeSils");
			//ComputeSiLs(0, L2_0, eValLE[0], eVecLE[0]) ;
			//ComputeSiLs(1, L2_1, eValLE[0], eVecLE[0]) ;
			
			m_InternalVoxelData->Fractions[0] = id[2]*(MAX_F-MIN_F)/NB_F + MIN_F;
			m_InternalVoxelData->Fractions[1]=1-m_InternalVoxelData->Fractions[0];
			m_InternalVoxelData->normalizedFractions[0] = m_InternalVoxelData->Fractions[0];
			m_InternalVoxelData->normalizedFractions[1]=m_InternalVoxelData->Fractions[1];
		}
		else if (mode==1)
		{
			r1 = RotationMatrix(id[0]*(MAX_X-MIN_X)/NB_X + MIN_X, 0);

			double v1 = id[1]*(MAX_Y-MIN_Y)/NB_Y + MIN_Y;
			double v2 = id[2]*(MAX_F-MIN_F)/NB_F + MIN_F;

			r1[1][1]=v1*r1[1][1];
			r1[2][2]=v2*r1[2][2];

			crl::FastOps::AtBA_Bsym_3x3(L2_0, r1, eVal2_0) ;

			throw itk::ExceptionObject(__FILE__,__LINE__,"WARNING. CODE HAS TO BE MODIFIED SINCE THE MODIF OF ComputeSils");
			//ComputeSiLs(0, L2_0, eValLE[0], eVecLE[0]) ;
		}
		else if ( mode==2 )
		{
			m_InternalVoxelData->current_EigenValuesLE[0][2] = id[0]*(MAX_X-MIN_X)/NB_X + MIN_X;
			m_InternalVoxelData->current_EigenValuesLE[1][2] = id[1]*(MAX_Y-MIN_Y)/NB_Y + MIN_Y;

			m_InternalVoxelData->Fractions[0] = id[2]*(MAX_F-MIN_F)/NB_F + MIN_F;
			m_InternalVoxelData->Fractions[1]=1-m_InternalVoxelData->Fractions[0];
			m_InternalVoxelData->normalizedFractions[0] = m_InternalVoxelData->Fractions[0];
			m_InternalVoxelData->normalizedFractions[1]=m_InternalVoxelData->Fractions[1];

			ComputeSiLs();

		}

		// Now return the value of the function
		double value = 0;

		//---------------------------------------------------------
		// For EACH gradient direction
		//---------------------------------------------------------
		_ASSERT( m_CurrentMaxBValueIndexToUse!=0 );
		for( int i = 0; i< m_CurrentMaxBValueIndexToUse ; i++ )
		{
			double estSi = 0;
			for ( int k=0 ; k<m_InternalVoxelData->nbTensorAtVoxel ; k++ )	estSi += m_InternalVoxelData->normalizedFractions[k] * m_InternalVoxelData->SiLs(k,i);

			if ( m_UseIsotropicPart )
				estSi += m_InternalVoxelData->normalizedFractions[m_InternalVoxelData->normalizedFractions.size()-1] *  m_InternalVoxelData->SiL_Isotropic[i];

			double error = (m_InternalVoxelData->S0*estSi) - m_InternalVoxelData->Si[i];
			value+=error*error;
		}

		it.Set(value/m_CurrentMaxBValueIndexToUse);	//divide by m_CurrentMaxBValueIndexToUse to compute the mean error over the directions

		++it;
	}

	//-------------------------------------------
	// Save the image
	//-------------------------------------------
	char szName[1024];
	sprintf(szName,"testcostfunction_%d.nrrd",(int)m_AlgoParameters->m_BValues[m_AlgoParameters->m_BValues.size()-1]);

	typedef itk::ImageFileWriter<TestImageType > TestWriterType;
	TestWriterType::Pointer testWriter = TestWriterType::New();
	testWriter->SetFileName(szName );
	testWriter->SetInput( imgTest );
	testWriter->Update();

	//-------------------------------------------
	// Save the gradient of the image
	//-------------------------------------------
	typedef itk::GradientMagnitudeImageFilter< TestImageType,TestImageType > GradFilterType;
	GradFilterType::Pointer grad=GradFilterType::New();
	grad->SetInput(imgTest);
	grad->Update();
	sprintf(szName,"testcostfunction_grad%d.nrrd",(int)m_AlgoParameters->m_BValues[m_AlgoParameters->m_BValues.size()-1]);

	testWriter->SetFileName(szName );
	testWriter->SetInput( grad->GetOutput() );
	testWriter->Update();

	cout<<"bip bip"<<endl;
}

} // end of namespace DWI
} // end of namespace crl

#endif
