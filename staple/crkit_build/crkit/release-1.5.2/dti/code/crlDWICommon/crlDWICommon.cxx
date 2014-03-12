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

#include "crlDWICommon.h"

namespace crl {
namespace DWI {

void Tensor2VnlMatrix(vnl_matrix<double> &out, const TensorPixelType &tensor ) 
{
	const double *p = tensor.GetDataPointer ();
	out(0,0) = *p ;
	out(0,1) = *(p+1);  //tensor[1];
	out(0,2) = *(p+2);  //tensor[2];
	out(1,0) = *(p+1);  //tensor[1];
	out(1,1) = *(p+3);  //tensor[3];
	out(1,2) = *(p+4);  //tensor[4];
	out(2,0) = *(p+2);  //tensor[2];
	out(2,1) = *(p+4);  //tensor[4];
	out(2,2) = *(p+5);  //tensor[5];
}

void Tensor2VnlMatrix(vnl_matrix<double> &out, const FloatTensorPixelType &tensor ) 
{
	const float *p = tensor.GetDataPointer ();
	out(0,0) = *p ;
	out(0,1) = *(p+1);  //tensor[1];
	out(0,2) = *(p+2);  //tensor[2];
	out(1,0) = *(p+1);  //tensor[1];
	out(1,1) = *(p+3);  //tensor[3];
	out(1,2) = *(p+4);  //tensor[4];
	out(2,0) = *(p+2);  //tensor[2];
	out(2,1) = *(p+4);  //tensor[4];
	out(2,2) = *(p+5);  //tensor[5];
}


void VnlMatrix2Tensor( TensorPixelType &out , const vnl_matrix<double> &m ) 
{
	double *p = out.GetDataPointer ();
	*p = m(0,0);
	*(p+1) =  m(0,1);
	*(p+2) = m(0,2);
	*(p+3) = m(1,1);
	*(p+4) = m(1,2);
	*(p+5) = m(2,2);
}

void VnlMatrix2Tensor( FloatTensorPixelType &out , const vnl_matrix<double> &m ) 
{
	float *p = out.GetDataPointer ();
	*p = (float)m(0,0);
	*(p+1) =  (float)m(0,1);
	*(p+2) = (float)m(0,2);
	*(p+3) = (float)m(1,1);
	*(p+4) = (float)m(1,2);
	*(p+5) = (float)m(2,2);
}


vnl_matrix<double> GetVnlRotationMatrix(int ax, double alpha)
{
	vnl_matrix<double> m(3,3);

	if ( ax==0 )
	{
		m(0,0)=1;
		m(0,1)=0;
		m(0,2)=0;
		m(1,0)=0;
		m(1,1)=cos(alpha);
		m(1,2)=-sin(alpha);
		m(2,0)=0;
		m(2,1)=sin(alpha);
		m(2,2)=cos(alpha);
	}
	else if ( ax==1 )
	{
		m(0,0)=cos(alpha);
		m(0,1)=0;
		m(0,2)=sin(alpha);
		m(1,0)=0;
		m(1,1)=1;
		m(1,2)=0;
		m(2,0)=-sin(alpha);
		m(2,1)=0;
		m(2,2)=cos(alpha);
	}
	else
	{
		m(0,0)=cos(alpha);
		m(0,1)=-sin(alpha);
		m(0,2)=0;
		m(1,0)=sin(alpha);
		m(1,1)=cos(alpha);
		m(1,2)=0;
		m(2,0)=0;
		m(2,1)=0;
		m(2,2)=1;
	}

	return m;
}

vnl_matrix<double> GetVnlRotationMatrixPhiTheta(double theta, double phi)
{
	vnl_matrix<double> m(3,3);

	m(0,0)=sin(theta)*cos(phi);
	m(0,1)=sin(theta)*sin(phi);
	m(0,2)=cos(theta);

	m(1,0)=cos(theta)*cos(phi);
	m(1,1)=cos(theta)*sin(phi);
	m(1,2)=-sin(theta);

	m(2,0)=-sin(phi);
	m(2,1)=cos(phi);
	m(2,2)=0;

	return m;
}

void eucl2logeucl(SymmetricEigenAnalysisType &evAnalysis, vnl_matrix<double> &L, const vnl_matrix<double> &D )
{
	vnl_matrix<double>		eVec(3,3);
	vnl_diag_matrix<double> eVal(3) ;

	evAnalysis.ComputeEigenValuesAndVectors(D, eVal, eVec);
	
	// Put it in logeucl
	for ( int i=0 ; i<3 ; i++ )
		eVal[i] = ((eVal[i]>1e-15 ? log(eVal[i]):log(1e-15)));

	L = eVec.transpose() * eVal * eVec;
}


void logeucl2eucl(SymmetricEigenAnalysisType &evAnalysis, vnl_matrix<double> &D, const vnl_matrix<double> &L )
{
	vnl_matrix<double>		eVec(3,3);
	vnl_diag_matrix<double> eVal(3) ;

	evAnalysis.ComputeEigenValuesAndVectors(L, eVal, eVec);
	
	// Put it in logeucl
	for ( int i=0 ; i<3 ; i++ )
		eVal[i] = exp(eVal[i]);

	D = eVec.transpose() * eVal * eVec;
}


} // end namespace DWI.
} // end namespace crl.

//***********************************************************************
