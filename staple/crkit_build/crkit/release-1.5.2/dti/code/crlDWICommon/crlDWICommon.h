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

#ifndef CRL_DWI_COMMON_H
#define CRL_DWI_COMMON_H

#include "itkVectorContainer.h"
#include "itkVectorImage.h"
#include "itkDiffusionTensor3D.h"
#include "itkImage.h"
#include "itkSymmetricEigenAnalysis.h"

namespace crl {
namespace DWI {

typedef unsigned short											PixelType;
typedef unsigned short											ReferenceImagePixelType;
typedef unsigned short											GradientImagePixelType ;

typedef itk::VectorImage<unsigned short, 3>							RawDataImageType;

typedef vnl_vector_fixed< double,3 > 								GradientDirectionType;
typedef itk::VectorContainer< unsigned int, GradientDirectionType > GradientDirectionContainerType;
typedef itk::VectorImage< GradientImagePixelType, 3 > 				GradientImagesType;
typedef itk::Image< GradientImagePixelType, 3 > 					GradientImageType;

typedef vnl_matrix_fixed< double, 6, 6 > 							TensorBasisMatrixType;
typedef itk::DiffusionTensor3D< double >							TensorPixelType;
typedef itk::DiffusionTensor3D< float >								FloatTensorPixelType;
typedef itk::Image< TensorPixelType, 3 > 							TensorImageType;

typedef itk::Image< float, 4 > 										FractionsImageType;

typedef itk::SymmetricEigenAnalysis< vnl_matrix<double> , vnl_diag_matrix<double>, vnl_matrix<double> > SymmetricEigenAnalysisType;

// cf http://en.wikipedia.org/wiki/Mass_diffusivity#cite_note-0
const double g_MediumDiffCoeff = 0.003; //0.0016; //0.0007;

//--------------------------------------------
// Useful functions
//--------------------------------------------
void Tensor2VnlMatrix(vnl_matrix<double> &out, const TensorPixelType &tensor ) ;
void VnlMatrix2Tensor( TensorPixelType &out , const vnl_matrix<double> &m ) ;

void Tensor2VnlMatrix(vnl_matrix<double> &out, const FloatTensorPixelType &tensor ) ;
void VnlMatrix2Tensor( FloatTensorPixelType &out , const vnl_matrix<double> &m ) ;


vnl_matrix<double> GetVnlRotationMatrix(int ax, double dg);
vnl_matrix<double> GetVnlRotationMatrixPhiTheta(double theta, double phi);

void logeucl2eucl(SymmetricEigenAnalysisType &evAnalysis, vnl_matrix<double> &D, const vnl_matrix<double> &L );
void eucl2logeucl(SymmetricEigenAnalysisType &evAnalysis, vnl_matrix<double> &L, const vnl_matrix<double> &D );




} // end namespace DWI.
} // end namespace crl.

#endif


