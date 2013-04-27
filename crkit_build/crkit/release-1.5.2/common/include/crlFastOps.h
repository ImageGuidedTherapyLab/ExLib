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

#ifndef crl_fastops_h
#define crl_fastops_h

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_diag_matrix.h>

namespace crl {

/**********************************************************************************************//**
 * \class	FastOps
 *
 * \brief	Fast operations
 *
 * \author	Benoit Scherrer
 * \date	August 2010
*************************************************************************************************/
class FastOps
{
 public:
   static void AtBA_3x3(vnl_matrix<double>& out, const vnl_matrix<double>& A, const vnl_matrix<double>& B);
   static void AtBA_3x3(vnl_matrix_ref<double> out, const vnl_matrix<double>& A, const vnl_matrix<double>& B);

   static void AtBA(vnl_matrix<double>& out, const vnl_matrix<double>& A, const vnl_matrix<double>& B);

   static void ggt( vnl_matrix<double> &out, const vnl_vector<double>& g );
 
   static void AtBA_Bsym_3x3(vnl_matrix<double>& out, const vnl_matrix<double>& A, const vnl_matrix<double>& B);
   static void AtBA_Bsym_3x3(vnl_matrix_ref<double>& out, const vnl_matrix<double>& A, const vnl_matrix<double>& B);

   static void AtBA_Bdiag_3x3(vnl_matrix<double>& out, const vnl_matrix<double>& A, const vnl_diag_matrix<double>& B);
   static void AtBA_Bdiag_3x3(vnl_matrix_ref<double> out, const vnl_matrix<double>& A, const vnl_diag_matrix<double>& B);


   static double btAb_Asym_3x3 (const vnl_matrix<double>& A, const vnl_vector<double>& b);
};

} // end of namespace
#endif // vnl_fastops_h_
