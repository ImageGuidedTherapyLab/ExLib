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

#ifndef __crlCCSegmentFilter_h
#define __crlCCSegmentFilter_h

#include "itkImageToImageFilter.h"
#include "itkDiffusionTensor3D.h"
#include "itkImage.h"
#include "itkSymmetricEigenAnalysis.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_diag_matrix.h"

#include "crlDWIStudy.h"
#include "crlOptimizerUtils.h"

using namespace itk;

namespace crl { 


	/** \class CCSegmentFilter
	*/

	class ITK_EXPORT CCSegmentFilter :
		public ImageToImageFilter< itk::Image< itk::DiffusionTensor3D< float > , 3 > , itk::Image< unsigned short,3>  >
	{
	public:
		typedef CCSegmentFilter								Self;
		typedef SmartPointer<Self>                          Pointer;
		typedef SmartPointer<const Self>                    ConstPointer;
		typedef ImageToImageFilter< itk::Image< itk::DiffusionTensor3D< float > , 3 > , itk::Image< unsigned short , 3 >  > Superclass;

		/** Method for creation through the object factory. */
		itkNewMacro(Self);  

		/** Runtime information support. */
		itkTypeMacro(CCSegmentFilter,	ImageToImageFilter);

		typedef itk::DiffusionTensor3D< float >			TensorPixelType;
		typedef itk::Image< TensorPixelType, 3 >		TensorImageType;
		typedef itk::Image< unsigned short, 3 >			OutputImageType;

		typedef crl::DWIStudy<unsigned short>			StudyType;

		typedef StudyType::DWIGradientImageSetType		DWIGradientImageSetType;

		typedef TensorImageType			 				InputImageType;

		typedef vnl_matrix<double>						VnlTensorType;
		typedef vnl_diag_matrix<double>					VnlEigValueType;
		typedef vnl_matrix<double>						VnlEigVectorType;

		typedef itk::SymmetricEigenAnalysis< VnlTensorType , VnlEigValueType, VnlEigVectorType > SymmetricEigenAnalysisType;

	protected:
		CCSegmentFilter();
		~CCSegmentFilter();
		void PrintSelf(std::ostream& os, Indent indent) const;

			/** Multi-threaded pipeline methods.*/
		virtual void BeforeThreadedGenerateData();
		virtual void ThreadedGenerateData( const OutputImageRegionType &outputRegionForThread, int);
		virtual void AfterThreadedGenerateData ();

	private:
		CCSegmentFilter(const Self&);		//purposely not implemented
		void operator=(const Self&);	//purposely not implemented

	protected:
		TensorImageType*					m_InputTensor;
		OutputImageType::Pointer			m_OutputSegm;
	};
	
}   // end namespace crl

#endif
