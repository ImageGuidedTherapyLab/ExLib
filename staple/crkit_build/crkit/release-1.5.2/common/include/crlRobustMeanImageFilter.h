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

#ifndef CRL_ROBUSTMEAN_IMAGE_FILTER_H
#define CRL_ROBUSTMEAN_IMAGE_FILTER_H

#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "crlRobustMeanCostFunction.h"

using namespace itk;

namespace crl{

		/**********************************************************************************************//**
		 * \class		RobustMeanImageFilter
		 * \brief	 Compute the robust mean of a set of input images
		 *
		*************************************************************************************************/
	template <class TInputImage, class TOutputImage>
	class ITK_EXPORT RobustMeanImageFilter:
		public ImageToImageFilter < TInputImage, TOutputImage >
	{
	public:
		/** Standard class typedefs. */
		typedef RobustMeanImageFilter							Self;
		typedef ImageToImageFilter<TInputImage,TOutputImage>	Superclass;
		typedef SmartPointer <Self>								Pointer;
		typedef SmartPointer<const Self>						ConstPointer;


		/** Method for creation through the object factory. */
		itkNewMacro(Self);

		/** Run-time type information (and related methods). */
		itkTypeMacro(RobustMeanImageFilter,ImageToImageFilter);

		///** Capture the image dimension from the input template parameters. */
		itkStaticConstMacro(ImageDimension, unsigned int,
			TInputImage::ImageDimension);

		/** Convenient typedefs. */
		typedef			TInputImage						InputImageType;
		typedef			TOutputImage					OutputImageType;
		typedef typename InputImageType::ConstPointer	InputImageConstPointer;
		typedef typename OutputImageType::Pointer		OutputImagePointer;
		typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

		/** Set/Get if also compute the mean image \sa GetMeanImage() */
		itkSetMacro(OutputMeanImage, bool);
		itkGetMacro(OutputMeanImage, bool);

		/** Set/Get if also computes the difference between the mean image and the robust mean image \sa GetDiffImage() */
		itkSetMacro(OutputDiffImage, bool);
		itkGetMacro(OutputDiffImage, bool);

		/** Set/Get if it shows the task progress */
		itkSetMacro(ShowTaskProgress, bool);
		itkGetMacro(ShowTaskProgress, bool);


		/**********************************************************************************************//**
		 * \fn	OutputImageType* GetMeanImage()
		 *
		 * \brief	Gets the mean image. (Only valid if GetOutputMeanImage() is true)
		 *
		 * \author	Benoit Scherrer
		 * \date	April 2010
		 *
		 * \return	null if it fails, else the mean image. 
		*************************************************************************************************/
		OutputImageType*	GetMeanImage() { return m_MeanImage.GetPointer(); }

		/**********************************************************************************************//**
		 * \fn	OutputImageType* GetDiffImage()
		 *
		 * \brief	Gets the difference between the mean and the robust mean image (Only valid if
		 * 			GetOutputDiffImage() is true) 
		 *
		 * \author	Benoit Scherrer
		 * \date	April 2010
		 *
		 * \return	null if it fails, else the difference image. 
		*************************************************************************************************/
		OutputImageType*	GetDiffImage() { return m_DiffImage.GetPointer(); }
	
	protected:
		RobustMeanImageFilter();
		~RobustMeanImageFilter();
		void PrintSelf(std::ostream& os, Indent indent) const;

		/** Multi-threaded pipeline method.*/
		void BeforeThreadedGenerateData();
		void ThreadedGenerateData( const OutputImageRegionType &outputRegionForThread, int);

		std::vector<InputImageConstPointer>	m_InputImages;

		OutputImagePointer					m_RobustMeanEstimate; 
		OutputImagePointer					m_MeanImage;
		OutputImagePointer					m_DiffImage;

		bool								m_OutputMeanImage;
		bool								m_OutputDiffImage;
		bool								m_ShowTaskProgress;

	private:
		RobustMeanImageFilter(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented


	};


} // end namespace crl.

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlRobustMeanImageFilter.txx"
#endif

#endif
