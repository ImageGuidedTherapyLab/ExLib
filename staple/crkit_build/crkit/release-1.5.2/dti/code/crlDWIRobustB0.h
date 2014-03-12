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

#ifndef CRL_DWI_ROBUST_B0_H
#define CRL_DWI_ROBUST_B0_H

#include <itkLightObject.h>
#include <itkImageIOBase.h>

#include "crlGeneralizedInstantiatorMacros.h"
#include "crlDWIStudy.h"
#include "crlFileName.h"
#include "crlVectorImageUtils.h"
#include "crlNHDRFileUtils.h"
#include "crlConvertAndWriteScalarImage.h"
#include "crlRobustMeanImageFilter.h"

using namespace itk;

namespace crl{
	namespace DWI {

		/**********************************************************************************************//**
		 * \class	DWIRobustB0_Base
		 *
		 * \brief	Base class for the instantiator 
		 *
		 * \author	Benoit Scherrer
		 * \date	January 2011
		*************************************************************************************************/
		class ITK_EXPORT DWIRobustB0_Base : public itk::LightObject
		{
		CRL_INSTANCIATOR_BASE(DWIRobustB0_Base);

		public:
			// Algo parameters
			void	SetInputFileName (const std::string& fileName ) { m_InputFileName=fileName; }
			void	SetOutputFileName (const std::string& fileName ) { m_OutputFileName=fileName; }
			void	SetMeanFileName (const std::string& fileName ) { m_MeanFileName=fileName; }
			void	SetDiffFileName (const std::string& fileName ) { m_DiffFileName=fileName; }
			void	SetNewRAWFileName (const std::string& fileName ) { m_NewRAWFileName=fileName; }
			void	SetB0Threshold ( float b0 ) { m_B0Threshold = b0; }
			void	SetNumberOfThreads(int nb) { m_NbProcs = nb; }

			// Constructor/Destructor
			DWIRobustB0_Base():
				m_B0Threshold(10),
				m_NbProcs(1) {}

			~DWIRobustB0_Base() {}

			// Function to implement 
			virtual int Execute() = 0;

		protected:
			std::string		m_InputFileName;
			std::string		m_OutputFileName;
			std::string		m_MeanFileName;
			std::string		m_DiffFileName;
			std::string		m_NewRAWFileName;
			float			m_B0Threshold;
			int				m_NbProcs;
		};

		/**********************************************************************************************//**
		 * \class	DWIRobustB0
		 *
		 * \brief	Main class
		 *
		 * \author	Benoit Scherrer
		 * \date	January 2011
		*************************************************************************************************/
		template <typename TPixel, unsigned int VImageDimension  >
		class ITK_EXPORT DWIRobustB0 : public DWIRobustB0_Base
		{
		CRL_INSTANCIATOR_CLASS(DWIRobustB0_Base, DWIRobustB0);
		public:
			typedef crl::DWIStudy<TPixel> StudyType;

			DWIRobustB0() {}
			~DWIRobustB0() {}

			virtual int Execute()
			{
				// Use crl::FileName for filenames
				crl::FileName	input_image(m_InputFileName);
				crl::FileName	output_image(m_OutputFileName);

				/* --- If not output file specified, just add "_B0" suffix to input filename ---*/
				if ( m_OutputFileName=="" ) output_image.setCompleteFilePath(input_image.getCompleteFilePath_WithSufix("_B0"));
				if ( output_image.getExtension()=="nhdr" ) output_image.setExtension("nrrd");

				//---------------------------------------
				// Load the study
				//---------------------------------------
				StudyType inputStudy;
				inputStudy.LoadStudy(m_InputFileName); 

				//-------------------------------------------------------
				// Instantiate the robust mean filter
				//-------------------------------------------------------
				typedef typename crl::RobustMeanImageFilter< typename StudyType::DWIGradientImageType, typename StudyType::DWIGradientImageType >  RobustMeanImageFilterType;
				typename RobustMeanImageFilterType::Pointer filterRobustMean = RobustMeanImageFilterType::New();
				filterRobustMean->SetNumberOfThreads( m_NbProcs );

				//-------------------------------------------------------
				// Add all B0 images to the filter
				//-------------------------------------------------------
				cout<<"- Extract all B=0 images..."<<endl;
				int nbB0 = 0; 
				for ( int i=0; i<inputStudy.BValues.size(); i++ )
				{
					if ( inputStudy.BValues[i] < m_B0Threshold )
					{
						typename StudyType::DWIGradientImageType::Pointer oneImage = crl::ExtractComponentOfVectorImage< TPixel, 3>( inputStudy.DWIData, i);
						filterRobustMean->SetInput(nbB0, oneImage);
						nbB0++;
					}
				}

				//-------------------------------------------------------
				// Compute the robust mean
				//-------------------------------------------------------
				cout<<"- Compute robust mean..."<<endl;
				filterRobustMean->SetOutputMeanImage(m_MeanFileName!="");
				filterRobustMean->SetOutputDiffImage(m_DiffFileName!="");
				filterRobustMean->Update();

				/*-------------------------------------
				Save the results
				-------------------------------------*/
				typedef typename itk::ImageFileWriter< typename StudyType::DWIGradientImageType > WriterType;
				typename WriterType::Pointer writer = WriterType::New();
				writer->SetUseCompression(true);

				// I. Save the robust estimate
				cout<<"- Save robust mean in <"<<output_image.getCompleteFilePath().c_str()<<">..."<<endl;
				writer->SetFileName( output_image.getCompleteFilePath() );
				writer->SetInput(  filterRobustMean->GetOutput() );
				writer->Update();

				// II. Save the mean if needed
				if ( m_MeanFileName!="" )
				{
					cout<<"- Save non-robust mean in <"<<m_MeanFileName.c_str()<<">..."<<endl;
					writer->SetFileName( m_MeanFileName );
					writer->SetInput(  filterRobustMean->GetMeanImage() );
					writer->Update();
				}

				// III Save the absolute difference between the mean and the robust mean
				if ( m_DiffFileName!="" )
				{
					cout<<"- Save diff-image in <"<<m_DiffFileName.c_str()<<">..."<<endl;
					writer->SetFileName( m_DiffFileName );
					writer->SetInput(  filterRobustMean->GetDiffImage() );
					writer->Update();
				}

				// IV Save he new RAW file
				if ( m_NewRAWFileName!="" )
				{
					cout << "- Create new RAW DWI image in <"<<m_NewRAWFileName.c_str()<<">..."<<endl;

					//-------------------------------------------------------
					// Create the new set of gradient vectors, with ONLY one null
					//-------------------------------------------------------
					typename StudyType::GradientVectorSetType newGVecs;
					typename StudyType::BValueSetType newBVals;
					typename StudyType::GradientVectorType newvec;

					newvec.fill(0);
					newBVals.push_back(0);
					newGVecs.push_back(newvec);

					//-------------------------------------------------------
					// Copy the previous non-null gradients
					//-------------------------------------------------------
					std::vector<int> idx;	//copy the index of non-null gradients
					for ( int i=0; i<inputStudy.BValues.size(); i++ )
					{
						if ( inputStudy.BValues[i] >= m_B0Threshold )
						{
							newGVecs.push_back(inputStudy.Gradients[i]);
							newBVals.push_back(inputStudy.BValues[i]);
							idx.push_back(i);
						}
					}

					//-------------------------------------------------------
					// Create the new vector image
					//-------------------------------------------------------
					typename StudyType::DWIGradientImageSetType::Pointer newImage = StudyType::DWIGradientImageSetType::New();
					newImage->SetRegions(inputStudy.DWIData->GetLargestPossibleRegion());
					newImage->SetOrigin(inputStudy.DWIData->GetOrigin());
					newImage->SetSpacing(inputStudy.DWIData->GetSpacing());
					newImage->SetDirection(inputStudy.DWIData->GetDirection());
					newImage->SetNumberOfComponentsPerPixel(newGVecs.size());
					newImage->Allocate();

					//-------------------------------------------------------
					// Fill the new vector image
					//-------------------------------------------------------
					typedef typename itk::ImageRegionConstIterator<typename StudyType::DWIGradientImageSetType>	InIteratorType;
					typedef typename itk::ImageRegionIterator<typename StudyType::DWIGradientImageSetType>		OutIteratorType;
					typedef typename itk::ImageRegionConstIterator<typename StudyType::DWIGradientImageType>	InB0IteratorType;

					InIteratorType inItr(inputStudy.DWIData, newImage->GetLargestPossibleRegion() );
					OutIteratorType outItr(newImage, newImage->GetLargestPossibleRegion() );
					InB0IteratorType inB0Itr ( filterRobustMean->GetOutput(), newImage->GetLargestPossibleRegion()); 
					inItr.GoToBegin(); inB0Itr.GoToBegin();

					typename StudyType::DWIGradientImageSetType::PixelType inP(inputStudy.DWIData->GetNumberOfComponentsPerPixel() ) ;
					typename StudyType::DWIGradientImageSetType::PixelType outP(newImage->GetNumberOfComponentsPerPixel() ) ;

					//-------------------------------------------------------
					// For all pixels
					//-------------------------------------------------------
					for ( outItr.GoToBegin() ; !outItr.IsAtEnd() ; ++outItr )
					{
						// Copy from the robust B=0
						typename StudyType::DWIGradientImageType::PixelType b0P = inB0Itr.Value();
						outP[0] = b0P;

						// Copy only a part of the vector
						inP = inItr.Get();
						for ( unsigned int i=0; i<idx.size(); i++ )
							outP[i+1] = inP[idx[i]];

						outItr.Set(outP);

						++inItr;
						++inB0Itr;
					}

					//-------------------------------------------------------
					// Create all the new RAW NHDR file
					//-------------------------------------------------------
					inputStudy.BValues = newBVals;
					inputStudy.Gradients = newGVecs;
					inputStudy.WriteStudyWithNewData(m_NewRAWFileName, newImage, false);
				}		
				return 0;
			}
		};




	} // end namespace DWI
} // end namespace crl.

#endif
