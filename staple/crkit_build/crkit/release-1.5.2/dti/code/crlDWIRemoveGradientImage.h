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

#ifndef CRL_DWI_REMOVE_GRADIENT_IMAGE_H
#define CRL_DWI_REMOVE_GRADIENT_IMAGE_H

#include <vector>
#include <stdlib.h>

#include <itkLightObject.h>
#include <itkImageIOBase.h>

#include "crlGeneralizedInstantiatorMacros.h"
#include "crlFileName.h"
#include "crlVectorImageUtils.h"
#include "crlNHDRFileUtils.h"
#include "crlDWIStudy.h"
#include "crlDWICommon.h"

using namespace itk;

namespace crl{
	namespace DWI {

		/**********************************************************************************************//**
		 * \class	RemoveGradientImage_Base
		 *
		 * \brief	Base class for the instantiator 
		 *
		 * \author	Benoit Scherrer
		 * \date	January 2011
		*************************************************************************************************/
		class ITK_EXPORT RemoveGradientImage_Base : public itk::LightObject
		{
		CRL_INSTANCIATOR_BASE(RemoveGradientImage_Base);
		public:
			// Algo parameters
			void	SetInputFileName (const std::string& fileName ) { m_InputFileName=fileName; }
			void	SetOutputFileName (const std::string& fileName ) { m_OutputFileName=fileName; }
			void	SetBMax(double bmax) {m_BMax = bmax;}	// if -1, don't use bmax

			void	ClearListOfIndexes() { m_IndexListToRemove.clear(); }
			void	AddIndexToRemove(int id) { m_IndexListToRemove.push_back(id); }
			bool	IsIndexInList(int id) 
			{ 
				for ( unsigned int i=0; i<m_IndexListToRemove.size(); i++ )
					if ( m_IndexListToRemove[i]==id ) return true;
				return false;
			}

			// Constructor/Destructor
			RemoveGradientImage_Base():m_BMax(-1)
			{
				ClearListOfIndexes();

			}

			~RemoveGradientImage_Base() {}

			// Function to implement 
			virtual int Execute() = 0;

		protected:
			std::string		m_InputFileName;
			std::string		m_OutputFileName;
			std::vector<int> m_IndexListToRemove;
			double			m_BMax;
		};

		/**********************************************************************************************//**
		 * \class	RemoveGradientImage
		 *
		 * \brief	Main class
		 *
		 * \author	Benoit Scherrer
		 * \date	January 2011
		*************************************************************************************************/
		template <typename TPixel, unsigned int VImageDimension  >
		class ITK_EXPORT RemoveGradientImage : public RemoveGradientImage_Base
		{
		CRL_INSTANCIATOR_CLASS(RemoveGradientImage_Base, RemoveGradientImage);
		public:
			typedef crl::DWIStudy<TPixel> StudyType;

			RemoveGradientImage() {}
			~RemoveGradientImage() {}

			virtual int Execute()
			{
				// Use crl::FileName for filenames
				crl::FileName	input_image(m_InputFileName);
				crl::FileName	output_image(m_OutputFileName);

				//---------------------------------------
				// Load the study
				//---------------------------------------
				StudyType inputStudy;
				inputStudy.LoadStudy(m_InputFileName);

				//---------------------------------------
				// Creates the new one!
				//---------------------------------------
				return CreateNewStudy(inputStudy);
			}

			int CreateNewStudy( StudyType& dwi )
			{
				//-------------------------------------------------------
				// Creates the list of new vectors
				//-------------------------------------------------------
				typename StudyType::GradientVectorSetType newGVecs;
				typename StudyType::BValueSetType newBVals;

				for ( unsigned int i=0; i<dwi.Gradients.size() ; i++ )
				{
					if ( ! ( IsIndexInList(i) || ( m_BMax>=0 && dwi.BValues[i]>m_BMax) ))
					{
						newBVals.push_back( dwi.BValues[i] );
						newGVecs.push_back( dwi.Gradients[i] );
					}
				}

				//-------------------------------------------------------
				// Create the new vector image
				//-------------------------------------------------------
				typename StudyType::DWIGradientImageSetType::Pointer newImage = StudyType::DWIGradientImageSetType::New();
				newImage->SetRegions(dwi.DWIData->GetLargestPossibleRegion());
				newImage->SetOrigin(dwi.DWIData->GetOrigin());
				newImage->SetSpacing(dwi.DWIData->GetSpacing());
				newImage->SetDirection(dwi.DWIData->GetDirection());
				newImage->SetNumberOfComponentsPerPixel(newGVecs.size());
				newImage->Allocate();

				//-------------------------------------------------------
				// Now extract and copy if needed all images
				//-------------------------------------------------------
				std::cout<<"- Copy..."<<std::endl;
				int newId=0;
				for ( unsigned int g=0; g<dwi.Gradients.size() ; g++ )
				{
					if ( IsIndexInList(g) || ( m_BMax>=0 && dwi.BValues[g]>m_BMax) )
						std::cout << "   + Skip image index "<<g<<std::endl;
					else
					{
						typename StudyType::DWIGradientImageType::Pointer oneImage = crl::ExtractComponentOfVectorImage< TPixel, 3>( dwi.DWIData, g );
						SetComponentOfVectorImage<TPixel, VImageDimension>(newImage, newId, oneImage);
						newId++;
					}
				}

				//-------------------------------------------------------
				// Create all the new RAW NHDR file
				//-------------------------------------------------------
				std::cout << "- Create new RAW DWI image in <"<<m_OutputFileName.c_str()<<">..."<<std::endl;
				dwi.Gradients = newGVecs;
				dwi.BValues = newBVals;
				dwi.WriteStudyWithNewData(m_OutputFileName, newImage);

				return 0;
			}
		};




	} // end namespace DWI
} // end namespace crl.

#endif
