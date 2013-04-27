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


#ifndef CRL_DWI_STUDY_TXX
#define CRL_DWI_STUDY_TXX

#include "crlDWIStudy.h"
#include "crlNHDRFileUtils.h"
#include "crlRobustMeanImageFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkOrientImageFilter.h>


namespace crl{

/**********************************************************************************************//**
 * \fn	template <class TPixel> DWIStudy<TPixel>::DWIStudy()
 *
 * \brief	Constructor
 *
 * \author	Benoit Scherrer
 * \date	September 2010
 *
 * \typeparam	TPixel	Pixel type. 
*************************************************************************************************/
template <class TPixel>
DWIStudy<TPixel>::DWIStudy():
m_GradientCoordinatesMode(GRADIENTS_SCANNER_COORDINATES)
{
	NominalBValue = -1;
}

/**********************************************************************************************//**
 * \fn	DWIStudy::~DWIStudy()
 *
 * \brief	Destructor. 
 *
 * \author	Benoit Scherrer
 * \date	February 2011
*************************************************************************************************/
template <class TPixel>
DWIStudy<TPixel>::~DWIStudy()
{
	Gradients.clear();
	BValues.clear();
	OriginalGradients.clear();
}
		
template <class TPixel>
DWIStudy<TPixel>::DWIStudy(const DWIStudy<TPixel>&ref)
{
	this->DWIData						= ref.DWIData;
	this->Gradients						= ref.Gradients;
	this->BValues						= ref.BValues;
	this->NominalBValue					= ref.NominalBValue;
	this->OriginalGradients				= ref.OriginalGradients;	
	this->m_GradientCoordinatesMode		= ref.m_GradientCoordinatesMode;
}
		
template <class TPixel>
void DWIStudy<TPixel>::operator=(const DWIStudy<TPixel>& ref)
{
	this->DWIData						= ref.DWIData;
	this->Gradients						= ref.Gradients;
	this->BValues						= ref.BValues;
	this->NominalBValue					= ref.NominalBValue;
	this->OriginalGradients				= ref.OriginalGradients;	
	this->m_GradientCoordinatesMode		= ref.m_GradientCoordinatesMode;

}


/**********************************************************************************************//**
 * \fn	template <class TPixel> void DWIStudy<TPixel>::LoadStudy ( const std::string& fileName )
 *
 * \brief	Loads a DWI study. 
 *
 * \author	Benoit Scherrer
 * \date	September 2010
 *
 * \typeparam	TPixel	Pixel type.  
 * \param	fileName	Filename of the file. 
*************************************************************************************************/
template <class TPixel>
void DWIStudy<TPixel>::LoadStudy ( const std::string& fileName )
{
	crl::FileName file(fileName);
	std::string fileExt=file.getExtension();

	if ( fileExt=="nhdr" )
		LoadNHDRStudy(fileName);
	else if ( fileExt=="nii" || fileExt=="nii.gz" )
		LoadFSLStudy(fileName);
	else {
		std::string errMsg = "Error while opening "+fileName+".\nDon't know how to open a '"+fileExt+"' DWI study"; 
		throw itk::ExceptionObject(__FILE__,__LINE__,errMsg, "");
	}
}

template <class TPixel>
void DWIStudy<TPixel>::LoadFSLStudy ( const std::string& fileName )
{
	throw itk::ExceptionObject(__FILE__,__LINE__,"LoadFSLStudy: To DO", "");
}

template <class TPixel>
void DWIStudy<TPixel>::LoadNHDRStudy ( const std::string& fileName )
{
	//---------------------------------------
	// First load the DWI vector image
	//---------------------------------------
	std::cout<<"- Loading the DWI file <" << fileName << ">..." << std::endl;
	typedef itk::ImageFileReader<DWIGradientImageSetType> DWIStudyReaderType;
	typename DWIStudyReaderType::Pointer reader = DWIStudyReaderType::New();
	reader->SetFileName(fileName);
	reader->Update();
	DWIData = reader->GetOutput();

	//--------------------------------------------
	// Prepare to read the B-Values/B-Vectors
	//--------------------------------------------
	NominalBValue = 0;
	Gradients.clear();
	BValues.clear();

	GradientVectorType vect3d;
	bool m_ReadB0 = false;
	int nbNullGrad = 0;

	//--------------------------------------------
	// Read the keys in the meta dictionary
	//--------------------------------------------
	std::vector<std::string> imgMetaKeys = DWIData->GetMetaDataDictionary().GetKeys();
	std::vector<std::string>::const_iterator itKey = imgMetaKeys.begin();
	std::string metaString;
	for (; itKey != imgMetaKeys.end(); itKey ++)
	{
		itk::ExposeMetaData<std::string> (DWIData->GetMetaDataDictionary(), *itKey, metaString);

		//--------------------------------------------
		// Is it a gradient vector?
		//--------------------------------------------
		if (itKey->find("DWMRI_gradient") != std::string::npos)
		{ 
			double x,y,z;
			sscanf(metaString.c_str(), "%lf %lf %lf\n", &x, &y, &z);
			vect3d[0] = x; vect3d[1] = y; vect3d[2] = z;

			if ( vect3d.one_norm()<1e-8 ) nbNullGrad++;
			Gradients.push_back(vect3d);
		}

		//--------------------------------------------
		// Is it the nominal B-value?
		//--------------------------------------------
		else if (itKey->find("DWMRI_b-value") != std::string::npos)
		{
			m_ReadB0 = true;
			NominalBValue = atof(metaString.c_str());
		}
	}

	//--------------------------------------------
	// Check that the nominal b-value has been read
	//--------------------------------------------
	if(!NominalBValue)
		throw itk::ExceptionObject(__FILE__,__LINE__,"Nominal B-Value not found in the NHDR file", "");

	std::cout << "  - Found " << Gradients.size()-nbNullGrad << " gradient vectors, "<<nbNullGrad<<" B=0 and a nominal B-value of "<<NominalBValue<<"."<<std::endl;

	//--------------------------------------------
	// Now normalize the gradients and compute the b-values
	//--------------------------------------------
	for ( unsigned int i=0; i<Gradients.size(); i++ )
	{
		vect3d = Gradients[i];

		if ( vect3d.one_norm() <= 1e-8 )
		{
			vect3d[0]=vect3d[1]=vect3d[2]=0;
			Gradients[i] = vect3d;
			BValues.push_back(0);
		}
		else
		{
			double nn = vect3d.two_norm();
			double b = NominalBValue * nn * nn;
			Gradients[i] = ( vect3d / nn );
			BValues.push_back(b);
		}
	}

	//----------------------------------------
	// Keep the original gradients in case OrientStudyInAxial is called
	//----------------------------------------
	OriginalGradients.clear();
	for ( unsigned int i=0; i<Gradients.size(); i++ )
		OriginalGradients.push_back(Gradients[i]);
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> void DWIStudy<TPixel>::WriteStudy( const std::string& fileName )
 *
 * \brief	Writes a DWI study. 
 *
 * \author	Benoit Scherrer
 * \date	September 2010
 *
 * \typeparam	TPixel	Pixel type.  
 * \param	fileName	Filename of the file. 
*************************************************************************************************/
template <class TPixel>
void DWIStudy<TPixel>::WriteStudy( const std::string& fileName, bool useCompression )
{

	crl::NHDR_WriteFromVectorImage<TPixel, 3> ( DWIData,
		Gradients,
		BValues,
		fileName,
		useCompression
		);
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> void DWIStudy<TPixel>::WriteStudyWithNewData( const std::string& fileName,
 * 		typename DWIGradientImageSetType::Pointer data )
 *
 * \brief	Writes the study with new data but by keeping the same gradient vectors/b-values. 
 *
 * \author	Benoit Scherrer
 * \date	September 2010
 *
 * \typeparam	TPixel	. 
 * \param	fileName	Filename of the file. 
 * \param	data		The data. 
*************************************************************************************************/
template <class TPixel>
void DWIStudy<TPixel>::WriteStudyWithNewData( const std::string& fileName, typename DWIGradientImageSetType::Pointer data, bool useCompression )
{
	crl::NHDR_WriteFromVectorImage<TPixel, 3> ( data,
		Gradients,
		BValues,
		fileName,
		useCompression
		);
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> void DWIStudy<TPixel>::SetGradientCoordinates(GradientCoordinatesType gc)
 *
 * \brief	Sets the gradient coordinates mode. If GRADIENTS_MFRAME_COORDINATES, the
 * 			OrientStudyInAxial function will reorient the gradients to the scanner coordinates. 
 *
 * \author	Benoit Scherrer
 * \date	December 2010
 *
 * \typeparam	TPixel	. 
 * \param	gc	The gradient coordinates type. 
 *
 * \sa	OrientStudyInAxial
*************************************************************************************************/
template <class TPixel>
void 
DWIStudy<TPixel>::SetGradientCoordinates(GradientCoordinatesType gc)
{
	m_GradientCoordinatesMode = gc;
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> typename DWIStudy<TPixel>::GradientCoordinatesType DWIStudy<TPixel>::GetGradientCoordinates()
 *
 * \brief	Gets the current gradient coordinates. mode 
 *
 * \author	Benoit Scherrer
 * \date	December 2010
 *
 * \typeparam	TPixel	. 
 *
 * \return	The gradient coordinates mode. 
*************************************************************************************************/
template <class TPixel>
typename DWIStudy<TPixel>::GradientCoordinatesType
DWIStudy<TPixel>::GetGradientCoordinates()
{
	return m_GradientCoordinatesMode;
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> unsigned int DWIStudy<TPixel>::GetNumberOfB0() const
 *
 * \brief	Gets the number of b=0 images
 *
 * \author	Benoit Scherrer
 * \date	September 2010
 *
 * \typeparam	TPixel	. 
 *
 * \return	The number of b=0 images. 
*************************************************************************************************/
template <class TPixel>
unsigned int DWIStudy<TPixel>::GetNumberOfB0() const
{
	int nb=0;
	for ( unsigned int i=0; i<BValues.size(); i++ )
	{
		if ( BValues[i]==0 ) nb++;
	}
	return nb;
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> unsigned int DWIStudy<TPixel>::GetNumberOfNonNullVectors() const
 *
 * \brief	Gets the number of non null gradient vectors. 
 *
 * \author	Benoit Scherrer
 * \date	September 2010
 *
 * \typeparam	TPixel	. 
 *
 * \return	The number of non null gradient vectors. 
*************************************************************************************************/
template <class TPixel>
unsigned int DWIStudy<TPixel>::GetNumberOfNonNullVectors() const
{
	int nb=0;
	for ( unsigned int i=0; i<BValues.size(); i++ )
	{
		if ( BValues[i]!=0 ) nb++;
	}
	return nb;
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> unsigned int DWIStudy<TPixel>::GetTotalNumberOfImages() const
 *
 * \brief	Gets the total number of images.
 *
 * \author	Benoit Scherrer
 * \date	June 2011
 *
 * \tparam	TPixel	Type of the pixel.
 *
 * \return	The total number of images.
 **************************************************************************************************/
template <class TPixel>
unsigned int DWIStudy<TPixel>::GetTotalNumberOfImages() const
{
	return BValues.size();
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> double DWIStudy<TPixel>::GetMaximumBValue() const
 *
 * \brief	Gets the maximum b value.
 *
 * \author	Benoit Scherrer
 * \date	June 2011
 *
 * \exception	itk::ExceptionObject	Thrown if error
 *
 * \tparam	TPixel	Type of the pixel.
 *
 * \return	The maximum b value.
 **************************************************************************************************/
template <class TPixel>
double DWIStudy<TPixel>::GetMaximumBValue() const
{
	if (BValues.size()<1 )
		throw itk::ExceptionObject(__FILE__,__LINE__,"ERROR. The array of b-values is empty.", "DWIStudy::GetMaximumBValue");

	double max = BValues[0];
	for ( unsigned int i=0; i<BValues.size(); i++ )
		if (BValues[i]>max ) max = BValues[i];

	return max;
}

template <class TPixel>
void DWIStudy<TPixel>::ComputeScaledGradients( GradientVectorSetType& scaledGradients ) const
{
	if (BValues.size()!=Gradients.size())
		throw itk::ExceptionObject(__FILE__,__LINE__,"ERROR. The number of b-values does not match the number of gradients", "");

	double maxB = GetMaximumBValue();

	scaledGradients.clear();
	GradientVectorType v;
	for ( unsigned int i=0; i<Gradients.size(); i++ )
	{
		v = Gradients[i] * sqrt(BValues[i]/maxB);
		scaledGradients.push_back(v);
	}
}


/**********************************************************************************************//**
 * \fn	template <class TPixel> void DWIStudy<TPixel>::OrientStudyInAxial ()
 *
 * \brief	Orients the current study in axial. If the gradient coordinates mode is set to
 * 			GRADIENTS_MFRAME_COORDINATES, the gradient vectors are modified accordingly to the
 * 			axes flip and permutations performed by the orienter. 
 *
 * \author	Benoit Scherrer
 * \date	September 2010
 *
 * \typeparam	TPixel	. 
 *
 * \sa	SetGradientCoordinates
*************************************************************************************************/
template <class TPixel>
void DWIStudy<TPixel>::OrientStudyInAxial ()
{
	//-----------------------------------------
	// Orient the input image
	//-----------------------------------------
	std::cout<<"  - Orient the study to axial..."<<std::endl;
	typedef typename itk::OrientImageFilter<DWIGradientImageSetType, DWIGradientImageSetType> OrientVectorImageFilterType;
	typename OrientVectorImageFilterType::Pointer orienter = OrientVectorImageFilterType::New();
	orienter->SetInput(DWIData);
	orienter->UseImageDirectionOn();
	orienter->SetDesiredCoordinateOrientationToAxial();
	orienter->Update();
	DWIData = orienter->GetOutput();
	orienter->GetOutput()->DisconnectPipeline();
	DWIData->Update();

	//-------------------------------------------------------
	// Apply the permutation and then the flip performed by the orienter to the vectors
	//-------------------------------------------------------
	if ( m_GradientCoordinatesMode==GRADIENTS_MFRAME_COORDINATES )
	{
		std::cout<<"  - Apply the permutations/flips to the gradients..."<<std::endl;
		const typename OrientVectorImageFilterType::PermuteOrderArrayType & permuteOrder = orienter->GetPermuteOrder();
		const typename OrientVectorImageFilterType::FlipAxesArrayType & flippedAxes = orienter->GetFlipAxes();
		GradientVectorType outVector;
		for ( unsigned int v=0; v<Gradients.size(); v++ )
		{
			// Permutations
			for ( unsigned int i=0;i<3 ; i++ )
			{
				GradientVectorType inVector = Gradients[v];
				outVector[i] = inVector[permuteOrder[i]];
			}
			// Flips
			for ( unsigned int i=0;i<3 ; i++ )
			{
				if ( flippedAxes[i] && outVector[i]!=0 ) outVector[i] = -outVector[i];
			}
			Gradients[v] = outVector;
		}
	}
}

template <class TPixel>
typename DWIStudy<TPixel>::DWIGradientImageType::Pointer DWIStudy<TPixel>::ComputeMeanB0()
{
	//-----------------------------------------
	// Extract images
	//-----------------------------------------
	std::cout<<"- Extract the B=0 images..."<<std::endl;
	typedef crl::RobustMeanImageFilter< DWIGradientImageType, DWIGradientImageType >  RobustMeanImageFilterType;
	typename RobustMeanImageFilterType::Pointer meanFilter = RobustMeanImageFilterType::New();

	//-------------------------------------------------------
	// Add all B!=0 images and connect them to the robust mean filter
	//-------------------------------------------------------
	int numImg = 0; 
	for( unsigned int b=0; b<BValues.size(); b++ )
	{
		// Take only zero-norm gradients
		if( BValues[b]==0 ) 
		{
			typename DWIGradientImageType::Pointer oneImage = crl::ExtractComponentOfVectorImage< TPixel, 3 >( DWIData, b );
			meanFilter->SetInput(numImg, oneImage);
			numImg++;
		}
	}

	//-------------------------------------------------------
	// Checking...
	//-------------------------------------------------------
	if ( numImg==0 )
		throw itk::ExceptionObject(__FILE__,__LINE__,"ERROR. No null b-value found", "");

	//-------------------------------------------------------
	// Computes!
	//-------------------------------------------------------
	meanFilter->Update();
	std::cout<<std::endl;
	return meanFilter->GetOutput();
}

template <class TPixel>
typename DWIStudy<TPixel>::DWIGradientImageType::Pointer DWIStudy<TPixel>::ComputeMeanBNonNull()
{
	//-----------------------------------------
	// Extract images
	//-----------------------------------------
	std::cout<<"- Extract the B!=0 images..."<<std::endl;
	typedef crl::RobustMeanImageFilter< DWIGradientImageType, DWIGradientImageType >  RobustMeanImageFilterType;
	typename RobustMeanImageFilterType::Pointer meanFilter = RobustMeanImageFilterType::New();

	//-------------------------------------------------------
	// Add all B!=0 images and connect them to the robust mean filter
	//-------------------------------------------------------
	int numImg = 0; 
	for( unsigned int b=0; b<BValues.size(); b++ )
	{
		// Take only zero-norm gradients
		if( BValues[b]!=0 ) 
		{
			typename DWIGradientImageType::Pointer oneImage = crl::ExtractComponentOfVectorImage< TPixel, 3 >( DWIData, b );
			meanFilter->SetInput(numImg, oneImage);
			numImg++;
		}
	}

	//-------------------------------------------------------
	// Checking...
	//-------------------------------------------------------
	if ( numImg==0 )
		throw itk::ExceptionObject(__FILE__,__LINE__,"ERROR. No non-null b-value found", "");

	//-------------------------------------------------------
	// Computes!
	//-------------------------------------------------------
	meanFilter->Update();
	std::cout<<std::endl;
	return meanFilter->GetOutput();
}

} /* end namespace crl. */

#endif





  




