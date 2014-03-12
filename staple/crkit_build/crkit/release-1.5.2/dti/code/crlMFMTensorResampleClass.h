#ifndef CRL_MFM_TENSORRESAMPLE_H
#define CRL_MFM_TENSORRESAMPLE_H

#ifdef _MSC_VER
  #pragma warning( disable: 4996 )
#endif

#include <iostream>
#include <stdlib.h>

#include <itkLightObject.h>
#include <itkImageIOBase.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkDiffusionTensor3D.h>
#include <itkOrientImageFilter.h>
#include <itkTransformFileReader.h>
#include <itkAffineTransform.h>

#include "crlGeneralizedInstantiatorMacros.h"
#include "crlMFMTensorResampleImageFilter.h"
#include "crlFileName.h"
#include "crlMFMUtils.h"
#include "crlTaskProgress.h"

namespace crl {

/**********************************************************************************************//**
 * \class	crlMFMTensorResample_Base
 *
 * \brief	The base class for the instantiator 
 *
 * \author	Benoit Scherrer
 * \date	June 2010
*************************************************************************************************/
class ITK_EXPORT crlMFMTensorResample_Base : public itk::LightObject
{
	CRL_INSTANCIATOR_BASE(crlMFMTensorResample_Base);

public:
	//-----------------------------------
	// Some typedefs
	//-----------------------------------
	typedef enum {
		  INTERP_NONE = 0,
		  INTERP_EUCLIDEAN = 1,
		  INTERP_LOGEUCLIDEAN = 2
	} InterpolationModeValue;
	typedef unsigned int InterpolationModeType;

	//-----------------------------------
	// Constructor/Destructor
	//-----------------------------------
	crlMFMTensorResample_Base() {}
	~crlMFMTensorResample_Base() {}

	//-----------------------------------
	// Algo parameters
	//-----------------------------------
	void SetInputFileNames (const std::vector<std::string>& fileNames ) { m_InputFileNames=fileNames; }
	void SetInputFractionsFileName (const std::string& fileName ) { m_InputFractionsFileName=fileName; }
	void SetOutputBaseFileName (const std::string& fileName ) { m_OutputBaseFileName=fileName; }
	void SetGeometryFixedImage (const std::string& fileName ) { m_GeometryFixedImageFileName=fileName; }
	void SetVoxelSizeMode (const std::string& voxelSize ) { m_VoxelSize=voxelSize; }
	void SetTransformFile (const std::string& fileName ) { m_TransformFileName=fileName; }
	void SetTensorInterpolationMode(InterpolationModeType m) { m_TensorInterpolationMode = m; }

	//-----------------------------------
	// Function to implement 
	//-----------------------------------
	virtual int Execute() = 0;

protected:
	std::vector<std::string>		m_InputFileNames;
	std::string						m_InputFractionsFileName;
	std::string						m_OutputBaseFileName;
	std::string						m_GeometryFixedImageFileName;
	std::string						m_TransformFileName;
	std::string						m_VoxelSize;
	InterpolationModeType			m_TensorInterpolationMode; /// Tensor interpolation mode

};

/**********************************************************************************************//**
 * \class	crlMFMTensorResample
 *
 * \brief	The actual class for the instantiator
 *
 * \author	Benoit Scherrer
 * \date	June 2010
*************************************************************************************************/
template <typename TPixel, unsigned int VImageDimension  >
class ITK_EXPORT crlMFMTensorResample : public crlMFMTensorResample_Base
{
	CRL_INSTANCIATOR_CLASS(crlMFMTensorResample_Base, crlMFMTensorResample);

public:
	//-----------------------------------
	// Some typedefs
	//-----------------------------------
	typedef typename itk::DiffusionTensor3D<TPixel>						TensorPixelType;
	typedef typename itk::VectorImage< float, 3 >						FractionsImageType;
	typedef typename itk::Image< TensorPixelType, VImageDimension >		TensorImageType;
	typedef std::vector<typename TensorImageType::Pointer>				TensorImageVectorType;
	typedef MFMTensorResampleImageFilter<TPixel, VImageDimension> MTensorResampleFilterType;

  typedef  itk::AffineTransform<double,VImageDimension> 	TransformType;

	//-----------------------------------
	// Constructor/Destructor
	//-----------------------------------
	crlMFMTensorResample() {}
	~crlMFMTensorResample() {}

	//-----------------------------------
	// Execute function!
	//-----------------------------------
	virtual int Execute()
	{
		if ( (m_GeometryFixedImageFileName=="" && m_VoxelSize=="")
			|| (m_GeometryFixedImageFileName!="" && m_VoxelSize!="")  )
			throw itk::ExceptionObject(__FILE__, __LINE__, "ERROR. A geometry fixed image OR an output voxel size must be set", "");

		//-------------------------------------
		// Read the voxel size (and check if error before loading files)
		//-------------------------------------
		float sX, sY, sZ;
		sX = sY = sZ = 0;
		if ( m_VoxelSize!="" )
		{
			if ( sscanf(m_VoxelSize.c_str(),"%f,%f,%f", &sX, &sY, &sZ)!=3 )
				throw itk::ExceptionObject(__FILE__, __LINE__, "ERROR. Invalid voxel size. It should be given as sX,sY,sZ (comma-separated values)." , "");

			if ( sX<=0 || sY<=0 || sZ<=0 )
			{
				std::stringstream err;
				err << "ERROR. The specified spacing is : ("<<sX<<","<<sY<<","<<sZ<<") : a spacing component has to be strictely positive.";
				throw itk::ExceptionObject(__FILE__, __LINE__, err.str(), "");
			}
		}

		//-------------------------------------
		//  Creates the multi tensor resample filter
		//  and set properties
		//-------------------------------------
		typename MTensorResampleFilterType::Pointer resampleFilter = MTensorResampleFilterType::New();

		//-------------------------------------
		//  Read the input tensor files
		//-------------------------------------
		std::cout << "- Loading input moving tensor images..."<<std::endl;
		typedef typename itk::ImageFileReader<TensorImageType> TensorImageReaderType;
		typename TensorImageReaderType::Pointer reader = TensorImageReaderType::New();
		for ( unsigned int i=0; i<m_InputFileNames.size() ; i++ )
		{
			reader->SetFileName( m_InputFileNames[i] );
			reader->Update();
			resampleFilter->SetInput(i, reader->GetOutput());
			reader->GetOutput()->DisconnectPipeline();
		}
  
		//-------------------------------------
		//  If defined, open the fractions image
		//-------------------------------------
		if ( m_InputFractionsFileName!="" )
		{
			std::cout << "- Loading input fractions file..."<<std::endl;
			FractionsImageType::Pointer fractionImage = MFM_OpenFractionsImageAsVectorImage<float>(m_InputFractionsFileName);

			if ( fractionImage->GetNumberOfComponentsPerPixel()<m_InputFileNames.size() )
			{
				std::cout << "ERROR. The fraction file has less component than the number of tensors"<<std::endl;
				std::cout << "IGNORE THE FRACTIONS"<<std::endl;
			}
			else
				resampleFilter->SetFractionsImage(fractionImage);
		}

		//-------------------------------------
		// Open the geometry fixed image 
		//-------------------------------------
		typedef typename itk::Image<unsigned short, VImageDimension> GeometryImageType;
		typename GeometryImageType::Pointer geometryImage;
		if ( m_GeometryFixedImageFileName!="" )
		{
			std::cout << "- Loading input fixed geometry image..."<<std::endl;
			typedef typename itk::ImageFileReader<GeometryImageType> GeometryImageReaderType;
			typename GeometryImageReaderType::Pointer greader = GeometryImageReaderType::New();
			greader->SetFileName( m_GeometryFixedImageFileName );
			greader->Update();

			typename itk::OrientImageFilter<GeometryImageType,GeometryImageType>::Pointer orienter = itk::OrientImageFilter<GeometryImageType,GeometryImageType>::New();
			orienter->UseImageDirectionOn();
			orienter->SetDesiredCoordinateOrientationToAxial();
			orienter->SetInput( greader->GetOutput() );
			orienter->Update();

			geometryImage = orienter->GetOutput();

			/*
			FIX FIX FIX
			check for tensor data the implications of changing the space and
			space directions on the measurement frame
			*/
			/* TO FIX: 
			Compare the Direction Cosines of the m_InputMovingImage to that of the
			Direction Cosines of the m_InputFixedImage.
			If these are different, then the resampling operation will modify the
			direction cosines, and that will impact the measurement frame of the
			tensors.  
			Any difference in ordering of the axes or in direction of an axes, must
			be considered for the potential impact upon the measurement frame of
			the tensors.
			*/

			if (resampleFilter->GetInput()->GetDirection() != geometryImage->GetDirection() ) {
				std::cerr << "  Warning : changing the orientation of the input to match" <<
					" that of the existing image. May require a change in the " <<
					" measurement frame." << std::endl;
			}
		}

		//-------------------------------------
		// Read the transform file (or identity if no file)
		//-------------------------------------
		typename TransformType::Pointer transform = ReadTransformFromFile();

		//-------------------------------------
		// Set the resampler parameters
		//-------------------------------------
		resampleFilter->SetTensorInterpolationMode((typename MTensorResampleFilterType::InterpolationModeType)m_TensorInterpolationMode);
		resampleFilter->SetTransform( transform );
		
		if ( geometryImage.GetPointer()!=NULL )
		{
			resampleFilter->SetOutputOrigin( geometryImage->GetOrigin() );
			resampleFilter->SetOutputDirection ( geometryImage->GetDirection() );
			resampleFilter->SetOutputSpacing( geometryImage->GetSpacing() );
			resampleFilter->SetSize( geometryImage->GetLargestPossibleRegion().GetSize() );
			resampleFilter->SetOutputStartIndex( geometryImage->GetLargestPossibleRegion().GetIndex() );
		}
		else
		{
			resampleFilter->SetOutputOrigin( resampleFilter->GetInput(0)->GetOrigin() );
			resampleFilter->SetOutputDirection ( resampleFilter->GetInput(0)->GetDirection() );
			resampleFilter->SetOutputStartIndex( resampleFilter->GetInput(0)->GetLargestPossibleRegion().GetIndex() );

			typename GeometryImageType::SpacingType outputSpacing;
			typename GeometryImageType::SizeType outputSize = resampleFilter->GetInput(0)->GetLargestPossibleRegion().GetSize();
			outputSpacing[0]=sX; 
			outputSpacing[1]=sY; 
			outputSpacing[2]=sZ; //!!tmp not compatible with dim!=3
			for (unsigned int i = 0; i < 3; i++) 
			{
				outputSize[i] = static_cast<unsigned long int>(
					(resampleFilter->GetInput(0)->GetSpacing()[i] *	outputSize[i]) / outputSpacing[i] );
			}
			resampleFilter->SetOutputSpacing( outputSpacing );
			resampleFilter->SetSize( outputSize );
		}

		//-------------------------------------
		// Set the resamping filter properties:
		// Size, spacing, origin and direction will be taken from the fixed geometry image
		//-------------------------------------
		std::cout << "- Resample..."<<std::endl;

		//resampleFilter->AddObserver(itk::ProgressEvent(), callbackProgress);
		//ItkFilterProgressObserver::Pointer progressObs = ItkFilterProgressObserver::New();
		//progressObs->SetProgressLinePrefix("  Processing... ");
		//resampleFilter->AddObserver(ProgressEvent(), progressObs );
		//resampleFilter->AddObserver(StartEvent(), progressObs );
		resampleFilter->SetShowProgress(true);

		//-------------------------------------
		//  Compute
		//-------------------------------------
		resampleFilter->Update();

		//-------------------------------------
		// Save the output!
		//-------------------------------------
		std::cout << "- Writing the resulting tensor images..."<<std::endl;

		typedef typename itk::ImageFileWriter<TensorImageType> TensorImageWriterType;
		typename TensorImageWriterType::Pointer writer = TensorImageWriterType::New();

		if ( resampleFilter->GetNumberOfOutputs()==1 )
		{
			writer->SetFileName( m_OutputBaseFileName );
			writer->SetInput(resampleFilter->GetOutput());
			writer->Update();
		}
		else
		{
			for ( unsigned int i=0; i<resampleFilter->GetNumberOfOutputs() ; i++ )
			{
				char szBuffer[128];
				sprintf(szBuffer, "_t%d", i);

				crl::FileName fn(m_OutputBaseFileName);
				writer->SetFileName( fn.getCompleteFilePath_WithSufix(std::string(szBuffer)) );
				writer->SetInput(resampleFilter->GetOutput(i));
				writer->Update();
			}
		}

		//-------------------------------------
		// If necessary, output the resampled fractions
		//-------------------------------------
		FractionsImageType::Pointer outputFractions = resampleFilter->GetResampledFractionsImage();
		if ( outputFractions.GetPointer()!=NULL )
		{
			typedef typename itk::ImageFileWriter<FractionsImageType> FractionsImageWriterType;
			typename FractionsImageWriterType::Pointer writer = FractionsImageWriterType::New();

			crl::FileName fn(m_OutputBaseFileName);
			writer->SetFileName( fn.getCompleteFilePath_WithSufix("_fractions") );
			writer->SetInput(outputFractions);
			writer->Update();
		}

		return 0;
	}

	typename TransformType::Pointer ReadTransformFromFile()
	{

		typename TransformType::Pointer outputTransform = TransformType::New();
		outputTransform->SetIdentity();

		// Returns the identity if no file
		if ( m_TransformFileName=="" )
			return outputTransform;

		std::cout << "- Loading input transform file..."<<std::endl;

		typedef typename TransformType::ParametersType TransformParametersType;
		TransformParametersType parms = outputTransform->GetParameters();

		//---------------------------------------------
		// Let's see if we can read it with an ITK transform reader
		//---------------------------------------------
		bool readThroughITKTransformReader = true;
		itk::TransformFileReader::Pointer trsfreader;
		trsfreader = itk::TransformFileReader::New();
		trsfreader->SetFileName( m_TransformFileName );

		try {
			trsfreader->Update();
		} catch ( itk::ExceptionObject & excp )
		{
			(void)excp;
			readThroughITKTransformReader = false;
		}

		//---------------------------------------------
		// If it was successfull
		//---------------------------------------------
		if (readThroughITKTransformReader) 
		{
			typedef itk::TransformFileReader::TransformListType * TransformListType;
			TransformListType transforms = trsfreader->GetTransformList();
			itk::TransformFileReader::TransformListType::const_iterator it = transforms->begin();
			if (transforms->size() <= 0 || transforms->size() > 1) 
			{
				std::cerr << "Read " << transforms->size() << 
					" transforms but want 1." << std::endl;
				exit(1);
			}
			if (strcmp((*it)->GetNameOfClass(), "AffineTransform")) {
				std::cerr << "The transform must be an AffineTransform." << std::endl;
				exit(1);
			}
			typename TransformType::Pointer affine_read = 
				static_cast<TransformType*>((*it).GetPointer());
			outputTransform->SetTranslation( affine_read->GetTranslation() );
			outputTransform->SetCenter( affine_read->GetCenter() );
			outputTransform->SetMatrix( affine_read->GetMatrix() );

			return outputTransform; // The transform has been successfully read.
		}

		// Continue here if the ITK transform reader didn't work out.


		int expectedcount = parms.Size();

		int count = 0;
		std::ifstream parmfile(m_TransformFileName.c_str());
		std::string line;
		std::string delims(" \t:");
		std::string::size_type idx;
		while (getline(parmfile, line)) {
			// Skip lines that are blank or start with a #
			if (line.length() == 0) continue;
			if (line[0] == '#') continue;
			idx = line.find("parameters");
			if (idx == std::string::npos) {
				std::cerr << "Failed finding parameters line in " << m_TransformFileName.c_str() << std::endl;
				std::cerr << "line length is " << line.length() << std::endl;
				std::cerr << "line is " << line << std::endl;
				exit(2);
			}
			// Skip the colon and leading white space in the answer.
			idx = line.rfind(":");
			idx = line.find_first_not_of(delims, idx);
			if (idx == std::string::npos) {
				std::cerr << "Failed parsing transform - missing tag value."
					<< std::endl;
				exit(2);
			}
			if (expectedcount == 7) {
				count = sscanf(line.substr(idx).c_str(),
					"%lg %lg %lg %lg %lg %lg %lg",
					&parms[0], &parms[1], &parms[2],
					&parms[3], &parms[4], &parms[5],
					&parms[6]
				);
			} else if (expectedcount == 12) {
				count = sscanf(line.substr(idx).c_str(),
					"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
					&parms[0], &parms[1], &parms[2],
					&parms[3], &parms[4], &parms[5],
					&parms[6], &parms[7], &parms[8],
					&parms[9], &parms[10], &parms[11]
				);
			} else {
				std::cerr << "Unexpected number of transform parameters." <<
					std::endl;
			}

			if (count != expectedcount) {
				std::cerr << "Read " << count << " instead of " <<
					expectedcount << " doubles." <<
					std::endl;
			}
		}

		if (count == expectedcount) {
			outputTransform->SetParameters(parms);
		} else {
			std::cerr << "Failed to find transform parameters." << std::endl;
			exit(3);
		}

		return outputTransform;
	}

	

};

} //end of namespace crl

#endif
