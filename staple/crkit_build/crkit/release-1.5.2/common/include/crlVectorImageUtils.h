/**********************************************************************************************//**
 * \file	crlVectorImageUtils.h
 *
 * \brief	Declare useful functions for vector image manipulations
*************************************************************************************************/

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


#ifndef CRL_VECTOR_IMAGE_UTILS
#define CRL_VECTOR_IMAGE_UTILS

#ifdef _MSC_VER
  #pragma warning( disable: 4996 )
#endif

#include <string>

#include <itkVectorImage.h>
#include <itkImage.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageToVectorImageFilter.h>

#include "crlFileName.h"


namespace crl {

	/**********************************************************************************************//**
	 * \brief	Extracts a component of a vector image in a itk::Image. 
	 *
	 * \author	Benoit Scherrer
	 * \date	December 2009
	 *
	 * \exception	Thrown when error. 
	 *
	 * \param	inputImage	The input vector image. 
	 * \param	id			The image index to extract. 
	 *
	 * \return	the extracted itk::Image.
	*************************************************************************************************/
	
	template  < typename TPixel, unsigned int VImageDimension  >
	typename itk::Image< TPixel, VImageDimension >::Pointer ExtractComponentOfVectorImage( const typename itk::VectorImage<TPixel, VImageDimension>::Pointer &inputImage ,  int id )
	{
		typedef typename itk::VectorImage<TPixel, VImageDimension>	InputVectorImage;
		typedef typename itk::Image<TPixel, VImageDimension>		OutputOneImage;

		typedef typename itk::ImageRegionConstIterator<InputVectorImage>	VecImageInIteratorType;
		typedef typename itk::ImageRegionIterator<OutputOneImage>			OneImageOutIteratorType;

		if ( (id<0 ) || (id >= (int)inputImage->GetNumberOfComponentsPerPixel ()) )
			throw itk::ExceptionObject(__FILE__,__LINE__,"ExtractComponentOfVectorImage: the image id to extract is invalid", "");

		/*-------------------------------------------------
		Create the output image
		/-------------------------------------------------*/
		typename OutputOneImage::Pointer oneImage = OutputOneImage::New();
		oneImage->CopyInformation(inputImage);
		oneImage->SetRegions(inputImage->GetLargestPossibleRegion().GetSize());
		oneImage->SetOrigin(inputImage->GetOrigin());
		oneImage->SetSpacing(inputImage->GetSpacing());
		oneImage->SetDirection(inputImage->GetDirection());
		oneImage->Allocate();

		/*-------------------------------------------------
		Copy one image of the vector image to oneImage
		/-------------------------------------------------*/
		VecImageInIteratorType inItr(inputImage, inputImage->GetLargestPossibleRegion());
		OneImageOutIteratorType outItr(oneImage, oneImage->GetLargestPossibleRegion() );
		for ( inItr.GoToBegin(), outItr.GoToBegin() ; !outItr.IsAtEnd() ; ++inItr, ++outItr )
			outItr.Set(inItr.Get()[id]);

		return oneImage;
	}

	/**********************************************************************************************//**
	 * \brief	Sets the component of a vector image from a source scalar image. 
	 *
	 * \author	Benoit Scherrer
	 * \date	April 2010
	 *
	 * \exception	itk::ExceptionObject	Thrown when error. 
	 *
	 * \param [in,out]	destImage	Destination vector image. 
	 * \param	id					The component index. 
	 * \param	srcImage			Source scalar image. 
	 *
	*************************************************************************************************/
	template  < typename TPixel, unsigned int VImageDimension  >
	void SetComponentOfVectorImage( typename itk::VectorImage<TPixel, VImageDimension>::Pointer &destImage ,  unsigned int id, const typename itk::Image<TPixel, VImageDimension>::Pointer& srcImage )
	{
		//-------------------------------------------------------
		// Defines some useful types
		//-------------------------------------------------------
		typedef typename itk::VectorImage<TPixel, VImageDimension>	VectorImageType;
		typedef typename itk::Image<TPixel, VImageDimension>		OneImageType;

		typedef typename itk::ImageRegionIterator<VectorImageType>	VecImageOutIteratorType;
		typedef typename itk::ImageRegionConstIterator<OneImageType> OneImageInIteratorType;

		//-------------------------------------------------------
		// Check....
		//-------------------------------------------------------
		if ( id>=destImage->GetNumberOfComponentsPerPixel() )
			throw itk::ExceptionObject("SetComponentOfVectorImage: The image id to set is invalid.");

		//-------------------------------------------------------
		// Creates the iterators
		//-------------------------------------------------------
		VecImageOutIteratorType outItr(destImage, srcImage->GetLargestPossibleRegion() );
		OneImageInIteratorType inItr ( srcImage, srcImage->GetLargestPossibleRegion()); 
		inItr.GoToBegin(); outItr.GoToBegin();

		typename VectorImageType::PixelType inP(destImage->GetNumberOfComponentsPerPixel() ) ;

		//-------------------------------------------------------
		// Copy all pixels of the src scalar image
		//-------------------------------------------------------
		for ( outItr.GoToBegin() ; !outItr.IsAtEnd() ; ++outItr, ++inItr )
		{
			inP = outItr.Get();
			inP[id] = inItr.Get();
			outItr.Set(inP);
		}
	}

	/**********************************************************************************************//**
	 * \brief	Extracts and writes out one component of a vector image. 
	 *
	 * \author	Benoit Scherrer
	 * \date	December 2009
	 *
	 * \param	inputImage		The input vector image. 
	 * \param	id				The image index to extract. 
	 * \param	ouputFileName	Filename of the ouput file. 
	*************************************************************************************************/
	template  <typename TPixel, unsigned int VImageDimension  >
	void WriteComponentOfVectorImage( const typename  itk::VectorImage<TPixel, VImageDimension>::Pointer &inputImage ,  int id, const std::string& ouputFileName, bool useCompression=false )
	{
		typename itk::Image<TPixel, VImageDimension>::Pointer oneImage = ExtractComponentOfVectorImage<TPixel, VImageDimension>( inputImage, id );

		/*-------------------------------------------------
		Create the writer
		/-------------------------------------------------*/
		typedef typename itk::ImageFileWriter< itk::Image<TPixel, VImageDimension> > WriterType;
		typename WriterType::Pointer writer = WriterType::New();

		/*-------------------------------------------------
		Save the image
		/-------------------------------------------------*/
		writer->SetInput(oneImage);
		writer->SetFileName( ouputFileName );
		writer->SetUseCompression(useCompression);
		writer->Update();

	}

	/**********************************************************************************************//**
     * \brief	Writes out all components of a vector image to a set of nrrd files. 
	 *
	 * \author	Benoit Scherrer
	 * \date	December 2009
	 *
	 * \param	inputImage					The input image. 
	 * \param	outputFileBaseName			Base name of the output files (created file names will be
	 * 										[outputFileBaseName]_0000, [outputFileBaseName]_0001, ...). 
	 * \param [in,out]	createdFileNames	list of names of the created files. 
	*************************************************************************************************/
	template  <typename TPixel, unsigned int VImageDimension  >
	void WriteAllComponentsOfVectorImage( const typename itk::VectorImage<TPixel, VImageDimension>::Pointer &inputImage,  const std::string& outputFileBaseName, std::vector<std::string> &createdFileNames, bool useCompression=false)
	{
		typedef typename itk::VectorImage<TPixel, VImageDimension>	InputVectorImage;
		typedef typename itk::Image<TPixel, VImageDimension>		OutputOneImage;

		/*-------------------------------------------------
		Create the writer
		/-------------------------------------------------*/
		typedef typename itk::ImageFileWriter<OutputOneImage> WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetUseCompression(useCompression);

		crl::FileName fileBase(outputFileBaseName);

		// If the output type was a nhdr, change it to nrrd
		if ( fileBase.getExtension()=="nhdr" ) 
			fileBase.setExtension("nrrd");

		/*-------------------------------------------------
		For each vector component...
		/-------------------------------------------------*/
		for ( unsigned int i=0; i<inputImage->GetNumberOfComponentsPerPixel () ; i++ )
		{
			typename OutputOneImage::Pointer oneImage = ExtractComponentOfVectorImage<TPixel, VImageDimension> ( inputImage, i );

			/*-------------------------------------------------
			Save the image
			/-------------------------------------------------*/
			char szBuffer[4096];
			sprintf(szBuffer,"_%04d", i);
			writer->SetInput(oneImage);
			std::string fn = fileBase.getCompleteFilePath_WithSufix(std::string(szBuffer));
			writer->SetFileName( fn );
			writer->Update();

			createdFileNames.push_back( fn );
		}
	}

	/**********************************************************************************************//**
	 * \brief	Export all components of a vector image (given by its file name) to a set of nrrd
	 * 			files. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \param	inputFileName				Filename of the input file. 
	 * \param	outputFileBaseName			Base filename of the output file. 
	 * \param [in,out]	createdFileNames	List of create filenames. 
	*************************************************************************************************/
	template  <typename TPixel, unsigned int VImageDimension  >
	void WriteAllComponentsOfVectorImage( const std::string& inputFileName,  const std::string& outputFileBaseName, std::vector<std::string> &createdFileNames, bool useCompression=false)
	{
		typedef typename itk::VectorImage<TPixel, VImageDimension>	InputVectorImage;

		/*-------------------------------------------------
		Open the file
		/-------------------------------------------------*/
		typedef typename itk::ImageFileReader<InputVectorImage> ReaderType;
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( inputFileName );
		reader->Update();
		
		WriteAllComponentsOfVectorImage<TPixel, VImageDimension>(reader->GetOutput(), outputFileBaseName, createdFileNames, useCompression);
	}

	/**********************************************************************************************//**
	 * \brief	Loads n images (of same type) as a single vector image. 
	 *
	 * \author	Benoit Scherrer
	 * \date	August 2010
	 *
	 * \param	fileNames	List of file names.
	 *
	 * \return	The n images as one vector image. 
	*************************************************************************************************/
	template  <typename TPixel, unsigned int VImageDimension  >
	typename itk::VectorImage<TPixel, VImageDimension>::Pointer LoadImageSetAsVectorImage ( const std::vector<std::string> & fileNames  )
	{
		typedef itk::Image<TPixel, VImageDimension> OneImageType;
		typedef itk::VectorImage<TPixel, VImageDimension> VectorImageType;
		
		//-------------------------------------------------------
		// Construct the filter to gather all the single images
		//-------------------------------------------------------
		typedef itk::ImageToVectorImageFilter< OneImageType > ImageToVectorImageFilterType;
		typename ImageToVectorImageFilterType::Pointer filter = ImageToVectorImageFilterType::New();

		//-------------------------------------------------------
		// Read each separate image and connect them to the filter
		//-------------------------------------------------------
		typedef itk::ImageFileReader<OneImageType> OneImageReaderType;
		typename OneImageReaderType::Pointer reader = OneImageReaderType::New();

		for ( unsigned int i=0; i<fileNames.size() ; i++ )
		{
			reader->SetFileName( fileNames[i] );
			reader->Update();
			filter->SetInput(i, reader->GetOutput());
			reader->GetOutput()->DisconnectPipeline();
		}

		//-------------------------------------------------------
		// Execute the filter!	
		//-------------------------------------------------------
		filter->Update();

		return filter->GetOutput();
	}

	/**********************************************************************************************//**
	 *
	 * \brief	Converts a (N-1)-D itk::VectorImage to a N-D itk::Image. For example a 3D vector
	 * 			image can be converted to a 4D FSL scalar image. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \param	vectorImage	The input (N-1)-D vector image. 
	 *
	 * \return	The output N-D image. 
	*************************************************************************************************/
	template  <typename TPixel, unsigned int VInputImageDimension >
	typename itk::Image< TPixel, (VInputImageDimension+1) >::Pointer
	ConvertVectorImageToNDImage( const typename itk::VectorImage<TPixel, VInputImageDimension>::Pointer &vectorImage )
	{
		// Define some types
		typedef itk::VectorImage<TPixel, VInputImageDimension>	InputImageType;
		typedef itk::Image< TPixel, (VInputImageDimension+1) > 	OutputImageType;
		typedef itk::ImageRegionIterator<InputImageType>		InIteratorType;
		typedef itk::ImageRegionIterator<OutputImageType>		OutIteratorType;
	
		// Define some variables
		typename OutputImageType::SizeType 		outSize;
		typename OutputImageType::PointType 	outOrigin;
		typename OutputImageType::SpacingType 	outSpacing;
		typename OutputImageType::DirectionType outDir;

		typename InputImageType::SizeType 		inRegion = vectorImage->GetLargestPossibleRegion().GetSize();
		typename OutputImageType::Pointer 		outImage = OutputImageType::New();

		//-------------------------------------------------
		// Prepare the output size/origin/spacing/direction
		//-------------------------------------------------
		outSpacing.Fill(1); outDir.Fill(0); outOrigin.Fill(0);
		for ( unsigned int i=0 ; i<VInputImageDimension ; i++ )
		{
			outSize[i]    = inRegion[i];
			outSpacing[i] = vectorImage->GetSpacing()[i];
			outOrigin[i]  = vectorImage->GetOrigin()[i];
		}
		outSize.SetElement(VInputImageDimension, vectorImage->GetNumberOfComponentsPerPixel() );

		for ( unsigned int i=0; i<VInputImageDimension ; i++ ) {
			for ( unsigned int j=0; j<VInputImageDimension ; j++ )
			{
				outDir(i,j) = vectorImage->GetDirection()(i,j);
			}
		}
		outDir(VInputImageDimension,VInputImageDimension)=1.0;

		//-------------------------------------------------
		// Allocate the output image
		//-------------------------------------------------
		outImage->SetRegions(outSize);
		outImage->SetOrigin(outOrigin);
		outImage->SetSpacing(outSpacing);
		outImage->SetDirection(outDir);
		outImage->Allocate();

		//-------------------------------------------------
		// Prepare the output region for the iterator
		//-------------------------------------------------
		typename OutputImageType::RegionType outRegion;
		outRegion.SetSize(outSize); outRegion.SetSize(VInputImageDimension,1);
		for ( unsigned int i=0; i<=VInputImageDimension ; i++ )
			outRegion.SetIndex(i,0); 
	
		//-------------------------------------------------
		// Copy the VectorImage to the N-D image
		//-------------------------------------------------
		for (unsigned int i=0 ; i<outSize[VInputImageDimension] ; i++)
		{ 
			outRegion.SetIndex(VInputImageDimension,i);

			// Copy one image of the vector image to one dimension of the 4D image
			InIteratorType inItr(vectorImage, vectorImage->GetLargestPossibleRegion());
			OutIteratorType outItr(outImage, outRegion );
			for ( inItr.GoToBegin(), outItr.GoToBegin() ; !outItr.IsAtEnd() ; ++inItr, ++outItr )
			{
				outItr.Set(inItr.Get()[i]);
			}
		}

		return outImage;
	}

	/**********************************************************************************************//**
	 * \brief	Converts a (N-1)-D itk::VectorImage to a N-D image and writes the result in a file. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \param	inputVectorImage	The input (N-1)-D vector image. 
	 * \param	outputFileName		Filename of the output file. 
	*************************************************************************************************/
	template  <typename TPixel, unsigned int VInputImageDimension >
	void ConvertVectorImageToNDImage( const typename itk::VectorImage<TPixel, VInputImageDimension>::Pointer & inputVectorImage, const std::string& outputFileName, bool useCompression=false )
	{
		typedef typename itk::Image<TPixel, (VInputImageDimension+1) >	OutputImageType;

		/*-------------------------------------------------
		Convert
		/-------------------------------------------------*/
		typename OutputImageType::Pointer outImage = crl::ConvertVectorImageToNDImage<TPixel, VInputImageDimension> ( inputVectorImage );

		/*-------------------------------------------------
		Save the result
		/-------------------------------------------------*/
		typedef typename itk::ImageFileWriter<OutputImageType> WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( outputFileName );
		writer->SetInput(outImage);
		writer->SetUseCompression(useCompression);
		writer->Update();
	}

	/**********************************************************************************************//**
	 * \brief	Converts a (N-1)-D itk::VectorImage given by its filename to a N-D image and writes
	 * 			the result in a file. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \param	inputFileName	Filename of the input file. 
	 * \param	outputFileName	Filename of the output file. 
	*************************************************************************************************/
	template  <typename TPixel, unsigned int VInputImageDimension >
	void ConvertVectorImageToNDImage( const std::string& inputFileName, const std::string& outputFileName, bool useCompression=false  )
	{
		typedef typename itk::VectorImage<TPixel, VInputImageDimension>	InputVectorImageType;

		/*-------------------------------------------------
		Open the file
		/-------------------------------------------------*/
		typedef typename itk::ImageFileReader<InputVectorImageType> ReaderType;
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( inputFileName );
		reader->Update();
		
		/*-------------------------------------------------
		Convert
		/-------------------------------------------------*/
		ConvertVectorImageToNDImage<TPixel, VInputImageDimension> ( reader->GetOutput(), outputFileName, useCompression ) ;
	}

	/**********************************************************************************************//**
	 *
	 * \brief	Convert a N-D itk::Image to a (N-1)-D itk::VectorImage. For example a 4D FSL scalar
	 * 			image can be converted to a 3D vector image. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \param	ndImage	The nd image. 
	 *
	 * \return	. 
	*************************************************************************************************/
	template  <typename TPixel, unsigned int VInputImageDimension >
	typename itk::VectorImage< TPixel, (VInputImageDimension-1) >::Pointer
	ConvertNDToVectorImage( const typename itk::Image<TPixel, VInputImageDimension>::Pointer &ndImage )
	{
		// Define some types
		typedef itk::VectorImage<TPixel, (VInputImageDimension-1) >	OutputImageType;
		typedef itk::Image< TPixel, VInputImageDimension > 			InputImageType;
		typedef itk::ImageRegionIteratorWithIndex<OutputImageType>	OutIteratorType;
		typedef itk::ImageRegionConstIterator<InputImageType>		InIteratorType;
	
		// Define some variables
		typename OutputImageType::SizeType 		outSize;
		typename OutputImageType::PointType 	outOrigin;
		typename OutputImageType::SpacingType 	outSpacing;
		typename OutputImageType::DirectionType outDir;

		typename InputImageType::SizeType 		inSize = ndImage->GetLargestPossibleRegion().GetSize();
		typename OutputImageType::Pointer 		outImage = OutputImageType::New();

		//-------------------------------------------------
		// Prepare the output size/origin/spacing/direction
		//-------------------------------------------------
		for ( unsigned int i=0 ; i<(VInputImageDimension-1) ; i++ )
		{
			outSize[i]    = inSize[i];
			outSpacing[i] = ndImage->GetSpacing()[i];
			outOrigin[i]  = ndImage->GetOrigin()[i];
		}
		outImage->SetNumberOfComponentsPerPixel( inSize[VInputImageDimension-1] ) ;

		for ( unsigned int i=0; i<(VInputImageDimension-1) ; i++ ) {
			for ( unsigned int j=0; j<(VInputImageDimension-1) ; j++ )
			{
				outDir(i,j) = ndImage->GetDirection()(i,j);
			}
		}

		//-------------------------------------------------
		// Allocate the output image
		//-------------------------------------------------
		outImage->SetRegions(outSize);
		outImage->SetOrigin(outOrigin);
		outImage->SetSpacing(outSpacing);
		outImage->SetDirection(outDir);
		outImage->Allocate();

		//-------------------------------------------------
		// Prepare the output region for the iterator
		//-------------------------------------------------
		typename OutputImageType::RegionType outRegion;
		outRegion.SetSize(outSize);
		for ( unsigned int i=0; i<(VInputImageDimension-1) ; i++ )
			outRegion.SetIndex(i,0); 
	
		//-------------------------------------------------
		// Prepare the input region for the iterator
		//-------------------------------------------------
		typename InputImageType::RegionType inRegion;
		inRegion.SetSize((VInputImageDimension-1), inSize[(VInputImageDimension-1)] );
		for ( unsigned int i=0; i<(VInputImageDimension-1) ; i++ )
			inRegion.SetSize(i, 1);

		//-------------------------------------------------
		// Copy the N-D image to the VectorImage image
		//-------------------------------------------------
		OutIteratorType outItr(outImage, outRegion );
		for ( outItr.GoToBegin() ; !outItr.IsAtEnd() ; ++outItr )
		{
			typename OutputImageType::IndexType idxOutput = outItr.GetIndex ();
			for ( unsigned int i=0; i<(VInputImageDimension-1) ; i++ )
				inRegion.SetIndex(i,idxOutput[i]);

			// Copy the N-th dim of the N-D imagine to one vector 
			typename OutputImageType::PixelType p(inSize[(VInputImageDimension-1)]) ;
			
			InIteratorType inItr(ndImage, inRegion );
			int id=0;
			for ( inItr.GoToBegin() ; !inItr.IsAtEnd() ; ++inItr )
			{
				p[id] = inItr.Value();
				id++;
			}

			outItr.Set(p);
		}

		return outImage;
	}

	/**********************************************************************************************//**
	 * \brief	Convert a N-D itk::Image to a (N-1)-D itk::VectorImage and writes the result in a
	 * 			file. 
	 *
	 * \author	Benoit Scherrer
	 * \date	June 2010
	 *
	 * \param	ndImage			The nd image. 
	 * \param	outputFileName	Filename of the output file. 
	*************************************************************************************************/
	template  <typename TPixel, unsigned int VInputImageDimension >
	void ConvertNDToVectorImage( const typename itk::Image<TPixel, VInputImageDimension>::Pointer &ndImage, const std::string& outputFileName, bool useCompression=false )
	{
		typedef typename itk::VectorImage<TPixel, (VInputImageDimension-1) >	OutputImageType;

		/*-------------------------------------------------
		Convert
		/-------------------------------------------------*/
		typename OutputImageType::Pointer outImage = crl::ConvertNDToVectorImage<TPixel, VInputImageDimension> ( ndImage );

		/*-------------------------------------------------
		Save the result
		/-------------------------------------------------*/
		typedef typename itk::ImageFileWriter<OutputImageType> WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( outputFileName );
		writer->SetInput(outImage);
		writer->SetUseCompression(useCompression);
		writer->Update();
	}

	template  <typename TPixel, unsigned int VInputImageDimension >
	void ConvertNDToVectorImage( const std::string& inputFileName, const std::string& outputFileName, bool useCompression=false )
	{
		typedef typename itk::Image<TPixel, VInputImageDimension >				InputImageType ;

		/*-------------------------------------------------
		Open the file
		/-------------------------------------------------*/
		typedef typename itk::ImageFileReader<InputImageType> ReaderType;
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( inputFileName );
		reader->Update();
		
		/*-------------------------------------------------
		Convert
		/-------------------------------------------------*/
		ConvertNDToVectorImage<TPixel, VInputImageDimension> ( reader->GetOutput(), outputFileName, useCompression ) ;
	}
}

#endif
