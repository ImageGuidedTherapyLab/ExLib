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

#ifndef CRL_DWI_UTILS_H
#define CRL_DWI_UTILS_H

#include "crlFileName.h"
#include "crlVectorImageUtils.h"

#include "itkDiffusionTensor3D.h"
#include "itkVectorImage.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include <vnl/vnl_cross.h>

namespace crl {
namespace DWI {

	template  < typename TPixel, unsigned int VImageDimension  >
	typename itk::Image< typename itk::DiffusionTensor3D< TPixel >, VImageDimension >::Pointer FSLTensorToDiff3DTensor( const std::string& inputBaseFileName, bool verbose=false )
	{
		crl::FileName inputFile(inputBaseFileName);

		typedef typename itk::DiffusionTensor3D< TPixel >					TensorPixelType;
		typedef typename itk::Image< TensorPixelType, VImageDimension > 	TensorImageType;
		typedef typename itk::Image< TPixel, VImageDimension > 				EigenValueImageType;
		typedef typename itk::Image< TPixel, VImageDimension+1 >		 	FSLEigenVectorImageType;
		typedef typename itk::VectorImage< TPixel, VImageDimension >		EigenVectorImageType;

		//---------------------------------------------------
		// Read the eigen values
		//---------------------------------------------------
		if ( verbose ) std::cout<<"- Read eigen value from <"<<inputFile.getCompleteFilePath_WithSufix("_L1" ).c_str()<<">..." << std::endl;
		typedef typename itk::ImageFileReader< EigenValueImageType > EigenValueReaderType;
		typename EigenValueReaderType::Pointer fslL1 = EigenValueReaderType::New();
		fslL1->SetFileName( inputFile.getCompleteFilePath_WithSufix("_L1" ));
		fslL1->Update();

		if ( verbose ) std::cout<<"- Read eigen value from <"<<inputFile.getCompleteFilePath_WithSufix("_L2" ).c_str()<<">..." << std::endl;
		typename EigenValueReaderType::Pointer fslL2 = EigenValueReaderType::New();
		fslL2->SetFileName( inputFile.getCompleteFilePath_WithSufix("_L2" ));
		fslL2->Update();

		if ( verbose ) std::cout<<"- Read eigen value from <"<<inputFile.getCompleteFilePath_WithSufix("_L3" ).c_str()<<">..." << std::endl;
		typename EigenValueReaderType::Pointer fslL3 = EigenValueReaderType::New();
		fslL3->SetFileName( inputFile.getCompleteFilePath_WithSufix("_L3" ));
		fslL3->Update();
		
		//---------------------------------------------------
		// Read the eigen vectors
		//---------------------------------------------------
		if ( verbose ) std::cout<<"- Read eigen vector from <"<<inputFile.getCompleteFilePath_WithSufix("_V1" ).c_str()<<">..." << std::endl;
		typedef typename itk::ImageFileReader< FSLEigenVectorImageType > EigenVectorReaderType;
		typename EigenVectorReaderType::Pointer fslV1 = EigenVectorReaderType::New();
		fslV1->SetFileName( inputFile.getCompleteFilePath_WithSufix("_V1" ));
		fslV1->Update();
		typename EigenVectorImageType::Pointer V1 = ConvertNDToVectorImage<TPixel, VImageDimension+1>(fslV1->GetOutput());

		if ( verbose ) std::cout<<"- Read eigen vector from <"<<inputFile.getCompleteFilePath_WithSufix("_V2" ).c_str()<<">..." << std::endl;
		typename EigenVectorReaderType::Pointer fslV2 = EigenVectorReaderType::New();
		fslV2->SetFileName( inputFile.getCompleteFilePath_WithSufix("_V2" ));
		fslV2->Update();
		typename EigenVectorImageType::Pointer V2 = ConvertNDToVectorImage<TPixel, VImageDimension+1>(fslV2->GetOutput());

		if ( verbose ) std::cout<<"- Read eigen vector from <"<<inputFile.getCompleteFilePath_WithSufix("_V3" ).c_str()<<">..." << std::endl;
		typename EigenVectorReaderType::Pointer fslV3 = EigenVectorReaderType::New();
		fslV3->SetFileName( inputFile.getCompleteFilePath_WithSufix("_V3" ));
		fslV3->Update();
		typename EigenVectorImageType::Pointer V3 = ConvertNDToVectorImage<TPixel, VImageDimension+1>(fslV3->GetOutput());

		//---------------------------------------------------
		// Create the itk::DiffusionTensor3D image
		//---------------------------------------------------
		if ( verbose ) std::cout << "- Convert to a Diffusion Tensor 3D Image..."  << std::endl;
		typename TensorImageType::Pointer outputImage = TensorImageType::New();
		outputImage->SetRegions(fslL1->GetOutput()->GetLargestPossibleRegion());
		outputImage->SetOrigin(fslL1->GetOutput()->GetOrigin());
		outputImage->SetSpacing(fslL1->GetOutput()->GetSpacing());
		outputImage->SetDirection(fslL1->GetOutput()->GetDirection());
		outputImage->Allocate();

		//---------------------------------------------------
		// Create the iterators
		//---------------------------------------------------
		typedef typename itk::ImageRegionConstIterator<EigenValueImageType>		EigenValueIteratorType;
		typedef typename itk::ImageRegionConstIterator<EigenVectorImageType>	EigenVectorIteratorType;
		typedef typename itk::ImageRegionIterator<TensorImageType>				OutputIteratorType;

		EigenValueIteratorType itL1(fslL1->GetOutput(), fslL1->GetOutput()->GetLargestPossibleRegion() );
		EigenValueIteratorType itL2(fslL2->GetOutput(), fslL2->GetOutput()->GetLargestPossibleRegion() );
		EigenValueIteratorType itL3(fslL3->GetOutput(), fslL3->GetOutput()->GetLargestPossibleRegion() );

		EigenVectorIteratorType itV1(V1, V1->GetLargestPossibleRegion() );
		EigenVectorIteratorType itV2(V2, V2->GetLargestPossibleRegion() );
		EigenVectorIteratorType itV3(V3, V3->GetLargestPossibleRegion() );

		OutputIteratorType		itOut(outputImage, outputImage->GetLargestPossibleRegion());

		//---------------------------------------------------
		// And now convert!
		//---------------------------------------------------
		vnl_matrix<TPixel>		matrixOut(3,3);
		vnl_matrix<TPixel>		matrixV (3,3);
		vnl_diag_matrix<TPixel> matrixL(3) ;

		for ( itOut.GoToBegin(), 
			  itL1.GoToBegin(), itL2.GoToBegin(), itL3.GoToBegin(), 
			  itV1.GoToBegin(), itV2.GoToBegin(), itV3.GoToBegin() ; 
			  !itOut.IsAtEnd() ; 
			  ++itOut, ++itL1, ++itL2, ++itL3, ++itV1, ++itV2, ++itV3 )
		{
			TensorPixelType diffTensor;

			// Get all informations
			TPixel L1 = itL1.Value();
			TPixel L2 = itL2.Value();
			TPixel L3 = itL3.Value();
			typename EigenVectorImageType::PixelType vecV1  = itV1.Get();
			typename EigenVectorImageType::PixelType vecV2  = itV2.Get();
			typename EigenVectorImageType::PixelType vecV3  = itV3.Get();

			// Transform into a vnl_matrix and a vnl_diag_matrix
			for ( int i=0; i<3; i++ )
			{
				matrixV(i,0) = vecV1[i];
				matrixV(i,1) = vecV2[i];
				matrixV(i,2) = vecV3[i];
			}
			matrixL[0] = L1;
			matrixL[1] = L2;
			matrixL[2] = L3;

			// Compute the diffusion tensor matrix
			matrixOut = matrixV * matrixL * matrixV.transpose();

			// Transform into a diffusion tensor 
			TPixel *p = diffTensor.GetDataPointer ();
			*p = matrixOut(0,0);
			*(p+1) =  matrixOut(0,1);
			*(p+2) = matrixOut(0,2);
			*(p+3) = matrixOut(1,1);
			*(p+4) = matrixOut(1,2);
			*(p+5) = matrixOut(2,2);

			// And set the result in the output image!
			itOut.Set(diffTensor);
		}

		return outputImage;
	}

	template  < typename TPixel, unsigned int VImageDimension  >
	void FSLTensorToDiff3DTensor( const std::string& inputBaseFileName, const std::string& outputFileName, bool verbose=false )
	{
		typedef typename itk::DiffusionTensor3D< TPixel >					TensorPixelType;
		typedef typename itk::Image< TensorPixelType, VImageDimension > 	TensorImageType;
		typedef typename itk::ImageFileWriter<TensorImageType>				TensorWriterType;

		typename TensorImageType::Pointer outImage = FSLTensorToDiff3DTensor<TPixel, VImageDimension> ( inputBaseFileName, verbose );

		if ( verbose ) std::cout<<"- Save the result in <"<<outputFileName.c_str()<<">..." << std::endl;
		typename TensorWriterType::Pointer tensorWriter = TensorWriterType::New();
		tensorWriter->SetFileName(outputFileName);
		tensorWriter->SetInput(outImage);
		tensorWriter->SetUseCompression(true);
		tensorWriter->Update();
	}





	template  < typename TPixel, unsigned int VImageDimension  >
	typename itk::Image< typename itk::DiffusionTensor3D< TPixel >, VImageDimension >::Pointer FSLStickToDiff3DTensor( const std::string& inputFileName, bool verbose=false )
	{
		crl::FileName inputFile(inputFileName);

		typedef typename itk::DiffusionTensor3D< TPixel >					TensorPixelType;
		typedef typename itk::Image< TensorPixelType, VImageDimension > 	TensorImageType;
		typedef typename itk::Image< TPixel, VImageDimension+1 >		 	FSLEigenVectorImageType;
		typedef typename itk::VectorImage< TPixel, VImageDimension >		EigenVectorImageType;
		
		//---------------------------------------------------
		// Read the eigen vectors
		//---------------------------------------------------
		if ( verbose ) std::cout<<"- Read eigen vector from <"<<inputFileName<<">..." << std::endl;
		typedef typename itk::ImageFileReader< FSLEigenVectorImageType > EigenVectorReaderType;
		typename EigenVectorReaderType::Pointer fslV1 = EigenVectorReaderType::New();
		fslV1->SetFileName( inputFileName );
		fslV1->Update();
		typename EigenVectorImageType::Pointer V1 = ConvertNDToVectorImage<TPixel, VImageDimension+1>(fslV1->GetOutput());

		//---------------------------------------------------
		// Create the itk::DiffusionTensor3D image
		//---------------------------------------------------
		if ( verbose ) std::cout << "- Convert to a Diffusion Tensor 3D Image..."  << std::endl;
		typename TensorImageType::Pointer outputImage = TensorImageType::New();
		outputImage->SetRegions(V1->GetLargestPossibleRegion());
		outputImage->SetOrigin(V1->GetOrigin());
		outputImage->SetSpacing(V1->GetSpacing());
		outputImage->SetDirection(V1->GetDirection());
		outputImage->Allocate();

		//---------------------------------------------------
		// Create the iterators
		//---------------------------------------------------
		typedef typename itk::ImageRegionConstIterator<EigenVectorImageType>	EigenVectorIteratorType;
		typedef typename itk::ImageRegionIterator<TensorImageType>				OutputIteratorType;

		EigenVectorIteratorType itV1(V1, V1->GetLargestPossibleRegion() );
		OutputIteratorType		itOut(outputImage, outputImage->GetLargestPossibleRegion());

		//---------------------------------------------------
		// And now convert!
		//---------------------------------------------------
		vnl_matrix<TPixel>		matrixOut(3,3);
		vnl_matrix<TPixel>		matrixV (3,3);
		vnl_diag_matrix<TPixel> matrixL(3) ;

		for ( itOut.GoToBegin(), itV1.GoToBegin()  ; 
			  !itOut.IsAtEnd() ; 
			  ++itOut, ++itV1 )
		{
			TensorPixelType diffTensor;

			// Get the vector
			typename EigenVectorImageType::PixelType vecV1  = itV1.Get();

			// Tips to get an orthogonal basis!
			vnl_vector_fixed<double,3> v1,v2,v3;
			for ( int i=0; i<3; i++ )
			{
				v1[i] = vecV1[i];
				v2[i] = v1[i] + ((float)rand()) / (float)(RAND_MAX);
			}
			v3 = vnl_cross_3d<double>(v1, v2);
			v2 = vnl_cross_3d<double>(v1, v3);

			// Transform into a vnl_matrix and a vnl_diag_matrix
			matrixV.fill(0);
			for ( int i=0; i<3; i++ )
			{
				matrixV(i,0) = v1[i];
				matrixV(i,1) = v2[i];
				matrixV(i,2) = v3[i];
			}

			matrixL[0] = 1e-3;
			matrixL[1] = 1e-6;
			matrixL[2] = 1e-6;

			// Compute the diffusion tensor matrix
			matrixOut = matrixV * matrixL * matrixV.transpose();

			// Transform into a diffusion tensor 
			TPixel *p = diffTensor.GetDataPointer ();
			*p = matrixOut(0,0);
			*(p+1) =  matrixOut(0,1);
			*(p+2) = matrixOut(0,2);
			*(p+3) = matrixOut(1,1);
			*(p+4) = matrixOut(1,2);
			*(p+5) = matrixOut(2,2);

			// And set the result in the output image!
			itOut.Set(diffTensor);
		}

		return outputImage;
	}

	template  < typename TPixel, unsigned int VImageDimension  >
	void FSLStickToDiff3DTensor( const std::string& inputBaseFileName, const std::string& outputFileName, bool verbose=false )
	{
		typedef typename itk::DiffusionTensor3D< TPixel >					TensorPixelType;
		typedef typename itk::Image< TensorPixelType, VImageDimension > 	TensorImageType;
		typedef typename itk::ImageFileWriter<TensorImageType>				TensorWriterType;

		typename TensorImageType::Pointer outImage = FSLStickToDiff3DTensor<TPixel, VImageDimension> ( inputBaseFileName, verbose );

		if ( verbose ) std::cout<<"- Save the result in <"<<outputFileName.c_str()<<">..." << std::endl;
		typename TensorWriterType::Pointer tensorWriter = TensorWriterType::New();
		tensorWriter->SetFileName(outputFileName);
		tensorWriter->SetInput(outImage);
		tensorWriter->SetUseCompression(true);
		tensorWriter->Update();
	}


} // end namespace DWI.
} // end namespace crl.

#endif


