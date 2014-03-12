#ifndef __crlMultiVectorResampleImageFilter_txx
#define __crlMultiVectorResampleImageFilter_txx

#include "itkVectorResampleImageFilter.h"
#include "itkObjectFactory.h"
#include "itkIdentityTransform.h"
#include "itkProgressReporter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMath.h"

#include "TensorLogFunctor.h"
#include "TensorExpFunctor.h"

#include "crlFastOps.h"
#include "crlTaskProgress.h"

namespace crl
{
	/**********************************************************************************************//**
	 * \fn	template <class TInputImage, class TOutputImage,
	 * 		class TInterpolatorPrecisionType> MFMTensorResampleImageFilter<TInputImage,
	 * 		TOutputImage, TInterpolatorPrecisionType> ::MFMTensorResampleImageFilter()
	 *
	 * \brief	Initialize new instance
	 *
	 * \author	Benoit Scherrer
	 * \date	June 2010
	 *
	 * \return	. 
	*************************************************************************************************/
	template <typename TPixel, unsigned int VImageDimension >
	MFMTensorResampleImageFilter<TPixel, VImageDimension>
		::MFMTensorResampleImageFilter()
	{
		m_OutputSpacing.Fill(1.0);
		m_OutputOrigin.Fill(0.0);
		m_OutputDirection.SetIdentity();
		m_Size.Fill( 0 );
		m_OutputStartIndex.Fill( 0 );
		m_ShowProgress = false;

		m_Transform = IdentityTransform<double, VImageDimension>::New();
		m_DefaultPixelValue.Fill(0);
		m_TensorInterpolationMode = INTERP_EUCLIDEAN;

		// Borland does not produce reliable results for multi processors
#ifdef __BORLANDC__
		this->SetNumberOfThreads(1);
#endif
	}

	/**********************************************************************************************//**
	 * \fn	template <class TInputImage, class TOutputImage,
	 * 		class TInterpolatorPrecisionType> void MFMTensorResampleImageFilter<TInputImage,
	 * 		TOutputImage,TInterpolatorPrecisionType> ::PrintSelf(std::ostream& os,
	 * 		Indent indent) const
	 *
	 * \brief	Print out a description of self
	 *
	 * \author	Benoit Scherrer
	 * \date	June 2010
	 *
	 * \param [in,out]	os	The operating system. 
	 * \param	indent		The indent. 
	 *
	 * \return	. 
	*************************************************************************************************/
	template <typename TPixel, unsigned int VImageDimension >
	void 
		MFMTensorResampleImageFilter<TPixel, VImageDimension>
		::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os,indent);

		os << indent << "Number of inputs vectors: "<<this->GetNumberOfInputs () << std::endl;
		os << indent << "DefaultPixelValue: "
			<< static_cast<typename NumericTraits<PixelType>::PrintType>(m_DefaultPixelValue)
			<< std::endl;
		os << indent << "Size: " << m_Size << std::endl;
		os << indent << "OutputStartIndex: " << m_OutputStartIndex << std::endl;
		os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
		os << indent << "OutputOrigin: " << m_OutputOrigin << std::endl;
		os << indent << "OutputDirection: " << m_OutputDirection << std::endl;
		os << indent << "Transform: " << m_Transform.GetPointer() << std::endl;

		if ( m_TensorInterpolationMode==INTERP_NONE )
			os << indent << "Interpolation: none" << std::endl;
		else if ( m_TensorInterpolationMode==INTERP_EUCLIDEAN )
			os << indent << "Interpolation: Euclidean" << std::endl;
		else if ( m_TensorInterpolationMode==INTERP_LOGEUCLIDEAN )
			os << indent << "Interpolation: log-Euclidean" << std::endl;
		else
			os << indent << "Interpolation: Unknown (error?)" << std::endl;

		return;
	}

	/**********************************************************************************************//**
	 * \fn	template <class TInputImage, class TOutputImage,
	 * 		class TInterpolatorPrecisionType> void MFMTensorResampleImageFilter<TInputImage,
	 * 		TOutputImage,TInterpolatorPrecisionType> ::SetOutputSpacing(const double* spacing)
	 *
	 * \brief	Set the output image spacing.
	 *
	 * \author	Benoit Scherrer
	 * \date	June 2010
	 *
	 * \param	spacing	The spacing. 
	 *
	 * \return	. 
	*************************************************************************************************/
	template <typename TPixel, unsigned int VImageDimension >
	void 
		MFMTensorResampleImageFilter<TPixel, VImageDimension>
		::SetOutputSpacing(const double* spacing)
	{
		SpacingType s(spacing);
		this->SetOutputSpacing( s );
	}

	/**********************************************************************************************//**
	 * \fn	template <class TInputImage, class TOutputImage,
	 * 		class TInterpolatorPrecisionType> void MFMTensorResampleImageFilter<TInputImage,
	 * 		TOutputImage,TInterpolatorPrecisionType> ::SetOutputOrigin(const double* origin)
	 *
	 * \brief	Set the output image origin.
	 *
	 * \author	Benoit Scherrer
	 * \date	June 2010
	 *
	 * \param	origin	The origin. 
	 *
	 * \return	. 
	*************************************************************************************************/
	template <typename TPixel, unsigned int VImageDimension >
	void 
		MFMTensorResampleImageFilter<TPixel, VImageDimension>
		::SetOutputOrigin(const double* origin)
	{
		PointType p(origin);
		this->SetOutputOrigin( p );
	}

	/**********************************************************************************************//**
	 * \fn	template <class TInputImage, class TOutputImage,
	 * 		class TInterpolatorPrecisionType> void MFMTensorResampleImageFilter<TInputImage,
	 * 		TOutputImage,TInterpolatorPrecisionType> ::BeforeThreadedGenerateData()
	 *
	 * \brief	Set up state of filter before multi-threading.
	 *
	 * \author	Benoit Scherrer
	 * \date	June 2010
	 *
	 * \return	. 
	*************************************************************************************************/
	template <typename TPixel, unsigned int VImageDimension >
	void 
		MFMTensorResampleImageFilter<TPixel, VImageDimension>
		::BeforeThreadedGenerateData()
	{
		TaskProgress::InitTaskProgress();
	}

	/**********************************************************************************************//**
	 * \fn	template <class TInputImage, class TOutputImage,
	 * 		class TInterpolatorPrecisionType> void MFMTensorResampleImageFilter<TInputImage,
	 * 		TOutputImage,TInterpolatorPrecisionType> ::AfterThreadedGenerateData()
	 *
	 * \brief	Set up state of filter after multi-threading.
	 *
	 * \author	Benoit Scherrer
	 * \date	June 2010
	 *
	 * \return	. 
	*************************************************************************************************/
	template <typename TPixel, unsigned int VImageDimension >
	void 
		MFMTensorResampleImageFilter<TPixel, VImageDimension>
		::AfterThreadedGenerateData()
	{

	}

	/**********************************************************************************************//**
	 * \fn	template <class TInputImage, class TOutputImage,
	 * 		class TInterpolatorPrecisionType> void MFMTensorResampleImageFilter<TInputImage,
	 * 		TOutputImage,
	 * 		TInterpolatorPrecisionType> ::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread,
	 * 		int threadId)
	 *
	 * \brief	Threaded generate data. 
	 *
	 * \author	Benoit Scherrer
	 * \date	June 2010
	 *
	 * \param	outputRegionForThread	The output region for thread. 
	 * \param	threadId				Identifier for the thread. 
	 *
	 * \return	. 
	*************************************************************************************************/
	template <typename TPixel, unsigned int VImageDimension >
	void 
		MFMTensorResampleImageFilter<TPixel, VImageDimension>
		::ThreadedGenerateData(
		const OutputImageRegionType& outputRegionForThread,
		int threadId)
	{
		itkDebugMacro(<<"Actually executing");
		unsigned int nbTensors = this->GetNumberOfInputs();

		//--------------------------------------------------
		// Get the output pointers
		//--------------------------------------------------
		std::vector<OutputImagePointer> outputPtrs(nbTensors);
		for ( unsigned int i=0; i<this->GetNumberOfOutputs() ; i++ )
			outputPtrs[i] = this->GetOutput(i);

		//--------------------------------------------------
		// Get ths input pointers
		//--------------------------------------------------
		std::vector<InputImageConstPointer> inputPtrs(nbTensors);
		for ( unsigned int i=0; i<nbTensors ; i++ )
			inputPtrs[i] = (this->GetInput(i));

		//--------------------------------------------------
		// Create the iterators that will walk the output region for this thread.
		//--------------------------------------------------
		typedef ImageRegionIteratorWithIndex<TOutputImage> OutputIterator;
		std::vector<OutputIterator> outIts;
		for ( unsigned int i=0; i<outputPtrs.size() ; i++ )
			outIts.push_back( OutputIterator(outputPtrs[i], outputRegionForThread));

		//-------------------------------------------------
		// Creates the rotation matrix from the transform parameters
		//-------------------------------------------------
		typename TransformType::ParametersType p;
		p = m_Transform->GetParameters();
		
		vnl_matrix<double> m_AMatrix(3,3);
		m_AMatrix(0,0) = p[0]; m_AMatrix(0,1) = p[1]; m_AMatrix(0,2) = p[2];
		m_AMatrix(1,0) = p[3]; m_AMatrix(1,1) = p[4]; m_AMatrix(1,2) = p[5];
		m_AMatrix(2,0) = p[6]; m_AMatrix(2,1) = p[7]; m_AMatrix(2,2) = p[8];

		// Computes the SVD (protected by a mutex because not thread safe)
		mutexAccessSVD.Lock();
		vnl_svd<double> svd( m_AMatrix, 1e-8 );
		mutexAccessSVD.Unlock();
		vnl_matrix<double> m_VMatrix = svd.V();
		vnl_matrix<double> m_UMatrix = svd.U();
		vnl_matrix<double> RotationMatrix = m_UMatrix*m_VMatrix.transpose();

		//--------------------------------------------------
		// Define a few indices that will be used to translate from an input pixel
		// to an output pixel
		//--------------------------------------------------
		PointType outputPoint;         // Coordinates of current output pixel
		PointType inputPoint;          // Coordinates of current input pixel
		ContinuousIndexType cindex;

		//--------------------------------------------------
		// Create the output tensors and fractions
		//--------------------------------------------------
		std::vector<PixelType> tensors(nbTensors);
		std::vector<itk::Vector< double, 3> > tensorVector(nbTensors);

		typename FractionsImageType::PixelType fractions;
		bool resampleFractions = ( m_FractionsImage.GetPointer()!=NULL);
		if ( resampleFractions )
			fractions.SetSize(m_FractionsImage->GetNumberOfComponentsPerPixel());

		//--------------------------------------------------
		// Create the iterators on the output fractions
		//--------------------------------------------------
		typedef ImageRegionIteratorWithIndex<FractionsImageType> OutputFractionsIterator;
		OutputFractionsIterator outFractionsIt;
		if ( resampleFractions )
			outFractionsIt = OutputFractionsIterator(m_ResampledFractionsImage, outputRegionForThread);

		//--------------------------------------------------
		// For the log-Euclidean interpolation
		//--------------------------------------------------
		itk::Functor::TensorExp<PixelType, PixelType> tensorExp;
		itk::Functor::TensorLog<PixelType,	PixelType> tensorLog;

		//--------------------------------------------------
		// Support for progress methods/callbacks
		//--------------------------------------------------
		// ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
		// ProgressReporter is really bad with multi threaded filters! (cf its source code)
		TaskProgress taskProgress(threadId, this->GetNumberOfThreads(), outputRegionForThread.GetNumberOfPixels(), true, false);
		taskProgress.SetLinePrefix("    ");
		if (m_ShowProgress)  taskProgress.InitShowProgress();

		//--------------------------------------------------
		// Walk the output region of all iterators
		//--------------------------------------------------
		outFractionsIt.GoToBegin();
		for ( unsigned int i=0; i<outIts.size() ; i++ )
			outIts[i].GoToBegin();

		while ( !outIts[0].IsAtEnd() )
		{
			if (m_ShowProgress) taskProgress.Update();

			//---------------------------------------------
			// Determine the index of the current output pixel
			//---------------------------------------------
			outputPtrs[0]->TransformIndexToPhysicalPoint( outIts[0].GetIndex(), outputPoint );

			//---------------------------------------------
			// Initialize (for both inside or outside cases)
			//---------------------------------------------
			fractions.Fill(0);
			for ( unsigned int i=0; i<nbTensors ; i++ )
				tensors[i].Fill(0);

			//---------------------------------------------
			// Compute the corresponding input pixel position
			//---------------------------------------------
			inputPoint = m_Transform->TransformPoint(outputPoint);
			bool isInside = inputPtrs[0]->TransformPhysicalPointToContinuousIndex(inputPoint, cindex);

			if( isInside )
			{
				//---------------------------------------------
				// Convert the input continuous index to an input index
				//---------------------------------------------
				IndexType index;
				index.CopyWithRound( cindex );

				//---------------------------------------------
				// Get the fractions for the input index if required
				//---------------------------------------------
				if ( resampleFractions )
					fractions = m_FractionsImage->GetPixel(index);







				//------------------------------------------------------
				// Compute the maj eigen vector of each input tensor
				// to feed SelectMatchingTensor
				// (only needed if more than one input)
				//
				// RQ: matching tensors are compared in term of similar
				// primary eigen vector. They could also be compared
				// with the log-euclidean metric
				//------------------------------------------------------
				if ( nbTensors>1 )
				{
					for ( unsigned int i=0; i<nbTensors ; i++ )
					{
						const PixelType& tensor = inputPtrs[i]->GetPixel(index);
						EstimatePrimaryEigenvector(tensor, tensorVector[i] );
					}
				}

				//------------------------------------------------------
				// If interpolation, uses the continuous index
				//------------------------------------------------------
				if ( m_TensorInterpolationMode!=INTERP_NONE )
				{
					//------------------------------------------------------
					// here we implement a Euclidean/logEuclidean interpolation.
					// The tensors at the corners of the voxel are selected based on similarity
					// of their primary eigenvector to the input tensor direction.
					//------------------------------------------------------

					//------------------------------------------------------
					// Compute base index = closest index below point
					// Compute distance from point to base index (from the point to the voxel center)
					//------------------------------------------------------
					unsigned int dim;						// index over dimension
					signed long baseIndex[VImageDimension];
					double distance[VImageDimension];

					for( dim = 0; dim < VImageDimension; dim++ )
					{
						baseIndex[dim] = Math::Floor< typename IndexType::IndexValueType >( cindex[dim] );
						distance[dim] = cindex[dim] - static_cast< double >( baseIndex[dim] );

						//				baseIndex[dim] = long (vnl_math_rnd( cindex[dim] ));
						//				distance[dim] = cindex[dim] - double( baseIndex[dim] ) + 0.5;
					}

					//------------------------------------------------------
					// Interpolated value is the weight some of each of the surrounding
					// neighbors. The weight for each neighbour is the fraction overlap
					// of the neighbor pixel with respect to a pixel centered on point.
					//------------------------------------------------------
					double totalOverlap = 0.0;
					for( unsigned int counter = 0; counter < (1<<VImageDimension); counter++ )
					{

						double			overlap = 1.0;          // fraction overlap
						unsigned int	upper = counter;		// each bit indicates upper/lower neighbour
						IndexType		neighIndex;

						//-------------------------------------------------
						// get neighbor index and overlap fraction
						//-------------------------------------------------
						for( dim = 0; dim < VImageDimension; dim++ )
						{
							if ( upper & 1 )
							{
								neighIndex[dim] = baseIndex[dim] + 1;
#ifdef ITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY
								// Take care of the case where the pixel is just
								// in the outer upper boundary of the image grid.
								/*if( neighIndex[dim] > this->m_EndIndex[dim] )
								{
								neighIndex[dim] = this->m_EndIndex[dim];
								}*/
#endif
								overlap *= distance[dim];
							}
							else
							{
								neighIndex[dim] = baseIndex[dim];
#ifdef ITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY
								// Take care of the case where the pixel is just
								// in the outer lower boundary of the image grid.
								/*if( neighIndex[dim] < this->m_StartIndex[dim] )
								{
								neighIndex[dim] = this->m_StartIndex[dim];
								}*/
#endif
								overlap *= 1.0 - distance[dim];
							}

							upper >>= 1;
						}

						//-------------------------------------------------
						// get neighbor value only if overlap is not zero
						//-------------------------------------------------
						if( overlap )
						{
							// Skip illegal memory accesses - same as setting outside to be zero.
							if ( this->GetInput(0)->GetLargestPossibleRegion().IsInside( neighIndex )) 
							{
								PixelType cornerTensor;
								typename FractionsImageType::PixelType cornerFraction;
								if ( resampleFractions )
									cornerFraction = m_FractionsImage->GetPixel(neighIndex);

								//-------------------------------------------------
								// For each tensor, update 
								//-------------------------------------------------
								for ( unsigned int t=0; t<nbTensors ; t++ )
								{ 
									int idTensorMatch = SelectMatchingTensor(tensorVector[t], neighIndex, cornerTensor);

									//---------------------------------------
									// Tensor interpolation: log-Euclidean:
									//---------------------------------------
									if ( m_TensorInterpolationMode==INTERP_LOGEUCLIDEAN ) 
									{
										PixelType logCornerTensor = tensorLog( cornerTensor );
										for(unsigned int k = 0; k < PixelType::InternalDimension; k++ ) 
											tensors[t][k] += overlap * static_cast<double>( logCornerTensor[k] );
									}
									//---------------------------------------
									// Tensor interpolation: Euclidean:
									//---------------------------------------
									else {
										for(unsigned int k = 0; k < PixelType::InternalDimension; k++ ) 
											tensors[t][k] += overlap * static_cast<double>( cornerTensor[k] );		  
									}

									//---------------------------------------
									// Update the tensor fractions
									//---------------------------------------
									if ( resampleFractions ) 
										fractions[t] += overlap * cornerFraction[idTensorMatch];
								}

								//-------------------------------------------------
								// Update also for the other fractions if needed
								//-------------------------------------------------
								if ( resampleFractions && ( fractions.GetSize() > nbTensors ) ) 
								{
									for ( unsigned int i=nbTensors ; i< fractions.GetSize() ; i++ )
										fractions[i] += overlap * cornerFraction[i];
								}
							}
							totalOverlap += overlap;
						}

						if( totalOverlap == 1.0 )
						{
							// finished
							break;
						}
					} //end for( unsigned int counter = 0; counter < (1<<VImageDimension); counter++ )

					//----------------------------------------------------
					// Converts back from LE if needed
					//----------------------------------------------------
					if ( m_TensorInterpolationMode==INTERP_LOGEUCLIDEAN )
					{
						for ( unsigned int t=0; t<nbTensors ; t++ )
						{
							if ( tensors[t].GetTrace()!=0 )
								tensors[t] = tensorExp( tensors[t] );
						}
					}
				}

				//------------------------------------------------------------
				// NO INTERPOLATION : m_TensorInterpolationMode == INTERP_NONE
				//------------------------------------------------------------
				else
				{
					for ( unsigned int t=0; t<nbTensors ; t++ )
					{
						// Select the matching tensor for that index (nearest neighbor interp)
						tensors[t] = inputPtrs[t]->GetPixel(index);
					}

					if ( resampleFractions ) 
						fractions = m_FractionsImage->GetPixel(index);

				}




				//---------------------------------------------
				// Now writes the resulting tensors/fractions
				//---------------------------------------------
				PixelType rotatedTensor; 
				for ( unsigned int t=0; t<nbTensors ; t++ )
				{
					// Rotate the tensor
					RotateOneTensor(tensors[t], RotationMatrix, rotatedTensor);
					outIts[t].Set( rotatedTensor );

				}
				
				//---------------------------------------------
				// Normalize and set fractions if required
				//---------------------------------------------
				if ( resampleFractions )
				{
					double sumF = 0.0;

					for ( unsigned int i=0; i<fractions.GetSize() ; i++ )
						sumF += fractions[i];		// to renormalize the fractions?

					if (sumF!=0)
					{
						for ( unsigned int i=0; i<fractions.GetSize() ; i++ )
							fractions[i] = fractions[i]/sumF;
					}

					outFractionsIt.Set(fractions);
				}
			}
			//if( isInside )
			else
			{
				// the fractions were set to NULL before
				if ( resampleFractions )
					outFractionsIt.Set(fractions);

				for ( unsigned int i=0; i<nbTensors ; i++ )
					outIts[i].Set(m_DefaultPixelValue); // default background value
			}

			//if (threadId==0) progress.CompletedPixel();

			if ( resampleFractions )
				++outFractionsIt;
			
			for ( unsigned int i=0; i<outputPtrs.size() ; i++ )
				++(outIts[i]);
		}

		//--------------------------------------------------
		// Transfer the meta data dictionary from the input
		// (This aims to transfer the tensor measurement frame from the NRRD data)
		//--------------------------------------------------
		for ( unsigned int i=0; i<this->GetNumberOfOutputs() ; i++ )
			outputPtrs[i]->SetMetaDataDictionary( inputPtrs[i]->GetMetaDataDictionary() );

		tensors.clear();
		outIts.clear();
		outputPtrs.clear();
		inputPtrs.clear();
		return;
	}

	/**********************************************************************************************//**
	 *
	 * \brief	Select in all inputs the best tensor matching refVector at the index index (based on
	 * 			the primary eigen vector). 
	 *
	 * \author	Benoit Scherrer
	 * \date	June 2010
	 *
	 * \param	refVector	The reference vector. 
	 * \param	index		Zero-based index in the input images. 
	 * \param [in,out]	out	out	The output tensor best matching refVector. 
	 *
	 * \return	The 0-based index of the selected tensor 
	*************************************************************************************************/
	template <typename TPixel, unsigned int VImageDimension >
	int MFMTensorResampleImageFilter<TPixel, VImageDimension>
		::SelectMatchingTensor(  const itk::Vector< double, 3> &refVector, const IndexType index, PixelType &out ) 
	{
		int idTensor = -1;
		double projection = -1;
		itk::Vector< double, 3> tensorVector;

		//------------------------------------------------------------------
		// If only one input, just return the tensor
		//------------------------------------------------------------------
		if ( this->GetNumberOfInputs()==1 )
		{
			out = this->GetInput(0)->GetPixel(index);
			return 0;
		}

		//------------------------------------------------------------------
		// Looks in all the input at the index <index>
		//------------------------------------------------------------------
		for ( unsigned int i=0; i<this->GetNumberOfInputs() ; i++ )
		{
			//-----------------------------------------------------
			// Take the tensor and computes its major eigen value
			//-----------------------------------------------------
			PixelType oneTensor = this->GetInput(i)->GetPixel(index);
			EstimatePrimaryEigenvector(oneTensor, tensorVector );

			//-----------------------------------------------------
			// Compare the projections on the reference vector
			//-----------------------------------------------------
			double proj2 = fabs(refVector * tensorVector);
			if (proj2 > projection) 
			{
				projection = proj2;
				idTensor = i;
				out = oneTensor;
			}
		}

		if ( idTensor<0 ) std::cout<<"Strange error. Cannot match the tensor"<<std::endl;
		return idTensor;
	}

		/**********************************************************************************************//**
		 * \fn	template <typename TPixel,
		 * 		unsigned int VImageDimension > void MFMTensorResampleImageFilter<TPixel,
		 * 		VImageDimension> ::EstimatePrimaryEigenvector( PixelType &in, itk::Vector< double,
		 * 		3> &evec)
		 *
		 * \brief	Computes the primary eigenvector of a tensor. 
		 *
		 * \author	Benoit Scherrer
		 * \date	June 2010
		 *
		 * \param [in,out]	in		The input tensor. 
		 * \param [in,out]	evec	The primary eigen vector. 
		 *
		 * \return	. 
		*************************************************************************************************/
	template <typename TPixel, unsigned int VImageDimension >
	void 
		MFMTensorResampleImageFilter<TPixel, VImageDimension>
		::EstimatePrimaryEigenvector( const PixelType &in, itk::Vector< double, 3> &evec) const
	{
		typedef typename PixelType::ValueType				ValueType;
		typedef typename PixelType::EigenVectorsMatrixType	MatrixType;
		typedef typename PixelType::EigenValuesArrayType	EigenValuesType;

		MatrixType V;
		EigenValuesType e;
		in.ComputeEigenAnalysis(e, V);

		// The vnl convention is for the eigenvalues to be ordered from
		// the smallest to the largest, and for the eigenvectors to be in the rows.
		evec[0] = V[2][0];
		evec[1] = V[2][1];
		evec[2] = V[2][2];
	}

	/**********************************************************************************************//**
	 * \fn	template <typename TPixel,
	 * 		unsigned int VImageDimension > void MFMTensorResampleImageFilter<TPixel,
	 * 		VImageDimension> ::InterpolateMatchingTensor( const TensorImageType::PixelType &in,
	 * 		const ContinuousIndexType &cindex, TensorImageType::PixelType &out)
	 *
	 * \brief	Interpolate matching tensor. 
	 *
	 * \author	Benoit Scherrer
	 * \date	June 2010
	 *
	 * \param	in			The in. 
	 * \param	cindex		The cindex. 
	 * \param [in,out]	out	The out. 
	 *
	 * \return	. 
	*************************************************************************************************/
	template <typename TPixel, unsigned int VImageDimension >
	void 
		MFMTensorResampleImageFilter<TPixel, VImageDimension>
		::InterpolateMatchingTensor(
			int idTensorIn, 
			const PixelType &in,
			const ContinuousIndexType &cindex,
			PixelType &out ,										// output tensor.
			typename FractionsImageType::PixelType &outFractions )  // output fractions.
	{
		//------------------------------------------------------
		// Compute the maj eigen vector of the input tensor
		// to feed SelectMatchingTensor
		// (only needed if more than one input)
		//
		// RQ: matching tensors are compared in term of similar
		// primary eigen vector. They could also be compared
		// with the log-euclidean metric
		//------------------------------------------------------
		itk::Vector< double, 3>  inVector;
		if ( this->GetNumberOfInputs()>1 ) EstimatePrimaryEigenvector(in, inVector);

		//------------------------------------------------------
		// Initialize the fractions resampling
		//------------------------------------------------------
		bool resampleFractions = (outFractions.GetSize()>0);
		// The fractions were already initialized to 0 before the call of InterpolateMatchingTensor

		//------------------------------------------------------
		// If interpolation, uses the continuous index
		//------------------------------------------------------
		if ( m_TensorInterpolationMode!=INTERP_NONE )
		{
			//------------------------------------------------------
			// here we implement a Euclidean/logEuclidean interpolation.
			// The tensors at the corners of the voxel are selected based on similarity
			// of their primary eigenvector to the input tensor direction.
			//------------------------------------------------------
			itk::Functor::TensorExp<PixelType, PixelType> tensorExp;
			itk::Functor::TensorLog<PixelType,	PixelType> tensorLog;

			//------------------------------------------------------
			// Compute base index = closest index below point
			// Compute distance from point to base index (from the point to the voxel center)
			//------------------------------------------------------
			unsigned int dim;						// index over dimension
			signed long baseIndex[VImageDimension];
			double distance[VImageDimension];

			for( dim = 0; dim < VImageDimension; dim++ )
			{
				baseIndex[dim] = Math::Floor< typename IndexType::IndexValueType >( cindex[dim] );
				distance[dim] = cindex[dim] - static_cast< double >( baseIndex[dim] );

//				baseIndex[dim] = long (vnl_math_rnd( cindex[dim] ));
//				distance[dim] = cindex[dim] - double( baseIndex[dim] ) + 0.5;
			}

			//------------------------------------------------------
			// Interpolated value is the weight some of each of the surrounding
			// neighbors. The weight for each neighbour is the fraction overlap
			// of the neighbor pixel with respect to a pixel centered on point.
			//------------------------------------------------------
			PixelType output;
			output.Fill( 0.0 );

			double totalOverlap = 0.0;
			for( unsigned int counter = 0; counter < (1<<VImageDimension); counter++ )
			{

				double			overlap = 1.0;          // fraction overlap
				unsigned int	upper = counter;		// each bit indicates upper/lower neighbour
				IndexType		neighIndex;

				//-------------------------------------------------
				// get neighbor index and overlap fraction
				//-------------------------------------------------
				for( dim = 0; dim < VImageDimension; dim++ )
				{
					if ( upper & 1 )
					{
						neighIndex[dim] = baseIndex[dim] + 1;
#ifdef ITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY
						// Take care of the case where the pixel is just
						// in the outer upper boundary of the image grid.
						/*if( neighIndex[dim] > this->m_EndIndex[dim] )
						{
							neighIndex[dim] = this->m_EndIndex[dim];
						}*/
#endif
						overlap *= distance[dim];
					}
					else
					{
						neighIndex[dim] = baseIndex[dim];
#ifdef ITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY
						// Take care of the case where the pixel is just
						// in the outer lower boundary of the image grid.
						/*if( neighIndex[dim] < this->m_StartIndex[dim] )
						{
							neighIndex[dim] = this->m_StartIndex[dim];
						}*/
#endif
						overlap *= 1.0 - distance[dim];
					}

					upper >>= 1;
				}

				//-------------------------------------------------
				// get neighbor value only if overlap is not zero
				//-------------------------------------------------
				if( overlap )
				{
					// Skip illegal memory accesses - same as setting outside to be zero.
					if ( this->GetInput(0)->GetLargestPossibleRegion().IsInside( neighIndex )) 
					{
						PixelType cornerTensor;
						int idTensorMatch = SelectMatchingTensor(inVector, neighIndex, cornerTensor);

						// Tensor interpolation: log-Euclidean:
						if ( m_TensorInterpolationMode==INTERP_LOGEUCLIDEAN ) 
						{
							PixelType logCornerTensor = tensorLog( cornerTensor );
							for(unsigned int k = 0; k < PixelType::InternalDimension; k++ ) 
								output[k] += overlap * static_cast<double>( logCornerTensor[k] );
						}
						// Tensor interpolation: Euclidean:
						else {
							for(unsigned int k = 0; k < PixelType::InternalDimension; k++ ) 
								output[k] += overlap * static_cast<double>( cornerTensor[k] );		  
						}

						// Fractions interpolation
						if ( resampleFractions ) 
						{
							const typename FractionsImageType::PixelType& cornerFraction = m_FractionsImage->GetPixel(neighIndex);
							outFractions[idTensorIn] += overlap * cornerFraction[idTensorMatch];

							// Trick: for the idTensorIn==0, resample also all the other fractions
							// (if there are other fractions...)
							if ( idTensorIn==0 && (this->GetNumberOfInputs()>outFractions.GetSize()) )
							{
								for ( int i=this->GetNumberOfInputs(); i< outFractions.GetSize() ; i++ )
									outFractions[i] += overlap * cornerFraction[i];
							}
						}
					}
					totalOverlap += overlap;
				}

				if( totalOverlap == 1.0 )
				{
					// finished
					break;
				}
			}

			// Now set the output
			if ( m_TensorInterpolationMode==INTERP_LOGEUCLIDEAN )
			{
				if ( output.GetTrace()==0 )
					out=output;
				else
					out = tensorExp( output );
			}
			else
				out = output;
		}

		// NO INTERPOLATION : m_TensorInterpolationMode == INTERP_NONE
		else
		{
			// Convert the continuous index to an index
			IndexType index;
			index.CopyWithRound( cindex );

			// Select the matching tensor for that index (nearest neighbor interp)
			SelectMatchingTensor(inVector, index, out); 

			if ( resampleFractions ) 
			{
				const typename FractionsImageType::PixelType& fractions = m_FractionsImage->GetPixel(index);
				outFractions[idTensorIn] = fractions[idTensorIn];

							// Trick: for the idTensorIn==0, resample also all the other fractions
							// (if there are other fractions...)
							if ( idTensorIn==0 && (this->GetNumberOfInputs()>outFractions.GetSize()) )
							{
								for ( int i=this->GetNumberOfInputs(); i< outFractions.GetSize() ; i++ )
									outFractions[i] += fractions[i];
							}
			}
		}


	}

	/**********************************************************************************************//**
	 *
	 * \brief	Rotate tensors. Rotate the tensors at each pixel according to the transformation:
	 * 			Tnew = Rt*T*R the rotation matrix is extracted from the affine transform parameters,
	 * 			using the method described by Raimundo Sierra, in his Master's Thesis. See also the
	 * 			paper by Juan Ruiz-Alzola. The matrix A of the affine transform is decomposed into a
	 * 			pure rotation part and another part. The rotation is used to reorient the rensors.
	 * 			
	 * 			Using svd: A = R*W = UWVt = UVt * VWVt, where UVt is the rotation matrix. 
	 *
	 * \author	Benoit Scherrer
	 * \date	June 2010
	 *
	 * \return	. 
	*************************************************************************************************/
	template <typename TPixel, unsigned int VImageDimension >
	void 
		MFMTensorResampleImageFilter<TPixel, VImageDimension>
		::RotateOneTensor(const PixelType& inTensor, const vnl_matrix<double>& R, PixelType& outTensor )
	{
		//------------------------------------------
		// Transform the PixelType to a vnl_matrix
		//------------------------------------------
		vnl_matrix<double> vnlTensorMatrix(3,3), vnlOut(3,3);
		double *o = vnlTensorMatrix.data_block();
		const typename PixelType::ComponentType *p = inTensor.GetDataPointer ();
		*(o) = *(p); *(o+1) = *(p+1); *(o+2) = *(p+2);
		*(o+3) = *(p+1); *(o+4) = *(p+3); *(o+5) = *(p+4);
		*(o+6) = *(p+2); *(o+7) = *(p+4); *(o+8) = *(p+5);

		//------------------------------------------
		// Apply the rotation (fast)
		//------------------------------------------		
		crl::FastOps::AtBA_3x3( vnlOut, R, vnlTensorMatrix );

		//------------------------------------------
		// Converts back to a PixelType
		//------------------------------------------
		outTensor[0] = vnlOut(0,0);
		outTensor[1] = vnlOut(0,1);
		outTensor[2] = vnlOut(0,2);
		outTensor[3] = vnlOut(1,1);
		outTensor[4] = vnlOut(1,2);
		outTensor[5] = vnlOut(2,2);
	}

	/** 
	* Inform pipeline of necessary input image region
	*
	* Determining the actual input region is non-trivial, especially
	* when we cannot assume anything about the transform being used.
	* So we do the easy thing and request the entire input image.
	*/
	template <typename TPixel, unsigned int VImageDimension >
	void 
		MFMTensorResampleImageFilter<TPixel, VImageDimension>
		::GenerateInputRequestedRegion()
	{
		// call the superclass's implementation of this method
		Superclass::GenerateInputRequestedRegion();

		// Prepare inputs
		for ( unsigned int i=0; i<this->GetNumberOfInputs() ; i++ )
		{
			if ( this->GetInput(i)!=NULL )
			{
				InputImagePointer  inputPtr  = 
					const_cast< TInputImage *>( this->GetInput(i) );

				// Request the entire input image
				InputImageRegionType inputRegion;
				inputRegion = inputPtr->GetLargestPossibleRegion();
				inputPtr->SetRequestedRegion(inputRegion);
			}
		}
		return;
	}

	/** 
	* Inform pipeline of required output region
	*/
	template <typename TPixel, unsigned int VImageDimension >
	void 
		MFMTensorResampleImageFilter<TPixel, VImageDimension>
		::GenerateOutputInformation()
	{
		// Creates the additional outputs if required
		unsigned int prevNum = this->GetNumberOfOutputs();
		unsigned int nbOutputs = this->GetNumberOfInputs();
		this->SetNumberOfOutputs( nbOutputs );
		for ( unsigned int i = prevNum; i < nbOutputs; ++i )
		{
			this->SetNthOutput( i, this->MakeOutput( i ).GetPointer() );
		 }

		// call the superclass' implementation of this method
		Superclass::GenerateOutputInformation();

		// get pointers to the input and output
		typename TOutputImage::RegionType outputLargestPossibleRegion;
		outputLargestPossibleRegion.SetSize( m_Size );
		outputLargestPossibleRegion.SetIndex( m_OutputStartIndex );

		for ( unsigned int i=0; i<this->GetNumberOfOutputs() ; i++ )
		{
			OutputImagePointer outputPtr = this->GetOutput(i);
			if ( outputPtr )
			{
				// Set the size of the output region
				outputPtr->SetLargestPossibleRegion( outputLargestPossibleRegion );

				// Set spacing and origin
				outputPtr->SetSpacing( m_OutputSpacing );
				outputPtr->SetOrigin( m_OutputOrigin );
				outputPtr->SetDirection( m_OutputDirection );

				//outputPtr->Allocate();
			}
		}

		//---------------------------------------------
		// If the fractions image is defined, allocate
		// the output resampled fractions image
		//---------------------------------------------
		if ( m_FractionsImage.GetPointer()!=NULL )
		{
			m_ResampledFractionsImage = FractionsImageType::New();
			m_ResampledFractionsImage->SetRegions(outputLargestPossibleRegion);
			m_ResampledFractionsImage->SetNumberOfComponentsPerPixel(m_FractionsImage->GetNumberOfComponentsPerPixel());
			m_ResampledFractionsImage->SetSpacing( m_OutputSpacing );
			m_ResampledFractionsImage->SetOrigin( m_OutputOrigin );
			m_ResampledFractionsImage->SetDirection( m_OutputDirection );
			m_ResampledFractionsImage->Allocate();
		}
		return;
	}

	/** 
	* Verify if any of the components has been modified.
	*/
	template <typename TPixel, unsigned int VImageDimension >
	unsigned long 
		MFMTensorResampleImageFilter<TPixel, VImageDimension>
		::GetMTime( void ) const
	{
		unsigned long latestTime = Object::GetMTime(); 

		if( m_Transform )
		{
			if( latestTime < m_Transform->GetMTime() )
			{
				latestTime = m_Transform->GetMTime();
			}
		}

		return latestTime;
	}

} // end namespace itk

#endif
