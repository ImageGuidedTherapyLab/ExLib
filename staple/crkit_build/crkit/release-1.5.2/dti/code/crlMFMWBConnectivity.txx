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

#ifndef _CRL_MFMWBCONNECTIVITY_TXX_INCLUDED
#define _CRL_MFMWBCONNECTIVITY_TXX_INCLUDED 1

#ifdef _MSC_VER
  #pragma warning( disable: 4996 )
#endif

#include "crlDilateLabelImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "crlSingleLineMultipleCout.h"

/**********************************************************************************************//**
 * \fn	ITK_THREAD_RETURN_TYPE ThreadedAnalysis( void * arg )
 *
 * \brief	Threaded analysis subroutine. 
 *
 * \author	Benoit Scherrer
 * \date	March 2011
 *
 * \param [in,out]	arg	If non-null, the argument. 
 *
 * \return	. 
*************************************************************************************************/
ITK_THREAD_RETURN_TYPE ThreadedAnalysis( void * arg )
{
	typedef crl::MFMWBConnectivity<3, float> ConnectivityAnalysisFloatType;
	typedef crl::MFMWBConnectivity<3, double> ConnectivityAnalysisDoubleType;

	typedef itk::MultiThreader::ThreadInfoStruct  ThreadInfoType;
	ThreadInfoType * infoStruct = static_cast< ThreadInfoType * >( arg );
	crl::MFMTractGeneratorBase *algo = static_cast< crl::MFMTractGeneratorBase * >(infoStruct->UserData);

	ConnectivityAnalysisFloatType *algoFloat = dynamic_cast<ConnectivityAnalysisFloatType *>(algo);
	ConnectivityAnalysisDoubleType *algoDouble = dynamic_cast<ConnectivityAnalysisDoubleType *>(algo);

	if ( algoFloat!=NULL ) algoFloat->ThreadedTractAnalysis(arg);
	else if ( algoDouble!=NULL ) algoDouble->ThreadedTractAnalysis(arg);
	else
		std::cout<<"ERROR in ThreadedAnalysis"<<std::endl;

	return 0;
}

/**********************************************************************************************//**
 * \fn	template <const unsigned int TImageDimension,
 * 		class TComponentType> int crl::MFMWBConnectivity<TImageDimension,
 * 		TComponentType> ::LoadParcellationImage()
 *
 * \brief	Load the parcellation image used to compute the connectivity matrix. Dilate the
 * 			labels if necessary.Compute the Index->Label and Label->Index tables. Allocate the
 * 			connectivity matrix and the color matrix (random colors). 
 *
 * \author	Benoit Scherrer
 * \date	December 2010
 *
 * \return	The calculated seeding region. 
*************************************************************************************************/
template <const unsigned int TImageDimension, class TComponentType>
int crl::MFMWBConnectivity<TImageDimension, TComponentType>
::LoadParcellationImage()
{
	//-------------------------------------------
	// No parcellation to open
	//-------------------------------------------
	if (this->m_ROIFileName=="" ) return 0;

	//-------------------------------------------
	// Read in the parcellation label file
	//-------------------------------------------
	typename LabelImageReaderType::Pointer labelReader = LabelImageReaderType::New();
	std::cout << "   * "<<this->m_ROIFileName.c_str()<<std::endl;
	labelReader->SetFileName( this->m_ROIFileName );
	try {
		labelReader->Update();
	} catch (itk::ExceptionObject &excp)
	{
		std::cerr << "Failed to read from file " << this->m_ROIFileName << std::endl;    
		std::cerr << "Exception error message is : " << excp << std::endl;
		return EXIT_FAILURE;
	}
	this->m_LabelImage = labelReader->GetOutput();


	//-------------------------------------------
	// Determines the different label values
	//-------------------------------------------
	typedef itk::hash_set<LabelValueType> HashSetType;
	HashSetType tmpTable;
	itk::ImageRegionConstIterator<LabelImageType> it(this->m_LabelImage, this->m_LabelImage->GetRequestedRegion());
	for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
	{
		LabelValueType  v = it.Get();

		// if the input pixel is not the background
		if ( v != 0 )
		{
			// Since we use hash_set, we don t need to check if
			// the label is already in the set!
			tmpTable.insert( (unsigned long)v );
		}
	}

	//-----------------------------------
	// Create the structuring element for the morphological operator
	//-----------------------------------
	if ( m_NbVoxelDilatation>0 )
	{
		std::cout<<"- Dilate labels image (Nb Voxel="<<m_NbVoxelDilatation<<")"<<std::endl;

		typedef itk::BinaryBallStructuringElement<unsigned short,3> KernelType;
		KernelType element;
		element.SetRadius(2);
		element.CreateStructuringElement();

		//-------------------------------------------
		// Dilate
		//-------------------------------------------
		typedef crlDilateLabelImageFilter<LabelImageType,LabelImageType, KernelType> DilateFilterType;
		typename DilateFilterType::Pointer filter = DilateFilterType::New();
		filter->SetInput(m_LabelImage);
		filter->SetKernel(element);
		filter->Update();

		m_LabelImage = filter->GetOutput();
		m_LabelImage->DisconnectPipeline();

		if ( m_OutputDilatedLabelImage!="" )
		{
			typedef itk::ImageFileWriter<LabelImageType> LabelWriterType;
			typename LabelWriterType::Pointer writer=LabelWriterType::New();
			writer->SetInput(m_LabelImage);
			writer->SetFileName(m_OutputDilatedLabelImage);
			writer->Update();
		}
	}

	//-------------------------------------------
	// Compute the Index to Label table
	//-------------------------------------------
	typedef HashSetType::iterator HashSetIteratorType;
	this->m_TableIndexToLabel.set_size(tmpTable.size());
	LabelValueType maxLabel = 0;
	int i=0;
	for ( HashSetIteratorType it=tmpTable.begin(); it!=tmpTable.end(); ++it, i++ )
	{
		this->m_TableIndexToLabel[i] = (*it);

		if ( *it>maxLabel ) maxLabel=(*it);
	}

	//-------------------------------------------
	// Compute the Label to Index table
	//-------------------------------------------
	this->m_TableLabelToIndex.set_size(maxLabel+1);
	this->m_TableLabelToIndex.fill(0);
	for ( unsigned int i=0; i<this->m_TableIndexToLabel.size(); i++ )
		this->m_TableLabelToIndex[this->m_TableIndexToLabel[i]]=i;

	//-------------------------------------------
	// Allocate the connectivity matrix
	//-------------------------------------------
	int nb = this->m_TableIndexToLabel.size();
	this->m_ConnectivityMatrix.set_size(nb,nb);
	this->m_ConnectivityMatrix.fill(0);

	//-------------------------------------------
	// Allocate the color matrix
	//-------------------------------------------
	this->m_ColorMatrix.resize(nb);
	for ( int i=0; i<nb; i++ ) this->m_ColorMatrix[i].resize(nb);

	for ( int i=0; i<nb; i++ )
		for ( int j=i; j<nb; j++ )
		{
			vnl_vector_fixed<float, 3> c;
			float sum=0;
			while ( sum<0.5 ) {		// Don 't select a too dark color
				c[0] = ((float)rand() / (float)RAND_MAX);
				c[1] = ((float)rand() / (float)RAND_MAX);
				c[2] = ((float)rand() / (float)RAND_MAX);
				sum=c[0]+c[1]+c[2];
			}

			this->m_ColorMatrix[i][j] = c;
			this->m_ColorMatrix[j][i] = c;

		}

	return 0;
}

/**********************************************************************************************//**
 * \fn	template <const unsigned int TImageDimension,
 * 		class TComponentType> int crl::MFMWBConnectivity<TImageDimension,
 * 		TComponentType> ::LoadSeedingImageFromFile()
 *
 * \brief	Loads the seeding image from a file (If this file is not define do nothing. the
 * 			seeding region will be computed by thresholding the FA image) 
 *
 * \author	Benoit Scherrer
 * \date	January 2011
 *
 * \return	Error code. 
*************************************************************************************************/
template <const unsigned int TImageDimension, class TComponentType>
int crl::MFMWBConnectivity<TImageDimension, TComponentType>
::LoadSeedingImageFromFile()
{
	if ( m_SeedingImageFileName=="" ) return 0;

	//------------------------------------
	// Load the file
	//------------------------------------
	typename LabelImageReaderType::Pointer imageReader = LabelImageReaderType::New();
	std::cout << "   * "<<this->m_SeedingImageFileName.c_str()<<std::endl;
	imageReader->SetFileName( this->m_SeedingImageFileName );
	try {
		imageReader->Update();
	} catch (itk::ExceptionObject &excp)
	{
		std::cerr << "Failed to read from file " << this->m_ROIFileName << std::endl;    
		std::cerr << "Exception error message is : " << excp << std::endl;
		return EXIT_FAILURE;
	}
	this->m_ROIImage = imageReader->GetOutput();

	//------------------------------------
	// Binarize the seeding region and computes m_NbVoxelsInROI
	//------------------------------------
	this->m_NbVoxelsInROI = 0;
	typedef itk::ImageRegionIterator<LabelImageType> LabelIteratorType;
	LabelIteratorType it(this->m_ROIImage, this->m_ROIImage->GetRequestedRegion());
	for ( it.GoToBegin(); !it.IsAtEnd() ; ++it )
	{
		if ( it.Value() !=0 ) { it.Set(1);this->m_NbVoxelsInROI++; }
		else it.Set(0);
	}
	return 0;
}

/**********************************************************************************************//**
 * \fn	template <const unsigned int TImageDimension,
 * 		class TComponentType> void crl::MFMWBConnectivity<TImageDimension,
 * 		TComponentType> ::ComputeSeedingImage( )
 *
 * \brief	Calculates the seeding region based on a FA threshold.
 *
 * \author	Benoit Scherrer
 * \date	December 2010
 *
 * \return	The calculated seeding region. 
 *
 * \sa	SetThresholdToGenerateROI
 * \sa	GetThresholdToGenerateROI. 
*************************************************************************************************/
template <const unsigned int TImageDimension, class TComponentType>
void crl::MFMWBConnectivity<TImageDimension, TComponentType>
::ComputeSeedingImage( )
{
	std::cout << "- Compute the seeding region (voxels with FA>"<<m_ThresholdToGenerateROI<<")"<<std::endl;

	//------------------------------------
	// Iterators on the tensor images
	//------------------------------------
    typedef   itk::ImageRegionConstIterator< ImageType> TensorImageIteratorType;
	std::vector<TensorImageIteratorType> its;
	for ( unsigned int i=0; i<this->m_TensorImageVector.size(); i++ )
		its.push_back(TensorImageIteratorType(this->m_TensorImageVector[i], this->m_TensorImageVector[i]->GetLargestPossibleRegion() ));

	//------------------------------------
	// Allocate the seeding region
	//------------------------------------ 
	this->m_ROIImage = ROIImageType::New();
	this->m_ROIImage->SetRegions(this->m_TensorImageVector[0]->GetLargestPossibleRegion());
	this->m_ROIImage->SetOrigin(this->m_TensorImageVector[0]->GetOrigin());
	this->m_ROIImage->SetSpacing(this->m_TensorImageVector[0]->GetSpacing());
	this->m_ROIImage->SetDirection(this->m_TensorImageVector[0]->GetDirection());
	this->m_ROIImage->Allocate();

	//------------------------------------
	// Create the iterator on the seeding region
	//------------------------------------
	typedef itk::ImageRegionIterator<ROIImageType> ROIIteratorType;
	ROIIteratorType itROI(this->m_ROIImage, this->m_ROIImage->GetLargestPossibleRegion());

	//------------------------------------
	// Computes the seeding region
	//------------------------------------
	typename ImageType::PixelType tensor;
	this->m_NbVoxelsInROI = 0;
	for ( itROI.GoToBegin(); !itROI.IsAtEnd() ; ++itROI )
	{
		//--------------------------------
		// Computes the max FA for that voxel
		//--------------------------------
		double maxFA = 0;
		for ( unsigned int i=0; i<this->m_TensorImageVector.size(); i++ ) {
			tensor = its[i].Get();
			double fa = tensor.GetFractionalAnisotropy();
			if ( fa>maxFA ) 
				maxFA=fa;
		}

		//--------------------------------
		// Check if larger than the treshold
		//--------------------------------
		if ( maxFA>this->m_ThresholdToGenerateROI ) {
			itROI.Set(1);
			this->m_NbVoxelsInROI++;
		}
		else
			itROI.Set(0);

		for ( unsigned int i=0; i<its.size(); i++ )
			++(its[i]);
	}

	if (m_OutputSeedingImageFileName!="" )
	{
		std::cout << "  Write the seeding region in "<<m_OutputSeedingImageFileName<<std::endl;

		typedef itk::ImageFileWriter<ROIImageType> ROIWriterType;
		typename ROIWriterType::Pointer writer=ROIWriterType::New();
		writer->SetInput(this->m_ROIImage);
		writer->SetFileName(m_OutputSeedingImageFileName);
		writer->Update();
	}
}

/**********************************************************************************************//**
 *
 * \brief	Updates the connectivity matrix and returns the two connected rois. 
 *
 * \author	Benoit Scherrer
 * \date	December 2010
 *
 * \param	tract			The tract. 
 * \param [in,out]	roi1	The first roi. 
 * \param [in,out]	roi2	The second roi. 
 *
 * \return	. 
*************************************************************************************************/
template <const unsigned int TImageDimension, class TComponentType>
bool crl::MFMWBConnectivity<TImageDimension, TComponentType>
::UpdateConnectivityMatrix( typename TractType::Pointer tract, int& roi1, int& roi2)
{
	if ( this->m_LabelImage.GetPointer()==NULL )
		throw itk::ExceptionObject(__FILE__, __LINE__, "INTERNAL ERROR. MFMWBConnectivity::UpdateConnectivityMatrix was called but the ConnectivityMatrix is OFF.", "MFMWBConnectivity::UpdateConnectivityMatrix");

	//-------------------------------------
	// Take the two extreme physical points of the tract
	// and converts to an index for the label image
	//-------------------------------------
	typename TractType::ContinuousIndexType cindex;
	typename ImageType::PointType point1, point2;
	typename ImageType::IndexType index1, index2;

	point1 = tract->GetVertexList()->GetElement(0);
	//this->m_TensorImageVector[0]->TransformContinuousIndexToPhysicalPoint(cindex, point1);
	if ( !this->m_LabelImage->TransformPhysicalPointToIndex(point1, index1))
	{
		std::cout<<point1<< " is outside"<<std::endl;
		return false;
	}

	point2 = tract->GetVertexList()->GetElement(tract->GetVertexList()->Size()-1);
	//this->m_TensorImageVector[0]->TransformContinuousIndexToPhysicalPoint(cindex, point2);
	if ( !this->m_LabelImage->TransformPhysicalPointToIndex(point2, index2))
	{
		std::cout<<point2<< " is outside"<<std::endl;
		return false;
	}

	//-------------------------------------
	// Get the labels
	//-------------------------------------
	LabelValueType l1 = this->m_LabelImage->GetPixel(index1);
	LabelValueType l2 = this->m_LabelImage->GetPixel(index2);
	if ( l1>=this->m_TableLabelToIndex.size() || l2>=this->m_TableLabelToIndex.size() ) {
		std::cout<<"STRANGE ERROR. The label value is bigger than the label->index table size"<<std::endl;
		return false;
	}

	//-------------------------------------
	// If the two labels not null, update the connectivity matrix
	//-------------------------------------
	if ( l1!=0 && l2!=0 )
	{
		LabelValueType i1 = this->m_TableLabelToIndex[l1];
		LabelValueType i2 = this->m_TableLabelToIndex[l2];
		if ( i1<this->m_ConnectivityMatrix.cols() && i2<this->m_ConnectivityMatrix.cols() )
		{
			if ( i1!=i2 ) {
				this->m_ConnectivityMatrix[i1][i2]++;
				this->m_ConnectivityMatrix[i2][i1]++;

				roi1 = i1;
				roi2 = i2;
				return true;
			}
		}
		else
		{
			std::cout<<"STRANGE ERROR. The index associated with the label value is bigger than the number of matrix elements"<<std::endl; 
			return false;
		}
		
	}
	return false;
//ConvertContinuousIndexToNearestIndex (const ContinuousIndexType &cindex, IndexType &index)
}

/**********************************************************************************************//**
 * \brief	Compute the whole brain connectivity matrix 
 *
 * \author	Benoit Scherrer
 * \date	December 2010
 *
 * \return	Error code. 
*************************************************************************************************/
template <const unsigned int TImageDimension, class TComponentType>
int crl::MFMWBConnectivity<TImageDimension, TComponentType>
::Execute()
{
	typedef typename PixelType::ValueType ValueType;
	typedef typename PixelType::EigenVectorsMatrixType MatrixType;
	typedef typename PixelType::EigenValuesArrayType EigenValuesType;

	if ( this->m_ROIFileName!="" )	
		this->m_ConnectivityMatrixMode=true;
	else							
		this->m_ConnectivityMatrixMode=false;

	//-------------------------------------
	// Load the multi fiber model
	//-------------------------------------
	int r;
	if ( (r=this->LoadMFM())!=0 ) return r;
	if ( (r=this->LoadParcellationImage())!=0 ) return r;

	//-------------------------------------
	// If defined load the seeding image
	//-------------------------------------
	if ( m_SeedingImageFileName!="" )
	{
		if ( (r=this->LoadSeedingImageFromFile())!=0 ) return r;
	}

	//-------------------------------------
	// Else compute the seeding region based on
	// a threshold on the max FA of the MFM
	//-------------------------------------
	else this->ComputeSeedingImage();

	this->InitializeAlgo(m_NbThreads);


	//-------------------------------------------
	// If ConnectivityMatrix MODE, then
	// allocate the tracts container for each connected regions
	//-------------------------------------------
	if (m_ConnectivityMatrixMode)
	{
		m_ConnectedRegionsTracts.resize(m_ConnectivityMatrix.cols());
		for ( unsigned int i=0; i<m_ConnectedRegionsTracts.size();i++ )
		{
			for ( unsigned int j=0; j<m_ConnectedRegionsTracts.size(); j++ )
				m_ConnectedRegionsTracts[i].push_back(TractContainerType::New());
		}
	}
	//-------------------------------------------
	// If ConnectivityMatrix is OFF, only allocate a single container
	//-------------------------------------------
	else
		m_Tracts = TractContainerType::New();


	std::cout << "- Compute the whole brain connectivity..." << std::endl;


	m_Iterator = TensorImageIteratorType(this->m_TensorImageVector[0], this->m_TensorImageVector[0]->GetLargestPossibleRegion() );
	m_Iterator.GoToBegin();

#ifdef MULTITHREAD_VERSION

	//-------------------------------------------
	// Create the barrier
	//-------------------------------------------
	m_BarrierSync = itk::Barrier::New();
	m_BarrierSync->Initialize(m_NbThreads);

	m_NbVoxelProcessed = 0;

	m_MultiThreader = ThreaderType::New();
	m_MultiThreader->SetNumberOfThreads( m_NbThreads );
	m_MultiThreader->SetGlobalDefaultNumberOfThreads( m_NbThreads );
	m_MultiThreader->SetGlobalMaximumNumberOfThreads( m_NbThreads );
	m_MultiThreader->SetSingleMethod( ThreadedAnalysis, (void *)this );
	m_MultiThreader->SingleMethodExecute();	

#else
	//-------------------------------------
	// Iterate over the tensor data, and construct 
	// a tract at every location in the ROI.
	//-------------------------------------
	PixelType in;


	crl::TaskProgress generatorProgress(0, 1, this->m_NbVoxelsInROI );
	generatorProgress.SetLinePrefix("  ");
	if ( this->m_ShowProgress ) generatorProgress.InitShowProgress();

	unsigned long int tractCount = 0;
	typename TractType::Pointer tract = 0;
	typename ImageType::PointType point;
	typename ImageType::IndexType index;
	for (m_Iterator.GoToBegin(); !m_Iterator.IsAtEnd(); ++m_Iterator) 
	{
		in = m_Iterator.Get();
		index = m_Iterator.GetIndex();

		//-------------------------------------
		// Check that the current voxel is labelled of our seeding label
		//-------------------------------------
		if (this->m_ROIImage->GetPixel( index ) == 1 ) {

			if ( this->m_ShowProgress ) generatorProgress.Update();

			//-------------------------------------
			// Estimate the tract that touches this location : generate 
			// m_TractsPerVoxel tracts from this voxel
			//-------------------------------------
			for (unsigned int i = 0; i < this->m_TractsPerVoxel; i++) {

				//-------------------------------------
				// Select a random subvoxel coordinate inside the current voxel
				//-------------------------------------
				this->m_TensorImageVector[0]->TransformIndexToPhysicalPoint(index, point);
				for (unsigned int j = 0; j < ImageType::ImageDimension; j++) {
					point[j] += ( (this->m_RandGen->GetVariate()-0.5)*
						this->m_TensorImageVector[0]->GetSpacing()[j] );
				}

				//-------------------------------------
				// Estimate the tract from this subvoxel coordinate
				//-------------------------------------
				this->EstimateTract(point, tract);

				//---------------------------------------------------------
				// If required, simplify the tracts
				//---------------------------------------------------------
				if ( m_TractSimplificationMinAngle>0 )
					tract = SimplifyTract(m_TractSimplificationMinAngle, tract);

				//-------------------------------------
				// If we found a valid tract
				//-------------------------------------
				if (tract->GetVertexList()->Size() > 0) 
				{
					// Connectivity analysis
					if (m_ConnectivityMatrixMode)
					{
						int roi1, roi2;
						if ( this->UpdateConnectivityMatrix( tract, roi1, roi2 ) && this->m_OutputFileName!="" ) 
						{
							if ( roi1>=roi2 ) m_ConnectedRegionsTracts[roi1][roi2]->InsertElement( m_ConnectedRegionsTracts[roi1][roi2]->Size(), tract );
							else m_ConnectedRegionsTracts[roi2][roi1]->InsertElement( m_ConnectedRegionsTracts[roi2][roi1]->Size(), tract );
							++tractCount;
						}
					}
					else
					{
							m_Tracts->InsertElement( m_Tracts->Size(), tract );
					}				
				} 
				else {
					; // std::cout << "Tract has length 0" << std::endl;
				}
			} // for (unsigned int i = 0; i < m_TractsPerVoxel; i++)
		} // if (m_ROIImage->GetPixel( index ) == m_ROILabels.at(labelID) )

	}
#endif

	if (m_ConnectivityMatrixMode
		&&   ( m_OutputMatrixFileName!=""
		    || m_OutputFAMatrixFileName!="" 
			|| m_OutputMDMatrixFileName!="") )
		ComputeWeightedStatisticMatrixes();

	//-------------------------------------
	// Write out the tracts with an ASCII format VTK file.
	//-------------------------------------
	if ( this->m_OutputFileName!="" ) 
	{
		if (m_ConnectivityMatrixMode) 
		{
			//-----------------------------------------------------
			// Copy all the tracts in a single output tract structure
			//-----------------------------------------------------
			typename TractContainerType::Pointer		tracts = TractContainerType::New();
			typename TractColorContainerType::Pointer	tractColors = TractColorContainerType::New();

			for ( unsigned int i=0; i<m_ConnectedRegionsTracts.size(); i++ )
			{
				for ( unsigned int j=0; j<m_ConnectedRegionsTracts[i].size(); j++ )
				{
					for (unsigned int k=0; k<m_ConnectedRegionsTracts[i][j]->Size() ; k++ ) 
					{
						tracts->InsertElement( tracts->Size(), m_ConnectedRegionsTracts[i][j]->GetElement(k) );
						tractColors->InsertElement( tractColors->Size(), this->m_ColorMatrix[i][j] );
					}
				}
				m_ConnectedRegionsTracts[i].clear();				// Free memory
			}
			m_ConnectedRegionsTracts.clear();						// Free memory

			//-----------------------------------------------------
			// Save the result
			//-----------------------------------------------------
			std::cout << "- Save tracts in <"<<this->m_OutputFileName <<">."<< std::endl;
			WriteTracts(tracts, tractColors, this->m_OutputFileName);
		}
		else
		{
			std::cout << "- Save tracts in <"<<this->m_OutputFileName <<">."<< std::endl;
			typename TractColorContainerType::Pointer noColor;
			WriteTracts(m_Tracts, noColor, this->m_OutputFileName);
		}
	}

	//-------------------------------------
	// Write out the connectivity matrix
	//-------------------------------------
	if (m_ConnectivityMatrixMode && m_OutputMatrixFileName!="")
		WriteConnectivityMatrix(m_ConnectivityMatrix, this->m_OutputMatrixFileName);

	if (m_GraphCoordinatesFileName!="")
		WriteGraphCoordinates(m_GraphCoordinatesFileName);

	//tmp
	/*std::cout<<std::endl;
	for ( unsigned int i=0; i<this->m_ConnectivityMatrix.rows(); i++ )
	{
		for ( unsigned int j=0; j<this->m_ConnectivityMatrix.cols(); j++ )
			std::cout<<std::setw(5)<<this->m_ConnectivityMatrix[i][j]<<" ";
		std::cout<<std::endl;
	}*/
	return 0; // Success
}



template <const unsigned int TImageDimension, class TComponentType>
void crl::MFMWBConnectivity<TImageDimension, TComponentType>
::ThreadedTractAnalysis( void * arg )
{
	//-------------------------------------------
	// Get the thread info
	//-------------------------------------------
	typedef itk::MultiThreader::ThreadInfoStruct  ThreadInfoType;
	ThreadInfoType * infoStruct = static_cast< ThreadInfoType * >( arg );
	const unsigned int threadId = infoStruct->ThreadID;
	const unsigned int threadCount = this->m_MultiThreader->GetNumberOfThreads();

	//-------------------------------------------
	// Init the task progress for this thread
	//-------------------------------------------
	crl::TaskProgress generatorProgress(threadId, threadCount, this->m_NbVoxelsInROI );
	typename TractType::Pointer tract = 0;
	typename ImageType::PointType point;
	typename ImageType::IndexType index;

	bool keepTractsInMemory = (this->m_OutputFileName!="") 
		|| (m_OutputFAMatrixFileName!="")
		|| (m_OutputMDMatrixFileName!="");


	//-------------------------------------------
	// Loop !
	//-------------------------------------------
	while ( true )
	{
		//-------------------------------------------
		// Take the next voxel
		//-------------------------------------------
		bool nextFound=false;
		m_NextSeedingPointMutex.Lock();
		while ( !this->m_Iterator.IsAtEnd() && !nextFound )
		{
			index=m_Iterator.GetIndex() ;
			if (this->m_ROIImage->GetPixel( index ) !=0 ) nextFound=true;
			++(this->m_Iterator);
		}
		if ( nextFound ) m_NbVoxelProcessed++;
		m_NextSeedingPointMutex.Unlock();

		// Update the progress
		if ( this->m_ShowProgress ) generatorProgress.SetCurrentCounter(m_NbVoxelProcessed); //generatorProgress.Update(); //SetCurrentCounter(m_NbVoxelProcessed-1);

		//-------------------------------------------
		// If nextFound==false and we exit the loop,
		// we have reached the end of the image
		//-------------------------------------------
		if ( !nextFound )
			break;

		//-------------------------------------
		// Estimate the tract that touches this location : generate 
		// m_TractsPerVoxel tracts from this voxel
		//-------------------------------------
		for (unsigned int i = 0; i < this->m_TractsPerVoxel; i++) 
		{
			//-------------------------------------
			// Select a random subvoxel coordinate inside the current voxel
			//-------------------------------------
			this->m_TensorImageVector[0]->TransformIndexToPhysicalPoint(index, point);
			for (unsigned int j = 0; j < ImageType::ImageDimension; j++) {
				point[j] += ( (this->m_RandGen->GetVariate()-0.5)*
					this->m_TensorImageVector[0]->GetSpacing()[j] );
			}

			//-------------------------------------
			// Estimate the tract from this subvoxel coordinate
			// (return list of physical coordinates)
			//-------------------------------------
			this->EstimateTract(point, tract);

			//---------------------------------------------------------
			// If required, simplify the tracts
			//---------------------------------------------------------
			if ( this->m_TractSimplificationMinAngle>0 )
				tract = SimplifyTract(this->m_TractSimplificationMinAngle, tract);

			//-------------------------------------
			// If we found a valid tract
			//-------------------------------------
			if (tract->GetVertexList()->Size() > 0) 
			{
				m_UpdateResultsMutex.Lock();

				if (m_ConnectivityMatrixMode)
				{
					// Connectivity analysis
					int roi1, roi2;
					if ( this->UpdateConnectivityMatrix( tract, roi1, roi2 ) && keepTractsInMemory ) 
					{
						if ( roi1>=roi2 ) m_ConnectedRegionsTracts[roi1][roi2]->InsertElement( m_ConnectedRegionsTracts[roi1][roi2]->Size(), tract );
						else m_ConnectedRegionsTracts[roi2][roi1]->InsertElement( m_ConnectedRegionsTracts[roi2][roi1]->Size(), tract );
					}
				}
				else
					m_Tracts->InsertElement( m_Tracts->Size(), tract );

				m_UpdateResultsMutex.Unlock();
			} 
			else {
				; // std::cout << "Tract has length 0" << std::endl;
			}
		} // for (unsigned int i = 0; i < m_TractsPerVoxel; i++)
	} // while ( true )

	m_BarrierSync->Wait();


	if ( (int)threadId == this->m_MultiThreader->GetNumberOfThreads()-1 )
	{
		generatorProgress.Update();
		std::cout<<std::endl;
	}
}

template <const unsigned int TImageDimension, class TComponentType>
void crl::MFMWBConnectivity<TImageDimension, TComponentType>
::ComputeStatisticImage(MatrixScalarType type, typename StatisticsImageType::Pointer outImage )
{
	//--------------------------------------------
	// Prepare the iterators
	//--------------------------------------------
	std::vector<TensorImageIteratorType> itTensor;
	std::vector<StatisticsImageIteratorType> itCount;
	for ( unsigned int i=0; i<this->m_TensorImageVector.size(); i++ )
	{
		itTensor.push_back(TensorImageIteratorType(this->m_TensorImageVector[i], this->m_TensorImageVector[i]->GetLargestPossibleRegion()) );
		itCount.push_back(StatisticsImageIteratorType(m_CountConsistentFiberDirections[i], m_CountConsistentFiberDirections[i]->GetLargestPossibleRegion()));
	}
	StatisticsImageIteratorType itOut(outImage, outImage->GetLargestPossibleRegion());

	FractionsImageIteratorType itFractions;
	if ( !this->m_FractionImage.IsNull() ) itFractions = FractionsImageIteratorType(this->m_FractionImage, this->m_FractionImage->GetLargestPossibleRegion());

	std::vector<float> w(itCount.size());

	//--------------------------------------------
	// For every voxel
	//--------------------------------------------
	while ( !itOut.IsAtEnd() )
	{
		//--------------------------------------------
		// Compute the sum of the fiber counts
		//--------------------------------------------
		float sum=0;
		for ( unsigned int i=0; i<itCount.size(); i++ ) 
		{
			w[i] = itCount[i].Value();
			sum  += w[i];
		}
		if ( sum!=0 )
			for ( unsigned int i=0; i<itCount.size(); i++ )  w[i]=w[i]/sum;

		//--------------------------------------------
		// Compute the sum of the fractions
		// (not the isotropic part!)
		//--------------------------------------------
		if ( !this->m_FractionImage.IsNull() ) 
		{
			sum = 0;
			typename FractionImageType::PixelType f = itFractions.Get();
			for ( unsigned int i=0; i<itCount.size(); i++ ) 
			{
				w[i]*=f[i];
				sum += w[i];
			}
			if ( sum!=0 )
				for ( unsigned int i=0; i<itCount.size(); i++ )  w[i]=w[i]/sum;
		}

		//--------------------------------------------
		// Compute the mixture 
		//--------------------------------------------
		float value = 0;
		for ( unsigned int i=0; i<itTensor.size(); i++ )
		{
			float v;

			if ( type==FA ) v = itTensor[i].Get().GetFractionalAnisotropy();
			else v = itTensor[i].Get().GetTrace()/3;

			value += w[i] * v;
		}

		itOut.Set(value);

		//--------------------------------------------
		// Increment the iterators
		//--------------------------------------------
		if ( !this->m_FractionImage.IsNull() ) ++itFractions;
		for ( unsigned int i=0; i<itTensor.size(); i++ )
		{
			++(itTensor[i]);
			++(itCount[i]);
		}
		++itOut;
	}

	/*typedef   itk::ImageFileWriter< StatisticsImageType > StatisticsImageWriterType;
	typename StatisticsImageWriterType::Pointer writer = StatisticsImageWriterType::New();
	writer->SetInput(outImage);
	writer->SetFileName( "testFA.nrrd" );
	writer->Update();*/
}


template <const unsigned int TImageDimension, class TComponentType>
void crl::MFMWBConnectivity<TImageDimension, TComponentType>
::ComputeWeightedStatistics(typename StatisticsImageType::Pointer statImage, typename StatisticsImageType::Pointer densityImage, double& mean, double& variance )
{
	StatisticsImageIteratorType itd( densityImage, densityImage->GetLargestPossibleRegion() );
	StatisticsImageIteratorType its( statImage, statImage->GetLargestPossibleRegion());

	//----------------------------------------
	// Compute the weighted mean
	//----------------------------------------
	mean = 0.0;
	double totalweight = 0.0;
	for (itd.GoToBegin(), its.GoToBegin(); !itd.IsAtEnd(), !its.IsAtEnd();  ++its, ++itd) 
	{
		double density = itd.Value();
		double scalar = its.Value();
		mean += density * scalar;
		totalweight += density;
	}
	if (totalweight != 0.0) mean /= totalweight;

	//----------------------------------------
	// Now loop over the data again to measure the weighted variance.
	//----------------------------------------
	double weightsquared = 0.0;
	variance = 0.0;
	if ( totalweight!=0 )
	{
		for (itd.GoToBegin(), its.GoToBegin(); !itd.IsAtEnd(), !its.IsAtEnd(); ++its, ++itd) 
		{
			double density = itd.Value();
			double scalar = its.Value();
			density /= totalweight; // Normalize the weights.
			variance += density*(scalar - mean)*(scalar - mean);
			weightsquared += (density*density);
		}
		if (weightsquared != 1.0) variance /= (1.0 - weightsquared);
	}

}

template <const unsigned int TImageDimension, class TComponentType>
void crl::MFMWBConnectivity<TImageDimension, TComponentType>
::ComputeWeightedStatisticMatrixes()
{
	std::cout<<"- Compute connectivity matrices..."<<std::endl;

	//------------------------------------------
	// Allocate the tract density image
	//------------------------------------------
	typename StatisticsImageType::Pointer tractDensity = StatisticsImageType::New();
	tractDensity->SetRegions(this->m_TensorImageVector[0]->GetLargestPossibleRegion());
	tractDensity->SetOrigin(this->m_TensorImageVector[0]->GetOrigin());
	tractDensity->SetSpacing(this->m_TensorImageVector[0]->GetSpacing());
	tractDensity->SetDirection(this->m_TensorImageVector[0]->GetDirection());
	tractDensity->Allocate();

	//------------------------------------------
	// Statistic image for one diffusion parameter
	// (FA, MD, RD...)
	//------------------------------------------
	typename StatisticsImageType::Pointer statisticImage = StatisticsImageType::New();
	statisticImage->SetRegions(this->m_TensorImageVector[0]->GetLargestPossibleRegion());
	statisticImage->SetOrigin(this->m_TensorImageVector[0]->GetOrigin());
	statisticImage->SetSpacing(this->m_TensorImageVector[0]->GetSpacing());
	statisticImage->SetDirection(this->m_TensorImageVector[0]->GetDirection());
	statisticImage->Allocate();

	//------------------------------------------
	// Count the number of fiber directions corresponding to each
	// tensor main orientation 
	//------------------------------------------
	m_CountConsistentFiberDirections.clear();
	for ( unsigned int i=0;i<this->m_TensorImageVector.size(); i++ )
	{
		m_CountConsistentFiberDirections.push_back( StatisticsImageType::New() );
		m_CountConsistentFiberDirections[i]->SetRegions(this->m_TensorImageVector[0]->GetLargestPossibleRegion());
		m_CountConsistentFiberDirections[i]->SetOrigin(this->m_TensorImageVector[0]->GetOrigin());
		m_CountConsistentFiberDirections[i]->SetSpacing(this->m_TensorImageVector[0]->GetSpacing());
		m_CountConsistentFiberDirections[i]->SetDirection(this->m_TensorImageVector[0]->GetDirection());
		m_CountConsistentFiberDirections[i]->Allocate();
	}

	typename TractType::ContinuousIndexType pt;	//use ContinuousIndex but it is actually the physical position
	typename StatisticsImageType::IndexType index, previndex;  
	typename StatisticsImageType::PointType p, prevpoint;  
	typename ImageType::PixelType compareTensor;
	typename TractType::Pointer tract = 0;
	itk::Vector< TComponentType, 3> tractVector, compareVector;
	prevpoint.Fill(0);


	ConnectivityMatrixType	matrixFAMean(m_ConnectivityMatrix.rows(), m_ConnectivityMatrix.cols());
	ConnectivityMatrixType	matrixFAVar(m_ConnectivityMatrix.rows(), m_ConnectivityMatrix.cols());
	ConnectivityMatrixType	matrixMDMean(m_ConnectivityMatrix.rows(), m_ConnectivityMatrix.cols());
	ConnectivityMatrixType	matrixMDVar(m_ConnectivityMatrix.rows(), m_ConnectivityMatrix.cols());

	matrixFAMean.fill(0);
	matrixFAVar.fill(0);
	matrixMDMean.fill(0);
	matrixMDVar.fill(0);

	crl::SingleLineMultipleCout msg;

	for ( unsigned int i=0; i<m_ConnectedRegionsTracts.size(); i++ )
	{
		for ( unsigned int j=0; j<m_ConnectedRegionsTracts[i].size(); j++ )
		{
			if ( i<=j ) continue;


			std::ostringstream str;
			str<<"  - Regions "<<this->m_TableIndexToLabel[i]<<"-"<<this->m_TableIndexToLabel[j]<<"           ";
			msg.cout(str.str());

			//--------------------------------------------
			// Initialize the computation of a new density image
			// for the tracts linking 'i' to 'j'
			//--------------------------------------------
			tractDensity->FillBuffer(0);
			float nbTracts = (float)m_ConnectedRegionsTracts[i][j]->Size();
			if ( nbTracts==0 ) continue;

			for ( unsigned int k=0; k<m_CountConsistentFiberDirections.size(); k++ )
				m_CountConsistentFiberDirections[k]->FillBuffer(0);

			//--------------------------------------------
			// For each fiber bundle connecting the regions i and j
			//--------------------------------------------
			for ( unsigned int k=0; k<m_ConnectedRegionsTracts[i][j]->Size(); k++ )
			{
				previndex.Fill(-1);
				//----------------------------------
				// For each point of the tract, compute the tract
				// density and update countConsistentFiberDirections in
				// each voxel
				//----------------------------------
				tract = m_ConnectedRegionsTracts[i][j]->GetElement(k);
				for (unsigned int l=0 ; l<tract->GetVertexList()->Size() ; l++) 
				{
					pt = tract->GetVertexList()->GetElement(l);
					p[0]=pt[0]; p[1]=pt[1]; p[2]=pt[2];

					tractDensity->TransformPhysicalPointToIndex( p, index );

					//--------------------------------------------
					// Increase the density count each time the tract enters a voxel.
					//--------------------------------------------
					if (previndex != index) {
						tractDensity->SetPixel(index, tractDensity->GetPixel(index) + 1.0f/nbTracts );
						previndex = index; // Record the fact that we entered this voxel.
					
						//--------------------------------------------
						// If not the first tract point, compute the tract vector
						// and determine the most aligned tensor
						//--------------------------------------------
						if ( k>0 ) 
						{
							tractVector[0]=p[0]-prevpoint[0];
							tractVector[1]=p[1]-prevpoint[1];
							tractVector[2]=p[2]-prevpoint[2];
						
							double	projectionLength = -1.0;
							double	compareProjectionLength = 0.0;
							int		idTensor=-1;
							for ( unsigned int t=0; t<this->m_TensorImageVector.size() ; t++) 
							{
								compareTensor = this->m_TensorImageVector[t]->GetPixel(index);
								EstimatePrimaryEigenvector(compareTensor, compareVector);
								compareProjectionLength = fabs(compareVector * tractVector);
								if (compareProjectionLength > projectionLength) 
								{
									projectionLength = compareProjectionLength;
									idTensor = (int)t;
								}
							}

							m_CountConsistentFiberDirections[idTensor]->SetPixel(index, m_CountConsistentFiberDirections[idTensor]->GetPixel(index) + 1.0f );
						}
					}		
					prevpoint[0]=p[0]; prevpoint[1]=p[1]; prevpoint[2]=p[2]; 
				}
			}

			//----------------------------------
			// Write the density map if necessary
			//----------------------------------
			if ( this->m_OutputTractDensityBaseFileName!="" )
			{
				typedef   itk::ImageFileWriter< StatisticsImageType > StatisticsImageWriterType;
				typename StatisticsImageWriterType::Pointer writer = StatisticsImageWriterType::New();

				char szBuffer[1024];
				crl::FileName fn(this->m_OutputTractDensityBaseFileName);
				fn.setExtension("nrrd");
				sprintf(szBuffer, "_l%03d_l%03d", this->m_TableIndexToLabel[i], this->m_TableIndexToLabel[j] );

				writer->SetInput(tractDensity);
				writer->SetFileName( fn.getCompleteFilePath_WithSufix(szBuffer) );
				writer->Update();
			}

			//----------------------------------
			//----------------------------------
			double mean, variance;
			if ( m_OutputFAMatrixFileName!="" ) 
			{
				ComputeStatisticImage(FA, statisticImage);
				ComputeWeightedStatistics(statisticImage, tractDensity, mean, variance);
				matrixFAMean[i][j]=mean;
				matrixFAMean[j][i]=mean;
				matrixFAVar[i][j]=variance;
				matrixFAVar[j][i]=variance;
			}
			if ( m_OutputMDMatrixFileName!="" ) 
			{
				ComputeStatisticImage(MD, statisticImage);
				ComputeWeightedStatistics(statisticImage, tractDensity, mean, variance);
				matrixMDMean[i][j]=mean;
				matrixMDMean[j][i]=mean;
				matrixMDVar[i][j]=variance;
				matrixMDVar[j][i]=variance;
			}
		}
	}

	//------------------------------------------
	// Save the matrices
	//------------------------------------------
	if ( m_OutputFAMatrixFileName!="" ) 
	{
		crl::FileName fn(m_OutputFAMatrixFileName);
		WriteConnectivityMatrix(matrixFAMean, fn.getCompleteFilePath());
		WriteConnectivityMatrix(matrixFAVar, fn.getCompleteFilePath_WithSufix("_var"));		
	}
	
	if ( m_OutputMDMatrixFileName!="" ) 
	{
		crl::FileName fn(m_OutputMDMatrixFileName);
		WriteConnectivityMatrix(matrixMDMean, fn.getCompleteFilePath());
		WriteConnectivityMatrix(matrixMDVar, fn.getCompleteFilePath_WithSufix("_var"));		
	}


	//------------------------------------------
	// Create the iterator
	//------------------------------------------
	//typedef itk::ImageRegionIterator<StatisticsImageType> StatisticsImageIteratorType
	//StatisticsImageIteratorType it(tractDensity, tractDensity->GetLargestPossibleRegion() );
	//
	//for ( it.GoToBegin(); !it.IsAtEnd() ; ++it )
	//{

	//}
}


template <const unsigned int TImageDimension, class TComponentType>
int crl::MFMWBConnectivity<TImageDimension, TComponentType>
::WriteGraphCoordinates(const std::string& fileName) const
{
	//-------------------------------------
	// Open the file for writing
	//-------------------------------------
	std::ofstream ofile(fileName.c_str(), std::ios::out);
	if ( !ofile.is_open())
	{
		std::cout << "Failed writing the graph coordinates in "<< fileName.c_str()<<"."<< std::endl;
		return EXIT_FAILURE;
	}

	//-------------------------------------
	// For each matrix index
	//-------------------------------------
	for ( unsigned int i=0; i<this->m_TableIndexToLabel.size(); i++ )
	{
		int label = this->m_TableIndexToLabel[i];
		//-------------------------------------
		// Compute the barycenter for that label
		//-------------------------------------
		unsigned int nbVoxels=0;
		vnl_vector<double> barycenter(3);
		barycenter.fill(0);
		typedef itk::ImageRegionIteratorWithIndex< LabelImageType > LabelImageIteratorWithIndexType;
		LabelImageIteratorWithIndexType it(m_LabelImage, m_LabelImage->GetLargestPossibleRegion());
		for (it.GoToBegin(); !it.IsAtEnd(); ++it )
		{
			if (it.Value()==label)
			{
				typename LabelImageType::IndexType id = it.GetIndex();
				typename LabelImageType::PointType p;
				m_LabelImage->TransformIndexToPhysicalPoint(id, p);

				for ( int j=0; j<3 ; j++ )
				{
					if (nbVoxels==0) barycenter[j]=p[j];
					else barycenter[j]=(barycenter[j]*nbVoxels+p[j])/(nbVoxels+1);
				}
				nbVoxels++;
			}
		}

		//-------------------------------------
		// Output in the file
		//-------------------------------------
		ofile << i << "\t" << label << "\t" << barycenter[0]<< "\t" << barycenter[1]<< "\t" << barycenter[2] << std::endl;
	}

	//-------------------------------------
	// Close the file
	//-------------------------------------
	ofile.close();
	return 0;
}


template <const unsigned int TImageDimension, class TComponentType>
int crl::MFMWBConnectivity<TImageDimension, TComponentType>
::WriteConnectivityMatrix(const ConnectivityMatrixType &matrix, const std::string& fileName) const
{
	//-------------------------------------
	// Check...
	//-------------------------------------
	if ( matrix.rows()==0 || matrix.cols()==0 )
	{
		std::cout << "Failed writing the connectivity matrix in "<< fileName.c_str()<<". The matrix has not been computed yet."<< std::endl;
		return EXIT_FAILURE;
	}

	//-------------------------------------
	// Open the file for writing
	//-------------------------------------
	std::ofstream ofile(fileName.c_str(), std::ios::out);
	if ( !ofile.is_open())
	{
		std::cout << "Failed writing the connectivity matrix in "<< fileName.c_str()<<"."<< std::endl;
		return EXIT_FAILURE;
	}

	//-------------------------------------
	// Output the matrix
	//-------------------------------------
	ofile<<std::setprecision(8);
	for ( unsigned int i=0; i<matrix.rows(); i++ )
	{
		for ( unsigned int j=0; j<matrix.cols(); j++ )
		{
			ofile << matrix[i][j] << "\t";
		}
		ofile << std::endl;

	}

	//-------------------------------------
	// Close the file
	//-------------------------------------
	ofile.close();
	return 0;
}

#endif

