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

#ifndef _CRL_MFMWBCONNECTIVITY_H_INCLUDED 
#define _CRL_MFMWBCONNECTIVITY_H_INCLUDED 1

#include <itk_hash_set.h>
#include <itkMultiThreader.h>
#include <itkBarrier.h>

#include "crlMFMTractGenerator.h"

#ifdef _MSC_VER
  #pragma warning( disable: 4996 )
#endif

namespace crl {

template <const unsigned int TImageDimension, class TComponentType >
class ITK_EXPORT MFMWBConnectivity : 
                    public crl::MFMTractGenerator<TImageDimension, TComponentType>
{
  public: 
  /** Standard class typedefs. */
  typedef MFMWBConnectivity Self;
  typedef MFMTractGenerator<TImageDimension, TComponentType> Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef vnl_matrix<double> ConnectivityMatrixType;
  typedef std::vector< std::vector< vnl_vector_fixed<float,3> > > ColorMatrixType;

  typedef  unsigned short LabelValueType;
  typedef   itk::OrientedImage< LabelValueType, TImageDimension> LabelImageType;
  typedef   itk::ImageFileReader<LabelImageType> LabelImageReaderType;

  typedef   itk::Image< float, TImageDimension> StatisticsImageType;
  typedef   typename Superclass::FractionImageType FractionImageType;

  typedef	vnl_vector<LabelValueType> TableConversionType;

  typedef typename Superclass::TractType TractType;
  typedef typename	Superclass::ROIImageType ROIImageType;
  typedef typename	Superclass::ImageType ImageType;
  typedef typename	Superclass::PixelType PixelType;
  typedef typename	Superclass::TractContainerType TractContainerType;
  typedef typename  Superclass::TractColorContainerType TractColorContainerType;

  typedef   itk::ImageRegionConstIterator< ImageType> TensorImageIteratorType;
  typedef   itk::ImageRegionIterator< StatisticsImageType> StatisticsImageIteratorType;
  typedef   itk::ImageRegionIterator< FractionImageType > FractionsImageIteratorType;
  typedef   itk::ImageRegionIterator< LabelImageType > LabelImageIteratorType;

  typedef itk::MultiThreader  ThreaderType;

  typedef enum {
	FA,
	MD
  } MatrixScalarValue;
  typedef unsigned int MatrixScalarType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  itkSetMacro( ThresholdToGenerateROI, double);
  itkGetMacro( ThresholdToGenerateROI, double);

  itkSetMacro( NbThreads, int );
  itkGetMacro( NbThreads, int );

  itkSetMacro( OutputMatrixFileName, std::string);
  itkGetMacro( OutputMatrixFileName, std::string);

  itkSetMacro( OutputFAMatrixFileName, std::string);
  itkGetMacro( OutputFAMatrixFileName, std::string);

  itkSetMacro( OutputMDMatrixFileName, std::string);
  itkGetMacro( OutputMDMatrixFileName, std::string);

  itkSetMacro( OutputTractDensityBaseFileName, std::string);
  itkGetMacro( OutputTractDensityBaseFileName, std::string);

  itkSetMacro( SeedingImageFileName, std::string);
  itkGetMacro( SeedingImageFileName, std::string);

  itkSetMacro( OutputSeedingImageFileName, std::string);
  itkGetMacro( OutputSeedingImageFileName, std::string);

  itkSetMacro( GraphCoordinatesFileName, std::string);
  itkGetMacro( GraphCoordinatesFileName, std::string);

  itkSetMacro( OutputDilatedLabelImage, std::string);
  itkGetMacro( OutputDilatedLabelImage, std::string);

  itkSetMacro( NbVoxelDilatation, unsigned int);
  itkGetMacro( NbVoxelDilatation, unsigned int);

 
  MFMWBConnectivity() {
	  m_ThresholdToGenerateROI = 0.7;
	  m_NbThreads = 1;
  };

  ~MFMWBConnectivity() {};

  int Execute();

	void ThreadedTractAnalysis( void * arg ) ;

protected:
	int LoadParcellationImage();

	int LoadSeedingImageFromFile();
	void ComputeSeedingImage();

	bool UpdateConnectivityMatrix( typename TractType::Pointer tract, int& roi1, int& roi2);

	void ComputeWeightedStatisticMatrixes();
	void ComputeStatisticImage(MatrixScalarType type, typename StatisticsImageType::Pointer outImage );
	void ComputeWeightedStatistics(typename StatisticsImageType::Pointer statImage, typename StatisticsImageType::Pointer densityImage, double& mean, double& variance );

	int WriteConnectivityMatrix(const ConnectivityMatrixType &matrix, const std::string& fileName) const;
	int WriteGraphCoordinates(const std::string& fileName) const;

protected:
	typename LabelImageType::Pointer	m_LabelImage;
	ConnectivityMatrixType				m_ConnectivityMatrix;
	ColorMatrixType						m_ColorMatrix;


	TableConversionType					m_TableIndexToLabel;
	TableConversionType					m_TableLabelToIndex;

	double								m_ThresholdToGenerateROI;
	std::string							m_OutputMatrixFileName;
	std::string							m_OutputFAMatrixFileName;
	std::string							m_OutputMDMatrixFileName;
	std::string							m_SeedingImageFileName;
	std::string							m_OutputSeedingImageFileName;
	std::string							m_OutputTractDensityBaseFileName;
	std::string							m_GraphCoordinatesFileName;

	std::string							m_OutputDilatedLabelImage;
	unsigned int						m_NbVoxelDilatation;


	int									m_NbThreads;

	std::vector<std::vector<typename TractContainerType::Pointer> > m_ConnectedRegionsTracts;	// If ConnectivityMatrix Mode
	typename TractContainerType::Pointer							m_Tracts ;					// If NOT ConnectivityMatrix Mode

	std::vector<typename StatisticsImageType::Pointer> m_CountConsistentFiberDirections;
	TensorImageIteratorType						m_Iterator;

	// For multi thread
	ThreaderType::Pointer				m_MultiThreader;
	itk::SimpleFastMutexLock			m_NextSeedingPointMutex;
	itk::SimpleFastMutexLock			m_UpdateResultsMutex;
	itk::Barrier::Pointer				m_BarrierSync;
	unsigned long int					m_NbVoxelProcessed;

	bool								m_ConnectivityMatrixMode;
};

} // end namespace crl

#endif

