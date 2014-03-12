
#ifndef _TENSORRESAMPLER_INCLUDED
#define _TENSORRESAMPLER_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>

/** itk classes for image type /pixel type */
#include "itkImage.h"
#include "itkOrientedImage.h"
#include "itkVector.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkDiffusionTensor3D.h"


/** itk classes for resampling */
#include "itkVectorResampleImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVectorNearestNeighborInterpolateImageFunction.h"
#include "itkAffineTransform.h"
#include "itkRotateTensorImageFilter.h"

#include "itkFlipImageFilter.h"

/** the classes to read / write itk compatible formats */
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"


#include "itkCastImageFilter.h"
#include "itkMetaDataObject.h"

/** classes needed to decompose the affine transform into a rotation and rest */
#include "vnl/vnl_matrix.h"
//#include "vnl/vnl_transpose.h"
#include "vnl/algo/vnl_svd.h"




class TensorResampler
{


  static int const NumberOfComponents = 6;
  static int const ImageDimension = 3;
  typedef   float ComponentType;

  typedef   itk::Vector<ComponentType,6>  VectorPixelType;
  //typedef   itk::SymmetricSecondRankTensor<ComponentType,3>  TensorPixelType;
  typedef   itk::DiffusionTensor3D<ComponentType>  TensorPixelType;
  
  typedef   VectorPixelType PixelType;


  typedef  itk::AffineTransform<double,ImageDimension> 	TransformType;

  int interptype; // 1 - linear , 2 - nearest neighbor

  typedef   itk::OrientedImage< PixelType, ImageDimension >  MovingImageType;
  typedef   itk::OrientedImage< TensorPixelType, ImageDimension >  
                                             OutputMovingImageType;
  typedef   itk::OrientedImage< PixelType, ImageDimension >  FixedImageType;
  typedef   itk::OrientedImage< ComponentType, ImageDimension >  GeometryImageType;


  typedef   TransformType::ParametersType TransformParametersType;

  typedef itk::VectorResampleImageFilter<
                            MovingImageType, MovingImageType > 
               ResampleFilterType;

  typedef itk::RotateTensorImageFilter<MovingImageType, OutputMovingImageType> RotateTensorFilterType;


  typedef itk::CastImageFilter<MovingImageType,MovingImageType>  CopierType;
  
  /// for nrrd, mhd, mha reading and writing
  typedef   itk::ImageFileReader<GeometryImageType>  GeometryImageReaderType;
  typedef   itk::ImageFileReader<FixedImageType>      itkFixedImageReaderType;

  typedef   itk::ImageFileReader<MovingImageType>   itkMovingImageReaderType;

  typedef   itk::ImageFileWriter<OutputMovingImageType>   itkMovingImageWriterType;
  ///-----------------------------------------------------------------


public:
  TensorResampler();
  ~TensorResampler() {}

  bool LoadGeometryImage( const char * filename );
  bool itkLoadFixedImage( const char * filename );
  void  itkLoadMovingImage( const char * filename );
  void  itkLoadMovingImage( const char * filename1, const char * filename2, const char * filename3, const char * filename4, const char * filename5, const char * filename6 );
  void  itkWriteMovingImage( const char * filename );
  void  itkWriteMovingImage7component( const char * filename );
  
  void  vtkLoadFixedImage( const char * filename );
  void  vtkLoadMovingImage( const char * filename );
  void  vtkWriteMovingImage( const char * filename );
  void  vtkWriteMovingImageASCII( const char * filename );
  void  vtkWriteMovingImageBinary( const char * filename );
  
  void  GenerateMovingImage( void );
  void  PrintTransformParameters( void );
  void  ParseParameterFile( const char *fname );
  void  SetInterpolationToNearestNeighbor( );
  void  SetInterpolationToLinear( );
  int     GetInterpolationType( void );
  void  RotateTensors( void );
  void  FlipImage( int axis );
  void  TestVectorImage(int x, int y, int z);

  void TransferMetaDataDictionary() {
    m_OutputMovingImage->SetMetaDataDictionary(
          m_InputMovingImage->GetMetaDataDictionary() );
  }
  
  bool isPlatformLittleEndian(void);
  void writeBigEndianFloat( std::ostream& s, float f );
  void writeLittleEndianFloat( std::ostream& s, float f );
  
protected:

  GeometryImageType::Pointer m_InputGeometryImage;
  FixedImageType::Pointer m_InputFixedImage;
  MovingImageType::Pointer m_InputMovingImage;
  OutputMovingImageType::Pointer m_OutputMovingImage;

  ResampleFilterType::Pointer m_ResampleMovingImageFilter;
  TransformType::Pointer m_Transform;
  RotateTensorFilterType::Pointer m_RotateTensorFilterType;

};

#endif

