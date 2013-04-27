
#ifndef _TENSORRESAMPLER_INCLUDED
#define _TENSORRESAMPLER_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>

/** itk classes for image type /pixel type */
#include "itkImage.h"
#include "itkVector.h"
#include "itkMatrix.h"
#include "itkDiffusionTensor3D.h"

/** itk classes for resampling */
#include "itkWarpVectorImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVectorNearestNeighborInterpolateImageFunction.h"

#include "itkDisplacementFieldJacobianRecursiveGaussianImageFilter.h"
#include "itkRotateTensorsPixelwiseImageFilter.h"

/** the classes to read / write itk compatible formats */
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"


#include "itkCastImageFilter.h"

/** classes needed to decompose the affine transform into a rotation and rest */
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_svd.h"



/****************************************************************************
 * Rotate the tensors at each pixel according to the transformation:
 *         Tnew = Rt*T*R
 *  the rotation matrix is extracted from the Jacobian A of the displacement 
 *  at each pixel.
 *  The method is described by Raimundo Sierra, Master's thesis 2001.
 * The matrix A is decomposed into a pure rotation part 
 * and another part. The rotation is used to reorient the rensors.
 *
 *  - using svd: A = R*W = UWVt = UVt * VWVt, where UVt is the rotation matrix.
 ****************************************************************************/
class TensorDeformableResampler
{


  static int const     NumberOfComponents = 6;

  static int const      ImageDimension = 3;

  /** This component type was chosen, because the svd in (itkJacobianRecursiveGaussian) was not working well with float */
  typedef   double   ComponentType;

  typedef   itk::Vector<float,6>  VectorPixelType;
  typedef   itk::Vector<ComponentType, 3>  DisplacementPixelType;
  typedef   itk::Matrix<ComponentType, 3>  MatrixPixelType;
  
  typedef   VectorPixelType PixelType;

  int interptype; // 1 - linear , 2 - nearest neighbor

  typedef   itk::Image< PixelType, ImageDimension >  MovingImageType;
  typedef   itk::Image< PixelType, ImageDimension >  OutputMovingImageType;

  typedef   itk::Image< DisplacementPixelType, ImageDimension >  DisplacementFieldImageType;
//  typedef   itk::Image< MatrixPixelType, ImageDimension >  MatrixImageType;

  typedef itk::DisplacementFieldJacobianRecursiveGaussianImageFilter<DisplacementFieldImageType>  JacobianFilterType;
  
  typedef itk::RotateTensorsPixelwiseImageFilter<MovingImageType, JacobianFilterType::OutputImageType, OutputMovingImageType >  RotateTensorImageFilterType;
  
  typedef itk::WarpVectorImageFilter< OutputMovingImageType, 
                                OutputMovingImageType, 
                                DisplacementFieldImageType  >  WarpFilterType;


  typedef itk::CastImageFilter<MovingImageType,MovingImageType>  CopierType;
  
  /// for nrrd, mhd, mha reading and writing
  typedef   itk::ImageFileReader<DisplacementFieldImageType>      itkDisplacementFieldImageReaderType;

  typedef   itk::ImageFileReader<MovingImageType>   itkMovingImageReaderType;

  typedef   itk::ImageFileWriter<OutputMovingImageType>   itkMovingImageWriterType;
  ///-----------------------------------------------------------------


public:
  TensorDeformableResampler();
  ~TensorDeformableResampler() {}

  void  itkLoadDisplacementField( const char * filename );
  void  itkLoadDisplacementField3ScalarImages( const char * filenamex,  const char * filenamey,  const char * filenamez );
  void  itkLoadMovingImage( const char * filename );
  void  itkLoadMovingImage( const char * filename1, const char * filename2, const char * filename3, const char * filename4, const char * filename5, const char * filename6 );
  void  itkWriteMovingImage( const char * filename );
  void  itkWriteMovingImage7component( const char * filename );
  
  void  RotateTensors( float Sigma );
  void  WarpMovingImage( void );
  void  SetInterpolationToNearestNeighbor( );
  void  SetInterpolationToLinear( );
  int     GetInterpolationType( void );
  
protected:
  DisplacementFieldImageType::Pointer       m_InputDisplacementFieldImage;
  MovingImageType::Pointer                            m_InputMovingImage;
  OutputMovingImageType::Pointer                m_OutputMovingImage;

  WarpFilterType::Pointer m_WarpMovingImageFilter;
  
  JacobianFilterType::Pointer  m_JacobianFilter;

  RotateTensorImageFilterType::Pointer m_RotateTensorFilter;

};

#endif

