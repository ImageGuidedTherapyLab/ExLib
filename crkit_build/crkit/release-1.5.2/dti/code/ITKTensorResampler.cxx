#include "ITKTensorResampler.h"

#include <itkTransformFileReader.h>
#include <itkAffineTransform.h>

#include <itkOrientImageFilter.h>

TensorResampler::TensorResampler()
{
  m_ResampleMovingImageFilter  = ResampleFilterType::New();

  m_Transform   = TransformType::New();
  m_Transform->SetIdentity();

  m_ResampleMovingImageFilter->SetTransform( m_Transform.GetPointer() );

  interptype = 1; //linear

  m_RotateTensorFilterType    = RotateTensorFilterType::New();
}


/************************************
 *  Load Fixed Image (from a file format like nrrd or mha, mhd , using ITK reader)
 * This provides the geometry for the output.
 ***********************************/
bool TensorResampler::itkLoadFixedImage( const char * filename )
{
  if( !filename ) {
    std::cerr << " filename of fixed image is not set !" << std::endl;
    return EXIT_FAILURE;
  }
  
  itkFixedImageReaderType::Pointer  reader = itkFixedImageReaderType::New();
   reader->SetFileName( filename );
   try {
     reader->Update();
   } catch( itk::ExceptionObject & excp )
    {
      std::cerr << "Error while reading fixed image file" << filename << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
  
  m_InputFixedImage  = reader->GetOutput();
  return EXIT_SUCCESS;
}

/************************************
 *  Load Geometry Image
 * This provides the geometry for the output.
 ***********************************/
bool TensorResampler::LoadGeometryImage( const char * filename )
{
  if( !filename ) {
    std::cerr << " filename of geometry image is not set !" << std::endl;
    return EXIT_FAILURE;
  }

  itk::OrientImageFilter<GeometryImageType,GeometryImageType>::Pointer orienter = itk::OrientImageFilter<GeometryImageType,GeometryImageType>::New();

  
  GeometryImageReaderType::Pointer  reader = GeometryImageReaderType::New();
   reader->SetFileName( filename );
   try {
     reader->Update();
   } catch( itk::ExceptionObject & excp )
    {
      std::cerr << "Error while reading geometry image file" << filename << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
  
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientationToAxial();
    orienter->SetInput( reader->GetOutput() );
    orienter->Update();

    m_InputGeometryImage  = orienter->GetOutput();

  return EXIT_SUCCESS;
}


/************************************
 *  Load Moving Image (from a file format like nrrd, mhd, using ITK reader. )
 * This is the image that will be resampled.
 ***********************************/
void
TensorResampler
::itkLoadMovingImage( const char * filename )
{
  if( !filename ) {
    std::cerr << " filename of moving image is not set !" << std::endl;
    exit(1);
  }
  
  itkMovingImageReaderType::Pointer  reader = itkMovingImageReaderType::New();
   reader->SetFileName( filename );

  try {
    reader->GenerateOutputInformation();
  } catch (itk::ExceptionObject &excp)
  {
    std::cerr << "Failed to read from file " << filename << std::endl;
    std::cerr << "Exception error message is : " << excp << std::endl;
    exit(1);
  }

  itk::ImageIOBase *imageIO = reader->GetImageIO();
  if (imageIO->GetComponentType() != itk::ImageIOBase::FLOAT) {
    std::cerr << "File component type does not match expected type FLOAT."
              << std::endl;
    std::cerr << "Use crlCastSymMatToFloat to change the type of the components."
              << std::endl;
    exit(1);
  }

  try {
    reader->Update();
  } catch (itk::ExceptionObject &excp)
  {
    std::cerr << "Failed to read from file " << filename << std::endl;    
    std::cerr << "Exception error message is : " << excp << std::endl;
    exit(1);
  }
  
  m_InputMovingImage  = reader->GetOutput();
}


/************************************
 *  Load Moving Image (from 6 files, in a format like nrrd, mhd, using ITK reader. )
 * This is the image that will be resampled.
 ***********************************/
void  
TensorResampler
::itkLoadMovingImage( const char * filename1, const char * filename2, const char * filename3, const char * filename4, const char * filename5, const char * filename6 )
{

  typedef itk::OrientedImage<ComponentType, ImageDimension> ScalarImageType;
  typedef itk::ImageFileReader<ScalarImageType>  ScalarReaderType;
  
  ScalarReaderType::Pointer  reader1                  = ScalarReaderType::New();
   reader1->SetFileName( filename1 );
   reader1->Update();
  ScalarReaderType::Pointer  reader2                  = ScalarReaderType::New();
   reader2->SetFileName( filename2 );
   reader2->Update();
  ScalarReaderType::Pointer  reader3                  = ScalarReaderType::New();
   reader3->SetFileName( filename3 );
   reader3->Update();
  ScalarReaderType::Pointer  reader4                  = ScalarReaderType::New();
   reader4->SetFileName( filename4 );
   reader4->Update();
  ScalarReaderType::Pointer  reader5                  = ScalarReaderType::New();
   reader5->SetFileName( filename5 );
   reader5->Update();
  ScalarReaderType::Pointer  reader6                  = ScalarReaderType::New();
   reader6->SetFileName( filename6 );
   reader6->Update();

  m_InputMovingImage = MovingImageType::New();
  m_InputMovingImage->SetRegions( reader1->GetOutput()->GetLargestPossibleRegion() );
  m_InputMovingImage->CopyInformation( reader1->GetOutput() );
  m_InputMovingImage->Allocate();
  PixelType testvalue;
  testvalue.Fill(3);
  m_InputMovingImage->FillBuffer( testvalue );
  
  itk::ImageRegionIterator<ScalarImageType> itc1(reader1->GetOutput(), reader1->GetOutput()->GetBufferedRegion() );
  itk::ImageRegionIterator<ScalarImageType> itc2(reader2->GetOutput(), reader2->GetOutput()->GetBufferedRegion() );
  itk::ImageRegionIterator<ScalarImageType> itc3(reader3->GetOutput(), reader3->GetOutput()->GetBufferedRegion() );
  itk::ImageRegionIterator<ScalarImageType> itc4(reader4->GetOutput(), reader4->GetOutput()->GetBufferedRegion() );
  itk::ImageRegionIterator<ScalarImageType> itc5(reader5->GetOutput(), reader5->GetOutput()->GetBufferedRegion() );
  itk::ImageRegionIterator<ScalarImageType> itc6(reader6->GetOutput(), reader6->GetOutput()->GetBufferedRegion() );
  
  itk::ImageRegionIterator<MovingImageType> ittensor(m_InputMovingImage, m_InputMovingImage->GetBufferedRegion() );

  itc1.GoToBegin();
  itc2.GoToBegin();
  itc3.GoToBegin();
  itc4.GoToBegin();
  itc5.GoToBegin();
  itc6.GoToBegin();
  ittensor.GoToBegin();
  
  for (; !itc1.IsAtEnd(); ++itc1, ++itc2, ++itc3, ++itc4, ++itc5, ++itc6, ++ittensor) {
    PixelType p;
    p[0] = itc1.Get();
    p[1] = itc2.Get();
    p[2] = itc3.Get();
    p[3] = itc4.Get();
    p[4] = itc5.Get();
    p[5] = itc6.Get();
    ittensor.Set( p );
// if (p[1]>0)
//   std::cout << p << "  " << ittensor.Get() <<  std::endl;    
  }
   
}


/************************************
 *  Save Moving Image 
 *  (in the a nrrd, mha, mhd or similar format, supported by ITK writer)
 ***********************************/
void
TensorResampler
::itkWriteMovingImage( const char * filename )
{
  if( !filename )
  {
    return;
  }

  itkMovingImageWriterType::Pointer writer = itkMovingImageWriterType::New();
   writer->SetInput( m_OutputMovingImage );
   writer->SetFileName( filename );
   writer->UseCompressionOff(); // Currently off due to bug writing compressed
                                // files larger than 2GB.
   writer->Write();

}


/************************************
 *  Save Moving Image 
 *  (in the a nrrd, mha, mhd or similar format, supported by ITK writer, adds a component (ones at the first position))
 ***********************************/
void
TensorResampler
::itkWriteMovingImage7component( const char * filename )
{
  if( !filename )
  {
    return;
  }
  typedef   itk::Vector<ComponentType,7>  Vector7PixelType;
  typedef itk::OrientedImage<Vector7PixelType, ImageDimension>  OutputImageType;
  typedef itk::ImageFileWriter<OutputImageType>  OutputImageWriterType;
  
  OutputImageType::Pointer output = OutputImageType::New();
   output->SetRegions( m_OutputMovingImage->GetLargestPossibleRegion() );
   output->CopyInformation( m_OutputMovingImage );
   output->Allocate();
   
  itk::ImageRegionIterator<OutputMovingImageType> ittensor(
             m_OutputMovingImage, 
             m_OutputMovingImage->GetBufferedRegion() );
  itk::ImageRegionIterator<OutputImageType> ittensorE(output, m_OutputMovingImage->GetBufferedRegion() );
  
  for (ittensor.GoToBegin(), ittensorE.GoToBegin(); !ittensor.IsAtEnd() && !ittensorE.IsAtEnd();  ++ittensor, ++ittensorE) {
    OutputMovingImageType::PixelType pin = ittensor.Get();
    OutputImageType::PixelType pout;
    pout[0] = 1;
    for (unsigned int kk=0; kk<6; kk++)  {
      pout[kk+1] = pin[kk];
    }
    ittensorE.Set( pout );
  }
  
  OutputImageWriterType::Pointer writer = OutputImageWriterType::New();
   writer->SetInput( output );
   writer->UseCompressionOff(); // Currently off due to bug writing compressed
                                // files larger than 2GB.
   writer->SetFileName( filename );
   writer->Write();

}


/******************************************************
 *
 * Test function
 *
 ******************************************************/
void
TensorResampler
::TestVectorImage(int x, int y, int z) 
{

  MovingImageType::IndexType pixelIndex;
  pixelIndex[0] = x;
  pixelIndex[1] = y;
  pixelIndex[2] = z;

  std::cout << "got passed index init" << std::endl;

  m_RotateTensorFilterType->Update();

  OutputMovingImageType::PixelType pixelValue = m_RotateTensorFilterType->GetOutput()->GetPixel(pixelIndex);

  std::cout << "got passed getPixel" << std::endl;

  for (int i=0; i<NumberOfComponents; i++) {
    std::cout << i << ":  " << pixelValue[i] << std::endl;
  }
}

/************************************
 *
 *  Make the resampled moving image
 *
 ***********************************/
void
TensorResampler
::FlipImage( int axis )
{

  // The appropriate transform to apply should have been set already by
  // ParseParameterFile().
  typedef itk::FlipImageFilter< MovingImageType >  FilterType;

  FilterType::Pointer filter = FilterType::New();

  typedef FilterType::FlipAxesArrayType     FlipAxesArrayType;

  FlipAxesArrayType flipArray;
  flipArray.Fill( 0 );
  
  if (axis==0) {
    flipArray[0] = 1;
  } else if (axis==1) {
    flipArray[1] = 1;
  } else if (axis==2) {
    flipArray[2] = 1;
  } else  {
    std::cerr << "Wrong argument: " << axis << " is not a recognised axis name." << std::endl;
    exit(1);
  }

  filter->SetFlipAxes( flipArray );

  filter->SetInput( m_InputMovingImage );

  filter->Update();
  
  m_InputMovingImage = filter->GetOutput();

}



/************************************
 *
 *  Make the resampled moving image
 *
 ***********************************/
void
TensorResampler
::GenerateMovingImage( void )
{

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

  if (m_InputMovingImage->GetDirection() != 
      m_InputGeometryImage->GetDirection() ) {
    std::cerr << "Warning : changing the orientation of the input to match" <<
    " that of the existing image. May require a change in the " <<
    " measurement frame." << std::endl;
  }

  // The appropriate transform to apply should have been set already by
  // ParseParameterFile().
  m_ResampleMovingImageFilter->SetInput( m_InputMovingImage );

  // Size, spacing, origin and direction will be taken from the fixed image
  m_ResampleMovingImageFilter->SetOutputSpacing( m_InputGeometryImage->GetSpacing() );
  m_ResampleMovingImageFilter->SetOutputOrigin(m_InputGeometryImage->GetOrigin() );
  m_ResampleMovingImageFilter->SetSize(m_InputGeometryImage->GetLargestPossibleRegion().GetSize() );
  m_ResampleMovingImageFilter->SetOutputStartIndex(m_InputGeometryImage->GetLargestPossibleRegion().GetIndex() );
  m_ResampleMovingImageFilter->SetOutputDirection (
                                  m_InputGeometryImage->GetDirection() );

  if (interptype == 2) {
    m_ResampleMovingImageFilter->SetInterpolator(
      itk::VectorNearestNeighborInterpolateImageFunction<MovingImageType, double >::New() );
  } else if (interptype == 1) {
    m_ResampleMovingImageFilter->SetInterpolator(
      itk::VectorLinearInterpolateImageFunction<MovingImageType, double >::New() );
  } else {
    std::cerr << "No interpolation type set." << std::endl;
  }

  m_ResampleMovingImageFilter->Update();

  /*
    std::cout << "Parameters : \n" <<
      m_ResampleMovingImageFilter->GetTransform()->GetParameters() 
      << std::endl;
   */
}



/****************************************************************************
 * Rotate the tensors at each pixel according to the transformation:
 *         Tnew = Rt*T*R
 *  the rotation matrix is extracted from the affine transform parameters, 
 *  using the method described by Raimundo Sierra, in his Master's Thesis.
 *   See also the paper by Juan Ruiz-Alzola.
 * The matrix A of the affine transform is decomposed into a pure rotation part 
 * and another part. The rotation is used to reorient the rensors.
 *
 *  - using svd: A = R*W = UWVt = UVt * VWVt, where UVt is the rotation matrix.
 ****************************************************************************/
void
TensorResampler
::RotateTensors()
{
  m_RotateTensorFilterType->SetInput( m_ResampleMovingImageFilter->GetOutput() );
  
  TransformParametersType p;
  p = m_Transform->GetParameters();
  
  typedef vnl_matrix<ComponentType> vnlMatrixType;
  vnlMatrixType m_AMatrix;
  m_AMatrix.set_size( 3, 3 );
   m_AMatrix(0,0) = p[0]; m_AMatrix(0,1) = p[1]; m_AMatrix(0,2) = p[2];
   m_AMatrix(1,0) = p[3]; m_AMatrix(1,1) = p[4]; m_AMatrix(1,2) = p[5];
   m_AMatrix(2,0) = p[6]; m_AMatrix(2,1) = p[7]; m_AMatrix(2,2) = p[8];
    
   typedef vnl_svd<ComponentType>  SVDSolverType;

   SVDSolverType svd( m_AMatrix, 1e-8 );

   vnlMatrixType m_VMatrix = svd.V();
   vnlMatrixType m_UMatrix = svd.U();
   vnlMatrixType m_RMatrix = m_UMatrix*m_VMatrix.transpose();

  itk::Matrix<ComponentType,3,3> rotationmatrix;
  rotationmatrix.SetIdentity();
  rotationmatrix *= m_RMatrix;

std::cout << "Rotation matrix is: " << std::endl;
std::cout << rotationmatrix << std::endl;

  m_RotateTensorFilterType->SetRotationMatrix( rotationmatrix );

  m_RotateTensorFilterType->Update();

  m_OutputMovingImage = m_RotateTensorFilterType->GetOutput();


/*
  itk::ImageRegionIterator<OutputMovingImageType> ittensor(
          m_RotateTensorFilterType->GetOutput(), 
          m_RotateTensorFilterType->GetOutput()->GetBufferedRegion() );
  itk::ImageRegionIterator<MovingImageType> ittensor2(
          m_InputMovingImage, 
          m_InputMovingImage->GetBufferedRegion() );

  ittensor.GoToBegin();
  ittensor2.GoToBegin();
*/

/* This is just for debugging:
  for (; !ittensor.IsAtEnd();  ++ittensor, ++ittensor2) {
    PixelType p;
// turned off 2006.03.18    p = ittensor.Get();
//     if (p[0]>0 || p[1]>0)
//       std::cout << ittensor2.Get() << "\t" << ittensor.Get() << std::endl;    
  }
  */

}


void
TensorResampler
::PrintTransformParameters( void )
{
  TransformParametersType p;
  p = m_Transform->GetParameters();

  for (unsigned int i = 0; i <  m_Transform->GetNumberOfParameters(); i++) {
    std::cout << "parm i " << p[i] << std::endl;
  }

}


void
TensorResampler
::ParseParameterFile( const  char *fname )
{

  TransformParametersType parms = m_Transform->GetParameters();

  // Let's see if we can read it with an ITK transform reader
  bool readThroughITKTransformReader = true;
  itk::TransformFileReader::Pointer trsfreader;
  trsfreader = itk::TransformFileReader::New();
  trsfreader->SetFileName( fname );

  try {
    trsfreader->Update();
  } catch ( itk::ExceptionObject & excp )
  {
    readThroughITKTransformReader = false;
  }

  if (readThroughITKTransformReader) {
    typedef itk::TransformFileReader::TransformListType * TransformListType;
    TransformListType transforms = trsfreader->GetTransformList();
    itk::TransformFileReader::TransformListType::const_iterator it =
               transforms->begin();
    if (transforms->size() <= 0 || transforms->size() > 1) {
      std::cerr << "Read " << transforms->size() << 
        " transforms but want 1." << std::endl;
      exit(1);
    }
    if (strcmp((*it)->GetNameOfClass(), "AffineTransform")) {
      std::cerr << "The transform must be an AffineTransform." << std::endl;
      exit(1);
    }
    TransformType::Pointer affine_read = 
                 static_cast<TransformType*>((*it).GetPointer());
    m_Transform->SetTranslation( affine_read->GetTranslation() );
    m_Transform->SetCenter( affine_read->GetCenter() );
    m_Transform->SetMatrix( affine_read->GetMatrix() );

    return; // The transform has been successfully read.
  }

  // Continue here if the ITK transform reader didn't work out.
  

  int expectedcount = parms.Size();

  int count = 0;
  std::ifstream parmfile(fname);
  std::string line;
  std::string delims(" \t:");
  std::string::size_type idx;
  while (getline(parmfile, line)) {
    // Skip lines that are blank or start with a #
    if (line.length() == 0) continue;
    if (line[0] == '#') continue;
    idx = line.find("parameters");
    if (idx == std::string::npos) {
      std::cerr << "Failed finding parameters line in " << fname << std::endl;
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
    m_Transform->SetParameters(parms);
  } else {
    std::cerr << "Failed to find transform parameters." << std::endl;
    exit(3);
  }
  
}

void TensorResampler::SetInterpolationToLinear( )
{
  interptype = 1;
}

void TensorResampler::SetInterpolationToNearestNeighbor( )
{
  interptype = 2;
}

int TensorResampler::GetInterpolationType( void )
{
  return interptype;
}


bool
TensorResampler
::isPlatformLittleEndian(void)
{
     union mybytes {
         unsigned short int l;
         unsigned char b[sizeof(unsigned short int)];
     } p;

     p.l = 0xffU;

     return (p.b[0] != 0);
}


void
TensorResampler
::writeBigEndianFloat( std::ostream& s, float f )
{
     union {
         float f;
         unsigned char b[sizeof(float)];
     } u;


     if ( isPlatformLittleEndian() ) {
        for (int i = sizeof(f) - 1; i >= 0; i-- )
           s << u.b[i];
     } else {
        for (unsigned int i = 0; i < sizeof(f); i++)
           s << u.b[i];
     }

}


void
TensorResampler
::writeLittleEndianFloat( std::ostream& s, float f )
{
     union {
         float f;
         unsigned char b[sizeof(float)];
     } u;


     if ( !isPlatformLittleEndian() ) {
        for (int i = sizeof(f) - 1; i >= 0; i-- )
           s << u.b[i];
     } else {
        for (unsigned int i = 0; i < sizeof(f); i++)
           s << u.b[i];
     }

}
