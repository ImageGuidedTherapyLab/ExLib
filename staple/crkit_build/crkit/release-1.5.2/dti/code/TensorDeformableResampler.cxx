
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>

#include "TensorDeformableResampler.h"

#include "itkVectorResampleImageFilter.h"
#include "itkTranslationTransform.h"

TensorDeformableResampler::TensorDeformableResampler()
{
  m_WarpMovingImageFilter  = WarpFilterType::New();

  interptype = 1; //linear

  m_JacobianFilter    = JacobianFilterType::New();
  
  m_RotateTensorFilter = RotateTensorImageFilterType::New();
}


/************************************
 *  Load Fixed Image (from a file format like nrrd or mha, mhd , using ITK reader)
 * This provides the geometry for the output.
 ***********************************/
void
TensorDeformableResampler
::itkLoadDisplacementField( const char * filename )
{
  if( !filename ) {
    std::cerr << " filename of displacement field image is not set !" << std::endl;
    exit(1);
  }
  
  itkDisplacementFieldImageReaderType::Pointer  reader = itkDisplacementFieldImageReaderType::New();
   reader->SetFileName( filename );
   reader->Update();
  
  m_InputDisplacementFieldImage  = reader->GetOutput();
}


void  
TensorDeformableResampler
::itkLoadDisplacementField3ScalarImages( const char * filenamex,  const char * filenamey,  const char * filenamez )
{
  typedef itk::Image<ComponentType, ImageDimension> ScalarImageType;
  typedef itk::ImageFileReader<ScalarImageType>   ReaderType;
  ReaderType::Pointer readerX = ReaderType::New();
  ReaderType::Pointer readerY = ReaderType::New();
  ReaderType::Pointer readerZ = ReaderType::New();
   readerX->SetFileName( filenamex );
   readerY->SetFileName( filenamey );
   readerZ->SetFileName( filenamez );
   readerX->Update();
   readerY->Update();
   readerZ->Update();
  
  m_InputDisplacementFieldImage = DisplacementFieldImageType::New();
  m_InputDisplacementFieldImage->SetRegions( readerX->GetOutput()->GetLargestPossibleRegion() );
  m_InputDisplacementFieldImage->CopyInformation( readerX->GetOutput() );
  m_InputDisplacementFieldImage->Allocate();
  
  itk::ImageRegionIterator<ScalarImageType> itX(readerX->GetOutput(), readerX->GetOutput()->GetBufferedRegion() );
  itk::ImageRegionIterator<ScalarImageType> itY(readerY->GetOutput(), readerY->GetOutput()->GetBufferedRegion() );
  itk::ImageRegionIterator<ScalarImageType> itZ(readerZ->GetOutput(), readerZ->GetOutput()->GetBufferedRegion() );

  itk::ImageRegionIterator<DisplacementFieldImageType> itDisplacementField(m_InputDisplacementFieldImage, m_InputDisplacementFieldImage->GetBufferedRegion() );

  itX.GoToBegin();
  itY.GoToBegin();
  itZ.GoToBegin();
  itDisplacementField.GoToBegin();
  
  for (; !itX.IsAtEnd(); ++itX, ++itY, ++itZ, ++itDisplacementField) {
    DisplacementFieldImageType::PixelType p;
    p[0] = itX.Get();
    p[1] = itY.Get();
    p[2] = itZ.Get();
    itDisplacementField.Set( p );
  }
}

/************************************
 *  Load Moving Image (from a file format like nrrd, mhd, using ITK reader. )
 * This is the image that will be resampled.
 ***********************************/
void
TensorDeformableResampler
::itkLoadMovingImage( const char * filename )
{
  if( !filename ) {
    std::cerr << " filename of moving image is not set !" << std::endl;
    exit(1);
  }
  
  itkMovingImageReaderType::Pointer reader = itkMovingImageReaderType::New();
   reader->SetFileName( filename );
   reader->Update();
  
  m_InputMovingImage  = reader->GetOutput();
}


/************************************
 *  Load Moving Image (from 6 files, in a format like nrrd, mhd, using ITK reader. )
 * This is the image that will be resampled.
 ***********************************/
void  
TensorDeformableResampler
::itkLoadMovingImage( const char * filename1, const char * filename2, const char * filename3, const char * filename4, const char * filename5, const char * filename6 )
{

  typedef itk::Image<ComponentType, ImageDimension> ScalarImageType;
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
TensorDeformableResampler
::itkWriteMovingImage( const char * filename )
{
  if( !filename )
  {
    return;
  }
  itkMovingImageWriterType::Pointer writer = itkMovingImageWriterType::New();

   writer->SetInput( m_OutputMovingImage );
   writer->SetFileName( filename );
   writer->UseCompressionOn();

  // Transfer the measurement frame from the input data to the output data
  // Ensure the measurement frame is transferred to the output as well.
  itk::MetaDataDictionary &inDic = m_InputMovingImage->GetMetaDataDictionary();
  std::vector<std::vector<double> > msrFrame(ImageDimension);
  for (signed int j = 0; j < ImageDimension; j++) {
    msrFrame[j].resize(ImageDimension);
  }
  itk::ExposeMetaData(inDic, "NRRD_measurement frame", msrFrame);

 itk::MetaDataDictionary &outDic = m_OutputMovingImage->GetMetaDataDictionary();
  itk::EncapsulateMetaData<std::vector<std::vector<double> > >(outDic,
     std::string("NRRD_measurement frame"), msrFrame);

  try {
   writer->Update();
   writer->Write();
  } catch ( itk::ExceptionObject & excp ) {
    std::cerr << "Failed to write out the image." << std::endl;
    std::cerr << "Exception error message is : " << excp << std::endl;
  }


}


/************************************
 *  Save Moving Image 
 *  (in the a nrrd, mha, mhd or similar format, supported by ITK writer, adds a component (ones at the first position))
 ***********************************/
void
TensorDeformableResampler
::itkWriteMovingImage7component( const char * filename )
{
  if( !filename )
  {
    return;
  }
  typedef   itk::Vector<ComponentType,7>  Vector7PixelType;
  typedef itk::Image<Vector7PixelType, ImageDimension>  OutputImageType;
  typedef itk::ImageFileWriter<OutputImageType>  OutputImageWriterType;
  
  OutputImageType::Pointer output = OutputImageType::New();
   output->SetRegions( m_OutputMovingImage->GetLargestPossibleRegion() );
   output->CopyInformation( m_OutputMovingImage );
   output->Allocate();
   
  itk::ImageRegionIterator<MovingImageType> ittensor(m_OutputMovingImage, m_OutputMovingImage->GetBufferedRegion() );
  itk::ImageRegionIterator<OutputImageType> ittensorE(output, m_OutputMovingImage->GetBufferedRegion() );
  
  for (ittensor.GoToBegin(), ittensorE.GoToBegin(); !ittensor.IsAtEnd() && !ittensorE.IsAtEnd();  ++ittensor, ++ittensorE) {
    PixelType pin = ittensor.Get();
    OutputImageType::PixelType pout;
    pout[0] = 1;
    for (unsigned int kk=0; kk<6; kk++)  {
      pout[kk+1] = pin[kk];
    }
    ittensorE.Set( pout );
  }
  
  OutputImageWriterType::Pointer writer = OutputImageWriterType::New();
   writer->SetInput( output );
   writer->SetFileName( filename );
   writer->Write();

}

/************************************
 *
 ***********************************/
void
TensorDeformableResampler
::WarpMovingImage( void )
{
    
  m_WarpMovingImageFilter->SetInput( m_RotateTensorFilter->GetOutput() );

  m_WarpMovingImageFilter->SetDeformationField( m_InputDisplacementFieldImage );
  m_WarpMovingImageFilter->SetOutputSpacing( m_InputDisplacementFieldImage->GetSpacing() );
  m_WarpMovingImageFilter->SetOutputDirection( m_InputDisplacementFieldImage->GetDirection() );
  m_WarpMovingImageFilter->SetOutputOrigin(m_InputDisplacementFieldImage->GetOrigin() );

  m_WarpMovingImageFilter->Update();

  m_OutputMovingImage = m_WarpMovingImageFilter->GetOutput();

}

void
TensorDeformableResampler
::RotateTensors( float Sigma )
{
  
std::cout << "starting rotate tensors " << std::endl;

  m_JacobianFilter->SetInput( m_InputDisplacementFieldImage );
  m_JacobianFilter->SetSigma( Sigma );
  m_JacobianFilter->SetNormalizeAcrossScale( true );
  m_JacobianFilter->Update();
std::cout << " done estimate rotation matrices " << std::endl;
    
  typedef itk::VectorResampleImageFilter<
                            MovingImageType, MovingImageType >    ResampleFilterType;
  ResampleFilterType::Pointer m_ResampleMovingImageFilter  = ResampleFilterType::New();
  typedef  itk::TranslationTransform<double,ImageDimension> 	TransformType;
  TransformType::Pointer m_Transform   = TransformType::New();
  m_Transform->SetIdentity();
  m_ResampleMovingImageFilter->SetTransform( m_Transform.GetPointer() );
  m_ResampleMovingImageFilter->SetInput( m_InputMovingImage );

  m_ResampleMovingImageFilter->SetOutputSpacing( m_InputDisplacementFieldImage->GetSpacing() );
  m_ResampleMovingImageFilter->SetOutputOrigin(m_InputDisplacementFieldImage->GetOrigin() );
  m_ResampleMovingImageFilter->SetOutputDirection(m_InputDisplacementFieldImage->GetDirection() );
  m_ResampleMovingImageFilter->SetSize(m_InputDisplacementFieldImage->GetLargestPossibleRegion().GetSize() );
  m_ResampleMovingImageFilter->SetOutputStartIndex(m_InputDisplacementFieldImage->GetLargestPossibleRegion().GetIndex() );

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
std::cout << " done initial resample to get same locations as displacement field " << std::endl;
  
  m_RotateTensorFilter->SetInput1( m_ResampleMovingImageFilter->GetOutput() );
  m_RotateTensorFilter->SetInput2( m_JacobianFilter->GetOutput() );
  m_RotateTensorFilter->Update();

std::cout << " done rotate tensors " << std::endl;

}

void TensorDeformableResampler::SetInterpolationToLinear( )
{
  interptype = 1;
}

void TensorDeformableResampler::SetInterpolationToNearestNeighbor( )
{
  interptype = 2;
}

int TensorDeformableResampler::GetInterpolationType( void )
{
  return interptype;
}

