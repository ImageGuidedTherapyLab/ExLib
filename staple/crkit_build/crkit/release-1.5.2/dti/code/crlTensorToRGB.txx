
#ifndef _CRLTENSORTORGB_TXX_INCLUDED
#define _CRLTENSORTORGB_TXX_INCLUDED 1

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkRGBPixel.h>
#include <math.h>

template <const unsigned int TImageDimension, class TComponentType>
void crl::TensorToRGB<TImageDimension, TComponentType>
::EstimateRGB(PixelType &in, ColorPixelType &color)
{
  typedef typename PixelType::ValueType ValueType;
  typedef typename PixelType::EigenVectorsMatrixType MatrixType;
  typedef typename PixelType::EigenValuesArrayType EigenValuesType;
  typedef typename itk::Vector< ValueType, 3> VectorType;

  MatrixType V;
  EigenValuesType e;
  MatrixType xin;
  MatrixType xout;
  MatrixType diag;
  VectorType v3;
  VectorType v2;
  VectorType v1;
  float minnonzero = FLT_EPSILON;   // Approximate precision of type float.
  bool isbad[3] = {false, false, false};

    in.ComputeEigenAnalysis(e, V);
    // ComputeEigenAnalysis returns the eigenvectors in the rows.
    // The vnl convention is for the eigenvalues to be ordered from 
    // the smallest to the largest.

    for (unsigned int i = 0; i < 3; i++) {
      isbad[i] = false;
    }
    for (unsigned int i = 0; i < 3; i++) {
      if (e[i] < minnonzero) {
        isbad[i] = true;
      }
    }
    if ((isbad[0] || isbad[1] || isbad[2])) {
      std::cout << "Warning: Some invalid tensors were encountered." <<
        std::endl;
      std::cout << "Use crlTensorClean to correct the tensors." << std::endl;
    }
    // Let's go ahead and estimate the color encoding.

    // Let's convert the eigenvectors to be unit magnitude.
    for (unsigned int i = 0; i < 3; i++) {
      double norm = V[i][0]*V[i][0] + V[i][1]*V[i][1]+ V[i][2]*V[i][2];
      norm = sqrt(norm);
      V[i][0] /= norm;
      V[i][1] /= norm;
      V[i][2] /= norm;
    }

    color[0] = static_cast<unsigned char>( 
                  floor(fabs(V[2][0]*255.0*in.GetFractionalAnisotropy()/
                           m_MaxFractionalAnisotropy))
               );
    color[1] = static_cast<unsigned char>(
                  floor(fabs(V[2][1]*255.0*in.GetFractionalAnisotropy()/
                           m_MaxFractionalAnisotropy))
               );
    color[2] = static_cast<unsigned char>(
                   floor(fabs(V[2][2]*255.0*in.GetFractionalAnisotropy()/
                           m_MaxFractionalAnisotropy))
               );
}

template <const unsigned int TImageDimension, class TComponentType>
int crl::TensorToRGB<TImageDimension, TComponentType>
::Execute()
{
  typedef   itk::DiffusionTensor3D<TComponentType> PixelType;
  typedef   itk::Image< PixelType, TImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  typedef   itk::ImageRegionConstIteratorWithIndex<ImageType> ImageIterator;

  typedef typename PixelType::ValueType ValueType;
  typedef typename PixelType::EigenVectorsMatrixType MatrixType;
  typedef typename PixelType::EigenValuesArrayType EigenValuesType;
  
  typename ImageReaderType::Pointer  reader = ImageReaderType::New();

  reader->SetFileName( m_InputFileName );

  try {
    reader->GenerateOutputInformation();
  } catch (itk::ExceptionObject &excp)
  {
    std::cerr << "Failed to read from file " << m_InputFileName << std::endl;    
    std::cerr << "Exception error message is : " << excp << std::endl;
    return EXIT_FAILURE;
  }

  itk::ImageIOBase *imageIO = reader->GetImageIO();
  if (imageIO->GetComponentTypeInfo() != typeid(TComponentType) ) {
    std::cerr << "File component type " <<
         imageIO->GetComponentTypeAsString(imageIO->GetComponentType()) << 
         " does not match expected type " << std::endl;
    std::cerr << "Use crlCastSymMat to change the type of the components."
              << std::endl;
    return EXIT_FAILURE;
  }

  reader->Update();

  // Allocate the new output image:
  typename ColorImageType::Pointer outputImage = ColorImageType::New();
  outputImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  outputImage->CopyInformation( reader->GetOutput() );
  outputImage->Allocate();

  /* Iterate over the tensor data , computing the color coding. */

  ImageIterator it(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());
  PixelType in;
  ColorPixelType color;

  typename ImageType::IndexType index;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    in = it.Get();
    index = it.GetIndex();
    EstimateRGB(in, color);
    outputImage->SetPixel(index, color);
  }

  typename ImageWriterType::Pointer writer = ImageWriterType::New();

  writer->SetInput ( outputImage );
  writer->SetFileName( m_OutputFileName );
  try {
    writer->Write(); // 
  } catch ( itk::ExceptionObject & excp ) {
    std::cerr << "Failed to write out the image." << std::endl;
    std::cerr << "Exception error message is : " << excp << std::endl;
  }

  exit(0); // success
}

#endif

