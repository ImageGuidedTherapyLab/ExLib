
#ifndef _CRLTENSORTOEIGENVECTOR_TXX_INCLUDED
#define _CRLTENSORTOEIGENVECTOR_TXX_INCLUDED 1

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <math.h>

template <const unsigned int TImageDimension, class TComponentType>
void crl::TensorToEigenvector<TImageDimension, TComponentType>
::EstimateEigenvector(PixelType &in, VectorPixelType &v)
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
    // Let's go ahead and estimate the encoding.

    for (unsigned int i = 0; i < 3; i++) {
      v[i] = V[m_EigenvectorIndex][i];
    }
}

template <const unsigned int TImageDimension, class TComponentType>
int crl::TensorToEigenvector<TImageDimension, TComponentType>
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
  typename VectorImageType::Pointer outputImage = VectorImageType::New();
  outputImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  outputImage->CopyInformation( reader->GetOutput() );
  outputImage->Allocate();

  /* Iterate over the tensor data , computing the coding. */

  ImageIterator it(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());
  PixelType in;
  VectorPixelType v;

  typename ImageType::IndexType index;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    in = it.Get();
    index = it.GetIndex();
    EstimateEigenvector(in, v);
    outputImage->SetPixel(index, v);
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

