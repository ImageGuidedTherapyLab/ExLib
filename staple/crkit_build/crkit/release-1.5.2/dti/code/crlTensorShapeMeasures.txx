
#ifndef _CRLTENSORSHAPEMEASURES_TXX_INCLUDED
#define _CRLTENSORSHAPEMEASURES_TXX_INCLUDED 1

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkVector.h>

#include "crlTensorShapeMeasures.h"

template <const unsigned int TImageDimension, class TComponentType>
void crl::TensorShapeMeasures<TImageDimension, TComponentType>
::EstimateMeasures(PixelType &in, TComponentType &Cs,
                        TComponentType &Cp, TComponentType &Cl)
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
    V = V.GetTranspose(); // make the columns be the eigenvectors
    // The vnl convention is for the eigenvalues to be ordered from 
    // the smallest to the largest.
    // Convert from symmetric tensor class to matrix
    // Indexes are matrix(row index, column index)
    xin(0,0) = in[0];  xin(0,1) = in[1];  xin(0,2) = in[2];
    xin(1,0) = in[1];  xin(1,1) = in[3];  xin(1,2) = in[4];
    xin(2,0) = in[2];  xin(2,1) = in[4];  xin(2,2) = in[5];
    diag = V.GetInverse(); // diag <- V^{-1}
    diag = diag * xin;     // diag <- V^{-1}* x
    diag = diag * V;       // diag <- V^{-1}* x * V
    for (unsigned int i = 0; i < 3; i++) {
      isbad[i] = false;
    }
    for (signed int i = 2; i >= 0; i--) {
      if (diag(i,i) < minnonzero) {
        isbad[i] = true;
      }
    }
    if ((isbad[0] || isbad[1] || isbad[2])) {
      std::cout << "Warning: Some invalid tensors were encountered." <<
        std::endl;
      std::cout << "Use crlTensorClean to correct the tensors." << std::endl;
    }
    // Let's go ahead and estimate the three parameters.
    // Note the notation:
    // Let \lambda_1 represent the largest eigenvalue
    // Let \lambda_2 represent the middle eigenvalue
    // Let \lambda_3 represent the smallest eigenvalue
    // The ordering returned by x.ComputeEigenAnalysis is smallest to largest,
    // so e[0] == \lambda_3, e[1] == \lambda_2, e[2] == \lambda_1

    // We need to choose a norm of one of these three: 
    double norm1 = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
    double norm2 = e[0] + e[1] + e[2];
    // double norm3 = e[2];

    double norm = norm1;

    // We will adaptively correct for noise to ensure that Cs + Cp + Cl = 1.
    // We will compute a different correction for each tensor.
    // See Westin et al. Med Imag Anal 2002 page 99.
    //   We can see that:
    // enforcing the requirement that the sum of the shape values is 1
    // at every voxel is equivalent to normalizing by the trace.
    double sigma = norm2 - norm;

    Cs = 3.0 * (e[0]       ) / (norm + sigma);
    Cp = 2.0 * (e[1] - e[0]) / (norm + sigma);
    Cl = 1.0 * (e[2] - e[1]) / (norm + sigma);
    // Ca = Cl + Cp;   // A measure of anisotropy, analogous to RA.

}

template <const unsigned int TImageDimension, class TComponentType>
int crl::TensorShapeMeasures<TImageDimension, TComponentType>
::Execute()
{
  typedef   itk::SymmetricSecondRankTensor<TComponentType, TImageDimension>  
                                                      PixelType;
  typedef   itk::Image< PixelType, TImageDimension >  ImageType;
  typedef   itk::Image< TComponentType, TImageDimension >  OutImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  typedef   itk::ImageFileWriter<OutImageType> ImageWriterType;
  typedef   itk::ImageRegionConstIteratorWithIndex<ImageType> ImageIterator;

  typedef typename PixelType::ValueType ValueType;
  typedef typename PixelType::EigenVectorsMatrixType MatrixType;
  typedef typename PixelType::EigenValuesArrayType EigenValuesType;
  typedef typename itk::Vector< ValueType, 3> VectorType;
  
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
  typename OutImageType::Pointer outputImageCs = OutImageType::New();
  typename OutImageType::Pointer outputImageCp = OutImageType::New();
  typename OutImageType::Pointer outputImageCl = OutImageType::New();
  outputImageCs->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  outputImageCs->CopyInformation( reader->GetOutput() );
  outputImageCs->Allocate();

  outputImageCp->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  outputImageCp->CopyInformation( reader->GetOutput() );
  outputImageCp->Allocate();

  outputImageCl->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  outputImageCl->CopyInformation( reader->GetOutput() );
  outputImageCl->Allocate();

  /* Iterate over the tensor data , computing the shape measures. */

  ImageIterator it(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());
  PixelType in;
  TComponentType Cs;
  TComponentType Cp;
  TComponentType Cl;

  typename ImageType::IndexType index;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    in = it.Get();
    index = it.GetIndex();
    EstimateMeasures(in, Cs, Cp, Cl);
    outputImageCs->SetPixel(index, Cs);
    outputImageCl->SetPixel(index, Cl);
    outputImageCp->SetPixel(index, Cp);
  }

  typename ImageWriterType::Pointer writer = ImageWriterType::New();

  writer->SetInput ( outputImageCs );
  writer->SetFileName( m_OutputCsFileName );
  try {
    writer->Write(); // 
  } catch ( itk::ExceptionObject & excp ) {
    std::cerr << "Failed to write out the image." << std::endl;
    std::cerr << "Exception error message is : " << excp << std::endl;
  }

  writer->SetInput ( outputImageCp );
  writer->SetFileName( m_OutputCpFileName );
  try {
    writer->Write(); // 
  } catch ( itk::ExceptionObject & excp ) {
    std::cerr << "Failed to write out the image." << std::endl;
    std::cerr << "Exception error message is : " << excp << std::endl;
  }

  writer->SetInput ( outputImageCl );
  writer->SetFileName( m_OutputClFileName );
  try {
    writer->Write(); // 
  } catch ( itk::ExceptionObject & excp ) {
    std::cerr << "Failed to write out the image." << std::endl;
    std::cerr << "Exception error message is : " << excp << std::endl;
  }

  exit(0); // success
}

#endif

