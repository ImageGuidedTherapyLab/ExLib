
#ifndef _TENSORCHECK_TXX_INCLUDED
#define _TENSORCHECK_TXX_INCLUDED 1

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkVector.h>

#include "crlTensorClean.h"

template <const unsigned int TImageDimension, class TComponentType>
bool TensorCheck<TImageDimension, TComponentType>
::CorrectTensor(PixelType &in, PixelType &out)
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
  float twiceminnonzero = 2.*FLT_EPSILON;   // Approximate precision of type float.
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
    if (!(isbad[0] || isbad[1] || isbad[2])) {
      out = in;
      return true; // No more processing to be done
    }

    // We only get to here if there is a bad eigenvalue and we need to fix it.  
    // Correct the tensors by projecting.
    if (isbad[2]) { // Largest eigenvalue is bad -> all are bad
      // Reset the tensor to a small isotropic tensor
      diag(2,2) = twiceminnonzero; 
      diag(0,2) = 0.0; diag(1,2) = 0.0;
      diag(1,1) = twiceminnonzero;
      diag(0,1) = 0.0; diag(2,1) = 0.0;
      diag(0,0) = twiceminnonzero;
      diag(1,0) = 0.0; diag(2,0) = 0.0;
      v3.Fill(0);
      v3[0] = 1.;
      v2.Fill(0);
      v2[1] = 1.;
      v1.Fill(0);
      v1[2] = 1.;
    } else if (isbad[1]) { // Two smallest evals are bad.
      // Reset the small eigenvalues and orthogonalize the tensor
      diag(0,2) = 0.0; diag(1,2) = 0.0;
      diag(1,1) = twiceminnonzero;
      diag(0,1) = 0.0; diag(2,1) = 0.0;
      diag(0,0) = twiceminnonzero;
      diag(1,0) = 0.0; diag(2,0) = 0.0;
      // Assign the eigenvectors to separate vector data structures.
      for (unsigned int i = 0; i < 3; i++) {
        v3[i] = V(i,2);
        v2[i] = V(i,1);
        v1[i] = V(i,0);
      }
      // Compute a new v2, perpendicular to v3.
      v2 = v2 - v3*(v2*v3);
      v2.Normalize();
      // Now compute a new v1, perpendicular to v3 and to v2.
      // Compute the component 0 of the cross product
      v1[0] = v2[1]*v3[2] - v2[2]*v3[1];
      // Compute the component 1 of the cross product
      v1[1] = v2[2]*v3[0] - v2[0]*v3[2];
      // Compute the component 2 of the cross product
      v1[2] = v2[0]*v3[1] - v2[1]*v3[0];
      v1.Normalize();
    } else if (isbad[0]) { // One smallest eval is bad.
      diag(0,0) = twiceminnonzero;
      diag(1,0) = 0.0; diag(2,0) = 0.0;
      diag(0,1) = 0.0; diag(2,1) = 0.0;
      diag(0,2) = 0.0; diag(1,2) = 0.0;
      // Assign the eigenvectors to separate vector data structures.
      for (unsigned int i = 0; i < 3; i++) {
        v3[i] = V(i,2);
        v2[i] = V(i,1);
        v1[i] = V(i,0);
      }
      // Now compute a new v1, perpendicular to v3 and to v2.
      // Compute the component 0 of the cross product
      v1[0] = v2[1]*v3[2] - v2[2]*v3[1];
      // Compute the component 1 of the cross product
      v1[1] = v2[2]*v3[0] - v2[0]*v3[2];
      // Compute the component 2 of the cross product
      v1[2] = v2[0]*v3[1] - v2[1]*v3[0];
      v1.Normalize();
    }

    /* Copy the new data back into the V matrix */
    for (unsigned int i = 0; i < 3; i++) {
      V(i,2) = v3[i];
      V(i,1) = v2[i];
      V(i,0) = v1[i];
    }
    xout = V * diag;
    xout = xout * V.GetInverse();

    // Convert from matrix type to symmetric tensor output
    // The output should already be perfectly symmetrical, but due to 
    // small numerical differences, it might not be.
    out[0] = xout(0,0); out[1] = 0.5*(xout(0,1)+xout(1,0));  
    out[2] = 0.5*(xout(0,2)+xout(2,0));
    out[3] = xout(1,1); out[4] = 0.5*(xout(1,2)+xout(2,1));
    out[5] = xout(2,2);

    return false; // We modified the tensor
}


template <const unsigned int TImageDimension, class TComponentType>
int TensorCheck<TImageDimension, TComponentType>
::Execute()
{
  typedef   itk::SymmetricSecondRankTensor<TComponentType, TImageDimension>  
                                                      PixelType;
  typedef   itk::Image< PixelType, TImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  typedef   itk::ImageFileWriter<ImageType> ImageWriterType;
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
  typename ImageType::Pointer outputImage = ImageType::New();
  outputImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  outputImage->CopyInformation( reader->GetOutput() );
  outputImage->Allocate();

  /* Iterate over the tensor data here, checking the quality of the tensors. */

  ImageIterator it(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());
  PixelType out;
  PixelType in;
  unsigned long goodcount = 0;
  unsigned long badcount = 0;
  typename ImageType::IndexType index;
  bool correct;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    in = it.Get();
    index = it.GetIndex();

    correct = CorrectTensor(in, out);
    if (correct) {
      goodcount++;
    } else {
      badcount++;
    }
    outputImage->SetPixel(index, out);

  }

  std::cout << "Total tensor count : " << (goodcount + badcount) << std::endl;
  std::cout << "Good tensors : " << (goodcount) << std::endl;
  std::cout << "Bad tensors : " << (badcount) << std::endl;

  // Transfer the measurement frame from the input data to the output data

  // Ensure the measurement frame is transferred to the output as well.
  itk::MetaDataDictionary &inDic = reader->GetOutput()->GetMetaDataDictionary();
  std::vector<std::vector<double> > msrFrame(TImageDimension);
  for (unsigned int j = 0; j < TImageDimension; j++) {
    msrFrame[j].resize(TImageDimension);
  }
  itk::ExposeMetaData(inDic, "NRRD_measurement frame", msrFrame);

  itk::MetaDataDictionary &outDic = 
                outputImage->GetMetaDataDictionary();
  itk::EncapsulateMetaData<std::vector<std::vector<double> > >(outDic,
     std::string("NRRD_measurement frame"), msrFrame);

  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput ( outputImage );
  writer->SetFileName( m_OutputFileName );

  writer->SetUseCompression( m_UseCompression ); 
  // Really slow writing compressed tensors. Discard previous comment. 
  // Believe me, compared to writing and reading time, the compression is fast.
  // Compression fails for files for which the internal buffer needs to be 
  // larger than 2GB, but uncompressed writing of such files is fine. 
  // Therefore, by default we have compression writing disabled.

  try {
    writer->Write(); // 
  } catch ( itk::ExceptionObject & excp ) {
    std::cerr << "Failed to write out the image." << std::endl;
    std::cerr << "Exception error message is : " << excp << std::endl;
  }

  exit(0); // success
}

#endif

