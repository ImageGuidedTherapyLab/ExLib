
#ifndef _CRLTRACTGENERATOR_TXX_INCLUDED
#define _CRLTRACTGENERATOR_TXX_INCLUDED 1

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkContinuousIndex.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "TensorLogFunctor.h"
#include "TensorExpFunctor.h"
#include "crlTensorLinearInterpolateImageFunction.h"

template <const unsigned int TImageDimension, class TComponentType>
void crl::TractGenerator<TImageDimension, TComponentType>
::EstimatePrimaryEigenvector(typename ImageType::PixelType &in, 
                               itk::Vector< TComponentType, 3> &evec)
{
  typedef typename ImageType::PixelType::ValueType ValueType;
  typedef typename ImageType::PixelType::EigenVectorsMatrixType MatrixType;
  typedef typename ImageType::PixelType::EigenValuesArrayType EigenValuesType;

  MatrixType V;
  EigenValuesType e;
    // The vnl convention is for the eigenvalues to be ordered from
    // the smallest to the largest, and for the eigenvectors to be in the rows.
    in.ComputeEigenAnalysis(e, V);

  evec[0] = V[0][2];
  evec[1] = V[1][2];
  evec[2] = V[2][2];
}

template <const unsigned int TImageDimension, class TComponentType>
void crl::TractGenerator<TImageDimension, TComponentType>
::EstimateTract(typename ImageType::PointType &point, 
    typename TractType::Pointer &tract)
{
  // Initialize the tract container.
  tract = TractType::New();
  
  // Store the forward and reverse directions of the tract
  typename TractType::Pointer forwardTract = TractType::New();
  typename TractType::Pointer backwardTract = TractType::New();

  // Start from point, find first part of the tract, 
  //    find second part of the tract, merge the parts and return.
  // std::cout << "Estimating a tract starting from " << point <<  "." <<std::endl;

  // Do Euler integration with a fixed step size.

  // This seems to be a large step size.
  /*
  double stepSize = m_ROIImage->GetSpacing()[0]*m_ROIImage->GetSpacing()[0] +
                    m_ROIImage->GetSpacing()[1]*m_ROIImage->GetSpacing()[1] +
                    m_ROIImage->GetSpacing()[2]*m_ROIImage->GetSpacing()[2];
    stepSize = sqrt(stepSize);
  */
  // It might be better to aim to use about half the minimum size of any edge of
  // any voxel, so we get two or more steps inside each voxel.
  double stepSize = fmin(fmin(m_ROIImage->GetSpacing()[0], 
                              m_ROIImage->GetSpacing()[1]), 
                           m_ROIImage->GetSpacing()[2]);
  stepSize /= static_cast<double>(m_StepsPerVoxel);

  // std::cout << "Point to point step size is " << stepSize << std::endl;

  unsigned int stepCount = 0;
  typename ImageType::PointType newpoint = point;
  typename ImageType::PixelType tensor;
  typename ImageType::PixelType logtensor;
  itk::Vector< TComponentType, 3> currentEigenvector;
  itk::Vector< TComponentType, 3> previousEigenvector;
  itk::Vector< TComponentType, 3> tendDeflection;
  itk::Vector< TComponentType, 3> tendDeflection2;
  itk::Vector< TComponentType, 3> initialDirection;

  previousEigenvector.Fill(0.0);
  tendDeflection.Fill(0.0);
  tendDeflection2.Fill(0.0);
  initialDirection.Fill(0.0);
  itk::Functor::TensorExp<typename ImageType::PixelType, 
                          typename ImageType::PixelType> exponentiator;

  typename TractType::ContinuousIndexType cindex;
  // Used for NN interpolation: typename ImageType::IndexType dindex;
  bool isInside;
  bool forwardDirection = true;
  double projectionLength = 0.0;
  double previousAngleCriterion = 0.0;
  double fractionalAnisotropyM = 1.0;
  // We first step along in one direction,
  // then restart along the opposite direction.
  // First we work on the forward direction.
  while (stepCount < m_MaxTractLength) {
    // std::cout << "Current point is " << newpoint << std::endl;
    isInside = m_TensorImage->TransformPhysicalPointToContinuousIndex(
       newpoint, cindex
       );
    if (!isInside) {
      break;
    };
    // We could compute the tensor and the logtensor to see how big a difference
    // we get. We exponentiate so we can easily get the FA of the new tensor.
    logtensor = m_TensorInterpolator->EvaluateAtContinuousIndex(cindex);
    tensor = exponentiator(logtensor);
    
    fractionalAnisotropyM = 
              fractionalAnisotropyM * m_FractionalAnisotropyMomentum + 
        (1 - m_FractionalAnisotropyMomentum)*tensor.GetFractionalAnisotropy();
    if (fractionalAnisotropyM < m_MinimumFA) {
      // std::cout << "Terminating tract due to low FA. " << std::endl;
      break;
    }

    // Compute the largest eigenvector of the tensor
    EstimatePrimaryEigenvector(tensor, currentEigenvector);
    currentEigenvector  /= (currentEigenvector.GetNorm());

    // Select the direction of the step, and step it.
    // Update the point we are working on.
    if ((stepCount == 0) && (forwardDirection)) {
      // Save the first step direction, so later we can go the other way.
      initialDirection = currentEigenvector;
      // To initialize the tensor deflection calculation:
      previousEigenvector = currentEigenvector; 
    } else {
      // If the eigenvector is opposite to our current direction, reverse it,
      // since its orientation is arbitrary.
      // This could be generalized to be a test on angle.
        projectionLength = (currentEigenvector * previousEigenvector);
        if ((projectionLength) < 0.0) {
          currentEigenvector *= -1.0;
        }
        // Utilize the maximum angle test, based on a low pass filtered
        // set of angle changes.
        previousAngleCriterion = 
                    previousAngleCriterion*m_AngleThresholdMomentum +
                    (1.0 - m_AngleThresholdMomentum)*fabs(projectionLength);
        if (previousAngleCriterion < m_MaximumAngle) {
          // std::cout << "Terminating tract since step exceeds maximum angle." << std::endl;
          break;
        }
    }

    // The current point is inside the image and has acceptable FA, so 
    // we will add it to the tract. The tract stores the continuous index of
    // each of the points that we add.
    if (forwardDirection) {
      forwardTract->AddVertex(cindex);
    } else {
      backwardTract->AddVertex(cindex);
    }
      
    // Record our direction.
    // Compute the new direction indicated by Tensor Deflection.
    // TEND is based on the incoming direction, not the current eigenvector.
    if (m_TendFraction > 0.0) {
      // alpha = 1.0
      tendDeflection.Fill(0.0);
      for (unsigned int i = 0; i < 3; i++) {
        tendDeflection[0] += tensor(0,i) * previousEigenvector[i];
        tendDeflection[1] += tensor(1,i) * previousEigenvector[i];
        tendDeflection[2] += tensor(2,i) * previousEigenvector[i];
      }
      tendDeflection /= (tendDeflection.GetNorm());
      // alpha = 2.0
      tendDeflection2.Fill(0.0);
      for (unsigned int i = 0; i < 3; i++) {
        tendDeflection2[0] += tensor(0,i) * tendDeflection[i];
        tendDeflection2[1] += tensor(1,i) * tendDeflection[i];
        tendDeflection2[2] += tensor(2,i) * tendDeflection[i];
      }
      tendDeflection = tendDeflection2 / tendDeflection2.GetNorm();
    }
    previousEigenvector = (previousEigenvector*m_DirectionMomentum +
      (tendDeflection*m_TendFraction + currentEigenvector*(1.0-m_TendFraction))
               *(1.0 - m_DirectionMomentum) );
    //Add it to the current position after scaling by the step size.
    // Utilize direction momentum in computing the new step.
    newpoint += ( previousEigenvector * stepSize );

    stepCount++;
  }

  // std::cout << "Working in the backward direction." << std::endl;

  // Let's do it again in the opposite direction.
  // The first point to consider is one step away from the initial start.
  forwardDirection = false;
  // If we could not start a tract in the forward direction,
  // and the initialDirection was not set, this will again stop at the
  // very first voxel, and the total tract length will be zero.
  previousEigenvector = initialDirection * -1.0;
  previousAngleCriterion = 0.0;
  fractionalAnisotropyM = 1.0;
  newpoint = point + previousEigenvector*stepSize;
  while (stepCount < m_MaxTractLength) {
    // std::cout << "Current point is " << newpoint << std::endl;
    isInside = m_TensorImage->TransformPhysicalPointToContinuousIndex(
       newpoint, cindex
       );
    if (!isInside) {
      break;
    };
    // We could compute the tensor and the logtensor to see how big a difference
    // we get. We exponentiate so we can easily get the FA of the new tensor.
    logtensor = m_TensorInterpolator->EvaluateAtContinuousIndex(cindex);
    tensor = exponentiator(logtensor);

    fractionalAnisotropyM = 
              fractionalAnisotropyM * m_FractionalAnisotropyMomentum + 
        (1 - m_FractionalAnisotropyMomentum)*tensor.GetFractionalAnisotropy();
    if (fractionalAnisotropyM < m_MinimumFA) {
      // std::cout << "Terminating tract due to low FA. " << std::endl;
      break;
    }

    // Select the direction of the step, and step it.
    // Update the point we are working on.
    // Compute the largest eigenvector of the tensor
    EstimatePrimaryEigenvector(tensor, currentEigenvector);
    currentEigenvector  /= (currentEigenvector.GetNorm());
    if ((stepCount == 0) && (forwardDirection)) {
      // Save the first step direction, so later we can go the other way.
      initialDirection = currentEigenvector;
      // To initialize the tensor deflection calculation:
      previousEigenvector = currentEigenvector; 
    } else {
      // If the eigenvector is opposite to our current direction, reverse it.
        projectionLength = (currentEigenvector * previousEigenvector);
        if ((projectionLength) < 0.0) {
          currentEigenvector *= -1.0;
        }
        // Utilize the maximum angle test, based on a low pass filtered
        // set of angle changes.
        previousAngleCriterion = 
                    previousAngleCriterion*m_AngleThresholdMomentum +
                    (1.0 - m_AngleThresholdMomentum)*fabs(projectionLength);
        if (previousAngleCriterion < m_MaximumAngle) {
          // std::cout << "Terminating tract since step exceeds maximum angle." << std::endl;
          break;
        }
    }

    // The current point is inside the image and has acceptable FA, so 
    // we will add it to the tract.
    if (forwardDirection) {
      forwardTract->AddVertex(cindex);
    } else {
      backwardTract->AddVertex(cindex);
    }

    // Record our direction.
    // Compute the new direction indicated by Tensor Deflection.
    // TEND is based on the incoming direction, not the current eigenvector.
    if (m_TendFraction > 0.0) {
      // alpha = 1.0
      tendDeflection.Fill(0.0);
      for (unsigned int i = 0; i < 3; i++) {
        tendDeflection[0] += tensor(0,i) * previousEigenvector[i];
        tendDeflection[1] += tensor(1,i) * previousEigenvector[i];
        tendDeflection[2] += tensor(2,i) * previousEigenvector[i];
      }
      tendDeflection /= (tendDeflection.GetNorm());
      // alpha = 2.0
      tendDeflection2.Fill(0.0);
      for (unsigned int i = 0; i < 3; i++) {
        tendDeflection2[0] += tensor(0,i) * tendDeflection[i];
        tendDeflection2[1] += tensor(1,i) * tendDeflection[i];
        tendDeflection2[2] += tensor(2,i) * tendDeflection[i];
      }
      tendDeflection = tendDeflection2 / tendDeflection2.GetNorm();
    }
    previousEigenvector = (previousEigenvector*m_DirectionMomentum +
      (tendDeflection*m_TendFraction + currentEigenvector*(1.0-m_TendFraction))
               *(1.0 - m_DirectionMomentum) );
    // Utilize direction momentum in computing the new step.
    newpoint += ( previousEigenvector * stepSize );
    stepCount++;
  }

  // Now we have the forward tract from the start point to its end,
  // and the backward tract from next to the start point to its end.
  // We need to iterate over the length of the paths to make one tract,
  // and return that.
  //   We go backwards from the end of the backward tract to its first point,
  // and then forwards over the forward tract from its start to its end.
  for (signed int i = backwardTract->GetVertexList()->Size() - 1; i >= 0; i--) {
    cindex = backwardTract->GetVertexList()->GetElement(i);
    tract->AddVertex( cindex );
  }
  for (unsigned int i = 0; i < forwardTract->GetVertexList()->Size(); i++) {
    cindex = forwardTract->GetVertexList()->GetElement(i);
    tract->AddVertex( cindex );
  }
  // We don't delete the forwardTract or reverseTract because the list
  // elements appear to be shared by reference.
}

template <const unsigned int TImageDimension, class TComponentType>
int crl::TractGenerator<TImageDimension, TComponentType>
::Execute()
{
  typedef typename PixelType::ValueType ValueType;
  typedef typename PixelType::EigenVectorsMatrixType MatrixType;
  typedef typename PixelType::EigenValuesArrayType EigenValuesType;

  typedef itk::ImageRegionConstIterator< ImageType> 
                           TensorImageIteratorType;
  
  typename TensorImageReaderType::Pointer reader = TensorImageReaderType::New();
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
  m_TensorImage = reader->GetOutput();

  typename ROIImageReaderType::Pointer ROIReader = ROIImageReaderType::New();
  ROIReader->SetFileName( m_ROIFileName );

  try {
    ROIReader->Update();
  } catch (itk::ExceptionObject &excp)
  {
    std::cerr << "Failed to read from file " << m_ROIFileName << std::endl;    
    std::cerr << "Exception error message is : " << excp << std::endl;
    return EXIT_FAILURE;
  }

  m_ROIImage = ROIReader->GetOutput();

  // Construct a log tensor vector image,
  typedef itk::TensorLogImageFilter<ImageType, ImageType> 
                                                   TensorLogImageFilterType;
  typename TensorLogImageFilterType::Pointer tensorLogFilter =
                                            TensorLogImageFilterType::New();
  tensorLogFilter->SetInput( m_TensorImage );
  tensorLogFilter->Update();
  m_LogTensorImage = tensorLogFilter->GetOutput();

  // connect it to tensor image linear interpolator.
  m_TensorInterpolator->SetInputImage(m_LogTensorImage);

  /*
   *   Use the linear interpolator on the log tensor data to achieve sub-voxel
   * interpolation of the tensor data.
   *   After doing the log-Euclidean interpolation, the largest eigenvector of
   * the log tensor points in the same direction as that of the Euclidean
   * tensor we would get if we exponentiated. 
   *   Therefore, for the purposes of eigenvector tractography, we do not
   * need to compute the exponentiated log-tensor.
   *   However, we do need to know the fractional anisotropy.
   */

  // The generated tracts are stored in the Tract Container
  typename TractContainerType::Pointer tracts = TractContainerType::New();

  /* Iterate over the tensor data, and construct a tract at every location
   * in the ROI.
   */
  TensorImageIteratorType it(m_TensorImage, 
                             m_TensorImage->GetLargestPossibleRegion() );
  PixelType in;

  unsigned long int tractCount = 0;
  unsigned long int tractPointCount = 0;
  unsigned long int tractTopologyCount = 0;
  typename TractType::Pointer tract = 0;
  typename ImageType::PointType point;
  typename ImageType::IndexType index;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    in = it.Get();
    index = it.GetIndex();
    if (m_ROIImage->GetPixel( index ) == m_ROILabel) {
      // Estimate the tract that touches this location.
      for (unsigned int i = 0; i < m_TractsPerVoxel; i++) {
        reader->GetOutput()->TransformIndexToPhysicalPoint(index, point);
        for (unsigned int j = 0; j < ImageType::ImageDimension; j++) {
          point[j] += ( (m_RandGen->GetVariate()-0.5)*
                            m_TensorImage->GetSpacing()[j] );
        }
        EstimateTract(point, tract);
        if (tract->GetVertexList()->Size() > 0) {
          tracts->InsertElement( tractCount, tract );
          ++tractCount;
          tractPointCount += tract->GetVertexList()->Size();
          tractTopologyCount += (tract->GetVertexList()->Size());
        } else {
          ; // std::cout << "Tract has length 0" << std::endl;
        }
      }
    }
  }
  tractTopologyCount += tractCount;
  // std::cout << "Number of tracts identified: " << tractCount << std::endl;

  /* Now we will write out the tracts. We will write an ASCII format VTK file.
   */
  std::ofstream outFile;
  outFile.open( m_OutputFileName.c_str(), std::ios::out);
  if (!outFile.is_open()) {
    std::cerr << "Failed to open file " << m_OutputFileName << std::endl;
    return 1; // Failure
  }
  typename TractType::ContinuousIndexType cindex;
  outFile << "# vtk DataFile Version 3.0" << std::endl;
  outFile << "vtk output" << std::endl;
  outFile << "ASCII" << std::endl << "DATASET POLYDATA" << std::endl;
  outFile << "POINTS " << tractPointCount << " float " << std::endl;
  for (unsigned int i = 0; i < tractCount; i++) {
    tract = tracts->GetElement(i);
    for (unsigned int j = 0; j < tract->GetVertexList()->Size(); j++) {
      cindex = tract->GetVertexList()->GetElement(j);
      m_TensorImage->TransformContinuousIndexToPhysicalPoint(cindex, point);
      for (unsigned int k = 0; k < ImageType::ImageDimension; k++) {
        outFile << point[k] << " " ;
      }
    }
    outFile << std::endl;
  }
  outFile << std::endl;
  // Now print out the poly line data of the tracts.
  outFile << "LINES " << tractCount << " " << tractTopologyCount << std::endl;
  unsigned long pointCounter = 0;
  for (unsigned int i = 0; i < tractCount; i++) {
    tract = tracts->GetElement(i);
    outFile << tract->GetVertexList()->Size() << " ";
    for (unsigned int j = 0; j < tract->GetVertexList()->Size(); j++) {
      outFile << pointCounter << " " ;
      pointCounter++;
    }
    outFile << std::endl;
  }
  outFile << std::endl;
  outFile.close();

  return 0; // Success
}

#endif

