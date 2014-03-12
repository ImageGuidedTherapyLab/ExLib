
/* In order to compute stream lines across the Laplace solution, we
 * will use a numerical integration scheme.
 *  We will use a Carp-Cash Runge-Kutta numerical integration scheme with
 * adaptive step size control.
 *
 *  Given values of n variables y[0,...,n-1] and their derivatives
 * dydx[0,...,n-1] known at x, use the fifth order Cash-Karp Runge-Kutta
 * method to advance the solution over an interval h and return the
 * incremented variables as yout[0,...,n-1]. Also return an estimate of
 * the local truncation error in yout using the embedded fourth-order method.
 * The routine derivates(x,y,dydx) is a user-supplied routine that
 *   returns derivatives dydx at x.
 */
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::a2 = 0.2;
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::a3 = 0.3;
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::a4 = 0.6;
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::a5 = 1.0;
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::a6 = 0.875;

template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::b21 = 0.2;
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::b31 = 3.0/40.0;
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::b32 = 9.0/40.0;
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::b41 = 0.3; 
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::b42 = -0.9;     
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::b43 = 1.2;
         
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::b51 = -11.0/54.0;
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::b52 = 2.5; 
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::b53 = -70.0/27.0; 
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::b54 = 35.0/27.0;

template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::b61 = 1631.0/55296.0; 
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::b62 = 175.0/512.0; 
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::b63 = 575.0/13824.0;
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::b64 = 44275.0/110592.0; 
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::b65 = 253.0/4096.0;
         
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::c1 = 37.0/378.0; 
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::c3 = 250.0/621.0; 
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::c4 = 125.0/594.0; 
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::c6 = 512.0/1771.0;

template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::dc5 = -277.0/14336.0;
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::dc1 = (37.0/378.0 - 2825.0/27648.0); // c1 -
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::dc3 = (250.0/621.0 - 18575.0/48384.0); // c3 -
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::dc4 = (125.0/594.0 - 13525.0/55296.0); // c4 -
template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::dc6 = (512.0/1771.0 - 0.25); // c6 -

template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::SAFETY  = 1.0;

template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::PGROW   = -0.2;

template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::PSHRINK = -0.25;

template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::ERRCON = 1.89e-04;
 // ERRCON is (5/SAFETY) to power (1/PGROW).


template <unsigned int Dimension, class PixelType>
int
crlImageLaplace<Dimension, PixelType >
::Execute()
{
  m_Reader->SetFileName( m_LabelImageFileName );
  m_Writer->SetFileName( m_OutputFileName);
  m_Writer->UseCompressionOn();

  try {
    m_Reader->Update();
  } catch( itk::ExceptionObject & exp )
  {
    std::cerr << "Failed to read image " << m_LabelImageFileName << std::endl;
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
  }

  // Set the m_LabelImage pointer to hold the label image:
  m_LabelImage = m_Reader->GetOutput();

  // Allocate the output image to the size of the input image.
  if (!m_Image) { return EXIT_FAILURE; }
  m_Image->SetSpacing( m_LabelImage->GetSpacing() );
  m_Image->SetOrigin ( m_LabelImage->GetOrigin() );
  m_Image->SetDirection( m_LabelImage->GetDirection() );

  typename ImageType::RegionType region;
  region.SetSize( m_LabelImage->GetLargestPossibleRegion().GetSize());
  region.SetIndex( m_LabelImage->GetLargestPossibleRegion().GetIndex());
  m_Image->SetRegions( region );
  m_Image->Allocate(); 

  typename ImageType::PixelType  initialValue = 0.0;
  m_Image->FillBuffer( initialValue );

  // Allocate the thickness image to the size of the input image.
  if (!m_ThicknessImage) { return EXIT_FAILURE; }
  m_ThicknessImage->SetSpacing( m_Image->GetSpacing() );
  m_ThicknessImage->SetOrigin ( m_Image->GetOrigin() );
  m_ThicknessImage->SetDirection( m_Image->GetDirection() );

  typename ImageType::RegionType thicknessRegion;
  thicknessRegion.SetSize( m_Image->GetLargestPossibleRegion().GetSize());
  thicknessRegion.SetIndex( m_Image->GetLargestPossibleRegion().GetIndex());
  m_ThicknessImage->SetRegions( thicknessRegion );
  m_ThicknessImage->Allocate(); 
  m_ThicknessImage->FillBuffer( initialValue );

  // Allocate the L0 image to the size of the input image.
  if (!m_L0Image) { return EXIT_FAILURE; }
  m_L0Image->SetSpacing( m_Image->GetSpacing() );
  m_L0Image->SetOrigin ( m_Image->GetOrigin() );
  m_L0Image->SetDirection( m_Image->GetDirection() );

  typename ImageType::RegionType L0Region;
  L0Region.SetSize( m_Image->GetLargestPossibleRegion().GetSize());
  L0Region.SetIndex( m_Image->GetLargestPossibleRegion().GetIndex());
  m_L0Image->SetRegions( L0Region );
  m_L0Image->Allocate(); 
  m_L0Image->FillBuffer( initialValue );

  // Allocate the L1 image to the size of the input image.
  if (!m_L1Image) { return EXIT_FAILURE; }
  m_L1Image->SetSpacing( m_Image->GetSpacing() );
  m_L1Image->SetOrigin ( m_Image->GetOrigin() );
  m_L1Image->SetDirection( m_Image->GetDirection() );

  typename ImageType::RegionType L1Region;
  L1Region.SetSize( m_Image->GetLargestPossibleRegion().GetSize());
  L1Region.SetIndex( m_Image->GetLargestPossibleRegion().GetIndex());
  m_L1Image->SetRegions( L1Region );
  m_L1Image->Allocate(); 
  m_L1Image->FillBuffer( initialValue );

  if (!m_OutputTangentImage) { return EXIT_FAILURE; }
  m_OutputTangentImage->SetSpacing( m_Image->GetSpacing() );
  m_OutputTangentImage->SetOrigin ( m_Image->GetOrigin() );
  m_OutputTangentImage->SetDirection( m_Image->GetDirection() );
  typename VectorImageType::RegionType vectorRegion;
  vectorRegion.SetSize( m_Image->GetLargestPossibleRegion().GetSize());
  vectorRegion.SetIndex( m_Image->GetLargestPossibleRegion().GetIndex());
  m_OutputTangentImage->SetRegions( vectorRegion );
  m_OutputTangentImage->Allocate(); 

  // Initializes to be all zero.
  typename VectorImageType::PixelType  initialVectorValue(0.0);
  m_OutputTangentImage->FillBuffer( initialVectorValue );
  // Ensure the vector interpolator is hooked up to the tangent image.
  m_VecInterpolator->SetInputImage( m_OutputTangentImage );

  if (!m_OutputDisplacementFieldImage) { return EXIT_FAILURE; }
  m_OutputDisplacementFieldImage->SetSpacing( m_Image->GetSpacing() );
  m_OutputDisplacementFieldImage->SetOrigin ( m_Image->GetOrigin() );
  m_OutputDisplacementFieldImage->SetDirection( m_Image->GetDirection() );
  typename VectorImageType::RegionType vectorDisplacementFieldRegion;
  vectorDisplacementFieldRegion.SetSize( 
                          m_Image->GetLargestPossibleRegion().GetSize());
  vectorDisplacementFieldRegion.SetIndex( 
                          m_Image->GetLargestPossibleRegion().GetIndex());
  m_OutputDisplacementFieldImage->SetRegions( vectorDisplacementFieldRegion );
  m_OutputDisplacementFieldImage->Allocate(); 
  m_OutputDisplacementFieldImage->FillBuffer( initialVectorValue );


  // Now we iterate over the input label volume, initializing the Laplace 
  // calculation.
  IteratorType it( m_LabelImage, 
                   m_LabelImage->GetLargestPossibleRegion() );
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    PixelType label = it.Get();
    if (label == m_InnerBorderLabel) {
      m_Image->SetPixel(it.GetIndex(), m_InnerBorderVoltage); 
    } else if (label == m_OuterBorderLabel) {
      m_Image->SetPixel(it.GetIndex(), m_OuterBorderVoltage); 
    } else if (label == m_SolveRegionLabel) {
      m_Image->SetPixel(it.GetIndex(), (m_OuterBorderVoltage+m_InnerBorderVoltage)/2.0);
    } else {
      m_Image->SetPixel(it.GetIndex(), 0.0); // Regions where we don't solve.
    }
  }

  // Here we aim to solve Laplace's Equation \grad \dot \grad \psi = 0
  // We use a finite difference discretization and an in place iterative 
  // solver.
  //
  // \nabla u \approx 
  //      \frac{1}{dx^2}\left[ f(x+dx,y,z) +f(d-dx,y,z) -2f(x,y,z)\right] +
  //      \frac{1}{dy^2}\left[ f(x,y+dy,z) +f(d,y-dy,z) -2f(x,y,z)\right] +
  //      \frac{1}{dz^2}\left[ f(x,y,z+dz) +f(x,y,z-dz) -2f(x,y,z)\right] 
  // 
  // Insert the above discretization into Laplace's Equation and solve for 0
  // This leads to an equation for the update value of the scalar field \psi.
  //
  std::vector<double> scalingForDerivative(Dimension);
  double scalingCoefficient = 0.0;
  for (unsigned int i = 0; i < Dimension; i++) {
    scalingForDerivative[i] = m_Image->GetSpacing()[i] *
                              m_Image->GetSpacing()[i];
    scalingForDerivative[i] = 1.0/scalingForDerivative[i];
    scalingCoefficient += 2.0 * scalingForDerivative[i];
  }
  // This is the multiplier needed for Laplace's Equation.
  scalingCoefficient = 1.0 / scalingCoefficient;

  double maxRelativeError = 1.0;
  double relativeError = 0.0;
  unsigned int iterationCount = 0;

  // We will use a Gauss-Seidel iterative scheme, and converge when we
  // have reached a small relative error or a maximum number of iterations.
  while ( (maxRelativeError > m_RelativeErrorTolerance) && 
          (iterationCount < m_MaxNumberOfIterations) ) {
    unsigned long int countSolveRegion = 0;
    IteratorType lit( m_Image, 
                    m_Image->GetLargestPossibleRegion() );
    maxRelativeError = 0.0;

    for (lit.GoToBegin(), it.GoToBegin(); !lit.IsAtEnd(); ++lit, ++it) {
      // We only calculate in the solve region.
      PixelType label = it.Get();
      if (label != m_SolveRegionLabel) continue;
      double center = lit.Get();
      double newLaplaceValue = 0.0;
      double v1 = 0.0;
      double v2 = 0.0;
      typename ImageType::IndexType ind = lit.GetIndex();
      for (unsigned int i = 0; i < Dimension; i++) {
        ind[i] += 1;
        if (m_Image->GetLargestPossibleRegion().IsInside( ind )) {
          v1 = m_Image->GetPixel( ind );
        } else {
          // Reflect the boundary, pretend outside exactly matches the edge
          v1 = center;
        }
        ind[i] -= 2;
        if (m_Image->GetLargestPossibleRegion().IsInside( ind )) {
          v2 = m_Image->GetPixel( ind );
        } else {
          // Reflect the boundary, pretend outside exactly matches the edge
          v2 = center;
        }
        newLaplaceValue += ( (v1 + v2 )*scalingForDerivative[i]); 
        ind[i] += 1;
      }

      newLaplaceValue = newLaplaceValue * scalingCoefficient;
      // Write the updated values into the original image, so that we use
      // the new values immediately.
      lit.Set(newLaplaceValue);

      if (newLaplaceValue != 0.0) {
        relativeError = fabs( (newLaplaceValue - center)/newLaplaceValue );
        if (relativeError > maxRelativeError) {
          maxRelativeError = relativeError;
        }
      }

      countSolveRegion++;
    }
    //std::cout<< "maximum relative error is " << maxRelativeError << std::endl;
    // std::cout << "number of iterations is " << iterationCount << std::endl;

    ++iterationCount;
  }

  // Write out m_Image since it has the updated values.
  m_Writer->SetInput( m_Image );
  try {
    m_Writer->Update();
  } catch (itk::ExceptionObject& err ) {
    std::cerr << "error writing output file: " << m_OutputFileName 
      << ": " << err << std::endl;
    return 1;
  }

  // Now compute the tangent vector field 

  // We iterate over the input label volume, and the tangent vector image
  // calculating the gradient at locations inside the solution volume.
  VectorIteratorType vectorImageIterator( m_OutputTangentImage, 
                    m_OutputTangentImage->GetLargestPossibleRegion() );
  typename VectorImageType::PixelType vPixel;
  for (it.GoToBegin(), vectorImageIterator.GoToBegin(); !it.IsAtEnd(); 
       ++it, ++vectorImageIterator) {
    PixelType label = it.Get();
    if (label != m_SolveRegionLabel) continue;
    // The location is a place where we need to compute the gradient.
    typename VectorImageType::IndexType vIndex = vectorImageIterator.GetIndex();
    typename VectorImageType::IndexType vIndexPlus = vIndex;
    typename VectorImageType::IndexType vIndexMinus = vIndex;

    // Compute the gradient using finite differences, to get the tangent field.
    double gtmp = 0.0;
    for (unsigned int i = 0; i < Dimension; i++) {
      vIndexPlus = vIndex;
      vIndexPlus[i] += 1;
      vIndexMinus = vIndex;
      vIndexMinus[i] -= 1;
      if (m_Image->GetLargestPossibleRegion().IsInside( vIndexPlus ) &&
          m_Image->GetLargestPossibleRegion().IsInside( vIndexMinus ) ) {
        // Central difference is OK
        gtmp = (m_Image->GetPixel(vIndexPlus) - m_Image->GetPixel(vIndexMinus))/
             (2.0 * m_Image->GetSpacing()[i]);
        // Insert partial derivative into i'th component of the vector image
        vPixel[i] = gtmp;
      } else if (m_Image->GetLargestPossibleRegion().IsInside( vIndexPlus ) ) {
        // Forward difference is OK
        gtmp = (m_Image->GetPixel(vIndexPlus) - m_Image->GetPixel(vIndex))/
               (m_Image->GetSpacing()[i]);
        // Insert partial derivative into i'th component of the vector image
        vPixel[i] = gtmp;
      } else if (m_Image->GetLargestPossibleRegion().IsInside( vIndexMinus ) ) {
        // Backward difference is OK
        gtmp = (m_Image->GetPixel(vIndexMinus) - m_Image->GetPixel(vIndex))/
               (m_Image->GetSpacing()[i]);
        // Insert partial derivative into i'th component of the vector image
        vPixel[i] = gtmp;
      } else {
        // No derivative can be taken - let's assume region is flat.
        vPixel[i] = 0.0;
      }
    }
    // Now ensure it is a unit vector
    gtmp = vPixel.GetNorm();
    if (gtmp != 0.0) {
      vPixel /= gtmp;
    }
    // Now save the tangent vector back to the location under consideration
    vectorImageIterator.Set( vPixel );
  }

  m_TangentFieldWriter->SetFileName( m_OutputTangentFieldFileName);
  m_TangentFieldWriter->UseCompressionOn();
  m_TangentFieldWriter->SetInput( m_OutputTangentImage );
  try {
    m_TangentFieldWriter->Update();
  } catch (itk::ExceptionObject& err ) {
    std::cerr << "error writing output file: " << m_OutputTangentFieldFileName 
      << ": " << err << std::endl;
    return 1;
  }

  // Next step: use the tangent field in the computation of a set of two PDEs
  // that will provide us with the value of the thickness at every location
  // inside the solve region.

  // Here we will calculate L0 and L1 as in Yezzi and Prince 2003 IEEE TMI
  maxRelativeError = 1.0;
  relativeError = 0.0;
  iterationCount = 0;

  typename VectorImageType::PixelType vPixelComp;
  
  // We will use a Gauss-Seidel iterative scheme, and converge when we
  // have reached a small relative error or a maximum number of iterations.
  while ( (maxRelativeError > m_RelativeErrorTolerance) && 
          (iterationCount < m_MaxNumberOfIterations) ) {
    IteratorType lit( m_L0Image, 
                      m_L0Image->GetLargestPossibleRegion() );
    maxRelativeError = 0.0;

    for (lit.GoToBegin(), it.GoToBegin(); !lit.IsAtEnd(); ++lit, ++it) {
      // We only calculate in the solve region.
      PixelType label = it.Get();
      if (label != m_SolveRegionLabel) continue;
      typename ImageType::IndexType ind = lit.GetIndex();
      typename ImageType::IndexType vIndexPlus = ind;
      typename ImageType::IndexType vIndexMinus = ind;

      // Compute the L0 PDE
      double coefficient = 0.0;
      double coefficienttmp = 0.0;
      double L0tmp = 1.0;
      double L1tmp = 1.0;
      vPixelComp = m_OutputTangentImage->GetPixel( ind );

      for (unsigned int i = 0; i < Dimension; i++) {
        coefficienttmp = fabs(vPixelComp[i])/m_L0Image->GetSpacing()[i];
        coefficient += coefficienttmp;
        if (vPixelComp[i] >= 0) {
          // For L0, we now do backward difference.
          vIndexMinus = ind;
          vIndexMinus[i] -= 1;
          if (m_L0Image->GetLargestPossibleRegion().IsInside( vIndexMinus ) ) {
            L0tmp += coefficienttmp * m_L0Image->GetPixel(vIndexMinus); 
          }
          // For L1, we now do forward difference.
          vIndexPlus = ind;
          vIndexPlus[i] += 1;
          if (m_L1Image->GetLargestPossibleRegion().IsInside( vIndexPlus ) ) {
            L1tmp += coefficienttmp * m_L1Image->GetPixel(vIndexPlus); 
          }
        } else {
          // For L0, we now do forward difference.
          vIndexPlus = ind;
          vIndexPlus[i] += 1;
          if (m_L0Image->GetLargestPossibleRegion().IsInside( vIndexPlus) ) {
            L0tmp += coefficienttmp * m_L0Image->GetPixel( vIndexPlus ); 
          }
          // For L1, we now do backward difference.
          vIndexMinus = ind;
          vIndexMinus[i] -= 1;
          if (m_L1Image->GetLargestPossibleRegion().IsInside( vIndexMinus ) ) {
            L1tmp += coefficienttmp * m_L1Image->GetPixel( vIndexMinus ); 
          }
        }
      }
      if (coefficient != 0.0) {
        L0tmp = L0tmp / coefficient;
      } else {
        L0tmp = 0.0;
      }
      if (L0tmp != 0.0) {
        relativeError = fabs((L0tmp - m_L0Image->GetPixel(ind))/L0tmp);
        if (relativeError > maxRelativeError) {
          maxRelativeError = relativeError;
        }
      }
      m_L0Image->SetPixel( ind, L0tmp );

      if (coefficient != 0.0) {
        L1tmp = L1tmp / coefficient;
      } else {
        L1tmp = 0.0;
      }
      if (L1tmp != 0.0) {
        relativeError = fabs((L1tmp - m_L1Image->GetPixel(ind))/L1tmp);
        if (relativeError > maxRelativeError) {
          maxRelativeError = relativeError;
        }
      }
      m_L1Image->SetPixel( ind, L1tmp );
    }

   //std::cout << "maximum relative error is " << maxRelativeError << std::endl;
    //std::cout << "number of iterations is " << iterationCount << std::endl;

    ++iterationCount;
  }
  /*
  std::cout << "Completed calculation of L0 and L1 in " 
            << iterationCount << " iterations."
            << " with maxRelativeError " << maxRelativeError <<  std::endl;
  */

  IteratorType lit( m_ThicknessImage, 
                      m_ThicknessImage->GetLargestPossibleRegion() );
  // I guess this could be done only in the solve region.
  double tmp = 0.0;
  for (lit.GoToBegin(); !lit.IsAtEnd(); ++lit) {
    tmp = m_L0Image->GetPixel( lit.GetIndex() )
        + m_L1Image->GetPixel( lit.GetIndex() );
    m_ThicknessImage->SetPixel( lit.GetIndex(), tmp );
  }

  m_ThicknessWriter->SetFileName( m_OutputThicknessFileName);
  m_ThicknessWriter->UseCompressionOn();
  m_ThicknessWriter->SetInput( m_ThicknessImage );
  try {
    m_ThicknessWriter->Update();
  } catch (itk::ExceptionObject& err ) {
    std::cerr << "error writing output file: " << m_OutputThicknessFileName 
      << ": " << err << std::endl;
    return 1;
  }

  /*
  m_ThicknessWriter->SetFileName( "L0Image.nrrd" );
  m_ThicknessWriter->UseCompressionOn();
  m_ThicknessWriter->SetInput( m_L0Image );
  try {
    m_ThicknessWriter->Update();
  } catch (itk::ExceptionObject& err ) {
    std::cerr << "error writing output file: " << "L0Image "
      << ": " << err << std::endl;
    return 1;
  }

  m_ThicknessWriter->SetFileName( "L1Image.nrrd" );
  m_ThicknessWriter->UseCompressionOn();
  m_ThicknessWriter->SetInput( m_L1Image );
  try {
    m_ThicknessWriter->Update();
  } catch (itk::ExceptionObject& err ) {
    std::cerr << "error writing output file: " << "L1Image "
      << ": " << err << std::endl;
    return 1;
  }
  */

  // Compute and store a displacement field representing the displacement to
  // move a particle from any point in the solution region to the end of the
  // streamline.
  for (lit.GoToBegin(), it.GoToBegin(); !lit.IsAtEnd(); ++lit, ++it) {
      // We only calculate in the solve region.
      PixelType label = it.Get();
      if (label != m_SolveRegionLabel) continue;
      VectorPixelType displacement;

      itk::Point<double, Dimension> start;
      std::vector< itk::Point<double, Dimension> > ysave(0);
      m_Image->TransformIndexToPhysicalPoint(it.GetIndex(), start);
      ComputeStreamLine(start, ysave);
      if (ysave.size() > 0) {
        for (unsigned int i = 0; i < Dimension; i++) {
          displacement[i] = ysave[ysave.size() - 1][i] - start[i];
        }
      }
      m_OutputDisplacementFieldImage->SetPixel(it.GetIndex(), displacement);
  }

  m_DisplacementFieldWriter->SetFileName( m_OutputDisplacementFieldFileName );
  m_DisplacementFieldWriter->UseCompressionOn();
  m_DisplacementFieldWriter->SetInput( m_OutputDisplacementFieldImage );
  try {
    m_DisplacementFieldWriter->Update();
  } catch (itk::ExceptionObject& err ) {
    std::cerr << "error writing output file: " << 
        m_OutputDisplacementFieldFileName << err << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

  // Initiate a stream line calculation at the specified start point,
  // follow it by numerical integration along the tangent field until
  // the path reaches out of the solution region, project the end point
  // back to the boundary of the solution region, and return the set of
  // steps from the start to the end representing the stream line.

template <unsigned int Dimension, class PixelType>
int
crlImageLaplace<Dimension, PixelType >
::ComputeStreamLine(
    itk::Point<double, Dimension> start,
    std::vector< itk::Point<double, Dimension> > &ysave )
{
    // We require that: 
    // 1. the tangent field has been computed.
    // 2. the label volume provides the inner, outer and solve region.

    // 3. start - will represent a physical start location in space.
    //        The stream line path will be computed and 
    //    stored in physical coordinates.

    const static unsigned int MAXSTEPS = 10000;
    const static double TINY = 1.0e-30;
    const static double eps = 1.0e-2;

    crlRungeKutta<Dimension, PixelType> rk(this);

    double s = 0.0; // Start integrating at s=0.0 on the streamline.

    itk::Point<double, Dimension> yscal;
    itk::Point<double, Dimension> y;
    itk::Point<double, Dimension> dydx;

    itk::Point<double, Dimension> vin;
    itk::Point<double, Dimension> vout;
    itk::Point<double, Dimension> vmid;

    // vectors along the cardinal axes for computing a last step in some cases.
    itk::Point<double, Dimension> delta[Dimension];
    for (unsigned int i = 0; i < Dimension; i++) {
      for (unsigned int j = 0; j < Dimension; j++) {
        if (i == j) {
          delta[i][j] = m_LabelImage->GetSpacing()[i];
        } else {
          delta[i][j] = 0.0;
        }
      }
    }

    ysave.resize(0);

    double hmin = 0.0;
    double hinitial = 0.0;

    // Compute the smallest step based on a quarter of the root mean square 
    // size of each voxel.
    for (unsigned int i = 0; i < Dimension; i++) {
      hmin += (m_LabelImage->GetSpacing()[i]*m_LabelImage->GetSpacing()[i]);
    }
    hmin = 0.25*sqrt(hmin);
    hinitial = 2.0 * hmin;

    double h = hinitial;
    double hdid = 0.0;
    double hnext = 0.0;

    int nok = 0;
    int nbad = 0;

    for (unsigned int i = 0; i < Dimension; i++) {
      y[i] = start[i];
    }
    unsigned int nsteps = 0;

    typename ImageType::IndexType index;
    double norm = 0.0;
    for (nsteps = 0; nsteps < MAXSTEPS; nsteps++) {
      if (!  this->Derivative(s,y,dydx) ) {
        // The derivative can fail to be calculatable when we step outside
        // the solution zone.
        // The norm of the tangent field gets smaller than one due to 
        // the linear interpolation when we approach close to the boundary.
        ; // Do nothing for now.
      }
      // Check on norm of dydx to see if we are close to solution boundary.
      norm = 0.0;
      for (unsigned int i = 0; i < Dimension; i++) {
        norm += (dydx[i]*dydx[i]);
      }
      if (norm < 0.64) { 
        // This is too small - terminate the search.
        // It indicates linear interpolation between a unit magnitude tangent
        // vector, and zero vectors outside of the solution region. Therefore,
        // it indicates we are at the edge of the region.
        ysave.push_back(y);  // Save the current step
        if (!OneVoxelStep(y, vout)) {
          // We can't even step one more voxel - we must be at the edge of
          // the image.
          ; // Nothing more to do.
        } else {
          // We will save the last step position as the final position.
          ysave.push_back(vout);
        }
        return EXIT_SUCCESS; // This is the only successful exit point
      }

      for (unsigned int i = 0; i < Dimension; i++) {
        // Initialize the scaling factors used to monitor accuracy.
        yscal[i] = fabs(y[i])+fabs(dydx[i]*h) + TINY;
      }
      ysave.push_back(y); // Save the current step
      for (unsigned int i = 0; i < Dimension; i++) {
        vin[i] = y[i];
      }
      // This will update y and s.
      rk.rkqs(y, dydx, s, h, eps, yscal, hdid, hnext);
      if (hdid == h) {
        ++nok;
      } else {
        ++nbad;
      }
      // Check if y is now outside the solve region, indicating we are done.
      for (unsigned int i = 0; i < Dimension; i++) {
        vout[i] = y[i];
      }
      bool isinside = m_LabelImage->TransformPhysicalPointToIndex(vout, index);
      if ((isinside) && (m_LabelImage->GetPixel(index) != m_SolveRegionLabel)) {
        // stop integrating now that we are out of the solve region.
        // binary search for the exact location of the boundary we match to.
        // search between y and yout for the location of the boundary.
        // push the boundary location as the last step.
        double diffnorm = 0.0;
        for (unsigned int i = 0; i < Dimension; i++) {
          diffnorm += ( (vin[i] - vout[i])*(vin[i] - vout[i]) );
        }
        diffnorm = sqrt(diffnorm);
        while (diffnorm > hmin) {
          for (unsigned int i = 0; i < Dimension; i++) {
            vmid[i] = 0.5*(vin[i] + vout[i]);
          }
          isinside = m_LabelImage->TransformPhysicalPointToIndex( vmid, index );
          if ( !isinside || 
               (m_LabelImage->GetPixel( index ) != m_SolveRegionLabel) ) {
            vout = vmid;
          } else {
            vin = vmid;
          }
          diffnorm = 0.0;
          for (unsigned int i = 0; i < Dimension; i++) {
            diffnorm += ( (vin[i] - vout[i])*(vin[i] - vout[i]) );
          }
          diffnorm = sqrt(diffnorm);
        }
        // We have now determined that vout is very close (within hmin) of the
        // border of the solve region, so we can declare it the end point.
        ysave.push_back(vout); // Save the current step
        return EXIT_SUCCESS; // This is the only successful exit point
      } else if (!isinside) {
        std::cerr << "Step is now outside the image. We need to handle that."
            << std::endl;
      }
      if (fabs(hnext) <= hmin) {
        //std::cerr << "Step size too small in ComputeStreamLine." << std::endl;
        if (!OneVoxelStep(y, vout)) {
          // We can't even step one more voxel - we must be at the edge of
          // the image.
          ; // Nothing more to do.
        } else {
          // We will save the last step position as the final position.
          ysave.push_back(vout);
        }
        return EXIT_SUCCESS;
      }
      h = hnext;
    }
    std::cerr << "Too many steps to compute the streamline." << std::endl;
    return EXIT_FAILURE;
};

template <unsigned int Dimension, class PixelType>
bool
crlImageLaplace<Dimension, PixelType>
::Derivative(double s, 
           itk::FixedArray<double, Dimension> &y,
           itk::FixedArray<double, Dimension> &dydx)
{
  // Compute the value of the derivative of the stream line function:
  // dy/ds = N(C(s))
  //  Here, s is the arc length parameterization of the streamline curve,
  // y is the position in space that we are integrating at, and
  // dy/dx is the derivative value to return.

  // The algorithm here is to:
  //   look up the location of y in the tangent field image, and return
  // the value of the interpolated tangent field.
  typename VectorInterpolatorType::PixelType value;
  itk::Point<double, Dimension> p;
  for (unsigned int i = 0; i < Dimension; i++) {
    p[i] = y[i];
    dydx[i] = 0.0; // Initialize the derivative to be zero.
  }
  typename ImageType::IndexType lindex;
  if (!m_LabelImage->TransformPhysicalPointToIndex( p, lindex )) {
    return false; // Cannot compute derivative.
  }
  if ( m_LabelImage->GetPixel( lindex ) != m_SolveRegionLabel ) {
    // std::cout << "Asking for derivative outside of the solve region." << std::endl;
    // std::cout << "Reduce the step size and try again..." << std::endl;
    return false; // Cannot compute derivative
  }
  if ( m_VecInterpolator->IsInsideBuffer( p ) ) {
     value = m_VecInterpolator->Evaluate( p );
  } else {
    // Can't interpolate a value that is outside the image.
    return false; // Cannot compute derivative
  }
  for (unsigned int i = 0; i < Dimension; i++) {
    dydx[i] = value[i];
  }
  return true; // Successful derivative calculation.

};

