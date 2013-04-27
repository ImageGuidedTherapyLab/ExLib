
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
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::SAFETY  = 0.999; //1.0;

template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::PGROW   = -0.2;

template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::PSHRINK = -0.25;

template <unsigned int Dimension, class PixelType>
const double ITK_EXPORT crlRungeKutta<Dimension, PixelType>::ERRCON = 1.89e-04;
 // ERRCON is (5/SAFETY) to power (1/PGROW).

template <unsigned int Dimension, class PixelType>
int crlImageLaplace<Dimension, PixelType>::SolveLaplace()
{
  // Now we iterate over the input label volume, initializing the Laplace 
  // calculation.

  m_LapSolImage->FillBuffer(0);

  std::cout << "    Initializing Solution to Laplace's Equation" << std::endl;
  IteratorType it( m_LabelImage, m_LabelImage->GetLargestPossibleRegion() );

  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    PixelType label = it.Get();
    if (label == m_InnerBorderLabel) {
      m_LapSolImage->SetPixel(it.GetIndex(), m_InnerBorderVoltage); 
    } else if (label == m_OuterBorderLabel) {
      m_LapSolImage->SetPixel(it.GetIndex(), m_OuterBorderVoltage); 
    } else if ( (label == m_CorticalLabelNew) || 
		( label == m_CorticalLabelChecked) || 
		(!m_ExcludeSulci && ((label==m_TempLabel) || ( label == m_SulcalLabel))) ) {
      m_LapSolImage->SetPixel(it.GetIndex(), (m_OuterBorderVoltage+m_InnerBorderVoltage)/2.0);
    } else {
      m_LapSolImage->SetPixel(it.GetIndex(), m_OuterBorderVoltage); // Regions where we don't solve.
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
	scalingForDerivative[i] = m_LapSolImage->GetSpacing()[i] *
				  m_LapSolImage->GetSpacing()[i];
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
      std::cout << "    Solving Using Gauss-Seidel" << std::endl;
      while ( (maxRelativeError > m_RelativeErrorTolerance) && 
	      (iterationCount < m_MaxNumberOfIterations) ) {
	unsigned long int countSolveRegion = 0;
	IteratorType lit( m_LapSolImage, m_LapSolImage->GetLargestPossibleRegion() );
	maxRelativeError = 0.0;

	for (lit.GoToBegin(), it.GoToBegin(); !lit.IsAtEnd(); ++lit, ++it) {
	  // We only calculate in the solve region.
	  PixelType label = it.Get();
	  if ((label != m_CorticalLabelNew)     && 
              (label != m_CorticalLabelChecked) && 
              (m_ExcludeSulci || ((label !=m_TempLabel) && (label != m_SulcalLabel)))) continue;
	  double center = lit.Get();
	  double newLaplaceValue = 0.0;
	  double v1 = 0.0;
	  double v2 = 0.0;
	  

	  typename ImageType::IndexType ind = lit.GetIndex();
	  for (unsigned int i = 0; i < Dimension; i++) {
	    ind[i] += 1;
	    if (m_LapSolImage->GetLargestPossibleRegion().IsInside( ind )) {
	      v1 = m_LapSolImage->GetPixel( ind );
	    } else {
	      // Reflect the boundary, pretend outside exactly matches the edge
	      v1 = center;
	    }
	    ind[i] -= 2;
	    if (m_LapSolImage->GetLargestPossibleRegion().IsInside( ind )) {
	      v2 = m_LapSolImage->GetPixel( ind );
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

	/*(
	std::cout << "Maximum relative error is "  << maxRelativeError ;
	std::cout << " after " << iterationCount << " iterations." << std::endl;
	*/

	iterationCount++;    
      }

      std::cout << "    Completed Gauss-Seidel Solution in " << iterationCount << " iterations." <<  std::endl;
  return 1;


}

template <unsigned int Dimension, class PixelType> 
int crlImageLaplace<Dimension, PixelType>::ComputeGradient()
{

  // Now compute the tangent vector field 
  std::cout << "    Reinitializing tangent image" << std::endl;
  typename VectorImageType::PixelType initialVectorValue(0.0);
  m_TangentImage->FillBuffer(initialVectorValue);

  // We iterate over the input label volume, and the tangent vector image
  // calculating the gradient at locations inside the solution volume.
  IteratorType       it(     m_LapSolImage,     m_LapSolImage->GetLargestPossibleRegion() );
  IteratorType label_it(m_LabelImage,m_LabelImage->GetLargestPossibleRegion() );
  VectorIteratorType vectorImageIterator( m_TangentImage, 
                    m_TangentImage->GetLargestPossibleRegion() );

  typename VectorImageType::PixelType vPixel;
  std::cout << "    Computing normalized tangent vectors" << std::endl;
  for (label_it.GoToBegin(), it.GoToBegin(), vectorImageIterator.GoToBegin(); !it.IsAtEnd(); 
       ++it, ++vectorImageIterator, ++label_it) {

    // Check if we're withing the relabeled solution area (Smaller than or equal to the total cortical area)
    PixelType label = label_it.Get();	 
    /*    
	  if (label != 0)
	{
	  std::cout << "Pixel is labeled: " << label << " compared to a New Cortical label of " << m_CorticalLabelNew 
                    << " and a Checked Cortical label " << m_CorticalLabelChecked << std::endl;
	}
    */
    if ((label != m_CorticalLabelNew) && 
        (label != m_CorticalLabelChecked) && 
        (m_ExcludeSulci || ((label !=m_TempLabel) && (label != m_SulcalLabel)))) continue;

    // The location is a place where we need to compute the gradient.
    typename VectorImageType::IndexType vIndex      = vectorImageIterator.GetIndex();
    typename VectorImageType::IndexType vIndexPlus  = vIndex;
    typename VectorImageType::IndexType vIndexMinus = vIndex;

    // Compute the gradient using finite differences, to get the tangent field.
    double gtmp = 0.0;
    for (unsigned int i = 0; i < Dimension; i++) {
      // Get Indices for Forward and Reverse Directions
      vIndexPlus = vIndex;
      vIndexPlus[i]  += 1;
      vIndexMinus = vIndex;
      vIndexMinus[i] -= 1;

      if (m_LapSolImage->GetLargestPossibleRegion().IsInside( vIndexPlus ) &&
          m_LapSolImage->GetLargestPossibleRegion().IsInside( vIndexMinus ) ) {
	//	std::cout << "Central Difference" << std::endl;
        // Central difference is OK
        gtmp = (m_LapSolImage->GetPixel(vIndexPlus) - m_LapSolImage->GetPixel(vIndexMinus))/
             (2.0 * m_LapSolImage->GetSpacing()[i]);
        // Insert partial derivative into i'th component of the vector image
        vPixel[i] = gtmp;
      } 
      else if (m_LapSolImage->GetLargestPossibleRegion().IsInside( vIndexPlus ) ) {
        // Forward difference is OK
	//std::cout << "Forward Difference" << std::endl;
        gtmp = (m_LapSolImage->GetPixel(vIndexPlus) - m_LapSolImage->GetPixel(vIndex))/
               (m_LapSolImage->GetSpacing()[i]);
        // Insert partial derivative into i'th component of the vector image
        vPixel[i] = gtmp;
      } 
      else if (m_LapSolImage->GetLargestPossibleRegion().IsInside( vIndexMinus ) ) {
        // Backward difference is OK
	//	std::cout << "Backward Difference" << std::endl;
	//        gtmp = (m_LapSolImage->GetPixel(vIndexMinus) - m_LapSolImage->GetPixel(vIndex))/
        //       (m_LapSolImage->GetSpacing()[i]);
	gtmp = (m_LapSolImage->GetPixel(vIndex) - m_LapSolImage->GetPixel(vIndex))/ (m_LapSolImage->GetSpacing()[i]);
        // Insert partial derivative into i'th component of the vector image
        vPixel[i] = gtmp;
      } 
      else {
	//	std::cout << " No Derivative! " << std::endl;
        // No derivative can be taken - let's assume region is flat.
        vPixel[i] = 0.0;
      }
    } // END Dimension Loop

    // Now ensure it is a unit vector
    gtmp = vPixel.GetNorm();
    if (gtmp != 0.0) {
      vPixel /= gtmp;
    }

    // Now save the tangent vector back to the location under consideration
    vectorImageIterator.Set( vPixel );
  } // END Voxel Loop

  bool NonZeroVal = false;
  double val = 0;
  unsigned int i;

  for (vectorImageIterator.GoToBegin(); !vectorImageIterator.IsAtEnd() ; ++vectorImageIterator)
    {
      vPixel = vectorImageIterator.Get();
      val = 0;
      for (i=0 ; i<Dimension ; i++)
	{
	  val += vPixel[i]*vPixel[i];
	}
      if (val>0) NonZeroVal = true;
    }

  if (!NonZeroVal)
    {
      std::cout << "All Tangent Vectors Have Magnitude Zero!!!!" << std::endl;
      return 0;
    }

  std::cout << "    Successfully completed computation of tangent vectors" << std::endl;

 return 1;

}

template <unsigned int Dimension, class PixelType>
int crlImageLaplace<Dimension, PixelType>::SolveYezziPrince()
{
 
  // Next step: use the tangent field in the computation of a set of two PDEs
  // that will provide us with the value of the thickness at every location
  // inside the solve region.

  double maxRelativeError;
  double relativeError;
  unsigned int iterationCount;

  // Here we will calculate L0 and L1 as in Yezzi and Prince 2003 IEEE TMI
  maxRelativeError = 1.0;
  relativeError = 0.0;
  iterationCount = 0;

 

  // Define Iterators
  IteratorType Label_it( m_LabelImage, m_LabelImage->GetLargestPossibleRegion() );
  IteratorType L0_it(    m_L0Image,       m_L0Image->GetLargestPossibleRegion() );
  
  // Reinitialize L0 and L1
  typename VectorImageType::PixelType vPixelComp;
  typename       ImageType::PixelType initialValue = 0.0;
  m_L0Image->FillBuffer(initialValue);
  m_L1Image->FillBuffer(initialValue);

  bool allZeros = true;
  for (L0_it.GoToBegin() ; !L0_it.IsAtEnd() ; ++L0_it)
    {
      if ( L0_it.Get() != 0 )
	{
	  allZeros = false;
	}
    }
  if (!allZeros) { std::cout << " L0 Image is Not All Zeros!!! " << std:: endl; }


  // We will use a Gauss-Seidel iterative scheme, and converge when we
  // have reached a small relative error or a maximum number of iterations.
  while ( ( maxRelativeError > m_RelativeErrorTolerance ) && 
          (   iterationCount < m_MaxNumberOfIterations  ) ) {
    unsigned long int countSolveRegion = 0;
   
    maxRelativeError = 0.0;

    for (Label_it.GoToBegin(), L0_it.GoToBegin(); !L0_it.IsAtEnd(); ++L0_it, ++Label_it) {
      // We only calculate in the solve region.
      PixelType label = Label_it.Get();
      if ((label != m_CorticalLabelNew) && 
          (label != m_CorticalLabelChecked) && 
          (m_ExcludeSulci || ((label != m_TempLabel) && (label != m_SulcalLabel))))  continue;
      typename ImageType::IndexType ind = L0_it.GetIndex();
      typename ImageType::IndexType vIndexPlus = ind;
      typename ImageType::IndexType vIndexMinus = ind;

      // Compute the L0 PDE
      double coefficient = 0.0;
      double coefficienttmp = 0.0;
      double L0tmp = 1.0;
      double L1tmp = 1.0;
      vPixelComp = m_TangentImage->GetPixel( ind );

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
      if (L0tmp != 0.0) {      countSolveRegion++;
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

    //std::cout << "At iteration " << iterationCount << " Maximum relative error is " << maxRelativeError << 
    // " Compared to threshold of " << m_RelativeErrorTolerance << std::endl;
    //std::cout << "number of iterations is " << iterationCount << std::endl;

    ++iterationCount;
  }

  std::cout << "    Completed Yezzi-Prince Computations in " << iterationCount << " iterations" << std::endl;

  std::cout << "    Computing tissue thicknesses" << std::endl;

  IteratorType Thickness_it( m_ThicknessImage, m_ThicknessImage->GetLargestPossibleRegion() );
  // I guess this could be done only in the solve region.
  double tmp = 0.0;
  for (Thickness_it.GoToBegin(); !Thickness_it.IsAtEnd(); ++Thickness_it) 
  {
    tmp = m_L0Image->GetPixel( Thickness_it.GetIndex() )
        + m_L1Image->GetPixel( Thickness_it.GetIndex() );
    m_ThicknessImage->SetPixel( Thickness_it.GetIndex(), tmp );
  }

  std::cout << "    Completed thickness computations" << std::endl; 

 return 1;

}

//
// BEGIN crlImageLaplace::CheckTopology6
//


template <unsigned int Dimension, class PixelType>
bool crlImageLaplace<Dimension, PixelType>::CheckLocalTopology6(typename ImageType::IndexType idxVox, double baseLabel,  double newLayerLabel)
{
  int i,j;
  bool valid_add = true;
  
  // Initialize Component Labels to Zero
  int CompLabels[125],CompCount=0;
  for (i=0;i<125;i++)
	CompLabels[i] = 0;
  int Que[125], Que_pos, Que_len;

  int nbrs6[6] = {37, 57, 61, 63, 67, 87};
  int nbrs18[18] = {32, 36, 37, 38, 42, 56, 57, 58, 61, 63, 66, 67, 68, 82, 86, 87, 88, 92};
  int localIdx,queIdx,nbrIdx;

  typename ImageType::PixelType nbrPixel;
  typename ImageType::PixelType nbrsNbrPixel;


  // Set up neighborhood iterators for the current voxel and its neighbors
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(2);
  NeighborhoodIteratorType itLocal(radius, m_LabelImageTmp,m_LabelImageTmp->GetRequestedRegion());

  itLocal.SetLocation(idxVox);
  // Construct the connected component(s) associated with the neighboring voxels
  // Loop across 18-Neighborhood
  for (i = 0; i<6 ; i++) {    
    localIdx = nbrs6[i];
    nbrPixel  = itLocal.GetPixel(localIdx);

    // If the voxel is the underlying type and we haven't assigned
    // a component to it yet.
    if ( (CompLabels[localIdx]==0) && (nbrPixel==baseLabel) ){
      // Assign a New Component Label
      CompLabels[localIdx] = ++CompCount;
      
      // Find the full connected component
      Que_pos=0; Que_len=1; Que[0] = localIdx;
      while (Que_pos<Que_len) {
	queIdx = Que[Que_pos]; 		
	for (j=1; j<=m_eListExpanded6[queIdx][0];j++){
	  nbrIdx = m_eListExpanded6[queIdx][j];      
	  nbrsNbrPixel = itLocal.GetPixel(nbrIdx);

	  if ((CompLabels[nbrIdx]==0) && (nbrsNbrPixel==baseLabel)) {
	    CompLabels[nbrIdx] = CompCount;
	    Que[Que_len++]=nbrIdx;
	  }

          if ((CompLabels[nbrIdx]==0) && (nbrsNbrPixel==newLayerLabel)) {
	    CompLabels[nbrIdx] = CompCount;
	  }
	} // END for
	Que_pos++;
      } // END while        
    } // END if
  } // END for

  if (CompCount>1)  // If more than one component was identified, there is a topology violation
    {
       valid_add = false;
    }
  else {
    // Check newly added voxels for adjacent components
    for (i = 0; i<18 ; i++) {
      localIdx = nbrs18[i];
      nbrPixel  = itLocal.GetPixel(localIdx);
      
      // If the voxel is the underlying type and we haven't assigned
      // a component to it yet.
      if ( (CompLabels[localIdx]==0) && (nbrPixel==newLayerLabel) ) {
	valid_add = false;
      } // END if        
    } // END for
  } // END if 

  //  std::cout << "6:"<<  idxVox[0] << " " << idxVox[1] << " " << idxVox[2] << std::endl;
  int k;
  if ((idxVox[0]==127)&(idxVox[1]==119)&(idxVox[2]==129)){
    std::cout << "Found Voxel " << valid_add << std::endl;
 
    //    if (!valid_add){
    std::cout<< "CompLabels" << std::endl;
    for (i=0;i<125;i=i+25){
      for (k=i;k<i+25;k=k+5){
	for (j=k;j<k+5;j++){
	  std::cout << std::setw(2) << CompLabels[j] << " " ;
	}
	std::cout << "    ";
	for (j=k;j<k+5;j++){
	  std::cout << std::setw(2) << itLocal.GetPixel(j) << " " ;
	}
	std::cout << std::endl;
      }
      std::cout<< std::endl;
    }
  }
  
  return(valid_add);
}


//
// BEGIN crlImageLaplace::CheckTopology18
//

template <unsigned int Dimension, class PixelType>
bool crlImageLaplace<Dimension, PixelType>::CheckLocalTopology18(typename ImageType::IndexType idxVox, double baseLabel,  double newLayerLabel)
{
  int i,j;
  bool valid_add = true;
  
  // Initialize Component Labels to Zero
  int CompLabels[125],CompCount=0;
  for (i=0;i<125;i++)
	CompLabels[i] = 0;
  int Que[125], Que_pos, Que_len;

  int nbrs18[18] = {32, 36, 37, 38, 42, 56, 57, 58, 61, 63, 66, 67, 68, 82, 86, 87, 88, 92};
  int localIdx,queIdx,nbrIdx;

  typename ImageType::PixelType nbrPixel;
  typename ImageType::PixelType nbrsNbrPixel;


  // Set up neighborhood iterators for the current voxel and its neighbors
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(2);
  NeighborhoodIteratorType itLocal(radius, m_LabelImageTmp,m_LabelImageTmp->GetRequestedRegion());

  itLocal.SetLocation(idxVox);
  // Construct the connected component(s) associated with the neighboring voxels
  // Loop across 18-Neighborhood
  for (i = 0; i<18 ; i++) {    
    localIdx = nbrs18[i];
    nbrPixel  = itLocal.GetPixel(localIdx);

    // If the voxel is the underlying type and we haven't assigned
    // a component to it yet.
    if ( (CompLabels[localIdx]==0) && (nbrPixel==baseLabel) ){
      // Assign a New Component Label
      CompLabels[localIdx] = ++CompCount;
      
      // Find the full connected component
      Que_pos=0; Que_len=1; Que[0] = localIdx;
      while (Que_pos<Que_len) {
	queIdx = Que[Que_pos]; 		
	for (j=1; j<=m_eListExpanded18[queIdx][0];j++){
	  nbrIdx = m_eListExpanded18[queIdx][j];      
	  nbrsNbrPixel = itLocal.GetPixel(nbrIdx);

	  if ((CompLabels[nbrIdx]==0) && (nbrsNbrPixel==baseLabel)) {
	    CompLabels[nbrIdx] = CompCount;
	    Que[Que_len++]=nbrIdx;
	  }

          if ((CompLabels[nbrIdx]==0) && (nbrsNbrPixel==newLayerLabel)) {
	    CompLabels[nbrIdx] = CompCount;
	  }
	} // END for
	Que_pos++;
      } // END while        
    } // END if
  } // END for

  if (CompCount>1)  // If more than one component was identified, there is a topology violation
    {
       valid_add = false;
    }
  else {
    // Check newly added voxels for adjacent components
    for (i = 0; i<18 ; i++) {
      localIdx = nbrs18[i];
      nbrPixel  = itLocal.GetPixel(localIdx);
      
      // If the voxel is the underlying type and we haven't assigned
      // a component to it yet.
      if ( (CompLabels[localIdx]==0) && (nbrPixel==newLayerLabel) ) {
	valid_add = false;
      } // END if        
    } // END for
  } // END if 
  
  // std::cout << idxVox[0] << " " << idxVox[1] << " " << idxVox[2] << std::endl;

  int k;
  if ((idxVox[0]==127)&(idxVox[1]==119)&(idxVox[2]==129)){
    std::cout << "Found Voxel " << valid_add << std::endl;
    //    if (!valid_add){

    std::cout<< "CompLabels" << std::endl;
    for (i=0;i<125;i=i+25){
      for (k=i;k<i+25;k=k+5){
	for (j=k;j<k+5;j++){
	  std::cout << std::setw(2) << CompLabels[j] << " " ;
	}
	std::cout << "    ";
	for (j=k;j<k+5;j++){
	  std::cout << std::setw(2) << itLocal.GetPixel(j) << " " ;
	}
	std::cout << std::endl;
      }
      std::cout<< std::endl;
    }
  }


  return(valid_add);
}



//
//  BEGIN crlImageLaplace::AddLayer
//


template <unsigned int Dimension, class PixelType>
bool crlImageLaplace<Dimension, PixelType>::AddLayer(double baseLabel, double addToLabel , double newLayerLabel, bool DiagExpansion, bool CheckTopology )
{
 //   std::cout << "Converting all voxels labeled " << addToLabel << " adjacent to " << baseLabel << " labeled voxels to new label " << newLayerLabel << std::endl;

    bool MadeChanges = false;

    // Construct Iterators for LabelImage and LabelImgTmp
    IteratorType label_it( m_LabelImage    ,    m_LabelImage->GetLargestPossibleRegion() );
    IteratorType   tmp_it( m_LabelImageTmp , m_LabelImageTmp->GetLargestPossibleRegion() );

    // Copy data from m_LabelImage to m_LabelImageTmp
    for (label_it.GoToBegin(), tmp_it.GoToBegin(); !label_it.IsAtEnd(); ++label_it, ++tmp_it)
    {
     tmp_it.Set(label_it.Get());
    }

    typename ImageType::PixelType  NewPixel = newLayerLabel;
    typename ImageType::PixelType  OldPixel = addToLabel;

    // Loop across voxels
    //for (label_it.GoToBegin(), tmp_it.GoToBegin(); !label_it.IsAtEnd(); ++label_it, ++tmp_it) {
    for (label_it.GoToBegin() ; !label_it.IsAtEnd(); ++label_it)
      {

	// Check that we're adding a layer onto the appropriate underlying type
	PixelType label = label_it.Get();
	if (label != baseLabel) continue;
	
	typename ImageType::IndexType ind = label_it.GetIndex();
	
	
	// Add immediately adjacent voxels
	for (unsigned int i = 0; i < Dimension; i++) {
	  ind[i] += 1;
	  if (m_LabelImage->GetLargestPossibleRegion().IsInside( ind )) {         
	    PixelType TestPix = m_LabelImage->GetPixel(ind);
	    
	    if ((TestPix ==  addToLabel)&(m_LabelImageTmp->GetPixel(ind) != newLayerLabel)) {
	      m_LabelImageTmp->SetPixel(ind,NewPixel);
	      // std::cout << ind << std::endl;
	      MadeChanges = true;
	    }
	  }
	  ind[i] -= 2;
	  if (m_LabelImage->GetLargestPossibleRegion().IsInside( ind )) {
	    PixelType TestPix = m_LabelImage->GetPixel(ind);
	    if ((TestPix ==  addToLabel)&(m_LabelImageTmp->GetPixel(ind) != newLayerLabel)) {
	      m_LabelImageTmp->SetPixel(ind,NewPixel);
	      // std::cout << ind << std::endl;
	      MadeChanges = true;
	    }         
	  }
	  ind[i] += 1;
	}
	
	if (DiagExpansion) { // We're also adding voxels with a taxicab distance of 2 or less.
	  for (unsigned int j = 0 ; j < Dimension ; j++) {

	    ind[j] += 1;
	    for (unsigned int i = 0; i < Dimension; i++) {
	      if ( i!=j){
		ind[i] += 1;
		if (m_LabelImage->GetLargestPossibleRegion().IsInside( ind )) {         
		  PixelType TestPix = m_LabelImage->GetPixel(ind);
		  
		  if ((TestPix ==  addToLabel)&(m_LabelImageTmp->GetPixel(ind) != newLayerLabel)) {
		    m_LabelImageTmp->SetPixel(ind,NewPixel);
		    //  std::cout << ind << std::endl;
		    MadeChanges = true;
		  }
		}
		ind[i] -= 2;
		if (m_LabelImage->GetLargestPossibleRegion().IsInside( ind )) {
		  PixelType TestPix = m_LabelImage->GetPixel(ind);
		  if ((TestPix ==  addToLabel)&(m_LabelImageTmp->GetPixel(ind) != newLayerLabel))  {
		    m_LabelImageTmp->SetPixel(ind,NewPixel);
		    // std::cout << ind << std::endl;
		    MadeChanges = true;
		  }         
		}
		ind[i] += 1;
	      }
	    } // END 1st FOR loop
		
	    ind[j] -=2;
	    for (unsigned int i = 0; i < Dimension; i++) {
	      if (i != j)   {
		ind[i] += 1;
		if (m_LabelImage->GetLargestPossibleRegion().IsInside( ind )) {         
		  PixelType TestPix = m_LabelImage->GetPixel(ind);
		  
		  if ((TestPix ==  addToLabel)&(m_LabelImageTmp->GetPixel(ind) != newLayerLabel))  {
		    m_LabelImageTmp->SetPixel(ind,NewPixel);
		    // std::cout << ind << std::endl;
			    MadeChanges = true;
		  }
		}
		ind[i] -= 2;
		if (m_LabelImage->GetLargestPossibleRegion().IsInside( ind )) {
		  PixelType TestPix = m_LabelImage->GetPixel(ind);
		  if ((TestPix ==  addToLabel)&(m_LabelImageTmp->GetPixel(ind) != newLayerLabel)) {
		    m_LabelImageTmp->SetPixel(ind,NewPixel);
		    //  std::cout << ind << std::endl;
		    MadeChanges = true;
		  }         
		}
		ind[i] += 1;
	      }
	    } // END 2nd FOR loop
		
	    ind[j] += 1;
	    
	  } // END outer FOR loop
	} // END if (diagExpansion)

      } // END loop across voxels



   // Update values in m_LabelImage
    for (label_it.GoToBegin(), tmp_it.GoToBegin(); !label_it.IsAtEnd(); ++label_it, ++tmp_it)
    {
     label_it.Set(tmp_it.Get());
    }

    // Fix labeling where the topology is violated
    bool test;
    if (CheckTopology&&false)
      {
	std::cout << "SHOULD NOT BE GETTING HERE RIGHT NOW!!!!!!!!" << std::endl;
	std::cout << "    Checking Topology" << std::endl;
	MadeChanges = false;
	int numViolations=0;
	int numValid = 0;
	for (label_it.GoToBegin(), tmp_it.GoToBegin(); !label_it.IsAtEnd(); ++label_it, ++tmp_it) {
	  typename ImageType::IndexType ind = tmp_it.GetIndex();   
	  if (label_it.Get()==newLayerLabel){
	    //std::cout << "Testing voxel " << ind << std::endl;
	    test = CheckLocalTopology18(ind,baseLabel,newLayerLabel);
	    //std::cout << "Voxel: " << ind << " of type " << label_it.Get() << " tested as " << test << std::endl;
	    //std::cin.get();
	    if (test) { 
	      MadeChanges = true; 
	      numValid++;
	    }
	    else {
	      label_it.Set(OldPixel);
	      numViolations++;
	    }
	  }
	}
	std::cout << "      Found " << numValid << " valid additions and " << numViolations << " topology violations" << std::endl;
      }


 return MadeChanges;   

}

template<unsigned int Dimension, class PixelType>
template< class ImgType1, class ImgType2, class PixType, class RegType> 
void crlImageLaplace<Dimension, PixelType>::InitializeImage(ImgType1 ImgToInitialize, ImgType2 ImgBase, PixType initialValue, RegType Region)
{
   // Allocate the thickness image to the size of the input image.

  ImgToInitialize->SetSpacing( ImgBase->GetSpacing() );
  ImgToInitialize->SetOrigin ( ImgBase->GetOrigin() );
  ImgToInitialize->SetDirection( ImgBase->GetDirection() );

  Region.SetSize( ImgBase->GetLargestPossibleRegion().GetSize());
  Region.SetIndex( ImgBase->GetLargestPossibleRegion().GetIndex());
  ImgToInitialize->SetRegions( Region );
  ImgToInitialize->Allocate(); 
  ImgToInitialize->FillBuffer( initialValue );
}

template<unsigned int Dimension, class PixelType>
void crlImageLaplace<Dimension, PixelType>::ComputeThicknessesStreamlines()
{
  IteratorType it( m_LabelImage    ,    m_LapSolImage->GetLargestPossibleRegion() );
  IteratorType lit( m_ThicknessImage, m_ThicknessImage->GetLargestPossibleRegion() );

  double tmpThick;

  typename ImageType::IndexType index;
  
  // Compute and store a displacement field representing the displacement to
  // move a particle from any point in the solution region to the end of the
  // streamline.
  for (lit.GoToBegin(), it.GoToBegin(); !lit.IsAtEnd(); ++lit, ++it) {
    // We only calculate in the new cortical region.
    PixelType label = it.Get();
    index = it.GetIndex();
    if ((label != m_CorticalLabelNew)&&( label != m_CorticalLabelChecked)) continue;
    // std::cout << "Computing Streamline For Pixel: " << index << std::endl;
    VectorPixelType displacement;
    
    itk::Point<double, Dimension> start;
    std::vector< itk::Point<double, Dimension> > ysave(0);
    m_LapSolImage->TransformIndexToPhysicalPoint(it.GetIndex(), start);
    ComputeStreamLine(start, ysave);

    tmpThick = 0;
    if (ysave.size() > 0) {
      for (unsigned int i = 0; i < Dimension; i++) {
	displacement[i] = ysave[ysave.size() - 1][i] - start[i];
	tmpThick += displacement[i]*displacement[i];
      }
      tmpThick = sqrt(tmpThick);
    }
    m_DisplacementFieldImage->SetPixel(it.GetIndex(), displacement);

    lit.Set(tmpThick);
  }
}


template<unsigned int Dimension, class PixelType>
void crlImageLaplace<Dimension, PixelType>::CheckThicknesses(double ExpectedThickness, double ExpectedThickness_old, bool diagExpansion)
{
  //  IteratorType thickness_it(m_ThicknessImage , m_ThicknessImage->GetLargestPossibleRegion() ); 
  IteratorType  thickness_it(m_ThicknessImage, m_DisplacementFieldImage->GetLargestPossibleRegion() ); 
  IteratorType      label_it(m_LabelImage     , m_LabelImage->GetLargestPossibleRegion()    );
  
  std::cout << "    Expected Thickness is: " << ExpectedThickness << std::endl;

  unsigned int numSulcal = 0;
  unsigned int numCortex = 0;
  unsigned int numChecked = 0;
  unsigned int numOverall = 0;
  double meanThickness = 0;
  double meanThicknessNew = 0;
  double meanThicknessChecked = 0;
  typename ImageType::PixelType TestPix;

  // Get Mean Thickness of Voxels that Have already been confirmed as grey matter
  for (thickness_it.GoToBegin(), label_it.GoToBegin(); !thickness_it.IsAtEnd() ; ++thickness_it, ++label_it)
    {
      if (label_it.Get() != m_CorticalLabelChecked) { continue; }

      meanThicknessChecked += thickness_it.Get();
      numChecked++;
    }
  meanThicknessChecked = meanThicknessChecked/ numChecked;



  for (thickness_it.GoToBegin(), label_it.GoToBegin() ;  !thickness_it.IsAtEnd() ; ++thickness_it, ++label_it)
    {
      
      TestPix = thickness_it.Get();
      meanThickness += TestPix;
      if (TestPix>0) {numOverall++;}

     
      if (label_it.Get() != m_CorticalLabelNew) { continue; }

      meanThicknessNew += TestPix;
      if ((TestPix>ExpectedThickness)) // ||(TestPix<ExpectedThickness_old))
	{ // If it's a newly added cortical voxel, and the thickness is higher than expected, label it sulcal
	  label_it.Set(m_SulcalLabel);
          numSulcal++;
	}
      else
        { // If it's newly added cortex, and the thickness is within the expected range, label it cortex.
	  if (m_CheckTopology)
	    label_it.Set(m_TempLabel);
	  else
	    label_it.Set(m_CorticalLabelChecked);
          numCortex++;
        }      
    }



   // Fix labeling where the topology is violated
    bool test;
    if (m_CheckTopology)
      {
	IteratorType   tmp_it( m_LabelImageTmp , m_LabelImageTmp->GetLargestPossibleRegion() );
	// Update values in m_LabelImageTmp
	for (label_it.GoToBegin(), tmp_it.GoToBegin(); !label_it.IsAtEnd(); ++label_it, ++tmp_it)
	  {
	    tmp_it.Set(label_it.Get());
	  }
	std::cout << "    Checking Topology" << std::endl;
	int numViolations=0;
	int numValid = 0;
	for (label_it.GoToBegin(),tmp_it.GoToBegin(); !label_it.IsAtEnd(); ++label_it,++tmp_it) {
	  typename ImageType::IndexType ind = label_it.GetIndex();   
	  if (label_it.Get()==m_TempLabel){
	    //std::cout << "Testing voxel " << ind << std::endl;
	    if (diagExpansion)
	      test = CheckLocalTopology18(ind,m_CorticalLabelChecked,m_TempLabel);
	    else
	      test = CheckLocalTopology6(ind,m_CorticalLabelChecked,m_TempLabel);
	    
	    //test = CheckLocalTopology18(ind,m_CorticalLabelChecked,m_TempLabel);

	    //std::cout << "Voxel: " << ind << " of type " << label_it.Get() << " tested as " << test << std::endl;
	    //std::cin.get();

            if (((int)ind[0]==127)&((int)ind[1]==119)&((int)ind[2]==129)){
	      std::cout << "Test Returned " << test << std::endl;
	    }

	    if (test) { 
	      label_it.Set(m_CorticalLabelChecked);
	      numValid++;
	    }
	    else {
	      if (m_ExpandSulci)
		// label_it.Set(m_SolveRegionLabel);
		label_it.Set(m_SulcalLabelTopCheck);
	      else
		label_it.Set(m_SulcalLabelTopCheck);
	      numViolations++;
	    }
	  }
	}
	std::cout << "      Found " << numValid << " valid additions and " << numViolations << " topology violations" << std::endl;
      }



  meanThickness    = meanThickness/numOverall;
  meanThicknessNew = meanThicknessNew/(numSulcal + numCortex);
  

  std::cout << "    Average cortical thickness at this stage: " << meanThickness << std::endl;
  std::cout << "    Average thickness for new cortical voxels: " << meanThicknessNew << std::endl;
  std::cout << "    Average thickness for tested cortical voxels: " << meanThicknessChecked << std::endl;
  std::cout << "    Confirmed " << numCortex << " cortical voxels and identified " << numSulcal << " sulcal voxels." << std::endl;

}


template<unsigned int Dimension, class PixelType>
void crlImageLaplace<Dimension, PixelType>::SaveSol(int whichFile, int numLayersAdded)
{
  // 1: Laplacian
  // 2: Segmentation
  // 3: Thickness
  // 4: Tangent Field
  std::string fName;
  switch (whichFile){
  case 1 : m_Writer->SetInput( m_LapSolImage );
    fName = "LapSol" + stringify(numLayersAdded) + ".nrrd";
    std::cout << "  Writing intermediate solution to Laplace's Equation as " << fName << std::endl;
    break;
  case 2:  m_Writer->SetInput( m_LabelImage );
    fName = "NewSeg" + stringify(numLayersAdded) + ".nrrd";
    std::cout << "  Writing intermediate segmentation as " << fName << std::endl;
    break;
  case 3:  m_Writer->SetInput( m_ThicknessImage );
    fName = "ThickSol" + stringify(numLayersAdded) + ".nrrd";    
    std::cout << "  Writing intermediate Thickness measures as " << fName << std::endl;
    break;
  case 4:  m_TangentFieldWriter->SetInput( m_TangentImage );
    fName = "TanSol" + stringify(numLayersAdded) + ".nrrd";    
    std::cout << "  Writing intermediate Tangent Field as " << fName << std::endl;
    break;
  case 5:  m_DisplacementFieldWriter->SetInput( m_DisplacementFieldImage );
    fName = "DispSol" + stringify(numLayersAdded) + ".nrrd";    
    std::cout << "  Writing intermediate Displacement Field as " << fName << std::endl;  
    break;
  }


  switch (whichFile){
  case 1:
  case 2:
  case 3:
      m_Writer->SetFileName(fName);
      m_Writer->UseCompressionOn();
      try {
	m_Writer->Update();
      } 
      catch (itk::ExceptionObject& err ) {
	std::cerr << "error writing output file: " << fName 
		  << ": " << err << std::endl;
      }
      break;
  case 4:
      m_TangentFieldWriter->SetFileName( fName);
      m_TangentFieldWriter->UseCompressionOn();
      
      try {
	m_TangentFieldWriter->Update();
      } catch (itk::ExceptionObject& err ) {
	std::cerr << "error writing output file: " << m_OutputTangentFieldFileName 
		  << ": " << err << std::endl;
      }
      break;

  case 5:      
      m_DisplacementFieldWriter->SetFileName( fName );
      m_DisplacementFieldWriter->UseCompressionOn();
      
      try {
	m_DisplacementFieldWriter->Update();
      } catch (itk::ExceptionObject& err ) {
	std::cerr << "error writing output file: " << 
	  m_OutputDisplacementFieldFileName << err << std::endl;
	}
      break;
  }
}



//
// BEGIN crlImageLaplace::Execute()
//

template <unsigned int Dimension, class PixelType>
int crlImageLaplace<Dimension, PixelType >::Execute()
{

  // Set Input and Output Files
  m_Reader->SetFileName( m_LabelImageFileName );
  m_Writer->SetFileName( m_OutputLaplaceFileName);
  m_Writer->UseCompressionOff();

  // Try Reading Segmentation File
  try {
    m_Reader->Update();
  } catch( itk::ExceptionObject & exp )
  {
    std::cerr << "Failed to read image " << m_LabelImageFileName << std::endl;
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Successfully read segmentation image" << std::endl;

  // Set the m_LabelImage pointer to hold the label image:
  m_LabelImage = m_Reader->GetOutput();

  // Assign New Cortical Label and Sulcal Voxel Label
  IteratorType label_it(m_LabelImage,m_LabelImage->GetLargestPossibleRegion() );
  double MaxLabel = 0;
  for (label_it.GoToBegin(); !label_it.IsAtEnd()   ; ++label_it)
    {
      if (label_it.Get()>MaxLabel)
	{ MaxLabel = label_it.Get();}
    }
  m_CorticalLabelNew     = MaxLabel + 1;
  m_CorticalLabelChecked = MaxLabel + 2;
  m_SulcalLabel          = MaxLabel + 3;
  m_SulcalLabelTopCheck  = MaxLabel + 5;
  m_TempLabel            = MaxLabel + 4;
  
  std::cout << "Setting New Segment Labels  "                             << std::endl;
  std::cout << "         New Cortical Layers: " << m_CorticalLabelNew     << std::endl;
  std::cout << "     Checked Cortical Layers: " << m_CorticalLabelChecked << std::endl;
  std::cout << "               Sulcal Layers: " << m_SulcalLabel          << std::endl;
  std::cout << "            Temporary Layers: " << m_TempLabel            << std::endl;


  // Set Initial Pixel Values
  typename ImageType::PixelType  initialValue = 0.0;
  typename VectorImageType::PixelType initialVectorValue(0.0);

  typename ImageType::RegionType Region;
  typename VectorImageType::RegionType VecRegion;

  // Allocate the output image to the size of the input image.
  if (!m_LapSolImage) { return EXIT_FAILURE; }
  std::cout<< "Initializing Laplace Solution Image" << std::endl;
  InitializeImage(m_LapSolImage,m_LabelImage,initialValue,Region);
 
  // Allocate the thickness image to the size of the input image.
  if (!m_ThicknessImage) { return EXIT_FAILURE; }
  std::cout<< "Initializing Thickness Image" << std::endl;
  InitializeImage(m_ThicknessImage,m_LapSolImage,initialValue,Region);


  // Allocate the LabelImageTmp Image to the size of the input image
  if (!m_LabelImageTmp) {  return EXIT_FAILURE; }
  std::cout<< "Initializing Temporary Label Image" << std::endl;
  InitializeImage(m_LabelImageTmp,m_LapSolImage,initialValue,Region);

  // Allocate the L0 image to the size of the input image.
  if (!m_L0Image) { return EXIT_FAILURE; }
  std::cout<< "Initializing L0 Image" << std::endl;
  InitializeImage(m_L0Image,m_LapSolImage,initialValue,Region);

  // Allocate the L1 image to the size of the input image.
  if (!m_L1Image) { return EXIT_FAILURE; }
  std::cout << "Initializing L1 Image" << std::endl;
  InitializeImage(m_L1Image,m_LapSolImage,initialValue,Region);

  // Allocate the Output Tangent Image to the size of the input image
  if (!m_TangentImage) 
   { return EXIT_FAILURE; }
  InitializeImage(m_TangentImage,m_LapSolImage,initialVectorValue,VecRegion);

  // Ensure the vector interpolator is hooked up to the tangent image.
  m_VecInterpolator->SetInputImage( m_TangentImage );

  // Allocate the Output Displacement Field Image to the size of the input image
  if (!m_DisplacementFieldImage) { return EXIT_FAILURE; }
  InitializeImage(m_DisplacementFieldImage,m_LapSolImage,initialVectorValue,VecRegion);

  // Compute Expected Thicknesses
  std::cout << "Computing expected thicknesses for each cortical layer" << std::endl;
  double BaseThickness  = 0.0;
  double BaseThickness2 = 0.0;

  unsigned int i;
  for (i=0 ; i<Dimension ; i++)
    {
      BaseThickness += m_LabelImage->GetSpacing()[i]*m_LabelImage->GetSpacing()[i];
      if (m_LabelImage->GetSpacing()[i] > BaseThickness2)
        { BaseThickness2 = m_LabelImage->GetSpacing()[i];
	}
    } 
  BaseThickness *= 0.6666; // We only ever expand in two dimensions at a time
  BaseThickness = sqrt(BaseThickness);
  //  BaseThickness2 = BaseThickness; // If we want both thickness to be the same

  std::cout << " Thickness #1: " << BaseThickness << " Thickness #2: " << BaseThickness2 << std::endl;


  // Add a 1 voxel thick layer off gray matter onto the white matter
  // Assume that this first layer is always gray matter (ie - don't do a thickness check).
  bool wasUpdated = 0;
  bool foo = 0;
  bool useDiagExpansion = 1;
  int numLayersAdded = 0;
  double ExpectedThickness = BaseThickness2; // BaseThickness2;
  double ExpectedThickness_old = 0;

  // Adding initial unchecked layer of cortex
  if (m_InitialCorticalLayers>0)
    {
      std::cout << "Adding Initial Cortical Gray Matter Layer... ";
      wasUpdated = AddLayer(m_InnerBorderLabel,m_SolveRegionLabel,m_CorticalLabelChecked,useDiagExpansion,false);
      if (wasUpdated) {
	std::cout << "... Layer Successfully Added" << std::endl;
	if (useDiagExpansion){
	  ExpectedThickness += BaseThickness;}
	else { 
	  ExpectedThickness += BaseThickness2;} 
      }  
      m_InitialCorticalLayers--;
      numLayersAdded++;
    }

  SaveSol(2,numLayersAdded);

  // Add additional, checked layers of graymatter
  while (m_InitialCorticalLayers>0)
      {
	std::cout << m_InitialCorticalLayers << " Layers Remain to be Added ...";
	useDiagExpansion = !useDiagExpansion;
	wasUpdated = AddLayer(m_CorticalLabelChecked,m_SolveRegionLabel,m_CorticalLabelChecked,useDiagExpansion,false);
	if (wasUpdated) {
	  std::cout << " ... Layer Successfully Added" << std::endl;
	  if (useDiagExpansion){
	    ExpectedThickness += BaseThickness;}
	  else { 
	    ExpectedThickness += BaseThickness2;} 	  
	}  
        m_InitialCorticalLayers--;
	numLayersAdded++;
	SaveSol(2,numLayersAdded);
      }
      
 

  /*
  // Solve Laplace's Equation on the cortical layer
  std::cout << "  Solving Laplace's Equation On New Layer" << std::endl;
  SolveLaplace();      
  SaveSol(1,numLayersAdded);
  
  // Compute Tangent Field
  std::cout << "  Computing Tangent Field to Streamlines" << std::endl;
  ComputeGradient();
  SaveSol(4,numLayersAdded);
  
  
  // Compute Thickness Using Yezzi-Prince Approach
  std::cout << "  Computing Thickness Using Yezzi-Prince Approach" << std::endl;
  SolveYezziPrince();
           
   // Compute Displacement Vectors By Tracing Streamlines
     std::cout << "Computing Displacement Image" << std::endl;
     ComputeThicknessesStreamlines();
     SaveSol(5,numLayersAdded);
  
  

  // Write out intermediate tangent thickness solution
  SaveSol(3,numLayersAdded);

         
  // Check the Thicknesses
  std::cout << "  Checking Thicknesses and labeling sulcal voxels" << std::endl;
  CheckThicknesses(ExpectedThickness);
  */

  std::cout << "Starting iterative addition of cortical layers" << std::endl;
  int numIterations = 0;
  while (wasUpdated&&(++numIterations<m_MaxIterations)){

    std::cout << std::endl << " Starting Iteration " << numLayersAdded << std::endl;
    
    // Add an additional new layer of graymatter to existing material.
    std::cout << "  Adding Graymatter Layer Number " << ++numLayersAdded << std::endl;

    useDiagExpansion = !useDiagExpansion; // Use diagonal expansion every other iteration
    wasUpdated = AddLayer(m_CorticalLabelChecked,m_SolveRegionLabel,m_CorticalLabelNew, useDiagExpansion,m_CheckTopology);

    if (m_ExpandSulci)
      foo = AddLayer(m_SulcalLabel,m_SolveRegionLabel,m_SulcalLabel,useDiagExpansion,false);
    
    if (wasUpdated) {
      ExpectedThickness_old = ExpectedThickness;
      if (useDiagExpansion){        
	ExpectedThickness += BaseThickness;}
      else { 
	ExpectedThickness += BaseThickness2;} 

      SaveSol(2,numLayersAdded);

      // Solve Laplace's Equation on the cortical layer
      std::cout << "  Solving Laplace's Equation On New Layer" << std::endl;
      SolveLaplace();
      SaveSol(1,numLayersAdded);

      
      // Compute Tangent Field
      std::cout << "  Computing Tangent Field to Streamlines" << std::endl;
      ComputeGradient();
      SaveSol(4,numLayersAdded);

      // Compute Thickness Using Yezzi-Prince Approach
      std::cout << "  Computing Thickness Using Yezzi-Prince Approach" << std::endl;
      SolveYezziPrince();
           
      /* // Compute Displacement Vectors By Tracing Streamlines
       std::cout << "Computing Displacement Image" << std::endl;
      ComputeThicknessesStreamlines();
      SaveSol(5,numLayersAdded);
       */

      SaveSol(3,numLayersAdded);

         
      // Check the Thicknesses
      std::cout << "  Checking Thicknesses and labeling sulcal voxels" << std::endl;
      CheckThicknesses(ExpectedThickness,ExpectedThickness_old,useDiagExpansion);

      
    }  
   


  }  // END WHILE LOOP


  for (label_it.GoToBegin(); !label_it.IsAtEnd()   ; ++label_it)
    {
      if (label_it.Get()==m_SolveRegionLabel)
	{ label_it.Set(m_SulcalLabel);}
    }


  // Solve things one last time to get the actual cortical thicknesses, excluding sulci
  m_ExcludeSulci = true;
  SolveLaplace();
  ComputeGradient();
  SolveYezziPrince();


  // Write out the final segmentation
    m_Writer->SetInput( m_LabelImage );
    m_Writer->SetFileName(m_OutputSegFileName);
    std::cout << "  Writing final segmentation as " << m_OutputSegFileName << std::endl;
    try {
      m_Writer->Update();
    } 
    catch (itk::ExceptionObject& err ) {
      std::cerr << "error writing output file: " << m_OutputSegFileName
		<< ": " << err << std::endl;
      return 1;
    }
 

    // Write out final thickness solution
    std::cout << "  Writing final Thickness measures as " << m_OutputThicknessFileName << std::endl;
    m_Writer->SetFileName( m_OutputThicknessFileName);
    m_Writer->SetInput( m_ThicknessImage );
    try {
      m_Writer->Update();
    } catch (itk::ExceptionObject& err ) {
      std::cerr << "error writing output file: " << m_OutputThicknessFileName 
		<< ": " << err << std::endl;
      return 1;
    }

  // Write out m_LapSolImage since it has the updated values.
  m_Writer->SetInput( m_LapSolImage );
  m_Writer->SetFileName( m_OutputLaplaceFileName);
  try {
    m_Writer->Update();
  } catch (itk::ExceptionObject& err ) {
    std::cerr << "error writing output file: " << m_OutputLaplaceFileName 
      << ": " << err << std::endl;
    return 1;
  }

  // Write out final tangent image
  m_TangentFieldWriter->SetFileName( m_OutputTangentFieldFileName);
  m_TangentFieldWriter->UseCompressionOff();
  m_TangentFieldWriter->SetInput( m_TangentImage );
  try {
    m_TangentFieldWriter->Update();
  } catch (itk::ExceptionObject& err ) {
    std::cerr << "error writing output file: " << m_OutputTangentFieldFileName 
      << ": " << err << std::endl;
    return 1;
  }
 


  /*
  m_ThicknessWriter->SetFileName( "L0Image.nrrd" );
  m_ThicknessWriter->UseCompressionOff();
  m_ThicknessWriter->SetInput( m_L0Image );

  try { m_ThicknessWriter->Update(); } 
  catch (itk::ExceptionObject& err ) 
  { std::cerr << "error writing output file: " << "L0Image "
      << ": " << err << std::endl;
    return 1;
  }

  m_ThicknessWriter->SetFileName( "L1Image.nrrd" );
  m_ThicknessWriter->UseCompressionOff();
  m_ThicknessWriter->SetInput( m_L1Image );
  try { m_ThicknessWriter->Update(); } 
  catch (itk::ExceptionObject& err ) 
  { std::cerr << "error writing output file: " << "L1Image "
      << ": " << err << std::endl;
    return 1;
  }
  */





  std::cout << "Successfully Completed Running crlImageLaplace_dh" << std::endl;
  return EXIT_SUCCESS;
}

//
//  END crlImageLaplace::Execute()
//



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
    std::vector< itk::Point<double, Dimension> > &ysave)
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
//    int count = 0;

    for (unsigned int i = 0; i < Dimension; i++) {
      y[i] = start[i];
    }
    unsigned int nsteps = 0;

    typename ImageType::IndexType index;
//    typename VectorImageType::IndexType vindex;
    double norm = 0.0;
    for (nsteps = 0; nsteps < MAXSTEPS; nsteps++) {
      if (!  this->Derivative(s,y,dydx) ) {
        // The derivative can fail to be calculatable when we step outside
        // the solution zone.
        // The norm of the tangent field gets smaller than one due to 
        // the linear interpolation when we approach close to the boundary.
	std::cout << "Derivative Incalcuable at " << y << std::endl;
	char hold;
	std::cin.get(hold);
        // Do nothing for now.
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
     if ((isinside) && (m_LabelImage->GetPixel(index) != m_CorticalLabelChecked) &&
          	  (m_LabelImage->GetPixel(index) != m_CorticalLabelNew)) {
        // stop integrating now that we are out of the solve region.
        // binary search for the exact location of the boundary we match to.
        // search between y and yout for the location of the boundary.
        // push the boundary location as the last step.
        double diffnorm = 0.0;
        for (unsigned int i = 0; i < Dimension; i++) {
          diffnorm += ( (vin[i] - vout[i])*(vin[i] - vout[i]) );
        }
        diffnorm = sqrt(diffnorm);
	int loopidx = 0;
        while (diffnorm > hmin) {
          for (unsigned int i = 0; i < Dimension; i++) {
            vmid[i] = 0.5*(vin[i] + vout[i]);
          }
          isinside = m_LabelImage->TransformPhysicalPointToIndex( vmid, index );
          if ( !isinside || 
               ((m_LabelImage->GetPixel( index ) != m_CorticalLabelChecked) &&
                (m_LabelImage->GetPixel( index ) != m_CorticalLabelNew ) &&
		(((m_LabelImage->GetPixel(index)!= m_TempLabel) && (m_LabelImage->GetPixel( index ) != m_SulcalLabel )) || m_ExcludeSulci))) {
            vout = vmid;
          } else {
            vin = vmid;
          }
          diffnorm = 0.0;
          for (unsigned int i = 0; i < Dimension; i++) {
            diffnorm += ( (vin[i] - vout[i])*(vin[i] - vout[i]) );
          }
          loopidx++;
          if (loopidx>100)
	    { std::cout<< "We seem to be stuck in a while loop with index: " << index << " and value " << m_LabelImage->GetPixel( index) << std::endl;}
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
//  bool retval = true;
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
    p[i] = y[i];   // Initialize Point Location
    dydx[i] = 0.0; // Initialize derivative to be zero
  }

  // Get the index corresponding to location p
  typename ImageType::IndexType lindex;
  typename ImageType::PixelType label;
  if (!m_LabelImage->TransformPhysicalPointToIndex( p, lindex )) {
    return false; // Cannot compute derivative.
  }

  // If we're not within the solution region, return false
  label = m_LabelImage->GetPixel( lindex );
  if (( m_LabelImage->GetPixel( lindex ) != m_CorticalLabelNew ) &&
      ( m_LabelImage->GetPixel( lindex ) != m_CorticalLabelChecked )&&
      ((( m_LabelImage->GetPixel( lindex ) != m_SulcalLabel) && (m_LabelImage->GetPixel(lindex)!= m_TempLabel)) || m_ExcludeSulci)) {
    // std::cout << "Asking for derivative outside of the solve region. Point index is: " << lindex << std::endl;
    // std::cout << " Label is:" << m_LabelImage->GetPixel(lindex) << " instead of " << m_CorticalLabelNew << " or " << m_CorticalLabelChecked << " or " << m_SulcalLabel << std::endl;
    // std::cout << "Reduce the step size and try again..." << std::endl;
    return false; // Cannot compute derivative
  }

  // Get Interpolated Gradient Vector
  if ( m_VecInterpolator->IsInsideBuffer( p ) ) {
     value = m_VecInterpolator->Evaluate( p );
  } else {
    // Can't interpolate a value that is outside the image.
    return false; // Cannot compute derivative
  }

  // Assign Output
  for (unsigned int i = 0; i < Dimension; i++) {
    dydx[i] = value[i];
  }

  return true; // Successful derivative calculation.

};

