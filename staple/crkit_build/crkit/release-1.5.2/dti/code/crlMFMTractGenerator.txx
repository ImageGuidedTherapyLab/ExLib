
#ifndef _CRLMFMTRACTGENERATOR_TXX_INCLUDED
#define _CRLMFMTRACTGENERATOR_TXX_INCLUDED 1

#include <vnl/vnl_vector.h>
#include <vnl/vnl_cross.h>
#include <vnl/vnl_math.h>

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkContinuousIndex.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "TensorLogFunctor.h"
#include "TensorExpFunctor.h"

#include "crlMFMUtils.h"
#include "crlVtkMeshIO.h"
#include "crlTaskProgress.h"

#include <vtkGIFTIWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPointData.h>

// #include "crlTensorLinearInterpolateImageFunction.h"

/**
 * This finds the tensor in the input vector that is oriented most
 * similarly to the input tensor.
 * Returns the index of the corresponding tensor.
 */
template <const unsigned int TImageDimension, class TComponentType>
int
crl::MFMTractGenerator<TImageDimension, TComponentType>
::SelectMatchingTensor(
        const typename ImageType::PixelType &in,
        const typename ImageType::IndexType &index,
        itk::Vector< TComponentType, 3> &evec, // output e vec.
        typename ImageType::PixelType &out) const   // output tensor.
{
  //unsigned int m_TensorIndex = 0; //not necessary
  itk::Vector< TComponentType, 3> refEvec;
  itk::Vector< TComponentType, 3> compareEvec;
  typename ImageType::PixelType compareTensor;
  double projectionLength = -1.0;
  double compareProjectionLength = 0.0;
  int idTensor=-1;
  /* Makes use of member variable m_TensorImageVector. */

  EstimatePrimaryEigenvector(in, refEvec);

  // rq: if one tensor is invalid the compareProjectionLength will be null, so OK.
  //     if all the tensors are invalid, return a null tensor and null vector.
  for (unsigned int i = 0; i < m_TensorImageVector.size(); i++) {
    compareTensor = m_TensorImageVector[i]->GetPixel(index);
    EstimatePrimaryEigenvector(compareTensor, compareEvec);
    compareProjectionLength = fabs(compareEvec * refEvec);
    if (compareProjectionLength > projectionLength) {
      projectionLength = compareProjectionLength;
      idTensor = (int)i;
      out = compareTensor;
      evec = compareEvec;
    }
  }

  return idTensor;
}

/**
 * This carries out Euclidean interpolation of the input tensor field,
 * selecting at each corner of the voxel the tensor in the input vector that 
 * is oriented most similarly to the input tensor.
 */
template <const unsigned int TImageDimension, class TComponentType>
void
crl::MFMTractGenerator<TImageDimension, TComponentType>
::InterpolateMatchingTensor(
							const typename ImageType::PixelType &in,
							const typename TractType::ContinuousIndexType &cindex,
							itk::Vector< TComponentType, 3> &evec, // output e vec.
							typename ImageType::PixelType &out) const   // output tensor.
{
	if ( m_TensorInterpolationMode!=INTERP_NONE )
	{

		// here we implement a Euclidean interpolation.
		// The tensors at the corners of the voxel are selected based on similarity
		// of their primary eigenvector to the input tensor direction.
		itk::Functor::TensorExp<typename ImageType::PixelType, 
			typename ImageType::PixelType> tensorExp;
		itk::Functor::TensorLog<typename ImageType::PixelType, 
			typename ImageType::PixelType> tensorLog;


		unsigned int dim;  // index over dimension

		/**
		* Compute base index = closest index below point
		* Compute distance from point to base index (from the point to the voxel center)
		*/
		signed long baseIndex[TImageDimension];
		double distance[TImageDimension];

		for( dim = 0; dim < TImageDimension; dim++ )
		{
			baseIndex[dim] = long (vnl_math_rnd( cindex[dim] ));
			distance[dim] = cindex[dim] - double( baseIndex[dim] ) + 0.5;
		}

		/**
		* Interpolated value is the weight some of each of the surrounding
		* neighbors. The weight for each neighbour is the fraction overlap
		* of the neighbor pixel with respect to a pixel centered on point.
		*/
		typename ImageType::PixelType output;
		output.Fill( 0.0 );

		double totalOverlap = 0.0;

		for( unsigned int counter = 0; counter < (1<<TImageDimension); counter++ )
		{

			double overlap = 1.0;          // fraction overlap
			unsigned int upper = counter;  // each bit indicates upper/lower neighbour
			typename ImageType::IndexType neighIndex;

			// get neighbor index and overlap fraction
			for( dim = 0; dim < TImageDimension; dim++ )
			{

				if ( upper & 1 )
				{
					neighIndex[dim] = baseIndex[dim] + 1;
					overlap *= distance[dim];
				}
				else
				{
					neighIndex[dim] = baseIndex[dim];
					overlap *= 1.0 - distance[dim];
				}

				upper >>= 1;

			}

			// get neighbor value only if overlap is not zero
			if( overlap )
			{
				// Skip illegal memory accesses - same as setting outside to be zero.
				if (m_TensorImageVector[0]->GetLargestPossibleRegion().IsInside( neighIndex )) {
					typename ImageType::PixelType cornerTensor;
					SelectMatchingTensor(in, neighIndex, evec, cornerTensor);

					// log-Euclidean:
					if ( m_TensorInterpolationMode==INTERP_LOGEUCLIDEAN ) {
						typename ImageType::PixelType logCornerTensor;
						logCornerTensor = tensorLog( cornerTensor );
						for(unsigned int k = 0; k < ImageType::PixelType::InternalDimension; k++ ) 
							output[k] += overlap * static_cast<double>( logCornerTensor[k] );
					}
					else {
						for(unsigned int k = 0; k < ImageType::PixelType::InternalDimension; k++ ) 
							output[k] += overlap * static_cast<double>( cornerTensor[k] );		  
					}
				}
				totalOverlap += overlap;
			}

			if( totalOverlap == 1.0 )
			{
				// finished
				break;
			}
		}

		// Now set the output
		if ( m_TensorInterpolationMode==INTERP_LOGEUCLIDEAN )
			out = tensorExp( output );
		else
			out = output;
	}

	// NO INTERPOLATION : m_TensorInterpolationMode == INTERP_NONE
	else
	{
		// Convert the continuous index to an index
		typename ImageType::IndexType index;
		index.CopyWithRound( cindex );

		// Select the matching tensor for that index
		SelectMatchingTensor(in, index, evec, out);
	}

	EstimatePrimaryEigenvector(out, evec);
}

template <const unsigned int TImageDimension, class TComponentType>
void crl::MFMTractGenerator<TImageDimension, TComponentType>
::EstimatePrimaryEigenvector(const typename ImageType::PixelType &in, 
                               itk::Vector< TComponentType, 3> &evec) const
{
  typedef typename ImageType::PixelType::ValueType ValueType;
  typedef typename ImageType::PixelType::EigenVectorsMatrixType MatrixType;
  typedef typename ImageType::PixelType::EigenValuesArrayType EigenValuesType;

  if ( in.GetTrace()<1e-6 )
  {
    evec[0] = 0;
    evec[1] = 0;
    evec[2] = 0;
  }
  else
  {
    MatrixType V;
    EigenValuesType e;
    
    // The vnl convention is for the eigenvalues to be ordered from
    // the smallest to the largest, and for the eigenvectors to be in the rows.
    in.ComputeEigenAnalysis(e, V);

    // std::cout << "Input tensor is " << in << std::endl;
    // std::cout << "eigenvector matrix is " << V << std::endl;

    evec[0] = V[2][0];
    evec[1] = V[2][1];
    evec[2] = V[2][2];
  }
}
 
template <const unsigned int TImageDimension, class TComponentType>
void crl::MFMTractGenerator<TImageDimension, TComponentType>
::GetFractionValues(const typename ImageType::IndexType &index, std::vector<float> &fractions ) const
{
	// Prepare the index for the fraction image
	//unsigned int tensorID;
	//FractionImageType::IndexType findex;
	//for ( tensorID=0; tensorID<TImageDimension ; tensorID++ ) findex[tensorID]=index[tensorID];

	// Grab all the tensor fraction for that index
	fractions.resize(m_TensorImageVector.size());
	for ( unsigned int t=0; t<fractions.size() ; t++ )
	{
		if ( m_FractionImage.IsNull() )
			fractions[t] = 1.0/m_TensorImageVector.size();
		else {
			
			//findex[tensorID] = t;
			fractions[t] = m_FractionImage->GetPixel(index)[t];
		}
	}
}

// Randomly select a tensor in proportion to the fraction of the voxel that each tensor occupies.
template <const unsigned int TImageDimension, class TComponentType>
unsigned int crl::MFMTractGenerator<TImageDimension, TComponentType>
::SelectTensorDependingOnFractions(const typename TractType::ContinuousIndexType& cindex) const
{
	// If the fraction image was not defined just select the tensor randomly
	if ( m_FractionImage.IsNull() ) {
		double tmp = m_RandGen->GetVariateWithOpenUpperRange() * m_TensorImageVector.size();
		return ( static_cast<unsigned int>( floor( tmp ) ) );
	}
	// If there are some fractions, use them to select the tensor;
	else
	{
		// Take the fractions from the cindex
		std::vector<float> fractions;
		typename ImageType::IndexType idx;
		idx.CopyWithRound( cindex );
		this->GetFractionValues(idx, fractions);

		// Compute the sum of the tensor fractions 
		// (only the tensor fractions, do not take into account a possible isotropic part fraction)
		float sumF = 0.0f;
		for ( unsigned int i=0; i<m_TensorImageVector.size() ; i++ )
			sumF += fractions[i];

		// If no tensor, just return the first index
		if ( sumF==0 ) return 0;

		// Select the tensor according to the fractions
		int outputID = -1;
		double tmp = m_RandGen->GetVariateWithOpenUpperRange() * sumF;
		for ( unsigned int i=0; i<m_TensorImageVector.size() ; i++ )
		{
			if ( tmp < fractions[i] ) {
				outputID = i;
				break;
			}
			else tmp -= fractions[i];
		}

		if ( outputID==-1 ) {
			std::cout << "Strange error. Cannot determine the tensor from the fractions. Select id=0" << std::endl;
			return 0;
		}
		else return static_cast<unsigned int>(outputID);

	}
}

/**********************************************************************************************//**
 * \brief	Estimate the tract from a seeding point. 
 *
 *			!!!!!!! Modification by benoit: the output tract
 * 			points are not anymore a list of continuous indexes but a list of PHYSICAL
 * 			COORDINATES !!!!!!!!! (necessary to simplify the tract in crlMFMWBConnectivity)
*************************************************************************************************/
template <const unsigned int TImageDimension, class TComponentType>
void crl::MFMTractGenerator<TImageDimension, TComponentType>
::EstimateTract( const typename ImageType::PointType &point, 
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

  // It might be better to aim to use about half the minimum size of any edge of
  // any voxel, so we get two or more steps inside each voxel.
  double stepSize = fmin(fmin(m_ROIImage->GetSpacing()[0], 
                              m_ROIImage->GetSpacing()[1]), 
                           m_ROIImage->GetSpacing()[2]);
  stepSize /= static_cast<double>(m_StepsPerVoxel);

  //std::cout << "Point to point step size is " << stepSize << std::endl;

  typename ImageType::PointType newpoint = point;
  typename ImageType::PixelType tensor;
  typename ImageType::PixelType logtensor;
  typename ImageType::PixelType initialTensor;
  itk::Vector< TComponentType, 3> currentEigenvector;
  itk::Vector< TComponentType, 3> previousEigenvector;
  itk::Vector< TComponentType, 3> tendDeflection;  // for: D v
  itk::Vector< TComponentType, 3> tendDeflection2; // for: D^2 v
  itk::Vector< TComponentType, 3> initialDirection;

  previousEigenvector.Fill(0.0);
  tendDeflection.Fill(0.0);
  tendDeflection2.Fill(0.0);
  initialDirection.Fill(0.0);

  typename TractType::ContinuousIndexType cindex;
  // Used for NN interpolation: 
  typename ImageType::IndexType dindex;
  bool isInside;
  bool forwardDirection = true;
  double projectionLength = 0.0;
  double previousAngleCriterion = 1.0;
  double fractionalAnisotropyM = 1.0;
  unsigned int stepCount = 0;

  // We first step along in one direction,
  // then restart along the opposite direction.

  // First we check the start point.

  isInside = m_TensorImageVector[0]->TransformPhysicalPointToContinuousIndex(
       newpoint, cindex );
  // We can't make tracts from a point not inside the volume.
  if (!isInside) return;

  //std::cout << "Estimating a tract with cindex " << cindex << "." <<std::endl;

  //MFM: Choose which tensor we interpolate.
  // Since we are selecting the starting tract, we will choose it at random,
  // in function of the fraction of occupancies.
  m_TensorIndex = SelectTensorDependingOnFractions(cindex);
 
  // std::cout << "Choosing from tensor model of size " << m_TensorImageVector.size() << std::endl;
  // std::cout << "Stochastically selected to start with tensor " << m_TensorIndex << std::endl;
  isInside = m_TensorImageVector[m_TensorIndex]->TransformPhysicalPointToIndex(
       newpoint, dindex );
  if (!isInside) {
    std::cout << "cindex is inside but dindex is outside." << std::endl;
  }
  tensor = m_TensorImageVector[m_TensorIndex]->GetPixel(dindex);

  // Now that we have selected the tensor at the middle that we want,
  // we interpolate the tensor to subvoxel accuracy.
  InterpolateMatchingTensor(tensor, cindex, currentEigenvector, initialTensor);
  tensor = initialTensor;

  //std::cout << "primary evec is " << currentEigenvector << std::endl;

  // Save the first step direction, so later we can go the other way.
  initialDirection = currentEigenvector;
  // To initialize the tensor deflection calculation:
  previousEigenvector = currentEigenvector; 

  // First we work on the forward direction.
  while (isInside && (stepCount < m_MaxTractLength)) {
    // Test the stopping criteria for the tensor at this location.
    fractionalAnisotropyM = 
              fractionalAnisotropyM * m_FractionalAnisotropyMomentum + 
        (1 - m_FractionalAnisotropyMomentum)*tensor.GetFractionalAnisotropy();
    if (fractionalAnisotropyM < m_MinimumFA) {
      /* std::cout << "Terminating tract due to low FA. " <<
       fractionalAnisotropyM << " < " << m_MinimumFA << std::endl; */
      break;
    }

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
        /* std::cout << "Terminating tract since step exceeds maximum angle." 
             << previousAngleCriterion << " < " << m_MaximumAngle
                     << std::endl; */
        //std::cout << "fabs(projectionLength) is " << fabs(projectionLength) << std::endl;
        break;
    }

    // The current point is inside the image and has acceptable FA, so 
    // we will add it to the tract. The tract stores the PHYSICAL COORDINATES of
    // each of the points that we add.
	typename TractType::ContinuousIndexType pp;
	for ( unsigned int i=0; i<TImageDimension; i++ ) pp[i]=newpoint[i];
    if (forwardDirection) {
      forwardTract->AddVertex(pp);
    } else {
      backwardTract->AddVertex(pp);
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
    // Compute the new direction of travel based on combining the 
    // principal eigenvector and the tensor deflection direction,
    // and save it as the previousEigenvector.
    previousEigenvector = (previousEigenvector*m_DirectionMomentum +
      (tendDeflection*m_TendFraction + currentEigenvector*(1.0-m_TendFraction))
               *(1.0 - m_DirectionMomentum) );
    //Add it to the current position after scaling by the step size.
    // Utilize direction momentum in computing the new step.
    newpoint += ( previousEigenvector * stepSize );
    stepCount++;

    // Check if the newpoint is still inside the volume.
    isInside = m_TensorImageVector[0]->TransformPhysicalPointToContinuousIndex(
                   newpoint, cindex );
    if (isInside) {
      InterpolateMatchingTensor(tensor, cindex, currentEigenvector, tensor);
    }
  }

  // std::cout << "Working in the backward direction." << std::endl;

  // Let's do it again in the opposite direction.
  // The first point to consider is one step away from the initial start.
  forwardDirection = false;
  // If we could not start a tract in the forward direction,
  // and the initialDirection was not set, this will again stop at the
  // very first voxel, and the total tract length will be zero.
  previousEigenvector = initialDirection * -1.0;
  previousAngleCriterion = 1.0;
  fractionalAnisotropyM = 1.0;
  newpoint = point + previousEigenvector*stepSize;
  isInside = m_TensorImageVector[0]->TransformPhysicalPointToContinuousIndex(
       newpoint, cindex
       );

  if (isInside) {
    InterpolateMatchingTensor(initialTensor, cindex, 
                              currentEigenvector, tensor);
  }

  // Now we work on the backward direction.
  while (isInside && (stepCount < m_MaxTractLength)) {
    // Test the stopping criteria for the tensor at this location.
    fractionalAnisotropyM = 
              fractionalAnisotropyM * m_FractionalAnisotropyMomentum + 
        (1 - m_FractionalAnisotropyMomentum)*tensor.GetFractionalAnisotropy();
    if (fractionalAnisotropyM < m_MinimumFA) {
      /* std::cout << "Terminating tract due to low FA. " <<
       fractionalAnisotropyM << " < " << m_MinimumFA << std::endl; */
      break;
    }

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
        /* std::cout << "Terminating tract since step exceeds maximum angle." 
             << previousAngleCriterion << " < " << m_MaximumAngle
                     << std::endl; */
        // std::cout << "fabs(projectionLength) is " << fabs(projectionLength) << std::endl;
        break;
    }

    // The current point is inside the image and has acceptable FA, so 
    // we will add it to the tract. The tract stores the PHYSICAL COORDINATES of
    // each of the points that we add.
	typename TractType::ContinuousIndexType pp;
	for ( unsigned int i=0; i<TImageDimension; i++ ) pp[i]=newpoint[i];
    if (forwardDirection) {
      forwardTract->AddVertex(pp);
    } else {
      backwardTract->AddVertex(pp);
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
    // Compute the new direction of travel based on combining the 
    // principal eigenvector and the tensor deflection direction,
    // and save it as the previousEigenvector.
    previousEigenvector = (previousEigenvector*m_DirectionMomentum +
      (tendDeflection*m_TendFraction + currentEigenvector*(1.0-m_TendFraction))
               *(1.0 - m_DirectionMomentum) );
    //Add it to the current position after scaling by the step size.
    // Utilize direction momentum in computing the new step.
    newpoint += ( previousEigenvector * stepSize );
    stepCount++;

    // Check if the newpoint is still inside the volume.
    isInside = m_TensorImageVector[0]->TransformPhysicalPointToContinuousIndex(
                   newpoint, cindex );
    if (isInside) {
      InterpolateMatchingTensor(tensor, cindex, currentEigenvector, tensor);
    }
  }

  // Now we have the forward tract from the start point to its end,
  // and the backward tract from next to the start point to its end.
  // We need to iterate over the length of the paths to make one tract,
  // and return that.
  //   We go backwards from the end of the backward tract to its first point,
  // and then forwards over the forward tract from its start to its end.
  typename TractType::ContinuousIndexType pp;
  for (signed int i = backwardTract->GetVertexList()->Size() - 1; i >= 0; i--) {
    pp = backwardTract->GetVertexList()->GetElement(i);
    tract->AddVertex( pp );
  }
  for (unsigned int i = 0; i < forwardTract->GetVertexList()->Size(); i++) {
    pp = forwardTract->GetVertexList()->GetElement(i);
    tract->AddVertex( pp );
  }
  // We don't delete the forwardTract or reverseTract because the list
  // elements appear to be shared by reference.
}

/*
 Simplify a tract (angleMin in degrees)
 WARNING. The input tract points MUST be in physical coordinates to have valid angles.
 The code of EstimateTracts has been modified for that.
*/
template <const unsigned int TImageDimension, class TComponentType>
typename crl::MFMTractGenerator<TImageDimension, TComponentType>::TractType::Pointer  
crl::MFMTractGenerator<TImageDimension, TComponentType>::SimplifyTract ( double angleMin, typename TractType::Pointer &inTract )
{
	// Convert to radians
	double angleMinRad = angleMin*vnl_math::pi/180;

	//---------------------------------------------------------
	// Now loop on all points
	//---------------------------------------------------------
	typename TractType::ContinuousIndexType currentPoint, nextPoint, lastInsertedPoint;
	typename TractType::Pointer outTract = TractType::New();
	for (unsigned int i = 0; i < inTract->GetVertexList()->Size(); i++) 
	{
		currentPoint = inTract->GetVertexList()->GetElement(i);
		
		bool insertThisPoint = true;

		//----------------------------------------------------
		//keep the first two and last two points
		//----------------------------------------------------
		if ( i>1 && i<inTract->GetVertexList()->Size()-2 )	
		{
			nextPoint = inTract->GetVertexList()->GetElement(i+1);

			//----------------------------------------------------
			// Considering the triangle (ptLastInserted,ptCurrent,ptNext) :
			// Compute sin(theta)=A^B/(|A| |B|)
			//----------------------------------------------------
			vnl_vector<double> v1(3), v2(3), crossprod(3);
			for ( int i=0; i<3; i++ ) { 
				v1[i] = lastInsertedPoint[i]-currentPoint[i] ;
				v2[i] = nextPoint[i]-currentPoint[i] ;
			}
			crossprod = vnl_cross_3d<double>(v1, v2);
			double angle = asin ( crossprod.two_norm() / (v1.two_norm()*v2.two_norm()) );
			
			//----------------------------------------------------
			// If the angle is too small, remove this point!
			//----------------------------------------------------
			if ( angle<angleMinRad ) insertThisPoint=false;
		}

		//----------------------------------------------------
		// Insert the point if needed
		//----------------------------------------------------
		if ( insertThisPoint )
		{
			lastInsertedPoint = currentPoint;
			outTract->AddVertex( currentPoint );
		}
	}

	return outTract;
}



  /*
    Load the multi fiber model
  */
  template <const unsigned int TImageDimension, class TComponentType>
  int crl::MFMTractGenerator<TImageDimension, TComponentType>
  ::LoadMFM()
  {
  // Read in the tensor images.
  typename TensorImageReaderType::Pointer reader = TensorImageReaderType::New();
  
  std::cout << "- Load data..." << std::endl;

  m_TensorImageVector.resize(    m_InputFileNames.size() );
  m_LogTensorImageVector.resize( m_InputFileNames.size() );
  m_TensorInterpolatorVector.resize( m_InputFileNames.size() );

  for (unsigned int i = 0; i < m_InputFileNames.size(); i++) {
	std::cout << "   * "<<m_InputFileNames[i].c_str()<<std::endl;
    reader->SetFileName( m_InputFileNames[i] );

    try {
      reader->GenerateOutputInformation();
    } catch (itk::ExceptionObject &excp)
    {
      std::cerr << "Failed to read from file " << 
             m_InputFileNames[i] << std::endl; 
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

    try {
      reader->Update();
    } catch (itk::ExceptionObject &excp)
    {
      (void)excp;		//avoid compiler warning
      std::cerr << "Failed to read from file " << 
             m_InputFileNames[i] << std::endl; 
    return EXIT_FAILURE;
    }

    m_TensorImageVector[i] = reader->GetOutput();
    m_TensorImageVector[i]->DisconnectPipeline();
  }

  // Read in the Fraction image
  // The Fraction image is a 4D file of floats, where the 4th dimension selects the corresponding tensor.
  // The first (second, etc...) fraction always refers to the first (second, etc...) tensor.
  // In the case of the estimation of an isotropic part, then the isotropic part fraction
  // is the (NbTensor+1)th fraction.
  if ( m_FractionFileName!="" ) {
   	try {
		std::cout << "   * "<<m_FractionFileName.c_str()<<std::endl;
		m_FractionImage = MFM_OpenFractionsImageAsVectorImage<float>(m_FractionFileName);
	} 
	catch(itk::ExceptionObject &excp)
    {
		std::cerr << "Failed to read from file " << m_FractionFileName << std::endl; 
		std::cerr << "Error: "<<excp.GetDescription() << std::endl;
		return EXIT_FAILURE;
    }
  }

  return 0;
  }

#include "crlMFMTensorEstimate/clock.h"

  template <const unsigned int TImageDimension, class TComponentType>
  void crl::MFMTractGenerator<TImageDimension, TComponentType>
  ::InitializeAlgo(int nbThreads)
  {
  std::cout << "- Initialize the algorithm..." <<std::endl;

  // Construct a log tensor vector image in m_LogTensorImageVector
  typedef itk::TensorLogImageFilter<ImageType, ImageType> 
                                                   TensorLogImageFilterType;
  typename TensorLogImageFilterType::Pointer tensorLogFilter =
                                            TensorLogImageFilterType::New();

  if ( nbThreads>0 ) tensorLogFilter->SetNumberOfThreads(nbThreads);

  for (unsigned int i = 0; i < m_TensorImageVector.size(); i++) {
    tensorLogFilter->SetInput( m_TensorImageVector[i] );
    tensorLogFilter->Update();
    m_LogTensorImageVector[i] = tensorLogFilter->GetOutput();
  //  tensorLogFilter->Update();
    m_LogTensorImageVector[i]->DisconnectPipeline();
  }

  // connect it to tensor image linear interpolator.
  for (unsigned int i = 0; i < m_TensorImageVector.size(); i++) {
    m_TensorInterpolatorVector[i] = 
      crl::TensorLinearInterpolateImageFunction< ImageType, double >::New();
    m_TensorInterpolatorVector[i]->SetInputImage(m_LogTensorImageVector[i]);
  }
  }


  /*
    Load the seeding region
  */
  template <const unsigned int TImageDimension, class TComponentType>
  int crl::MFMTractGenerator<TImageDimension, TComponentType>
  ::LoadROI()
  {
  // Read in the seeding ROI.
  typename ROIImageReaderType::Pointer ROIReader = ROIImageReaderType::New();
  std::cout << "   * "<<m_ROIFileName.c_str()<<std::endl;
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

  // Quickly compute the number of seed points to show the tract generator progress 
  // (should be really quick, doesn't add a significant processing time)
  typedef itk::ImageRegionConstIterator<ROIImageType> ROIIteratorType;
  ROIIteratorType itROI(m_ROIImage, m_ROIImage->GetLargestPossibleRegion());
  m_NbVoxelsInROI = 0;
  if ( m_ShowProgress )
  {
	for ( itROI.GoToBegin(); !itROI.IsAtEnd() ; ++itROI )
	{
	  unsigned short label = itROI.Value();
	  for ( unsigned int labelID=0; labelID<m_ROILabels.size() ; labelID++ )
	    if ( label==m_ROILabels[labelID] )  m_NbVoxelsInROI++;
	}
  }
  return 0;
 }

template <const unsigned int TImageDimension, class TComponentType>
vtkPolyData* crl::MFMTractGenerator<TImageDimension, TComponentType>
::createPolyDataFromITK( typename TractContainerType::Pointer itkTracts, typename TractColorContainerType::Pointer itkTractColors) const
{
	typename ImageType::PointType	itkPoint;
	typename TractType::Pointer		itkTract = 0;

	//---------------------------------------------------------
	// Allocate the new VTK Polydata 
	//---------------------------------------------------------
	vtkPolyData*	newTracts = vtkPolyData::New();
	vtkCellArray*	newTractArray = vtkCellArray::New();
	vtkPoints*		newPoints = vtkPoints::New();
	double			vtkPoint[3];

	//---------------------------------------------------------
	// Precompute the total number of points to avoid a memory 
	// error when calling newPoints->InsertNextPoint
	//---------------------------------------------------------
	vtkIdType numPoints=0;
	for (unsigned int i = 0; i < itkTracts->Size() ; i++) 
		numPoints += itkTracts->GetElement(i)->GetVertexList()->Size();
	newPoints->Allocate(numPoints);

	//---------------------------------------------------------
	// For each tract
	//---------------------------------------------------------
	for (unsigned int i = 0; i < itkTracts->Size() ; i++) 
	{
		//---------------------------------------------------------
		// Take the tract and loop on all points
		//---------------------------------------------------------
		itkTract = itkTracts->GetElement(i);
		vtkIdList *vtkOneTract = vtkIdList::New();
		for (unsigned int j = 0; j < itkTract->GetVertexList()->Size(); j++) 
		{
			typename ImageType::PointType itkPoint = itkTract->GetVertexList()->GetElement(j);
			vtkPoint[0] = itkPoint[0];
			vtkPoint[1] = itkPoint[1];
			vtkPoint[2] = itkPoint[2];
			
			vtkIdType newID = newPoints->InsertNextPoint(vtkPoint);
			if (newID==-1) std::cout<<"MEMORY ERROR. Cannot add an additional tract point."<<std::endl;
			vtkOneTract->InsertNextId(newID);
		}
		newTractArray->InsertNextCell(vtkOneTract);
	}
	newPoints->Squeeze();
	newTractArray->Squeeze();

	newTracts->SetPoints( newPoints );
	newTracts->SetLines( newTractArray );

	//---------------------------------------------------------
	// Current color mode: one RGB color for the whole tract
	//---------------------------------------------------------
	if ( itkTractColors.GetPointer()!=NULL && itkTracts->Size()==itkTractColors->Size() ) 
	{
		vtkUnsignedCharArray *pointColor = vtkUnsignedCharArray::New();
		pointColor->SetNumberOfComponents(3);
		pointColor->SetNumberOfTuples(newPoints->GetNumberOfPoints());

		int id=0;
		for (unsigned int i = 0; i < itkTractColors->Size() ; i++) 
		{
			vnl_vector_fixed<float,3> rgb = itkTractColors->GetElement(i);
			for (unsigned int j = 0; j < itkTracts->GetElement(i)->GetVertexList()->Size(); j++)
			{
				pointColor->SetTuple3(id, (unsigned char)(unsigned int)(255*rgb[0]),(unsigned char)(unsigned int)(255*rgb[1]),(unsigned char)(unsigned int)(255*rgb[2]));
				id++;
			}
		}

		newTracts->GetPointData()->SetScalars(pointColor);
	}

	return newTracts;
}

// Write the tracts
// Depending on the file extension, write in vtp, gii, or vtk (default) file format
// If defined, itkTractColors can contain the color of each FULL tract (not one color per point)
template <const unsigned int TImageDimension, class TComponentType>
int crl::MFMTractGenerator<TImageDimension, TComponentType>
::WriteTracts( typename TractContainerType::Pointer tracts, typename TractColorContainerType::Pointer tractColors, const std::string& fileName) const
{
	//-------------------------------------------------
	// For vtp and gii extensions, writes using vtk
	//-------------------------------------------------
	crl::FileName fn(fileName);
	if ( fn.getExtension()=="vtp" || fn.getExtension()=="gii" )
	{
		// Creates the poly data
		vtkPolyData* polyData = createPolyDataFromITK( tracts, NULL); //tractColors );
		if ( !crlVtkMeshIO::WriteMesh(polyData, fileName, true)) return -1;
	}
	//-------------------------------------------------
	// Else, for vtk extension, use previous code
	// (could use VTK writer, but to be sure for compatibility)
	//-------------------------------------------------
	else
	{
		typename TractType::Pointer tract = 0;
		typename ImageType::PointType point;

		//-------------------------------------
		// Open the file for writing
		//-------------------------------------
		std::ofstream outFile;
		outFile.open( fileName.c_str(), std::ios::out);
		if (!outFile.is_open()) {
			std::cerr << "Failed to open file " << fileName << std::endl;
			return 1; // Failure
		}

		//-------------------------------------
		// Prepare the vtk header
		//-------------------------------------
		typename TractType::ContinuousIndexType cindex;
		outFile << "# vtk DataFile Version 3.0" << std::endl;
		outFile << "vtk output" << std::endl;
		outFile << "ASCII" << std::endl << "DATASET POLYDATA" << std::endl;

		//-------------------------------------
		// Output the number of points
		//-------------------------------------
		unsigned long int tractPointCount = 0;
		for (unsigned int i = 0; i < tracts->Size() ; i++)
			tractPointCount += tracts->GetElement(i)->GetVertexList()->Size(); 
		outFile << "POINTS " << tractPointCount << " float " << std::endl;

		//-------------------------------------
		// Output the point list
		//-------------------------------------
		for (unsigned int i = 0; i < tracts->Size() ; i++) {
			tract = tracts->GetElement(i);
			for (unsigned int j = 0; j < tract->GetVertexList()->Size(); j++) {
				point = tract->GetVertexList()->GetElement(j);
				//this->m_TensorImageVector[0]->TransformContinuousIndexToPhysicalPoint(cindex, point);
				for (unsigned int k = 0; k < ImageType::ImageDimension; k++) {
					outFile << point[k] << " " ;
				}
			}
			outFile << std::endl;
		}
		outFile << std::endl;

		//-------------------------------------
		// Now print out the poly line data of the tracts.
		//-------------------------------------
		outFile << "LINES " << tracts->Size() << " " << tractPointCount+tracts->Size() << std::endl;
		unsigned long pointCounter = 0;
		for (unsigned int i = 0; i < tracts->Size(); i++) {
			tract = tracts->GetElement(i);
			outFile << tract->GetVertexList()->Size() << " ";
			for (unsigned int j = 0; j < tract->GetVertexList()->Size(); j++) {
				outFile << pointCounter << " " ;
				pointCounter++;
			}
			outFile << std::endl;
		}

		if ( tractColors.GetPointer()!=NULL && tracts->Size()==tractColors->Size() ) 
		{
			outFile << std::endl;
			outFile << "POINT_DATA " << tractPointCount << std::endl;
			outFile << "COLOR_SCALARS scalars 3" << std::endl;

			for (unsigned int i = 0; i < tractColors->Size() ; i++) {
				vnl_vector_fixed<float,3> color = tractColors->GetElement(i);
				for (unsigned int j = 0; j < tracts->GetElement(i)->GetVertexList()->Size(); j++) {
					outFile << color[0] << " " << color[1] << " " << color[2] << std::endl;
				}
			}
		}
		outFile << std::endl;

		//-------------------------------------
		// Close the file
		//-------------------------------------
		outFile << std::endl;
		outFile.close();
	}

	return 0;
}

template <const unsigned int TImageDimension, class TComponentType>
int crl::MFMTractGenerator<TImageDimension, TComponentType>
::Execute()
{
  typedef typename PixelType::ValueType ValueType;
  typedef typename PixelType::EigenVectorsMatrixType MatrixType;
  typedef typename PixelType::EigenValuesArrayType EigenValuesType;


  
  // Load the data and initialize the algorithm
  int r;  
  if ( (r=LoadMFM())!=0 ) return r;
  if ( (r=LoadROI())!=0 ) return r;
  InitializeAlgo();

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
  TensorImageIteratorType it(m_TensorImageVector[0], 
                      m_TensorImageVector[0]->GetLargestPossibleRegion() );
  PixelType in;

  std::cout << "- Generate tracts..." << std::endl;

  crl::TaskProgress generatorProgress(0, 1, m_NbVoxelsInROI );
  generatorProgress.SetLinePrefix("  ");
  if ( m_ShowProgress )
	  generatorProgress.InitShowProgress();

  unsigned long int tractCount = 0;
  unsigned long int tractPointCount = 0;
  unsigned long int tractTopologyCount = 0;
  typename TractType::Pointer tract = 0;
  typename ImageType::PointType point;
  typename ImageType::IndexType index;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    in = it.Get();
    index = it.GetIndex();

	// For all seeding label values
	for ( unsigned int labelID=0; labelID<m_ROILabels.size() ; labelID++ ) {

	  // Check that the current voxel is labelled one of our seeding label
	  if (m_ROIImage->GetPixel( index ) == m_ROILabels[labelID] ) {

        if ( m_ShowProgress ) generatorProgress.Update();

        // Estimate the tract that touches this location : generate 
		// m_TractsPerVoxel tracts from this voxel
		for (unsigned int i = 0; i < m_TractsPerVoxel; i++) {
          
		  // Select a random subvoxel coordinate inside the current voxel
		  m_TensorImageVector[0]->TransformIndexToPhysicalPoint(index, point);
          for (unsigned int j = 0; j < ImageType::ImageDimension; j++) {
            point[j] += ( (m_RandGen->GetVariate()-0.5)*
                            m_TensorImageVector[0]->GetSpacing()[j] );
          }

		  // Estimate the tract from this subvoxel coordinate
          EstimateTract(point, tract);

		  // Record the generated tract in the tract list
          if (tract->GetVertexList()->Size() > 0) {
            tracts->InsertElement( tractCount, tract );
            ++tractCount;
            tractPointCount += tract->GetVertexList()->Size();
            tractTopologyCount += (tract->GetVertexList()->Size());
          } else {
            ; // std::cout << "Tract has length 0" << std::endl;
		  }
        } // for (unsigned int i = 0; i < m_TractsPerVoxel; i++)
      } // if (m_ROIImage->GetPixel( index ) == m_ROILabels.at(labelID) )
	}  // for ( int labelID=0; labelID<m_ROILabels.size() ; labelID++ )
  }
  tractTopologyCount += tractCount;
  // std::cout << "Number of tracts identified: " << tractCount << std::endl;

  /* Now we will write out the tracts. We will write an ASCII format VTK file.
   */
  std::cout << "- Save tracts in <"<<m_OutputFileName <<">."<< std::endl;
  TractColorContainerType::Pointer noColor;
  WriteTracts(tracts, noColor, this->m_OutputFileName);

 /* std::ofstream outFile;
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
	  // Code of EstimateTract modified to directly output physical coordinates
	  //     cindex = tract->GetVertexList()->GetElement(j);	
	  //     m_TensorImageVector[0]->TransformContinuousIndexToPhysicalPoint(cindex, point);
      point = tract->GetVertexList()->GetElement(j);
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
  outFile.close();*/

  return 0; // Success
}

#endif

