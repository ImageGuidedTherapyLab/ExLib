
#ifndef _itkDisplacementFieldJacobianRecursiveGaussianImageFilter_txx
#define _itkDisplacementFieldJacobianRecursiveGaussianImageFilter_txx

#include "itkDisplacementFieldJacobianRecursiveGaussianImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_svd.h"

namespace itk
{


/**
 * Constructor
 */
template <typename TInputImage, typename TOutputImage >
DisplacementFieldJacobianRecursiveGaussianImageFilter<TInputImage,TOutputImage>
::DisplacementFieldJacobianRecursiveGaussianImageFilter()
{

  m_NormalizeAcrossScale = false;
  m_ComputeRotationMatrix = true;

  // Create a process accumulator for tracking the progress of this
  // minipipeline
  m_Progress = ProgressAccumulator::New();
  m_Progress->SetMiniPipelineFilter(this);
 
  // Compute the contribution of each filter to the total progress.
  const double weight = 1.0 / ( ImageDimension );


  m_ImageAdaptorU = VectorScalarSelectionFilterType::New();
  m_ImageAdaptorV = VectorScalarSelectionFilterType::New();
  m_ImageAdaptorW = VectorScalarSelectionFilterType::New();

  m_GradientFilterU = GradientFilterType::New();
  m_GradientFilterV = GradientFilterType::New();
  m_GradientFilterW = GradientFilterType::New();

  m_GradientFilterU->SetNormalizeAcrossScale( m_NormalizeAcrossScale );
  m_GradientFilterV->SetNormalizeAcrossScale( m_NormalizeAcrossScale );
  m_GradientFilterW->SetNormalizeAcrossScale( m_NormalizeAcrossScale );
  

  m_Progress->RegisterInternalFilter( m_GradientFilterU, weight );
  m_Progress->RegisterInternalFilter( m_GradientFilterV, weight );
  m_Progress->RegisterInternalFilter( m_GradientFilterW, weight );

  this->SetSigma( 1.0 );

}



/**
 * Set value of Sigma
 */
template <typename TInputImage, typename TOutputImage>
void 
DisplacementFieldJacobianRecursiveGaussianImageFilter<TInputImage,TOutputImage>
::SetSigma( ComponentType sigma )
{

  m_GradientFilterU->SetSigma( sigma );
  m_GradientFilterV->SetSigma( sigma );
  m_GradientFilterW->SetSigma( sigma );

  this->Modified();

}



/**
 * Set Normalize Across Scale Space
 */
template <typename TInputImage, typename TOutputImage>
void 
DisplacementFieldJacobianRecursiveGaussianImageFilter<TInputImage,TOutputImage>
::SetNormalizeAcrossScale( bool normalize )
{

  if (m_NormalizeAcrossScale != normalize)  {
    m_NormalizeAcrossScale = normalize;

    m_GradientFilterU->SetNormalizeAcrossScale( normalize );
    m_GradientFilterV->SetNormalizeAcrossScale( normalize );
    m_GradientFilterW->SetNormalizeAcrossScale( normalize );

    this->Modified();
  }

}


/**
 * Set Compute Rotation Matrix. Boolean to control if the pure rotation matrix should be computed from the Jacobian.
 */
template <typename TInputImage, typename TOutputImage>
void 
DisplacementFieldJacobianRecursiveGaussianImageFilter<TInputImage,TOutputImage>
::SetComputeRotationMatrix( bool compute )
{

  if (m_ComputeRotationMatrix != compute) {
    m_ComputeRotationMatrix = compute;
    this->Modified();
    }

}


//
//
//
template <typename TInputImage, typename TOutputImage>
void
DisplacementFieldJacobianRecursiveGaussianImageFilter<TInputImage,TOutputImage>
::GenerateInputRequestedRegion() throw(InvalidRequestedRegionError)
{
  /// This function is the same as HessianRecursiveGaussianImageFilter, do we need this?

  // call the superclass' implementation of this method. this should
  // copy the output requested region to the input requested region
  Superclass::GenerateInputRequestedRegion();

  // This filter needs all of the input
  typename DisplacementFieldJacobianRecursiveGaussianImageFilter<TInputImage,TOutputImage>::InputImagePointer image = const_cast<InputImageType *>( this->GetInput() );
  image->SetRequestedRegion( this->GetInput()->GetLargestPossibleRegion() );
}


//
//
//
template <typename TInputImage, typename TOutputImage>
void
DisplacementFieldJacobianRecursiveGaussianImageFilter<TInputImage,TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *output)
{
  TOutputImage *out = dynamic_cast<TOutputImage*>(output);

  if (out)
    {
    out->SetRequestedRegion( out->GetLargestPossibleRegion() );
    }
}

/**
 * Compute filter for Gaussian kernel
 */
template <typename TInputImage, typename TOutputImage >
void
DisplacementFieldJacobianRecursiveGaussianImageFilter<TInputImage,TOutputImage >
::GenerateData(void)
{

  itkDebugMacro(<< "DisplacementFieldJacobianRecursiveGaussianImageFilter generating data ");
//std::cout << "entering GenerateData " << std::endl;

  m_Progress->ResetProgress();

  const typename TInputImage::ConstPointer   inputImage( this->GetInput() );

  OutputImagePointer outputimage = this->GetOutput();
  
  outputimage->SetLargestPossibleRegion( 
    inputImage->GetLargestPossibleRegion() );

  outputimage->SetBufferedRegion( 
    inputImage->GetBufferedRegion() );

  outputimage->SetRequestedRegion( 
    inputImage->GetRequestedRegion() );
  
  outputimage->Allocate();


  m_ImageAdaptorU->SetInput(  inputImage );
  m_ImageAdaptorU->SetIndex( 0 );
  m_ImageAdaptorV->SetInput(  inputImage );
  m_ImageAdaptorV->SetIndex( 1 );
  m_ImageAdaptorW->SetInput(  inputImage );
  m_ImageAdaptorW->SetIndex( 2 ); 
  
  m_GradientFilterU->SetInput( m_ImageAdaptorU->GetOutput() );
  m_GradientFilterV->SetInput( m_ImageAdaptorV->GetOutput() );
  m_GradientFilterW->SetInput( m_ImageAdaptorW->GetOutput() );

  m_GradientFilterU->Update();
  m_GradientFilterV->Update();
  m_GradientFilterW->Update();

  /// What does this do??
//  m_Progress->ResetFilterProgressAndKeepAccumulatedProgress();

   // Copy the results to the corresponding component
   // on the output image of vectors

   ImageRegionIteratorWithIndex< InputImageType > itU( 
              m_GradientFilterU->GetOutput(), 
              m_GradientFilterU->GetOutput()->GetRequestedRegion() );

   ImageRegionIteratorWithIndex< InputImageType > itV( 
              m_GradientFilterV->GetOutput(), 
              m_GradientFilterV->GetOutput()->GetRequestedRegion() );

   ImageRegionIteratorWithIndex< InputImageType > itW( 
              m_GradientFilterW->GetOutput(), 
              m_GradientFilterW->GetOutput()->GetRequestedRegion() );

   ImageRegionIteratorWithIndex< OutputImageType > ot( 
        outputimage, 
        outputimage->GetRequestedRegion() );
    

    itU.GoToBegin();
    itV.GoToBegin();
    itW.GoToBegin();
    ot.GoToBegin();
    typedef vnl_matrix<ComponentType> vnlMatrixType;
    typedef vnl_svd<ComponentType>  SVDSolverType;
    
    while( !itU.IsAtEnd() )  {
      vnlMatrixType tmp;
      tmp.set_size( 3, 3 );
      tmp(0,0) = 1 + itU.Get()[0];
      tmp(0,1) = itU.Get()[1];
      tmp(0,2) = itU.Get()[2];
      
      tmp(1,0) = itV.Get()[0];
      tmp(1,1) = 1 + itV.Get()[1];
      tmp(1,2) = itV.Get()[2];
      
      tmp(2,0) = itW.Get()[0];
      tmp(2,1) = itW.Get()[1];
      tmp(2,2) = 1 + itW.Get()[2];
      
//std::cout << "1\n" << tmp << std::endl;
      if (m_ComputeRotationMatrix)   {
        SVDSolverType svd( tmp, 1e-6 );

        vnlMatrixType m_VMatrix = svd.V();
        vnlMatrixType m_UMatrix = svd.U();
        tmp = m_UMatrix*m_VMatrix.transpose();               
      }
//std::cout << tmp << std::endl;      
      OutputPixelType outputpixel;
      outputpixel.SetIdentity();
      outputpixel *= tmp;

      ot.Set( outputpixel );
      ++itU;
      ++itV;
      ++itW;
      ++ot;
      }

}


template <typename TInputImage, typename TOutputImage>
void
DisplacementFieldJacobianRecursiveGaussianImageFilter<TInputImage,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << "NormalizeAcrossScale: " << m_NormalizeAcrossScale << std::endl;
}


} // end namespace itk

#endif
