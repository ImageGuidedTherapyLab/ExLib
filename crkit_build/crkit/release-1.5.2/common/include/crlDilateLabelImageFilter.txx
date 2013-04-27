
#ifndef __crlDilateLabelImageFilter_txx
#define __crlDilateLabelImageFilter_txx

#include "crlDilateLabelImageFilter.h"


template<class TInputImage, class TOutputImage, class TKernel>
crlDilateLabelImageFilter<TInputImage, TOutputImage, TKernel>
::crlDilateLabelImageFilter()
{
  //m_DilateBoundaryCondition.SetConstant( NumericTraits<PixelType>::NonpositiveMin() );
  //this->OverrideBoundaryCondition( &m_DilateBoundaryCondition );
}

template<class TInputImage, class TOutputImage, class TKernel>
typename crlDilateLabelImageFilter<TInputImage, TOutputImage, TKernel>::PixelType
crlDilateLabelImageFilter<TInputImage, TOutputImage, TKernel>
::Evaluate(const NeighborhoodIteratorType &nit,
           const KernelIteratorType kernelBegin,
           const KernelIteratorType kernelEnd)
{
  unsigned int i;
  PixelType max = NumericTraits<PixelType>::NonpositiveMin();
  PixelType temp;

  KernelIteratorType kernel_it;
	
  PixelType nCenPxVal=nit.GetCenterPixel();
	if ( nCenPxVal!=0 ) return nCenPxVal;

  for( i=0, kernel_it=kernelBegin; kernel_it<kernelEnd; ++kernel_it, ++i )
    {
    // if structuring element is positive, use the pixel under that element
    // in the image
    if( *kernel_it > NumericTraits<KernelPixelType>::Zero )
      {
      // note we use GetPixel() on the SmartNeighborhoodIterator to
      // respect boundary conditions
      temp = nit.GetPixel(i);

     /* if( temp > max )
        {
        max = temp;
        }*/
		if( temp > max )
        {
        max = temp;
        }
      }
    }
  
  return max;
} 


#endif

