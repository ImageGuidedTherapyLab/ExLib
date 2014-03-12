
#ifndef _itkRotateTensorImageFilter_txx
#define _itkRotateTensorImageFilter_txx

#include "itkRotateTensorImageFilter.h"

namespace itk
{

/**
 *
 */
template <class TInputImage, class TOutputImage>
RotateTensorImageFilter<TInputImage, TOutputImage>
::RotateTensorImageFilter()
{
  //Identity Matrix
  m_RotationMatrix.SetIdentity();
}


/**
 *
 */
template <class TInputImage, class TOutputImage>
void
RotateTensorImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Internal Rotation Matrix : " << std::endl;
  for (unsigned int i =0; i<3;i++)
    {
    for (unsigned int j=0; j<3; j++)
      {
      os << m_RotationMatrix(i,j) << "  ";
      }
    os << std::endl;
    }
}



template <class TInputImage, class TOutputImage>
void
RotateTensorImageFilter<TInputImage, TOutputImage>
::PrintParameters()
{
for (unsigned int i =0; i<3;i++)
  {
  for (unsigned int j=0; j<3; j++)
    {
    std::cout << m_RotationMatrix(i,j) << "  ";
    }
    std::cout << std::endl;
  }
}


/**
 *
 */
template <class TInputImage, class TOutputImage>
void
RotateTensorImageFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  // set up the functor rotation matrix
  this->GetFunctor().SetRotationMatrix( this->m_RotationMatrix );
}


} // end namespace itk

#endif
