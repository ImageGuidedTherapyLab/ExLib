
#ifndef _TensorScalarMultiplyFunctorIncluded
#define _TensorScalarMultiplyFunctorIncluded 1

#include "TensorScalarMultiplyFunctor.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
void
TensorScalarMultiplyImageFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  // set up the functor scalar value
  this->GetFunctor().SetScalar( this->scalar );
}

} // end namespace itk

#endif
