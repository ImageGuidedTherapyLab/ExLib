/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkNeighborhoodSubsetThresholdImageFunction_hxx
#define itkNeighborhoodSubsetThresholdImageFunction_hxx

#include "itkNeighborhoodSubsetThresholdImageFunction.h"
#include "itkNumericTraits.h"
#include "itkConstNeighborhoodIterator.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TInputImage, typename TCoordRep>
NeighborhoodSubsetThresholdImageFunction<TInputImage,
TCoordRep>::NeighborhoodSubsetThresholdImageFunction()
{
  m_Radius.Fill(1);
  m_PercentInside=0.5;
}

/**
 *
 */
template <typename TInputImage, typename TCoordRep>
void
NeighborhoodSubsetThresholdImageFunction<TInputImage, TCoordRep>::PrintSelf(std::ostream & os, Indent indent) const
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "Radius: " << m_Radius << std::endl;
  os << indent << "Percent Inside: " << m_PercentInside << std::endl;
}

/**
 *
 */
template <typename TInputImage, typename TCoordRep>
bool
NeighborhoodSubsetThresholdImageFunction<TInputImage, TCoordRep>::EvaluateAtIndex(const IndexType & index) const
{
  if (!this->GetInputImage())
  {
    return (false);
  }

  if (!this->IsInsideBuffer(index))
  {
    return (false);
  }

  // Create an N-d neighborhood kernel, using a zeroflux boundary condition
  ConstNeighborhoodIterator<InputImageType> it(
    m_Radius, this->GetInputImage(), this->GetInputImage()->GetBufferedRegion());

  // Set the iterator at the desired location
  it.SetLocation(index);

  // Walk the neighborhood
  bool               mostInside = true;
  int                countOutside = 0;
  PixelType          lower = this->GetLower();
  PixelType          upper = this->GetUpper();
  PixelType          value;
  const unsigned int size = it.Size();
  for (unsigned int i = 0; i < size; ++i)
  {
    value = it.GetPixel(i);
    if (lower > value || value > upper)
    {
      countOutside++;
    }
  }
  if ( static_cast< float > (countOutside)/static_cast< float >(size) > m_PercentInside )
    {
      mostInside = false;
    }
  return (mostInside);
}
} // end namespace itk

#endif
