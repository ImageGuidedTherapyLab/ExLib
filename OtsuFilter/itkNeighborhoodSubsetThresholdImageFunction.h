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
#ifndef itkNeighborhoodSubsetThresholdImageFunction_h
#define itkNeighborhoodSubsetThresholdImageFunction_h

#include "itkBinaryThresholdImageFunction.h"

namespace itk
{
/**
 * \class NeighborhoodSubsetThresholdImageFunction
 * \brief Determine whether all the pixels in the specified neighborhood meet a threshold criteria
 *
 * Determine whether all the pixels in the specified neighborhood meet
 * a threshold criteria.
 *
 * If called with a ContinuousIndex or Point, the calculation is performed
 * at the nearest neighbor.
 *
 * This class is templated over the input image type and the coordinate
 * representation type (e.g. float or double).
 *
 * \ingroup ImageFunctions
 * \ingroup ITKImageFunction
 */
template <typename TInputImage, typename TCoordRep = float>
class ITK_TEMPLATE_EXPORT NeighborhoodSubsetThresholdImageFunction : public BinaryThresholdImageFunction<TInputImage, TCoordRep>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(NeighborhoodSubsetThresholdImageFunction);

  /** Standard class type aliases. */
  using Self = NeighborhoodSubsetThresholdImageFunction;
  using Superclass = NeighborhoodSubsetThresholdImageFunction<TInputImage, TCoordRep>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(NeighborhoodSubsetThresholdImageFunction, SubsetThresholdImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** InputImageType type alias support */
  using InputImageType = TInputImage;

  /** OutputType typdef support. */
  using OutputType = typename Superclass::OutputType;

  /** Index type alias support */
  using IndexType = typename Superclass::IndexType;

  /** ContinuousIndex type alias support */
  using ContinuousIndexType = typename Superclass::ContinuousIndexType;

  /** Point type alias support */
  using PointType = typename Superclass::PointType;

  /** PixelType type alias support */
  using PixelType = typename Superclass::PixelType;

  /** Dimension of the underlying image. */
  static constexpr unsigned int ImageDimension = InputImageType::ImageDimension;

  /** SizeType of the input image */
  using InputSizeType = typename InputImageType::SizeType;

  /** Set the radius of the neighborhood used in computation. */
  itkSetMacro(Radius, InputSizeType);

  /** Set the percentage of the neighborhood inside */
  itkSetMacro(PercentInside, float);

  /** Get the radius of the neighborhood used in computation */
  itkGetConstReferenceMacro(Radius, InputSizeType);

  /** Evalulate the function at specified index */
  bool
  EvaluateAtIndex(const IndexType & index) const override;

  /** Evaluate the function at non-integer positions */
  bool
  Evaluate(const PointType & point) const override
  {
    IndexType index;

    this->ConvertPointToNearestIndex(point, index);
    return this->EvaluateAtIndex(index);
  }

  bool
  EvaluateAtContinuousIndex(const ContinuousIndexType & cindex) const override
  {
    IndexType index;

    this->ConvertContinuousIndexToNearestIndex(cindex, index);
    return this->EvaluateAtIndex(index);
  }

protected:
  NeighborhoodSubsetThresholdImageFunction();
  ~NeighborhoodSubsetThresholdImageFunction() override = default;
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  InputSizeType m_Radius;
  float         m_PercentInside;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkNeighborhoodSubsetThresholdImageFunction.hxx"
#endif

/*
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNeighborhoodSubsetThresholdImageFunction.hxx"
#endif
*/

#endif
