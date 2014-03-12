/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkHausdorffDistanceImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-04-25 12:27:26 $
  Version:   $Revision: 1.12 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __crlHausdorffDistanceImageFilter_h
#define __crlHausdorffDistanceImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkNumericTraits.h"

/** \class HausdorffDistanceImageFilter
 * \brief Computes the Hausdorff distance between the set of
 * non-zero pixels of two images.
 *
 *
 * HausdorffDistanceImageFilter computes the distance between the set
 * non-zero pixels of two images using the following formula:
 * \f[ H(A,B) = \max(h(A,B),h(B,A)) \f]
 * where
 * \f[ h(A,B) = \max_{a \in A} \min_{b \in B} \| a - b\| \f] is the directed
 * Hausdorff distance
 * and \f$A\f$ and \f$B\f$ are respectively the set of non-zero pixels
 * in the first and second input images.
 *
 * In particular, this filter uses the DirectedHausdorffImageFilter inside to
 * compute the two directed distances and then select the largest of the two.
 *
 * The Hausdorff distance measures the degree of mismatch between two sets and
 * behaves like a metric over the set of all closedm bounded sets -
 * with properties of identity, symmetry and triangle inequality.
 *
 * This filter requires the largest possible region of the first image
 * and the same corresponding region in the second image.
 * It behaves as filter with
 * two input and one output. Thus it can be inserted in a pipeline with
 * other filters. The filter passes the first input through unmodified.
 *
 * This filter is templated over the two input image type. It assume
 * both image have the same number of dimensions.
 *
 * \sa DirectedHausdorffDistanceImageFilter
 *
 * \ingroup MultiThreaded
 */
template<class TInputImage1, class TInputImage2>
class crlHausdorffDistanceImageFilter :
  public itk::ImageToImageFilter<TInputImage1, TInputImage1>
{
public:
  /** Standard Self typedef */
  typedef crlHausdorffDistanceImageFilter                   Self;
  typedef itk::ImageToImageFilter<TInputImage1,TInputImage1>  Superclass;
  typedef itk::SmartPointer<Self>                             Pointer;
  typedef itk::SmartPointer<const Self>                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(crlHausdorffDistanceImageFilter, itk::ImageToImageFilter);

  /** Image related typedefs. */
  typedef TInputImage1                        InputImage1Type;
  typedef TInputImage2                        InputImage2Type;
  typedef typename TInputImage1::Pointer      InputImage1Pointer;
  typedef typename TInputImage2::Pointer      InputImage2Pointer;
  typedef typename TInputImage1::ConstPointer InputImage1ConstPointer;
  typedef typename TInputImage2::ConstPointer InputImage2ConstPointer;

  typedef typename TInputImage1::RegionType RegionType;
  typedef typename TInputImage1::SizeType   SizeType;
  typedef typename TInputImage1::IndexType  IndexType;

  typedef typename TInputImage1::PixelType InputImage1PixelType;
  typedef typename TInputImage2::PixelType InputImage2PixelType;

  /** Image related typedefs. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage1::ImageDimension);

  /** Type to use form computations. */
  typedef typename itk::NumericTraits<InputImage1PixelType>::RealType RealType;

  /** Set the first input. */
  void SetInput1( const InputImage1Type * image )
    { this->SetInput( image ); }

  /** Set the second input. */
  void SetInput2( const InputImage2Type * image );

  /** Get the first input. */
  const InputImage1Type * GetInput1(void)
    { return this->GetInput(); }

  /** Get the second input. */
  const InputImage2Type * GetInput2(void);

  /** Return the computed Hausdorff distance. */
  itkGetConstMacro(HausdorffDistance,RealType);
  itkGetConstMacro(AverageHausdorffDistance,RealType);

  double GetHausdorffDistance(int percentage);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(Input1HasNumericTraitsCheck,
		  (itk::Concept::HasNumericTraits<InputImage1PixelType>));
  /** End concept checking */
#endif

protected:
  crlHausdorffDistanceImageFilter();
  ~crlHausdorffDistanceImageFilter(){};
  void PrintSelf(std::ostream& os, itk::Indent indent) const;

  /** GenerateData. */
  void  GenerateData();

  // Override since the filter needs all the data for the algorithm
  void GenerateInputRequestedRegion();

  // Override since the filter produces all of its output
  void EnlargeOutputRequestedRegion(itk::DataObject *data);

private:
  crlHausdorffDistanceImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  RealType                            m_HausdorffDistance;
  RealType                            m_AverageHausdorffDistance;
  RealType m_HausdorffDistances[101];

}; // end of class


#ifndef ITK_MANUAL_INSTANTIATION
#include "crlHausdorffDistanceImageFilter.txx"
#endif

#endif
