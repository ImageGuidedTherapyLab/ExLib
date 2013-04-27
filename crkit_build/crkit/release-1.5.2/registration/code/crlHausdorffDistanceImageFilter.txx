/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkHausdorffDistanceImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2008-10-16 16:45:09 $
  Version:   $Revision: 1.13 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __crlHausdorffDistanceImageFilter_txx
#define __crlHausdorffDistanceImageFilter_txx
#include "crlHausdorffDistanceImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNumericTraits.h"
#include "itkProgressAccumulator.h"
#include "crlDirectedHausdorffDistanceImageFilter.h"



template<class TInputImage1, class TInputImage2>
crlHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::crlHausdorffDistanceImageFilter()
{

  // this filter requires two input images
  this->SetNumberOfRequiredInputs( 2 );

  m_HausdorffDistance = itk::NumericTraits<RealType>::Zero;
  m_AverageHausdorffDistance = itk::NumericTraits<RealType>::Zero;
}


template<class TInputImage1, class TInputImage2>
void
crlHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::SetInput2( const TInputImage2 * image )
{
  this->SetNthInput(1, const_cast<TInputImage2 *>( image ) );
}


template<class TInputImage1, class TInputImage2>
const typename crlHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::InputImage2Type *
crlHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::GetInput2()
{
  return static_cast< const TInputImage2 * >
    (this->itk::ProcessObject::GetInput(1));
}

template<class TInputImage1, class TInputImage2>
void
crlHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  // this filter requires:
  // - the largeset possible region of the first image
  // - the corresponding region of the second image
  if ( this->GetInput1() )
    {
    InputImage1Pointer image1 =
      const_cast< InputImage1Type * >( this->GetInput1() );
    image1->SetRequestedRegionToLargestPossibleRegion();

    if ( this->GetInput2() )
      {
      InputImage2Pointer image2 =
        const_cast< InputImage2Type * >( this->GetInput2() );
      image2->SetRequestedRegion(
        this->GetInput1()->GetRequestedRegion() );
      }

    }
}


template<class TInputImage1, class TInputImage2>
void
crlHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::EnlargeOutputRequestedRegion(itk::DataObject *data)
{
  Superclass::EnlargeOutputRequestedRegion(data);
  data->SetRequestedRegionToLargestPossibleRegion();
}


template<class TInputImage1, class TInputImage2>
void
crlHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::GenerateData()
{

  // Pass the first input through as the output
  InputImage1Pointer image =
    const_cast< TInputImage1 * >( this->GetInput1() );
  this->GraftOutput( image );

  RealType distance12, distance21;

  // Create a process accumulator for tracking the progress of this minipipeline
  itk::ProgressAccumulator::Pointer progress = itk::ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  typedef crlDirectedHausdorffDistanceImageFilter<InputImage1Type,InputImage2Type>
    Filter12Type;

  typename Filter12Type::Pointer filter12 = Filter12Type::New();

  filter12->SetInput1( this->GetInput1() );
  filter12->SetInput2( this->GetInput2() );

  typedef crlDirectedHausdorffDistanceImageFilter<InputImage2Type,InputImage1Type>
    Filter21Type;

  typename Filter21Type::Pointer filter21 = Filter21Type::New();

  filter21->SetInput1( this->GetInput2() );
  filter21->SetInput2( this->GetInput1() );

  // Register the filter with the with progress accumulator using
  // equal weight proportion
  progress->RegisterInternalFilter(filter12,.5f);
  progress->RegisterInternalFilter(filter21,.5f);

  filter12->Update();
  distance12 = filter12->GetDirectedHausdorffDistance();
  filter21->Update();
  distance21 = filter21->GetDirectedHausdorffDistance();

  if ( distance12 > distance21 )
    {
    m_HausdorffDistance = distance12;
    }
  else
    {
    m_HausdorffDistance = distance21;
    }
  m_AverageHausdorffDistance = ( filter12->GetAverageHausdorffDistance() + filter21->GetAverageHausdorffDistance() ) / 2.0;

  for (int i = 0; i < 101; i++)
    {
      distance12 = filter12->GetDirectedHausdorffDistance((double)(i)/100);
      distance21 = filter21->GetDirectedHausdorffDistance((double)(i)/100);
      if (distance12 > distance21)
	m_HausdorffDistances[i] = distance12;
      else
	m_HausdorffDistances[i] = distance21;
    }
}

template<class TInputImage1, class TInputImage2>
void
crlHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "HausdorffDistance: "
     << m_HausdorffDistance << std::endl;
  os << indent << "AverageHausdorffDistance: "
     << m_AverageHausdorffDistance << std::endl;
}


template<class TInputImage1, class TInputImage2>
double crlHausdorffDistanceImageFilter<TInputImage1, TInputImage2>::GetHausdorffDistance(int percentage)
{
  return m_HausdorffDistances[percentage];
}

#endif
