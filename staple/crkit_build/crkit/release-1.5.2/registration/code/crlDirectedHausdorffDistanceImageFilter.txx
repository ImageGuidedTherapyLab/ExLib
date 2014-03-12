/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkDirectedHausdorffDistanceImageFilter.txx,v $
Language:  C++
Date:      $Date: 2008-10-14 19:20:33 $
Version:   $Revision: 1.14 $

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __crlDirectedHausdorffDistanceImageFilter_txx
#define __crlDirectedHausdorffDistanceImageFilter_txx
#include "crlDirectedHausdorffDistanceImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNumericTraits.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkProgressReporter.h"

#include <vector>
#include <algorithm>


template<class TInputImage1, class TInputImage2>
crlDirectedHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::crlDirectedHausdorffDistanceImageFilter(): m_MaxDistance(1)
{

  // this filter requires two input images
  this->SetNumberOfRequiredInputs( 2 );

  m_DistanceMap = NULL;
  m_DirectedHausdorffDistance = itk::NumericTraits<RealType>::Zero;
  m_AverageHausdorffDistance = itk::NumericTraits<RealType>::Zero;
}


template<class TInputImage1, class TInputImage2>
void
crlDirectedHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::SetInput2( const TInputImage2 * image )
{
  this->SetNthInput(1, const_cast<TInputImage2 *>( image ) );
}


template<class TInputImage1, class TInputImage2>
const typename crlDirectedHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::InputImage2Type *
crlDirectedHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::GetInput2()
{
  return static_cast< const TInputImage2 * >
    (this->itk::ProcessObject::GetInput(1));
}

template<class TInputImage1, class TInputImage2>
void
crlDirectedHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
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
crlDirectedHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::EnlargeOutputRequestedRegion(itk::DataObject *data)
{
  Superclass::EnlargeOutputRequestedRegion(data);
  data->SetRequestedRegionToLargestPossibleRegion();
}


template<class TInputImage1, class TInputImage2>
void
crlDirectedHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::AllocateOutputs()
{
  // Pass the first input through as the output
  InputImage1Pointer image =
    const_cast< TInputImage1 * >( this->GetInput1() );
  this->GraftOutput( image );
}


template<class TInputImage1, class TInputImage2>
void
crlDirectedHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::BeforeThreadedGenerateData()
{
  int numberOfThreads = this->GetNumberOfThreads();

  // Resize the thread temporaries
  m_MaxDistance.SetSize(numberOfThreads);
  m_PixelCount.SetSize(numberOfThreads);
  m_Sum.SetSize(numberOfThreads);
  m_Distances.resize(numberOfThreads);

  // Initialize the temporaries
  m_MaxDistance.Fill(itk::NumericTraits<RealType>::Zero);
  m_PixelCount.Fill(0);
  m_Sum.Fill(itk::NumericTraits<RealType>::Zero);

  // Compute Danielsson distance from non-zero pixels in the second image
  typedef itk::DanielssonDistanceMapImageFilter<InputImage2Type,DistanceMapType>
    FilterType;

  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput( this->GetInput2() );
  filter->UseImageSpacingOn();
  filter->Update();

  m_DistanceMap = filter->GetOutput();

}


template<class TInputImage1, class TInputImage2>
void
crlDirectedHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::AfterThreadedGenerateData()
{
  int i;

  int numberOfThreads = this->GetNumberOfThreads();

  m_DirectedHausdorffDistance = itk::NumericTraits<RealType>::Zero;
  RealType sum = itk::NumericTraits<RealType>::Zero;
  unsigned int pixelcount = 0;

  // find max over all threads
  for( i = 0; i < numberOfThreads; i++)
    {
      if ( m_MaxDistance[i] > m_DirectedHausdorffDistance )
	{
	  m_DirectedHausdorffDistance = m_MaxDistance[i];
	}
      pixelcount += m_PixelCount[i];
      sum += m_Sum[i];
      for (unsigned long j = 0; j < m_Distances.at(i).size(); j++)
	{
	  m_SortedDistances.push_back(m_Distances.at(i).at(j));
	}
    }

  m_AverageHausdorffDistance = sum / (RealType) pixelcount;

  sort(m_SortedDistances.begin(),m_SortedDistances.end());

  // clean up
  m_DistanceMap = NULL;

}

template<class TInputImage1, class TInputImage2>
void
crlDirectedHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::ThreadedGenerateData(const RegionType& regionForThread,
                       int threadId)
{

  itk::ImageRegionConstIterator<TInputImage1> it1 (this->GetInput1(), regionForThread);
  itk::ImageRegionConstIterator<DistanceMapType> it2 (m_DistanceMap, regionForThread);

  // support progress methods/callbacks
  itk::ProgressReporter progress(this, threadId, regionForThread.GetNumberOfPixels());

  // do the work
  while (!it1.IsAtEnd())
    {

      if( it1.Get() != itk::NumericTraits<InputImage1PixelType>::Zero )
	{
	  if ( it2.Get() > m_MaxDistance[threadId] )
	    {
	      m_MaxDistance[threadId] = it2.Get();
	    }
	  m_PixelCount[threadId]++;
	  m_Sum[threadId] += it2.Get();
	  m_Distances.at(threadId).push_back(it2.Get());
	}

      ++it1;
      ++it2;

      progress.CompletedPixel();
    }
}


template<class TInputImage1, class TInputImage2>
void
crlDirectedHausdorffDistanceImageFilter<TInputImage1, TInputImage2>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "DirectedHausdorffDistance: "
     << m_DirectedHausdorffDistance << std::endl;
  os << indent << "AverageHausdorffDistance: "
     << m_AverageHausdorffDistance << std::endl;
}


template<class TInputImage1, class TInputImage2>
double crlDirectedHausdorffDistanceImageFilter<TInputImage1, TInputImage2>::GetDirectedHausdorffDistance(double percentage)
{
  if ((percentage < 0.0) || (percentage > 1.0))
    return 0;

  unsigned long element = (unsigned long)ceil((m_SortedDistances.size()-1)*percentage);
  if (element < 0)
    element = 0;
  if (element > m_SortedDistances.size()-1)
    element = m_SortedDistances.size()-1;

  return m_SortedDistances[element];
}


#endif
