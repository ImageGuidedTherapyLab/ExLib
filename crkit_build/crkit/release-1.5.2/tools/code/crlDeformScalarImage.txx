/*
 * Copyright 2007 Simon Warfield simon.warfield@childrens.harvard.edu
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#ifndef _CRL_DEFORM_SCALAR_IMAGE_TXX_INCLUDED
#define _CRL_DEFORM_SCALAR_IMAGE_TXX_INCLUDED 1

#define USE_BROKEN_WARPIMAGE_FILTER 0

template <const unsigned int TImageDimension, class TPixelType >
int crlDeformScalarImage< TImageDimension, TPixelType >
::Execute()
{

  typename ReaderType::Pointer reader = ReaderType::New();
  typename WriterType::Pointer writer = WriterType::New();
  reader->SetFileName( m_InputFileName );
  writer->SetFileName( m_OutputFileName );

  try {
    reader->Update();
  } catch (itk::ExceptionObject &ex) {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }

  typename FieldReaderType::Pointer fieldReader = FieldReaderType::New();
  fieldReader->SetFileName( m_DeformationFieldFileName );
  try {
    fieldReader->Update();
  }  catch (itk::ExceptionObject &ex) {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }

  typename DeformationFieldType::ConstPointer deformationField = fieldReader->GetOutput();
  typedef itk::InterpolateImageFunction<ImageType, double>   InterpolatorType;
  typename InterpolatorType::Pointer interp;
  


  // set up the proper interpolator
  if ( !strncmp( m_InterpolationType.c_str(), "linear",6 ) )
    {  // no arg specified or set to linear
    typename LInterpolatorType::Pointer interpolator = LInterpolatorType::New();
    interp = interpolator;
    }
  else if ( !strncmp( m_InterpolationType.c_str(), "nearest",7 ) )
    {  // nearest neighbor interpolation - good for segmented data
    typename NInterpolatorType::Pointer interpolator = NInterpolatorType::New();
    interp = interpolator;
    }
  else if ( !strncmp( m_InterpolationType.c_str(), "bspline",6 ) )
    {  // spline interpolation - watch for ringing.
    typename BSplineInterpolatorType::Pointer interpolator = 
                                 BSplineInterpolatorType::New();
    interp = interpolator;
    }
  else if ( !strncmp( m_InterpolationType.c_str(), "sinc",4 ) )
    {  // windowed sinc interpolation.
    typename SincInterpolatorType::Pointer interpolator = 
                                 SincInterpolatorType::New();
    interp = interpolator;
    }
  else
    {  // unknown interpolation type
    std::cerr << "unknown interpolation type: " << m_InterpolationType
              << std::endl;
    return EXIT_FAILURE;
    }

  /* Since we may be working with oriented images, ensure that the output
   * direction is set properly.
   */
#if USE_BROKEN_WARPIMAGE_FILTER
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInterpolator(interp);
  filter->SetOutputSpacing( deformationField->GetSpacing() );
  filter->SetOutputOrigin( deformationField->GetOrigin() );
  filter->SetDeformationField( deformationField );

  filter->SetInput( reader->GetOutput() );
  writer->SetInput( filter->GetOutput() );
  

  // Ensure the necessary output image fields are set.
  filter->GenerateOutputInformation();
  filter->GetOutput()->SetDirection( deformationField->GetDirection() );
  std::cerr << "output directions BEFORE UPDATE: " << deformationField->GetDirection()
  << std::endl;
  filter->Update();
#warning "using possibly broken warpimage filter."
#else

  typename ImageType::Pointer output = ImageType::New();
  output->CopyInformation(deformationField);
  output->SetRegions( deformationField->GetLargestPossibleRegion() );
  output->Allocate();
  
  itk::ImageRegionConstIterator<DeformationFieldType> defIter(deformationField, deformationField->GetLargestPossibleRegion());
  itk::ImageRegionIteratorWithIndex<ImageType> outIter(output, output->GetLargestPossibleRegion());
  typename ImageType::PointType point;
  

  interp->SetInputImage(reader->GetOutput());
  
  for ( outIter.GoToBegin(), defIter.GoToBegin(); !outIter.IsAtEnd(); ++outIter, ++defIter )
    {
    output->TransformIndexToPhysicalPoint( outIter.GetIndex(), point );

    point += defIter.Get();
    if( interp->IsInsideBuffer( point ) )
        {
        TPixelType value = static_cast<TPixelType>( 
          interp->Evaluate( point ) );
        outIter.Set( value );
        }
      else
        {
        outIter.Set(0);
        }   
    }
  
  
  writer->SetInput(output);
#endif

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown " << std::endl;
    std::cerr << excp << std::endl;
    }

  return EXIT_SUCCESS;

};

#endif
