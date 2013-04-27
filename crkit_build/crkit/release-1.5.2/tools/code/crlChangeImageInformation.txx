/*
 * Copyright 2008 Simon Warfield simon.warfield@childrens.harvard.edu
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 * http://crl.med.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#ifndef _CRL_CHANGEIMAGEINFORMATION_TXX_INCLUDED
#define _CRL_CHANGEIMAGEINFORMATION_TXX_INCLUDED

template < const unsigned int TDimension, class TPixelType >
int ChangeImageInformation< TDimension, TPixelType >
::Execute()
{

  unsigned int exitStatus = EXIT_SUCCESS;
  typename ImageReaderType::Pointer reader;
  typename ImageReaderType::Pointer refreader;
  typename ImageWriterType::Pointer writer;
  typename FilterType::Pointer filter;

  reader = ImageReaderType::New();
  reader->SetFileName(m_InputImageFileName.c_str());

  try {
    reader->Update();
  } catch (itk::ExceptionObject &e) {
    std::cerr << "Failed to read image named " << m_InputImageFileName << 
      std::endl;
        std::cerr << "Caught ITK exception: " << e << std::endl;
        exitStatus = EXIT_FAILURE;
  }


  filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  typename ImageTypeIn::RegionType inputRegion = 
                   reader->GetOutput()->GetLargestPossibleRegion();
  typename ImageTypeIn::SizeType size = inputRegion.GetSize();
  typename ImageTypeIn::IndexType start = inputRegion.GetIndex();

  if (m_UseReferenceImage) {
    refreader = ImageReaderType::New();
    refreader->SetFileName(m_ReferenceImageFileName.c_str());
    try {
      refreader->Update();
    } catch (itk::ExceptionObject &e) {
        std::cerr << "Failed to read reference image called " << 
        m_ReferenceImageFileName << std::endl;
        std::cerr << "Caught ITK exception: " << e << std::endl;
        exitStatus = EXIT_FAILURE;
    }
    filter->SetUseReferenceImage(true);
    filter->SetReferenceImage( refreader->GetOutput() );
  }
  filter->SetChangeSpacing( m_ChangeSpacing );
  filter->SetChangeDirection( m_ChangeDirection );
  filter->SetChangeOrigin( m_ChangeOrigin );
  filter->SetChangeRegion( m_ChangeRegion );

  writer = ImageWriterType::New();
  writer->SetFileName(m_OutputImageFileName.c_str());
  writer->SetInput( filter->GetOutput() );

  try {
      writer->Update();
      writer->UseCompressionOn();
      writer->Write();
  } catch (itk::ExceptionObject &e) {
        std::cerr << "Failed to write file called " 
                  << m_OutputImageFileName
                  << std::endl;
        std::cerr << "Caught ITK exception: " << e << std::endl;
        exitStatus = EXIT_FAILURE;
  }

  return exitStatus;
}

#endif 

