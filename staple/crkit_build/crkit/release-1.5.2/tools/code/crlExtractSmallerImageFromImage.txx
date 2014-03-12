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

#ifndef _CRL_EXTRACTSMALLERIMAGERFROMIMAGE_TXX_INCLUDED
#define _CRL_EXTRACTSMALLERIMAGERFROMIMAGE_TXX_INCLUDED

template < const unsigned int TDimension, class TPixelType >
int ExtractSmallerImageFromImage< TDimension, TPixelType >
::Execute()
{

  unsigned int exitStatus = EXIT_SUCCESS;
  typename ImageReaderType::Pointer r;
  typename ImageWriterType::Pointer w;
  typename FilterType::Pointer f;

  r = ImageReaderType::New();
  r->SetFileName(m_InputImageFileName.c_str());
  r->Update();
  w = ImageWriterType::New();
  w->SetFileName(m_OutputImageFileName.c_str());

  f = FilterType::New();
  f->SetInput( r->GetOutput() );
  typename ImageTypeIn::RegionType inputRegion = 
                   r->GetOutput()->GetLargestPossibleRegion();
  typename ImageTypeIn::SizeType size = inputRegion.GetSize();
  typename ImageTypeIn::IndexType start = inputRegion.GetIndex();

  size[m_Axis] = 0;
  int bufferlen = 1 + static_cast<int>(vcl_ceil(1 + log10(fabs((float)m_HighIndex))));
  bufferlen += strlen(m_OutputImageFileName.c_str());
  char *buffer = new char[bufferlen];
  for (unsigned int i = m_LowIndex; i < m_HighIndex; i++) {
    memset(buffer, '\0', bufferlen);
    snprintf(buffer, bufferlen, m_OutputImageFileName.c_str(), i);
    w->SetFileName(buffer);

    start[m_Axis] = i;

    typename ImageTypeIn::RegionType desiredRegion;
    desiredRegion.SetSize( size );
    desiredRegion.SetIndex( start );

    f->SetExtractionRegion( desiredRegion );
    f->Update();
    w->SetInput(f->GetOutput());

    try {
      w->Update();
      w->UseCompressionOn();
      w->Write();
    } catch (itk::ExceptionObject &e) {
        std::cerr << "Caught ITK exception: " << e << std::endl;
        exitStatus = EXIT_FAILURE;
    }
  }
  delete [] buffer;

  return exitStatus;
}

#endif 

