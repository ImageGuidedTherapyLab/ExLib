
#ifndef _CRL_GAUSSIANSMOOTH_TXX_INCLUDED
#define _CRL_GAUSSIANSMOOTH_TXX_INCLUDED

template < const unsigned int TDimension, class TPixelType >
int GaussianSmooth< TDimension, TPixelType >
::Execute()
{
  unsigned int exitStatus = EXIT_SUCCESS;
  typename ImageReaderType::Pointer reader;
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
  filter->SetNormalizeAcrossScale( m_NormalizeAcrossScale );
  filter->SetSigma( m_Sigma );
  filter->Update();

  if (m_Clamp) {
    typename ImageTypeOut::Pointer image = filter->GetOutput();
    itk::ImageRegionIterator< ImageTypeOut > iter(image,
                                   image->GetLargestPossibleRegion() );
    for ( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
        {
        if ( iter.Value() > m_MaxValue)
          iter.Value() = m_MaxValue;
        else if ( iter.Value() < m_MinValue)
          iter.Value() = m_MinValue;
        }
  }

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

