
template <unsigned int Dimension, class PixelType>
int
crlBinaryMorphology<Dimension, PixelType>
::Execute()
{

  typename ImageReaderType::Pointer reader = ImageReaderType::New();
  typename ImageWriterType::Pointer writer = ImageWriterType::New();

  typename ErodeFilterType::Pointer   binaryErode   = ErodeFilterType::New();
  typename DilateFilterType::Pointer  binaryDilate  = DilateFilterType::New();
  // ClosingFilterType::Pointer binaryClosing = ClosingFilterType::New();
  // OpeningFilterType::Pointer binaryOpening = OpeningFilterType::New();

  binaryErode->SetKernel(   structuringElement );
  binaryDilate->SetKernel(  structuringElement );
  // binaryClosing->SetKernel( structuringElement );
  // binaryOpening->SetKernel( structuringElement );

  reader->SetFileName( m_InputFileName );
  writer->SetFileName( m_OutputFileName );

  if (m_Operation == 1) {        // dilate
    binaryDilate->SetInput( reader->GetOutput() );
    binaryDilate->SetDilateValue( static_cast<PixelType>(m_ForegroundValue) );
    writer->SetInput( binaryDilate->GetOutput() );
  } else if (m_Operation == 2) { // erode
    binaryErode->SetInput( reader->GetOutput() );
    binaryErode->SetErodeValue( static_cast<PixelType>(m_ForegroundValue) );
    writer->SetInput( binaryErode->GetOutput() );
  } else if (m_Operation == 3) { // closing
    binaryDilate->SetInput( reader->GetOutput() );
    binaryDilate->SetDilateValue( static_cast<PixelType>(m_ForegroundValue) );
    binaryErode->SetInput( binaryDilate->GetOutput() );
    binaryErode->SetErodeValue( static_cast<PixelType>(m_ForegroundValue) );
    writer->SetInput( binaryErode->GetOutput() );
  } else if (m_Operation == 4) { // opening
    binaryErode->SetInput( reader->GetOutput() );
    binaryErode->SetErodeValue( static_cast<PixelType>(m_ForegroundValue) );
    binaryDilate->SetInput( binaryErode->GetOutput() );
    binaryDilate->SetDilateValue( static_cast<PixelType>(m_ForegroundValue) );
    writer->SetInput( binaryDilate->GetOutput() );
  }

  try
  {
    writer->Update();
  }
  catch( itk::ExceptionObject & exp )
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

template <unsigned int Dimension, class PixelType>
void 
crlBinaryMorphology<Dimension, PixelType>
::SetRadius(unsigned int radius)
{
    structuringElement.SetRadius( radius ); // ball of given radius
    structuringElement.CreateStructuringElement(); // create the ball
}

