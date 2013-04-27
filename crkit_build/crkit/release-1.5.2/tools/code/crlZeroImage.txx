
template <unsigned int Dimension, class PixelType>
int
crlZeroImage<Dimension, PixelType >
::Execute()
{

  typename ImageWriterType::Pointer w;
  w = ImageWriterType::New();
  if (!w) return EXIT_FAILURE;
  w->SetFileName( m_OutputFileName );

  if (!m_Image) {
    return EXIT_FAILURE;
  }

  m_Image->SetSpacing( m_Spacing );
  m_Image->SetOrigin( m_Origin);

  typename ImageType::RegionType region;
  region.SetSize( m_Size );
  region.SetIndex( m_Start );  // Note that regions assume a start index of
                               // zero if none are supplied.
  m_Image->SetRegions( region );
  m_Image->Allocate();

  //typename ImageType::PixelType  initialValue = 0.0;
  //m_Image->FillBuffer( initialValue );
  m_Image->FillBuffer(m_Value);

  w->SetInput (m_Image);

  try
  {
    w->UseCompressionOn();
    w->Update();
  }
  catch( itk::ExceptionObject & exp )
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
