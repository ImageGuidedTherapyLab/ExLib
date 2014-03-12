
template <unsigned int Dimension, class PixelType>
int
crlImageAlgebra<Dimension, PixelType >
::Execute()
  // It is a bit strange that I can't just get a handle on an opaque
  // BinaryFunctorImageFilter and then point that at the appropriate sub-type
  // of filter.
{

    m_Reader1->SetFileName( m_Input1FileName );
    m_Reader2->SetFileName( m_Input2FileName );
    m_Writer->SetFileName( m_OutputFileName);
    m_Writer->UseCompressionOn();

    // Try to determine what type of filter is to be executed and then run it
    if ( m_OperationName.find("add") != std::string::npos ) {

      m_AddFilter = AddFilterType::New();
      m_AddFilter->SetInput1( m_Reader1->GetOutput() );
      m_AddFilter->SetInput2( m_Reader2->GetOutput() );

      try
      {
        m_AddFilter->Update();
        m_Writer->SetInput( m_AddFilter->GetOutput() );
        m_Writer->Update();
      } catch ( itk::ExceptionObject &err) {
        std::cout << "Exception caught writing file." << std::endl;
        std::cout << err << std::endl;
        return 1;
      }

    } else
    if ( m_OperationName.find("subtract") != std::string::npos ) {
      m_SubtractFilter = SubtractFilterType::New();
      m_SubtractFilter->SetInput1( m_Reader1->GetOutput() );
      m_SubtractFilter->SetInput2( m_Reader2->GetOutput() );

      try
      {
        m_SubtractFilter->Update();
        m_Writer->SetInput( m_SubtractFilter->GetOutput() );
        m_Writer->Update();
      } catch ( itk::ExceptionObject &err) {
        std::cout << "Exception caught writing file." << std::endl;
        std::cout << err << std::endl;
        return 1;
      }
    } else
    if ( m_OperationName.find("multiply") != std::string::npos ) {
      m_MultiplyFilter = MultiplyFilterType::New();
      m_MultiplyFilter->SetInput1( m_Reader1->GetOutput() );
      m_MultiplyFilter->SetInput2( m_Reader2->GetOutput() );

      try
      {
        m_MultiplyFilter->Update();
        m_Writer->SetInput( m_MultiplyFilter->GetOutput() );
        m_Writer->Update();
      } catch ( itk::ExceptionObject &err) {
        std::cout << "Exception caught writing file." << std::endl;
        std::cout << err << std::endl;
        return 1;
      }
    } else
#if 0 // Turning off filters that don't work for vectors.
    if ( m_OperationName.find("divide") != std::string::npos ) {
      m_DivideFilter = DivideFilterType::New();
      m_DivideFilter->SetInput1( m_Reader1->GetOutput() );
      m_DivideFilter->SetInput2( m_Reader2->GetOutput() );

      try
      {
        m_DivideFilter->Update();
        m_Writer->SetInput( m_DivideFilter->GetOutput() );
        m_Writer->Update();
      } catch ( itk::ExceptionObject &err) {
        std::cout << "Exception caught writing file." << std::endl;
        std::cout << err << std::endl;
        return 1;
      }
    } else
    if ( m_OperationName.find("absolutevaluedifference") != std::string::npos ) {
      m_AbsoluteValueDifferenceFilter = AbsoluteValueDifferenceFilterType::New();
      m_AbsoluteValueDifferenceFilter->SetInput1( m_Reader1->GetOutput() );
      m_AbsoluteValueDifferenceFilter->SetInput2( m_Reader2->GetOutput() );

      try
      {
        m_AbsoluteValueDifferenceFilter->Update();
        m_Writer->SetInput( m_AbsoluteValueDifferenceFilter->GetOutput() );
        m_Writer->Update();
      } catch ( itk::ExceptionObject &err) {
        std::cout << "Exception caught writing file." << std::endl;
        std::cout << err << std::endl;
        return 1;
      }
    } else
    if ( m_OperationName.find("squareddifference") != std::string::npos ) {
      m_SquaredDifferenceFilter = SquaredDifferenceFilterType::New();
      m_SquaredDifferenceFilter->SetInput1( m_Reader1->GetOutput() );
      m_SquaredDifferenceFilter->SetInput2( m_Reader2->GetOutput() );

      try
      {
        m_SquaredDifferenceFilter->Update();
        m_Writer->SetInput( m_SquaredDifferenceFilter->GetOutput() );
        m_Writer->Update();
      } catch ( itk::ExceptionObject &err) {
        std::cout << "Exception caught writing file." << std::endl;
        std::cout << err << std::endl;
        return 1;
      }
    } else 
#endif
    {
      std::cerr << "Don't know how to carry out " << m_OperationName.c_str() << std::endl;
      std::cerr << "Valid options are : " ;
      std::cerr << "add " ;
      std::cerr << "subtract " ;
      std::cerr << "multiply " ;
      std::cerr << std::endl;
    }

  return 0;
}

