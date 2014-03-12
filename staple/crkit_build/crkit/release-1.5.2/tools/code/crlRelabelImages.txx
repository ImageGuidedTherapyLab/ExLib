
#ifndef _CRL_RELABELIMAGES_TXX_INCLUDED
#define _CRL_RELABELIMAGES_TXX_INCLUDED 1

int RelabelImagesBase::countFloatsInString(const char *string)
/* return -1 on malformed string. 
 * return the count of number of floating point numbers.
 */
{
  char *end;
  const char *start = string;
  double d;
  int count = 0;
  while ((static_cast<unsigned char>(*start) != '\0') && isspace(static_cast<unsigned char>(*start))) { start++; }
  if (static_cast<unsigned char>(*start) == '\0') return -1; /* Don't ask to convert empty strings */  do {
    d = strtod(start, (char **)&end);
    if (end == start) {
      /* I want to parse strings of numbers with comments at the end */
      /* This return is executed when the next thing along can't be parsed */
      return count;
    }
    count++;   /* Count the number of floats */
    start = end;  /* Keep converting from the returned position. */
    while ((static_cast<unsigned char>(*start) != '\0') && isspace(static_cast<unsigned char>(*start))) { start++; }
  } while (static_cast<unsigned char>(*start) != '\0');
  return count; /* Success */
};

std::vector<double> * 
RelabelImagesBase::getFloatString(const char *string)
{
  std::vector<double> * v = new std::vector<double>();

  int numFloats = countFloatsInString(string);
  if (numFloats == -1) return NULL;
  if (numFloats == 0) return v;

  char *end;
  const char *start = string;
  double d;
  int count = 0;
  while ((static_cast<unsigned char>(*start) != '\0') && isspace(static_cast<unsigned char>(*start))) { start++; }
  do {
    d = strtod(start, (char **)&end);
    if (end == start) {
      /* Can't do any more conversions on this line */
      return v;
    }
    v->push_back(d);
    start = end;  /* Keep converting from the returned position. */
    while ((static_cast<unsigned char>(*start) != '\0') && isspace(static_cast<unsigned char>(*start))) { start++; }
    if (count == numFloats) return v; /* I don't care if there are more on
                        the line as long as I got what I wanted */
  } while (static_cast<unsigned char>(*start) != '\0');
  return v;
};

template < const unsigned int TImageDimension, class TPixelType >
int RelabelImages< TImageDimension, TPixelType > 
::Execute()
{
  // The image type is largely defined by the template parameters;
  typedef itk::Image< TPixelType, TImageDimension >         ImageType;

  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  typename ImageReaderType::Pointer reader1 = ImageReaderType::New();
  typename ImageReaderType::Pointer reader2 = ImageReaderType::New();

  if (m_Image1FileName.size() <= 0) {
    std::cerr << "Input 1 file name must be set." << std::endl;
  }
  if (m_Image2FileName.size() <= 0) {
    std::cerr << "Input 2 file name must be set." << std::endl;
  }
  if (m_OutputImageFileName.size() <= 0) {
    std::cerr << "Output image file name must be set." << std::endl;
  }
  reader1->SetFileName( m_Image1FileName );
  reader2->SetFileName( m_Image2FileName );

  try {
    reader1->Update();
    reader2->Update();
  } catch (itk::ExceptionObject &ex) {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }

  if (m_InLabels->size() <= 0) {
    std::cout << "The input labels must be set." << std::endl;
  }
  if (m_InLabels->size() != m_OutLabels->size()) {
    std::cout << "The number of input labels must equal the number of output labels." << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::ImageRegionIteratorWithIndex< ImageType >        IteratorType;

  typename ImageType::Pointer output = ImageType::New();

  /* Set the output to be the values of reader2 */
  output = reader2->GetOutput();
  output->Update();

  /* Now iterator over reader1 and compare the pixel values in the
   * tentative output.
   */
  IteratorType it( reader1->GetOutput(),
                               reader1->GetOutput()->GetRequestedRegion() );

  IteratorType out(output, reader2->GetOutput()->GetRequestedRegion());

  // This assumes the same iterators go over the same geometry the same way.
  for (it.GoToBegin(), out.GoToBegin(); !it.IsAtEnd(); ++it, ++out)
    {
      TPixelType start = it.Value();
      TPixelType outval = static_cast<TPixelType>( 0.0 );
      if (m_DefaultOutputValueSupplied) {
        outval = static_cast<TPixelType>( m_DefaultOutputValue );
      } else {
        outval = out.Value();  // preserve the value in the second input
                                   // unless the first input matches one
                                   // of the given labels.
      }
      for (unsigned int i = 0; i < m_InLabels->size(); i++) {
        if (start == (*m_InLabels)[i]) { 
          outval = static_cast<TPixelType>( (*m_OutLabels)[i] );
        }
      }
      out.Set( outval );
    }

  typedef itk::ImageFileWriter< ImageType > WriterType;

  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( m_OutputImageFileName );

  try
    {
    writer->SetInput(output);
    writer->UseCompressionOn( );
    writer->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;

};

#endif

