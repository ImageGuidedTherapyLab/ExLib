/*
 * Copyright (c) 2008-2011 Children's Hospital Boston.
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 * http://www.crl.med.harvard.edu
 * For more information contact: simon.warfield@childrens.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/


template< unsigned int VDimension>
int MSTAPLE<VDimension>::Execute()
{
  int i;

  typename itk::ImageFileReader<InputImageType>::Pointer  reader;
  typename itk::ImageFileWriter<OutputImageType>::Pointer writer
    = itk::ImageFileWriter<OutputImageType>::New();

  int numberOfFiles = m_Files.size();

  // Set the inputs
  bool sizeInitialized = false;
  typename InputImageType::SizeType initialSize;
  initialSize.Fill( 0 );
  typename InputImageType::SizeType currentSize;
  for (i = 0; i < numberOfFiles; i++)
    {
      try {
	// Check that the geometry of all the segmentations is the same.
	reader = itk::ImageFileReader<InputImageType>::New();
	// std::cout << "Reading file " << m_Files[i] << std::endl;
	reader->SetFileName( m_Files[i].c_str() );
	reader->Update();
	currentSize = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
	if (!sizeInitialized) {
	  initialSize = currentSize;
	  sizeInitialized = true;
	}
	if (currentSize != initialSize) {
	  std::cerr << "All input segmentations must have the same size."
		    << std::endl;
	  std::cerr << "File " << m_Files[i] << " has a different size."
		    << std::endl;
	  return -1;
	}
	m_MSTAPLE->SetInput(i, reader->GetOutput());
      } catch (itk::ExceptionObject &e) {
	std::cerr << e << std::endl;
	return -1;
      }
    }

  try
    {
      writer->SetFileName( m_OutputFile.c_str() );
      writer->SetInput(m_MSTAPLE->GetOutput());
      if ( m_MSTAPLE->GetUseWriteCompression() ) {
	writer->UseCompressionOn();
      }
      writer->Update();
    }
  catch( itk::ExceptionObject &e )
    {
      std::cerr << e << std::endl;
      return -2;
    }

  return 0;
}
