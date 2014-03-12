/*
 * Copyright (c) 2008, Children's Hospital Boston.
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

#ifndef _CRL_IMAGEADDMULTIPLYADD_TXX_INCLUDED
#define _CRL_IMAGEADDMULTIPLYADD_TXX_INCLUDED 1

template <const unsigned int TImageDimension, class TPixelType >
int crlImageAddMultiplyAdd< TImageDimension, TPixelType >
::Execute()
{
  // The image type is largely defined by the template parameters;
  typedef itk::Image< TPixelType, TImageDimension >         ImageType;

  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  typedef   itk::ImageRegionIterator<ImageType> ImageIterator;
  typedef   itk::ImageFileWriter< ImageType >    ImageWriterType;

  typename ImageReaderType::Pointer inputImageReader = ImageReaderType::New();

  inputImageReader->SetFileName( m_InputImageFileName );

  try {
    inputImageReader->Update();
  } catch (itk::ExceptionObject &ex) {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }

  // Modify this to be scalar values:
  double add1 = 0.0;
  double mult1 = 0.0;
  double add2 = 0.0;

  // Convert from string representation to TPixelType.
  // This needs to be done carefully for scalar, vector and tensor types.
  add1 = atof(m_Add1.c_str());
  mult1 = atof(m_Multiply.c_str());
  add2 = atof(m_Add2.c_str());

std::cout << "add1 " << add1 << std::endl;
std::cout << "mult1 " << mult1 << std::endl;
std::cout << "add2 " << add2 << std::endl;

  typename ImageType::Pointer image1;
  image1 = ImageType::New();
  image1->SetSpacing( inputImageReader->GetOutput()->GetSpacing() );
  image1->SetOrigin( inputImageReader->GetOutput()->GetOrigin() );
  image1->SetDirection( inputImageReader->GetOutput()->GetDirection() );
  typename ImageType::RegionType region;
  region.SetSize( inputImageReader->GetOutput()->GetLargestPossibleRegion().GetSize());
  region.SetIndex( inputImageReader->GetOutput()->GetLargestPossibleRegion().GetIndex());
  image1->SetRegions( region );
  image1->Allocate();

  typename ImageWriterType::Pointer writer;
  writer = ImageWriterType::New();
  writer->SetFileName( m_OutputImageFileName );

  // using itk iterators. Fast and beautiful.
  // We rely on these two iterators traversing the same sized volumes in
  // the same order.
  TPixelType outval;
  ImageIterator iter1(image1,  image1->GetLargestPossibleRegion());
  ImageIterator iterInput(inputImageReader->GetOutput(),  
       inputImageReader->GetOutput()->GetLargestPossibleRegion());
  for (iter1.GoToBegin(),iterInput.GoToBegin() ; 
          !iter1.IsAtEnd(); ++iter1, ++iterInput) {
    TPixelType p1 = iterInput.Get();
    // The calculation is : 
    // outval = (p1 + add1)*mult1 + add2;
    // but splitting it out in the following manner allows the compiler to
    // instantiate different operators for scalars, vector and tensors.
    outval = p1 + add1;
    outval = outval * mult1;
    outval = outval + add2;
    iter1.Set(outval);
  }
 
  writer->SetInput (image1);
  writer->UseCompressionOn( );
  
  try
  {
    writer->Update();
  }
  catch( itk::ExceptionObject & exp )
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
  }

  return 0; // success
}

#endif

