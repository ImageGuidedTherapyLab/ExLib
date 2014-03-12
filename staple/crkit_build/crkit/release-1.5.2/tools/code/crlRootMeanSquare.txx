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


#ifndef _CRL_ROOTMEANSQUARE_TXX_INCLUDED
#define _CRL_ROOTMEANSQUARE_TXX_INCLUDED 1

template <const unsigned int TImageDimension, class TPixelType >
int crlRootMeanSquare< TImageDimension, TPixelType >
::Execute()
{
  // The image type is largely defined by the template parameters;
  typedef itk::Image< TPixelType, TImageDimension >         ImageType;

  typedef itk::ImageFileReader< ImageType > VolumeReaderType;
  typename VolumeReaderType::Pointer volumereader1 = VolumeReaderType::New();
  typename VolumeReaderType::Pointer volumereader2 = VolumeReaderType::New();

  volumereader1->SetFileName( m_Image1FileName );
  volumereader2->SetFileName( m_Image2FileName );

  try {
    volumereader1->Update();
    volumereader2->Update();
  } catch (itk::ExceptionObject &ex) {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }

  typedef   itk::ImageRegionConstIteratorWithIndex<ImageType> ImageIterator;

  typename ImageType::Pointer image1 = volumereader1->GetOutput();
  typename ImageType::Pointer image2 = volumereader2->GetOutput();

  m_RootMeanSquare = 0.0;
  unsigned int count = 0;

  ImageIterator iter1(image1,  image1->GetLargestPossibleRegion());
  for (iter1.GoToBegin(); !iter1.IsAtEnd(); ++iter1) {
    typename ImageType::IndexType index = iter1.GetIndex();
    TPixelType p1 = image1->GetPixel(index);
    TPixelType p2 = image2->GetPixel(index);
    m_RootMeanSquare += (p1 - p2)*(p1-p2);
    count++;
  }
  m_RootMeanSquare = sqrt(m_RootMeanSquare / static_cast<float>(count));

  return 0; // success

};

// Implementation of partial specialization for itk::DiffusionTensor3D<double>
template < const unsigned int TImageDimension >
int crlRootMeanSquare< TImageDimension , itk::DiffusionTensor3D<double> >
::Execute()
{
  typedef itk::DiffusionTensor3D<double> TPixelType;

  // The image type is largely defined by the template parameters;
  typedef itk::Image< TPixelType, TImageDimension >         ImageType;

  typedef itk::ImageFileReader< ImageType > VolumeReaderType;
  typename VolumeReaderType::Pointer volumereader1 = VolumeReaderType::New();
  typename VolumeReaderType::Pointer volumereader2 = VolumeReaderType::New();

  volumereader1->SetFileName( m_Image1FileName );
  volumereader2->SetFileName( m_Image2FileName );

  try {
    volumereader1->Update();
    volumereader2->Update();
  } catch (itk::ExceptionObject &ex) {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }

  typedef   itk::ImageRegionConstIteratorWithIndex<ImageType> ImageIterator;

  typename ImageType::Pointer image1 = volumereader1->GetOutput();
  typename ImageType::Pointer image2 = volumereader2->GetOutput();

  m_RootMeanSquare = 0.0;
  ImageIterator iter1(image1,  image1->GetLargestPossibleRegion());

  double norm = 0.0;
  unsigned int count = 0;
  double v1,v2,diff;
  for (iter1.GoToBegin(); !iter1.IsAtEnd(); ++iter1) {
    typename ImageType::IndexType index = iter1.GetIndex();
    TPixelType p1 = image1->GetPixel(index);
    TPixelType p2 = image2->GetPixel(index);
    for (unsigned int i = 0; i < TImageDimension; i++) {
      for (unsigned int j = 0; j < TImageDimension; j++) {
        v1 = p1(i,j);
        v2 = p2(i,j);
        diff = v1 - v2;
        norm += (diff*diff);
      }
    }
    count++;
    m_RootMeanSquare += norm;
  }
  m_RootMeanSquare = sqrt(m_RootMeanSquare / static_cast<float>(count));

  return 0; // success

};

// Implementation of partial specialization for itk::DiffusionTensor3D<float>
template < const unsigned int TImageDimension >
int crlRootMeanSquare< TImageDimension , itk::DiffusionTensor3D<float> >
::Execute()
{
  typedef itk::DiffusionTensor3D<float> TPixelType;

  // The image type is largely defined by the template parameters;
  typedef itk::Image< TPixelType, TImageDimension >         ImageType;

  typedef itk::ImageFileReader< ImageType > VolumeReaderType;
  typename VolumeReaderType::Pointer volumereader1 = VolumeReaderType::New();
  typename VolumeReaderType::Pointer volumereader2 = VolumeReaderType::New();

  volumereader1->SetFileName( m_Image1FileName );
  volumereader2->SetFileName( m_Image2FileName );

  try {
    volumereader1->Update();
    volumereader2->Update();
  } catch (itk::ExceptionObject &ex) {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }

  typedef   itk::ImageRegionConstIteratorWithIndex<ImageType> ImageIterator;

  typename ImageType::Pointer image1 = volumereader1->GetOutput();
  typename ImageType::Pointer image2 = volumereader2->GetOutput();

  m_RootMeanSquare = 0.0;
  ImageIterator iter1(image1,  image1->GetLargestPossibleRegion());

  double norm = 0.0;
  unsigned int count = 0;
  double v1,v2,diff;
  for (iter1.GoToBegin(); !iter1.IsAtEnd(); ++iter1) {
    typename ImageType::IndexType index = iter1.GetIndex();
    TPixelType p1 = image1->GetPixel(index);
    TPixelType p2 = image2->GetPixel(index);
    for (unsigned int i = 0; i < TImageDimension; i++) {
      for (unsigned int j = 0; j < TImageDimension; j++) {
        v1 = p1(i,j);
        v2 = p2(i,j);
        diff = v1 - v2;
        norm += (diff*diff);
      }
    }
    count++;
    m_RootMeanSquare += norm;
  }
  m_RootMeanSquare = sqrt(m_RootMeanSquare / static_cast<float>(count));

  return 0; // success

};

#endif

