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

#ifndef _CRL_CHANGEIMAGEINFORMATION_H
#define _CRL_CHANGEIMAGEINFORMATION_H 1

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkChangeInformationImageFilter.h>

#include <iostream>

#include "tclap/CmdLine.h"
#include "configuration.h"

#ifndef PATH_MAX
#define PATH_MAX 32767
#endif

#ifdef WIN32
#define snprintf _snprintf
#endif

int GetImageType (std::string fileName,
                     itk::ImageIOBase::IOPixelType &pixelType,
                     itk::ImageIOBase::IOComponentType &componentType,
                     itk::ImageIOBase::Pointer &imageIO,
                     unsigned int &dimensionality )
{
  typedef itk::Image<unsigned char, 2> ImageType;
  itk::ImageFileReader<ImageType>::Pointer imageReader =
                                itk::ImageFileReader<ImageType>::New();
  imageReader->SetFileName(fileName.c_str());
  try {
    imageReader->UpdateOutputInformation();
  } catch ( itk::ExceptionObject &err )
  {
    std::cout << "ExceptionObject caught!" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }

  pixelType = imageReader->GetImageIO()->GetPixelType();
  componentType = imageReader->GetImageIO()->GetComponentType();
  imageIO = imageReader->GetImageIO();
  dimensionality = imageIO->GetNumberOfDimensions();
  return EXIT_SUCCESS;
};

namespace crl
{

class ITK_EXPORT ChangeImageInformationBase : public itk::LightObject
{
  // Defines the public API for running ChangeImageInformation.

  public:
  typedef ChangeImageInformationBase Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ChangeImageInformationBase, LightObject);

  ChangeImageInformationBase() { 
    m_UseReferenceImage = false;
  };

  ~ChangeImageInformationBase() { };

  void SetInputImageFileName(std::string name) { m_InputImageFileName = name; }
  void SetReferenceImageFileName(std::string name) { 
     m_ReferenceImageFileName = name; 
     m_UseReferenceImage = true;
  }
  void SetOutputImageFileName(std::string name) { m_OutputImageFileName = name; }
  void SetChangeSpacing(bool changeSpacing) { m_ChangeSpacing = changeSpacing; }
  void SetChangeDirection(bool changeDirection) { m_ChangeDirection = changeDirection; }
  void SetChangeOrigin(bool changeOrigin) { m_ChangeOrigin = changeOrigin; }
  void SetChangeRegion(bool changeRegion) { m_ChangeRegion = changeRegion; }

  virtual int Execute() = 0;

  protected:
  std::string m_InputImageFileName;
  std::string m_ReferenceImageFileName;
  bool m_UseReferenceImage;
  std::string m_OutputImageFileName;
  bool m_ChangeSpacing;
  bool m_ChangeDirection;
  bool m_ChangeOrigin;
  bool m_ChangeRegion;

  private:
  // purposely not implemented
  ChangeImageInformationBase(const ChangeImageInformationBase &);
  void operator=(const ChangeImageInformationBase &);
};

template <const unsigned int TDimension, class TPixelType >
class ITK_EXPORT ChangeImageInformation : 
      public ChangeImageInformationBase
{
public:
  /** Standard class typedefs. */
  typedef   ChangeImageInformation Self;
  typedef   ChangeImageInformationBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::Image< TPixelType, TDimension >  ImageTypeIn;
  typedef   itk::Image< TPixelType, TDimension >  ImageTypeOut;
  typedef   itk::ImageFileReader<ImageTypeIn> ImageReaderType;
  typedef   itk::ImageFileWriter<ImageTypeOut> ImageWriterType;
  typedef   itk::ChangeInformationImageFilter< ImageTypeIn > FilterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  ChangeImageInformation() { };

  ~ChangeImageInformation() {};

  int Execute();

private:
  // purposely not implemented
  ChangeImageInformation(const Self &); 
  void operator=(const Self &);

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlChangeImageInformation.txx"
#endif

} // end namespace crl

#endif

