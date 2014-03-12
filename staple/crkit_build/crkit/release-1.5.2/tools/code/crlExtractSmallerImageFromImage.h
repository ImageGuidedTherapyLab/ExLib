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

#ifndef _CRL_EXTRACTSMALLERIMAGEFROMIMAGE_H
#define _CRL_EXTRACTSMALLERIMAGEFROMIMAGE_H 1

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>

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

class ITK_EXPORT ExtractSmallerImageFromImageBase : public itk::LightObject
{
  // Defines the public API for running ExtractSmallerImageFromImage .

  public:
  typedef ExtractSmallerImageFromImageBase Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ExtractSmallerImageFromImageBase, LightObject);

  ExtractSmallerImageFromImageBase() { };
  ~ExtractSmallerImageFromImageBase() { };

  void SetInputImageFileName(std::string name) { m_InputImageFileName = name; }
  void SetOutputImageFileName(std::string name) { m_OutputImageFileName = name; }
  void SetAxis(unsigned int axis) { m_Axis = axis; }
  void SetLowIndex(unsigned int index) { m_LowIndex = index; }
  void SetHighIndex(unsigned int index) { m_HighIndex = index; }

  virtual int Execute() = 0;

  protected:
  std::string m_InputImageFileName;
  std::string m_OutputImageFileName;
  unsigned int m_Axis;
  unsigned int m_LowIndex;
  unsigned int m_HighIndex;

  private:
  // purposely not implemented
  ExtractSmallerImageFromImageBase(const ExtractSmallerImageFromImageBase &);
  void operator=(const ExtractSmallerImageFromImageBase &);
};

template <const unsigned int TDimension, class TPixelType >
class ITK_EXPORT ExtractSmallerImageFromImage : 
      public ExtractSmallerImageFromImageBase
{
public:
  /** Standard class typedefs. */
  typedef   ExtractSmallerImageFromImage Self;
  typedef   ExtractSmallerImageFromImageBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::Image< TPixelType, TDimension >  ImageTypeIn;
  typedef   itk::Image< TPixelType, TDimension - 1 >  ImageTypeOut;
  typedef   itk::ImageFileReader<ImageTypeIn> ImageReaderType;
  typedef   itk::ImageFileWriter<ImageTypeOut> ImageWriterType;
  typedef   itk::ExtractImageFilter< ImageTypeIn, ImageTypeOut > FilterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  ExtractSmallerImageFromImage() { m_Axis = TDimension - 1; };

  ~ExtractSmallerImageFromImage() {};

  int Execute();

private:
  // purposely not implemented
  ExtractSmallerImageFromImage(const Self &); 
  void operator=(const Self &);

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlExtractSmallerImageFromImage.txx"
#endif

} // end namespace crl

#endif

