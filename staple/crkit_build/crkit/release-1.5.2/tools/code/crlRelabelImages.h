/*
 * Copyright 2007,2008 Simon Warfield simon.warfield@childrens.harvard.edu
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
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#ifndef _CRLRELABELIMAGES_H
#define _CRLRELABELIMAGES_H

#include <vector>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkImageRegionIteratorWithIndex.h>

namespace crl
{

class ITK_EXPORT RelabelImagesBase : public itk::LightObject
{
  // Defines the public API for running crlRelabelImages .
  // The templated arguments are provided by classes that inherit this API.
  //  This allows us to create a pointer to a class of this type, and
  // then to change it at run time to match the types we actually need to
  // deal with.

  public: 
  typedef RelabelImagesBase Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(RelabelImagesBase, LightObject);

  RelabelImagesBase() { m_DefaultOutputValueSupplied = false; };
  ~RelabelImagesBase() { };

  void SetImage1FileName(std::string name) { m_Image1FileName = name; }
  void SetImage2FileName(std::string name) { m_Image2FileName = name; }
  void SetOutputImageFileName(std::string name) { m_OutputImageFileName = name; }
  void SetDefaultOutputValue(double val) { 
    m_DefaultOutputValue = val; 
    m_DefaultOutputValueSupplied = true; 
  };

  int countFloatsInString(const char *string);
  std::vector<double> * getFloatString(const char *string);

  void SetInputLabels(std::string labels) { 
    m_InLabels = getFloatString(labels.c_str());
  }
  void SetOutputLabels(std::string labels) { 
    m_OutLabels = getFloatString(labels.c_str());
  }

  virtual int Execute() = 0;

  protected:
  std::string m_Image1FileName;
  std::string m_Image2FileName;
  std::string m_OutputImageFileName;
  std::vector<double> *m_InLabels;
  std::vector<double> *m_OutLabels;
  double m_DefaultOutputValue;
  bool m_DefaultOutputValueSupplied;

private:
  RelabelImagesBase(const RelabelImagesBase &); // purposely not implemented
  void operator=(const RelabelImagesBase &); // purposely not implemented.

};

template <const unsigned int TImageDimension, class TPixelType >
class ITK_EXPORT RelabelImages : public RelabelImagesBase
{
public:
  /** Standard class typedefs. */
  typedef   RelabelImages Self;
  typedef   RelabelImagesBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::Image< TPixelType, TImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  typedef   itk::ImageFileWriter<ImageType> ImageWriterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  RelabelImages() { };

  ~RelabelImages() {};

  int Execute();

private:
  RelabelImages(const Self&); // purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlRelabelImages.txx"
#endif

} // end namespace crl

#endif

