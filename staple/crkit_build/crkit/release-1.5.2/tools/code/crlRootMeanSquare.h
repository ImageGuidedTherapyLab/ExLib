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

#ifndef _CRL_ROOTMEANSQUARE_H_
#define _CRL_ROOTMEANSQUARE_H_ 1

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <math.h>

class ITK_EXPORT crlRootMeanSquareBase : public itk::LightObject
{
  // Defines the public API for running
  //   crlRootMeanSquareBase
  // but has no templated arguments.  The templated arguments are provided by
  // classes that inherit this API.
  //  This allows us to create a pointer to a class of this type, and
  // then to change it at run time to match the types we actually need to
  // deal with.

public:
  typedef crlRootMeanSquareBase Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(crlRootMeanSquareBase, LightObject);

  crlRootMeanSquareBase() {
  };

  ~crlRootMeanSquareBase() {
  };

  void SetImage1FileName(std::string name) { m_Image1FileName = name; }
  void SetImage2FileName(std::string name) { m_Image2FileName = name; }

  virtual int Execute() = 0;

  virtual double GetRootMeanSquare() = 0;

protected:
  double m_RootMeanSquare;
  std::string m_Image1FileName;
  std::string m_Image2FileName;

private:
  crlRootMeanSquareBase(const
			    crlRootMeanSquareBase &); // purposely not implemented
  void operator=(const crlRootMeanSquareBase &); // purposely
  // not implemented.

};

template <const unsigned int TImageDimension, class TPixelType >
class ITK_EXPORT crlRootMeanSquare :
  public crlRootMeanSquareBase
{
public:
  /** Standard class typedefs. */
  typedef   crlRootMeanSquare Self;
  typedef   crlRootMeanSquareBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::Image< TPixelType, TImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  typedef   itk::ImageFileWriter<ImageType> ImageWriterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  crlRootMeanSquare() {
  };

  ~crlRootMeanSquare() {};

  // This won't make sense unless the class has been initialized and
  // executed properly.
  double GetRootMeanSquare() { 
    return m_RootMeanSquare;
  };

  int Execute();

private:
  crlRootMeanSquare(const Self&); // purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

// Partial specialization for DiffusionTensor3D< double >
template <const unsigned int TImageDimension >
class ITK_EXPORT crlRootMeanSquare<TImageDimension, 
                         itk::DiffusionTensor3D<double> > :
  public crlRootMeanSquareBase
{
public:

  typedef itk::DiffusionTensor3D< double > TPixelType;

  /** Standard class typedefs. */
  typedef   crlRootMeanSquare Self;
  typedef   crlRootMeanSquareBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::Image< TPixelType, TImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  typedef   itk::ImageFileWriter<ImageType> ImageWriterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  crlRootMeanSquare() {
  };

  ~crlRootMeanSquare() {};

  // This won't make sense unless the class has been initialized and
  // executed properly.
  double GetRootMeanSquare() { 
    return m_RootMeanSquare;
  };

  int Execute();

private:
  crlRootMeanSquare(const Self&); // purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

// Partial specialization for DiffusionTensor3D< float >
template <const unsigned int TImageDimension >
class ITK_EXPORT crlRootMeanSquare<TImageDimension, 
                         itk::DiffusionTensor3D<float> > :
  public crlRootMeanSquareBase
{
public:

  typedef itk::DiffusionTensor3D< double > TPixelType;

  /** Standard class typedefs. */
  typedef   crlRootMeanSquare Self;
  typedef   crlRootMeanSquareBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::Image< TPixelType, TImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  typedef   itk::ImageFileWriter<ImageType> ImageWriterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  crlRootMeanSquare() {
  };

  ~crlRootMeanSquare() {};

  // This won't make sense unless the class has been initialized and
  // executed properly.
  double GetRootMeanSquare() { 
    return m_RootMeanSquare;
  };

  int Execute();

private:
  crlRootMeanSquare(const Self&); // purposely not implemented
  void operator=(const Self&); //purposely not implemented

};



#ifndef ITK_MANUAL_INSTANTIATION
#include "crlRootMeanSquare.txx"
#endif


#endif
