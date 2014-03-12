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

#ifndef _CRL_IMAGEADDMULTIPLYADD_
#define _CRL_IMAGEADDMULTIPLYADD_ 1

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>

#include <itkDiffusionTensor3D.h>
#include <itkNumericTraits.h>

class ITK_EXPORT crlImageAddMultiplyAddBase : public itk::LightObject
{
  // This defines the public API for running
  //   crlImageAddMultiplyAdd

  public:
    typedef crlImageAddMultiplyAddBase Self;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(crlImageAddMultiplyAddBase, LightObject);

    crlImageAddMultiplyAddBase() {
    };

    ~crlImageAddMultiplyAddBase() {
    };

    void SetInputImageFileName(std::string name) { 
      m_InputImageFileName = name;
    };

    void SetAdd1(std::string name) { m_Add1 = name; };
    void SetAdd2(std::string name) { m_Add2 = name; };
    void SetMultiply(std::string name) { m_Multiply = name; };

    void SetOutputImageFileName(std::string name) { 
      m_OutputImageFileName = name;
    };

    virtual int Execute() = 0;

  protected:
    std::string m_InputImageFileName;
    std::string m_OutputImageFileName;
    std::string m_Add1;
    std::string m_Add2;
    std::string m_Multiply;

  private:
    crlImageAddMultiplyAddBase(const crlImageAddMultiplyAddBase &); // not impl.
    void operator=(const crlImageAddMultiplyAddBase &); // not impl.

};

template <const unsigned int TImageDimension, class TPixelType >
class ITK_EXPORT crlImageAddMultiplyAdd : public crlImageAddMultiplyAddBase
{
  public:
  /** Standard class typedefs. */
  typedef   crlImageAddMultiplyAdd Self;
  typedef   crlImageAddMultiplyAddBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::Image< TPixelType, TImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  typedef   itk::ImageFileWriter<ImageType> ImageWriterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  crlImageAddMultiplyAdd() {
  };

  ~crlImageAddMultiplyAdd() {};

  int Execute();

  void GetFirstAdditionValue(TPixelType a1) { m_FirstAdditionValue = a1; };
  void GetSecondAdditionValue(TPixelType a1) { m_SecondAdditionValue = a1; };
  void GetMultiplicationValue(TPixelType a1) { m_MultiplicationValue = a1; };

private:
  crlImageAddMultiplyAdd(const Self&); // purposely not implemented
  void operator=(const Self&); //purposely not implemented

  TPixelType m_FirstAdditionValue;
  TPixelType m_SecondAdditionValue;
  TPixelType m_MultiplicationValue;

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlImageAddMultiplyAdd.txx"
#endif

#endif

