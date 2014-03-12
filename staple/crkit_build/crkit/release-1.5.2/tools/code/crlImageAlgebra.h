
#ifndef _CRL_IMAGEALGEBRA_INCLUDED
#define _CRL_IMAGEALGEBRA_INCLUDED 1

#include <itkImage.h>
#include <itkVector.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkSubtractImageFilter.h>

#include <iostream>

class ITK_EXPORT crlImageAlgebraBase : public itk::LightObject
{
public:
  typedef crlImageAlgebraBase Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(crlImageAlgebraBase, LightObject);

  crlImageAlgebraBase() {};
  ~crlImageAlgebraBase() {};

  void SetInput1FileName(std::string name) {
    m_Input1FileName = name;
  };
  void SetInput2FileName(std::string name) {
    m_Input2FileName = name;
  };
  void SetOutputFileName(std::string name) {
    m_OutputFileName = name;
  };
  void SetOperationName(std::string name) {
    m_OperationName = name;
  };

  virtual int Execute() = 0;

protected:
  std::string m_Input1FileName;
  std::string m_Input2FileName;
  std::string m_OutputFileName;
  std::string m_OperationName;

private:
  crlImageAlgebraBase(const crlImageAlgebraBase&); // purposely not implemented
  void operator=(const crlImageAlgebraBase &); // purposely not implemented
};

template <unsigned int Dimension, class PixelType>
class ITK_EXPORT crlImageAlgebra : public crlImageAlgebraBase
{
public:
  /** Standard class typedefs. */
  typedef crlImageAlgebra Self;
  typedef crlImageAlgebraBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

    typedef itk::Image< PixelType, Dimension> ImageType;
    typedef   itk::ImageFileReader< ImageType >    ImageReaderType;
    typedef   itk::ImageRegionConstIterator<ImageType> ImageIterator;
    typedef   itk::ImageFileWriter< ImageType >    ImageWriterType;

    // Filters that this program knows how to use:
    typedef itk::AddImageFilter<ImageType, ImageType, ImageType> AddFilterType;
    typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> 
                 MultiplyFilterType;
    typedef itk::SubtractImageFilter<ImageType, ImageType, ImageType> 
                 SubtractFilterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  crlImageAlgebra()
  {

    m_Reader1 = ImageReaderType::New();
    m_Reader2 = ImageReaderType::New();
    m_Writer = ImageWriterType::New();

    m_AddFilter = 0;
    m_MultiplyFilter = 0;
    m_SubtractFilter = 0;
    
  }

  ~crlImageAlgebra()
  {
  }

  // It is a bit strange that I can't just get a handle on an opaque
  // BinaryFunctorImageFilter and then point that at the appropriate sub-type
  // of filter.
  int Execute();

  protected :

    typename ImageReaderType::Pointer m_Reader1;
    typename ImageReaderType::Pointer m_Reader2;
    typename ImageWriterType::Pointer m_Writer;

    typename ImageType::Pointer m_InImage1;
    typename ImageType::Pointer  m_InImage2;
    typename ImageType::Pointer  m_OutImage;

    typename AddFilterType::Pointer m_AddFilter;
    typename MultiplyFilterType::Pointer m_MultiplyFilter;
    typename SubtractFilterType::Pointer m_SubtractFilter;

private:
  crlImageAlgebra(const Self&); // purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlImageAlgebra.txx"
#endif

#endif

