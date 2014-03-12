
#ifndef _CRL_ZEROIMAGE_INCLUDED
#define _CRL_ZEROIMAGE_INCLUDED 1

#include <itkImage.h>
#include <itkImageFileWriter.h>

#include <iostream>

class ITK_EXPORT crlZeroImageBase : public itk::LightObject
{
  public:
  typedef crlZeroImageBase Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(crlZeroImageBase, LightObject);

  crlZeroImageBase() {};
  ~crlZeroImageBase() {};

  void SetOutputFileName(std::string name) {
    m_OutputFileName = name;
  };

  virtual void SetStartIndex(unsigned int i, unsigned int index) = 0;
  virtual void SetSize(unsigned int i, unsigned int size) = 0;
  virtual void SetSpacing(unsigned int i, float spacing) = 0;
  virtual void SetOrigin(unsigned int i, float origin) = 0;
  virtual void SetValue(double value) = 0;

  virtual int Execute() = 0;

protected:
  std::string m_OutputFileName;

private:
  crlZeroImageBase(const crlZeroImageBase &); // purposely not implemented
  void operator=(const crlZeroImageBase &); // purposely not implemented

};

template <unsigned int Dimension, class PixelType>
class ITK_EXPORT crlZeroImage : public crlZeroImageBase
{
public:
  /** Standard class typedefs. */
  typedef crlZeroImage Self;
  typedef crlZeroImageBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef itk::Image< PixelType, Dimension> ImageType;
  typedef   itk::ImageFileWriter< ImageType >    ImageWriterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  crlZeroImage()  { m_Image = ImageType::New();  m_Value = 0.0;};
  ~crlZeroImage() { m_Image = 0; };

  void SetStartIndex(unsigned int i, unsigned int index) {
    m_Start[i] = index;
  };

  void SetSize(unsigned int i, unsigned int size) {
    m_Size[i] = size;
  };

  void SetSpacing(unsigned int i, float spacing) {
    m_Spacing[i] = spacing;
  };

  void SetOrigin(unsigned int i, float origin) {
    m_Origin[i] = origin;
  };

  void SetValue(double value) {
    m_Value = value;
  };

  int Execute();

private:
  crlZeroImage(const Self&); // purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typename ImageType::Pointer m_Image;
  typename ImageType::IndexType m_Start;
  typename ImageType::SizeType  m_Size;
  typename ImageType::SpacingType m_Spacing;
  typename ImageType::PointType m_Origin;
  typename ImageType::PixelType m_Value;

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlZeroImage.txx"
#endif

#endif
