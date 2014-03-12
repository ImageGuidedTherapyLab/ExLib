/* 
 * Neil Weisenfeld neil.weisenfeld@childrens.harvard.edu
 */

#ifndef _CRL_GAUSSIANSMOOTH_H
#define _CRL_GAUSSIANSMOOTH_H 1

#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>

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

class ITK_EXPORT GaussianSmoothBase : public itk::LightObject
{
  // Defines the public API for running GaussianSmooth.

  public:
  typedef GaussianSmoothBase Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, LightObject);

  GaussianSmoothBase() {
    m_NormalizeAcrossScale = false;
    m_Sigma = 1.0;
    m_Clamp = false;
    m_MinValue = 0.0;
    m_MaxValue = 0.0;
  };

  ~GaussianSmoothBase() { };

  void SetInputImageFileName(std::string name) { m_InputImageFileName = name; }
  void SetOutputImageFileName(std::string name) { m_OutputImageFileName = name;}
  void SetNormalizeAcrossScale(bool normalize) { 
    m_NormalizeAcrossScale = normalize;
  }

  void SetSigma(double sigma) { m_Sigma = sigma; }
  void SetClamp(bool doClamp) { m_Clamp = doClamp; }
  void SetMaxValue(double maxValue) { m_MaxValue = maxValue; }
  void SetMinValue(double minValue) { m_MinValue = minValue; }

  virtual int Execute() = 0;

  protected:
  std::string m_InputImageFileName;
  std::string m_OutputImageFileName;
  bool m_NormalizeAcrossScale;
  double m_Sigma;
  bool m_Clamp;
  double m_MinValue;
  double m_MaxValue;

  private:
  // purposely not implemented
  GaussianSmoothBase(const Self &);
  void operator=(const Self &);
};

template <const unsigned int TDimension, class TPixelType >
class ITK_EXPORT GaussianSmooth : public GaussianSmoothBase 
{
public:
  /** Standard class typedefs. */
  typedef   GaussianSmooth Self;
  typedef   GaussianSmoothBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::OrientedImage< TPixelType, TDimension >  ImageTypeIn;
  typedef   itk::OrientedImage< TPixelType, TDimension >  ImageTypeOut;
  typedef   itk::ImageFileReader<ImageTypeIn> ImageReaderType;
  typedef   itk::ImageFileWriter<ImageTypeOut> ImageWriterType;
  typedef   itk::SmoothingRecursiveGaussianImageFilter< ImageTypeIn, 
                      ImageTypeOut > FilterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  GaussianSmooth() { };
  ~GaussianSmooth() { };

  int Execute();

private:
  // purposely not implemented
  GaussianSmooth(const Self &);
  void operator=(const Self &);

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlGaussianSmooth.txx"
#endif

} // end namespace crl

#endif

