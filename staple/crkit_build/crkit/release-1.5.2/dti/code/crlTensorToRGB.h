
#ifndef _CRL_TENSORTORGB_H_INCLUDED 
#define _CRL_TENSORTORGB_H_INCLUDED 1

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkDiffusionTensor3D.h>
#include <itkRGBPixel.h>
#include <itkMatrix.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkSmartPointer.h>
#include <itkLightObject.h>

namespace crl {

class ITK_EXPORT TensorToRGBBase : public itk::LightObject
{
  // Defines the public API for this class.
  // but has no templated arguments.  These are provided by classes that
  // inherit this API.
  //  This allows us to create a pointer to a class of this type, 
  // and then change it at run time to match the types we actually need to 
  // deal with.

  public:
    typedef TensorToRGBBase Self;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(TensorToRGBBase, LightObject);

    TensorToRGBBase() : m_MaxFractionalAnisotropy(1.0) {
      m_InputFileName = "";
      m_OutputFileName = "";
    }

    virtual ~TensorToRGBBase() {
    }

    void SetInputFileName( const char *s ) { m_InputFileName = s; }
    void SetOutputFileName( const char *s ) { m_OutputFileName = s; }

    virtual int Execute() = 0;

  private:
    TensorToRGBBase(const TensorToRGBBase &); // purposely not implemented
    void operator=(const TensorToRGBBase &); //purposely not implemented

  protected:
    std::string m_OutputFileName;
    std::string m_InputFileName;
    double m_MaxFractionalAnisotropy;

};

template <const unsigned int TImageDimension, class TComponentType >
class ITK_EXPORT TensorToRGB : 
                    public crl::TensorToRGBBase
{
  public: 
  /** Standard class typedefs. */
  typedef TensorToRGB Self;
  typedef TensorToRGBBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::DiffusionTensor3D<TComponentType> PixelType;
  typedef   itk::Image< PixelType, TImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  typedef   itk::RGBPixel< unsigned char > ColorPixelType;
  typedef   itk::Image< ColorPixelType, TImageDimension> ColorImageType;
  typedef   itk::ImageFileWriter<ColorImageType> ImageWriterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  TensorToRGB() {
  };

  ~TensorToRGB() {};

  int Execute();
  void EstimateRGB(PixelType &in, ColorPixelType &color);

  private:
    TensorToRGB(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

};

} // end namespace crl

#endif

