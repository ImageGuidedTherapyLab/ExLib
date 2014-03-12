
#ifndef _CRL_VECTOR3TORGB_H_INCLUDED 
#define _CRL_VECTOR3TORGB_H_INCLUDED 1

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkRGBPixel.h>
#include <itkVector.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkSmartPointer.h>
#include <itkLightObject.h>

namespace crl {

class ITK_EXPORT Vector3ToRGBBase : public itk::LightObject
{
  // Defines the public API for this class.
  // but has no templated arguments.  These are provided by classes that
  // inherit this API.
  //  This allows us to create a pointer to a class of this type, 
  // and then change it at run time to match the types we actually need to 
  // deal with.

  public:
    typedef Vector3ToRGBBase Self;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(Vector3ToRGBBase, LightObject);

    Vector3ToRGBBase() {
      m_InputFileName = "";
      m_OutputFileName = "";
      m_ScaleFactor = 1.0;
    }

    virtual ~Vector3ToRGBBase() {
    }

    void SetInputFileName( const char *s ) { m_InputFileName = s; }
    void SetOutputFileName( const char *s ) { m_OutputFileName = s; }

    virtual int Execute() = 0;

  private:
    Vector3ToRGBBase(const Vector3ToRGBBase &); // purposely not implemented
    void operator=(const Vector3ToRGBBase &); //purposely not implemented

  protected:
    std::string m_OutputFileName;
    std::string m_InputFileName;
    double m_ScaleFactor;

};

template <const unsigned int TImageDimension, class TComponentType >
class ITK_EXPORT Vector3ToRGB : 
                    public crl::Vector3ToRGBBase
{
  public: 
  /** Standard class typedefs. */
  typedef Vector3ToRGB Self;
  typedef Vector3ToRGBBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::Vector<TComponentType, 3> PixelType;
  typedef   itk::Image< PixelType, TImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  typedef   itk::RGBPixel< unsigned char > ColorPixelType;
  typedef   itk::Image< ColorPixelType, TImageDimension> ColorImageType;
  typedef   itk::ImageFileWriter<ColorImageType> ImageWriterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  Vector3ToRGB() {
  };

  ~Vector3ToRGB() {};

  int Execute();
  void EstimateRGB(PixelType &in, ColorPixelType &color);

  private:
    Vector3ToRGB(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

};

} // end namespace crl

#endif

