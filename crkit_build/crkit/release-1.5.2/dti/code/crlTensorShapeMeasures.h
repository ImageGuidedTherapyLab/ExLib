
#ifndef _CRL_TENSORSHAPEMEASURES_H_INCLUDED 
#define _CRL_TENSORSHAPEMEASURES_H_INCLUDED 1

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSymmetricSecondRankTensor.h>
#include <itkDiffusionTensor3D.h>
#include <itkMatrix.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkSmartPointer.h>
#include <itkLightObject.h>

namespace crl {

class ITK_EXPORT TensorShapeMeasuresBase : public itk::LightObject
{
  // Defines the public API for this class.
  // but has no templated arguments.  These are provided by classes that
  // inherit this API.
  //  This allows us to create a pointer to a class of this type, 
  // and then change it at run time to match the types we actually need to 
  // deal with.

  public:
    typedef TensorShapeMeasuresBase Self;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(TensorShapeMeasuresBase, LightObject);

    TensorShapeMeasuresBase() {
      m_InputFileName = "";
      m_OutputCsFileName = "";
      m_OutputCpFileName = "";
      m_OutputClFileName = "";
    }

    virtual ~TensorShapeMeasuresBase() {
    }

    void SetInputFileName( const char *s ) { m_InputFileName = s; }
    void SetOutputCsFileName( const char *s ) { m_OutputCsFileName = s; }
    void SetOutputCpFileName( const char *s ) { m_OutputCpFileName = s; }
    void SetOutputClFileName( const char *s ) { m_OutputClFileName = s; }

    virtual int Execute() = 0;

  private:
    TensorShapeMeasuresBase(const TensorShapeMeasuresBase &); // purposely not implemented
    void operator=(const TensorShapeMeasuresBase &); //purposely not implemented

  protected:
    std::string m_OutputCsFileName;
    std::string m_OutputCpFileName;
    std::string m_OutputClFileName;
    std::string m_InputFileName;

};

template <const unsigned int TImageDimension, class TComponentType >
class ITK_EXPORT TensorShapeMeasures : 
                    public crl::TensorShapeMeasuresBase
{
  public: 
  /** Standard class typedefs. */
  typedef TensorShapeMeasures Self;
  typedef TensorShapeMeasuresBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::SymmetricSecondRankTensor<TComponentType, TImageDimension>  PixelType;
  typedef   itk::Image< PixelType, TImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  typedef   itk::ImageFileWriter<ImageType> ImageWriterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  TensorShapeMeasures() {
  };

  ~TensorShapeMeasures() {};

  int Execute();
  void EstimateMeasures(PixelType &in, TComponentType &Cs, 
                        TComponentType &Cp, TComponentType &Cl);

  private:
    TensorShapeMeasures(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

};

} // end namespace crl

#endif

