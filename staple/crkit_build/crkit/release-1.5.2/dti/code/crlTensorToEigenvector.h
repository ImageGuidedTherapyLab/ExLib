
#ifndef _CRL_TENSORTOEIGENVECTOR_H_INCLUDED 
#define _CRL_TENSORTOEIGENVECTOR_H_INCLUDED 1

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkDiffusionTensor3D.h>
#include <itkVector.h>
#include <itkMatrix.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkSmartPointer.h>
#include <itkLightObject.h>

namespace crl {

class ITK_EXPORT TensorToEigenvectorBase : public itk::LightObject
{
  // Defines the public API for this class.
  // but has no templated arguments.  These are provided by classes that
  // inherit this API.
  //  This allows us to create a pointer to a class of this type, 
  // and then change it at run time to match the types we actually need to 
  // deal with.

  public:
    typedef TensorToEigenvectorBase Self;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(TensorToEigenvectorBase, LightObject);

    TensorToEigenvectorBase() {
      m_InputFileName = "";
      m_OutputFileName = "";
      m_EigenvectorIndex = 2;
    }

    virtual ~TensorToEigenvectorBase() {
    }

    void SetInputFileName( const char *s ) { m_InputFileName = s; }
    void SetOutputFileName( const char *s ) { m_OutputFileName = s; }

    virtual int Execute() = 0;

  private:
    TensorToEigenvectorBase(const TensorToEigenvectorBase &); // purposely not implemented
    void operator=(const TensorToEigenvectorBase &); //purposely not implemented

  protected:
    std::string m_OutputFileName;
    std::string m_InputFileName;
    // 2 is eigenvector with largest eigenvalue, 1 is middle, 0 is smallest.
    unsigned int m_EigenvectorIndex; 

};

template <const unsigned int TImageDimension, class TComponentType >
class ITK_EXPORT TensorToEigenvector : 
                    public crl::TensorToEigenvectorBase
{
  public: 
  /** Standard class typedefs. */
  typedef TensorToEigenvector Self;
  typedef TensorToEigenvectorBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::DiffusionTensor3D<TComponentType> PixelType;
  typedef   itk::Image< PixelType, TImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;

  typedef   itk::Vector<float, 3> VectorPixelType;
  typedef   itk::Image< VectorPixelType, TImageDimension> VectorImageType;
  typedef   itk::ImageFileWriter<VectorImageType> ImageWriterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  TensorToEigenvector() {
  };

  ~TensorToEigenvector() {};

  int Execute();
  void EstimateEigenvector(PixelType &in, VectorPixelType &v);

  private:
    TensorToEigenvector(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

};

} // end namespace crl

#endif

