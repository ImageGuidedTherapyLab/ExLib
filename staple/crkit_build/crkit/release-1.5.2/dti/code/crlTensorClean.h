
#ifndef _TENSORCHECK_H_INCLUDED
#define _TENSORCHECK_H_INCLUDED 1

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSymmetricSecondRankTensor.h>
#include <itkDiffusionTensor3D.h>
#include <itkMatrix.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkSmartPointer.h>
#include <itkLightObject.h>

class ITK_EXPORT TensorCheckBase : public itk::LightObject
{
  // Defines the public API for running TensorCheck
  // but has no templated arguments.  These are provided by classes that
  // inherit this API.
  //  This allows us to create a pointer to a class of this type, 
  // and then change it at run time to match the types we actually need to 
  // deal with.

  public:
    typedef TensorCheckBase Self;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(TensorCheckBase, LightObject);

    TensorCheckBase() {
      m_InputFileName = "";
      m_OutputFileName = "";
      m_CountTensorsChanged = 0;
    }

    virtual ~TensorCheckBase() {
    }

    void SetInputFileName( const char *s ) { m_InputFileName = s; }
    void SetOutputFileName( const char *s ) { m_OutputFileName = s; }
    void SetUseCompression( bool useCompression ) { 
                                  m_UseCompression = useCompression; }
    void SetCountTensorsChanged( int s) { m_CountTensorsChanged = s; }
    int GetCountTensorsChanged( ) { return m_CountTensorsChanged; }
    void IncrementCountTensorsChanged( ) { m_CountTensorsChanged++; }

    virtual int Execute() = 0;

  private:
    TensorCheckBase(const TensorCheckBase &); // purposely not implemented
    void operator=(const TensorCheckBase &); //purposely not implemented

  protected:
    std::string m_OutputFileName;
    std::string m_InputFileName;
    int m_CountTensorsChanged;
    bool m_UseCompression;

};

template <const unsigned int TImageDimension, class TComponentType >
class ITK_EXPORT TensorCheck : 
                    public TensorCheckBase
{
  public: 
  /** Standard class typedefs. */
  typedef TensorCheck Self;
  typedef TensorCheckBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::SymmetricSecondRankTensor<TComponentType, TImageDimension>  PixelType;
  typedef   itk::Image< PixelType, TImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  typedef   itk::ImageFileWriter<ImageType> ImageWriterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  TensorCheck() {
  };

  ~TensorCheck() {};

  int Execute();
  bool CorrectTensor(PixelType &in, PixelType &out);

  private:
    TensorCheck(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

};

#endif

