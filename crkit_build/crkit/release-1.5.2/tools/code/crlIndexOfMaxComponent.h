
#ifndef _CRL_INDEXOFMAXCOMPONENT_INCLUDED
#define _CRL_INDEXOFMAXCOMPONENT_INCLUDED 1

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>

#include <iostream>

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

class ITK_EXPORT IndexOfMaxComponentBase : public itk::LightObject
{
  public:
    typedef IndexOfMaxComponentBase Self;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(IndexOfMaxComponentBase, LightObject);

    IndexOfMaxComponentBase() {};
    ~IndexOfMaxComponentBase() {};

    void SetInputFileName(std::string name) {
      m_InputFileName = name;
    };

    void SetOutputFileName(std::string name) {
      m_OutputFileName = name;
    };

    virtual int Execute() = 0;

  protected:
    std::string m_InputFileName;
    std::string m_OutputFileName;

  private:
    // purposely not implemented
    IndexOfMaxComponentBase(const IndexOfMaxComponentBase &); 
    void operator=(const IndexOfMaxComponentBase &); 

};

template <unsigned int Dimension, class PixelType>
class ITK_EXPORT IndexOfMaxComponent : public IndexOfMaxComponentBase
{
  public:
    /** Standard class typedefs. */
    typedef IndexOfMaxComponent Self;
    typedef IndexOfMaxComponentBase Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    static int const OutImageDimension = Dimension - 1;

    typedef   itk::Image<PixelType, Dimension>  InImageType;
    typedef   itk::Image<PixelType, OutImageDimension>  OutImageType;
    typedef   itk::ImageFileReader< InImageType >    ImageReaderType;
    typedef   itk::ImageFileWriter< OutImageType >    ImageWriterType;

    typedef typename OutImageType::RegionType OutRegionType;
    typedef typename OutImageType::SizeType   OutSizeType;
    typedef typename OutImageType::IndexType  OutIndexType;
    typedef typename OutImageType::Pointer    OutImagePointer;
    typedef typename itk::ImageRegionIteratorWithIndex< OutImageType >
                     OutRegionIteratorType;

    typedef typename InImageType::RegionType InRegionType;
    typedef typename InImageType::SizeType   InSizeType;
    typedef typename InImageType::IndexType  InIndexType;
    typedef typename InImageType::Pointer    InImagePointer;
    typedef typename itk::ImageRegionIteratorWithIndex< InImageType >
                     InRegionIteratorType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    IndexOfMaxComponent()  { };
    ~IndexOfMaxComponent() { };

    int Execute();

  private:
    IndexOfMaxComponent(const Self&); // purposely not implemented
    void operator=(const Self&); //purposely not implemented

};

}  // end namespace crl

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlIndexOfMaxComponent.txx"
#endif

#endif

