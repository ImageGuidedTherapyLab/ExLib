
#ifndef _CRL_CONVERTDICOMMOSAIC_H_
#define _CRL_CONVERTDICOMMOSAIC_H_ 1

#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include <itkImage.h>
#include <itkImageSeriesReader.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <ctype.h>

struct DWIData
{
  std::string filename;
  unsigned int b_value;
  std::string g_directions;
  std::string series_number;
  std::string sizes;
  std::string space_directions;
  std::string space_origin;
};

class ITK_EXPORT crlConvertDICOMMOSAICBase : public itk::LightObject
{
  // Defines the public API for running
  //   crlConvertDICOMMOSAIC
  // but has no templated arguments.  The templated arguments are provided by
  // classes that inherit this API.
  //  This allows us to create a pointer to a class of this type, and
  // then to change it at run time to match the types we actually need to
  // deal with.

public:
  typedef crlConvertDICOMMOSAICBase Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(crlConvertDICOMMOSAICBase, LightObject);

  crlConvertDICOMMOSAICBase() {
  };

  ~crlConvertDICOMMOSAICBase() {
  };

  virtual int Execute() = 0;

  void SetInputDICOMDirectoryName( std::string dicomDir ) {
    m_dicomDir = dicomDir;
  };
  void SetOutputPrefix( std::string outputPrefix) {
    m_outputPrefix = outputPrefix;
  };
  void SetOutputSuffix( std::string outputSuffix) {
    m_outputSuffix = outputSuffix;
  };
  void SetSeriesID( std::string seriesID) {
    m_SeriesID = seriesID;
  };
  void AddRestrictionTag( std::string restrictionTag) {
    m_RestrictionTag = restrictionTag;
  };
  void AddRestrictionTags(std::vector<std::string> &restrictionTags) {
    m_RestrictionTags = restrictionTags;
  }
  void SetIgnoreDuplicateData( bool ignoreDuplicateData ) {
    m_IgnoreDuplicateData = ignoreDuplicateData;
  };
  void SetProcessNonMOSAICData(bool processNonMOSAICData) {
    m_processNonMOSAICData = processNonMOSAICData;
  };
  void SetRecursiveProcessing(bool recursiveProcessing) {
    m_recursiveProcessing = recursiveProcessing;
  };
  void SetFMRIProcessing(bool fMRIProcessing) {
    m_fMRIProcessing = fMRIProcessing;
  };

  void SetMOSAICOrientation(std::string mosaicOrientation) {
    m_mosaicOrientation = mosaicOrientation;
  };

private:
  crlConvertDICOMMOSAICBase(const crlConvertDICOMMOSAICBase &); // purposely not implemented
  void operator=(const crlConvertDICOMMOSAICBase &); // purposely not implemented.

protected:
  bool m_IgnoreDuplicateData;
  std::string m_SeriesID;
  std::string m_outputSuffix;
  std::string m_outputPrefix;
  std::string m_dicomDir;
  std::string m_RestrictionTag;
  std::vector<std::string> m_RestrictionTags;
  bool m_processNonMOSAICData;
  bool m_recursiveProcessing;
  bool m_fMRIProcessing;
  std::string m_mosaicOrientation;
  DWIData dwiData;

public:
  bool dwiImage;
  std::vector<DWIData> dwiDataVector;

};

template <const unsigned int TImageDimension, class TComponentType >
class ITK_EXPORT crlConvertDICOMMOSAIC :
  public crlConvertDICOMMOSAICBase
{
public:
  /** Standard class typedefs. */
  typedef   crlConvertDICOMMOSAIC Self;
  typedef   crlConvertDICOMMOSAICBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::Image< TComponentType, TImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  typedef   itk::ImageFileWriter<ImageType> ImageWriterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  crlConvertDICOMMOSAIC() {
  };

  ~crlConvertDICOMMOSAIC() {};

  int Execute();

private:
  crlConvertDICOMMOSAIC(const Self&); // purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

#endif
