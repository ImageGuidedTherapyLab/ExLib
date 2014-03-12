
#ifndef _CRL_BINARYMORPH_INCLUDED
#define _CRL_BINARYMORPH_INCLUDED 1

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <iostream>

#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
// #include <itkBinaryMorphologicalClosingImageFilter.h>
// #include <itkBinaryMorphologicalOpeningImageFilter.h>
#include <itkBinaryBallStructuringElement.h>

class ITK_EXPORT crlBinaryMorphologyBase : public itk::LightObject
{
  public: 
  typedef crlBinaryMorphologyBase Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(crlBinaryMorphologyBase, LightObject);

  crlBinaryMorphologyBase() {};
  ~crlBinaryMorphologyBase() {};

  void SetInputFileName( std::string name) { m_InputFileName  = name;  };
  void SetOutputFileName(std::string name) { m_OutputFileName = name; };

  void SetOperationToDilation() { m_Operation = 1; /* Dilation */ };
  void SetOperationToErosion()  { m_Operation = 2; /* Erosion  */ };
  void SetOperationToClosing()  { m_Operation = 3; /* Closing  */ };
  void SetOperationToOpening()  { m_Operation = 4; /* Opening  */ };

  virtual void SetRadius(unsigned int radius) = 0;
  virtual int Execute() = 0;

  void SetForegroundValue(double value) {
    m_ForegroundValue = value;
  };

protected:
  std::string m_InputFileName;
  std::string m_OutputFileName;
  unsigned int m_Operation;     // Stores the type of operation to apply.
  double m_ForegroundValue;

private:
  crlBinaryMorphologyBase(const crlBinaryMorphologyBase &); // purposely not implemented
  void operator=(const crlBinaryMorphologyBase &); // purposely not implemented

};

template <unsigned int Dimension, class PixelType>
class ITK_EXPORT crlBinaryMorphology : public crlBinaryMorphologyBase
{
public:
  /** Standard class typedefs. */
  typedef crlBinaryMorphology Self;
  typedef crlBinaryMorphologyBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::Image< PixelType, Dimension> ImageType;
  typedef   itk::ImageFileWriter< ImageType >    ImageWriterType;
  typedef   itk::ImageFileReader< ImageType >    ImageReaderType;

  typedef itk::BinaryBallStructuringElement<
               PixelType, Dimension > StructuringElementType;

  typedef itk::BinaryErodeImageFilter< ImageType, ImageType,
               StructuringElementType > ErodeFilterType;

  typedef itk::BinaryDilateImageFilter< ImageType, ImageType,
               StructuringElementType > DilateFilterType;

  /* Future expansion:
  typedef itk::BinaryMorphologicalClosingImageFilter< ImageType, ImageType,
               StructuringElementType > ClosingFilterType;

  typedef itk::BinaryMorphologicalOpeningImageFilter< ImageType, ImageType,
               StructuringElementType > OpeningFilterType;
  */

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  crlBinaryMorphology()  { 
    structuringElement.SetRadius( 1 ); // ball of dimension 3
    structuringElement.CreateStructuringElement(); // ball of dimension 3
  };

  ~crlBinaryMorphology() { m_Operation = 0; };

  void SetRadius(unsigned int radius);

  int Execute();

private:
  crlBinaryMorphology(const Self&); // purposely not implemented
  void operator=(const Self&); //purposely not implemented

  StructuringElementType structuringElement;

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlBinaryMorphology.txx"
#endif

#endif

