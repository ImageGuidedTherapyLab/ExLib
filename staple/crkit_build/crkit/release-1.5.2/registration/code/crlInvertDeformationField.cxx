#include <iostream>
#include <fstream>
#include <itkNumericTraits.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkFixedPointInverseDeformationFieldImageFilter.h"




int main(int argc, char **argv)
{
  if( argc < 3 )
    {
      std::cerr << "usage: " << argv[0] << " inputImage outputImage" << std::endl;
    return 1;
    }

  const unsigned int Dimension = 3;

  typedef itk::Vector< float, Dimension >    VectorPixelType;
  typedef  itk::Image< VectorPixelType, Dimension > InputDFType;
  typedef  itk::Image< VectorPixelType, Dimension > OutputDFType;

  typedef itk::ImageFileReader< InputDFType >      ReaderType;
  typedef itk::ImageFileWriter< OutputDFType >     WriterType;

  typedef itk::FixedPointInverseDeformationFieldImageFilter<InputDFType, OutputDFType> FPInverseType;

  // read the file
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[1]);

  try
    {
      reader->UpdateLargestPossibleRegion();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception detected while reading " << argv[1];
    std::cerr << " : "  << e.GetDescription();
    return 1;
    }

  // invert the deformationfield
  InputDFType::Pointer inputDf = reader->GetOutput();

  FPInverseType::Pointer inverter = FPInverseType::New();
  inverter->SetInput(inputDf);
  inverter->SetOutputOrigin(inputDf->GetOrigin());
  inverter->SetSize(inputDf->GetLargestPossibleRegion().GetSize());
  inverter->SetOutputSpacing(inputDf->GetSpacing());
  inverter->SetNumberOfIterations(20);
  inverter->Update();
  OutputDFType::Pointer outputDf = inverter->GetOutput();


  // write the file
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(inverter->GetOutput());
  try
    {
    writer->SetFileName( argv[2] );
    writer->Update();
    }
  catch (...)
    {
    std::cerr << "Error during write of " << argv[2] << std::endl;
    return 1;
    }

  return 0;
}
