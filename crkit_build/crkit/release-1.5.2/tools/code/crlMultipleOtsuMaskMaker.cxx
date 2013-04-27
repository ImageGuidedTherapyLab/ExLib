/*
 * Copyright 2007 Simon Warfield simon.warfield@childrens.harvard.edu
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/



#include <tclap/CmdLine.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkOtsuMultipleThresholdsImageFilter.h>

#include "configuration.h"

typedef float PixelType;
typedef itk::Image<PixelType,3> ImageType;

int main(int argc, char* argv[])
{ 
  TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ', CRKIT_VERSION_STRING );
  TCLAP::UnlabeledValueArg<std::string> inputImageArg("inputImage","pathname to an ITK format input image",true,"","input image",cmd);
  TCLAP::ValueArg<std::string> outputImageArg("o","out","Output Image Argument",true,"","output image",cmd);
  TCLAP::ValueArg<unsigned int> binsArg("b","bins","number of bins in histogram", false, 256, "number of bins",cmd);
  TCLAP::UnlabeledValueArg<unsigned int> nthreshArg("numThreshes","number of desired thresholds",true,0,"num threshes",cmd);
  try
    {
    cmd.parse(argc,argv);
    }
  catch (TCLAP::ArgException& e)
    {
    std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
    exit(1);
    }

  std::string& infilen = inputImageArg.getValue();
  std::string& outfilen = outputImageArg.getValue();
  unsigned int bins = binsArg.getValue();
  unsigned int nthreshes = nthreshArg.getValue();
 
  if (nthreshes <= 0)
    {
    std::cerr << "number of threshes must be greater than zero" << std::endl;
    return 1;
    }
  
  std::cout << "input file=" << infilen << std::endl;
  std::cout << "output file=" << outfilen << std::endl;
  std::cout << "bins=" << bins << std::endl;
  std::cout << "nthreshes=" << nthreshes << std::endl;
  
  // read image
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer reader = ImageReaderType::New();
  
  try
    {
    reader->SetFileName(infilen);
    reader->Update();
    }
  catch (itk::ExceptionObject& err)
    {
    std::cerr << "error reading input file " << infilen << ": " << err << std::endl;
    return 1;
    }
  
  // set up otsu
  typedef itk::OtsuMultipleThresholdsImageFilter<ImageType,ImageType> OtsuType;
  OtsuType::Pointer otsu = OtsuType::New();
  otsu->SetInput(reader->GetOutput());
  otsu->SetNumberOfHistogramBins(bins);
  otsu->SetNumberOfThresholds(nthreshes);
  

  
  // write answer
  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  ImageWriterType::Pointer writer = ImageWriterType::New();
  
  try
    {
    writer->SetInput( otsu->GetOutput() );
    writer->SetFileName(outfilen);
    writer->Update();	
    }
  catch( itk::ExceptionObject& err )
    {
    std::cerr << "Failed either running Otsu or writing image " << outfilen << ": " << err << std::endl;
    return 1;
    }  
  
  return 0;
}
