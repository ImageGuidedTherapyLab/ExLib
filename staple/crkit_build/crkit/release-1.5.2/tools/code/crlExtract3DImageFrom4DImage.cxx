/*
 * Copyright 2008 Simon Warfield simon.warfield@childrens.harvard.edu
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 * http://crl.med.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>

#include <iostream>

#include "tclap/CmdLine.h"
#include "configuration.h"

#ifdef WIN32
#define snprintf _snprintf
#define PATH_MAX 65536
#endif


int main(int argc, char *argv[])
{

  static int const ImageDimensionIn = 4;
  static int const ImageDimensionOut = 3;
  typedef   float PixelType;

  typedef   itk::Image<PixelType,ImageDimensionIn>  ImageTypeIn;
  typedef   itk::Image<PixelType,ImageDimensionOut>  ImageTypeOut;
  typedef   itk::ImageFileReader< ImageTypeIn >    ImageReaderType;
  typedef   itk::ImageFileWriter< ImageTypeOut >    ImageWriterType;
  typedef   itk::ExtractImageFilter< ImageTypeIn, ImageTypeOut> 
                                                    FilterType;

  int lowIndexComponent = -1;
  int highIndexComponent = -1;
  int axisToExtractFrom = ImageDimensionIn - 1;
  std::string *inputImage = new std::string("");
  std::string *outputImage = new std::string("");
  int exitStatus = EXIT_SUCCESS;

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',
              CRKIT_VERSION_STRING);
      TCLAP::ValueArg<std::string> outputImageArg("o", "outputImage",
        "Output Image File Name", 
        true,"",
        "Output image file name", cmd);

      TCLAP::ValueArg<std::string> inputImageArg("i", "inputImage",
          "Input Image File Name",
          true,"",
          "Input image file name", cmd);

      TCLAP::ValueArg<int> highIndexComponentArg("m", "highIndexComponent",
          "Highest Component Index",
           true,-1,
           "integer value of highest component index",cmd);

      TCLAP::ValueArg<int> lowIndexComponentArg("l", "lowIndexComponent",
           "Lowest Component Index",
           true,-1,
           "integer value of lowest component index", cmd);
      TCLAP::ValueArg<int> axisToExtractFromArg("a","axis",
         "Axis to extract from",
          true,0,
          "integer value of axis to extract from", cmd);

      cmd.parse(argc,argv);
      if (inputImageArg.isSet()) inputImage = new std::string(inputImageArg.getValue());
      if (outputImageArg.isSet()) outputImage = new std::string(outputImageArg.getValue());
      if (lowIndexComponentArg.isSet()) lowIndexComponent = lowIndexComponentArg.getValue();
      if (highIndexComponentArg.isSet()) highIndexComponent = highIndexComponentArg.getValue();
      if (axisToExtractFromArg.isSet()) {
        axisToExtractFrom = axisToExtractFromArg.getValue();
      };
    }
  catch (TCLAP::ArgException &e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(1);
    }

  if ( (lowIndexComponent == -1) && (highIndexComponent == -1) ) {
    std::cout << "Need a component index or range from zero up." << std::endl;
    exit(1);
  }
  if ( (lowIndexComponent == -1) && (highIndexComponent != -1) ||
       (lowIndexComponent != -1) && (highIndexComponent == -1) ) {
    if (lowIndexComponent == -1) lowIndexComponent = highIndexComponent;
    if (highIndexComponent == -1) highIndexComponent = lowIndexComponent;
  }

  if (lowIndexComponent > highIndexComponent) {
    int tmp = highIndexComponent;
    highIndexComponent = lowIndexComponent;
    lowIndexComponent = tmp;
  }

  ImageReaderType::Pointer r;
  ImageWriterType::Pointer w;
  FilterType::Pointer f;

  r = ImageReaderType::New();
  r->SetFileName(inputImage->c_str());
  r->Update();
  w = ImageWriterType::New();
  w->SetFileName(outputImage->c_str());

  f = FilterType::New();
  f->SetInput( r->GetOutput() );
  ImageTypeIn::RegionType inputRegion = 
                   r->GetOutput()->GetLargestPossibleRegion();
  ImageTypeIn::SizeType size = inputRegion.GetSize();
  ImageTypeIn::IndexType start = inputRegion.GetIndex();

  /* Each tool should operate quietly in the case of correct function.
    for (int i = 0; i < ImageDimensionIn; i++) {
    std::cout << "size["<< i<<"] " << size[i] << std::endl;
    std::cout << "start["<< i<<"] " << start[i] << std::endl;
    }
  */

  size[axisToExtractFrom] = 0;
  int bufferlen = 1 + static_cast<int>(vcl_ceil(1 + log10(fabs((float)highIndexComponent))));
  bufferlen += strlen(outputImage->c_str());
  char *buffer = new char[bufferlen];
  for (int i = lowIndexComponent; i < highIndexComponent; i++) {
    memset(buffer, '\0', bufferlen);
    snprintf(buffer, bufferlen, outputImage->c_str(), i);
    w->SetFileName(buffer);

    start[axisToExtractFrom] = i;

    ImageTypeIn::RegionType desiredRegion;
    desiredRegion.SetSize( size );
    desiredRegion.SetIndex( start );

    f->SetExtractionRegion( desiredRegion );
    f->Update();
    w->SetInput(f->GetOutput());

    try {
      w->Update();
      w->UseCompressionOn();
      w->Write();
    } catch (itk::ExceptionObject &e) {
        std::cerr << "Caught ITK exception: " << e << std::endl;
        exitStatus = EXIT_FAILURE;
    }
  }
  delete [] buffer;
  return exitStatus;
}

