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

typedef float PixelType;	// using same scalar pixel type for mask image and gray image
typedef itk::Image<PixelType,3> ImageType;

int main(int argc, char* argv[])
{ 
  TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ', CRKIT_VERSION_STRING );
  TCLAP::UnlabeledValueArg<std::string> inputImageArg("inputImage","pathname to an ITK format input image",true,"","input image",cmd);
  TCLAP::UnlabeledValueArg<std::string> maskImageArg("maskImage","pathname to an ITK format segmentation or mask",true,"","mask image",cmd);
  TCLAP::ValueArg<std::string> outputImageArg("o","out","Output Image Argument",true,"","output image",cmd);
  TCLAP::UnlabeledMultiArg<unsigned int > tissueClassArg("classes","values for each class to select",true,"classes",cmd);
  TCLAP::ValueArg<unsigned int> backArg("b","background","value to assign to background", false, 0, "background",cmd);
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
  std::string& maskfilen = maskImageArg.getValue();
  std::string& outfilen = outputImageArg.getValue();
  PixelType background = backArg.getValue();
  const std::vector<unsigned int>& classes = tissueClassArg.getValue();
  
  std::cout << "input file=" << infilen << std::endl;
  std::cout << "mask file=" << maskfilen << std::endl;
  std::cout << "output file=" << outfilen << std::endl;
  std::cout << "background value=" << background << std::endl;
  std::cout << "classes to select=";
  for ( unsigned int i = 0; i < classes.size(); i++)
    {
    std::cout << " " << classes[i];
    }
  std::cout << std::endl;
  
  // read images
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer reader = ImageReaderType::New();
  ImageReaderType::Pointer maskreader = ImageReaderType::New();
  
  try
    {
    reader->SetFileName(infilen);
    maskreader->SetFileName(maskfilen);
    reader->Update();
    maskreader->Update();
    }
  catch (itk::ExceptionObject& err)
    {
    std::cerr << "error reading input file: " << err << std::endl;
    return 1;
    }
  
  ImageType::Pointer input = reader->GetOutput();
  input->DisconnectPipeline();	//strictly necessary here because we'll modify the image
  ImageType::Pointer mask = maskreader->GetOutput();
  mask->DisconnectPipeline();
  
  typedef itk::ImageRegionIteratorWithIndex<ImageType> IterType;
  const ImageType::SizeType& t1 = input->GetLargestPossibleRegion().GetSize();
  const ImageType::SizeType& t2 = mask->GetLargestPossibleRegion().GetSize();
  
  if ( t1 != t2 )
    {
    std::cerr << "mask size: " << t2 << " is not the same as image size " << t1 << std::endl;
    std::cerr << "was this mask from this image?" << std::endl;
    return 1;
    }
  
  IterType iter( input, input->GetLargestPossibleRegion() );
  PixelType maskval;
  PixelType oldval;
  for ( iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
    {
    maskval = mask->GetPixel( iter.GetIndex() );
    oldval = iter.Value();
    iter.Value() = background;	// start off assuming background
    for ( unsigned int i = 0; i < classes.size(); ++i )
      {
      if ( classes[i] == maskval )
	{
	iter.Value() = oldval;	// if found, restore the grey value
	}
      }
    }

  
  // write answer
  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  ImageWriterType::Pointer writer = ImageWriterType::New();
  
  try
    {
    writer->SetInput( input );
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
