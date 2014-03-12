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


#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>

#include <iostream>

int main(int argc, char *argv[])
{

  static int const ImageDimensionIn = 3;
  static int const ImageDimensionOut = 3;
  typedef   unsigned long InPixelType;
  typedef   unsigned long OutPixelType;

  typedef   itk::Image<InPixelType,ImageDimensionIn>  ImageTypeIn;
  typedef   itk::Image<OutPixelType,ImageDimensionOut>  ImageTypeOut;
  typedef   itk::ImageFileReader< ImageTypeIn >    ImageReaderType;
  typedef   itk::ImageFileWriter< ImageTypeOut >    ImageWriterType;
  typedef   itk::ConnectedComponentImageFilter< ImageTypeIn, ImageTypeOut> 
                                                    CCFilterType;
  typedef   itk::RelabelComponentImageFilter< ImageTypeIn, ImageTypeOut> 
                                                    RCFilterType;

  if (argc != 4) {
    std::cout << 
      "Usage: crlConnectedComponentFilter inimage outimage minsize" 
              << std::endl;
    exit(1);
  }

  unsigned long int minsize = atol(argv[3]);

  ImageReaderType::Pointer r;
  ImageWriterType::Pointer w;
  CCFilterType::Pointer ccf;
  RCFilterType::Pointer rcf;

  r = ImageReaderType::New();
  r->SetFileName(argv[1]);
  r->Update();
  w = ImageWriterType::New();
  w->SetFileName(argv[2]);

  ccf = CCFilterType::New();
  ccf->SetInput( r->GetOutput() );
  ccf->Update();

  rcf = RCFilterType::New();
  rcf->SetInput(ccf->GetOutput());
  rcf->SetMinimumObjectSize (minsize);

  w->SetInput(rcf->GetOutput());

  try {
    w->Update();
    w->UseCompressionOn();
    w->Write();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  exit(0); // success
}

