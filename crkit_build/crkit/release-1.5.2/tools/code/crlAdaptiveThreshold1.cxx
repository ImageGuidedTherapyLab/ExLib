/*
 * Copyright 2007-2009 Children's Hospital Boston
 * Contact: Simon Warfield simon.warfield@childrens.harvard.edu
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
#include <itkStatisticsImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkScalarImageToHistogramGenerator.h>
#include "configuration.h"

typedef float PixelType;
typedef itk::Image<PixelType,3> ImageType;

int main(int argc, char* argv[])
{ 
  TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ', CRKIT_VERSION_STRING );
  TCLAP::UnlabeledValueArg<std::string> inputImageArg("inputImage","pathname to an ITK format input image",true,"","input image",cmd);
  TCLAP::ValueArg<std::string> outputImageArg("o","out","Output Image Argument",true,"","output image",cmd);
  TCLAP::ValueArg<unsigned int> binsArg("b","bins","number of bins in histogram", false, 256, "number of bins",cmd);
  TCLAP::UnlabeledValueArg<float> loThreshArg("backThresh","fraction of the mean to threshold below and call background",true,0.0,"mean fraction",cmd);
  TCLAP::UnlabeledValueArg<float> cumSumArg("cumSumFraction","fraction of the total above-threshold pixels to include",true,0.0,"histogram fraction",cmd);
  TCLAP::SwitchArg clipArg("c","clip","if set will cause pixels outside the threshold to be set to zero", cmd, false);	
  TCLAP::SwitchArg semiArg("s","semiclip","if set all pixels at lower threshold will go to zero, upper will be clamped", cmd, false);
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
  float lothresh = loThreshArg.getValue();
  float cumsum = cumSumArg.getValue();
  unsigned int bins = binsArg.getValue();
  bool clip = clipArg.getValue();
  bool semi = semiArg.getValue();

  if ( clip && semi )
    {
    std::cerr << "must choose *either* clip or semi, but not both (-c and -s)" << std::endl;
    return 1;
    }
  
  std::cout << "input file=" << infilen << std::endl;
  std::cout << "output file=" << outfilen << std::endl;
  std::cout << "background fraction=" << lothresh << std::endl;
  std::cout << "cum sum fraction=" << cumsum << std::endl;
  std::cout << "bins=" << bins << std::endl;
  std::cout << "clip=" << clip << std::endl;
  
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
  
  
  // generate mean
  typedef itk::StatisticsImageFilter<ImageType> StatsFilter;
  StatsFilter::Pointer stats = StatsFilter::New();
  stats->SetInput(reader->GetOutput());
  stats->Update();
  const float mean = stats->GetMean();
  const float max = stats->GetMaximum();
  const float min = stats->GetMinimum();
  std::cout << "input mean is " << mean << std::endl;
  
  // threshold background
  float thresh = lothresh*mean;
  std::cout << "will threshold for background at " << lothresh << "*" << mean << "=" << thresh << std::endl;

  // histogram image (excluding background)
  typedef itk::Statistics::ScalarImageToHistogramGenerator<ImageType> Histogrammer;
  Histogrammer::Pointer histogrammer = Histogrammer::New();

  histogrammer->SetInput(reader->GetOutput() );

  histogrammer->SetHistogramMin(thresh);
  histogrammer->SetHistogramMax(max);
  histogrammer->SetNumberOfBins(bins);
  histogrammer->Compute();
  
  typedef Histogrammer::HistogramType HistogramType;
  const HistogramType* histogram = histogrammer->GetOutput();
  
  HistogramType::ConstIterator itr = histogram->Begin();
  HistogramType::ConstIterator end = histogram->End();
  
  double count = 0.0;
  for ( unsigned int i=0; itr != end; ++itr,++i)
    {
    count += static_cast<double>(itr.GetFrequency());
    std::cout << "histogram[" << itr.GetMeasurementVector() << "]=" << itr.GetFrequency() << std::endl;
//    std::cout << "mv[" << i << "]=" << itr.GetMeasurementVector() << std::endl;
    }
  std::cout << "TOTAL=" << count << std::endl;
  
  // find adaptive threshold and do it
  typedef itk::ThresholdImageFilter<ImageType> ThreshFilter;
  ThreshFilter::Pointer thresher = ThreshFilter::New();
  double hithresh = min;
  double cutoff = cumsum*count;
  
  std::cout << "cutoff is " << cutoff << std::endl;
  
  count = 0.0;
  itr=histogram->Begin();
  end=histogram->End();
  while ( count < cutoff && itr != end)
    {
    count += itr.GetFrequency();
    hithresh = itr.GetMeasurementVector()[0];
    ++itr;
    }
  
  std::cout << "thresholding outside of [" << thresh << "," << hithresh << "]" << std::endl;
  if (clip)
    {
    std::cout << "using clipping (setting to zero)" << std::endl;
    thresher->SetInput(reader->GetOutput());
    thresher->SetLower(thresh);
    thresher->SetUpper(hithresh);
    thresher->SetOutsideValue(0.0);
    }
  else
    {
    if (semi )
      {
      std::cout << "using semi-clamping (below the bottom->zero, above the top gets clamped)" << std::endl;
      }
    else
      {
      std::cout << "using clamping (setting to the value of the threshold boundary exceeded)" << std::endl;
      }

    thresher->SetInput(reader->GetOutput());
    thresher->SetOutsideValue((semi ? 0.0 : thresh));
    thresher->ThresholdBelow(thresh);
    thresher->Update();
    ImageType::Pointer temp = thresher->GetOutput();
    temp->DisconnectPipeline();
    thresher->SetInput(temp);
    thresher->SetOutsideValue(hithresh);
    thresher->ThresholdAbove(hithresh);
    thresher->Update();
    }
  
  // write answer
  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  ImageWriterType::Pointer writer = ImageWriterType::New();
  
  try
    {
    writer->SetInput( thresher->GetOutput() );
    writer->SetFileName(outfilen);
    writer->Update();	
    }
  catch( itk::ExceptionObject& err )
    {
    std::cerr << "Failed either thresholding or writing image " << outfilen << ": " << err << std::endl;
    return 1;
    }  
  
  return 0;
}
