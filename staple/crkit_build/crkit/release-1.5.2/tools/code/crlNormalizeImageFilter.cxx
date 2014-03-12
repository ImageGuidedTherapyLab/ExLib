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
#include <itkNormalizeImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include "tclap/CmdLine.h"

int main( int argc, char * argv[] )
{

  double desiredMean = 0.0;
  double desiredVariance = 1.0;
  std::string inputFileName;
  std::string outputFileName;

  try {
    // Define the command line object, and insert a message
    // that describes the program.
    TCLAP::CmdLine cmd("Computational Radiology Laboratory", ' ',
        "1.0");

    // Define a value argument for the input image on the command line.
    // This can be specified without any flag. It is simply taken from the
    // first unflagged position.
    TCLAP::UnlabeledValueArg<std::string> inputImageArg("input",
            "Input Image", true,
            "/dev/null", // default value
            "input image file name as string",
            cmd);

    // Define a value argument for the output image on the command line.
    // This can be specified without any flag. It is simply taken from the
    // next unflagged position.
    TCLAP::UnlabeledValueArg<std::string> outputImageArg("output",
            "Output Image", true,
            "/dev/null",
            "output image file name as string",
            cmd);

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg< double > meanArg("m", "mean",
            "Mean value image will have", false,
            0.0,
            "double precision floating point number",
            cmd);

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg< double > varianceArg("v", "variance",
            "Image will have specified variance", false,
            1.0, 
            "double precision floating point number",
            cmd);

   // Now parse the command line and assign variable values as required.
   cmd.parse( argc, argv );

   inputFileName = inputImageArg.getValue();
   outputFileName = outputImageArg.getValue();
   if (meanArg.isSet()) { desiredMean = meanArg.getValue(); }
   if (varianceArg.isSet()) { desiredVariance = varianceArg.getValue(); }

  } catch  (TCLAP::ArgException &e)
  {
    std::cerr << "error: " << e.error() << " for arg " <<
                  e.argId() << std::endl;
    std::cerr << "This program computes normalized images."
        << std::endl;
    return EXIT_FAILURE;
  }

  // When dealing with parametric images such as fractional anisotropy, it
  // is a good idea to preserve as much precision as possible.
  typedef   float PixelType;
  static int const Dimension = 3;

  typedef itk::Image< PixelType,  Dimension >   ImageType;

  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( inputFileName );
  writer->SetFileName( outputFileName );

  typedef itk::StatisticsImageFilter< ImageType >  StatisticsFilterType;
  StatisticsFilterType::Pointer statFilter = StatisticsFilterType::New();

  typedef itk::ShiftScaleImageFilter< ImageType, ImageType >  ShiftScaleFilterType;
  ShiftScaleFilterType::Pointer shiftScaleFilter1 = ShiftScaleFilterType::New();
  ShiftScaleFilterType::Pointer shiftScaleFilter2 = ShiftScaleFilterType::New();

  statFilter->SetInput( reader->GetOutput() );
  statFilter->Update();

  shiftScaleFilter1->SetInput( statFilter->GetOutput() );
  // Shift the distribution to have a mean of zero
  shiftScaleFilter1->SetShift( - statFilter->GetMean() );
  // And a unit variance 
  shiftScaleFilter1->SetScale( 1.0/sqrt( statFilter->GetVariance() ) );

  shiftScaleFilter2->SetInput( shiftScaleFilter1->GetOutput() );
  shiftScaleFilter2->SetShift( desiredMean/sqrt(desiredVariance) );
  shiftScaleFilter2->SetScale( sqrt( desiredVariance ) );

  try
    {
    writer->SetInput( shiftScaleFilter2->GetOutput() );
    writer->UseCompressionOn();
    writer->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "Failed writing output file." << std::endl;
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
