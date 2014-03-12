/*
 * Copyright 2009 Simon Warfield simon.warfield@childrens.harvard.edu
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

// The purpose of this program is to resample a volume to have isotropic voxels,
// and to be represented on a cubic lattice.

// IO mechanisms
#include <itkOrientedImage.h>
#include <itkOrientImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

// Interpolators supported
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkWindowedSincInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>

// TCLAP
#include <tclap/CmdLine.h>
#include "configuration.h"

// Macro to allow for programmatic generation of other dimensionality data.
#define DIMENSIONALITY  (3)

int main( int argc, char * argv[] )
{
  std::string *inputImageFile = new std::string();
  std::string *interpolationStyle = new std::string();
  std::string *outputImageFile = new std::string();

  bool AxisSpacingSet[DIMENSIONALITY];
  double AxisSpacingValue[DIMENSIONALITY];

  for (unsigned int i = 0; i < DIMENSIONALITY; i++) {
    AxisSpacingSet[i] = false;
    AxisSpacingValue[i] = 0.0;
  }

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',
                          CRKIT_VERSION_STRING  );
      TCLAP::ValueArg<double> zAxisSpacingArg("z","zAxisSpacing",
    "desired z axis spacing.", false, 0.0,
    "double floating point value of the desired z axis spacing.",cmd);
      TCLAP::ValueArg<double> yAxisSpacingArg("y","yAxisSpacing",
    "desired y axis spacing.", false, 0.0,
    "double floating point value of the desired y axis spacing.",cmd);
      TCLAP::ValueArg<double> xAxisSpacingArg("x","xAxisSpacing",
    "desired x axis spacing.", false, 0.0,
    "double floating point value of the desired x axis spacing.",cmd);
      TCLAP::UnlabeledValueArg<std::string> inputImageFileArg("inputImageFile","may be in any ITK format.",true,"","input image file",cmd);
      TCLAP::UnlabeledValueArg<std::string> interpolationStyleArg("interpolationStyle","one of linear, nearest, bspline, sinc.",true,"","interpolation style",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputImageFileArg("outputImageFile","result is written here.",true,"","output image file",cmd);

      cmd.parse(argc,argv);

      if (inputImageFileArg.isSet()) inputImageFile = new std::string(inputImageFileArg.getValue());
      if (interpolationStyleArg.isSet()) interpolationStyle = new std::string(interpolationStyleArg.getValue());
      if (outputImageFileArg.isSet()) outputImageFile = new std::string(outputImageFileArg.getValue());

      if (xAxisSpacingArg.isSet()) {
        AxisSpacingSet[0] = true;
        AxisSpacingValue[0] = xAxisSpacingArg.getValue();
      }
      if (yAxisSpacingArg.isSet()) {
        AxisSpacingSet[1] = true;
        AxisSpacingValue[1] = yAxisSpacingArg.getValue();
      }
      if (zAxisSpacingArg.isSet()) {
        AxisSpacingSet[2] = true;
        AxisSpacingValue[2] = zAxisSpacingArg.getValue();
      }
    }
  catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      return EXIT_FAILURE;
    }

    // Sometimes it is helpful to turn off parallel resampling when debugging:
    // itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);

  const unsigned int Dimension = DIMENSIONALITY;
  typedef float PixelType;

  typedef itk::OrientedImage< PixelType,  Dimension >   ImageType;
  typedef itk::OrientedImage< PixelType,  Dimension >   OutImageType;

  typedef itk::ImageFileReader< ImageType  >  ReaderType;
  typedef itk::ImageFileWriter< OutImageType  >  WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( inputImageFile->c_str() );
  writer->SetFileName( outputImageFile->c_str() );

  itk::OrientImageFilter<ImageType,OutImageType>::Pointer imageOrienter =
            itk::OrientImageFilter<ImageType,OutImageType>::New();
  imageOrienter->UseImageDirectionOn();
  imageOrienter->SetDesiredCoordinateOrientationToAxial();
  imageOrienter->SetInput( reader->GetOutput() );
  try {
    imageOrienter->Update();
  } catch ( itk::ExceptionObject & excp )
  {
    std::cerr << "Error while orienting the input file " <<
                   inputImageFile << std::endl;
    std::cerr << excp << std::endl;
    std::cerr << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::ResampleImageFilter<
                  ImageType, OutImageType >  FilterType;
  FilterType::Pointer filter = FilterType::New();

  typedef itk::AffineTransform< double, Dimension >  AffineTransformType;
  AffineTransformType::Pointer identityTransform = AffineTransformType::New();
  identityTransform->SetIdentity();


  typedef itk::ConstantBoundaryCondition< ImageType >  BoundaryConditionType;
  const unsigned int WindowRadius = 5;
  typedef itk::Function::HammingWindowFunction<WindowRadius>  WindowFunctionType;

  typedef itk::WindowedSincInterpolateImageFunction<
                                          ImageType,
                                          WindowRadius,
                                          WindowFunctionType,
                                          BoundaryConditionType,
                                          double  >    SincInterpolatorType;

  SincInterpolatorType::Pointer sincinterpolator = SincInterpolatorType::New();

  typedef itk::BSplineInterpolateImageFunction<
                       ImageType, double >  BSplineInterpolatorType;
  BSplineInterpolatorType::Pointer bsplineInterpolator =
                       BSplineInterpolatorType::New();
  typedef itk::LinearInterpolateImageFunction<
                       ImageType, double >  LinearInterpolatorType;
  LinearInterpolatorType::Pointer linearInterpolator =
                       LinearInterpolatorType::New();
  typedef itk::NearestNeighborInterpolateImageFunction<
                       ImageType, double >  NearestNeighborInterpolatorType;
  NearestNeighborInterpolatorType::Pointer nninterpolator =
                       NearestNeighborInterpolatorType::New();

  // Set to the requested type of interpolator
  if (*interpolationStyle == "linear") {
    filter->SetInterpolator ( linearInterpolator ) ;
  } else if (*interpolationStyle == "nearest") {
    filter->SetInterpolator ( nninterpolator ) ;
  } else if (*interpolationStyle == "bspline") {
    filter->SetInterpolator ( bsplineInterpolator ) ;
  } else if (*interpolationStyle == "sinc") {
    filter->SetInterpolator ( sincinterpolator ) ;
  } else {
    std::cerr << "Don't know how to interpolate using " << (*interpolationStyle) << std::endl;
    std::cerr << "Options are: linear nearestneighbor bspline sinc"<< std::endl;
    return EXIT_FAILURE;
  }

  OutImageType::IndexType outIndex;
  OutImageType::SizeType outSize;
  OutImageType::SpacingType outSpacing;

  for (unsigned int i = 0; i < Dimension; i++) {
    outIndex[i] = imageOrienter->GetOutput()->GetLargestPossibleRegion().GetIndex()[i];
    outSize[i] = imageOrienter->GetOutput()->GetLargestPossibleRegion().GetSize()[i];
    outSpacing[i] = imageOrienter->GetOutput()->GetSpacing()[i];
  }

  // If no spacings are requested, set the spacing to be isotropic and as
  // small as possible.
  bool spacingSet = false;
  for (unsigned int i = 0; i < DIMENSIONALITY; i++) {
    if (AxisSpacingSet[i]) {
      spacingSet = true;
      outSpacing[i] = AxisSpacingValue[i];
    }
  }
  if (!spacingSet) {
    double minSpacing = outSpacing[0];
    for (unsigned int i = 1; i < DIMENSIONALITY; i++) {
      if (minSpacing > outSpacing[i]) {
        minSpacing = outSpacing[i];
      }
    }
    for (unsigned int i = 0; i < DIMENSIONALITY; i++) {
      outSpacing[i] = minSpacing;
    }
  }

  // Set the size to be a rectangle sufficiently large to hold the data.
  // Depending on the values, rounding could occur here.
  for (unsigned int i = 0; i < Dimension; i++) {
    outSize[i] = static_cast<unsigned long int>(
                   (imageOrienter->GetOutput()->GetSpacing()[i] * 
                     outSize[i])/ outSpacing[i] );
  }

  /* Use this if you want to make a cube.
  if ( (outSize[0] >= outSize[1]) && (outSize[0] >= outSize[2]) ) {
    outSize[2] = outSize[1] = outSize[0];
  } else if ( (outSize[1] >= outSize[0]) && (outSize[1] >= outSize[2]) ) {
    outSize[2] = outSize[0] = outSize[1];
  } else if ( (outSize[2] >= outSize[0]) && (outSize[2] >= outSize[1]) ) {
    outSize[0] = outSize[1] = outSize[2];
  }
  */

  filter->SetDefaultPixelValue ( 0.0 );

  filter->SetOutputOrigin( imageOrienter->GetOutput()->GetOrigin() );
  filter->SetOutputDirection( imageOrienter->GetOutput()->GetDirection() );
  filter->SetOutputStartIndex( outIndex );
  filter->SetOutputSpacing( outSpacing );
  filter->SetSize( outSize );

  filter->SetInput( imageOrienter->GetOutput() );
  writer->SetInput( filter->GetOutput() );
  writer->UseCompressionOn( );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Failed trying to write file " <<
               writer->GetFileName() << std::endl;
    std::cerr << "Exception occurred!" << std::endl;
    std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;

}
