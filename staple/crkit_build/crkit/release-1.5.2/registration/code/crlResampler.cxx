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

// This program will read in a transform, convert it
// to an affine representation and then run a resample filter with it
// on the input image data.

// IO mechanisms
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkOrientedImage.h>
#include <itkOrientImageFilter.h>

// Transforms supported
#include <itkTranslationTransform.h>
#include <itkVersorRigid3DTransform.h>
#include <itkSimilarity3DTransform.h>
#include <itkScaleSkewVersor3DTransform.h>
#include <itkAffineTransform.h>

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
  std::string *inputTransformFile = new std::string();
  std::string *outputGeometryExampleFile = new std::string();
  std::string *interpolationStyle = new std::string();
  std::string *outputImageFile = new std::string();
  std::string *imageProvidingOrigin = new std::string();
  bool useImageProvidingOrigin = false;
  bool computeOutputOriginFromInput = false;

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',
                          CRKIT_VERSION_STRING  );
      TCLAP::UnlabeledValueArg<std::string> inputImageFileArg("inputImageFile","may be in any ITK format.",true,"","input image file",cmd);
      TCLAP::UnlabeledValueArg<std::string> inputTransformFileArg("inputTransformFile","an ITK transform file.",true,"","input transform file",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputGeometryExampleFileArg("outputGeometryExampleFile","an image file with the desired geometry.",true,"","output geometry example file",cmd);
      TCLAP::UnlabeledValueArg<std::string> interpolationStyleArg("interpolationStyle","one of linear, nearest, bspline, sinc.",true,"","interpolation style",cmd);
      TCLAP::SwitchArg computeOriginFromInputArg("c",
              "computeOriginFromInput",
              "Compute the output origin by transforming the input origin",
              false);
      cmd.add(computeOriginFromInputArg);
      TCLAP::ValueArg<std::string> alternativeOriginImageArg("a","alternativeOriginImage","Image to take origin from",false,"","image providing origin (string file name)",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputImageFileArg("outputImageFile","result is written here.",true,"","output image file",cmd);

      cmd.parse(argc,argv);

      if (inputImageFileArg.isSet()) inputImageFile = new std::string(inputImageFileArg.getValue());
      if (inputTransformFileArg.isSet()) inputTransformFile = new std::string(inputTransformFileArg.getValue());
      if (outputGeometryExampleFileArg.isSet()) outputGeometryExampleFile = new std::string(outputGeometryExampleFileArg.getValue());
      if (interpolationStyleArg.isSet()) interpolationStyle = new std::string(interpolationStyleArg.getValue());
      if (outputImageFileArg.isSet()) outputImageFile = new std::string(outputImageFileArg.getValue());
      if (alternativeOriginImageArg.isSet()) {
        imageProvidingOrigin = new 
                             std::string(alternativeOriginImageArg.getValue());
        useImageProvidingOrigin = true;
      }
      if (computeOriginFromInputArg.getValue()) {
        computeOutputOriginFromInput = true;
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

  typedef itk::ImageFileReader< ImageType  >  ReaderType;
  typedef itk::ImageFileWriter< ImageType  >  WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  ReaderType::Pointer geomreader = ReaderType::New();

  reader->SetFileName( inputImageFile->c_str() );
  writer->SetFileName( outputImageFile->c_str() );

  itk::OrientImageFilter<ImageType,ImageType>::Pointer imageOrienter = 
            itk::OrientImageFilter<ImageType,ImageType>::New();
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

  geomreader->SetFileName( outputGeometryExampleFile->c_str() );
  itk::OrientImageFilter<ImageType,ImageType>::Pointer geomimageOrienter = 
            itk::OrientImageFilter<ImageType,ImageType>::New();
  geomimageOrienter->UseImageDirectionOn();
  geomimageOrienter->SetDesiredCoordinateOrientationToAxial();
  geomimageOrienter->SetInput( geomreader->GetOutput() );
  try {
    geomimageOrienter->Update();
  } catch ( itk::ExceptionObject & excp )
  {
    std::cerr << "Error while reading and orienting the geometry file " <<
                  (*outputGeometryExampleFile) << std::endl;
    std::cerr << excp << std::endl;
    std::cerr << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::ResampleImageFilter<
                  ImageType, ImageType >  FilterType;
  FilterType::Pointer filter = FilterType::New();

  typedef itk::AffineTransform< double, Dimension >  AffineTransformType;
  typedef itk::ScaleSkewVersor3DTransform< double >
                   ScaleSkewVersor3DTransformType;
  typedef itk::Similarity3DTransform< double >
                   Similarity3DTransformType;
  typedef itk::VersorRigid3DTransform< double >
                   VersorRigid3DTransformType;
  typedef itk::TranslationTransform< double, Dimension >
                   TranslationTransformType;


  itk::TransformFileReader::Pointer trsfreader;
  trsfreader = itk::TransformFileReader::New();
  trsfreader->SetFileName( inputTransformFile->c_str() );

  try {
    trsfreader->Update();
  } catch ( itk::ExceptionObject & excp )
  {
    // Display error from reading the transform file.
    std::cerr << "Error while reading the transform file " <<
                  (*inputTransformFile) << std::endl;
    std::cerr << excp << std::endl;
    std::cerr << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
  }

  // Now try to work out how many and what type of transforms were read.
  // We only want to get one transform.
  typedef itk::TransformFileReader::TransformListType * TransformListType;
  TransformListType transforms = trsfreader->GetTransformList();
  std::cout << "Number of transforms = " << transforms->size() << std::endl;
  itk::TransformFileReader::TransformListType::const_iterator it =
               transforms->begin();
  if (transforms->size() <= 0 || transforms->size() > 1) {
    std::cerr << "Read " << transforms->size() << " transforms but want 1." << std::endl;
    return EXIT_FAILURE;
  }

  // Let's see if we can keep around an inverse transform too
  AffineTransformType::Pointer affine_inverse = 0;

  if (!strcmp((*it)->GetNameOfClass(), "AffineTransform"))
  {
    AffineTransformType::Pointer affine_read = static_cast<AffineTransformType*>((*it).GetPointer());
    affine_read->Print(std::cout);
    affine_inverse = AffineTransformType::New();
    filter->SetTransform( affine_read );
    if(!affine_read->GetInverse(affine_inverse)) {
      std::cerr << "Input transform is not invertible." << std::endl;
    }
  } else if (!strcmp((*it)->GetNameOfClass(), "TranslationTransform")) {
    TranslationTransformType::Pointer  translation_read =
        static_cast<TranslationTransformType*>((*it).GetPointer());
    filter->SetTransform( translation_read );
  } else if (!strcmp((*it)->GetNameOfClass(), "VersorRigid3DTransform")) {
    VersorRigid3DTransformType::Pointer  versorrigid_read =
        static_cast<VersorRigid3DTransformType*>((*it).GetPointer());
    filter->SetTransform( versorrigid_read );
  } else if (!strcmp((*it)->GetNameOfClass(), "Similarity3DTransform")) {
    Similarity3DTransformType::Pointer  similarity_read =
        static_cast<Similarity3DTransformType*>((*it).GetPointer());
    filter->SetTransform( similarity_read );
  } else if (!strcmp((*it)->GetNameOfClass(), "ScaleSkewVersor3DTransform")) {
    ScaleSkewVersor3DTransformType::Pointer  scaleskewversor_read =
        static_cast<ScaleSkewVersor3DTransformType*>((*it).GetPointer());
    filter->SetTransform( scaleskewversor_read );
  } else {
    std::cerr << "Don't know how to resample using a " <<
            (*it)->GetNameOfClass() << " transform." << std::endl;
    return EXIT_FAILURE;
  }

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
  // argv[4]
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

  filter->SetDefaultPixelValue ( 0.0 );

  /* Since we may be working with oriented images, ensure that the output
   * direction is set properly. Recent releases of ITK (e.g. 3.4.0) now 
   * correctly set the output direction when setting the output parameters
   * from an image.
   */
  //filter->SetOutputDirection(geomimageOrienter->GetOutput()->GetDirection());
  filter->SetOutputParametersFromImage( geomimageOrienter->GetOutput() );

  if (useImageProvidingOrigin) {
    /* This provides a mechanism to specify the spacing and the size of an
     * output image using one image, but to select the origin by taking it
     * from a different image.
     */
    ReaderType::Pointer originImageReader = ReaderType::New();
    originImageReader->SetFileName( imageProvidingOrigin->c_str() );
    try {
      originImageReader->Update();
    } catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Error while reading alternative origin input file " <<
                   imageProvidingOrigin << std::endl;
      std::cerr << excp << std::endl;
      std::cerr << "[FAILED]" << std::endl;
      return EXIT_FAILURE;
    }
    // The origin is supposed to be provided always in LPS coordinates, so
    // orientation should play no role here.
    filter->SetOutputOrigin( originImageReader->GetOutput()->GetOrigin() );
  }

  if (computeOutputOriginFromInput) {
    // Here we allow the input origin shifted by the transform to define the
    // new origin. This is useful if some or all of the original geometry is
    // to be preserved, but in the new location indicated by the transform.
    FilterType::OriginPointType origin = filter->GetOutputOrigin(); 
    FilterType::OriginPointType neworigin; 
    AffineTransformType::Pointer transformInverse = AffineTransformType::New();
    if (!affine_inverse) {
      std::cerr << "Origin computation can only be done with an affine input transform." << std::endl;
      return EXIT_FAILURE;
    }
    neworigin = affine_inverse->TransformPoint( origin );
    filter->SetOutputOrigin( neworigin );
std::cout << "Initial output origin is " << origin << std::endl;
std::cout << "New output origin is " << neworigin << std::endl;
  }

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
