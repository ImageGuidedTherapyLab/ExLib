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
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#ifdef WIN32
#include <missingFunctionsWindows.h>
#else
#include <sys/time.h>
#endif

#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkAffineTransform.h>

#include <itkOrientedImage.h>
#include <itkOrientImageFilter.h>

#include <vnl/vnl_inverse.h>
#include <itkImageMomentsCalculator.h>

#include <itkImageFileReader.h>

#include <tclap/CmdLine.h>
#include "configuration.h"

int main( int argc, char *argv[] )
{
  std::string *fixedImageFile = new std::string("");
  std::string *movingImageFile = new std::string("");
  std::string *outputTransformFile = new std::string("");

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',
                         CRKIT_VERSION_STRING );
      TCLAP::UnlabeledValueArg<std::string> fixedImageArg("fixedImage","Fixed Image File Name",true,"","fixed image file name",cmd);
      TCLAP::UnlabeledValueArg<std::string> movingImageArg("movingImage","Moving Image File Name",true,"","moving image file name",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputTransformArg("outputTransform","Output File for Transform",true,"","transform file name",cmd);

      cmd.parse(argc,argv);

      if (fixedImageArg.isSet()) fixedImageFile = new std::string(fixedImageArg.getValue());
      if (movingImageArg.isSet()) movingImageFile = new std::string(movingImageArg.getValue());
      if (outputTransformArg.isSet()) outputTransformFile = new std::string(outputTransformArg.getValue());
    }
  catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(1);
    }

  const    unsigned int    Dimension = 3;
  typedef  float           PixelType;

  typedef itk::OrientedImage< PixelType, Dimension >  FixedImageType;
  typedef itk::OrientedImage< PixelType, Dimension >  MovingImageType;

  typedef   float     InternalPixelType;
  typedef itk::OrientedImage< InternalPixelType, Dimension > InternalImageType;

  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  itk::OrientImageFilter<FixedImageType,FixedImageType>::Pointer fixedImageOrienter = itk::OrientImageFilter<FixedImageType,FixedImageType>::New();

  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
  itk::OrientImageFilter<MovingImageType,MovingImageType>::Pointer movingImageOrienter = itk::OrientImageFilter<MovingImageType,MovingImageType>::New();

  fixedImageReader->SetFileName(  fixedImageFile->c_str() );
  fixedImageOrienter->UseImageDirectionOn();
  fixedImageOrienter->SetDesiredCoordinateOrientationToAxial();
  fixedImageOrienter->SetInput( fixedImageReader->GetOutput() );
  fixedImageOrienter->Update();

  movingImageReader->SetFileName( movingImageFile->c_str() );
  movingImageOrienter->UseImageDirectionOn();
  movingImageOrienter->SetDesiredCoordinateOrientationToAxial();
  movingImageOrienter->SetInput( movingImageReader->GetOutput() );
  movingImageOrienter->Update();

  typedef itk::AffineTransform< double, 3> AffineTransformType;
  AffineTransformType::Pointer affine_read = AffineTransformType::New();

  /* Now initialize the rotation based on the moments of the data. */
  typedef itk::ImageMomentsCalculator< FixedImageType >   FixedImageCalculatorType;
  typedef itk::ImageMomentsCalculator< MovingImageType >  MovingImageCalculatorType;
  FixedImageCalculatorType::Pointer fixedImageCalculator = FixedImageCalculatorType::New();
  MovingImageCalculatorType::Pointer movingImageCalculator = MovingImageCalculatorType::New();
    fixedImageCalculator->SetImage( fixedImageOrienter->GetOutput() );
    fixedImageCalculator->Compute();

    FixedImageCalculatorType::VectorType fixedCenter = fixedImageCalculator->GetCenterOfGravity();

    movingImageCalculator->SetImage( movingImageOrienter->GetOutput() );
    movingImageCalculator->Compute();
    MovingImageCalculatorType::VectorType movingCenter = movingImageCalculator->GetCenterOfGravity();

    // Now we need to compute the initial alignment transform.

  AffineTransformType::Pointer fixedAffine = AffineTransformType::New();
  fixedAffine = fixedImageCalculator->GetPrincipalAxesToPhysicalAxesTransform();

  AffineTransformType::Pointer movingAffine = AffineTransformType::New();
  movingAffine = movingImageCalculator->GetPrincipalAxesToPhysicalAxesTransform();

  std::cout << "Fixed affine transform : " << std::endl 
            << fixedAffine << std::endl;
  std::cout << "Moving affine transform : " << std::endl
            << movingAffine << std::endl;

  // The principal axes of the moving image give us a transform between the
  // intrinsic coordinate system of the moving image and the default space.
  //    This transform is Pm, and is exactly the principal axes matrix.
  // The principal axes of the fixed image give us a transform between the
  // intrinsic coordinate system of the fixed image and the default space.
  //    This transform is Pf, and is exactly the principal axes matrix.
  //   In order to align the principal axes of the moving image to the fixed 
  // image, we need to form: (Pf)^{-1}(Pm)
  //
  // We can then estimate the rotation by decomposing this change of basis
  // matrix into a versor.  
  //
  AffineTransformType::Pointer combinedTransform = AffineTransformType::New();
  fixedAffine->GetInverse(combinedTransform);
  combinedTransform->Compose(movingAffine, true);

  // Recall that we are estimating the transform that tells us where in the
  // moving data the points in the fixed data set come from. Therefore we
  // initialize with the inverse of the rotation from the moving to fixed data.
  combinedTransform->GetInverse(combinedTransform);
  std::cout << "Initial Affine Matrix = " << std::endl;
  combinedTransform->Print(std::cout);
  std::cout << std::endl;

  // Now write the transform to a file.
  itk::TransformFileWriter::Pointer trsfWriter;
  trsfWriter = itk::TransformFileWriter::New();
  trsfWriter->SetInput( combinedTransform );
  trsfWriter->SetFileName( outputTransformFile->c_str() );

  try
  {
    trsfWriter->Update();
  }
  catch( itk::ExceptionObject & excp)
  {
    std::cerr << "Error while saving the transform to file " 
              << (*outputTransformFile) << std::endl;
    std::cerr << excp << std::endl;
    return 0;
  }

  return 0;
}
