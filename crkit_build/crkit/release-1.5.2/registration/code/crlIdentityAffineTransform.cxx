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

// This program will write out an identity transform
// to an affine representation.
//
// It would benefit from being made to work with arbitrary dimensionality 
// transforms, and to be able to generate different types of identity 
// transform.

// IO mechanisms
#include <itkTransformFileWriter.h>

// Transforms supported
#include <itkAffineTransform.h>
#include "configuration.h"

int main( int argc, char * argv[] )
{
  if( (argc != 2 ) && (argc != 3) )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " outputTransformFile ?dimensionality?" << std::endl;
    std::cerr 
	<< std::endl 
	<< "Write out an identity transform.  Useful for in-place resampling of images to a new geometry." 
	<< std::endl;
    return 1;
    }


  const char* outfile = argv[1];

  int dim = 3;
  if (argc == 3) {
    dim = ::atoi( argv[2] );
  }

  std::cout << "Dimensionality is " << dim << std::endl;
 

  itk::TransformFileWriter::Pointer trsfwriter;
  trsfwriter = itk::TransformFileWriter::New();
  trsfwriter->SetFileName( outfile );

  if (dim == 1) {
    typedef itk::AffineTransform< double, 1>  AffineTransformType;
    AffineTransformType::Pointer outputTransform = AffineTransformType::New();
    outputTransform->SetIdentity();
    trsfwriter->SetInput ( outputTransform );
    try { trsfwriter->Update(); }
    catch( itk::ExceptionObject & excep )
    {
      std::cerr << "Exception occurred writing out the transform to file " 
                << outfile << std::endl;
      return 1;
    }
  } else if (dim == 2) {
    typedef itk::AffineTransform< double, 2>  AffineTransformType;
    AffineTransformType::Pointer outputTransform = AffineTransformType::New();
    outputTransform->SetIdentity();
    trsfwriter->SetFileName( outfile );
    trsfwriter->SetInput ( outputTransform );
    try { trsfwriter->Update(); }
    catch( itk::ExceptionObject & excep )
    {
      std::cerr << "Exception occurred writing out the transform to file " 
                << outfile << std::endl;
      return 1;
    }
  } else if (dim == 3) {
    typedef itk::AffineTransform< double, 3>  AffineTransformType;
    AffineTransformType::Pointer outputTransform = AffineTransformType::New();
    outputTransform->SetIdentity();
    trsfwriter->SetFileName( outfile );
    trsfwriter->SetInput ( outputTransform );
    try { trsfwriter->Update(); }
    catch( itk::ExceptionObject & excep )
    {
      std::cerr << "Exception occurred writing out the transform to file " 
                << outfile << std::endl;
      return 1;
    }
  } else if (dim == 4) {
    typedef itk::AffineTransform< double, 4>  AffineTransformType;
    AffineTransformType::Pointer outputTransform = AffineTransformType::New();
    outputTransform->SetIdentity();
    trsfwriter->SetFileName( outfile );
    trsfwriter->SetInput ( outputTransform );
    try { trsfwriter->Update(); }
    catch( itk::ExceptionObject & excep )
    {
      std::cerr << "Exception occurred writing out the transform to file " 
                << outfile << std::endl;
      return 1;
    }
  }

  return 0;
}
