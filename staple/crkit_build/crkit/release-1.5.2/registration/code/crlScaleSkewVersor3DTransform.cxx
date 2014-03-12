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

// This program will write out a transform
// with the provided parameter values, for a
// ScaleSkewVersor3D transform representation.
//
// It would benefit from being made to work with arbitrary dimensionality 
// transforms, and to be able to generate different types of identity 
// transform.

// IO mechanisms
#include <itkTransformFileWriter.h>

// Transforms supported
#include <itkScaleSkewVersor3DTransform.h>
#include "configuration.h"

int main( int argc, char * argv[] )
{
  if (( argc != 2 ) && (argc != 3))
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " outputTransformFile" << std::endl;
    std::cerr << "?optional y skew angle(gantry tilt, radians)?" << std::endl;
    std::cerr 
	<< std::endl 
	<< "Write out a ScaleSkewVersor3DTransform transform.  Useful for in-place resampling of images to a new geometry." 
	<< std::endl;

    std::cerr << "The serialization of the optimizable parameters is an array of 15 elements. The first 3 elements are the components of the versor representation of 3D rotation. The next 3 parameters defines the translation in each dimension. The next 3 parameters defines scaling in each dimension. The last 6 parameters defines the skew." << std::endl;

    return 1;
    }


  const char* outfile = argv[1];

  typedef itk::ScaleSkewVersor3DTransform< double >  AffineTransformType;

  AffineTransformType::Pointer outputTransform = AffineTransformType::New();
  outputTransform->SetIdentity();

  if (argc == 3) {
    // Attempt to set the skew using the specified gantry tilt
    double angle = std::atof( argv[2] );
    itk::Vector<double, 6> skews;
    /* Parameters map to : 
  newMatrix[0][1] = 2.0 * ( xy - zw )  + ( m_Skew[0] );
  newMatrix[0][2] = 2.0 * ( xz + yw )  + ( m_Skew[1] );
  newMatrix[1][0] = 2.0 * ( xy + zw )  + ( m_Skew[2] );
  newMatrix[1][2] = 2.0 * ( yz - xw )  + ( m_Skew[3] );
  newMatrix[2][0] = 2.0 * ( xz - yw )  + ( m_Skew[4] );
  newMatrix[2][1] = 2.0 * ( yz + xw )  + ( m_Skew[5] );
     */

    skews[0] = 0.0; skews[1] = 0.0;
    skews[2] = 0.0;        // shx 
    skews[3] = tan(angle); // shz
    skews[4] = 0.0; skews[5] = 0.0;
    outputTransform->SetSkew( skews );
  }

  itk::TransformFileWriter::Pointer trsfwriter;
  trsfwriter = itk::TransformFileWriter::New();
  trsfwriter->SetFileName( outfile );

  trsfwriter->SetInput ( outputTransform );

  try
    {
    trsfwriter->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception occurred writing out the transform to file " 
      << outfile << std::endl;
    std::cerr << excep << std::endl;
    return 1;
    }

  return 0;
}
