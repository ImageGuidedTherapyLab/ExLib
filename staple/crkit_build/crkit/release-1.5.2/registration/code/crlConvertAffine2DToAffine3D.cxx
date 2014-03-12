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


// This program will read in two affine transforms and compose them.

// IO mechanisms
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkTransformFactory.h>

// Transforms supported
#include <itkTranslationTransform.h>
#include <itkCenteredRigid2DTransform.h>
#include <itkSimilarity2DTransform.h>
#include <itkAffineTransform.h>
#include "configuration.h"

int main( int argc, char * argv[] )
{
  if( argc != 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputTransformFile" ;
      std::cerr << " outputTransformFile" << std::endl;
      return 1;
    }

  const unsigned int InDimension = 2;
  const unsigned int OutDimension = 3;

  typedef itk::AffineTransform< double, InDimension >  InAffineTransformType;
  typedef itk::CenteredRigid2DTransform< double >  InCenteredRigid2DTransformType;
  typedef itk::AffineTransform< double, OutDimension >  OutAffineTransformType;

  // In order to do IO with these transforms, they need to be registered
  // with the transform factory.  Some of these are done by default, but some
  // are not.
  itk::TransformFactory<InAffineTransformType>::RegisterDefaultTransforms();
  itk::TransformFactory<InCenteredRigid2DTransformType>::RegisterTransform();

  OutAffineTransformType::Pointer outputTransform = OutAffineTransformType::New();
  outputTransform->SetIdentity();

  InAffineTransformType::Pointer inputTransform = InAffineTransformType::New();
  inputTransform->SetIdentity();

  itk::TransformFileReader::Pointer trsfreader;
  trsfreader = itk::TransformFileReader::New();
  trsfreader->SetFileName( argv[1] );

  try {
    trsfreader->Update();
  } catch ( itk::ExceptionObject & excp )
  {
    std::cerr << "Error while reading the transform file " << 
                  argv[1] << std::endl;
    std::cerr << excp << std::endl;
    std::cerr << "[FAILED]" << std::endl;
    return 1;
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
    return 1;
  }
  if (!strcmp((*it)->GetNameOfClass(), "AffineTransform"))
  {
    InAffineTransformType::Pointer affine_read = static_cast<InAffineTransformType*>((*it).GetPointer());
    affine_read->Print(std::cout);

    // Setting the matrix and offset is sufficient to duplicate the transform 
    // performance operating on points, 
    // but does not set the center appropriately for follow on registration.
    
    InAffineTransformType::OutputVectorType intrans = 
                                             affine_read->GetTranslation();
    OutAffineTransformType::OutputVectorType outtrans;
          outtrans[0] = intrans[0];
          outtrans[1] = intrans[1];
          outtrans[2] = 0.0;
    outputTransform->SetTranslation(outtrans);

    InAffineTransformType::InputPointType inCenter = affine_read->GetCenter();
    OutAffineTransformType::InputPointType outCenter ;
    outCenter[0] = inCenter[0];
    outCenter[1] = inCenter[1];
    outCenter[2] = 0.0;
    InAffineTransformType::MatrixType inMatrix = affine_read->GetMatrix();
    OutAffineTransformType::MatrixType outMatrix = outputTransform->GetMatrix();
    for (unsigned int i = 0; i < InDimension; i++) {
      for (unsigned int j = 0; j < InDimension; j++) {
        outMatrix(i,j) = inMatrix(i,j);
      }
    }
    outputTransform->SetMatrix(outMatrix);

  } else if (!strcmp((*it)->GetNameOfClass(), "CenteredRigid2DTransform")) {
    InCenteredRigid2DTransformType::Pointer t_read = static_cast<InCenteredRigid2DTransformType*>((*it).GetPointer());
    t_read->Print(std::cout);

    InAffineTransformType::OutputVectorType intrans = 
                                             t_read->GetTranslation();
    OutAffineTransformType::OutputVectorType outtrans;
          outtrans[0] = intrans[0];
          outtrans[1] = intrans[1];
          outtrans[2] = 0.0;
    outputTransform->SetTranslation(outtrans);

    InAffineTransformType::InputPointType inCenter = t_read->GetCenter();
    OutAffineTransformType::InputPointType outCenter ;
    outCenter[0] = inCenter[0];
    outCenter[1] = inCenter[1];
    outCenter[2] = 0.0;
    InAffineTransformType::MatrixType inMatrix = t_read->GetMatrix();
    OutAffineTransformType::MatrixType outMatrix = outputTransform->GetMatrix();
    for (unsigned int i = 0; i < InDimension; i++) {
      for (unsigned int j = 0; j < InDimension; j++) {
        outMatrix(i,j) = inMatrix(i,j);
      }
    }
    outputTransform->SetMatrix(outMatrix);

  } else {
    std::cerr << "Don't know how to convert a " << 
            (*it)->GetNameOfClass() << " transform." << std::endl;
    return 1;
  }

  itk::TransformFileWriter::Pointer trsfwriter;
  trsfwriter = itk::TransformFileWriter::New();
  trsfwriter->SetFileName( argv[2] );

  trsfwriter->SetInput ( outputTransform );

  try
    {
    trsfwriter->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception occurred writing the transform!" << std::endl;
    std::cerr << excep << std::endl;
    }

  return EXIT_SUCCESS;
}
