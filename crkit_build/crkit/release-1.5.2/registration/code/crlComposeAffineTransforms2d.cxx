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
  if( argc != 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputTransformFile1" ;
      std::cerr << " inputTransformFile2 outputTransformFile" << std::endl;
      return 1;
    }

  const unsigned int Dimension = 2;

  typedef itk::AffineTransform< double, Dimension >  AffineTransformType;
  typedef itk::CenteredRigid2DTransform< double >  CenteredRigid2DTransformType;
  typedef itk::Similarity2DTransform< double >
                   Similarity2DTransformType;
  typedef itk::TranslationTransform< double, Dimension >
                   TranslationTransformType;

  typedef itk::AffineTransform< double, Dimension >  AffineTransformType;

  // In order to do IO with these transforms, they need to be registered
  // with the transform factory.  Some of these are done by default, but some
  // are not.
  itk::TransformFactory<Similarity2DTransformType>::RegisterDefaultTransforms();
  itk::TransformFactory<Similarity2DTransformType>::RegisterTransform();
  itk::TransformFactory<CenteredRigid2DTransformType>::RegisterTransform();

  AffineTransformType::Pointer outputTransform = AffineTransformType::New();
  outputTransform->SetIdentity();

  AffineTransformType::Pointer input2Transform = AffineTransformType::New();
  input2Transform->SetIdentity();

  itk::TransformFileReader::Pointer trsfreader1;
  itk::TransformFileReader::Pointer trsfreader2;
  trsfreader1 = itk::TransformFileReader::New();
  trsfreader2 = itk::TransformFileReader::New();
  trsfreader1->SetFileName( argv[1] );
  trsfreader2->SetFileName( argv[2] );

  try {
    trsfreader1->Update();
  } catch ( itk::ExceptionObject & excp )
  {
    std::cerr << "Error while reading the transform file " << 
                  argv[1] << std::endl;
    std::cerr << excp << std::endl;
    std::cerr << "[FAILED]" << std::endl;
    return 1;
  }
  try {
    trsfreader2->Update();
  } catch ( itk::ExceptionObject & excp )
  {
    std::cerr << "Error while reading the transform file " << 
                  argv[2] << std::endl;
    std::cerr << excp << std::endl;
    std::cerr << "[FAILED]" << std::endl;
    return 1;
  }

  // Now try to work out how many and what type of transforms were read.
  // We only want to get one transform.
  typedef itk::TransformFileReader::TransformListType * TransformListType;

  TransformListType transforms = trsfreader1->GetTransformList();
  std::cout << "Number of transforms = " << transforms->size() << std::endl;
  itk::TransformFileReader::TransformListType::const_iterator it = 
               transforms->begin();
  if (transforms->size() <= 0 || transforms->size() > 1) {
    std::cerr << "Read " << transforms->size() << " transforms but want 1." << std::endl;
    return 1;
  }
  if (!strcmp((*it)->GetNameOfClass(), "AffineTransform"))
  {
    AffineTransformType::Pointer affine_read = static_cast<AffineTransformType*>((*it).GetPointer());
    affine_read->Print(std::cout);

    // Setting the matrix and offset is sufficient to duplicate the transform 
    // performance operating on points, 
    // but does not set the center appropriately for follow on registration.
    outputTransform->SetTranslation(affine_read->GetTranslation());
    outputTransform->SetCenter(affine_read->GetCenter());
    outputTransform->SetMatrix(affine_read->GetMatrix());

  } else if (!strcmp((*it)->GetNameOfClass(), "TranslationTransform")) {
    TranslationTransformType::Pointer  translation_read = 
        static_cast<TranslationTransformType*>((*it).GetPointer());

    outputTransform->SetOffset(translation_read->GetOffset());

  } else if (!strcmp((*it)->GetNameOfClass(), "Similarity2DTransform")) {
    Similarity2DTransformType::Pointer  similarity_read = 
        static_cast<Similarity2DTransformType*>((*it).GetPointer());
    outputTransform->SetCenter(similarity_read->GetCenter());
    outputTransform->SetTranslation(similarity_read->GetTranslation());
    outputTransform->SetMatrix(similarity_read->GetMatrix());
  } else if (!strcmp((*it)->GetNameOfClass(), "CenteredRigid2DTransform")) {
    CenteredRigid2DTransformType::Pointer t_read = static_cast<CenteredRigid2DTransformType*>((*it).GetPointer());
    t_read->Print(std::cout);

    // Setting the matrix and offset is sufficient to duplicate the transform 
    // performance operating on points, 
    // but does not set the center appropriately for follow on registration.
    outputTransform->SetTranslation(t_read->GetTranslation());
    outputTransform->SetCenter(t_read->GetCenter());
    outputTransform->SetMatrix(t_read->GetMatrix());

  } else {
    std::cerr << "Don't know how to convert a " << 
            (*it)->GetNameOfClass() << " transform." << std::endl;
    return 1;
  }

  transforms = trsfreader2->GetTransformList();
  std::cout << "Number of transforms = " << transforms->size() << std::endl;
  it = transforms->begin();
  if (transforms->size() <= 0 || transforms->size() > 1) {
    std::cerr << "Read " << transforms->size() << " transforms but want 1." << std::endl;
    return 1;
  }
  if (!strcmp((*it)->GetNameOfClass(), "AffineTransform"))
  {
    AffineTransformType::Pointer affine_read = static_cast<AffineTransformType*>((*it).GetPointer());
    affine_read->Print(std::cout);

    // Setting the matrix and offset is sufficient to duplicate the transform 
    // performance operating on points, 
    // but does not set the center appropriately for follow on registration.
    input2Transform->SetTranslation(affine_read->GetTranslation());
    input2Transform->SetCenter(affine_read->GetCenter());
    input2Transform->SetMatrix(affine_read->GetMatrix());

  } else if (!strcmp((*it)->GetNameOfClass(), "TranslationTransform")) {
    TranslationTransformType::Pointer  translation_read = 
        static_cast<TranslationTransformType*>((*it).GetPointer());

    input2Transform->SetOffset(translation_read->GetOffset());

  } else if (!strcmp((*it)->GetNameOfClass(), "Similarity2DTransform")) {
    Similarity2DTransformType::Pointer  similarity_read = 
        static_cast<Similarity2DTransformType*>((*it).GetPointer());
    input2Transform->SetCenter(similarity_read->GetCenter());
    input2Transform->SetTranslation(similarity_read->GetTranslation());
    input2Transform->SetMatrix(similarity_read->GetMatrix());
  } else if (!strcmp((*it)->GetNameOfClass(), "CenteredRigid2DTransform")) {
    CenteredRigid2DTransformType::Pointer t_read = static_cast<CenteredRigid2DTransformType*>((*it).GetPointer());
    t_read->Print(std::cout);

    // Setting the matrix and offset is sufficient to duplicate the transform 
    // performance operating on points, 
    // but does not set the center appropriately for follow on registration.
    input2Transform->SetTranslation(t_read->GetTranslation());
    input2Transform->SetCenter(t_read->GetCenter());
    input2Transform->SetMatrix(t_read->GetMatrix());
  } else {
    std::cerr << "Don't know how to convert a " << 
            (*it)->GetNameOfClass() << " transform." << std::endl;
    return 1;
  }


  itk::TransformFileWriter::Pointer trsfwriter;
  trsfwriter = itk::TransformFileWriter::New();
  trsfwriter->SetFileName( argv[3] );

  // Can be a pre or post composition...
  // outputTransform is initialized with input1Transform.
  //   Here it is composed with input2Transform.
  // This method composes self with another MatrixOffsetTransformBase
  // of the same dimension, modifying self to be the composition of self
  // and other. If the argument pre is true, then other is precomposed with
  // self; that is, the resulting transformation consists of first applying
  // other to the source, followed by self. If pre is false or omitted, then
  // other is post-composed with self; that is the resulting transformation
  // consists of first applying self to the source, followed by other. This
  // updates the Translation based on current center.
  outputTransform->Compose( input2Transform );

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
