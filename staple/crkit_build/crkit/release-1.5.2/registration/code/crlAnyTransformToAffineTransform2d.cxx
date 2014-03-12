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


// This program will read in one of the supported transforms and convert it
// to an affine representation.

// IO mechanisms
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>

// Transforms supported
#include <itkTranslationTransform.h>
#include <itkRigid2DTransform.h>
#include <itkSimilarity2DTransform.h>
#include <itkAffineTransform.h>

#include "configuration.h"

int main( int argc, char * argv[] )
{
  if( ( argc != 3 ) && (argc != 4)) 
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputTransformFile " << std::endl;
    std::cerr << " outputTransformFile " ;
    std::cerr << " ?inverse?" << std::endl;
    std::cerr << "If a third argument is specified, the inverse transform "
              << "will be computed." << std::endl;
    return 1;
    }

  const unsigned int Dimension = 2;

  typedef itk::AffineTransform< double, Dimension >  AffineTransformType;
  typedef itk::TranslationTransform< double, Dimension >  
                   TranslationTransformType;
  typedef itk::Similarity2DTransform< double >  Similarity2DTransformType;
  typedef itk::Rigid2DTransform< double >  Rigid2DTransformType;

  AffineTransformType::Pointer outputTransform = AffineTransformType::New();
  outputTransform->SetIdentity();

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

  } else if (!strcmp((*it)->GetNameOfClass(), "Rigid2DTransform")) {
    Rigid2DTransformType::Pointer  rigid_read = 
        static_cast<Rigid2DTransformType*>((*it).GetPointer());
    outputTransform->SetCenter(rigid_read->GetCenter());
    outputTransform->SetTranslation(rigid_read->GetTranslation());
    outputTransform->SetMatrix(rigid_read->GetMatrix());
    // The above is sufficient to duplicate the transform performance
    // but does not set the center appropriately for follow on registration.

  } else if (!strcmp((*it)->GetNameOfClass(), "Similarity2DTransform")) {
    Similarity2DTransformType::Pointer  similarity_read = 
        static_cast<Similarity2DTransformType*>((*it).GetPointer());
    outputTransform->SetCenter(similarity_read->GetCenter());
    outputTransform->SetTranslation(similarity_read->GetTranslation());
    outputTransform->SetMatrix(similarity_read->GetMatrix());
  } else {
    std::cerr << "Don't know how to convert a " << 
            (*it)->GetNameOfClass() << " transform." << std::endl;
    return 1;
  }

  itk::TransformFileWriter::Pointer trsfwriter;
  trsfwriter = itk::TransformFileWriter::New();
  trsfwriter->SetFileName( argv[2] );

  AffineTransformType::Pointer inverse = AffineTransformType::New();

  try
    {
    outputTransform->GetInverse(inverse);
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception occurred computing the inverse transform" 
      << std::endl;
    std::cerr << excep << std::endl;
    }

  if (argc == 4) {
    std::cout << "Writing the inverse affine transform." << std::endl;
    trsfwriter->SetInput ( inverse );
  } else {
    trsfwriter->SetInput ( outputTransform );
  }

  try
    {
    trsfwriter->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception occurred writing out the transform to file " 
      << argv[2] << std::endl;
    std::cerr << excep << std::endl;
    }

}
