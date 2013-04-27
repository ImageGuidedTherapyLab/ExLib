/*
  Copyright 2007 Simon Warfield simon.warfield@childrens.harvard.edu
  This software is made available under the terms of the Open Software 
  License version 3.0: http://www.opensource.org/licenses/osl-3.0.php

Attribution Notice.

This research was carried out in the Computational Radiology Laboratory of
Children's Hospital, Boston and Harvard Medical School.

This research work was made possible by Grant Number R01 RR021885 (Principal
Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
from the National Center for Research Resources (NCRR), a component of the
National Institutes of Health (NIH).
*/


#include "itkDiffusionTensor3DReconstructionImageFilter.h"
#include "itkVectorImage.h"
#include "itkNrrdImageIO.h"
#include "itkMetaDataDictionary.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <iostream>

int main( int argc, char *argv[] )
{
  if(argc < 3)
    {
    std::cerr << "Usage: " << argv[0] 
     << " DWIfile.nhdr threshold(on B0) tensorOutputFileName " << std::endl;
    return EXIT_FAILURE;
    }
  
  const unsigned int ImageDimension = 3;
  unsigned int numberOfImages = 0;
  bool readb0 = false;
  double b0 = 0;
  
  typedef unsigned short                      InPixelType;
  typedef itk::VectorImage<InPixelType, 3> ImageType;
  typedef float OutPixelType;
  
  itk::ImageFileReader<ImageType>::Pointer reader 
    = itk::ImageFileReader<ImageType>::New();
  
  ImageType::Pointer img;
  
  // Set the properties for NrrdReader
  reader->SetFileName(argv[1]);

  // Read in the DWI format data.
  try
    {
    reader->Update();
    img = reader->GetOutput();
    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << ex << std::endl;
    std::cout << "Failed to read image called " << argv[1] << std::endl;
    return EXIT_FAILURE;
    }

  // This is the filter that fits the tensor model, using the method of 
  // Westin et al. 2002.
  typedef itk::DiffusionTensor3DReconstructionImageFilter< 
    InPixelType, InPixelType, OutPixelType> TensorReconstructionImageFilterType;

  /* 
  // Read the Nrrd headers to get the B value and the gradient directions used
  // for diffusion weighting. 
  // 
  // DWMRI_b-value:=1000
  // DWMRI_gradient_0000:= 1 -1 0
  // DWMRI_gradient_0001:= 0 0 0
  //           ...
  */ 
  itk::MetaDataDictionary imgMetaDictionary = img->GetMetaDataDictionary();    
  std::vector<std::string> imgMetaKeys = imgMetaDictionary.GetKeys();
  std::vector<std::string>::const_iterator itKey = imgMetaKeys.begin();
  std::string metaString;
  
  TensorReconstructionImageFilterType::GradientDirectionType vect3d;
  TensorReconstructionImageFilterType::GradientDirectionContainerType::Pointer 
    DiffusionVectors = 
    TensorReconstructionImageFilterType::GradientDirectionContainerType::New();
 

  for (; itKey != imgMetaKeys.end(); itKey ++)
    {
    double x,y,z;

    itk::ExposeMetaData<std::string> (imgMetaDictionary, *itKey, metaString);
    if (itKey->find("DWMRI_gradient") != std::string::npos)
      { 
      std::cout << *itKey << " ---> " << metaString << std::endl;      
      sscanf(metaString.c_str(), "%lf %lf %lf\n", &x, &y, &z);
      vect3d[0] = x; vect3d[1] = y; vect3d[2] = z;
      
      DiffusionVectors->InsertElement( numberOfImages, vect3d );
      ++numberOfImages;
      }
    else if (itKey->find("DWMRI_b-value") != std::string::npos)
      {
      std::cout << *itKey << " ---> " << metaString << std::endl;      
      readb0 = true;
      b0 = atof(metaString.c_str());
      }
    }
  std::cout << "Number Of Gradient images: " << numberOfImages << std::endl;
  if(!readb0)
    {
    std::cerr << "BValue not specified in header file" << std::endl;
    return EXIT_FAILURE;
    }

  TensorReconstructionImageFilterType::Pointer tensorReconstructionFilter = 
    TensorReconstructionImageFilterType::New();
  
  tensorReconstructionFilter->SetGradientImage( DiffusionVectors, reader->GetOutput() );

  // This is necessary until netlib/dsvdc.c becomes thread safe
  tensorReconstructionFilter->SetNumberOfThreads( 1 );
 
  tensorReconstructionFilter->SetBValue(b0);  
  tensorReconstructionFilter->SetThreshold( static_cast< 
      TensorReconstructionImageFilterType::ReferencePixelType >( 
                                                    atof(argv[2])));
  tensorReconstructionFilter->Update();

  // Manipulate the Meta Data so that the file written on disk has it correct
  // Ensure the measurement frame is transferred to the output as well.
  itk::MetaDataDictionary &inDic = reader->GetOutput()->GetMetaDataDictionary();
  std::vector<std::vector<double> > msrFrame(ImageDimension);
  for (unsigned int j = 0; j < ImageDimension; j++) {
    msrFrame[j].resize(ImageDimension);
  }
  itk::ExposeMetaData(inDic, "NRRD_measurement frame", msrFrame);

  itk::MetaDataDictionary &outDic = 
             tensorReconstructionFilter->GetOutput()->GetMetaDataDictionary();
  itk::EncapsulateMetaData<std::vector<std::vector<double> > >(outDic,
     std::string("NRRD_measurement frame"), msrFrame);

  typedef itk::ImageFileWriter< 
    TensorReconstructionImageFilterType::OutputImageType > TensorWriterType;
  TensorWriterType::Pointer tensorWriter = TensorWriterType::New();
  tensorWriter->SetFileName( argv[3] );
  tensorWriter->SetInput( tensorReconstructionFilter->GetOutput() );
  try {
    tensorWriter->UseCompressionOff();
    tensorWriter->Update();
  } catch  (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      std::cerr << "Failed trying to write to " << argv[3] << std::endl;
      exit(1);
  }


  return EXIT_SUCCESS;
}
