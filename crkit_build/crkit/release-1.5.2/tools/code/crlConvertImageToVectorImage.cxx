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


/* Reduce the dimensionality by 1 and make the most major axis a vector index
 *
 */


#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageIterator.h>
#include <itkVectorImage.h>
#include <itkImageAdaptor.h>
#include <itkVectorImageToImageAdaptor.h>

#include <iostream>

int main(int argc, char *argv[])
{

  static int const ImageDimensionIn = 4;
  static int const ImageDimensionOut = 3;
  typedef   float PixelType;

  typedef   itk::VectorImage< PixelType , ImageDimensionOut>  ImageTypeOut;
  typedef   itk::Image<PixelType,ImageDimensionIn>  ImageTypeIn;
  typedef   itk::ImageFileReader< ImageTypeIn,  
               itk::DefaultConvertPixelTraits<PixelType> > ImageReaderType;

  typedef   itk::ImageFileWriter< ImageTypeOut >    ImageWriterType;

  typedef itk::ImageRegionConstIteratorWithIndex< ImageTypeIn > InIteratorType;
  typedef itk::ImageRegionIteratorWithIndex< ImageTypeOut > OutIteratorType;

  if (argc != 3) {
    std::cout << 
      "Usage: " << argv[0] << " inimage outimage" << std::endl;
    std::cout << "Convert 4 dimensional image to a vector image." << std::endl;
    exit(1);
  }

  ImageReaderType::Pointer r;
  ImageWriterType::Pointer w;

  r = ImageReaderType::New();
  r->SetFileName(argv[1]);
  r->Update();
  
  ImageTypeIn::Pointer inImage = r->GetOutput();
  inImage->DisconnectPipeline();	// for fun 

  ImageTypeOut::Pointer outImage = ImageTypeOut::New();
  ImageTypeOut::RegionType outRegion;

  ImageTypeOut::SpacingType outSpacing;
  ImageTypeOut::PointType outOrigin;
  
  for (signed int i = 0; i < ImageDimensionOut; i++)
    {
    outRegion.SetSize(i, inImage->GetLargestPossibleRegion().GetSize(i));
    outRegion.SetIndex(i, inImage->GetLargestPossibleRegion().GetIndex(i));
    outSpacing[i] = inImage->GetSpacing()[i];
    outOrigin[i] = inImage->GetOrigin()[i];
    }
  
  
  outImage->SetRegions(outRegion);
  outImage->SetSpacing( outSpacing );
  outImage->SetOrigin( outOrigin );
  // size of last dimension of input image is the vector length

  unsigned int ncomponents = inImage->GetLargestPossibleRegion().GetSize()[ ImageDimensionIn - 1 ];
  std::cout << "making a " << ncomponents << " vector image." << std::endl;
  outImage->SetNumberOfComponentsPerPixel( ncomponents );
  outImage->Allocate();
  
  
  ImageTypeIn::RegionType inRegion = inImage->GetLargestPossibleRegion();
  
  for (unsigned int i = 0; i < ncomponents; i++ )
    {
    OutIteratorType outItr(outImage, outImage->GetLargestPossibleRegion());
    inRegion.SetSize(ImageDimensionIn-1, 1);
    inRegion.SetIndex(ImageDimensionIn-1, i);
    InIteratorType inItr(inImage, inRegion);
    
    for ( outItr.GoToBegin(), inItr.GoToBegin(); !inItr.IsAtEnd() && !outItr.IsAtEnd(); ++outItr, ++inItr)
      {
      itk::VariableLengthVector<double> pixel = outItr.Get();
      pixel[i] = inItr.Get();
      outItr.Set(pixel);
      }
    }

  w = ImageWriterType::New();
  w->SetFileName(argv[2]);
  w->SetInput( outImage );
  w->UseCompressionOn();

  try {
    w->Update();
    w->Write();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  exit(0); // success
}

