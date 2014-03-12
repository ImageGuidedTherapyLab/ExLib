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



/* Raise the dimensionality by 1 and insert the vector components as the
 * elements of the Nth dimension.
 *
 * WARNING: The geometry of the Nth dimension is not set.
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

  static int const ImageDimensionIn = 2;
  static int const ImageDimensionOut = 3;
  typedef   float PixelType;

  typedef   itk::VectorImage< PixelType , ImageDimensionIn>  ImageTypeIn;
  typedef   itk::Image<PixelType,ImageDimensionOut>  ImageTypeOut;
  typedef   itk::ImageFileReader< ImageTypeIn,
               itk::DefaultConvertPixelTraits<PixelType> > ImageReaderType;

  typedef   itk::ImageFileWriter< ImageTypeOut >    ImageWriterType;

  typedef itk::ImageRegionIteratorWithIndex< ImageTypeIn > InIteratorType;
  typedef itk::ImageRegionIteratorWithIndex< ImageTypeOut > OutIteratorType;

  if (argc != 3) {
    std::cout <<
      "Usage: " << argv[0] << " inimage outimage" << std::endl;
    std::cout << "Convert vector image to 4 dimensional image." << std::endl;
    exit(1);
  }

  ImageReaderType::Pointer r;
  ImageWriterType::Pointer w;

  r = ImageReaderType::New();
  r->SetFileName(argv[1]);
  r->Update();

  ImageTypeOut::Pointer outImage = ImageTypeOut::New();
  ImageTypeOut::RegionType region;

  ImageTypeIn::SpacingType inSpacing;
  ImageTypeOut::SpacingType outSpacing;
  ImageTypeIn::PointType inOrigin;
  ImageTypeOut::PointType outOrigin;
  ImageTypeIn::DirectionType inDirection;
  ImageTypeOut::DirectionType outDirection;

  inSpacing = r->GetOutput()->GetSpacing();
  inOrigin  = r->GetOutput()->GetOrigin();
  inDirection = r->GetOutput()->GetDirection();
  for (unsigned int i = 0; i < region.GetImageDimension(); i++) {
    region.SetSize(i, r->GetOutput()->GetLargestPossibleRegion().GetSize(i));
    region.SetIndex(i, r->GetOutput()->GetLargestPossibleRegion().GetIndex(i));
    outSpacing[i] = inSpacing[i];
    outOrigin[i] = inOrigin[i];
    for (unsigned int j = 0; j < region.GetImageDimension(); j++)
      {
	outDirection[i][j] = inDirection[i][j];
      }
  }
  region.SetSize(region.GetImageDimension()-1,
                         r->GetOutput()->GetVectorLength());
  region.SetIndex(region.GetImageDimension()-1, 0);
  outImage->SetRegions(region);
  outSpacing[region.GetImageDimension()-1] = 0.0;
  outOrigin[region.GetImageDimension()-1] = 0.0;
  for (unsigned int j = 0; j < region.GetImageDimension()-1; j++)
    {
      outDirection[region.GetImageDimension()-1][j] = 0.0;
      outDirection[j][region.GetImageDimension()-1] = 0.0;
    }
  outDirection[region.GetImageDimension()-1][region.GetImageDimension()-1] = 1.0;

  outImage->SetSpacing( outSpacing );
  outImage->SetOrigin( outOrigin );
  outImage->SetDirection( outDirection );
  outImage->Allocate();

  InIteratorType inItr(r->GetOutput(),
                       r->GetOutput()->GetLargestPossibleRegion());
  for (inItr.GoToBegin(); !inItr.IsAtEnd(); ++inItr ) {
    ImageTypeIn::PixelType pin = inItr.Get();
    ImageTypeIn::IndexType inIndex = inItr.GetIndex();
    ImageTypeOut::IndexType outIndex;
    for (signed int c = 0; c < ImageDimensionIn; c++) {
      outIndex[c] = inIndex[c];
    }
    for (unsigned int v = 0; v < pin.Size(); v++) {
      outIndex[ImageDimensionIn] = v;
      outImage->SetPixel(outIndex, pin[v]);
    }
  }

  w = ImageWriterType::New();
  w->SetFileName(argv[2]);
  w->UseCompressionOn();
  w->SetInput( outImage );

  try {
    w->Update();
    w->Write();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  exit(0); // success
}
