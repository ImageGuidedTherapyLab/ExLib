// Type definitions used throughout the display software.

#ifndef _ImageTypeDefinitions_Included
#define _ImageTypeDefinitions_Included 1

#include <itkImage.h>
#include <itkOrientedImage.h>
#include <itkOrientImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageToVTKImageFilter.h>
#include <itkVTKImageToImageFilter.h>


class ImageTypeDefinitions
{
public:
  static const unsigned int Dimension = 3;
  typedef float InternalPixelType;
  typedef itk::OrientedImage<InternalPixelType, Dimension> ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;
  typedef itk::ImageToVTKImageFilter< ImageType > ITK2VTKConnectorFilterType;
  typedef itk::VTKImageToImageFilter< ImageType > VTK2ITKConnectorFilterType;

  enum DataType {greyscale,probmap,segmentation,fmri,tensor,rgb,polydata,scalarpolydata};
  enum OrientationType {axial=0, sagittal, coronal};
};

#endif
