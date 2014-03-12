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

#ifndef _CRL_DEFORM_SCALAR_IMAGE_H
#define _CRL_DEFORM_SCALAR_IMAGE_H 1

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkOrientedImage.h>

#include <itkWarpImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkWindowedSincInterpolateImageFunction.h>

class ITK_EXPORT crlDeformScalarImageBase : public itk::LightObject
{
  // Define public API for deformation of a scalar image.
  // Templated arguments are added by derived classes that inherit this API.
  //

  // The deformation field is represented as an image of vector pixel types. The
  // dimension of the vectors is the same as the dimension of the input image.
  // Each vector in the deformation field represents the distance between a 
  // geometric point in the input space and a point in the output space 
  // such that:
  // p_{in} = p_{out} + distance

public:
  typedef crlDeformScalarImageBase Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(crlDeformScalarImageBase, LightObject);

  crlDeformScalarImageBase() {};
  ~crlDeformScalarImageBase() {};

  void SetDeformationFieldFileName(std::string name) {
    m_DeformationFieldFileName = name;
  };

  void SetInputFileName(std::string name) {
    m_InputFileName = name;
  };

  void SetOutputFileName(std::string name) {
    m_OutputFileName = name;
  };

  void SetInterpolationType(std::string name) {
    m_InterpolationType = name;
  };

  virtual int Execute() = 0;
 
protected:
  std::string m_DeformationFieldFileName;
  std::string m_InputFileName;
  std::string m_OutputFileName;
  std::string m_InterpolationType;

private:
  crlDeformScalarImageBase(const crlDeformScalarImageBase&); // purposely not implemented 
  void operator=(const crlDeformScalarImageBase&); // purposely not implemented

};


template <const unsigned int TImageDimension, class TPixelType >
class ITK_EXPORT crlDeformScalarImage : public crlDeformScalarImageBase
{
public:
  /** Standard class typedefs. */
  typedef crlDeformScalarImage Self;
  typedef crlDeformScalarImageBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   TPixelType VectorComponentType;
  typedef   itk::Vector< VectorComponentType, TImageDimension > VectorPixelType;
  typedef   itk::OrientedImage< VectorPixelType, TImageDimension > DeformationFieldType;
  
  typedef   itk::OrientedImage< TPixelType,  TImageDimension >   ImageType;

  typedef   itk::ImageFileReader< ImageType >  ReaderType;
  typedef   itk::ImageFileWriter< ImageType >  WriterType;
  
  // A reader for the deformation field:
  typedef   itk::ImageFileReader< DeformationFieldType >  FieldReaderType;

  // A deformation field warping filter
  typedef itk::WarpImageFilter< ImageType, 
                                ImageType, 
                                DeformationFieldType  >  FilterType;

  // Values are interpolated to the non-integer locations of the warped points.
  typedef itk::LinearInterpolateImageFunction< 
                       ImageType, double >  LInterpolatorType;
  typedef itk::NearestNeighborInterpolateImageFunction< ImageType,
          double> NInterpolatorType;
  typedef itk::BSplineInterpolateImageFunction<
                       ImageType, double >  BSplineInterpolatorType;
  typedef itk::ConstantBoundaryCondition< ImageType >  BoundaryConditionType;
  static const unsigned int WindowRadius = 5;
  typedef itk::Function::HammingWindowFunction<WindowRadius>  WindowFunctionType;

  typedef itk::WindowedSincInterpolateImageFunction<
                                          ImageType,
                                          WindowRadius,
                                          WindowFunctionType,
                                          BoundaryConditionType,
                                          double  >    SincInterpolatorType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  crlDeformScalarImage() {
    m_InterpolationType = std::string("linear");
  };

  ~crlDeformScalarImage() {};

  int Execute();

private:
  crlDeformScalarImage(const Self&); // purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlDeformScalarImage.txx"
#endif

#endif
