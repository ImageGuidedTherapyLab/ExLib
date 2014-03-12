/*
 * Copyright 2010 Simon Warfield simon.warfield@childrens.harvard.edu
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 * http://www.crl.med.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#ifndef _CRLDENSITYWEIGHTEDSTATS_H
#define _CRLDENSITYWEIGHTEDSTATS_H

#include <itkImage.h>
#include <itkImageFileReader.h>

#include <itkImageRegionIteratorWithIndex.h>

namespace crl
{

class ITK_EXPORT DensityWeightedStatsBase : public itk::LightObject
{
  // Defines the public API for running crlDensityWeightedStats .
  // The templated arguments are provided by classes that inherit this API.
  //  This allows us to create a pointer to a class of this type, and
  // then to change it at run time to match the types we actually need to
  // deal with.

  public: 
  typedef DensityWeightedStatsBase Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(DensityWeightedStatsBase, LightObject);

  DensityWeightedStatsBase() { };
  ~DensityWeightedStatsBase() { };

  void SetDensityImageFileName(std::string name) { 
         m_DensityImageFileName = name; 
  }
  void SetScalarImageFileName(std::string name) { 
         m_ScalarImageFileName = name; 
  }

  virtual int Execute() = 0;

  protected:
  std::string m_DensityImageFileName;
  std::string m_ScalarImageFileName;

private:
  DensityWeightedStatsBase(const DensityWeightedStatsBase &); // purposely not implemented
  void operator=(const DensityWeightedStatsBase &); // purposely not implemented.

};

template <const unsigned int TImageDimension, class TPixelType >
class ITK_EXPORT DensityWeightedStats : public DensityWeightedStatsBase
{
public:
  /** Standard class typedefs. */
  typedef   DensityWeightedStats Self;
  typedef   DensityWeightedStatsBase Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef   itk::Image< TPixelType, TImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  DensityWeightedStats() { };

  ~DensityWeightedStats() {};

  int Execute();

private:
  DensityWeightedStats(const Self&); // purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlDensityWeightedStats.txx"
#endif

} // end namespace crl

#endif

