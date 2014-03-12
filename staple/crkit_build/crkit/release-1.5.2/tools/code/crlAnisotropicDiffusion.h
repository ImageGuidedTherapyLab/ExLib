
#ifndef _ANDIFF_INCLUDED
#define _ANDIFF_INCLUDED 1

#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkResampleImageFilter.h>

#include "itkCurvatureAnisotropicDiffusionImageFilter.h"

#include <iostream>
#include <fstream>
#include <sstream>


class AnisotropicDiffusion
{

  static int const ImageDimension = 3;
  typedef   float PixelType;

  typedef   itk::OrientedImage<PixelType,ImageDimension>  ImageType;
  typedef   itk::ImageFileReader< ImageType >    ImageReaderType;
  typedef   itk::ImageFileWriter< ImageType >    ImageWriterType;
  typedef   itk::CurvatureAnisotropicDiffusionImageFilter< ImageType, ImageType>  DiffusionFilterType;

public:
  AnisotropicDiffusion() {
    // Default parameter values are reasonable for MRI
    m_Diffusion = DiffusionFilterType::New();
    m_ImageReader = ImageReaderType::New();
    m_ImageWriter = ImageWriterType::New();
    m_Diffusion->SetInput(m_ImageReader->GetOutput());
    m_ImageWriter->SetInput(m_Diffusion->GetOutput());

    m_Diffusion->SetNumberOfIterations( 5 );
    m_Diffusion->SetConductanceParameter( 7.0 );
  }

  ~AnisotropicDiffusion() { } ;

  void LoadImage( const char * filename )
  {
    if (!filename) return;
    m_ImageReader->SetFileName(filename);
    m_ImageReader->Update();
  };

  void WriteImage( const char * filename )
  {
    if (!filename) return;

    m_ImageWriter->SetFileName(filename);
    try {
      m_ImageWriter->Update();
      m_ImageWriter->UseCompressionOn();
      m_ImageWriter->Write();
    } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
    }
  }

  void SmoothImage( ) {
    m_Diffusion->Update();
  };

  void SetIterations( int iterations ) {
    m_Diffusion->SetNumberOfIterations(iterations);
  };

  void SetConductance( float conductance ) {
    m_Diffusion->SetConductanceParameter(conductance);
  };

  int GetIterations( ) {
    return m_Diffusion->GetNumberOfIterations();
  };

  float GetConductance( ) {
    return m_Diffusion->GetConductanceParameter();
  };

  void Execute(char *inname, char *outname) {
    LoadImage( inname );
    SmoothImage( );
    WriteImage( outname );
  };

protected:

  ImageReaderType::Pointer           m_ImageReader;
  ImageWriterType::Pointer           m_ImageWriter;
  DiffusionFilterType::Pointer       m_Diffusion;

};

#endif 
