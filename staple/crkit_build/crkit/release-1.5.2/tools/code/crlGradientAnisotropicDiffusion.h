#ifndef CRLGRADIENTANISOTROPICDIFFUSION_H_
#define CRLGRADIENTANISOTROPICDIFFUSION_H_

#include "itkMacro.h"
#include "itkImage.h"
#include "itkCovariantVector.h"
#include "itkVectorGradientAnisotropicDiffusionImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorImageToImageAdaptor.h"

namespace crl 
{


class GradientAnisotropicDiffusion 
{
public:
  // use the following static function (defined at the end of the file) in order to interrogate the input image and allocate the 
  // proper templated algorithm based on the vector size and image dimension.  You can then use the pointer to the base class generically.
  // 
  // So you do something like this:
  // crl::GradientAnisotropicDiffusion* diff = crl::GradientAnisotropicDiffusion::ProbeFile( inputfile );
  // diff->ReadFile(inputfile);
  // diff->Execute();
  // diff->WriteFile(outputfile);
  // delete diff;
  static GradientAnisotropicDiffusion* ProbeFile( const std::string& inputfile );

    
  virtual void ReadFile(const std::string&) = 0;
  virtual void Execute() = 0;
  virtual void WriteFile(const std::string&) = 0;
  virtual void SetIterations(unsigned int& i)  = 0;
  virtual void SetConductanceParameter( double& c ) = 0;
  virtual void SetTimeStep( double& t) = 0;
  virtual void SetThreads(unsigned int&) = 0;

  virtual ~GradientAnisotropicDiffusion() {};

};

template <const unsigned int VectorDimension = 2, const unsigned int ImageDimension = 3>
class GradientAnisotropicDiffusionSpecial : public GradientAnisotropicDiffusion
{
public:
  typedef float ElementType;
  typedef itk::CovariantVector<ElementType, VectorDimension > VectorType;
  typedef typename itk::Image<VectorType, ImageDimension >	ImageType;
  
  // recommended timestep <= ( 0.5 / 2.0^ImageDimension )
  // for 2D: 0.125
  // for 3D: 0.0625
  GradientAnisotropicDiffusionSpecial() : m_InputImage(0),m_OutputImage(0),m_Iterations(10),m_Threads(1),m_ConductanceParameter(1.0), m_TimeStep(0.05) {};
  

  void ReadFile(const std::string& inputfile)
    {
    typedef itk::ImageFileReader<ImageType > ImageReaderType;
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(inputfile);
    reader->Update();
    m_InputImage = reader->GetOutput();
    m_InputImage->DisconnectPipeline();	
    };
    
  void Execute()
    {
    typedef itk::VectorGradientAnisotropicDiffusionImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(m_InputImage);
    filter->SetConductanceParameter(m_ConductanceParameter);
    filter->SetNumberOfIterations(m_Iterations);
    filter->SetTimeStep(m_TimeStep);
    filter->SetNumberOfThreads(m_Threads);
    filter->Update();
    m_OutputImage = filter->GetOutput();
    m_OutputImage->DisconnectPipeline();
    };
    
  void WriteFile(const std::string& outputfile)
    {
    if ( VectorDimension != 1 )
      {
      typedef itk::ImageFileWriter<ImageType > ImageWriterType;
      typename ImageWriterType::Pointer writer = ImageWriterType::New();
      writer->SetInput( m_OutputImage );
      writer->SetFileName( outputfile );
      writer->Update();
      }
    else
      {
      typedef itk::Image<ElementType, ImageDimension > ScalarImageType;
      typedef itk::ImageFileWriter< ScalarImageType > ImageWriterType;
      typedef itk::ImageRegionConstIterator< ImageType > VectorIteratorType;
      typedef itk::ImageRegionIterator<ScalarImageType > ScalarIteratorType;
      
      typename ScalarImageType::Pointer scalar = ScalarImageType::New();
      scalar->CopyInformation( m_OutputImage );
      scalar->SetRegions( m_OutputImage->GetLargestPossibleRegion() );
      scalar->Allocate();
      
      VectorIteratorType viter(m_OutputImage, m_OutputImage->GetBufferedRegion() );
      ScalarIteratorType siter(scalar, scalar->GetBufferedRegion() );
      
      viter.GoToBegin();
      siter.GoToBegin();
      while ( !viter.IsAtEnd() && !siter.IsAtEnd() )
	{
	siter.Set(viter.Value()[0]);
	++siter;
	++viter;
	}
      
      typename ImageWriterType::Pointer writer = ImageWriterType::New();
      writer->SetInput( scalar );
      writer->SetFileName(outputfile);
      writer->Update();
      }
    };
  
  void SetIterations(unsigned int& i) { m_Iterations = i; };
  void SetConductanceParameter( double& c ) { m_ConductanceParameter = c;  };
  void SetTimeStep( double& t) { m_TimeStep = t; };
  void SetThreads(unsigned int& i) {m_Threads=i;};
  
  void GetVectorDimension() { return VectorDimension; };
  void GetImageDimension() { return ImageDimension; };

  virtual ~GradientAnisotropicDiffusionSpecial() {};
  
private:
  typename ImageType::Pointer m_InputImage;
  typename ImageType::Pointer m_OutputImage;
  unsigned int m_Iterations;
  unsigned int m_Threads;
  double m_ConductanceParameter;
  double m_TimeStep;


};

GradientAnisotropicDiffusion* GradientAnisotropicDiffusion::ProbeFile(const std::string& inputfile)
  {
  typedef itk::ImageFileReader< itk::VectorImage<float> > ImageReaderType;
  ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName(inputfile);
  reader->GenerateOutputInformation();
  unsigned int imageDimension = reader->GetImageIO()->GetNumberOfDimensions();
  unsigned int vectorDimension = reader->GetImageIO()->GetNumberOfComponents();
  
  std::cout << "Image Dimension = " << imageDimension << ", Vector Size = " << vectorDimension << std::endl;
  
  if      ( vectorDimension == 1 && imageDimension == 2)
    return new GradientAnisotropicDiffusionSpecial<1,2>();
  else if ( vectorDimension == 2 && imageDimension == 2)
    return new GradientAnisotropicDiffusionSpecial<2,2>();
  else if ( vectorDimension == 3 && imageDimension == 2)
    return new GradientAnisotropicDiffusionSpecial<3,2>();
  else if ( vectorDimension == 1 && imageDimension == 3 )
    return new GradientAnisotropicDiffusionSpecial<1,3>();
  else if ( vectorDimension == 2 && imageDimension == 3)
    return new GradientAnisotropicDiffusionSpecial<2,3>();
  else if ( vectorDimension == 3 && imageDimension == 3)
    return new GradientAnisotropicDiffusionSpecial<3,3>();
  else
    throw std::runtime_error("we only handle 2D or 3D images with 1 or 2 vector components");
  
  return 0;	// NEVER REACHED
  }

};


#endif /*CRLGRADIENTANISOTROPICDIFFUSION_H_*/

