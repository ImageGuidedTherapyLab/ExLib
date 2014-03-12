#include <itkMaskImageFilter.h>
#include <itkVector.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

using namespace std;

int main(int argc, char **argv)
{
  if (argc < 4)
  {
    cout << "Usage: " << argv[0] << " inputFileName maskFileName outputFileName" << endl;
    return 1;
  }

  string refName(argv[1]), maskName(argv[2]), resName(argv[3]);
  
  typedef itk::Image <float,3> FloatImageType;
  typedef itk::ImageFileReader <FloatImageType> itkFloatReader;
  typedef itk::ImageFileWriter <FloatImageType> itkFloatWriter;
  
  typedef itk::Image <itk::Vector <float,6>,3> VectorFloatImageType;
  typedef itk::ImageFileReader <VectorFloatImageType> itkVectorFloatReader;
  typedef itk::ImageFileWriter <VectorFloatImageType> itkVectorFloatWriter;
  
  typedef itk::Image <unsigned short,3> UShortImageType;
  typedef itk::ImageFileReader <UShortImageType> itkUShortReader;
  
  typedef itk::MaskImageFilter <FloatImageType,UShortImageType> itkMaskFilterType;
  typedef itk::MaskImageFilter <VectorFloatImageType,UShortImageType> itkVectorMaskFilterType;

  itkUShortReader::Pointer maskInput = itkUShortReader::New();
  maskInput->SetFileName(maskName.c_str());
  maskInput->Update();
  
  itkFloatReader::Pointer preInput = itkFloatReader::New();
  preInput->SetFileName(refName.c_str());
  preInput->Update();
  
  preInput->GetImageIO()->ReadImageInformation();
  
  bool isImageVector = (preInput->GetImageIO()->GetNumberOfComponents() > 1);
  
  if (isImageVector)
  {
    itkVectorFloatReader::Pointer vecInput = itkVectorFloatReader::New();
    vecInput->SetFileName(refName.c_str());
    vecInput->Update();
    
    itkVectorMaskFilterType::Pointer imageMasker = itkVectorMaskFilterType::New();
    imageMasker->SetInput(vecInput->GetOutput());
    imageMasker->SetInput2(maskInput->GetOutput());
    imageMasker->Update();
    
    itkVectorFloatWriter::Pointer vecOutput = itkVectorFloatWriter::New();
    vecOutput->SetFileName(resName.c_str());
    vecOutput->SetInput(imageMasker->GetOutput());
    vecOutput->SetUseCompression(true);
    
    vecOutput->Update();
  }
  else
  {
    itkMaskFilterType::Pointer imageMasker = itkMaskFilterType::New();
    imageMasker->SetInput(preInput->GetOutput());
    imageMasker->SetInput2(maskInput->GetOutput());
    imageMasker->Update();
    
    itkFloatWriter::Pointer scalarOutput = itkFloatWriter::New();
    scalarOutput->SetFileName(resName.c_str());
    scalarOutput->SetInput(imageMasker->GetOutput());
    scalarOutput->SetUseCompression(true);
    
    scalarOutput->Update();
  }
  
  return 0;
}
