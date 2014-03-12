
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSymmetricSecondRankTensor.h>
#include <itkDiffusionTensor3D.h>
#include <itkMatrix.h>
#include <itkUnaryFunctorImageFilter.h>
#include <itkMetaDataObject.h>
#include <itkMetaDataDictionary.h>

#include "TensorExpFunctor.h"

int main(int argc, const char *argv[])
{

  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " intensor outtensor" << std::endl;
    std::cerr << "This program computes the tensor exponential." << std::endl;
    exit(1);
  }

  static int const ImageDimension = 3;
  typedef float ComponentType;
  typedef   itk::DiffusionTensor3D<ComponentType>  PixelType;
  typedef   itk::Image< PixelType, ImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  typedef   itk::ImageFileWriter<ImageType> ImageWriterType;
  
  ImageReaderType::Pointer  reader = ImageReaderType::New();
  reader->SetFileName( argv[1] );

  try {
    reader->GenerateOutputInformation();
  } catch (itk::ExceptionObject &excp) 
  {
    std::cerr << "Failed to read from file " << argv[1] << std::endl;
    std::cerr << "Exception error message is : " << excp << std::endl;
    return EXIT_FAILURE;
  }

  itk::ImageIOBase *imageIO = reader->GetImageIO();
  if (imageIO->GetComponentType() != itk::ImageIOBase::FLOAT) {
    std::cerr << "File component type does not match expected type FLOAT."
              << std::endl;
    std::cerr << "Use crlCastSymMatToFloat to change the type of the components." 
              << std::endl;
    return EXIT_FAILURE;
  }

  reader->Update();

  typedef itk::TensorExpImageFilter<ImageType, ImageType> TensorExpImageFilterType;
  TensorExpImageFilterType::Pointer tensorExpFilter = 
                                            TensorExpImageFilterType::New();
  tensorExpFilter->SetInput( reader->GetOutput() );
  tensorExpFilter->Update();

  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput ( tensorExpFilter->GetOutput() );
  writer->SetFileName( argv[2] );

  // Ensure the measurement frame is transferred to the output as well.
  itk::MetaDataDictionary &inDic = reader->GetOutput()->GetMetaDataDictionary();
  std::vector<std::vector<double> > msrFrame(ImageDimension);
  for (unsigned int j = 0; j < ImageDimension; j++) {
    msrFrame[j].resize(ImageDimension);
  }
  itk::ExposeMetaData(inDic, "NRRD_measurement frame", msrFrame);

  itk::MetaDataDictionary &outDic = 
               tensorExpFilter->GetOutput()->GetMetaDataDictionary();
  itk::EncapsulateMetaData<std::vector<std::vector<double> > >(outDic,
     std::string("NRRD_measurement frame"), msrFrame);

  try {
    // Write compression fails with large diffusion tensor images.
    writer->UseCompressionOff();
    writer->Write();
  } catch ( itk::ExceptionObject & excp ) {
    std::cerr << "Failed to write out the image." << std::endl;
    std::cerr << "Exception error message is : " << excp << std::endl;
  }

  exit(0); // success
}

