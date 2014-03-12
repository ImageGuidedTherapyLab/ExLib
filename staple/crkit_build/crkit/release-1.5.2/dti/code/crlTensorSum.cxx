
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSymmetricSecondRankTensor.h>
#include <itkDiffusionTensor3D.h>
#include <itkMatrix.h>
#include <itkUnaryFunctorImageFilter.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <typeinfo>

#include "TensorSumFunctor.h"

int main(int argc, const char *argv[])
{

  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " intensor1 intensor2 outtensor" << std::endl;
    std::cerr << "This program computes the sum of two tensor images." << std::endl;
    exit(1);
  }

  static unsigned int const ImageDimension = 3;
  typedef float ComponentType;
  typedef   itk::DiffusionTensor3D<ComponentType>  PixelType;
  typedef   itk::Image< PixelType, ImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  typedef   itk::ImageFileWriter<ImageType> ImageWriterType;
  
  ImageReaderType::Pointer  reader1 = ImageReaderType::New();
  reader1->SetFileName( argv[1] );

  try {
    reader1->GenerateOutputInformation();
  } catch (itk::ExceptionObject &excp) 
  {
    std::cerr << "Failed to read from file " << argv[1] << std::endl;
    std::cerr << "Exception error message is : " << excp << std::endl;
    return EXIT_FAILURE;
  }

  itk::ImageIOBase *imageIO = reader1->GetImageIO();
  if (imageIO->GetComponentType() != itk::ImageIOBase::FLOAT) {
    std::cerr << "File component type does not match expected type FLOAT."
              << std::endl;
    std::cerr << "Use crlCastSymMatToFloat to change the type of the components." 
              << std::endl;
    return EXIT_FAILURE;
  }

  reader1->Update();

  ImageReaderType::Pointer  reader2 = ImageReaderType::New();
  reader2->SetFileName( argv[2] );

  try {
    reader2->GenerateOutputInformation();
  } catch (itk::ExceptionObject &excp) 
  {
    std::cerr << "Failed to read from file " << argv[2] << std::endl;
    std::cerr << "Exception error message is : " << excp << std::endl;
    return EXIT_FAILURE;
  }

  imageIO = reader2->GetImageIO();
  if (imageIO->GetComponentType() != itk::ImageIOBase::FLOAT) {
    std::cerr << "File component type does not match expected type FLOAT."
              << std::endl;
    std::cerr << "Use crlCastSymMatToFloat to change the type of the components." 
              << std::endl;
    return EXIT_FAILURE;
  }

  reader2->Update();

  typedef itk::TensorSumImageFilter<ImageType, ImageType, ImageType> TensorSumImageFilterType;
  TensorSumImageFilterType::Pointer tensorSumFilter = 
                                            TensorSumImageFilterType::New();
  tensorSumFilter->SetInput1( reader1->GetOutput() );
  tensorSumFilter->SetInput2( reader2->GetOutput() );
  tensorSumFilter->Update();

  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput ( tensorSumFilter->GetOutput() );
  writer->SetFileName( argv[3] );

  // Transfer the measurement frame from the input data to the output data
  // Ensure the measurement frame is transferred to the output as well.
  itk::MetaDataDictionary &inDic = reader1->GetOutput()->GetMetaDataDictionary();  std::vector<std::vector<double> > msrFrame(ImageDimension);
  for (unsigned int j = 0; j < ImageDimension; j++) {
    msrFrame[j].resize(ImageDimension);
  }
  itk::ExposeMetaData(inDic, "NRRD_measurement frame", msrFrame);

  itk::MetaDataDictionary &outDic =
             tensorSumFilter->GetOutput()->GetMetaDataDictionary();
  itk::EncapsulateMetaData<std::vector<std::vector<double> > >(outDic,
     std::string("NRRD_measurement frame"), msrFrame);

  try {
    // Compression fails with NRRD format and large images.
    writer->UseCompressionOff();
    writer->Write(); // 
  } catch ( itk::ExceptionObject & excp ) {
    std::cerr << "Failed to write out the image." << std::endl;
    std::cerr << "Exception error message is : " << excp << std::endl;
  }

  exit(0); // success
}

