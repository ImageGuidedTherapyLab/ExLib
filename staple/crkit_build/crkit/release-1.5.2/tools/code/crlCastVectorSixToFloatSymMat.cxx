
#include <itkImage.h>
#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkSymmetricSecondRankTensor.h>
#include <itkImageRegionIterator.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkVector.h>

#include <iostream>

int main(int argc, char *argv[])
{

  static int const ImageDimension = 3;
  typedef   itk::SymmetricSecondRankTensor< float >    FloatTensorPixelType;
  typedef   itk::SymmetricSecondRankTensor< double >    DoubleTensorPixelType;

  typedef itk::Vector<float,6>   FloatVectorPixelType;
  typedef itk::Vector<double,6>  DoubleVectorPixelType;

 // typedef   itk::Image< FloatVectorPixelType, ImageDimension > ImageTypeIn;
  typedef   itk::Image< FloatTensorPixelType, ImageDimension > ImageTypeOut;
  typedef   itk::Image< FloatTensorPixelType, ImageDimension > ImageTypeIn;
  typedef   itk::ImageFileReader< ImageTypeIn >    ImageReaderType;
  typedef   itk::ImageFileWriter< ImageTypeOut >    ImageWriterType;

  if (argc != 3) {
    std::cout << 
      "Usage: " << argv[0] << " inimage outimage " << std::endl;
    std::cout << "This purpose of this program is to take a " 
     << "3D second order symmetric tensor with components stored as float "
     << " in a vector, and modify it to have the components stored in a "
     << "symmetric matrix file type." << std::endl;
    exit(1);
  }

  ImageReaderType::Pointer reader;
  ImageWriterType::Pointer writer;

  reader = ImageReaderType::New();
  reader->SetFileName(argv[1]);

  try {
    reader->Update();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  writer = ImageWriterType::New();
  writer->SetFileName(argv[2]);
  writer->SetInput( reader->GetOutput() );

/*
  // Ensure the measurement frame is transferred to the output as well.
  itk::MetaDataDictionary &inDic = reader->GetOutput()->GetMetaDataDictionary();
  std::vector<std::vector<double> > msrFrame(ImageDimension);
  for (unsigned int j = 0; j < ImageDimension; j++) {
    msrFrame[j].resize(ImageDimension);
  }
  itk::ExposeMetaData(inDic, "NRRD_measurement frame", msrFrame);

  itk::MetaDataDictionary &outDic = outImage->GetMetaDataDictionary();
  itk::EncapsulateMetaData<std::vector<std::vector<double> > >(outDic,
     std::string("NRRD_measurement frame"), msrFrame);
*/

  try {
    writer->Update();
    // writer->UseCompressionOn();
    writer->Write();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  exit(0); // success
}

