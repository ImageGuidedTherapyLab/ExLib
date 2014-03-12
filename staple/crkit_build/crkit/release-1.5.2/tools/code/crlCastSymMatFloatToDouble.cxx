
#include <itkImage.h>
#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkSymmetricSecondRankTensor.h>
#include <itkImageRegionIterator.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>

#include <iostream>

int main(int argc, char *argv[])
{

  static int const ImageDimension = 3;
  typedef   itk::SymmetricSecondRankTensor< float >    FloatTensorPixelType;
  typedef   itk::SymmetricSecondRankTensor< double >    DoubleTensorPixelType;

  typedef  itk::OrientedImage<FloatTensorPixelType,ImageDimension>  ImageTypeIn;
  typedef  itk::OrientedImage<DoubleTensorPixelType,ImageDimension>  ImageTypeOut;
  typedef   itk::ImageFileReader< ImageTypeIn >    ImageReaderType;
  typedef   itk::ImageFileWriter< ImageTypeOut >    ImageWriterType;
  typedef   itk::ImageRegionConstIterator< ImageTypeIn > InIteratorType;
  typedef   itk::ImageRegionIterator< ImageTypeOut > OutIteratorType;


  if (argc != 3) {
    std::cout << 
      "Usage: castsymmattodouble inimage outimage " << std::endl;
    std::cout << "This purpose of this program is to take a " 
     << "3D second order symmetric tensor with components stored as float "
     << "and modify it to have the components stored with a precision of "
     << "double." << std::endl;
    exit(1);
  }

  ImageReaderType::Pointer reader;
  ImageWriterType::Pointer writer;


  reader = ImageReaderType::New();
  reader->SetFileName(argv[1]);

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

  try {
    reader->Update();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  ImageTypeOut::Pointer outImage = ImageTypeOut::New();
  outImage->SetOrigin( reader->GetOutput()->GetOrigin());
  outImage->SetSpacing( reader->GetOutput()->GetSpacing() );
  outImage->SetDirection( reader->GetOutput()->GetDirection() );
  outImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  outImage->Allocate();

  writer = ImageWriterType::New();
  writer->SetFileName(argv[2]);

  /* Iterate over the input image for every pixel, and construct a new
   * pixel of the desired type by explicit casting.
   */
  InIteratorType it(reader->GetOutput(), 
                    reader->GetOutput()->GetLargestPossibleRegion());
  OutIteratorType ot(outImage, outImage->GetLargestPossibleRegion());

  for (it.GoToBegin(), ot.GoToBegin(); !it.IsAtEnd(); ++it, ++ot) {
    FloatTensorPixelType pv = it.Value();
    DoubleTensorPixelType pw;
    for (unsigned int i = 0; i < pv.GetNumberOfComponents(); i++) {
      pw[i] = static_cast<float>(pv[i]);
    }
    // Now set the output voxel to have the correct values.
    ot.Set(pw);
  }

  // Ensure the measurement frame is transferred to the output as well.
  itk::MetaDataDictionary &inDic = reader->GetOutput()->GetMetaDataDictionary();
  std::vector<std::vector<double> > msrFrame(ImageDimension);
  for (signed int j = 0; j < ImageDimension; j++) {
    msrFrame[j].resize(ImageDimension);
  }
  itk::ExposeMetaData(inDic, "NRRD_measurement frame", msrFrame);

  itk::MetaDataDictionary &outDic = outImage->GetMetaDataDictionary();
  itk::EncapsulateMetaData<std::vector<std::vector<double> > >(outDic,
     std::string("NRRD_measurement frame"), msrFrame);

  writer->SetInput( outImage );

  try {
    writer->Update();
    writer->UseCompressionOn();
    writer->Write();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  exit(0); // success
}

