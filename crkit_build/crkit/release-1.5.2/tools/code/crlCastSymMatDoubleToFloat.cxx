
#include <itkImage.h>
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

  typedef   itk::Image<DoubleTensorPixelType,ImageDimension>  ImageTypeIn;
  typedef   itk::Image<FloatTensorPixelType,ImageDimension>  ImageTypeOut;
  typedef   itk::ImageFileReader< ImageTypeIn >    ImageReaderType;
  typedef   itk::ImageFileWriter< ImageTypeOut >    ImageWriterType;
  typedef   itk::ImageRegionConstIterator< ImageTypeIn > InIteratorType;
  typedef   itk::ImageRegionIterator< ImageTypeOut > OutIteratorType;


  if (argc != 3) {
    std::cout << 
      "Usage: castsymmattofloat inimage outimage " << std::endl;
    std::cout << "This purpose of this program is to take a " 
     << "3D second order symmetric tensor with components stored as double "
     << "and modify it to have the components stored with a precision of "
     << "float." << std::endl;
    exit(1);
  }

  ImageReaderType::Pointer r;
  ImageWriterType::Pointer w;

  r = ImageReaderType::New();
  r->SetFileName(argv[1]);

  try {
    r->GenerateOutputInformation();
  } catch (itk::ExceptionObject &excp)
  {
    std::cerr << "Failed to read from file " << argv[1] << std::endl;
    std::cerr << "Exception error message is : " << excp << std::endl;
    return EXIT_FAILURE;
  }

  itk::ImageIOBase *imageIO = r->GetImageIO();
  if (imageIO->GetComponentType() != itk::ImageIOBase::DOUBLE) {
    std::cerr << "File component type does not match expected type DOUBLE."
              << std::endl;
    std::cerr << "Use crlCastSymMatToFloat to change the type of the components."
              << std::endl;
    return EXIT_FAILURE;
  }

  try {
    r->Update();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  ImageTypeOut::Pointer outImage = ImageTypeOut::New();
  outImage->SetSpacing( r->GetOutput()->GetSpacing() );
  outImage->SetOrigin( r->GetOutput()->GetOrigin());
  outImage->SetRegions( r->GetOutput()->GetLargestPossibleRegion().GetSize() );
  outImage->Allocate();
  outImage->SetDirection( r->GetOutput()->GetDirection());

  w = ImageWriterType::New();
  w->SetFileName(argv[2]);

  /* Iterate over the input image for every pixel, and construct a new
   * pixel of the desired type by explicit casting.
   */
  InIteratorType it(r->GetOutput(), r->GetOutput()->GetLargestPossibleRegion());
  OutIteratorType ot(outImage, outImage->GetLargestPossibleRegion());

  for (it.GoToBegin(), ot.GoToBegin(); !it.IsAtEnd(); ++it, ++ot) {
    DoubleTensorPixelType pv = it.Value();
    FloatTensorPixelType pw;
    for (unsigned int i = 0; i < pv.GetNumberOfComponents(); i++) {
      pw[i] = static_cast<float>(pv[i]);
    }
    // Now set the output voxel to have the correct values.
    ot.Set(pw);
  }

  // Ensure the measurement frame is transferred to the output as well.
  itk::MetaDataDictionary &inDic = r->GetOutput()->GetMetaDataDictionary();
  std::vector<std::vector<double> > msrFrame(ImageDimension);
  for (signed int j = 0; j < ImageDimension; j++) {
    msrFrame[j].resize(ImageDimension);
  }
  itk::ExposeMetaData(inDic, "NRRD_measurement frame", msrFrame);

  itk::MetaDataDictionary &outDic = outImage->GetMetaDataDictionary();
  itk::EncapsulateMetaData<std::vector<std::vector<double> > >(outDic,
     std::string("NRRD_measurement frame"), msrFrame);

  w->SetInput( outImage );

  try {
    w->Update();
    w->UseCompressionOn();
    w->Write();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  exit(0); // success
}

