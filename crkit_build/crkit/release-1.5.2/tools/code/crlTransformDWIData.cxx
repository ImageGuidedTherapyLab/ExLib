#include <tclap/CmdLine.h>

#include <itkTransformFileReader.h>
#include <itkAffineTransform.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkVariableLengthVector.h>
#include <itkImageToVectorImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkWindowedSincInterpolateImageFunction.h>


int main( int argc, char * argv[] )
{
  std::string inputImageFile;
  std::string inputTransformFile;
  std::string outputGeometryExampleFile;
  std::string interpolationStyle;
  std::string outputImageFile;

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',"1.0");
      TCLAP::UnlabeledValueArg<std::string> inputImageFileArg("inputImageFile","may be in any ITK format.",true,"","input DWI image file",cmd);
      TCLAP::UnlabeledValueArg<std::string> inputTransformFileArg("inputTransformFile","an ITK transform file.",true,"","input transform file",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputGeometryExampleFileArg("outputGeometryExampleFile","an image file with the desired geometry.",true,"","output geometry example file",cmd);
      TCLAP::UnlabeledValueArg<std::string> interpolationStyleArg("interpolationStyle","one of linear, nearest, bspline, sinc.",true,"","interpolation style",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputImageFileArg("outputImageFile","result is written here.",true,"","output DWI image file",cmd);

      cmd.parse(argc,argv);

      if (inputImageFileArg.isSet()) inputImageFile = inputImageFileArg.getValue();
      if (inputTransformFileArg.isSet()) inputTransformFile = inputTransformFileArg.getValue();
      if (outputGeometryExampleFileArg.isSet()) outputGeometryExampleFile = outputGeometryExampleFileArg.getValue();
      if (interpolationStyleArg.isSet()) interpolationStyle = interpolationStyleArg.getValue();
      if (outputImageFileArg.isSet()) outputImageFile = outputImageFileArg.getValue();
    }
  catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(1);
    }

  itk::TransformFileReader::Pointer trsfreader;
  trsfreader = itk::TransformFileReader::New();
  trsfreader->SetFileName( inputTransformFile.c_str() );

  try
    {
      trsfreader->Update();
    }
  catch ( itk::ExceptionObject & excp )
    {
      // Display error from reading the transform file.
      std::cerr << "Error while reading the transform file " <<
	(inputTransformFile) << std::endl;
      std::cerr << excp << std::endl;
      std::cerr << "[FAILED]" << std::endl;
      return 1;
    }

  typedef itk::AffineTransform< double, 3 >  AffineTransformType;
  typedef itk::TransformFileReader::TransformListType * TransformListType;
  TransformListType transforms = trsfreader->GetTransformList();
  std::cout << "Number of transforms = " << transforms->size() << std::endl;
  itk::TransformFileReader::TransformListType::const_iterator it =
    transforms->begin();
  if (transforms->size() <= 0 || transforms->size() > 1) {
    std::cerr << "Read " << transforms->size() << " transforms but want 1." << std::endl;
    return 1;
  }
  AffineTransformType::Pointer affine_read;
  if (!strcmp((*it)->GetNameOfClass(), "AffineTransform"))
    {
      affine_read = static_cast<AffineTransformType*>((*it).GetPointer());
      affine_read->Print(std::cout);
    } else {
      std::cerr << "Don't know how to resample using a " <<
	(*it)->GetNameOfClass() << " transform." << std::endl;
      return 1;
    }

  typedef itk::VectorImage<signed short,  3 >   DWIImageType;
  typedef itk::OrientedImage<signed short, 3> ImageType;

  typedef itk::ImageFileReader< DWIImageType  >  DWIReaderType;
  typedef itk::ImageFileWriter< DWIImageType  >  DWIWriterType;
  typedef itk::ImageFileReader< ImageType > ImageReaderType;

  DWIReaderType::Pointer reader = DWIReaderType::New();
  DWIWriterType::Pointer writer = DWIWriterType::New();
  ImageReaderType::Pointer geometry = ImageReaderType::New();
  reader->SetFileName(inputImageFile.c_str());
  writer->SetFileName(outputImageFile.c_str());
  geometry->SetFileName(outputGeometryExampleFile.c_str());

  try
    {
      reader->Update();
    }
  catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Error while reading the input DWI file " << inputImageFile << std::endl;
      std::cerr << excp << std::endl;
      std::cerr << "[FAILED]" << std::endl;
      return 1;
    }

  try
    {
      geometry->Update();
    }
  catch (itk::ExceptionObject &excp)
    {
      std::cerr << "Error while reading the geometry example file " << outputGeometryExampleFile << std::endl;
      std::cerr << excp << std::endl;
      std::cerr << "[FAILED]" << std::endl;
      return 1;
    }

  itk::Matrix<double,3,3> measurementFrame;
  measurementFrame.SetIdentity();

  itk::MetaDataDictionary &mdd = reader->GetOutput()->GetMetaDataDictionary();
  if (mdd.HasKey("NRRD_measurement frame"))
    {
      std::cout << "Measurement frame found in NRRD file: (";
      std::vector<std::vector<double> > frame;
      itk::ExposeMetaData<std::vector<std::vector<double> > >(mdd,"NRRD_measurement frame",frame);
      std::cout << frame[0][0] << " " << frame[0][1] << " " << frame[0][2] << ") (";
      std::cout << frame[1][0] << " " << frame[1][1] << " " << frame[1][2] << ") (";
      std::cout << frame[2][0] << " " << frame[2][1] << " " << frame[2][2] << ")" << std::endl;
      measurementFrame[0][0] = frame[0][0];
      measurementFrame[0][1] = frame[0][1];
      measurementFrame[0][2] = frame[0][2];
      measurementFrame[1][0] = frame[1][0];
      measurementFrame[1][1] = frame[1][1];
      measurementFrame[1][2] = frame[1][2];
      measurementFrame[2][0] = frame[2][0];
      measurementFrame[2][1] = frame[2][1];
      measurementFrame[2][2] = frame[2][2];
    }
  else
    {
      std::cout << "No measurement frame found in NRRD file, assuming identity." << std::endl;
    }

  std::vector<std::vector<double> > frame;
  frame.resize(3);
  frame[0].resize(3);
  frame[1].resize(3);
  frame[2].resize(3);

  AffineTransformType::Pointer mmft = AffineTransformType::New();
  affine_read->GetInverse(mmft);

  itk::Vector<double,3> vec;
  vec = mmft->TransformVector(measurementFrame[0]);
  frame[0][0] = vec[0];
  frame[0][1] = vec[1];
  frame[0][2] = vec[2];

  vec = mmft->TransformVector(measurementFrame[1]);
  frame[1][0] = vec[0];
  frame[1][1] = vec[1];
  frame[1][2] = vec[2];

  vec = mmft->TransformVector(measurementFrame[2]);
  frame[2][0] = vec[0];
  frame[2][1] = vec[1];
  frame[2][2] = vec[2];

  writer->UseCompressionOn();

  // since there seems to be no working ITK resampler for variable length vector images, we have to resample each component separately
  itk::ImageToVectorImageFilter<ImageType>::Pointer composer = itk::ImageToVectorImageFilter<ImageType>::New();
  for (unsigned int c = 0; c < reader->GetOutput()->GetNumberOfComponentsPerPixel(); c++)
    {
      ImageType::Pointer currentComponent = ImageType::New();
      currentComponent->SetRegions(reader->GetOutput()->GetLargestPossibleRegion().GetSize());
      currentComponent->SetDirection(reader->GetOutput()->GetDirection());
      currentComponent->SetOrigin(reader->GetOutput()->GetOrigin());
      currentComponent->SetSpacing(reader->GetOutput()->GetSpacing());
      currentComponent->Allocate();

      itk::ImageRegionIterator<DWIImageType> imgIt(reader->GetOutput(),reader->GetOutput()->GetLargestPossibleRegion());
      itk::ImageRegionIterator<ImageType> compIt(currentComponent,currentComponent->GetLargestPossibleRegion());

      for (imgIt.GoToBegin(), compIt.GoToBegin(); !compIt.IsAtEnd(); ++imgIt, ++compIt)
	{
	  compIt.Set(imgIt.Get()[c]);
	}
      itk::ResampleImageFilter<ImageType,ImageType>::Pointer resampler = itk::ResampleImageFilter<ImageType,ImageType>::New();
      resampler->SetInput(currentComponent);
      if (interpolationStyle == "linear")
	resampler->SetInterpolator(itk::LinearInterpolateImageFunction<ImageType,double>::New());
      else if (interpolationStyle == "nearest")
	resampler->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<ImageType,double>::New());
      else if (interpolationStyle == "bspline")
	resampler->SetInterpolator(itk::BSplineInterpolateImageFunction<ImageType,double>::New());
      else if (interpolationStyle == "sinc")
	resampler->SetInterpolator(itk::WindowedSincInterpolateImageFunction<ImageType,5,itk::Function::HammingWindowFunction<5>,itk::ConstantBoundaryCondition<ImageType>,double>::New());
      else
	{
	  std::cerr << "Don't know how to interpolate using " << interpolationStyle  << std::endl;
	  std::cerr << "Options are: linear nearest bspline sinc" << std::endl;
	  return 1;
	}
      resampler->SetOutputDirection(geometry->GetOutput()->GetDirection());
      resampler->SetOutputParametersFromImage(geometry->GetOutput());
      resampler->SetDefaultPixelValue(0);
      resampler->SetTransform(affine_read);
      resampler->Update();
      composer->SetNthInput(c,resampler->GetOutput());
      composer->Update();
    }

  writer->SetInput(composer->GetOutput());
  writer->UseCompressionOn();

  itk::EncapsulateMetaData<std::vector<std::vector<double> > >(mdd,"NRRD_measurement frame",frame);
  composer->GetOutput()->SetMetaDataDictionary(mdd);

  /*mdd = writer->GetInput()->GetMetaDataDictionary();
  if (mdd.HasKey("NRRD_kinds[0]"))
    {
      for (unsigned int k = 0; k < 4; k++)
	{
	  std::string kind;
	  std::stringstream ss;
	  ss << "NRRD_kinds[" << k << "]";
	  itk::ExposeMetaData<std::string>(mdd,ss.str(),kind);
	  std::cout << kind << std::endl;
	}
     }*/

  try
    {
      writer->Write();
    }
  catch (itk::ExceptionObject &excp)
    {
      std::cerr << "Error while writing the output DWI file " << outputImageFile << std::endl;
      std::cerr << excp << std::endl;
      std::cerr << "[FAILED]" << std::endl;
      return 1;
    }
}
