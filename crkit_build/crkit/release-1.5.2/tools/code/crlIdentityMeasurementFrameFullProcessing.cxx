#include <tclap/CmdLine.h>

#include <itkAffineTransform.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

int main( int argc, char * argv[] )
{
  std::string inputImageFile;
  std::string outputImageFile;

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',"1.0");
      TCLAP::UnlabeledValueArg<std::string> inputImageFileArg("inputImageFile","has to be a 4D DWI image.",true,"","input DWI image file",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputImageFileArg("outputImageFile","result is written here.",true,"","output DWI image file",cmd);

      cmd.parse(argc,argv);

      if (inputImageFileArg.isSet()) inputImageFile = inputImageFileArg.getValue();
      if (outputImageFileArg.isSet()) outputImageFile = outputImageFileArg.getValue();
    }
  catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(1);
    }

  typedef itk::AffineTransform< double, 3 >  AffineTransformType;

  typedef itk::VectorImage<signed short,  3 >   DWIImageType;

  typedef itk::ImageFileReader< DWIImageType  >  DWIReaderType;
  typedef itk::ImageFileWriter< DWIImageType  >  DWIWriterType;

  DWIReaderType::Pointer reader = DWIReaderType::New();
  DWIWriterType::Pointer writer = DWIWriterType::New();
  reader->SetFileName(inputImageFile.c_str());
  writer->SetFileName(outputImageFile.c_str());

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

  itk::Matrix<double,3,3> measurementFrame;
  measurementFrame.SetIdentity();

  AffineTransformType::Pointer mfTransform = AffineTransformType::New();

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
      std::cout << "No measurement frame found in NRRD file." << std::endl;
      return 1;
    }

  mfTransform->SetMatrix(measurementFrame);
  //mfTransform->Print(std::cout);

  AffineTransformType::Pointer inverse = AffineTransformType::New();
  mfTransform->GetInverse(inverse);

  int k = -1;
  std::stringstream ss;

  do
    {
      k++;
      ss.str("");
      ss << "DWMRI_gradient_" << std::right << std::setw(4) << std::setfill('0') << k;
      if (mdd.HasKey(ss.str()))
	{
	  std::string strDirection;
	  itk::ExposeMetaData<std::string>(mdd,ss.str(),strDirection);

	  itk::Vector<double,3> direction;
	  sscanf(strDirection.c_str(),"%lf %lf %lf",&direction[0],&direction[1],&direction[2]);
	  //std::cout << direction[0] << " " << direction[1] << " " << direction[2] << std::endl;

	  itk::Vector<double,3> vec;
	  vec = inverse->TransformVector(direction);

	  std::stringstream output;
	  output << vec[0] << " " << vec[1] << " " << vec[2];
	  itk::EncapsulateMetaData<std::string>(mdd,ss.str(),output.str());
	  //std::cout << output.str() << std::endl << std::endl;
	}
    }
  while (mdd.HasKey(ss.str()));

  std::vector< std::vector<double> > idFrame;
  idFrame.resize(   3);
  idFrame[0].resize(3);
  idFrame[1].resize(3);
  idFrame[3].resize(3);
  idFrame[0][0] = 1;
  idFrame[0][1] = 0;
  idFrame[0][2] = 0;
  idFrame[1][0] = 0;
  idFrame[1][1] = 1;
  idFrame[1][2] = 0;
  idFrame[2][0] = 0;
  idFrame[2][1] = 0;
  idFrame[2][2] = 1;
  itk::EncapsulateMetaData<std::vector<std::vector<double> > >(mdd,"NRRD_measurement frame",idFrame);

  writer->UseCompressionOn();
  writer->SetInput(reader->GetOutput());
  reader->GetOutput()->SetMetaDataDictionary(mdd);

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
