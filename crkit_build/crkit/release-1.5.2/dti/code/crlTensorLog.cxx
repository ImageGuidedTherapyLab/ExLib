
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSymmetricSecondRankTensor.h>
#include <itkDiffusionTensor3D.h>
#include <itkMatrix.h>
#include <itkUnaryFunctorImageFilter.h>
#include <itkMultiThreader.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <typeinfo>

#include "TensorLogFunctor.h"

#include "tclap/CmdLine.h"
#include "configuration.h"
#include <limits.h>

#ifdef WIN32
#define snprintf _snprintf
#define PATH_MAX 512
#endif

int main(int argc, char *argv[])
{

  std::string inputTensorFileName;
  std::string outputLogTensorFileName;

  // Define the command line object, and insert a message
  // that describes the program.
  TCLAP::CmdLine cmd("Computational Radiology Laboratory", ' ',
        CRKIT_VERSION_STRING );

  TCLAP::ValueArg<std::string> outputArg("o", "outputFile",
            "Output Tensor File Name", true,
            "/dev/null", "output tensor file name",
            cmd);

  TCLAP::ValueArg<std::string> inputArg("i", "inputFile",
            "Input Tensor File Name", true,
            "/dev/null", "input tensor file name",
            cmd);

  TCLAP::SwitchArg compressOutputArg("z", "compressOutput",
       "compress data in output files", cmd, false);

  try {
    cmd.parse(argc, argv);
  } catch (TCLAP::ArgException &e) 
  {
    std::cerr << "Error: " << e.error() << "for argument " <<
      e.argId() << std::endl;
    exit (1);
  }

  const std::string& inputName = inputArg.getValue();
  const std::string& outputName = outputArg.getValue();
  bool compressOutput = compressOutputArg.getValue();

  static int const ImageDimension = 3;
  // typedef double ComponentType;
  typedef float ComponentType; 
             // Automatic component type conversion 
            // fails for NRRD data.
  typedef   itk::DiffusionTensor3D<ComponentType>  PixelType;
  typedef   itk::Image< PixelType, ImageDimension >  ImageType;
  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  typedef   itk::ImageFileWriter<ImageType> ImageWriterType;


// Some of the numerical code used here is not thread-safe
  itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);
  
  ImageReaderType::Pointer  reader = ImageReaderType::New();
  reader->SetFileName( inputName );

  try {
    reader->GenerateOutputInformation();
  } catch (itk::ExceptionObject &excp) 
  {
    std::cerr << "Failed to read from file " << argv[1] << std::endl;
    std::cerr << "Exception error message is : " << excp << std::endl;
    return EXIT_FAILURE;
  }

  itk::ImageIOBase *imageIO = reader->GetImageIO();
  if (imageIO->GetComponentTypeInfo() != typeid(ComponentType)) {
  // if (imageIO->GetComponentType() != itk::ImageIOBase::FLOAT) {
    std::cerr << "File component type does not match expected type FLOAT."
              << std::endl;
    std::cerr << "Use crlCastSymMatToFloat to change the type of the components." 
              << std::endl;
    return EXIT_FAILURE;
  }

  reader->Update();

  typedef itk::TensorLogImageFilter<ImageType, ImageType> TensorLogImageFilterType;
  TensorLogImageFilterType::Pointer tensorLogFilter = 
                                            TensorLogImageFilterType::New();
  tensorLogFilter->SetInput( reader->GetOutput() );
  tensorLogFilter->Update();

  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput ( tensorLogFilter->GetOutput() );
  writer->SetFileName( outputName );

  // Transfer the measurement frame from the input data to the output data 
  // Ensure the measurement frame is transferred to the output as well.
  itk::MetaDataDictionary &inDic = reader->GetOutput()->GetMetaDataDictionary();  std::vector<std::vector<double> > msrFrame(ImageDimension);
  for (signed int j = 0; j < ImageDimension; j++) {
    msrFrame[j].resize(ImageDimension);
  }
  itk::ExposeMetaData(inDic, "NRRD_measurement frame", msrFrame);

  itk::MetaDataDictionary &outDic = 
             tensorLogFilter->GetOutput()->GetMetaDataDictionary();
  itk::EncapsulateMetaData<std::vector<std::vector<double> > >(outDic,
     std::string("NRRD_measurement frame"), msrFrame);

  try {
    if (compressOutput) {
      writer->UseCompressionOn();
    }
    writer->Write();
  } catch ( itk::ExceptionObject & excp ) {
    std::cerr << "Failed to write out the image." << std::endl;
    std::cerr << "Exception error message is : " << excp << std::endl;
  }

  exit(0); // success
}

