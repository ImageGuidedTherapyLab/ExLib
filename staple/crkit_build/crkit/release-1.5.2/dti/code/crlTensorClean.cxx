
#include "crlTensorClean.h"
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "configuration.h"
#include "tclap/CmdLine.h"
#include <limits.h>

#ifdef WIN32
#define snprintf _snprintf
#define PATH_MAX 512
#endif

#include <math.h>
#include <float.h>

#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_det.h>

// Count the number of good tensors in a file.
// Count the number of bad tensors in a file.

int main(int argc, char *argv[])
{

  // Define the command line object, and insert a message
  // that describes the program.
  TCLAP::CmdLine cmd("Computational Radiology Laboratory", ' ',
        CRKIT_VERSION_STRING );

  // Define a value argument and add it to the command line.
  TCLAP::ValueArg<std::string> inputArg("i", "inputFile",
            "Input Tensor File Name", true,
            "/dev/null", "input tensor file name",
            cmd);

  TCLAP::ValueArg<std::string> outputArg("o", "outputFile",
            "Output Tensor File Name", true,
            "/dev/null", "output tensor file name",
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

  typedef double ComponentTypeDouble;
  typedef float  ComponentTypeFloat;
  const unsigned int ImageDimension = 3;
  typedef   itk::SymmetricSecondRankTensor<ComponentTypeDouble, ImageDimension>  PixelTypeDouble;
  typedef   itk::Image< PixelTypeDouble, ImageDimension >  ImageTypeDouble;
  typedef   itk::ImageFileReader<ImageTypeDouble> ImageReaderTypeDouble;

  ImageReaderTypeDouble::Pointer reader = ImageReaderTypeDouble::New();
  reader->SetFileName( inputName.c_str() );

  try {
    reader->GenerateOutputInformation();
  } catch (itk::ExceptionObject &excp)
  {
    std::cerr << "Failed to read from file " << inputName << std::endl;    
    std::cerr << "Exception error message is : " << excp << std::endl;
    return EXIT_FAILURE;
  }

  TensorCheckBase::Pointer tcb;

  itk::ImageIOBase *imageIO = reader->GetImageIO();
  if (imageIO->GetComponentType() == itk::ImageIOBase::DOUBLE) {
   // Component type of data on file is DOUBLE, so use filter templated
   // over double for this operation.
   tcb = TensorCheck<ImageDimension, double>::New();
   tcb->SetInputFileName( inputName.c_str() );
   tcb->SetOutputFileName( outputName.c_str() );
   tcb->SetUseCompression( compressOutput );
  } else if (imageIO->GetComponentType() == itk::ImageIOBase::FLOAT) {
   // Component type of data on file is FLOAT, so use filter templated
   // over float for this operation.
   tcb = TensorCheck<ImageDimension, float>::New();
   tcb->SetInputFileName( inputName.c_str() );
   tcb->SetOutputFileName( outputName.c_str() );
   tcb->SetUseCompression( compressOutput );
  } else {
    std::cerr << "File component type does not match expected type "
       << " of float or double." << std::endl;
    std::cerr << "Use castsymmat to change the type of the components."
              << std::endl;
    return EXIT_FAILURE;
  }

  int ret = tcb->Execute();

  if (ret != 0) {
    std::cerr << "TensorCheck failed with return error code " << ret <<
    std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

};

#include "crlTensorClean.txx"

