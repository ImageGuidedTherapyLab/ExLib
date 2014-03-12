
#include "crlTensorShapeMeasures.h"
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

// Estimate Cs, Cp and Cl for each tensor.
//
// By estimating all three together, an adaptive noise parameter can be
// computed to ensure that they meet certain constraints, such as 
// sum to 1.

int main(int argc, const char *argv[])
{

  if (argc != 5) {
    std::cerr << "Usage: " << argv[0] << " intensor csimage climage cpimage" << std::endl;
    std::cerr << "This program computes Cs, Cp and Cl." << std::endl;
    exit(1);
  }

  std::string inName(argv[1]);
  std::string outCsName(argv[2]);
  std::string outCpName(argv[3]);
  std::string outClName(argv[4]);

  typedef double ComponentTypeDouble;
  typedef float  ComponentTypeFloat;
  const unsigned int ImageDimension = 3;
  typedef   itk::SymmetricSecondRankTensor<ComponentTypeDouble, ImageDimension>  PixelTypeDouble;
  typedef   itk::Image< PixelTypeDouble, ImageDimension >  ImageTypeDouble;
  typedef   itk::ImageFileReader<ImageTypeDouble> ImageReaderTypeDouble;

  ImageReaderTypeDouble::Pointer reader = ImageReaderTypeDouble::New();
  reader->SetFileName( inName.c_str() );

  try {
    reader->GenerateOutputInformation();
  } catch (itk::ExceptionObject &excp)
  {
    std::cerr << "Failed to read from file " << inName << std::endl;    
    std::cerr << "Exception error message is : " << excp << std::endl;
    return EXIT_FAILURE;
  }

  crl::TensorShapeMeasuresBase::Pointer tsm;

  itk::ImageIOBase *imageIO = reader->GetImageIO();
  if (imageIO->GetComponentType() == itk::ImageIOBase::DOUBLE) {
   // Component type of data on file is DOUBLE, so use filter templated
   // over double for this operation.
   tsm = crl::TensorShapeMeasures<ImageDimension, double>::New();
  } else if (imageIO->GetComponentType() == itk::ImageIOBase::FLOAT) {
   // Component type of data on file is FLOAT, so use filter templated
   // over float for this operation.
   tsm = crl::TensorShapeMeasures<ImageDimension, float>::New();
  } else {
    std::cerr << "File component type does not match expected type "
       << " of float or double." << std::endl;
    std::cerr << "Use crlCastSymMat* to change the type of the components."
              << std::endl;
    return EXIT_FAILURE;
  }

  tsm->SetInputFileName( inName.c_str() );
  tsm->SetOutputCsFileName( outCsName.c_str() );
  tsm->SetOutputCpFileName( outCpName.c_str() );
  tsm->SetOutputClFileName( outClName.c_str() );

  int ret = tsm->Execute();

  if (ret != 0) {
    std::cerr << argv[0] << " failed with return error code " << ret <<
    std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

};

#include "crlTensorShapeMeasures.txx"

