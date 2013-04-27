
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

#include "crlVector3ToRGB.h"

// Compute an RGB image representing a color encoded 3 vector image.

int main(int argc, const char *argv[])
{

  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " vector3Image rgbImage" << std::endl;
    exit(1);
  }

  std::string inName(argv[1]);
  std::string outName(argv[2]);

  const unsigned int ImageDimension = 3;
  typedef   itk::Vector<double, ImageDimension>  VectorPixelType;
  typedef   itk::Image< VectorPixelType, ImageDimension >  VectorImageType;
  typedef   itk::ImageFileReader<VectorImageType> ImageReaderType;

  ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName( inName.c_str() );

  try {
    reader->GenerateOutputInformation();
  } catch (itk::ExceptionObject &excp)
  {
    std::cerr << "Failed to read from file " << inName << std::endl;    
    std::cerr << "Exception error message is : " << excp << std::endl;
    return EXIT_FAILURE;
  }

  crl::Vector3ToRGBBase::Pointer tsm;

  itk::ImageIOBase *imageIO = reader->GetImageIO();
  if (imageIO->GetComponentType() == itk::ImageIOBase::DOUBLE) {
   // Component type of data on file is DOUBLE, so use filter templated
   // over double for this operation.
   tsm = crl::Vector3ToRGB<ImageDimension, double>::New();
  } else {
   // Component type of data on file is FLOAT, so use filter templated
   // over float for this operation.
   tsm = crl::Vector3ToRGB<ImageDimension, float>::New();
  }

  tsm->SetInputFileName( inName.c_str() );
  tsm->SetOutputFileName( outName.c_str() );

  int ret = tsm->Execute();

  if (ret != 0) {
    std::cerr << argv[0] << " failed with return error code " << ret <<
    std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

};

#include "crlVector3ToRGB.txx"

