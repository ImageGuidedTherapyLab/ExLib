
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSymmetricSecondRankTensor.h>
#include <itkDiffusionTensor3D.h>
#include <itkMatrix.h>
#include <itkUnaryFunctorImageFilter.h>
#include <itkMultiThreader.h>

#include "TensorScalarParameterFunctor.h"

#include "tclap/CmdLine.h"
#include "configuration.h"
#include <limits.h>

#ifdef WIN32
#define snprintf _snprintf
#define PATH_MAX 512
#endif

int main(int argc, char *argv[])
{

  // for FA
  bool computeFA = false;
  std::string faFileName;

  // for mean diffusivity
  bool computeMeanDiffusivity = false;
  std::string meanDiffusivityFileName;

  // for axial diffusivity
  bool computeAxialDiffusivity = false;
  std::string axialDiffusivityFileName;

  // for radial diffusivity
  bool computeRadialDiffusivity = false;
  std::string radialDiffusivityFileName;

  // for eigenvalue 1
  bool computeEigenvalue1 = false;
  std::string eigenValue1Name;

  // for eigenvalue 2
  bool computeEigenvalue2 = false;
  std::string eigenValue2Name;

  // for eigenvalue 3
  bool computeEigenvalue3 = false;
  std::string eigenValue3Name;

  // for components
  bool computeComponents = false;
  std::string componentsPrefixFileName;

  // for Cs
  bool computeShapeSpherical = false;
  std::string shapeSphericalName;

  // for Cp
  bool computeShapePlanar = false;
  std::string shapePlanarName;

  // for Cl
  bool computeShapeLinear = false;
  std::string shapeLinearName;

  // for Tensor Frobenius Norm = Tensor Scalar Product
  bool computeFrobeniusNorm = false;
  std::string frobeniusNormFileName;

  // the input image
  std::string inputImageFileName;


  // Some of the numerical code used here is not thread-safe
  itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);

  // Exceptions are thrown by TCLAP for bad arguments in the CLI.
  try {
    // Define the command line object, and insert a message
    // that describes the program.
    TCLAP::CmdLine cmd("Computational Radiology Laboratory", ' ', 
        CRKIT_VERSION_STRING );

    // Define a value argument for the input tensor image on the command line.
    // This can be specified without any flag. It is simply taken from the 
    // first unflagged position.
    TCLAP::UnlabeledValueArg<std::string> inputImageArg("input", 
            "Input Tensor Image", true, 
            "input file name", 
            "input tensor image file name",
            cmd); 

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> faArg("f", "fractionalAnisotropy", 
            "Construct Fractional Anistropy Image", false, 
            "/dev/null", "fa file name", 
            cmd); 

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> mdArg("m", "meanDiffusivity", 
            "Construct mean diffusivity Image", false, 
            "/dev/null", "md file name", 
            cmd); 

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> adArg("a", "axialDiffusivity", 
            "Construct axial diffusivity Image", false, 
            "/dev/null", "axial diffusivity file name", 
            cmd); 

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> rdArg("r", "radialDiffusivity", 
            "Construct radial diffusivity Image", false, 
            "/dev/null", "radial diffusivity file name", 
            cmd); 

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> eval1Arg("1", "eval1", 
            "Construct eigenvalue 1 (largest) image ", false, 
            "/dev/null", "eig val 1 file name", 
            cmd); 

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> eval2Arg("2", "eval2", 
            "Construct eigenvalue 2 (middle) image ", false, 
            "/dev/null", "eig val 2 file name", 
            cmd); 

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> eval3Arg("3", "eval3", 
            "Construct eigenvalue 3 (smallest) image ", false, 
            "/dev/null", "eig val 3 file name", 
            cmd); 

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> shapeSphericalArg("s", "shapeSpherical", 
            "Construct spherical tensor shape image", false, 
            "/dev/null", "Cs file name", 
            cmd); 

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> shapeLinearArg("l", "shapeLinear", 
            "Construct linear tensor shape image", false, 
            "/dev/null", "Cl file name", 
            cmd); 

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> shapePlanarArg("p", "shapePlanar", 
            "Construct planar tensor shape image", false, 
            "/dev/null", "Cp file name", 
            cmd); 

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> frobeniusNormArg("n", "frobeniusNorm", 
            "Construct Frobenius Norm image", false, 
            "/dev/null", "norm file name", 
            cmd); 

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> componentsArg("c", "components", 
            "Construct set of separate images of components ", false, 
            "/dev/null", "prefix of component file names", 
            cmd); 

    // Now parse the command line and assign variable values as required.
    cmd.parse( argc, argv );

    // Was an FA File requested ?
    if (faArg.isSet()) {
      computeFA = true;
      faFileName = faArg.getValue();
    }

    // Was a mean diffusivity file requested ?
    if (mdArg.isSet()) {
      computeMeanDiffusivity = true;
      meanDiffusivityFileName = mdArg.getValue();
    }

    // Was a axial diffusivity file requested ?
    if (adArg.isSet()) {
      computeAxialDiffusivity = true;
      axialDiffusivityFileName = adArg.getValue();
    }

    // Was a radial diffusivity file requested ?
    if (rdArg.isSet()) {
      computeRadialDiffusivity = true;
      radialDiffusivityFileName = rdArg.getValue();
    }

    // Was the planar shape file requested ?
    if (shapePlanarArg.isSet()) {
      computeShapePlanar = true;
      shapePlanarName = shapePlanarArg.getValue();
    }

    // Was the spherical shape file requested ?
    if (shapeSphericalArg.isSet()) {
      computeShapeSpherical = true;
      shapeSphericalName = shapeSphericalArg.getValue();
    }

    // Was the linear shape file requested ?
    if (shapeLinearArg.isSet()) {
      computeShapeLinear = true;
      shapeLinearName = shapeLinearArg.getValue();
    }

    // Was the eigenvalue 1 file requested ?
    if (eval1Arg.isSet()) {
      computeEigenvalue1 = true;
      eigenValue1Name = eval1Arg.getValue();
    }
    // Was the eigenvalue 2 file requested ?
    if (eval2Arg.isSet()) {
      computeEigenvalue2 = true;
      eigenValue2Name = eval2Arg.getValue();
    }
    // Was the eigenvalue 3 file requested ?
    if (eval3Arg.isSet()) {
      computeEigenvalue3 = true;
      eigenValue3Name = eval3Arg.getValue();
    }

    // Was the components files requested ?
    if (componentsArg.isSet()) {
      computeComponents = true;
      componentsPrefixFileName = componentsArg.getValue();
    }

    // Was the frobenius norm file requested ?
    if (frobeniusNormArg.isSet()) {
      computeFrobeniusNorm = true;
      frobeniusNormFileName = frobeniusNormArg.getValue();
    }

    if (inputImageArg.isSet()) {
      inputImageFileName = inputImageArg.getValue();
    }

  } catch (TCLAP::ArgException &e) 
  {
    std::cerr << "error: " << e.error() << " for arg " << 
                  e.argId() << std::endl;
    std::cerr << "This program computes parametric images from tensors." 
        << std::endl;
    exit(1);
  }

  static int const ImageDimension = 3;
  // Note the type used here must match the type checked for in the file.
  // typedef double   ComponentType;
  // As of 2007.06.22 NRRD format has some bugs in dealing with double
  // as a component type, so we will do some testing to see if float
  // will work.
  typedef float ComponentType;
  typedef   itk::DiffusionTensor3D<ComponentType>  PixelType;
  typedef   itk::Image< PixelType, ImageDimension >  InImageType;
  typedef   itk::ImageFileReader<InImageType> ImageReaderType;
  // Choose to store the scalar output in a float
  typedef   itk::Image< float, ImageDimension >  OutImageType;
  typedef   itk::ImageFileWriter<OutImageType> ImageWriterType;
  
  ImageReaderType::Pointer  reader = ImageReaderType::New();
  reader->SetFileName( inputImageFileName );

  try {
    reader->GenerateOutputInformation();
  } catch (itk::ExceptionObject &excp) 
  {
    std::cerr << "Failed to read from file " << inputImageFileName << std::endl;
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

  reader->Update();

  // For each type of output that is desired, put the tensor data through
  // a filter that creates it, and then write out the resulting scalar image.

  typedef itk::TensorFractionalAnisotropyImageFilter<InImageType, OutImageType> TensorFractionalAnisotropyImageFilterType;
  typedef itk::TensorMeanDiffusivityImageFilter<InImageType, OutImageType> TensorMeanDiffusivityImageFilterType;
  typedef itk::TensorRadialDiffusivityImageFilter<InImageType, OutImageType> TensorRadialDiffusivityImageFilterType;
  typedef itk::TensorAxialDiffusivityImageFilter<InImageType, OutImageType> TensorAxialDiffusivityImageFilterType;
  typedef itk::TensorLargestEigenValueImageFilter<InImageType, OutImageType> TensorLargestEigenValueImageFilterType;
  typedef itk::TensorMiddleEigenValueImageFilter<InImageType, OutImageType> TensorMiddleEigenValueImageFilterType;
  typedef itk::TensorSmallestEigenValueImageFilter<InImageType, OutImageType> TensorSmallestEigenValueImageFilterType;
  typedef itk::TensorFrobeniusNormImageFilter<InImageType, OutImageType> TensorFrobeniusNormImageFilterType;
  typedef itk::TensorComponentImageFilter<InImageType, OutImageType> TensorComponentImageFilterType;
  typedef itk::TensorPlanarShapeImageFilter<InImageType, OutImageType> \
    TensorPlanarShapeImageFilterType;
  typedef itk::TensorLinearShapeImageFilter<InImageType, OutImageType> \
    TensorLinearShapeImageFilterType;
  typedef itk::TensorSphericalShapeImageFilter<InImageType, OutImageType> \
    TensorSphericalShapeImageFilterType;

  if (computeFA) {
    TensorFractionalAnisotropyImageFilterType::Pointer tensorFilter = 
                   TensorFractionalAnisotropyImageFilterType::New();

    tensorFilter->SetInput( reader->GetOutput() );
    tensorFilter->Update();

    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetInput ( tensorFilter->GetOutput() );
    writer->UseCompressionOn( );
    writer->SetFileName( faFileName );
    try {
      writer->Write(); // 
    } catch ( itk::ExceptionObject & excp ) {
      std::cerr << "Failed to write out the image." << std::endl;
      std::cerr << "Exception error message is : " << excp << std::endl;
    }
  }

  if (computeMeanDiffusivity) {
    TensorMeanDiffusivityImageFilterType::Pointer tensorFilter = 
                   TensorMeanDiffusivityImageFilterType::New();

    tensorFilter->SetInput( reader->GetOutput() );
    tensorFilter->Update();

    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetInput ( tensorFilter->GetOutput() );
    writer->SetFileName( meanDiffusivityFileName );
    writer->UseCompressionOn( );
    try {
      writer->Write(); // 
    } catch ( itk::ExceptionObject & excp ) {
      std::cerr << "Failed to write out the image." << std::endl;
      std::cerr << "Exception error message is : " << excp << std::endl;
    }
  }

  if (computeAxialDiffusivity) {
    TensorAxialDiffusivityImageFilterType::Pointer tensorFilter = 
                   TensorAxialDiffusivityImageFilterType::New();

    tensorFilter->SetInput( reader->GetOutput() );
    tensorFilter->Update();

    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetInput ( tensorFilter->GetOutput() );
    writer->UseCompressionOn( );
    writer->SetFileName( axialDiffusivityFileName );
    try {
      writer->Write(); // 
    } catch ( itk::ExceptionObject & excp ) {
      std::cerr << "Failed to write out the image." << std::endl;
      std::cerr << "Exception error message is : " << excp << std::endl;
    }
  }

  if (computeRadialDiffusivity) {
    TensorRadialDiffusivityImageFilterType::Pointer tensorFilter = 
                   TensorRadialDiffusivityImageFilterType::New();

    tensorFilter->SetInput( reader->GetOutput() );
    tensorFilter->Update();

    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetInput ( tensorFilter->GetOutput() );
    writer->UseCompressionOn( );
    writer->SetFileName( radialDiffusivityFileName );
    try {
      writer->Write(); // 
    } catch ( itk::ExceptionObject & excp ) {
      std::cerr << "Failed to write out the image." << std::endl;
      std::cerr << "Exception error message is : " << excp << std::endl;
    }
  }

  if (computeShapePlanar) {
    TensorPlanarShapeImageFilterType::Pointer tensorFilter = 
                   TensorPlanarShapeImageFilterType::New();

    tensorFilter->SetInput( reader->GetOutput() );
    tensorFilter->Update();

    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetInput ( tensorFilter->GetOutput() );
    writer->UseCompressionOn( );
    writer->SetFileName( shapePlanarName );
    try {
      writer->Write(); // 
    } catch ( itk::ExceptionObject & excp ) {
      std::cerr << "Failed to write out the image." << std::endl;
      std::cerr << "Exception error message is : " << excp << std::endl;
    }
  }

  if (computeShapeLinear) {
    TensorLinearShapeImageFilterType::Pointer tensorFilter = 
                   TensorLinearShapeImageFilterType::New();

    tensorFilter->SetInput( reader->GetOutput() );
    tensorFilter->Update();

    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetInput ( tensorFilter->GetOutput() );
    writer->UseCompressionOn( );
    writer->SetFileName( shapeLinearName );
    try {
      writer->Write(); // 
    } catch ( itk::ExceptionObject & excp ) {
      std::cerr << "Failed to write out the image." << std::endl;
      std::cerr << "Exception error message is : " << excp << std::endl;
    }
  }

  if (computeShapeSpherical) {
    TensorSphericalShapeImageFilterType::Pointer tensorFilter = 
                   TensorSphericalShapeImageFilterType::New();

    tensorFilter->SetInput( reader->GetOutput() );
    tensorFilter->Update();

    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetInput ( tensorFilter->GetOutput() );
    writer->UseCompressionOn( );
    writer->SetFileName( shapeSphericalName );
    try {
      writer->Write(); // 
    } catch ( itk::ExceptionObject & excp ) {
      std::cerr << "Failed to write out the image." << std::endl;
      std::cerr << "Exception error message is : " << excp << std::endl;
    }
  }

  if (computeEigenvalue1) {
    TensorLargestEigenValueImageFilterType::Pointer tensorFilter = 
                   TensorLargestEigenValueImageFilterType::New();

    tensorFilter->SetInput( reader->GetOutput() );
    tensorFilter->Update();

    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetInput ( tensorFilter->GetOutput() );
    writer->UseCompressionOn( );
    writer->SetFileName( eigenValue1Name );
    try {
      writer->Write(); // 
    } catch ( itk::ExceptionObject & excp ) {
      std::cerr << "Failed to write out the image." << std::endl;
      std::cerr << "Exception error message is : " << excp << std::endl;
    }
  }

  if (computeEigenvalue2) {
    TensorMiddleEigenValueImageFilterType::Pointer tensorFilter = 
                   TensorMiddleEigenValueImageFilterType::New();

    tensorFilter->SetInput( reader->GetOutput() );
    tensorFilter->Update();

    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetInput ( tensorFilter->GetOutput() );
    writer->UseCompressionOn( );
    writer->SetFileName( eigenValue2Name );
    try {
      writer->Write(); // 
    } catch ( itk::ExceptionObject & excp ) {
      std::cerr << "Failed to write out the image." << std::endl;
      std::cerr << "Exception error message is : " << excp << std::endl;
    }
  }

  if (computeEigenvalue3) {
    TensorSmallestEigenValueImageFilterType::Pointer tensorFilter =
                   TensorSmallestEigenValueImageFilterType::New();

    tensorFilter->SetInput( reader->GetOutput() );
    tensorFilter->Update();

    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetInput ( tensorFilter->GetOutput() );
    writer->UseCompressionOn( );
    writer->SetFileName( eigenValue3Name );
    try {
      writer->Write(); // 
    } catch ( itk::ExceptionObject & excp ) {
      std::cerr << "Failed to write out the image." << std::endl;
      std::cerr << "Exception error message is : " << excp << std::endl;
    }
  }

  if (computeFrobeniusNorm) {
    TensorFrobeniusNormImageFilterType::Pointer tensorFilter =
                   TensorFrobeniusNormImageFilterType::New();

    tensorFilter->SetInput( reader->GetOutput() );
    tensorFilter->Update();

    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetInput ( tensorFilter->GetOutput() );
    writer->UseCompressionOn( );
    writer->SetFileName( frobeniusNormFileName );
    try {
      writer->Write(); // 
    } catch ( itk::ExceptionObject & excp ) {
      std::cerr << "Failed to write out the image." << std::endl;
      std::cerr << "Exception error message is : " << excp << std::endl;
    }
  }

  if (computeComponents) {
    char fileName[PATH_MAX];
    componentsPrefixFileName += "-c%03d.nrrd";
  
    // For each component, make a file, write out the component
    for (unsigned int i = 0; i < 6; i++) {
      memset(fileName, 0, PATH_MAX);
      TensorComponentImageFilterType::Pointer tensorFilter =
                   TensorComponentImageFilterType::New();

      tensorFilter->SetComponent( i );
      tensorFilter->SetInput( reader->GetOutput() );
      tensorFilter->Update();

      ImageWriterType::Pointer writer = ImageWriterType::New();
      writer->SetInput ( tensorFilter->GetOutput() );
      writer->UseCompressionOn( );
      // modify the file name
      snprintf(fileName, PATH_MAX, componentsPrefixFileName.c_str(), i);
std::cout << "File name is now " << fileName << std::endl;
      writer->SetFileName( fileName );
      try {
        writer->Write(); // 
      } catch ( itk::ExceptionObject & excp ) {
        std::cerr << "Failed to write out the image." << std::endl;
        std::cerr << "Exception error message is : " << excp << std::endl;
      }
    }
  }




  exit(0); // success
}

