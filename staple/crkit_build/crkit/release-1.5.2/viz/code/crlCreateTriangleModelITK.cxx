//
// Create a triangle model from a segmentation
//

#include <itkBinaryMask3DMeshSource.h>
#include <itkMesh.h>
#include <itkVTKPolyDataWriter.h>

#include "ImageTypeDefinitions.h"

int main(int argc, char *argv[])
{

  if (argc != 4) {
    std::cout << "Usage: " << argv[0] << " segmentedImageFile "
      << "threshold triangleFileName.vtk" << std::endl;
    return EXIT_FAILURE;
  }

  // IMPROVEMENTS: provide capability to pass decimation and smoothing
  // parameters.
  //   Provide an option to use vtkDecimateMarchingCubes for generating
  // multiple surfaces from segmented objects.
  double threshold = std::atof(argv[2]);

  // Read in the segmentation.
  ImageTypeDefinitions::ImageType::Pointer inputImage =
    ImageTypeDefinitions::ImageType::New();

  // Declare and initialize an image orienter filter
  itk::OrientImageFilter< ImageTypeDefinitions::ImageType,
                ImageTypeDefinitions::ImageType >::Pointer imageOrienter;
  imageOrienter = itk::OrientImageFilter< ImageTypeDefinitions::ImageType,                ImageTypeDefinitions::ImageType >::New();

  // Now load the image off disk and initialize the image data.
  ImageTypeDefinitions::ReaderType::Pointer reader =
                ImageTypeDefinitions::ReaderType::New();

  reader->SetFileName( argv[1] );
  try {
    reader->Update();
  } catch (itk::ExceptionObject & err) {
                std::cout << "ExceptionObject caught !" << std::endl;
                std::cout << err << std::endl;
                return false;
  }
  imageOrienter->UseImageDirectionOn();
  imageOrienter->SetDesiredCoordinateOrientationToAxial();
  imageOrienter->SetInput( reader->GetOutput() );
  inputImage = imageOrienter->GetOutput();
  try {
    inputImage->Update();
  } catch (itk::ExceptionObject & err) {
                std::cout << "ExceptionObject caught !" << std::endl;
                std::cout << err << std::endl;
                return false;
  }

  // Disconnect the DataObject inputImage from the pipeline.
  inputImage->DisconnectPipeline();


  itk::BinaryMask3DMeshSource<ImageTypeDefinitions::ImageType,itk::Mesh<float> >::Pointer marchingCubes;

  marchingCubes = itk::BinaryMask3DMeshSource<ImageTypeDefinitions::ImageType,itk::Mesh<float> >::New();
  marchingCubes->SetInput(inputImage);
  marchingCubes->SetObjectValue(threshold);
  marchingCubes->Update();

  itk::VTKPolyDataWriter<itk::Mesh<float> >::Pointer surfaceWriter;
  surfaceWriter = itk::VTKPolyDataWriter<itk::Mesh<float> >::New();
  surfaceWriter->SetInput(marchingCubes->GetOutput());
  surfaceWriter->SetFileName(argv[3]);
  surfaceWriter->Write();

  return EXIT_SUCCESS;
}
