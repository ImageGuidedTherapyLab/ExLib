//
// Create scalars for a surface model
//

#include <vtkPolyData.h>
#include <vtkProbeFilter.h>
#include <vtkImageData.h>
#include <vtkIdTypeArray.h>

#include <crlVtkMeshIO.h>
#include "ImageTypeDefinitions.h"

int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cout << "Usage: " << argv[0] << " surfaceModelFile imageFile outputModelFile" << std::endl;
    return EXIT_FAILURE;
  }

  // Read in the image
  ImageTypeDefinitions::ImageType::Pointer inputImage =
    ImageTypeDefinitions::ImageType::New();

  // Declare and initialize an image orienter filter
  itk::OrientImageFilter< ImageTypeDefinitions::ImageType,
    ImageTypeDefinitions::ImageType >::Pointer imageOrienter;
  imageOrienter = itk::OrientImageFilter< ImageTypeDefinitions::ImageType, ImageTypeDefinitions::ImageType >::New();

  // Now load the image off disk and initialize the image data.
  ImageTypeDefinitions::ReaderType::Pointer reader =
    ImageTypeDefinitions::ReaderType::New();

  reader->SetFileName( argv[2] );
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

  ImageTypeDefinitions::ITK2VTKConnectorFilterType::Pointer
    ITK2VTKConnector =
    ImageTypeDefinitions::ITK2VTKConnectorFilterType::New();
  ITK2VTKConnector->SetInput( inputImage );

  // Read the surface
  vtkPolyData *surface = crlVtkMeshIO::ReadMesh(std::string(argv[1]));

  // Set up a pipeline to create scalars
  vtkProbeFilter *probe = vtkProbeFilter::New();
  probe->SetInputConnection(surface->GetProducerPort());
  probe->SetSource(ITK2VTKConnector->GetOutput());
  probe->Update();
  std::cout << probe->GetValidPoints()->GetNumberOfTuples() << std::endl;

  crlVtkMeshIO::WriteMesh(dynamic_cast<vtkPolyData*>(probe->GetOutput()),argv[3]);

  surface->Delete();
  probe->Delete();

  return EXIT_SUCCESS;

}
