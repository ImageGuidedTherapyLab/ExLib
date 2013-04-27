//
// Create a triangle model from a segmentation
//

#include <vtkImageData.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkMatrix4x4.h>

#include "ImageTypeDefinitions.h"
#include <crlVtkMeshIO.h>

int main(int argc, char *argv[])
{

  if (argc != 4) {
    std::cout << "Usage: " << argv[0] << " originalModelFile imageFile rotatedModelFile" << std::endl;
    return EXIT_FAILURE;
  }

  // Read in the model.
  vtkPolyData *surface = crlVtkMeshIO::ReadMesh(argv[1]);

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
  ImageTypeDefinitions::ImageType::DirectionType direction = reader->GetOutput()->GetDirection();

  double translation[3];
  vtkMatrix4x4 *rotation = vtkMatrix4x4::New();

  ImageTypeDefinitions::ImageType::IndexType idx;
  ImageTypeDefinitions::ImageType::PointType point;
  idx[0]=0; idx[1]=0; idx[2]=0;
  reader->GetOutput()->TransformIndexToPhysicalPoint(idx,point);
  translation[0] = point[0];
  translation[1] = point[1];
  translation[2] = point[2];
  rotation->SetElement(0,0,direction[0][0]);
  rotation->SetElement(0,1,direction[0][1]);
  rotation->SetElement(0,2,direction[0][2]);
  rotation->SetElement(1,0,direction[1][0]);
  rotation->SetElement(1,1,direction[1][1]);
  rotation->SetElement(1,2,direction[1][2]);
  rotation->SetElement(2,0,direction[2][0]);
  rotation->SetElement(2,1,direction[2][1]);
  rotation->SetElement(2,2,direction[2][2]);

  vtkTransform *transform = vtkTransform::New();
  vtkTransformPolyDataFilter *tpdf = vtkTransformPolyDataFilter::New();
  tpdf->SetInputConnection(surface->GetProducerPort());
  transform->PostMultiply();
  transform->Translate(-translation[0],-translation[1],-translation[2]);
  transform->Concatenate(rotation);
  transform->Translate(translation);
  tpdf->SetTransform(transform);

  crlVtkMeshIO::WriteMesh(tpdf->GetOutput(),argv[3]);

  return EXIT_SUCCESS;

}
