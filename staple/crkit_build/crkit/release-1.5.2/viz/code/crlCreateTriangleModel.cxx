//
// Create a triangle model from a segmentation
//

#include <vtkVolume16Reader.h>
#include <vtkStripper.h>
#include <vtkPolyDataNormals.h>
#include <vtkContourFilter.h>
#include <vtkImageData.h>
#include <vtkDecimatePro.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkDiscreteMarchingCubes.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkMatrix4x4.h>

#include "ImageTypeDefinitions.h"
#include <crlVtkMeshIO.h>

int main(int argc, char *argv[])
{
  
  if (argc != 4) {
    std::cout << "Usage: " << argv[0] << " segmentedImageFile "
	      << "threshold triangleFileName" << std::endl;
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
  ImageTypeDefinitions::ImageType::DirectionType direction = inputImage->GetDirection();
  
  ImageTypeDefinitions::ITK2VTKConnectorFilterType::Pointer
    ITK2VTKConnector =
    ImageTypeDefinitions::ITK2VTKConnectorFilterType::New();
  ITK2VTKConnector->SetInput( inputImage );
  
  // Set up a pipeline to extract a contour from the data
  vtkContourFilter *surfaceExtractor = vtkContourFilter::New();
  surfaceExtractor->SetInput( ITK2VTKConnector->GetOutput() );
  surfaceExtractor->SetValue(0, threshold );
  
  vtkDecimatePro *surfaceDecimation = vtkDecimatePro::New();
  surfaceDecimation->SetInputConnection( surfaceExtractor->GetOutputPort() );
  surfaceDecimation->SetTargetReduction( 0.1 );
  surfaceDecimation->PreserveTopologyOn( );
  
  vtkWindowedSincPolyDataFilter *smoother =vtkWindowedSincPolyDataFilter::New();
  smoother->SetInputConnection( surfaceDecimation->GetOutputPort() );
  // Smooth along sharp interior edges
  smoother->FeatureEdgeSmoothingOn();
  smoother->SetFeatureAngle( 60.0 );
  smoother->SetNumberOfIterations( 10 );
  
  vtkPolyDataNormals *surfaceNormals = vtkPolyDataNormals::New();
  surfaceNormals->SetInputConnection( smoother->GetOutputPort() );
  surfaceNormals->SetFeatureAngle( 60.0 );
  
  vtkStripper *surfaceStripper = vtkStripper::New();
  surfaceStripper->SetInputConnection( surfaceNormals->GetOutputPort() );
  
  //  vtkXMLPolyDataWriter *surfaceWriter = vtkXMLPolyDataWriter::New();
  
  // take care of rotations of the input image that were ignored by VTK
  double translation[3];
  vtkMatrix4x4 *rotation = vtkMatrix4x4::New();
  
  ImageTypeDefinitions::ImageType::IndexType idx;
  ImageTypeDefinitions::ImageType::PointType point;
  idx[0]=0; idx[1]=0; idx[2]=0;
  inputImage->TransformIndexToPhysicalPoint(idx,point);
  
  translation[0] = point[0];
  translation[1] = point[1];
  translation[2] = point[2];
  
  point = inputImage->GetOrigin();
  
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
  tpdf->SetInputConnection(surfaceStripper->GetOutputPort());
  transform->PostMultiply();
  transform->Translate(-translation[0],-translation[1],-translation[2]);
  transform->Concatenate(rotation);
  transform->Translate(translation);
  tpdf->SetTransform(transform);
  
  crlVtkMeshIO::WriteMesh(tpdf->GetOutput(),argv[3]);
  
  surfaceExtractor->Delete();
  surfaceDecimation->Delete();
  surfaceNormals->Delete();
  surfaceStripper->Delete();
  tpdf->Delete();
  
  return EXIT_SUCCESS;  
}
