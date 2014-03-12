#include <vtkDecimatePro.h>
#include <vtkPolyDataNormals.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkStripper.h>
#include <vtkTriangleFilter.h>
#include <crlVtkMeshIO.h>

int main(int argc, char *argv[])
{
  
  if ((argc < 3) || (argc > 4)) {
    std::cout << "Usage: " << argv[0] << " inputFileName outputFileName [reduction]" << std::endl;
    return EXIT_FAILURE;
  }
  
  vtkPolyData *surface = crlVtkMeshIO::ReadMesh(argv[1]);
  
  vtkTriangleFilter *triangulator = vtkTriangleFilter::New();
  triangulator->SetInput(surface);
  triangulator->Update();
  
  vtkDecimatePro *surfaceDecimation = vtkDecimatePro::New();
  surfaceDecimation->SetInput(triangulator->GetOutput());
  if (argc == 4)
    surfaceDecimation->SetTargetReduction(atof(argv[3]));
  else
    surfaceDecimation->SetTargetReduction( 0.1 );
  if (argc == 3)
    surfaceDecimation->PreserveTopologyOn( );
  surfaceDecimation->Update();
  
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
  
  crlVtkMeshIO::WriteMesh(surfaceStripper->GetOutput(),argv[2]);
  
  surfaceDecimation->Delete();
  surfaceNormals->Delete();
  surfaceStripper->Delete();
}
