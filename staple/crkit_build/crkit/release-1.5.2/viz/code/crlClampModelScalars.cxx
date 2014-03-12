// Clamp scalars in a VTK model

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>

#include <crlVtkMeshIO.h>

int main(int argc, char *argv[])
{
  if (argc != 5)
    {
      std::cout << "Usage: " << argv[0] << " inputModelFile lowerThreshold upperThreshold outputModelFile" << std::endl;
	return EXIT_FAILURE;
    }

  float min = atof(argv[2]);
  float max = atof(argv[3]);

  if (min > max)
    {
      std::cout << "Error: lowerThreshold has to be <= upperThreshold" << std::endl;
      return EXIT_FAILURE;
    }

  vtkPolyData *pd = crlVtkMeshIO::ReadMesh(argv[1]);
  vtkFloatArray *scalars = dynamic_cast<vtkFloatArray*>(pd->GetPointData()->GetScalars());

  if (scalars->GetNumberOfComponents() != 1)
    {
      std::cout << " Error: Scalar data cannot have more than one component." << std::endl;
      return EXIT_FAILURE;
    }

  float *s = scalars->GetPointer(0);

  for (long n = 0; n < scalars->GetNumberOfTuples(); n++)
    {
      if (*s < min)
	*s = min;
      else if (*s > max)
	*s = max;
      s++;
    }

  crlVtkMeshIO::WriteMesh(pd,argv[4]);

  return EXIT_SUCCESS;
}
