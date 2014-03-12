#ifndef __crlVtkMeshIO_h
#define __crlVtkMeshIO_h

#include <vtkPolyData.h>
#include <string>

class crlVtkMeshIO 
{
public:
  static bool WriteMesh(vtkPolyData *mesh, const std::string& fileName, bool binary = false);
	static vtkPolyData *ReadMesh(const std::string& fileName);

};

#endif
