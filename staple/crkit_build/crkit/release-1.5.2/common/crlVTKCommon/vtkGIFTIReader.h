#ifndef __vtkGIFTIReader_h
#define __vtkGIFTIReader_h

#include <vtkPolyDataAlgorithm.h>
#include <string>

class vtkGIFTIReader : public vtkPolyDataAlgorithm
{
public:
  vtkTypeRevisionMacro(vtkGIFTIReader,vtkPolyDataAlgorithm);

   static vtkGIFTIReader *New();

  void SetFileName(std::string name);
  virtual void SetFileName(const char* name);


protected:
  vtkGIFTIReader();
  ~vtkGIFTIReader();

  virtual int RequestData(vtkInformation*,vtkInformationVector**,vtkInformationVector*);

private:
  vtkGIFTIReader(const vtkGIFTIReader&);  // Not implemented.
  void operator=(const vtkGIFTIReader&);  // Not implemented.

  std::string fileName;
  bool fileNameSet;
};

#endif
