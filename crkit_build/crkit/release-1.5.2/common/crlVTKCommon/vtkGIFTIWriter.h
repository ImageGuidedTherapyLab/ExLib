#ifndef __vtkGIFTIWriter_h
#define __vtkGIFTIWriter_h

#include <vtkWriter.h>
#include <gifti_io.h>
#include <string>

class vtkGIFTIWriter : public vtkWriter
{
public:
  vtkTypeRevisionMacro(vtkGIFTIWriter,vtkWriter);

   static vtkGIFTIWriter *New();

  void SetFileName(std::string name, bool writeBinary = false);

  virtual void SetInput(vtkDataObject *input);
  virtual void SetWriteScalars(bool write);
  virtual void WriteData();
  virtual int Write();

protected:
  vtkGIFTIWriter();
  ~vtkGIFTIWriter();

private:
  vtkGIFTIWriter(const vtkGIFTIWriter&);  // Not implemented.
  void operator=(const vtkGIFTIWriter&);  // Not implemented.

  std::string fileName;
  std::string fileNameWithoutExt;
  bool fileNameSet;
  bool binary;
  gifti_image *gim;
  gifti_image *texgim;
  gifti_image *linegim;
  int scalarComponents;
  bool writeScalars;
  long numTriangles, numLines;
};

#endif
