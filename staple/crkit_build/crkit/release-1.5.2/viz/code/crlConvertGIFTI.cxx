#include <tclap/CmdLine.h>
#include "configuration.h"
#include <crlVtkMeshIO.h>

int main (int argc, char *argv[])
{
  TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',CRKIT_VERSION_STRING);
  TCLAP::UnlabeledValueArg<std::string> inputFileArg("input","Input file name",true,"","Input file name",cmd);
  TCLAP::UnlabeledValueArg<std::string> outputFileArg("output","Output file name",true,"","Output file name",cmd);
  TCLAP::SwitchArg binaryArg("e","external","Create external binary files",cmd,false);

  cmd.parse(argc,argv);

  vtkPolyData *pd = crlVtkMeshIO::ReadMesh(inputFileArg.getValue());
  crlVtkMeshIO::WriteMesh(pd,outputFileArg.getValue(), binaryArg.isSet());

  pd->Delete();
  return EXIT_SUCCESS;
}
