#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>

#include <itkMinimumMaximumImageCalculator.h>
#include "ImageTypeDefinitions.h"

#include <tclap/CmdLine.h>
#include "configuration.h"
#include <crlVtkMeshIO.h>

int main (int argc, char *argv[])
{
  TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',CRKIT_VERSION_STRING);
  TCLAP::UnlabeledValueArg<std::string> inputImageArg("input","Input file name",true,"","Input file name",cmd);
  TCLAP::ValueArg<unsigned int> binArg("b","bins","An integer describing the number of bins to be generated",false,0,"Number of bins",cmd);
  TCLAP::ValueArg<double> sizeArg("s","size","An integer describing the bin size",false,0.0,"Bin size",cmd);
  TCLAP::SwitchArg modelArg("m","model","If this switch is set, the input refers to a VTK model instead of an ITK image.",cmd,false);

  cmd.parse(argc,argv);

  unsigned int bins = binArg.getValue();
  double size = sizeArg.getValue();
  if ((bins == 0 && size == 0))
    {
      std::cout << "Either bin size or number of bins must be set." << std::endl;
      return EXIT_FAILURE;
    }

  std::vector<int> frequency;
  double min, max;

  if (modelArg.isSet())
    {
      vtkPolyData *pd = crlVtkMeshIO::ReadMesh(inputImageArg.getValue());
      vtkFloatArray *scalars = dynamic_cast<vtkFloatArray*>(pd->GetPointData()->GetScalars());

      min = pd->GetScalarRange()[0];
      max = pd->GetScalarRange()[1];
      if (size == 0)
	size = (max-min) / bins;
      else
	bins = (int)ceil((max-min)/size);

      std::cout << "Scalar range: " << min << " to " << max << std::endl;
      std::cout << "Bin size: " << size << std::endl;

      for (unsigned int b = 0; b < bins; b++)
	frequency.push_back(0);

      float *s = scalars->GetPointer(0);
      long size = scalars->GetNumberOfTuples();

      for (long n = 0; n < size; n++)
	{
	  frequency[(int)floor((*s-min)/size)]++;
	  s++;
	}
      pd->Delete();
    }
  else
    {
      ImageTypeDefinitions::ReaderType::Pointer reader = ImageTypeDefinitions::ReaderType::New();
      itk::MinimumMaximumImageCalculator<ImageTypeDefinitions::ImageType>::Pointer calc = itk::MinimumMaximumImageCalculator<ImageTypeDefinitions::ImageType>::New();

      reader->SetFileName(inputImageArg.getValue().c_str());
      reader->Update();

      calc->SetImage(reader->GetOutput());
      calc->Compute();
      min = calc->GetMinimum();
      max = calc->GetMaximum();
      if (size == 0)
	size = (max-min) / bins;
      else
	bins = (int)ceil((max-min)/size);

      std::cout << "Scalar range: " << min << " to " << max << std::endl;
      std::cout << "Bin size: " << size << std::endl;


      for (unsigned int b = 0; b < bins; b++)
	frequency.push_back(0);

      ImageTypeDefinitions::InternalPixelType *s = reader->GetOutput()->GetBufferPointer();
      unsigned long imageSize = reader->GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels();

      for (unsigned long n = 0; n < imageSize; n++)
	{
	  frequency[(int)floor((*s-min)/size)]++;
	  s++;
	}
    }

  for (unsigned int b = 0; b < bins; b++)
    {
      std::cout << setw(12) << std::left;
      std::cout << min+(b*size) << ": ";
      std::cout << frequency[b] << std::endl;
    }

  return EXIT_SUCCESS;
}
