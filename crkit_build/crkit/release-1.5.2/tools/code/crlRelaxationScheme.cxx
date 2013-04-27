#include "crlRelaxationScheme.h"
#include "tclap/CmdLine.h"
#include "configuration.h"



int main(int argc, char* argv[])
{
  try
    {
    // TCLAP command line parsing
    TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',CRKIT_VERSION_STRING );
    
    TCLAP::UnlabeledValueArg<std::string> inputImageFileArg("inputImage", "Filename of the vector-format segmentation to operate on.", true, "", "input image",cmd);
    TCLAP::ValueArg<std::string> outputImageFileArg("o", "outputImage", "Filename to write the output vector-format segmentation to.", true, "", "output image", cmd);
    TCLAP::ValueArg<unsigned int> iterationsArg("i","iterations","number of iterations to run", false, 5, "iterations", cmd);
    TCLAP::ValueArg<std::string> debugImageArg("d","debugImage","Filename to write a debugging label image to showing affected voxels.", false, "", "debug image", cmd);
    TCLAP::ValueArg<std::string> exampleImageArg("e","exampleImage", "Filename of example segmentation to compute statistics from", false, "", "segmentation", cmd);
    
    cmd.parse(argc,argv);
    
    typedef crl::RelaxationScheme RelaxationType;
    
    RelaxationType* relaxation = RelaxationType::ProbeFile(inputImageFileArg.getValue());

    relaxation->ReadFile(inputImageFileArg.getValue());
    relaxation->SetIterations(iterationsArg.getValue());

    if ( exampleImageArg.isSet() )
      {
      std::cout << "Example Image mode" << std::endl;
      relaxation->ReadExampleImage( exampleImageArg.getValue() );
      relaxation->Execute1();
      }
    else
      {
      relaxation->Execute();
      }

    if ( debugImageArg.isSet() )
	relaxation->WriteFile(outputImageFileArg.getValue(), debugImageArg.getValue());
    else
	relaxation->WriteFile(outputImageFileArg.getValue() );
    
    delete relaxation;
    
    }
  catch (TCLAP::ArgException& e)
    {
    std::cerr << "error parsing " << e.argId() << ": " << e.error() << std::endl;
    exit(1);
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "itk error: " << e << std::endl;
    exit(1);
    }
  catch (std::runtime_error& e)
    {
    std::cerr << "ERROR: " << e.what() << std::endl;
    exit(1);
    }
  
  
  exit(0);
}
