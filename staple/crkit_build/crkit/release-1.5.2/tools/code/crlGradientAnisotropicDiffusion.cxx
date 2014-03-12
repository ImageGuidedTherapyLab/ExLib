#include "crlGradientAnisotropicDiffusion.h"
#include "tclap/CmdLine.h"
#include "configuration.h"



int main(int argc, char* argv[])
{
  try
    {
    // TCLAP command line parsing
    TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',CRKIT_VERSION_STRING );
    
    TCLAP::UnlabeledValueArg<std::string> inputImageFileArg("inputImage", "Filename of the (vector) image to diffuse.", true, "", "input image",cmd);
    TCLAP::ValueArg<std::string> outputImageFileArg("o", "outputImage", "Filename to write the output (vector) image to.", true, "", "output image", cmd);
    TCLAP::ValueArg<double> conductanceArg("c", "conductance", "conductance value", false, 1.0, "conductance",cmd);
    TCLAP::ValueArg<unsigned int> iterationsArg("i","iterations","number of iterations to run", false, 5, "iterations", cmd);
    TCLAP::ValueArg<double> timestepArg("t","timestep","timestep to use", false, 0.05, "timestep", cmd);
    TCLAP::ValueArg<unsigned int> threadsArg("T","threads","number of threads to use (note: horribly inefficient)", false, 1, "num threads", cmd);
    
    cmd.parse(argc,argv);
    
    typedef crl::GradientAnisotropicDiffusion DiffusionType;
    
    DiffusionType* diffusion = DiffusionType::ProbeFile(inputImageFileArg.getValue());
    
    diffusion->ReadFile(inputImageFileArg.getValue());
    diffusion->SetConductanceParameter(conductanceArg.getValue());
    diffusion->SetIterations(iterationsArg.getValue());
    diffusion->SetTimeStep(timestepArg.getValue());
    diffusion->SetThreads(threadsArg.getValue());
    diffusion->Execute();
    diffusion->WriteFile(outputImageFileArg.getValue());
    
    delete diffusion;
    
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
