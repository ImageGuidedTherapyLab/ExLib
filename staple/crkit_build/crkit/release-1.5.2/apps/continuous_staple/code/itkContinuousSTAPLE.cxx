#include <iostream>
#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkContinuousSTAPLEImageFilter.h"
#include <itkTimeProbe.h>
#include "configuration.h"

using namespace std;

int main(int argc, char **argv)
{
  TCLAP::CmdLine cmd("Computational Radiology Laboratory", ' ',CRKIT_VERSION_STRING);

  TCLAP::ValueArg<std::string> lstArg("l","listname","File containing a list of images",true,"","image list",cmd);
  TCLAP::ValueArg<std::string> resArg("o","outputname","Result image",true,"","result image",cmd);
  TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask (highly recommended)",false,"","computation mask",cmd);
  TCLAP::ValueArg<std::string> resMaskArg("r","resmaskname","Outputs the computation mask",false,"","computation mask output",cmd);
  TCLAP::ValueArg<std::string> matlabArg("O","matlabname","Result parameters as a matlab file",false,"","result matlab file",cmd);

  TCLAP::ValueArg<float> biasArg("b","initialbias","Initial value of bias parameters (default: 0.0)",false,0.0,"initial bias",cmd);
  TCLAP::ValueArg<float> varArg("v","initialvar","Initial value of variance parameters (default: 2.0)",false,2.0,"initial variance",cmd);
  TCLAP::ValueArg<float> thrArg("t","relativethreshold","Relative stopping criterion for Staple (default : 1.0e-8)",false,1.0e-8,"stopping criterion",cmd);

  TCLAP::ValueArg<unsigned int> iterArg("i","iterations","Maximum number of iterations (default: 100)",false,100,"number of iterations",cmd);
  TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default : 1)",false,1,"number of threads",cmd);

  try
  {
    cmd.parse(argc,argv);
  }
  catch (TCLAP::ArgException& e)
  {
    std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
    return(1);
  }

  string lstName, resFile, maskName, matlabName;
  lstName = lstArg.getValue();
  resFile = resArg.getValue();
  maskName = maskArg.getValue();
  matlabName = matlabArg.getValue();

  unsigned int itMax = iterArg.getValue();

  float iniBias = biasArg.getValue();
  float iniVar = varArg.getValue();
  float relThr = thrArg.getValue();

  unsigned int nbProcs = nbpArg.getValue();

  unsigned int vdim = 6;
	
  itk::TimeProbe tmpTime;
  tmpTime.Start();

  typedef itk::ImageFileReader < itk::Image <unsigned char, 3> > itkMaskReader;
  typedef itk::ImageFileWriter < itk::Image <unsigned char, 3> > itkMaskWriter;

  typedef itk::VectorImage<double,3> LogTensorImageType;
  typedef itk::ImageFileReader <LogTensorImageType> itkInputReader;
  typedef itk::ImageFileWriter <LogTensorImageType> itkOutputWriter;
  typedef itk::ContinuousSTAPLEImageFilter StapleFilterType;
  
  StapleFilterType::Pointer contStapleFilter = StapleFilterType::New();
  contStapleFilter->SetMaximumIterations(itMax);
  contStapleFilter->SetOverrideComputationRegion(false);
  contStapleFilter->SetRelativeConvergenceThreshold(relThr);
  
  if (strcmp("",maskName.c_str()) != 0)
  {
    itkMaskReader::Pointer maskRead = itkMaskReader::New();
    maskRead->SetFileName(maskName.c_str());
    maskRead->Update();
    contStapleFilter->SetComputationMask(maskRead->GetOutput());
  }
  
  ifstream fileIn(lstName.c_str());
  
  int nbPats = 0;
  itkInputReader::Pointer imageReader;
  
  while (!fileIn.eof())
  {
    char tmpStr[2048];
    fileIn.getline(tmpStr,2048);
    
    if (strcmp(tmpStr,"") == 0)
      continue;
    
    std::cout << "Loading image " << nbPats << " " << tmpStr << "..." << std::endl;
    imageReader = itkInputReader::New();
    imageReader->SetFileName(tmpStr);
    imageReader->Update();
    
    contStapleFilter->SetInput(nbPats,imageReader->GetOutput());
    if (nbPats == 0)
    {
      vdim = imageReader->GetOutput()->GetVectorLength();
      contStapleFilter->SetNDim(vdim);
    }
    
    nbPats++;
  }
  
  fileIn.close();
  
  contStapleFilter->SetNumberOfThreads(nbProcs);
  contStapleFilter->InitializeExpertBias(iniBias);
  contStapleFilter->InitializeExpertCovariance(iniVar);
  
  contStapleFilter->Update();
  
  contStapleFilter->PrintPerformanceParameters(matlabName);
  
  itkOutputWriter::Pointer resultWriter = itkOutputWriter::New();
  resultWriter->SetFileName(resFile.c_str());
  resultWriter->SetUseCompression(true);
  resultWriter->SetInput(contStapleFilter->GetOutput());
  
  resultWriter->Update();
  
  if (strcmp(resMaskArg.getValue().c_str(),"") != 0)
  {
    itkMaskWriter::Pointer maskWr = itkMaskWriter::New();
    maskWr->SetInput(contStapleFilter->GetComputationMask());
    maskWr->SetFileName(resMaskArg.getValue());
    
    maskWr->Update();
  }

  tmpTime.Stop();

  cout << "Total computation time: " << tmpTime.GetTotal() << endl;

  return 0;
}
