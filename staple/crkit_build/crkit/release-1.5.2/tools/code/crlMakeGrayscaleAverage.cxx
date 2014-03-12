#include <tclap/CmdLine.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkNormalizeImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkExceptionObject.h>
#include <itkStatisticsImageFilter.h>
#include <itkShiftScaleImageFilter.h>

int main(int argc, char * argv[])
  {
  // TCLAP command line parsing
  //
  // usage: crlMakeGrayscaleAverage input1 [input2...] --out output.nrrd
  //
  // force use of an output command line argument --out to avoid the situation where someone does
  // crlMakeGrayscaleAverage *.nrrd 
  // and overwrites an input accidently.
  //
  TCLAP::CmdLine cmd("Computational Radiology Laboratory", ' ', "1.0");
  TCLAP::ValueArg<std::string> outputImageFileArg("o", "out",
      "may be in any ITK format.", true, "", "output image", cmd);
  TCLAP::UnlabeledMultiArg<std::string> inputImageFilesArg("in",
      "may be in any ITK format.", true, "input image", cmd);
  TCLAP::ValueArg<float> desiredMeanArg("m","mean","float value for output mean",false,500.0,"output mean",cmd);
  TCLAP::ValueArg<float> desiredVarArg("v","var","float value for desired var",false,100.0,"output variance",cmd);

  try
    {
    cmd.parse(argc,argv);
    }
  catch (TCLAP::ArgException& e)
    {
    std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
    exit(1);
    }

  const std::vector<std::string>& inputFiles = inputImageFilesArg.getValue();
  const std::string& outputFile = outputImageFileArg.getValue();
  float mean=desiredMeanArg.getValue();
  float var=desiredVarArg.getValue();

  std::cout << "Desired mean=" << mean << std::endl;
  std::cout << "Desired var=" << var << std::endl;
  std::cout << "Output file=" << outputFile << std::endl;
  for (unsigned int i = 0; i < inputFiles.size(); i++)
    {
    std::cout << "Input file[" << i << "]=" << inputFiles[i] << std::endl;
    }


  // itk typedefs
  typedef itk::Image<float,3>			ImageType;
  typedef itk::ImageFileReader<ImageType>	ImageReaderType;
  typedef itk::ImageFileWriter<ImageType>	ImageWriterType;
  typedef itk::NormalizeImageFilter<ImageType,ImageType>	NormalizerType;
  typedef itk::AddImageFilter<ImageType,ImageType,ImageType> 	AdderType;
  typedef itk::StatisticsImageFilter<ImageType> StatisticsFilterType;
  typedef itk::ShiftScaleImageFilter<ImageType,ImageType> ScalerType;

  // read in first image file and detach as our template for accumulating.
  ImageReaderType::Pointer reader = ImageReaderType::New();
  NormalizerType::Pointer normalizer = NormalizerType::New();
  
  reader->SetFileName(inputFiles[0]);
  normalizer->SetInput(reader->GetOutput());
  
  try
    {
    normalizer->Update();
    }
  catch (itk::ExceptionObject &err)
    {
    std::cerr << "Failed reading first file: " << inputFiles[0] << " error: "
    	<< err << std::endl;
    return EXIT_FAILURE;
    }
  
  ImageType::Pointer output = normalizer->GetOutput();
  output->DisconnectPipeline();
  output->FillBuffer(0.0);	// wasteful, but whatevs
  
  // read in each subsequent image file, normalize, and add to our template
  AdderType::Pointer adder = AdderType::New();
  StatisticsFilterType::Pointer stats = StatisticsFilterType::New();
  ScalerType::Pointer scaler = ScalerType::New();
  
  
  float n = static_cast<float>(inputFiles.size());
  for (unsigned int i = 0; i < inputFiles.size(); i++) 
    {
    reader->SetFileName(inputFiles[i]);
    stats->SetInput(reader->GetOutput());
    stats->Update();
//    std::cerr << "mean=" << stats->GetMean() << std::endl;
//    std::cerr << "variance=" << stats->GetVariance() << std::endl;
    
    scaler->SetInput(reader->GetOutput());
    scaler->SetShift(-stats->GetMean());
    scaler->SetScale(1/sqrt(stats->GetVariance()));
    scaler->Update();
    
//    stats->SetInput(scaler->GetOutput());
//    stats->Update();
//    std::cerr << "mean=" << stats->GetMean() << std::endl;
//    std::cerr << "variance=" << stats->GetVariance() << std::endl;
    adder->SetInput1(output);
    adder->SetInput2(scaler->GetOutput());
    try
      {
      std::cout << "reading " << inputFiles[i] << std::endl;
      adder->Update();
      }
    catch(itk::ExceptionObject &err)
      {
      std::cerr << "Failed reading file: " << inputFiles[i] << " error: " 
      	<< err << std::endl;
      return EXIT_FAILURE;
      }
    output = adder->GetOutput();
    output->DisconnectPipeline();
    }
  
  // scale the output
  stats->SetInput(output);
  stats->Update();
  
  scaler->SetInput(output);
  scaler->SetShift(-stats->GetMean()+mean);
  scaler->SetScale(sqrt(var)/sqrt(stats->GetVariance()));
  scaler->Update();
  
  // write out the output
  
  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName(outputFile);
  writer->SetInput(scaler->GetOutput());
  
  try
    {
    std::cout << "writing " << outputFile << std::endl;
    writer->Update();
    }
  catch(itk::ExceptionObject& err)
    {
    std::cerr << "Failed writing output file: " << outputFile << " error: " << err << std::endl;
    return EXIT_FAILURE;
    }
 
  
  return EXIT_SUCCESS;	  
  }
