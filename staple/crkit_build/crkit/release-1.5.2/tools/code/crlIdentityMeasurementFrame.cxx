#include <tclap/CmdLine.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <itkAffineTransform.h>

int main( int argc, char * argv[] )
{
  std::string inputImageFile;
  std::string outputImageFile;

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',"1.0");
      TCLAP::UnlabeledValueArg<std::string> inputImageFileArg("inputImageFile","has to be a 4D DWI image.",true,"","input DWI image file",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputImageFileArg("outputImageFile","result is written here.",true,"","output DWI image file",cmd);

      cmd.parse(argc,argv);

      if (inputImageFileArg.isSet()) inputImageFile = inputImageFileArg.getValue();
      if (outputImageFileArg.isSet()) outputImageFile = outputImageFileArg.getValue();
    }
  catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(1);
    }

  std::ifstream reader(inputImageFile.c_str());
  std::ofstream writer(outputImageFile.c_str());

  typedef itk::AffineTransform< double, 3 >  AffineTransformType;

  itk::Matrix<double,3,3> mf;
  AffineTransformType::Pointer mfTransform = AffineTransformType::New();
  AffineTransformType::Pointer inverse = AffineTransformType::New();
  mf.SetIdentity();

  std::string line;
  while (!reader.eof())
    {
      getline(reader,line);

      if (line.substr(0,17) == "measurement frame")
	{
	  sscanf(line.c_str(),"measurement frame: (%lf,%lf,%lf) (%lf,%lf,%lf) (%lf,%lf,%lf)",&mf[0][0],&mf[0][1],&mf[0][2],&mf[1][0],&mf[1][1],&mf[1][2],&mf[2][0],&mf[2][1],&mf[2][2]);
	  mfTransform->SetMatrix(mf);
	  mfTransform->GetInverse(inverse);

	  writer << "measurement frame: (1,0,0) (0,1,0) (0,0,1)" << std::endl;
	}
      else if (line.substr(0,14) == "DWMRI_gradient")
	{
	  int k;
	  itk::Vector<double,3> g;
	  std::stringstream ss;

	  sscanf(line.c_str(),"DWMRI_gradient_%d:= %lf %lf %lf",&k,&g[0],&g[1],&g[2]);

	  ss << "DWMRI_gradient_" << std::right << std::setw(4) << std::setfill('0') << k << ":= ";
	  writer << ss.str();

	  g = inverse->TransformVector(g);
	  writer << g[0] <<  " " << g[1] << " " << g[2] << std::endl;
	}
      else
	if (line != "")
	  writer << line << std::endl;
    }

  reader.close();
  writer.close();
}
