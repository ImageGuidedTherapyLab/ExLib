
#include "ITKTensorResampler.h"

#include "itkMetaDataObject.h"

const char *get_argument(const char **argv,
                         const char *keystr,
                         const char *deflt)
{
  int i;
  for (i=1; argv[i]; i++)
    if (strstr(argv[i],keystr))
      return argv[i+1];
  return deflt;
}

int exists_argument(int argc,const char **argv,const char *keystr)
{
  int i;

  for (i=1;i<argc;i++) if (strncmp(argv[i],keystr,strlen(keystr)) == 0) return 1;
  return 0;
}

int main(int argc, const char *argv[])
{

  
  const bool help = exists_argument(argc, argv, "-h");

  const bool sevencomponent = exists_argument(argc, argv, "-7");
  const bool image_must_be_flipped = exists_argument(argc, argv, "-F");
  const int axis = atoi(get_argument(argv, "-F", "10"));
  
  if (help)
   {
    std::cerr << "argc is " << argc << std::endl;
    std::cerr << "Usage: crlTensorResampler "<< std::endl;
    std::cerr << "  geometryExampleImage movingImage transformFileName interpolationType outputImage" << std::endl;
    std::cerr << std::endl;
    std::cerr << "OR:" <<std::endl;
    std::cerr << "  geometryExampleImage" << std::endl; 
    std::cerr << "  component(0,0)-Image" << std::endl; 
    std::cerr << "  component(0,1)-Image "<< std::endl; 
    std::cerr << "  component(0,2)-Image"<< std::endl; 
    std::cerr << "  component(1,1)-Image"<< std::endl; 
    std::cerr << "  component(1,2)-Image "<< std::endl; 
    std::cerr << "  component(2,3)-Image"<< std::endl; 
    std::cerr << "  transformFileName "<< std::endl; 
    std::cerr << "  interpolationType"<< std::endl; 
    std::cerr << "  outputImage" << std::endl << std::endl;
    std::cerr << "Option:" <<std::endl;
    std::cerr << "  -7 \t (add this option if you want write a 7 component image to use ''tend expand'':" << std::endl << " \t  do:  tend expand -i tensors.nrrd -o tensors.vtk" << std::endl;     
    std::cerr << "  -F  <axis> \t (add this option if you want the image to be flipped before resampling," << std::endl << " \t axis=0,1,2 for x,y,z)" << std::endl;     
    
    exit(1);
  }
  
  TensorResampler r;
  std::string  outputImageFile;
  
/**One input (vector) image */
if (argc==6 || (argc==7 && sevencomponent) || (argc==8 && image_must_be_flipped) || (argc==9 && image_must_be_flipped && sevencomponent))
{
  std::string  geometryImageFile(argv[1]);
  std::string  movingImageFile(argv[2]);
  std::string  transformFile(argv[3]);
  outputImageFile = std::string(argv[5]);
  
  // set interpolation type
  if (strncmp(argv[4], "linear",strlen("linear")) == 0) {
    r.SetInterpolationToLinear();
  } else if (strncmp(argv[4], "1",strlen("1")) == 0) {
    r.SetInterpolationToLinear();
  } else if (strncmp(argv[4], "nearestneighbor",strlen("nearest")) == 0) {
    r.SetInterpolationToNearestNeighbor();
  } else if (strncmp(argv[4], "2",strlen("2")) == 0) {
    r.SetInterpolationToNearestNeighbor();
  }

  r.ParseParameterFile(transformFile.c_str());
//  std::cout << "parsed parameter file" << std::endl;
//r.PrintTransformParameters();

  r.LoadGeometryImage(geometryImageFile.c_str());  // load image using ITK methods
  
//  std::cout << "loaded Fixed Image" << std::endl;

    r.itkLoadMovingImage(movingImageFile.c_str());  // load image using ITK methods
//  std::cout << "loaded moving image" << std::endl;
}


/** Six input (scalar) images: the components of the symetric tensor */
else if (argc==11 || (argc==12 && sevencomponent) || (argc==13 && image_must_be_flipped) || (argc==14 && image_must_be_flipped && sevencomponent))
{
// std::cout << "Loading scalar images containing tensor components." << std::endl;
  std::string  geometryImageFile(argv[1]);
  std::string  movingImageFile1(argv[2]);
  std::string  movingImageFile2(argv[3]);
  std::string  movingImageFile3(argv[4]);
  std::string  movingImageFile4(argv[5]);
  std::string  movingImageFile5(argv[6]);
  std::string  movingImageFile6(argv[7]);
  std::string  transformFile(argv[8]);
  outputImageFile = std::string(argv[10]);
  
  // set interpolation type
  if (strncmp(argv[9], "linear",strlen("linear")) == 0) {
    r.SetInterpolationToLinear();
  } else if (strncmp(argv[9], "1",strlen("1")) == 0) {
    r.SetInterpolationToLinear();
  } else if (strncmp(argv[9], "nearestneighbor",strlen("nearest")) == 0) {
    r.SetInterpolationToNearestNeighbor();
  } else if (strncmp(argv[9], "2",strlen("2")) == 0) {
    r.SetInterpolationToNearestNeighbor();
  }
  r.ParseParameterFile(transformFile.c_str());
  std::cout << "parsed parameter file" << std::endl;
//r.PrintTransformParameters();

  r.itkLoadFixedImage(geometryImageFile.c_str());  // load image using ITK methods
  std::cout << "loaded Fixed Image" << std::endl;

  // load image using ITK methods
  r.itkLoadMovingImage(movingImageFile1.c_str(), movingImageFile2.c_str(), movingImageFile3.c_str(), movingImageFile4.c_str(), movingImageFile5.c_str(), movingImageFile6.c_str()); 
//  std::cout << "loaded moving image" << std::endl;
}
else {
    std::cerr << "Usage: crlTensorResampler " ;
    std::cerr << "  geometryExampleImage movingImage  transformFileName interpolationType outputImage" << std::endl;
    std::cerr << std::endl;
    std::cerr << "OR:" <<std::endl << " crlTensorResampler " ;
    std::cerr << "  geometryExampleImage" << std::endl; 
    std::cerr << "  component(0,0)-Image" << std::endl; 
    std::cerr << "  component(0,1)-Image "<< std::endl; 
    std::cerr << "  component(0,2)-Image"<< std::endl; 
    std::cerr << "  component(1,1)-Image"<< std::endl; 
    std::cerr << "  component(1,2)-Image "<< std::endl; 
    std::cerr << "  component(2,3)-Image"<< std::endl; 
    std::cerr << "  transformFileName "<< std::endl; 
    std::cerr << "  interpolationType"<< std::endl; 
    std::cerr << "  outputImage" << std::endl << std::endl;
    std::cerr << "Option:" <<std::endl;
    std::cerr << "  -7 \t (add this option if you want write a 7 component image to use ''tend expand'':" << std::endl << " \t  do:  tend expand -i tensors.nrrd -o tensors.vtk" << std::endl;     
    std::cerr << "  -F  <axis> \t (add this option if you want the image to be flipped before resampling," << std::endl << " \t axis=0,1,2 for x,y,z)" << std::endl;     
    
    exit(1);
  }




/** doing the processing, resample, rotate tensors, write the output */

if (image_must_be_flipped) {
  r.FlipImage(axis);
}

r.GenerateMovingImage();
// std::cout << "finished resampling image" << std::endl;

r.RotateTensors();
// std::cout << "finished rotate tensors" << std::endl;

// This aims to transfer the tensor measurement frame from the NRRD data
r.TransferMetaDataDictionary( );

// std::cout << "writing image" << std::endl;
if (sevencomponent)  {
  r.itkWriteMovingImage7component(outputImageFile.c_str());
  exit(0);
}
else {
  r.itkWriteMovingImage(outputImageFile.c_str());
  }

  exit(0); // success
}

