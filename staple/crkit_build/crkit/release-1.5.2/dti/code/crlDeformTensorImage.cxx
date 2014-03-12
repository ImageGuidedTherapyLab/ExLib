
#include "TensorDeformableResampler.h"

int main(int argc, const char *argv[])
{

  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " movingTensorImageFile " <<
      "displacementVectorField " << "outputTensorImageField" << std::endl;
    std::cerr << "  \t This program resamples a tensor image. The image is resampled using the \n\t same spacings and size as the displacement field image.\n \t The image is warped by first estimating a local 'rotation field', then\n \t rotating and shifting the tensors." << std::endl;     
    exit(1);
  }
  
  TensorDeformableResampler r;
  std::string  movingImageFile( argv[1] );
  std::string  displacementFieldImageFile( argv[2] );
  std::string  outputImageFile( argv[3] );

  std::string  interpolationType("linear");
  float defaultRecursiveGaussianSmoothingBandwidth = 1.0;
      
  // set interpolation type
  if (strncmp(interpolationType.c_str(), "linear",strlen("linear")) == 0) {
    r.SetInterpolationToLinear();
  } else if (strncmp(interpolationType.c_str(), "1",strlen("1")) == 0) {
    r.SetInterpolationToLinear();
  } else if (strncmp(interpolationType.c_str(), "nearestneighbor",strlen("nearest")) == 0) {
    r.SetInterpolationToNearestNeighbor();
  } else if (strncmp(interpolationType.c_str(), "2",strlen("2")) == 0) {
    r.SetInterpolationToNearestNeighbor();
  }

  /** Load displacement field */
  r.itkLoadDisplacementField(displacementFieldImageFile.c_str());
  std::cout << "loaded Displacement Field " << std::endl;

  /** Load moving image (tensor image) */
  r.itkLoadMovingImage(movingImageFile.c_str());  // load image 
  std::cout << "loaded moving image" << std::endl;

  /** doing the processing: rotate tensors, resample, write the output */
  r.RotateTensors( defaultRecursiveGaussianSmoothingBandwidth );
  std::cout << "finished rotate tensors" << std::endl;
      
  r.WarpMovingImage();
  std::cout << "finished warping image" << std::endl;

  std::cout << "writing image ..." << std::endl;
  r.itkWriteMovingImage(outputImageFile.c_str());

  exit(0); // success
}

