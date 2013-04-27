#include "crlAnisotropicDiffusion.h"

int main(int argc, char *argv[])
{
  AnisotropicDiffusion a;

  if (argc != 5) {
    std::cout << 
      "Usage: crlAnisotropicDiffusion inimage iterations conductance outimage" 
             << std::endl;
    exit(1);
  }

  int iterations;
  ::sscanf(argv[2], "%d", &iterations);
  float conductance;
  ::sscanf(argv[3], "%f", &conductance);

  a.SetIterations(iterations);
  a.SetConductance(conductance);

  a.LoadImage(argv[1]);
  a.SmoothImage();
  a.WriteImage(argv[4]);
  
  exit(0); // success
}

