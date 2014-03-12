
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkNeighborhoodAlgorithm.h"


int main( int argc, char ** argv )
{
  if ( argc != 7 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << " 4dpdfInputImage classInteractionMatrixFile "
      << " convergenceTolerance mrfStrength maxIterationCount outputImageFile"
      << std::endl;
    std::cerr << "This is a work in progress." << std::endl;
      return -1;
    }

  typedef double PixelType;
  static int const ImageDimension = 4;  // needs to be 1 + spatial dimensions.
  static int const ImageDimensionMinusOne = ImageDimension - 1; 

  typedef itk::Image< PixelType, ImageDimension >  ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
  typedef itk::ImageLinearIteratorWithIndex< ImageType>        IteratorType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  typedef itk::NumericTraits< PixelType >::AccumulateType SumType;
  typedef itk::NumericTraits< SumType >::RealType MeanType;
  typedef ImageType::IndexType Index4DType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  double mrfboost = std::atof( argv[4] );
  signed int maxIterationCount = std::atoi( argv[5] );
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[6] );

  try
    {
    reader->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }

  ImageType::Pointer output = ImageType::New();
  output->SetRegions(reader->GetOutput()->GetLargestPossibleRegion());
  output->CopyInformation(reader->GetOutput());
  output->Allocate();

  ImageType::Pointer prevEstimate = ImageType::New();
  ImageType::Pointer curnEstimate = ImageType::New();
  ImageType::Pointer tmpImagePointer = ImageType::New();

  prevEstimate = reader->GetOutput();
  curnEstimate = output;

  ImageType::SizeType inSize;
  inSize = prevEstimate->GetLargestPossibleRegion().GetSize();

  unsigned int numberOfProbabilitiesStored = 1;
  for (unsigned int i = 0; i <= ImageDimensionMinusOne; i++) {
    numberOfProbabilitiesStored *= inSize[i];
  }
  std::cout << "Automatic initialiation of the interaction matrix " <<
       " with size " << 
                (inSize[ImageDimensionMinusOne] * 
                 inSize[ImageDimensionMinusOne]) << std::endl;
#ifdef WIN32
	double* interactionmatrix;
	interactionmatrix = (double*)malloc(sizeof(double)*inSize[ImageDimensionMinusOne]*inSize[ImageDimensionMinusOne]);
	assert(interactionmatrix);
#else
	double interactionmatrix[inSize[ImageDimensionMinusOne]*
                           inSize[ImageDimensionMinusOne]];
#endif
  for (unsigned int i = 0; i < inSize[ImageDimensionMinusOne]; i++) {
    for (unsigned int j = 0; j < inSize[ImageDimensionMinusOne]; j++) {
      if (i == j) {
        interactionmatrix[i*inSize[ImageDimensionMinusOne] + j] = 0.95;
      } else {
        interactionmatrix[i*inSize[ImageDimensionMinusOne] + j] = 
                                    0.05/(inSize[ImageDimensionMinusOne]-1);
      }
      std::cout << "intmat["<<i<<","<<j<<"] " << 
        interactionmatrix[i*inSize[ImageDimensionMinusOne] + j] << std::endl;
    }
  }

  // Repeatedly iterate over the calculation until we run out of iterations.
  do {

  IteratorType inIt( prevEstimate, prevEstimate->GetLargestPossibleRegion());

  // Lines will walk along the last axis.
  inIt.SetDirection( ImageDimensionMinusOne ); 
  inIt.GoToBegin();

  while( !inIt.IsAtEnd() ) { // Iterating over input space.
    inIt.GoToBeginOfLine();
    SumType sum = itk::NumericTraits< SumType >::Zero;
    while ( !inIt.IsAtEndOfLine() ) { // Iterating over space over line.
      // We initialize the images with the 
      //   product of the site prior and the likelihood.
      double tmp = inIt.Get();
      ImageType::IndexType inIdx = inIt.GetIndex();

      // inIt is going to iterate over all the voxels in 4D, 
      // and down each line in 4D, so s gets consecutively accessed.
      unsigned int s       = inIdx[ImageDimensionMinusOne];
      unsigned int numcomp = inSize[ImageDimensionMinusOne];
      double interaction = 0.0;
      for (int n = 0; n < inSize[ImageDimensionMinusOne]; n++) {
        ImageType::IndexType neighborIdx = inIt.GetIndex();
        // Hard coded for a particular neighborhood size, and for 3D.
        for (int rk = -1; rk <= 1; rk++) {
          for (int rj = -1; rj <= 1; rj++) {
            for (int ri = -1; ri <= 1; ri++) {
              if ( (inIdx[2] == 0) && (rk == -1) ) continue;
              if ( (inIdx[2] == (inSize[2]-1)) && (rk == +1) ) continue;
              if ( (inIdx[1] == 0) && (rj == -1) ) continue;
              if ( (inIdx[1] == (inSize[1]-1)) && (rj == +1) ) continue;
              if ( (inIdx[0] == 0) && (ri == -1) ) continue;
              if ( (inIdx[0] == (inSize[0]-1)) && (ri == +1) ) continue;
              neighborIdx[0] = inIdx[0] + ri;
              neighborIdx[1] = inIdx[1] + rj;
              neighborIdx[2] = inIdx[2] + rk;
              neighborIdx[3] = n;
              // mrfboost plays the role of inverse temperature, and modifies
              // the strength of the spatial prior.
              tmp *= exp( log(interactionmatrix[s*numcomp + n])*mrfboost * 
                    prevEstimate->GetPixel(neighborIdx) );
            }
          }
        }
      }
      curnEstimate->SetPixel( inIdx, tmp );
      sum += tmp; // Keep track of the sum of new values along this line.
      ++inIt; // move s onto the next label value in this spatial location.
    }
    
    inIt.GoToBeginOfLine(); // Now iterate over the same line, in order to
                            // compute the necessary normalization.
    if (sum <= 0.0) {
      std::cerr << "Found a sum of " << sum << " at " 
                << inIt.GetIndex() << std::endl;
      sum = 1.0;
    }
    while ( !inIt.IsAtEndOfLine() ) { // Iterating over space over line.
      ImageType::IndexType inIdx = inIt.GetIndex();
      curnEstimate->SetPixel( inIdx, curnEstimate->GetPixel( inIdx ) / sum );
      ++inIt;
    }
    inIt.NextLine(); // Move on to the next location in space.
  }

    // That completes the iteration. Let's prepare for the next.
    tmpImagePointer = prevEstimate;
    prevEstimate = curnEstimate;
    curnEstimate = tmpImagePointer;

    maxIterationCount--;
  } while (maxIterationCount >= 0);

#ifdef WIN32
	free(interactionmatrix);
#endif

  // Now prevEstimate holds the most up to date estimate of the mean field

  writer->SetInput( prevEstimate );
  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }

  return EXIT_SUCCESS;

}
