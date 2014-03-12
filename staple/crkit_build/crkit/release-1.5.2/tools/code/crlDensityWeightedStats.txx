
#ifndef _CRL_DENSITYWEIGHTEDSTATS_TXX_INCLUDED
#define _CRL_DENSITYWEIGHTEDSTATS_TXX_INCLUDED 1

template < const unsigned int TImageDimension, class TPixelType >
int DensityWeightedStats< TImageDimension, TPixelType > 
::Execute()
{
  typename ImageReaderType::Pointer readerDensity = ImageReaderType::New();
  typename ImageReaderType::Pointer readerScalar = ImageReaderType::New();

  if (m_DensityImageFileName.size() <= 0) {
    std::cerr << "Input 1 file name must be set." << std::endl;
  }
  if (m_ScalarImageFileName.size() <= 0) {
    std::cerr << "Input 2 file name must be set." << std::endl;
  }
  readerDensity->SetFileName( m_DensityImageFileName );
  readerScalar->SetFileName( m_ScalarImageFileName );

  try {
    readerDensity->Update();
    readerScalar->Update();
  } catch (itk::ExceptionObject &ex) {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;

  /* Now iterator over the images and compute the weighted mean and variance.
   */
  IteratorType itd( readerDensity->GetOutput(),
                    readerDensity->GetOutput()->GetRequestedRegion() );
  IteratorType its( readerScalar->GetOutput(), 
                    readerScalar->GetOutput()->GetRequestedRegion());

  double mean = 0.0;
  double totalweight = 0.0;
  // This assumes the same iterators go over the same geometry the same way.
  for (itd.GoToBegin(), its.GoToBegin(); !itd.IsAtEnd(), !its.IsAtEnd(); 
                                          ++its, ++itd) {
      TPixelType density = itd.Value();
      TPixelType scalar = its.Value();
      mean += density * scalar;
      totalweight += density;
  }
  if (totalweight != 0.0) {
    mean /= totalweight;
  }

  // Now loop over the data again to measure the weighted variance.
  double weightsquared = 0.0;
  double variance = 0.0;
  for (itd.GoToBegin(), its.GoToBegin(); !itd.IsAtEnd(), !its.IsAtEnd(); 
                                          ++its, ++itd) {
      TPixelType density = itd.Value();
      TPixelType scalar = its.Value();
      density /= totalweight; // Normalize the weights.
      variance += density*(scalar - mean)*(scalar - mean);
      weightsquared += (density*density);
  }
  if (weightsquared != 1.0) {
    variance /= (1.0 - weightsquared);
  }

  std::cout << "Mean: " << mean << " , Variance: " << variance <<
             " , Std. Dev: " << sqrt(variance) << std::endl;

  return EXIT_SUCCESS;

};

#endif

