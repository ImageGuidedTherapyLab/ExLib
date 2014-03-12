
template <unsigned int Dimension, class PixelType>
int
crl::IndexOfMaxComponent<Dimension, PixelType>
::Execute()
{   

  typename ImageReaderType::Pointer r;
  typename ImageWriterType::Pointer w;

  r = ImageReaderType::New();
  r->SetFileName(m_InputFileName);
  try {
    r->Update(); // read in the file , could catch errors here...
  } catch ( itk::ExceptionObject &exp )
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
  }

  InSizeType inSize;
  inSize = r->GetOutput()->GetLargestPossibleRegion().GetSize();

  w = ImageWriterType::New();
  w->SetFileName( m_OutputFileName );

  typename OutImageType::IndexType outIndex;
  typename OutImageType::SizeType outSize;
  typename OutImageType::SpacingType outSpacing;
  typename InImageType::SpacingType inSpacing;
  inSpacing = r->GetOutput()->GetSpacing();

  typename InImageType::PointType  inputOrigin;
  typename OutImageType::PointType  outputOrigin;
  inputOrigin = r->GetOutput()->GetOrigin();

  typename InImageType::DirectionType  inputDirection;
  typename OutImageType::DirectionType  outputDirection;
  inputDirection = r->GetOutput()->GetDirection();

  for (unsigned int i = 0; i < OutImageType::ImageDimension; i++) {
    outSize[i] = r->GetOutput()->GetLargestPossibleRegion().GetSize()[i];
    outIndex[i] = r->GetOutput()->GetLargestPossibleRegion().GetIndex()[i];
    outputOrigin[i] = inputOrigin[i];
    outSpacing[i] = inSpacing[i];
  }

  for (unsigned int r = 0; r < OutImageType::ImageDimension; r++) {
    for (unsigned int c = 0; c < OutImageType::ImageDimension; c++) {
      outputDirection(r,c) = inputDirection(r,c);
    }
  }

  typename OutImageType::RegionType outRegion;
  outRegion.SetIndex( outIndex );
  outRegion.SetSize( outSize );

  typename OutImageType::Pointer outImage = OutImageType::New();
  // Adjust for the change in size
  outImage->SetRegions( outRegion );
  outImage->Allocate();

  // Copy geometry information from the input without change
  outImage->SetSpacing( outSpacing );
  outImage->SetDirection( outputDirection );

  // This program is designed for taking an ND array of probability values
  // and constructing an (N-1)D array with the index of the maximum 
  // probability at each voxel. As such, the Nth axis does not contain
  // useful geometry such as spacing, or direction cosines. Here we use the
  // (N-1)D geometry information to set the values for the output.
  outImage->SetOrigin(outputOrigin);
  
  // Here in the code, we have the input image read from disk and the output
  // image of the desired size allocated.
  OutRegionIteratorType outItr(outImage, outRegion);

  typename InImageType::IndexType inIndex;

  for (outItr.Begin(); !outItr.IsAtEnd(); ++outItr) {
    outIndex = outItr.GetIndex();
    for (unsigned int i = 0; i < outIndex.GetIndexDimension(); i++) {
      inIndex[i] = outIndex[i];
    }
    PixelType max = itk::NumericTraits<PixelType>::min();
    unsigned int maxindex = 0;
    for (unsigned int i = 0; i < inSize[InImageType::ImageDimension-1]; i++) {
      inIndex[InImageType::ImageDimension-1] = i;
      if (r->GetOutput()->GetPixel(inIndex) > max) {
        max = r->GetOutput()->GetPixel(inIndex);
        maxindex = i;
      }
    }
    outImage->SetPixel(outIndex, maxindex);
  }

  w->SetInput( outImage );

  try {
    w->Update();
    w->UseCompressionOn();
    w->Write();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

