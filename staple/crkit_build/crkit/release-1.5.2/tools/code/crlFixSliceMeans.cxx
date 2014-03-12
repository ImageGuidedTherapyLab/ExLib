
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "vnl/vnl_math.h"

int main( int argc, char ** argv )
{
  if ( argc != 5 && argc != 6 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " " 
      << "inputImageFile axisDir1 axisDir2 [maskImageFile] ouptutImageFile"
      << std::endl
      << std::endl
      << "axisDir1 and axisDir2 are usually between 0 and 2 (for a 3D image)"
      << std::endl
      << "and define the axes WITHIN a slice."
      << std::endl
      << "For now, if a mask image is not offered, then flat 0.0 is the background "
      << std::endl
      << "and the image is assumed to have been previously masked."
      << std::endl
      << "Note that the output is always LAST" << std::endl;

    return 1;
    }

  const char* infile = argv[1];
  const unsigned int   axis0 = static_cast<unsigned int>(::atoi( argv[2]));
  const unsigned int   axis1 = static_cast<unsigned int>(::atoi( argv[3]));
  const unsigned int   axis2 = 3 - axis0 - axis1;        // yes, Dorothy, this works for 0..2
  const char* outfile = ( argc == 6 ) ? argv[5] : argv[4];
  const char* maskfile = ( argc == 6 ) ? argv[4] : 0;

  std::cerr << "infile=" << infile << std::endl;
  std::cerr << "axes = " << axis0 << "," << axis1 << std::endl;
  std::cerr << "slice axis = " << axis2 << std::endl;
  std::cerr << "outfile=" << outfile << std::endl;
  if ( maskfile ) std::cerr << "maskfile=" << maskfile << std::endl;


  typedef float PixelType;
  typedef itk::Image< PixelType, 3 >  ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  ImageType::RegionType::SizeType size;

  ReaderType::Pointer reader = ReaderType::New();
  ReaderType::Pointer maskreader = ReaderType::New();
  reader->SetFileName( infile );
  if ( maskfile) maskreader->SetFileName( maskfile );
  try
    {
    reader->Update();
    size = reader->GetOutput()->GetBufferedRegion().GetSize();
    if (maskfile) 
      {
      ImageType::RegionType::SizeType newsize;
      maskreader->Update();
      newsize = maskreader->GetOutput()->GetBufferedRegion().GetSize();
      for ( int i = 0; i < ImageType::ImageDimension; i++ )
        {
        if ( size[i] != newsize[i] )
          {
          std::cerr << "image size and mask size are not the same " << std::endl;
          return 1;
          }
        }
      }
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return 1;
    }

  ImageType::Pointer input = reader->GetOutput();
  ImageType::Pointer mask = (maskfile) ? maskreader->GetOutput() : 0;

  typedef itk::ImageSliceIteratorWithIndex<ImageType> SliceIterType;
  SliceIterType inputIter( input, input->GetBufferedRegion() );
  SliceIterType* maskIterp = 0;
  if ( mask ) maskIterp = new SliceIterType( mask, mask->GetBufferedRegion() );

  // Generate the slice means
  std::cerr << "generating slice means" << std::endl;
  std::vector<double> means(size[axis2],0.0) ;
  std::vector<long unsigned int> npixels(size[axis2],0U);

  inputIter.SetFirstDirection( axis0 );
  inputIter.SetSecondDirection( axis1 );
  if ( mask ) 
    {
    maskIterp->SetFirstDirection( axis0 );
    maskIterp->SetSecondDirection( axis1 );
    }

  inputIter.GoToBegin(); 
  if ( mask ) maskIterp->GoToBegin();

  unsigned int sliceno = 0;
  while( !inputIter.IsAtEnd() ) 
    { 
    std::cerr << "processing slice0 " << sliceno << std::endl;
    long unsigned int nslicepixels = 0;
    while ( !inputIter.IsAtEndOfSlice() ) 
      { 
      while ( !inputIter.IsAtEndOfLine() )
        {
        if ( ( mask && maskIterp->Value() != 0.0 ) || 
              ( !mask && inputIter.Value() != 0.0 ) ) 
          {
          ++nslicepixels;
          means[ sliceno ] += inputIter.Value();
          }
        ++inputIter; 
        if ( mask ) ++(*maskIterp);
        } 
      inputIter.NextLine(); 
      if ( mask ) maskIterp->NextLine();
      } 
    inputIter.NextSlice();
    if ( mask ) maskIterp->NextSlice();


    if ( nslicepixels > 0) means[ sliceno ] /= static_cast<double>( nslicepixels );
    npixels[sliceno] = nslicepixels;
    ++sliceno;
    }

  std::cerr << "slice means: ";
  for ( unsigned int i = 0; i < size[axis2]; i++ )
    {
    std::cerr << means[i] << " ";
    }
  std::cerr << std::endl;

  // find middle three slices and take mean
  long unsigned int s0 = size[axis2]/2;
  double normtarget=0.0;
  long unsigned int midpixels=0;
  for (long unsigned int i = vnl_math_max(s0-1, 0UL); i < vnl_math_min(s0+1,size[axis2]-1); i++)
    {
    normtarget += means[i]*static_cast<double>(npixels[i]);
    midpixels += npixels[i];
    }
  normtarget /= static_cast<double>(midpixels);

  std::cerr << "mean of middle three slices is " << normtarget << std::endl;


  // generate correction factors
  std::vector<double> corrections( size[axis2], 0.0 );
  for (unsigned int i = 0; i < size[axis2]; i++)
    {
    if ( means[i] > 0 )
      {
      corrections[i] = normtarget / means[i];
      }
    }
  for ( unsigned int i = 0; i < size[axis2]; i++ )
    {
    if ( corrections[i] == 0.0 )
      {
      // find closest neighbor and use it
      double corr = 1.0;        // this will be the default
      unsigned int dist = size[axis2];
      for ( unsigned int j = 0; j < i; j++ )
        {
        if ( corrections[j] > 0.0 )
          {
          corr = corrections[j];
          dist = i-j;
          }
        }
      for ( unsigned int j = i+1; j-i < dist && j < size[axis2]; j++ )
        {
        if ( corrections[j] > 0.0 )
          {
          corr = corrections[j];
          dist = j-i;
          }
        }
      corrections[i] = corr;
      }
    }

  std::cerr << "slice corrections: ";
  for (unsigned int i = 0; i < size[axis2]; i++)
    std::cerr << corrections[i] << ", ";
  std::cerr << std::endl;


  // apply correction factors to all pixels
  inputIter.SetFirstDirection( axis0 );
  inputIter.SetSecondDirection( axis1 );

  inputIter.GoToBegin(); 

  sliceno = 0;
  while( !inputIter.IsAtEnd() ) 
    { 
//    std::cerr << "processing slice0 " << sliceno << std::endl;
    while ( !inputIter.IsAtEndOfSlice() ) 
      { 
      while ( !inputIter.IsAtEndOfLine() )
        {
        inputIter.Value() *= corrections[sliceno];
        ++inputIter; 
        } 
      inputIter.NextLine(); 
      } 
    inputIter.NextSlice();

    ++sliceno;
    }


  // write the output

  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outfile );

  try
    {
    writer->SetInput(reader->GetOutput());
    writer->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return 1;
    }

  return 0;

}

