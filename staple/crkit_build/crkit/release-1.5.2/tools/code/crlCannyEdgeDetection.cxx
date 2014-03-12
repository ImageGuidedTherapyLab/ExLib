#include <iostream>
#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkCannyEdgeDetectionImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <tclap/CmdLine.h>
#include "configuration.h"


int main(int argc, char *argv[])
{

  TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ', CRKIT_VERSION_STRING );
  TCLAP::ValueArg<int> sliceArg("a", "axis",
      "Input axis to slice along for multi-slice 2D processing", false, -1,
      "slice axis", cmd);
  TCLAP::UnlabeledValueArg<std::string> imageArg("inimage",
      "Input filename of image", true, "", "input image", cmd);
  TCLAP::UnlabeledValueArg<std::string> outputArg("outimage",
      "Ouput edge image filename", true, "", "output image", cmd);
  TCLAP::UnlabeledMultiArg<float> multiArg("optional", "[variance maxerror [lthresh uthresh]]", false, "[variance maxerror [lthresh uthresh]]", cmd);


  std::vector<float> optArgs;
  try
    {
    cmd.parse(argc,argv);

    optArgs = multiArg.getValue();
    std::cout << "optArgs size is " << optArgs.size() << std::endl;

    if ( optArgs.size() == 1 )
      {
      throw TCLAP::CmdLineParseException("if the variance is set, you must also set the max error");
      }
    if ( optArgs.size() == 3 )
      {
      throw TCLAP::CmdLineParseException("if the lower threshold is set, you must also set the upper threshold");
      }

    }
  catch (TCLAP::ArgException& e)
    {
    std::cout << "argument error: " << e.error() << std::endl;
    return 1;
    }

  const int Dimension = 3;
  int axis = sliceArg.getValue();

  if ( axis >= Dimension )
    throw std::runtime_error( "axis specified to slice must be between [0,2]");

  typedef itk::Image<float,Dimension> ImageType;
  typedef itk::Image<float,Dimension-1> SliceType;
  typedef itk::ImageRegionConstIteratorWithIndex< SliceType > ConstSliceIterator;
  typedef itk::ImageRegionIteratorWithIndex<ImageType > ImageIterator;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName(imageArg.getValue());
  writer->SetFileName(outputArg.getValue());

  try
    {
      reader->Update();
    }
  catch (...)
    {
      std::cerr << "Failed to read input image: " << argv[1] << std::endl;
      return EXIT_FAILURE;
    }

  typedef itk::CannyEdgeDetectionImageFilter<ImageType,ImageType> EdgeFilter3DType;
  typedef itk::CannyEdgeDetectionImageFilter<SliceType,SliceType> EdgeFilter2DType;

  ImageType::Pointer output;

  if ( axis < 0 )
    {
    // running in natural 3D mode.
    EdgeFilter3DType::Pointer filter = EdgeFilter3DType::New();
    filter->SetInput( reader->GetOutput() );
    filter->SetVariance(5.0);
    if ( optArgs.size() > 0 )
      {
      filter->SetVariance(optArgs[0]);
      filter->SetMaximumError(optArgs[1]);
      }
    if (optArgs.size() > 2 )
      {
      filter->SetLowerThreshold(optArgs[2]);
      filter->SetUpperThreshold(optArgs[3]);
      }
    filter->Update();
    output = filter->GetOutput();
    }
  else
    {
    // running slice by slice, first set up filter
    EdgeFilter2DType::Pointer filter = EdgeFilter2DType::New();
    filter->SetVariance(5.0);
    if ( optArgs.size() > 0 )
      {
      filter->SetVariance(optArgs[0]);
      filter->SetMaximumError(optArgs[1]);
      }
    if (optArgs.size() > 2 )
      {
      filter->SetLowerThreshold(optArgs[2]);
      filter->SetUpperThreshold(optArgs[3]);
      }

    // allocate output image
    output = ImageType::New();
    output->SetRegions( reader->GetOutput()->GetBufferedRegion() );
    output->CopyInformation( reader->GetOutput() );
    output->Allocate();
    ImageType::RegionType region3D = reader->GetOutput()->GetBufferedRegion();

    typedef itk::ExtractImageFilter< ImageType, SliceType > ExtractType;
    ExtractType::Pointer extractor = ExtractType::New();

    unsigned int asize = region3D.GetSize(axis);
    unsigned int aindex = region3D.GetIndex(axis);

    extractor->SetInput(reader->GetOutput());
    for ( unsigned int i = 0; i < asize; i++ )
      {
      unsigned int slice = i + aindex;
      std::cout << "slice " << slice << std::endl;
      region3D.SetIndex(axis, slice);
      region3D.SetSize(axis, 0);		// zero is a flag for an axis to compress out in ExtractImageFilter
      extractor->SetExtractionRegion(region3D);
      filter->SetInput( extractor->GetOutput() );
      filter->Update();
      // okay, now copy the slice over to the output
      ConstSliceIterator srcIter( filter->GetOutput(), filter->GetOutput()->GetBufferedRegion() );
      region3D.SetSize(axis, 1);		// reset the flag value 0 back to 1 for a valid slice region
      ImageIterator dstIter( output, region3D );

      for ( srcIter.GoToBegin(), dstIter.GoToBegin(); !dstIter.IsAtEnd(); ++srcIter, ++dstIter )
	{
	dstIter.Set( srcIter.Value() );
	}
      if ( srcIter.IsAtEnd() != dstIter.IsAtEnd() )
	{
	throw std::runtime_error("BUG: one iterator ended, but not the other???");
	}
      }

    }

  writer->SetInput(output);

  try
    {
      writer->Write();;
    }
  catch (...)
    {
      std::cerr << "Failed to write image: " << argv[2] << std::endl;
      return EXIT_FAILURE;
    }


  return EXIT_SUCCESS;
}
