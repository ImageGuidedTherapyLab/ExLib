#include <iostream>
#include <itkOrientedImage.h>
#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <tclap/CmdLine.h>
#include <string>
#include <algorithm>
#include "configuration.h"

int main(int argc, char *argv[])
{
    std::string inImage;
    std::string outImage;
    std::string type;

    try
	{
	TCLAP::CmdLine cmd("Computational Radiology Laboratory", ' ',
		CRKIT_VERSION_STRING);

	TCLAP::ValueArg<std::string> typeArg("t", "type",
		"type of norm to compute, either: L2|Linf", false, "L2",
		"L2|Linf", cmd);
	TCLAP::UnlabeledValueArg<std::string> inputArg("input",
		"Input image filename", true, "", "input image", cmd);
	TCLAP::UnlabeledValueArg<std::string> outputArg("output",
		"Output image filename", true, "", "output image", cmd);

	cmd.parse(argc, argv);

	inImage = inputArg.getValue();
	outImage = outputArg.getValue();
	type = typeArg.getValue();
	std::transform(type.begin(), type.end(), type.begin(), ::toupper); // convert to upper case just to make typing easier

	if (type != "L2" && type != "LINF")
	    throw TCLAP::ArgException("type must be either L2 or Linf");
	}
    catch (TCLAP::ArgException& e)
	{
	std::cout << "argument error: " << e.error() << std::endl;
	return EXIT_FAILURE;
	}

    std::cout << "compute " << type << " norm of " << inImage << " into "
	    << outImage << std::endl;

    const unsigned int Dimension = 3;
    typedef itk::VectorImage<float, Dimension> VectorImageType;
    typedef itk::OrientedImage<float, Dimension> ScalarImageType;

    try
	{
	// read image
	typedef itk::ImageFileReader<VectorImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(inImage);
	reader->Update();
	VectorImageType::Pointer vimage = reader->GetOutput();
	vimage->DisconnectPipeline();

	// allocate scalar image
	ScalarImageType::Pointer norm = ScalarImageType::New();
	norm->SetRegions(vimage->GetBufferedRegion());
	norm->CopyInformation(vimage);
	norm->Allocate();

	// norm image
	typedef itk::ImageRegionConstIterator<VectorImageType>
		VectorConstIterator;
	typedef itk::ImageRegionIterator<ScalarImageType> ScalarIterator;

	VectorConstIterator viter(vimage, vimage->GetBufferedRegion());
	ScalarIterator niter(norm, norm->GetBufferedRegion());
	viter.GoToBegin();
	niter.GoToBegin();
	VectorImageType::PixelType vtemp;
	const unsigned int n = vimage->GetVectorLength();
	if (type == "L2")
	    {
	    // L2 norm
	    std::cout << "L2 norm..." << std::endl;
	    while (!niter.IsAtEnd())
		{
		vtemp = viter.Get();
		double sum = 0.0;
		for (unsigned int i = 0; i < n; ++i)
		    sum += vtemp[i] * vtemp[i];
		niter.Set(static_cast<float> (vcl_sqrt(sum)));
		++viter, ++niter;
		}
	    }
	else
	    {
	    // Linf norm
	    std::cout << "Linf norm..." << std::endl;
	    while (!niter.IsAtEnd())
		{
		vtemp = viter.Get();
		double max = -1;
		for ( unsigned int i = 0; i < n; ++i )
		    {
		    double temp = vcl_abs(vtemp[i]);
		    if ( temp > max ) max = temp;
		    }
		niter.Set(max);
		++viter, ++niter;
		}

	    }

	// write scalar image
	typedef itk::ImageFileWriter< ScalarImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(outImage);
	writer->SetInput(norm);
	writer->Update();
	}
    catch (itk::ExceptionObject& e)
	{
	std::cerr << "ITK error: " << e << std::endl;
	return EXIT_FAILURE;
	}

    return EXIT_SUCCESS;
}
