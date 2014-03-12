/*
 * Copyright 2007-2009 Children's Hospital Boston
 * Contact: Simon Warfield simon.warfield@childrens.harvard.edu
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
 */

#include "tclap/CmdLine.h"
#include "configuration.h"
#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkVectorImage.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkDanielssonDistanceMapImageFilter.h>


const unsigned int Dimension = 3;
typedef unsigned int SegPixelType;
typedef itk::OrientedImage<SegPixelType, Dimension> SegImageType;
typedef itk::ImageFileReader<SegImageType> ReaderType;
typedef itk::ImageRegionConstIterator<SegImageType> SegConstIter;
typedef itk::ImageRegionIterator<SegImageType> SegIterType;

typedef float AtlasPixelType;
typedef itk::VectorImage<AtlasPixelType, Dimension> AtlasImageType;
typedef itk::ImageFileWriter<AtlasImageType> WriterType;
typedef itk::ImageRegionIterator<AtlasImageType> AtlasIterType;

typedef itk::OrientedImage<float, Dimension> DistanceMapImageType;
typedef itk::ImageRegionConstIterator<DistanceMapImageType> MapIterType;

typedef itk::DanielssonDistanceMapImageFilter< SegImageType, DistanceMapImageType > DistanceMapType;

int main(int argc, char* argv[])
{
    float rho;
    std::string inFile;
    std::string outFile;

    try
	{
	TCLAP::CmdLine cmd("Computational Radiology Laboratory", ' ',
		CRKIT_VERSION_STRING);
	TCLAP::UnlabeledValueArg<std::string> inputFileArg("inputFile",
		"Input File Name", true, "", "input file name", cmd);
	TCLAP::ValueArg<std::string> outputFileArg("o","outputFile",
		"Output File Name", true, "", "output file name", cmd);
	TCLAP::ValueArg<float> rhoArg("r","rho", "smoothing value", false, 1.0, "smoothing value", cmd);

	cmd.parse(argc,argv);
	rho = rhoArg.getValue();
	inFile = inputFileArg.getValue();
	outFile = outputFileArg.getValue();
	}
    catch (TCLAP::ArgException& e)
	{
	std::cerr << "error: " << e.error() << " for argument " << e.argId()
		<< std::endl;
	exit(EXIT_FAILURE);
	}


    try
	{
	// read in the image
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( inFile );
	reader->Update();
	SegImageType::Pointer inImage = reader->GetOutput();
	inImage->DisconnectPipeline();

	// find the maximum label value
	unsigned int maxlabel = 0;
	SegConstIter iter( inImage, inImage->GetBufferedRegion());
	for ( iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
	    {
	    if (iter.Value() > maxlabel ) maxlabel = iter.Value();
	    }
	std::cout << "maximum image label is " << maxlabel << std::endl;
	std::vector<bool> foundlabels(maxlabel+1, false);
	for ( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
	    foundlabels[iter.Value()]=true;

	// allocate atlas image and temp seg image
    	AtlasImageType::Pointer atlas = AtlasImageType::New();
    	atlas->CopyInformation( inImage );
    	atlas->SetRegions( inImage->GetBufferedRegion() );
    	atlas->SetVectorLength( maxlabel+1 );
    	atlas->Allocate();
    	//ANONYMOUS
    	    {
    	    AtlasImageType::PixelType p(maxlabel+1);
    	    p.Fill(0.);
    	    atlas->FillBuffer(p);
    	    }
    	SegImageType::Pointer tempseg = SegImageType::New();
    	tempseg->CopyInformation(inImage);
    	tempseg->SetRegions(inImage->GetBufferedRegion());
    	tempseg->Allocate();


	// for each label value...

    	for ( unsigned int label = 0; label <= maxlabel; ++label )
    	    {
    	    if ( !foundlabels[label] )
    		continue;
    	    std::cout << "label " << label << std::endl;
	    // make a binary map of one particular label
    	    SegIterType dstIter( tempseg, tempseg->GetBufferedRegion() );
    	    SegConstIter srcIter( inImage, inImage->GetBufferedRegion() );
    	    tempseg->FillBuffer(0U);
    	    for ( dstIter.GoToBegin(), srcIter.GoToBegin(); !dstIter.IsAtEnd(); ++srcIter, ++dstIter)
    		{
    		if ( srcIter.Value() == label )
    		    {
    		    foundlabels[label]=true;
    		    dstIter.Set(1);
    		    }
    		}

    	    // distance map the temp image
    	    DistanceMapType::Pointer distmap = DistanceMapType::New();
    	    distmap->SetInput( tempseg );
    	    distmap->UseImageSpacingOn();
    	    distmap->Update();

    	    //DEBUGGING NIW
#if 0
    	    if ( foundlabels[label] )
    		{
    		std::ostringstream s;
    		s << std::string("distmap") << label << ".nrrd";
    		std::cout << "writing distance map to " << s.str() << std::endl;
    		typedef itk::ImageFileWriter<DistanceMapImageType> DistWriter;
    		DistWriter::Pointer writer = DistWriter::New();
    		writer->SetFileName(s.str());
    		writer->SetInput(distmap->GetDistanceMap());
    		writer->UseCompressionOn();
    		writer->Update();
    		}
#endif


	    // populate the atlas channel
    	    if ( foundlabels[label] )
    		{
    		AtlasIterType dstIter( atlas, atlas->GetBufferedRegion() );
    		MapIterType mapIter( distmap->GetDistanceMap(), distmap->GetDistanceMap()->GetBufferedRegion());
    		AtlasImageType::PixelType p;
    		for ( dstIter.GoToBegin(), mapIter.GoToBegin(); !dstIter.IsAtEnd(); ++dstIter, ++mapIter )
    		    {
    		    p = dstIter.Get();
    		    p[label] = vcl_exp(-1.0*rho*mapIter.Get() );
    		    dstIter.Set(p);
    		    }

    		}
    	    }

    	// normalize the atlas image
    	AtlasIterType dstIter( atlas, atlas->GetBufferedRegion() );
    	AtlasImageType::PixelType p;
    	for ( dstIter.GoToBegin(); !dstIter.IsAtEnd(); ++dstIter )
    	    {
    	    p = dstIter.Get();
    	    double norm = 0.0;
    	    for ( unsigned int i = 0; i < p.GetNumberOfElements(); ++i )
    		norm += p[i];
    	    if ( norm > 0 )
    		{
    		for ( unsigned int i = 0; i < p.GetNumberOfElements(); ++i )
    		    p[i] /= norm;
    		}
    	    else
    		{
    		p[0] = 1;
    		std::cerr << "warning: we found a zero voxel and set p(label=0)=1 for lack of anything better to do." << std::endl;
    		}
    	    dstIter.Set(p);
    	    }

	// write out the atlas image
    	std::cout << "writing output vector image atlas " << outFile << std::endl;
    	WriterType::Pointer writer = WriterType::New();
    	writer->SetFileName( outFile );
    	writer->SetInput( atlas );
    	writer->UseCompressionOn();
    	writer->Update();
	}
    catch ( itk::ExceptionObject& e )
	{
	std::cerr << "ITK Error: " << e << std::endl;
	return 1;
	}


    return 0;
}
