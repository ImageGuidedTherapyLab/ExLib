
// A class to provide access to simple calculations and 
// transformations of images.

#ifndef _IMAGE_OPERATIONS_INCLUDED
#define _IMAGE_OPERATIONS_INCLUDED 1

#include <ImageTypeDefinitions.h>

class ImageOperations
{
public:
	static void reportImageStatistics( ImageTypeDefinitions::ImageType::Pointer in);

	static void calculateImageStatistics( 
		ImageTypeDefinitions::ImageType::Pointer in, 
		ImageTypeDefinitions::InternalPixelType & min,
		ImageTypeDefinitions::InternalPixelType & max,
		ImageTypeDefinitions::InternalPixelType & mean,
		ImageTypeDefinitions::InternalPixelType & stddev
		);

	static void calculateOtsuThreshold( ImageTypeDefinitions::ImageType::Pointer in, ImageTypeDefinitions::InternalPixelType threshold );
};

#endif

