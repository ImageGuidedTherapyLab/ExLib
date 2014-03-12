#include <ImageOperations.h>

#include <itkMinimumMaximumImageCalculator.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkStatisticsImageFilter.h>

void ImageOperations::reportImageStatistics( ImageTypeDefinitions::ImageType::Pointer in) {

	itk::StatisticsImageFilter< ImageTypeDefinitions::ImageType >::Pointer calc =
		itk::StatisticsImageFilter< ImageTypeDefinitions::ImageType >::New();
	calc->SetInput( in );
	calc->Update();
	std::cout << "ITK image range is (" << calc->GetMinimum() << ", "
		<< calc->GetMaximum() << ")" << std::endl;
	std::cout << "ITK image mean and stddev is (" << calc->GetMean() << ", "
		<< calc->GetSigma() << ")" << std::endl;
};

void ImageOperations::calculateImageStatistics( 
	ImageTypeDefinitions::ImageType::Pointer in, 
	ImageTypeDefinitions::InternalPixelType & min,
	ImageTypeDefinitions::InternalPixelType & max,
	ImageTypeDefinitions::InternalPixelType & mean,
	ImageTypeDefinitions::InternalPixelType & stddev
	) {
		itk::StatisticsImageFilter< ImageTypeDefinitions::ImageType >::Pointer 
			calc = itk::StatisticsImageFilter< ImageTypeDefinitions::ImageType >::New();
		calc->SetInput( in );
		calc->Update();
		min = calc->GetMinimum();
		max = calc->GetMaximum();
		mean = calc->GetMean();
		stddev = calc->GetSigma();
};

void ImageOperations::calculateOtsuThreshold( ImageTypeDefinitions::ImageType::Pointer in, ImageTypeDefinitions::InternalPixelType threshold ) {
	typedef itk::OtsuThresholdImageFilter<
		ImageTypeDefinitions::ImageType, ImageTypeDefinitions::ImageType >
		OtsuFilter;
	OtsuFilter::Pointer otsuFilter = OtsuFilter::New();
	otsuFilter->SetInput( in );
	otsuFilter->SetNumberOfHistogramBins( 128 );
	otsuFilter->Update();
	threshold = otsuFilter->GetThreshold();
};
