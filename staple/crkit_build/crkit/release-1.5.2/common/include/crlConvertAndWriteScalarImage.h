/**********************************************************************************************//**
 * \file	crlCommon\crlConvertAndWriteImage.h
 *
 * \brief	Declares functions to write scalar images. 
*************************************************************************************************/

#include <string>

#include <itkCastImageFilter.h>
#include <itkImageFileWriter.h>

namespace crl {

/**********************************************************************************************//**
 * \brief	Converts an input scalar image of voxel type TInputVoxel and dimension VDim to an
 * 			image of voxel type TOutputVoxel and dimension VDim, and writes it out. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \param	outputFileName	Filename of the output file. 
 * \param	image			The image. 
 * 
*************************************************************************************************/
template <typename TInputVoxel, typename TOutputVoxel, unsigned int VDim>
void ConvertAndWriteScalarImage( const std::string& outputFileName, const  itk::Image<TInputVoxel,VDim> *image )
{
	/*-------------------------------------
	 First convert image component type
	-------------------------------------*/
	typedef itk::CastImageFilter< itk::Image<TInputVoxel,VDim>, itk::Image<TOutputVoxel,VDim> > ImageCaster;
	typename ImageCaster::Pointer caster = ImageCaster::New();
	caster->SetInput(image);

	/*-------------------------------------
	 And write the image
	-------------------------------------*/
	typedef itk::ImageFileWriter< itk::Image<TOutputVoxel,VDim> > WriterType;
	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( outputFileName  );
	writer->SetUseCompression(true);
	writer->SetInput(  caster->GetOutput() );
	writer->Update();
}

/**********************************************************************************************//**
 * \brief	Writes out a scalar image of voxel type TVoxel and dimension VDim with the new voxel
 * 			type described by outputType. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \exception	Thrown when error. 
 *
 * \param	outputFileName	Filename of the output file. 
 * \param [in,out]	image	If non-null, the image. 
 * \param	outputType		Type of the output. 
*************************************************************************************************/
template <typename TVoxel, unsigned int VDim>
void WriteScalarImage( const std::string& outputFileName,  itk::Image<TVoxel,VDim> *image, itk::ImageIOBase::IOComponentType outputType )
{
		switch(outputType)
		{
		case itk::ImageIOBase::UCHAR:
			ConvertAndWriteScalarImage<TVoxel, unsigned char, VDim> (outputFileName, image); break;
		case itk::ImageIOBase::CHAR:
			ConvertAndWriteScalarImage<TVoxel, char, VDim> (outputFileName, image); break;
		case itk::ImageIOBase::USHORT:
			ConvertAndWriteScalarImage<TVoxel, unsigned short, VDim> (outputFileName, image); break;
		case itk::ImageIOBase::SHORT:
			ConvertAndWriteScalarImage<TVoxel, short, VDim> (outputFileName, image); break;
		case itk::ImageIOBase::UINT:
			ConvertAndWriteScalarImage<TVoxel, unsigned int, VDim> (outputFileName, image); break;
		case itk::ImageIOBase::INT:
			ConvertAndWriteScalarImage<TVoxel, int, VDim> (outputFileName, image); break;
		case itk::ImageIOBase::ULONG:
			ConvertAndWriteScalarImage<TVoxel, unsigned long, VDim> (outputFileName, image); break;
		case itk::ImageIOBase::LONG:
			ConvertAndWriteScalarImage<TVoxel, long, VDim> (outputFileName, image); break;
		case itk::ImageIOBase::FLOAT:
			ConvertAndWriteScalarImage<TVoxel, float, VDim> (outputFileName, image); break;
		case itk::ImageIOBase::DOUBLE:  
			ConvertAndWriteScalarImage<TVoxel, double, VDim> (outputFileName, image); break;
		case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
		default:
			throw new itk::ExceptionObject("file");
			break;
		}
}

}
