/**********************************************************************************************//**
 * \file	crlGeneralizedInstantiator.h
 *
 * \brief	Declare the crl::Instantiator class. 
*************************************************************************************************/

/*
 * Copyright (c) 2008-2009 Children's Hospital Boston.
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 * http://www.crl.med.harvard.edu
 * For more information contact: simon.warfield@childrens.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/


#ifndef CRL_GENERALIZED_INSTANTIATOR_H
#define CRL_GENERALIZED_INSTANTIATOR_H

// Remove annoying warning with Visual Studio
#ifdef _MSC_VER
  #pragma warning( disable: 4996 )
#endif

#include <iostream>
#include <stdlib.h>

#include <itkImage.h>
#include <itkImageIOBase.h>
#include <itkImageFileReader.h>

//------------------------------------------------------
// Check that at least one of the macro is set
//------------------------------------------------------
#ifndef CRL_GA_SCALAR_CLASS
#ifndef CRL_GA_VECTORIMAGE_CLASS
#ifndef CRL_GA_DIFFTENSOR_CLASS  
#error At least one of these macros should be set: CRL_GA_SCALAR_CLASS, CRL_GA_VECTOR_CLASS, or CRL_GA_DIFFTENSOR_CLASS
#endif
#endif
#endif

//------------------------------------------------------
// If the following constants were not defined 
// before the inclusion of crlGeneralizedInstantiator.h,
// use default settings
//------------------------------------------------------
#ifndef MAX_IMAGE_DIM
  #define MAX_IMAGE_DIM 4
#endif
#ifndef MIN_IMAGE_DIM
  #define MIN_IMAGE_DIM 1
#endif

#ifndef MAX_VEC_COMPONENT_DIM 
  #define MAX_VEC_COMPONENT_DIM 6
#endif


namespace crl {

	/**********************************************************************************************//**
	 * \class	Instantiator
	 *
	 * \brief	Base class for the CRL Generalized Instantiator mecanism.
	 * 			
	 * 			The CRL Generalized Instantiator allows to create ITK tools that keep the input files
	 * 			type (pixel type, component type, number of dimensions, etc). It provides mecanisms
	 * 			to automatically instantiate an ITK template based on a Instantiator::ImageTypeInfos
	 * 			structure. 
	 *
	 *			It is based on a non-templatized algorithm base class (AlgoBaseType) overloaded by
	 *			(multiple) template versions of it.
	 *
	 * Here is some part of code to give an idea how it works. For further information,
	 * look at the source code of crlMaskImages or crlDWIConvertNHDRForFSL.
	 *
	 *<b> Convert_VectorImage_To_ND.h </b>
* \code
* #include "crlGeneralizedInstantiatorMacros.h"
*
* //-----------------------------------------------------------------------
* // DECLARE THE ALGORITHM BASE CLASS
* //-----------------------------------------------------------------------
* class ITK_EXPORT Convert_VectorImage_To_ND_Base : public itk::LightObject
* {
*   CRL_INSTANCIATOR_BASE(Convert_VectorImage_To_ND_Base);
*
* public:
* 	// Algo parameters
* 	void SetInputFileName (const std::string& fileName ) { m_InputFileName=fileName; }
* 	void SetOutputFileName (const std::string& fileName ) { m_OutputFileName=fileName; }
* 
* 	// Function to implement 
* 	virtual int Execute() = 0;
* 
* protected:
* 	std::string		m_InputFileName;
* 	std::string		m_OutputFileName;
* 
* (...)
* };
*
* //-----------------------------------------------------------------------
* // DECLARE A TEMPLATE VERSION OF THE ALGORITHM
* //-----------------------------------------------------------------------
* template <typename TPixel, unsigned int VImageDimension  >
* class ITK_EXPORT Convert_VectorImage_To_ND : public Convert_VectorImage_To_ND_Base
* {
*   CRL_INSTANCIATOR_CLASS(Convert_VectorImage_To_ND_Base, Convert_VectorImage_To_ND);
*
* public:
*   // Implement the pure-virtual function for that template
*	virtual int Execute()
*	{
*		ConvertVectorImageToNDImage<TPixel, VImageDimension>( m_InputFileName, m_OutputFileName  );
*		return 0;
*	}
*   
*   (...) 
* };
* \endcode
* 
*
* Then, to instantiate this class:
* <b> main.cxx </b>
* \code
* //-----------------------------------------------------------------------
* // Needs some declarations before including crlGeneralizedInstantiator.h
* //-----------------------------------------------------------------------
*
* // I. Include our template files (here only one)
* #include "Convert_VectorImage_To_ND.h"
*
* // II. Define some constants for crlGeneralizedInstantiator.h
* #undef       CRL_GA_SCALAR_CLASS			
* #define      CRL_GA_VECTORIMAGE_CLASS Convert_VectorImage_To_ND
* #undef       CRL_GA_VECTOR_CLASS
* #undef       CRL_GA_DIFFTENSOR_CLASS		
*
* typedef		crl::Convert_VectorImage_To_ND_Base AlgoBaseType;
*
* //-----------------------------------------------------------------------
* // Now includes the file crlGeneralizedInstantiator.h
* //-----------------------------------------------------------------------
* #include	"crlGeneralizedInstantiator.h"
*
*
* int main(int argc, char **argv)
* {
*    (...)
*
*    // Instantiate the algorithm for a given file (fileInput)
*    crl::Instantiator Instantiator;
*    AlgoBaseType::Pointer convertAlgo = Instantiator.Instantiate( fileInput );
*    
*    // Set the parameters for the algo
*    convertAlgo->SetInputFileName(fileInput);
*    convertAlgo->SetOutputFileName(fnData);
*    
*    // And execute!
*    convertAlgo->Execute();
*
*    (...)
* }
*
*\endcode
*
	 * \author	Benoit Scherrer
	 * \date	December 2009
	*************************************************************************************************/
	class Instantiator
	{
	public:

		/**********************************************************************************************//**
		 * \class	ImageTypeInfos
		 *
		 * \brief	Subclass containing the image informations. 
		 *
		 * \author	Benoit Scherrer
		 * \date	December 2009
		*************************************************************************************************/
		class ImageTypeInfos
		{
		public:
			std::string							filename;
			itk::ImageIOBase::IOPixelType		pixelType;
			itk::ImageIOBase::IOComponentType	componentType;
			unsigned int						numberOfComponents;
			itk::ImageIOBase::Pointer			imageIO;
			unsigned int						dimensionality ;
		};

		/**********************************************************************************************//**
		 * \fn	AlgoBaseType::Pointer Instantiate(const std::string& fileName)
		 *
		 * \brief	Creates an instance of the algorithm AlgoBaseType for a given file. 
		 *
		 * \author	Benoit Scherrer
		 * \date	December 2009
		 *
		 * \param	fileName	Filename of the file. 
		 *
		 * \return	. 
		*************************************************************************************************/
		AlgoBaseType::Pointer Instantiate(const std::string& fileName)
		{
			AlgoBaseType::Pointer algoPointer = 0;

			ImageTypeInfos infos = GetImageType( fileName );
			return Instantiate(infos);
		}

		/**********************************************************************************************//**
		 * \fn	AlgoBaseType::Pointer Instantiate(itk::ImageIOBase::IOPixelType pixelType,
		 * 		itk::ImageIOBase::IOComponentType componentType, unsigned int dimensionality)
		 *
		 * \brief	Creates an instance of the algorithm AlgoBaseType by forcing the pixel type,
		 * 			component type and dimension. 
		 *
		 * \author	Benoit Scherrer
		 * \date	November 2010
		 *
		 * \param	pixelType		Type of the pixel. 
		 * \param	componentType	Type of the component. 
		 * \param	dimensionality	The dimensionality. 
		 *
		 * \return	. 
		*************************************************************************************************/
		AlgoBaseType::Pointer Instantiate(itk::ImageIOBase::IOPixelType pixelType, itk::ImageIOBase::IOComponentType componentType, unsigned int dimensionality)
		{
			ImageTypeInfos infos ;
			infos.pixelType = pixelType;
			infos.componentType = componentType;
			infos.dimensionality = dimensionality;
			return Instantiate(infos);

		}

		/**********************************************************************************************//**
		 * \fn	AlgoBaseType::Pointer Instantiate (const ImageTypeInfos &infos)
		 *
		 * \brief	Creates an instance of the algorithm AlgoBaseType for a given Instantiator::
		 * 			ImageTypeInfos structure.
		 * 			
		 * 			This function can be overloaded to introduce specific instantiation mecanisms (see
		 * 			crlMaskImages). 
		 *
		 * \author	Benoit Scherrer
		 * \date	December 2009
		 *
		 * \param	infos	The infos. 
		 *
		 * \return	. 
		*************************************************************************************************/
		AlgoBaseType::Pointer Instantiate (const ImageTypeInfos &infos)
		{
#ifdef CRL_GA_SCALAR_CLASS
			if (infos.pixelType == itk::ImageIOBase::SCALAR)
			{
				switch(infos.componentType) 
				{
				case itk::ImageIOBase::UCHAR : 
					return InstantiateScalar<unsigned char, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::CHAR : 
					return InstantiateScalar<char, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::USHORT : 
					return InstantiateScalar<unsigned short, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::SHORT : 
					return InstantiateScalar<short, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::UINT :  
					return InstantiateScalar<unsigned int, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::INT :   
					return InstantiateScalar<int, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::ULONG : 
					return InstantiateScalar<unsigned long, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::LONG :  
					return InstantiateScalar<long, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::FLOAT : 
					return InstantiateScalar<float, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::DOUBLE : 
					return InstantiateScalar<double, MAX_IMAGE_DIM>::Instantiate (infos);
					break;	
				case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE : 
				default: 
					return NULL;
					break; 
				} 

			}
#endif


#ifdef CRL_GA_VECTORIMAGE_CLASS
			if (infos.pixelType == itk::ImageIOBase::VECTOR)
			{
				switch(infos.componentType) 
				{
				case itk::ImageIOBase::UCHAR : 
					return InstantiateVectorImage<unsigned char, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::CHAR : 
					return InstantiateVectorImage<char, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::USHORT : 
					return InstantiateVectorImage<unsigned short, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::SHORT : 
					return InstantiateVectorImage<short, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::UINT :  
					return InstantiateVectorImage<unsigned int, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::INT :   
					return InstantiateVectorImage<int, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::ULONG : 
					return InstantiateVectorImage<unsigned long, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::LONG :  
					return InstantiateVectorImage<long, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::FLOAT : 
					return InstantiateVectorImage<float, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::DOUBLE : 
					return InstantiateVectorImage<double, MAX_IMAGE_DIM>::Instantiate (infos);
					break;	
				case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE : 
				default: 
					return NULL;
					break; 
				} 

			}
#endif


#ifdef CRL_GA_VECTOR_CLASS
			if (infos.pixelType == itk::ImageIOBase::VECTOR) 
			{
				switch(infos.componentType) 
				{
				case itk::ImageIOBase::UCHAR : 
					return InstantiateVector<unsigned char, MAX_VEC_COMPONENT_DIM, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::CHAR : 
					return InstantiateVector<char, MAX_VEC_COMPONENT_DIM, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::USHORT : 
					return InstantiateVector<unsigned short, MAX_VEC_COMPONENT_DIM, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::SHORT : 
					return InstantiateVector<short, MAX_VEC_COMPONENT_DIM, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::UINT :  
					return InstantiateVector<unsigned int, MAX_VEC_COMPONENT_DIM, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::INT :   
					return InstantiateVector<int, MAX_VEC_COMPONENT_DIM, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::ULONG : 
					return InstantiateVector<unsigned long, MAX_VEC_COMPONENT_DIM, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::LONG :  
					return InstantiateVector<long, MAX_VEC_COMPONENT_DIM, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::FLOAT : 
					return InstantiateVector<float, MAX_VEC_COMPONENT_DIM, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::DOUBLE : 
					return InstantiateVector<double, MAX_VEC_COMPONENT_DIM, MAX_IMAGE_DIM>::Instantiate (infos);
					break;	
				case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE : 
				default: 
					break; 
				} 
			} 
#endif
#ifdef CRL_GA_DIFFTENSOR_CLASS
			if ( 	(infos.imageIO->GetPixelType() == itk::ImageIOBase::SYMMETRICSECONDRANKTENSOR)
				|| (infos.imageIO->GetPixelType() == itk::ImageIOBase::DIFFUSIONTENSOR3D) 	) 
			{
				switch(infos.componentType) 
				{
				case itk::ImageIOBase::FLOAT : 
					return InstantiateTensor<float, MAX_IMAGE_DIM>::Instantiate (infos);
					break;
				case itk::ImageIOBase::DOUBLE : 
					return InstantiateTensor<double, MAX_IMAGE_DIM>::Instantiate (infos);
					break;	
				default: 
					break; 
				}

			}
#endif

			return NULL;
		}


	


#ifdef CRL_GA_SCALAR_CLASS

		/**********************************************************************************************//**
		 * \class	InstantiateScalar
		 *
		 * \brief	Recursive template used to instantiate a scalar version of the algo. 
		 *
		 * \author	Benoit Scherrer.
		 * \date	December 2009.
		*************************************************************************************************/

		template <typename TScalar, unsigned int VImageDim >  class InstantiateScalar 
		{
		public:
			static inline AlgoBaseType::Pointer Instantiate( const ImageTypeInfos &imageInfos )
			{
				if ( imageInfos.dimensionality == VImageDim )
					return (CRL_GA_SCALAR_CLASS<TScalar,VImageDim>::New().GetPointer());
				else
					return InstantiateScalar<TScalar, VImageDim-1>::Instantiate(imageInfos);


			}
		};

		/**********************************************************************************************//**
		 * \class	0>
		 *
		 * \brief	Template partial specialization designed to stop the recursive template generation. 
		 *
		 * \author	Benoit Scherrer
		 * \date	December 2009
		*************************************************************************************************/

		template <typename TScalar> class InstantiateScalar<TScalar, (MIN_IMAGE_DIM-1) >
		{
		public:
			inline static AlgoBaseType::Pointer Instantiate(const ImageTypeInfos &imageInfos) { return NULL; }
		};

#endif


#ifdef CRL_GA_VECTORIMAGE_CLASS

		/**********************************************************************************************//**
		 * \class	InstantiateVectorImage
		 *
		 * \brief	Recursive template used to instantiate a scalar version of the algo. 
		 *
		 * \author	Benoit Scherrer.
		 * \date	December 2009.
		*************************************************************************************************/

		template <typename TScalar, unsigned int VImageDim >  class InstantiateVectorImage 
		{
		public:
			static inline AlgoBaseType::Pointer Instantiate( const ImageTypeInfos &imageInfos )
			{
				if ( imageInfos.dimensionality == VImageDim )
					return (CRL_GA_VECTORIMAGE_CLASS<TScalar,VImageDim>::New().GetPointer());
				else
					return InstantiateVectorImage<TScalar, VImageDim-1>::Instantiate(imageInfos);
			}
		};

		/**********************************************************************************************//**
		 * \class	0>
		 *
		 * \brief	Template partial specialization designed to stop the recursive template generation. 
		 *
		 * \author	Benoit Scherrer
		 * \date	December 2009
		*************************************************************************************************/

		template <typename TScalar> class InstantiateVectorImage<TScalar, (MIN_IMAGE_DIM-1) >
		{
		public:
			inline static AlgoBaseType::Pointer Instantiate(const ImageTypeInfos &imageInfos) { return NULL; }
		};
#endif


#ifdef CRL_GA_VECTOR_CLASS

		/**********************************************************************************************//**
		 * \class	InstantiateVector
		 *
		 * \brief	Recursive template used to instantiate a vector version of the algo. 
		 *
		 * \author	Benoit Scherrer.
		 * \date	December 2009.
		*************************************************************************************************/

		template <typename TScalar, unsigned int VComponentDim, unsigned int VImageDim >  class InstantiateVector 
		{
		public:
			static inline  void Instantiate( const ImageTypeInfos &imageInfos )
			{
				if ( imageInfos.dimensionality == VImageDim && imageInfos.numberOfComponents == VComponentDim)
					CRL_GA_VECTOR_CLASS< itk::Vector<TScalar, VComponentDim> ,VImageDim>::New().GetPointer();
				InstantiateVector<TScalar, VComponentDim, VImageDim-1>::Instantiate(imageInfos);
				InstantiateVector<TScalar, VComponentDim-1, VImageDim>::Instantiate(imageInfos);
			}
		};

		/**********************************************************************************************//**
		 * \class	>
		 *
		 * \brief	Template partial specialization designed to stop the recursive template generation. 
		 *
		 * \author	Benoit Scherrer
		 * \date	December 2009
		*************************************************************************************************/

		template <typename TScalar, unsigned int VComponentDim> class InstantiateVector<TScalar, VComponentDim, (MIN_IMAGE_DIM-1) >
		{
		public:
			inline static void Instantiate(const ImageTypeInfos &imageInfos) {}
		};
		template <typename TScalar, unsigned int VImageDim> class InstantiateVector<TScalar, 0, VImageDim>
		{
		public:
			inline static void Instantiate(const ImageTypeInfos &imageInfos) {}
		};
		template <typename TScalar> class InstantiateVector<TScalar, 0, 0>
		{
		public:
			inline static void Instantiate(const ImageTypeInfos &imageInfos) {}
		};
#endif


#ifdef CRL_GA_DIFFTENSOR_CLASS
		/**********************************************************************************************//**
		 *
		 * \brief	Recursive template used to instanciate the tensor version of the algo. 
		 *
		 * \author	Benoit Scherrer.
		 * \date	December 2009.
		*************************************************************************************************/

		template <typename TScalar, unsigned int VImageDim >  class InstantiateTensor
		{
		public:
			static inline AlgoBaseType::Pointer Instantiate( const ImageTypeInfos &imageInfos )
			{
				if ( imageInfos.dimensionality == VImageDim )
					return (CRL_GA_DIFFTENSOR_CLASS<TScalar,VImageDim>::New().GetPointer());
				else
					return InstantiateTensor<TScalar, VImageDim-1>::Instantiate(imageInfos);
			}
		};

		/**********************************************************************************************//**
		 *
		 * \brief	Template partial specialization designed to stop the recursive template generation. 
		 *
		 * \author	Benoit Scherrer
		 * \date	December 2009
		*************************************************************************************************/

		template <typename TScalar> class InstantiateTensor<TScalar, (MIN_IMAGE_DIM-1) >
		{
		public:
			inline static AlgoBaseType::Pointer Instantiate(const ImageTypeInfos &imageInfos) { return NULL; }
		};
#endif

public:
		/**********************************************************************************************//**
		 * \fn	ImageTypeInfos GetImageType (const std::string &fileName )
		 *
		 * \brief	Open the header of an image file and provides a bunch of informations about it
		 *			in a Instantiator::ImageTypeInfos structure, including the pixel type, the
		 *			component type, etc).
		 *
		 * \author	Benoit Scherrer
		 * \date	December 2009
		 *
		 * \exception	Thrown	when error. 
		 *
		 * \param	fileName	Filename of the file. 
		 *
		 * \return	The image type. 
		 *
		 * \see Instantiator::ImageTypeInfos
		*************************************************************************************************/
		ImageTypeInfos GetImageType (const std::string &fileName )
		{
			ImageTypeInfos infos;


			// arbitrary image type, just to read the image infos
			typedef itk::Image<unsigned char, 2> ImageType;

			// Read the image informations
			itk::ImageFileReader<ImageType>::Pointer imageReader = itk::ImageFileReader<ImageType>::New();
			imageReader->SetFileName(fileName.c_str());
			try 
			{
				imageReader->UpdateOutputInformation();
			} 
			catch ( itk::ExceptionObject &err )
			{
				std::cout << "ExceptionObject caught!" << std::endl;
				std::cout << err << std::endl;
				throw err;
			}

			// 
			infos.imageIO = imageReader->GetImageIO();
			infos.pixelType = infos.imageIO->GetPixelType();
			infos.componentType = infos.imageIO->GetComponentType();
			infos.numberOfComponents = infos.imageIO->GetNumberOfComponents();
			infos.dimensionality = infos.imageIO->GetNumberOfDimensions();
			infos.filename = fileName;

			return infos;
		}


	};

}

#endif
