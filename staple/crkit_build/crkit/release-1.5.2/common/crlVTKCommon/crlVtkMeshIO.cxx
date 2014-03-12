#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataReader.h>
#include <vtkGIFTIReader.h>

#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkGIFTIWriter.h>


#include "crlVtkMeshIO.h"
#include "crlFileName.h"

/**********************************************************************************************//**
 * \fn	vtkPolyData *crlVtkMeshIO::ReadMesh(const std::string& fileName)
 *
 * \brief	Reads a mesh.
 *
 * \author	Benoit Scherrer
 * \date	July 2011
 *
 * \param	fileName	Filename of the file to read.
 *
 * \return	null if it fails, else the mesh.
 **************************************************************************************************/
vtkPolyData *crlVtkMeshIO::ReadMesh(const std::string& fileName)
{
	vtkAlgorithm *reader;
	try 
	{
		crl::FileName fn(fileName);
		if ( fn.getExtension()=="vtp" )
		{
			vtkXMLPolyDataReader *vtpReader = vtkXMLPolyDataReader::New();
			vtpReader->SetFileName(fileName.c_str());
			reader = vtpReader;
		}
		else if ( fn.getExtension()=="vtk" )
		{
			vtkPolyDataReader *vtkReader = vtkPolyDataReader::New();
			vtkReader->SetFileName(fileName.c_str());
			reader = vtkReader;
		}
		else
		{
			vtkGIFTIReader *giiReader = vtkGIFTIReader::New();
			giiReader->SetFileName(fileName.c_str());
			reader = giiReader;
		}

		reader->Update();
		vtkPolyData *polyData = dynamic_cast<vtkPolyData*>(reader->GetOutputDataObject(0));
		polyData->Register(NULL);
		polyData->SetSource(NULL);
		reader->Delete();
		return polyData;
	}
	catch (...)
	{
		std::cout << "Writing model failed." << std::endl;
		reader->Delete();
	}
	return NULL;
}

/**********************************************************************************************//**
 * \fn	bool crlVtkMeshIO::WriteMesh(vtkPolyData *mesh, const std::string& fileName, bool binary)
 *
 * \brief	Writes a mesh.
 *
 * \author	Benoit Scherrer
 * \date	July 2011
 *
 * \param [in,out]	mesh	The mesh.
 * \param	fileName		Filename of the file.
 * \param	binary			ONLY used for the giftii mode, binary mode or not.
 *
 * \return	true if it succeeds, false if it fails.
 **************************************************************************************************/
bool crlVtkMeshIO::WriteMesh(vtkPolyData *mesh, const std::string& fileName, bool binary)
{
	vtkAlgorithm *writer;
	try 
	{
		crl::FileName fn(fileName);
		if ( fn.getExtension()=="vtp" )
		{
			vtkXMLPolyDataWriter *vtpWriter = vtkXMLPolyDataWriter::New();
			vtpWriter->SetFileName(fileName.c_str());
			vtpWriter->SetInput(mesh);
			//vtpWriter->EncodeAppendedDataOff();
			vtpWriter->SetDataModeToBinary();
			writer = vtpWriter;
			vtpWriter->Update();
		}
		else if ( fn.getExtension()=="vtk" )
		{
			vtkPolyDataWriter *vtkWriter = vtkPolyDataWriter::New();
			vtkWriter->SetFileName(fileName.c_str());
			vtkWriter->SetInput(mesh);
			writer = vtkWriter;
			vtkWriter->Update();
		}
		else
		{
			vtkGIFTIWriter *giiWriter = vtkGIFTIWriter::New();
			giiWriter->SetFileName(fileName.c_str(), binary);
			giiWriter->SetInput(mesh);
			writer = giiWriter;
			giiWriter->Write();
		}
	}
	catch (...)
	{
		std::cout << "Writing model failed." << std::endl;
		writer->Delete();
		return false;
	}
	return true;
}
