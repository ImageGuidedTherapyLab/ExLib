#include "vtkGIFTIReader.h"

#include <vtkObjectFactory.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkTransform.h>
#include <vtkPolyDataNormals.h>
#include <vtkStripper.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkCleanPolyData.h>

#include <gifti_io.h>


vtkCxxRevisionMacro(vtkGIFTIReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkGIFTIReader);

vtkGIFTIReader::vtkGIFTIReader()
{
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
  fileNameSet = false;
  fileName = "";
}

vtkGIFTIReader::~vtkGIFTIReader()
{
}

void vtkGIFTIReader::SetFileName(std::string name)
{
  fileName = name;
  fileNameSet = true;
}

void vtkGIFTIReader::SetFileName(const char* name)
{
  fileName = name;
  fileNameSet = true;
}

int vtkGIFTIReader::RequestData(
				vtkInformation *vtkNotUsed(request),
				vtkInformationVector **inputVector,
				vtkInformationVector *outputVector)
{
  if (!fileNameSet)
    return 0;
  
  bool texture = false;
  bool lines = false;
  
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  
  // get the output
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  
  gifti_image *gim = gifti_read_image(fileName.c_str(),1);
  gifti_image *linegim = NULL;
  gifti_image *texgim = NULL;
  
  vtkDataArray *scalars = NULL;
  int numcomp = 0;
  
  if (!gim)
    {
      std::cout << "Reading GIFTI model failed." << std::endl;
      return false;
    }
  std::cout << "GIFTI version: " << gim->version << std::endl;
  std::cout << "Number of data arrays in file: " << gim->numDA << std::endl;
  
  // try to find a file with lines
  std::string lineFileName = (fileName.substr(0,fileName.find_last_of("."))+".lines.gii");
  std::ifstream linefile(lineFileName.c_str());
  if (linefile)
    {
      lines = true;
      linefile.close();
    }

  // now, try to find a file with scalars
  std::string textureFileName = (fileName.substr(0,fileName.find_last_of(".")))+".rgba.gii";
  std::ifstream texfile(textureFileName.c_str());
  if (texfile)
    {
      texture = true;
      texfile.close();
    }
  else
    {
      textureFileName = (fileName.substr(0,fileName.find_last_of(".")))+".scalars.gii";
      texfile.open(textureFileName.c_str());
      if (texfile)
	{
	  texture = true;
	  texfile.close();
	}
    }
  if (lines)
    {
      linegim = gifti_read_image(lineFileName.c_str(),1);
      if (!linegim)
	{
	  std::cout << "Reading GIFTI lines file failed." << std::endl;
	  std::cout << "Continuing without lines." << std::endl;
	  lines = false;
	}
      else
	{
	  std::cout << "Lines GIFTI version: " << linegim->version << std::endl;
	  std::cout << "Number of data arrays in lines file: " << linegim->numDA << std::endl;
	}
    }

  if (texture)
    {
      texgim = gifti_read_image(textureFileName.c_str(),1);
      
      if (!texgim)
	{
	  std::cout << "Reading GIFTI texture file failed." << std::endl;
	  std::cout << "Continuing without scalars." << std::endl;
	  texture = false;
	}
      else
	{
	  std::cout << "Texture GIFTI version: " << texgim->version << std::endl;
	  std::cout << "Number of data arrays in texture file: " << texgim->numDA << std::endl;
	}
    }
  
  // check if there are two dimensions and if the second dimension has three values (i.e. three coordinates per point)
  if ((gim->darray[0]->num_dim != 2) || (gim->darray[0]->dims[1] != 3))
    {
      std::cout << "Surface points have to be in 3D space!" << std::endl;
      return false;
    }
  
  int giftiPointData = -1;
  int giftiTriangleData = -1;
  
  for (int n = 0; n < gim->numDA; n++)
    {
      if (gim->darray[n]->intent == NIFTI_INTENT_POINTSET)
	giftiPointData = n;
      else if (gim->darray[n]->intent == NIFTI_INTENT_TRIANGLE)
	giftiTriangleData = n;
    }
  if ((giftiPointData < 0 || giftiTriangleData < 0) && !lines)
    {
      std::cout << "GIFTI model is missing point, triangle or line data." << std::endl;
      gifti_free_image(gim);
      if (texture) gifti_free_image(texgim);
      return 0;
    }
  
  
  unsigned long totalPoints,totalTriangles,totalLines;
  totalPoints = gim->darray[giftiPointData]->dims[0];
  if (giftiTriangleData >= 0)
    totalTriangles = gim->darray[giftiTriangleData]->dims[0];
  else
    totalTriangles = 0;
  if (lines)
    totalLines = linegim->darray[0]->dims[0];
  else
    totalLines = 0;
  if (texture)
    {
      unsigned long totalTexturePoints = texgim->darray[0]->dims[0];
      if (totalTexturePoints != totalPoints)
	{
	  std::cout << "Number of points in texture file must match number of points in GIFTI model." << std::endl;
	  gifti_free_image(gim);
	  gifti_free_image(texgim);
	  return 0;
	}
    }
  
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *triangleCells = vtkCellArray::New();
  vtkCellArray * lineCells = vtkCellArray::New();
  triangleCells->EstimateSize(totalTriangles,3);
  lineCells->EstimateSize(totalLines,3);
  
  vtkTransform *transform = vtkTransform::New();
  bool useTransform = false;
  
  transform->Identity();
  if (gim->darray[giftiPointData]->numCS > 0)
    {
      // we have a coordinate transform
      transform->SetMatrix(*(gim->darray[giftiPointData]->coordsys[0]->xform));
      useTransform = true;
    }
  
  float *data = (float*)gim->darray[giftiPointData]->data;
  void *texdata = NULL;
  if (texture)
    {
      texdata = texgim->darray[0]->data;
      numcomp = texgim->darray[0]->dims[1];
      if (numcomp == 0)
	scalars = vtkFloatArray::New();
      else if (numcomp == 4)
	scalars = vtkUnsignedCharArray::New();
      else
	{
	  std::cout << "This number of components (" << numcomp << ") is not supported." << std::endl;
	  gifti_free_image(gim);
	  gifti_free_image(texgim);
	  return 0;
	}
      scalars->SetNumberOfComponents(numcomp);
    }
  
  for (unsigned long n = 0; n < totalPoints; n++)
    {
      float point[4];
      point[0] = data[3*n];
      point[1] = data[3*n+1];
      point[2] = data[3*n+2];
      point[3] = 1;
      
      if (useTransform)
	transform->MultiplyPoint(point,point);
      points->InsertNextPoint(point[0],point[1],point[2]);
      
      if (texture)
	{
	  if (numcomp == 0)
	    {
	      float* floatdata = (float*)texdata;
	      dynamic_cast<vtkFloatArray*>(scalars)->InsertNextValue(floatdata[n]);
	    }
	  else if (numcomp == 4)
	    {
	      float* floatdata = (float*)texdata;
	      float scalar[4];
	      scalar[0] = floatdata[4*n]*255;
	      scalar[1] = floatdata[4*n+1]*255;
	      scalar[2] = floatdata[4*n+2]*255;
	      scalar[3] = floatdata[4*n+3]*255;
	      dynamic_cast<vtkUnsignedCharArray*>(scalars)->InsertNextTuple(scalar);
	    }
	}
    }
  transform->Delete();
  
  for (unsigned long n = 0; n < totalTriangles; n++)
    {
      unsigned int *data = (unsigned int*)gim->darray[giftiTriangleData]->data;
      unsigned int p1,p2,p3;
      p1 = data[3*n];
      p2 = data[3*n+1];
      p3 = data[3*n+2];
      vtkIdList *triangle = vtkIdList::New();
      triangle->InsertNextId(p1);
      triangle->InsertNextId(p2);
      triangle->InsertNextId(p3);
      triangleCells->InsertNextCell(triangle);
      triangle->Delete();
    }
  
  if (lines)
    for (unsigned long n = 0; n < totalLines; n++)
      {
	unsigned int *data = (unsigned int*)linegim->darray[0]->data;
	unsigned int p1,p2;
	p1 = data[2*n];
	p2 = data[2*n+1];
	vtkIdList *line = vtkIdList::New();
	line->InsertNextId(p1);
	line->InsertNextId(p2);
	lineCells->InsertNextCell(line);
	line->Delete();
      }

  vtkPolyData *polydata = vtkPolyData::New();
  polydata->SetPoints(points);
  polydata->SetPolys(triangleCells);
  if (lines)
    polydata->SetLines(lineCells);
  if (texture)
    polydata->GetPointData()->SetScalars(scalars);
  points->Delete();
  triangleCells->Delete();
  lineCells->Delete();
  
  if (!lines)
    {
      vtkPolyDataNormals *pdn = vtkPolyDataNormals::New();
      pdn->SetInput(polydata);
      pdn->SetFeatureAngle(60.0);
      vtkCleanPolyData *cleaner = vtkCleanPolyData::New();
      cleaner->SetInput(pdn->GetOutput());
      cleaner->Update();
      vtkStripper *strip = vtkStripper::New();
      strip->SetInput(cleaner->GetOutput());
      strip->Update();
  
      output->DeepCopy(strip->GetOutput());
      strip->Delete();
      cleaner->Delete();
      pdn->Delete();
    }
  else
    output->DeepCopy(polydata);
  
  polydata->Delete();
  
  gifti_free_image(gim);
  if (texture)
    gifti_free_image(texgim);
  if (lines)
    gifti_free_image(linegim);
  return 1;
}
