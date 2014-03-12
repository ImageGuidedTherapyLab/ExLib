#include "vtkGIFTIWriter.h"

#include <vtkObjectFactory.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkTransform.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkTriangleFilter.h>

#include <gifti_io.h>


vtkCxxRevisionMacro(vtkGIFTIWriter, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkGIFTIWriter);

vtkGIFTIWriter::vtkGIFTIWriter()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(0);
  fileNameSet = false;
  fileName = "";
  gim = NULL;
  texgim = NULL;
  linegim = NULL;
  scalarComponents = 0;
  numTriangles = 0;
  numLines = 0;
  writeScalars = true;
}

vtkGIFTIWriter::~vtkGIFTIWriter()
{
  if (gim)
    delete gim;
  if (texgim)
    delete texgim;
  if (linegim)
    delete linegim;
}

void vtkGIFTIWriter::SetFileName(std::string name, bool writeBinary)
{
  fileName = name;
  fileNameWithoutExt = fileName.substr(0,fileName.find_last_of("."));
  binary = writeBinary;
  fileNameSet = true;
}

void vtkGIFTIWriter::SetWriteScalars(bool write)
{
  writeScalars = write;
}

void vtkGIFTIWriter::SetInput(vtkDataObject *input)
{
  vtkPolyData *orig = dynamic_cast<vtkPolyData*>(input);
  if (!orig)
    return;

  vtkTriangleFilter *triangulator = vtkTriangleFilter::New();
  triangulator->SetInput(orig);
  triangulator->Update();

  vtkPolyData *model = triangulator->GetOutput();
  model->BuildCells();
  model->BuildLinks();

  gim = gifti_create_image(2,0,NIFTI_TYPE_FLOAT32,0,0,0);

  gim->darray[0]->intent = NIFTI_INTENT_POINTSET;
  gim->darray[1]->intent = NIFTI_INTENT_TRIANGLE;

  gim->darray[0]->num_dim = 2;
  gim->darray[0]->dims[0] = model->GetNumberOfPoints();
  gim->darray[0]->dims[1] = 3;
  if (binary)
    {
      gim->darray[0]->encoding = GIFTI_ENCODING_EXTBIN;
      gim->darray[0]->endian = GIFTI_ENDIAN_LITTLE;
      gim->darray[0]->ext_offset = 0;
      gim->darray[0]->ext_fname = gifti_strdup((fileNameWithoutExt+".giibin").c_str());;
    }
  else
    gim->darray[0]->encoding = GIFTI_ENCODING_B64GZ;

  // check number of triangles and lines
  for (long n = 0; n < model->GetNumberOfCells(); n++)
    {
      if (model->GetCellType(n) == VTK_TRIANGLE)
	numTriangles++;
      else if (model->GetCellType(n) == VTK_LINE)
	numLines++;
    }
  std::cout << "Triangles: " << numTriangles << std::endl;
  std::cout << "Lines: " << numLines << std::endl;

  gim->darray[1]->num_dim = 2;
  gim->darray[1]->datatype = NIFTI_TYPE_INT32;
  gim->darray[1]->dims[0] = numTriangles;
  gim->darray[1]->dims[1] = 3;
  if (binary)
    {
      gim->darray[1]->encoding = GIFTI_ENCODING_EXTBIN;
      gim->darray[1]->endian = GIFTI_ENDIAN_LITTLE;
      gim->darray[1]->ext_offset = 0;
      gim->darray[1]->ext_fname = gifti_strdup((fileNameWithoutExt+".triangles.giibin").c_str());
    }
  else
    gim->darray[1]->encoding = GIFTI_ENCODING_B64GZ;

  if (numLines > 0)
    {
      linegim = gifti_create_image(1,0,NIFTI_TYPE_INT32,0,0,0);
      linegim->darray[0]->intent = NIFTI_INTENT_NONE;
      linegim->darray[0]->num_dim = 2;
      linegim->darray[0]->dims[0] = numLines;
      linegim->darray[0]->dims[1] = 2;
      if (binary)
	{
	  linegim->darray[0]->encoding = GIFTI_ENCODING_EXTBIN;
	  linegim->darray[0]->endian = GIFTI_ENDIAN_LITTLE;
	  linegim->darray[0]->ext_offset = 0;
	  linegim->darray[0]->ext_fname = gifti_strdup((fileNameWithoutExt+".lines.giibin").c_str());
	}
      else
	linegim->darray[0]->encoding = GIFTI_ENCODING_B64GZ;
    }
  if (numTriangles == 0)
    {
      gifti_clear_DataArray(gim->darray[1]);
      gifti_free_DataArray(gim->darray[1]);
      gim->numDA = 1;
    }

  // fill points
  double point[4];
  float *data = (float*)malloc(model->GetNumberOfPoints()*3*sizeof(float));
  gim->darray[0]->data = data;
  gim->darray[0]->nvals = model->GetNumberOfPoints()*3;
  gim->darray[0]->nbyper = sizeof(float);
  gim->darray[0]->numCS = 0;
  for (long n = 0; n < gim->darray[0]->dims[0]; n++)
    {
      model->GetPoint(n,point);
      data[3*n] = point[0];
      data[3*n+1] = point[1];
      data[3*n+2] = point[2];
    }
  // fill triangles and lines
  vtkIdType *p = 0,np;
  int *data2 = 0, *data3 = 0;
  if (numTriangles > 0)
    {

      data2 = (int*)malloc(numTriangles*3*sizeof(int));
      gim->darray[1]->data = data2;
      gim->darray[1]->nvals = numTriangles*3;
      gim->darray[1]->nbyper = sizeof(int);
      gim->darray[1]->numCS = 0;
    }

  if (linegim)
    {
      data3 = (int*)malloc(numLines*2*sizeof(int));
      linegim->darray[0]->data = data3;
      linegim->darray[0]->nvals = numLines*2;
      linegim->darray[0]->nbyper = sizeof(int);
      linegim->darray[0]->numCS = 0;
    }
  for (long n = 0; n < numTriangles; n++)
    {
      if (model->GetCellType(n) == VTK_TRIANGLE)
	{
	  model->GetCellPoints(n,np,p);
	  data2[3*n] = p[0];
	  data2[3*n+1] = p[1];
	  data2[3*n+2] = p[2];
	}
    }
  for (long n = 0; n < numLines; n++)
    {
      if (model->GetCellType(n) == VTK_LINE)
	{
	  model->GetCellPoints(n,np,p);
	  data3[2*n] = p[0];
	  data3[2*n+1] = p[1];
	}
    }

  if (model->GetPointData()->GetScalars() && model->GetPointData()->GetScalars()->GetNumberOfComponents() == 3)
    {
      vtkUnsignedCharArray *scalars = (vtkUnsignedCharArray*)(model->GetPointData()->GetScalars());
      texgim = gifti_create_image(1,0,NIFTI_TYPE_FLOAT32,0,0,0);
      texgim->darray[0]->intent = NIFTI_INTENT_RGBA_VECTOR;
      texgim->darray[0]->num_dim = 2;
      texgim->darray[0]->dims[0] = model->GetNumberOfPoints();
      texgim->darray[0]->dims[1] = 4;
      float *texdata = (float*)malloc(model->GetNumberOfPoints()*4*sizeof(float));
      texgim->darray[0]->data = texdata;
      texgim->darray[0]->nvals = model->GetNumberOfPoints()*4;
      texgim->darray[0]->nbyper = sizeof(float);
      texgim->darray[0]->numCS = 0;
      if (binary)
	{
	  texgim->darray[0]->encoding = GIFTI_ENCODING_EXTBIN;
	  texgim->darray[0]->endian = GIFTI_ENDIAN_LITTLE;
	  texgim->darray[0]->ext_offset = 0;
	  texgim->darray[0]->ext_fname = gifti_strdup((fileNameWithoutExt+".rgba.giibin").c_str());
	}
      else
	texgim->darray[0]->encoding = GIFTI_ENCODING_B64GZ;
      unsigned char rgb[4];
      for (long n = 0; n < texgim->darray[0]->dims[0]; n++)
	{
	  scalars->GetTupleValue(n,rgb);
	  texdata[4*n] = (float)(rgb[0])/255.0;
	  texdata[4*n+1] = (float)(rgb[1])/255.0;
	  texdata[4*n+2] = (float)(rgb[2])/255.0;
	  texdata[4*n+3] = 1.0;
	}
      scalarComponents = 4;
    }
  else if (model->GetPointData()->GetScalars() && model->GetPointData()->GetScalars()->GetNumberOfComponents() == 1)
    {
      vtkFloatArray *scalars = (vtkFloatArray*)(model->GetPointData())->GetScalars();
      texgim = gifti_create_image(1,0,NIFTI_TYPE_FLOAT32,0,0,0);
      texgim->darray[0]->intent = NIFTI_INTENT_NONE;
      texgim->darray[0]->num_dim = 1;
      texgim->darray[0]->dims[0] = model->GetNumberOfPoints();
      //texgim->darray[0]->dims[1] = 1;
      float *texdata = (float*)malloc(model->GetNumberOfPoints()*sizeof(float));
      texgim->darray[0]->data = texdata;
      texgim->darray[0]->nvals = model->GetNumberOfPoints();
      texgim->darray[0]->nbyper = sizeof(float);
      texgim->darray[0]->numCS = 0;
      if (binary)
	{
	  texgim->darray[0]->encoding = GIFTI_ENCODING_EXTBIN;
	  texgim->darray[0]->endian = GIFTI_ENDIAN_LITTLE;
	  texgim->darray[0]->ext_offset = 0;
	  texgim->darray[0]->ext_fname = gifti_strdup((fileNameWithoutExt+".scalars.giibin").c_str());
	}
      else
	texgim->darray[0]->encoding = GIFTI_ENCODING_B64GZ;
      float value;
      for (long n = 0; n < texgim->darray[0]->dims[0]; n++)
	{
	  value = scalars->GetValue(n);
	  texdata[n] = value;
	}
      scalarComponents = 1;
    }

}

void vtkGIFTIWriter::WriteData()
{
  if ((!gim) || (!fileNameSet))
    return;

  gifti_write_image(gim,fileName.c_str(),1);

  if (texgim && writeScalars && (scalarComponents == 4))
    gifti_write_image(texgim,(fileNameWithoutExt+".rgba.gii").c_str(),1);
  else if (texgim && writeScalars && (scalarComponents == 1))
    gifti_write_image(texgim,(fileNameWithoutExt+".scalars.gii").c_str(),1);

  if (linegim) gifti_write_image(linegim,(fileNameWithoutExt+".lines.gii").c_str(),1);
}

int vtkGIFTIWriter::Write()
{
  WriteData();
  return 0;
}
