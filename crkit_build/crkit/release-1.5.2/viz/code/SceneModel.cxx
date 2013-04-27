#include <SceneModel.h>

#include <SceneImageView.h>
#include <crlVtkMeshIO.h>

// VTK classes
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkPointData.h>
#include <vtkLookupTable.h>
#include <vtkIdList.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>


SceneModel::SceneModel()
{
  polyData = 0;
}

// copy constructor
SceneModel::SceneModel(const SceneModel &o) {
  name = o.GetName();
  polyData = o.polyData;
  lowerValue = 0;
  upperValue = 0;
  opacity = 1.0;
}

// Construct with a specified handle (name)
SceneModel::SceneModel(const QString newName)
{
  name = newName;
  polyData = 0;
  lowerValue = 0;
  upperValue = 0;
  opacity = 1.0;
}

// Construct with a specified handle (name) and load an image
SceneModel::SceneModel(QString newName, QString loadFileName)
{
  name = newName;
  polyData = 0;
  if (!LoadModel(loadFileName)) {
    fileName = "";
  }
  lowerValue = 0;
  upperValue = 0;
  opacity = 1.0;
}

// destructor
SceneModel::~SceneModel() {
  if (polyData)
    polyData->Delete();
  polyData = 0;

  std::vector<vtkActor*>::iterator it;
  for (it = modelActors.begin(); it != modelActors.end(); it++)
    (*it)->Delete();
  modelActors.clear();
}

void SceneModel::SetName(const char *n1)
{
  name = QString(n1);
}

void SceneModel::SetName(QString n1)
{
  name = n1;
}

QString SceneModel::GetName() const
{
  return name;
}

bool SceneModel::LoadModel(QString loadFileName, bool hasScalars)
{
  fileName = loadFileName;
  isFlippedX = false;
  isFlippedY = false;
  isFlippedZ = false;
  ignoreScalars = !hasScalars;
  logLookupTable = false;

  polyData = crlVtkMeshIO::ReadMesh(loadFileName.toStdString());
  /*vtkAlgorithm *reader;
  if (loadFileName.contains(".vtp",Qt::CaseInsensitive))
    {
      reader = vtkXMLPolyDataReader::New();
      dynamic_cast<vtkXMLPolyDataReader*>(reader)->SetFileName(loadFileName.toStdString().c_str());
    }
  else if (loadFileName.contains(".gii",Qt::CaseInsensitive))
    {
      reader = vtkGIFTIReader::New();
      dynamic_cast<vtkGIFTIReader*>(reader)->SetFileName(loadFileName.toStdString().c_str());
    }
  else
    {
      reader = vtkPolyDataReader::New();
      dynamic_cast<vtkPolyDataReader*>(reader)->SetFileName(loadFileName.toStdString().c_str());
    }
  try
    {
      reader->Update();
    }
  catch (...)
    {
      std::cout << "Reading model failed." << std::endl;
      reader->Delete();
      return false;
    }
    polyData = dynamic_cast<vtkPolyData*>(reader->GetOutputDataObject(0));*/
  // register the vtkPolyData object to itself to prevent it from being deleted
  polyData->Register(polyData);
  //reader->Delete();

  color = Qt::lightGray;
  lowerColor.setHsvF(0.0,1.0,1.0,1.0);
  upperColor.setHsvF(0.24,1.0,1.0,1.0);
  return true; // success
}

vtkPolyData *SceneModel::GetPolyData() {
  return polyData;
}

void SceneModel::SetPolyData(vtkPolyData *data)
{
  polyData->DeepCopy(data);
}

vtkActor *SceneModel::CreateNewActor()
{
  if (!polyData)
    return NULL;

  vtkActor *actor = vtkActor::New();
  vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
  if (ignoreScalars)
    polyData->GetPointData()->SetScalars(NULL);

  mapper->SetInput(polyData);

  if (!ignoreScalars)
    {
      vtkLookupTable *vtklut = vtkLookupTable::New();
      vtklut->SetTableRange(polyData->GetScalarRange()[0],polyData->GetScalarRange()[1]);
      std::cout << polyData->GetScalarRange()[0] << " " << polyData->GetScalarRange()[1] << std::endl;
      lowerValue = polyData->GetScalarRange()[0];
      upperValue = polyData->GetScalarRange()[1];

      if (logLookupTable)
	vtklut->SetScaleToLog10();
      vtklut->SetHueRange(lowerColor.hueF(),upperColor.hueF());
      vtklut->SetSaturationRange(1,1);
      vtklut->SetValueRange(1,1);
      vtklut->Build();
      //vtklut->SetTableValue(0,0.9,0.9,0.9);
      //if (polyData->GetScalarRange()[0] == 0)
      //vtklut->SetTableValue(0,1.0,1.0,1.0);
      mapper->UseLookupTableScalarRangeOn();
      mapper->SetLookupTable(vtklut);
      vtklut->Delete();
      //mapper->SetColorModeToMapScalars();
    }

  mapper->Update();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(color.redF(),color.greenF(),color.blueF());
  actor->GetProperty()->BackfaceCullingOff();
  actor->GetProperty()->FrontfaceCullingOff();
  actor->GetProperty()->SetOpacity(opacity);
  modelActors.push_back(actor);
  mapper->Delete();
  return actor;
}


vtkActor *SceneModel::CreateNewActor(SceneImageView *siv)
{
  if (!polyData)
    return NULL;

  vtkActor *actor = this->CreateNewActor();

  vtkPolyData *pd = dynamic_cast<vtkPolyData*>(actor->GetMapper()->GetInput());

  vtkTransform *transform = vtkTransform::New();
  vtkTransformPolyDataFilter *tpdf = vtkTransformPolyDataFilter::New();
  tpdf->SetInput(pd);
  transform->PostMultiply();
  transform->Translate(-(siv->GetTranslation()[0]),-(siv->GetTranslation()[1]),-(siv->GetTranslation()[2]));
  transform->Concatenate(siv->GetRotation());
  transform->Translate(siv->GetTranslation());

  tpdf->SetTransform(transform);
  dynamic_cast<vtkPolyDataMapper*>(actor->GetMapper())->SetInput(tpdf->GetOutput());

  transform->Delete();
  tpdf->Delete();

  return actor;
}


void SceneModel::AddTextureMap(vtkImageData* texture)
{
  /*bool del = false;
  // delete the texture if we get the same input again
  if (probe && (texture == probe->GetSource()))
  del = true;

  vtkProbeFilter* oldProbe = probe;

  if (!del)
  {
  probe = vtkProbeFilter::New();
  probe->SpatialMatchOff();
  probe->SetInput(polyData);
  probe->SetSource(texture);
  probe->SpatialMatchOn();
  probe->Update();
  std::cout << probe->GetValidPoints()->GetNumberOfTuples() << " out of " << probe->GetOutput()->GetNumberOfPoints() << std::endl;

  std::vector<vtkActor*>::iterator it;
  for (it = modelActors.begin(); it != modelActors.end(); it++)
  {
  vtkPolyDataMapper *mapper = dynamic_cast<vtkPolyDataMapper*>((*it)->GetMapper());
  mapper->SetInput(probe->GetPolyDataOutput());
  vtkLookupTable *vtklut = vtkLookupTable::New();
  vtklut->SetTableRange(0.001,5.0);
  vtklut->SetHueRange(0.0,0.18);
  vtklut->SetSaturationRange(1,1);
  vtklut->SetValueRange(1,1);
  vtklut->Build();
  vtklut->SetTableValue(0,0.9,0.9,0.9);
  mapper->SetLookupTable(vtklut);
  vtklut->Delete();
  }
  }
  else
  {
  probe = 0;
  std::vector<vtkActor*>::iterator it;
  for (it = modelActors.begin(); it != modelActors.end(); it++)
  dynamic_cast<vtkPolyDataMapper*>((*it)->GetMapper())->SetInput(polyData);
  }

  if (oldProbe)
  oldProbe->Delete();*/
}


void SceneModel::ToggleLookupTable()
{
  if (ignoreScalars)
    return;

  logLookupTable = !logLookupTable;
  std::vector<vtkActor*>::iterator it;

  for (it = modelActors.begin(); it != modelActors.end(); it++)
    {
      if (logLookupTable)
	dynamic_cast<vtkLookupTable*>((*it)->GetMapper()->GetLookupTable())->SetScaleToLog10();
      else
	dynamic_cast<vtkLookupTable*>((*it)->GetMapper()->GetLookupTable())->SetScaleToLinear();
    }
}

QColor SceneModel::GetColor(ColorType whichOne)
{
  if (whichOne == lower)
    return lowerColor;
  else if (whichOne == upper)
    return upperColor;
  else
    return color;
}

void SceneModel::SetColor(QColor newColor, ColorType whichOne)
{
  if (whichOne == noScalars)
    {
      color = newColor;

      std::vector<vtkActor*>::iterator it;
      for (it = modelActors.begin(); it != modelActors.end(); it++)
	{
	  (*it)->GetProperty()->SetColor(color.redF(),color.greenF(),color.blueF());
	}
    }
  else
    {
      if (whichOne == lower)
	lowerColor = newColor;
      else
	upperColor = newColor;

      std::vector<vtkActor*>::iterator it;
      for (it = modelActors.begin(); it != modelActors.end(); it++)
	{
	  vtkLookupTable *lut = vtkLookupTable::New();
	  lut->SetHueRange(lowerColor.hueF(),upperColor.hueF());
	  lut->SetSaturationRange(lowerColor.saturationF(),upperColor.saturationF());
	  lut->SetValueRange(lowerColor.valueF(),upperColor.valueF());
	  if (logLookupTable)
	    lut->SetScaleToLog10();
	  else
	    lut->SetScaleToLinear();
	  if ((lowerValue!=0) || (upperValue!=0))
	    lut->SetTableRange(lowerValue,upperValue);

	  lut->Build();
	  (*it)->GetMapper()->SetLookupTable(lut);
	  lut->Delete();
	}
    }
}

bool SceneModel::GetLogLookupTable()
{
  return logLookupTable;
}

void SceneModel::SetUpperValue(double value)
{
  upperValue = value;
  this->SetColor(upperColor,upper);
}

void SceneModel::SetLowerValue(double value)
{
  lowerValue = value;
  this->SetColor(lowerColor,lower);
}

double SceneModel::GetLowerValue()
{
  return lowerValue;
}

double SceneModel::GetUpperValue()
{
  return upperValue;
}

void SceneModel::SetOpacity(double value)
{
  opacity = value;
  std::vector<vtkActor*>::iterator it;
  for (it = modelActors.begin(); it != modelActors.end(); it++)
    {
      (*it)->GetProperty()->SetOpacity(value);
    }
}

double SceneModel::GetOpacity()
{
  return opacity;
}
