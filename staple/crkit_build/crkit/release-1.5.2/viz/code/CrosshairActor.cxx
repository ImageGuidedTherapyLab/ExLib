#include <CrosshairActor.h>

#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>

CrosshairActor::CrosshairActor()
{
  pos[0] = 0;
  pos[1] = 0;
  pos[2] = 0;
  pos[3] = 1;
  slice = 0;

  actor = vtkActor::New();

  sphere = vtkSphereSource::New();
  sphere->SetCenter(0,0,0);
  sphere->SetRadius(2.0);

  vtkPolyDataMapper *pdm = vtkPolyDataMapper::New();
  pdm->SetInput(sphere->GetOutput());

  actor->SetMapper(pdm);
  pdm->Delete();

  actor->GetProperty()->SetColor(1.0,0.0,0.0);
}

CrosshairActor::~CrosshairActor()
{
  actor->Delete();
  sphere->Delete();
}

vtkActor* CrosshairActor::GetActor()
{
  return actor;
}

void CrosshairActor::SetWorldPosition(double x,double y,double z)
{
  sphere->SetCenter(x,y,z);
  pos[0] = x;
  pos[1] = y;
  pos[2] = z;
}

void CrosshairActor::SetWorldPosition(double position[3])
{
  sphere->SetCenter(position);

  pos[0] = position[0];
  pos[1] = position[1];
  pos[2] = position[2];
}

double *CrosshairActor::GetWorldPosition()
{
  return pos;
}
