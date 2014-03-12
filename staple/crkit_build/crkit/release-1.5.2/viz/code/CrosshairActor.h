#ifndef _CROSSHAIRACTOR_H_
#define _CROSSHAIRACTOR_H_

class vtkCursor3D;
class vtkActor;
class vtkCursor2D;
class vtkActor2D;
class vtkSphereSource;

class CrosshairActor
{
public:
  CrosshairActor();
  ~CrosshairActor();

  vtkActor* GetActor();
  //vtkActor2D* GetActor2D();

  void SetWorldPosition(double x, double y, double z);
  void SetWorldPosition(double position[3]);
  double *GetWorldPosition();

  //void Set2DPosition(double x, double y);
  //int Get2DSlice();
  //void Set2DSlice(int s);

private:
  double pos[4];
  int slice;

  //vtkCursor3D* cursor;
  vtkActor* actor;
  //vtkCursor2D *cursor2D;
  //vtkActor2D* actor2D;

  vtkSphereSource *sphere;
};

#endif
