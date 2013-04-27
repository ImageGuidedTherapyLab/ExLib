
//
// Set up a QVTKWidget to render in 2D and 3D mode.
//
// Provide an API for
// 1. initialization of the rendering, the world coordinate system origin etc.
// 2. changing the data being rendered (setInput).

#ifndef _RENDER2D3D_INCLUDED
#define _RENDER2D3D_INCLUDED 1

#include <vtkObject.h>
#include <vtkCommand.h>
#include "ImageTypeDefinitions.h"
#include <set>

class vtkRenderer;
class vtkRenderWindow;
class QVTKInteractor;
class vtkImageActor;
class vtkImageMapToColors;
class vtkLookupTable;
class vtkImageBlend;
class QVTKWidget;
class SceneData;
class SceneImage;
class SceneImageView;
class SceneModel;
class VizViewerBase;
class CrosshairActor;
class SceneTensors;
class vtkActor;
class vtkAnnotatedCubeActor;
class vtkOrientationMarkerWidget;
class vtkBoxWidget;
class vtkActor2D;
class vtkCursor3D;
class vtkImageBlend;
class vtkMatrix4x4;

struct PointType
{
  int x,y,z;
  double pixelValue, overlayValue;
  unsigned int pointClass;
  CrosshairActor *actor;
  double world[4];
};

enum vizEventIDs
  {
    NoEvent = vtkCommand::UserEvent,
    Render2D3DWindowLevelEvent
  };

class Render2D3D : public vtkObject
{

public:
  vtkTypeRevisionMacro(Render2D3D,vtkObject);
  static Render2D3D *New();
  vtkRenderer *ren;
  QVTKWidget *qvtk;
  VizViewerBase *viewWin;
  bool stereo;

protected:
  // Data related to caching of properties of the ITK image data viewed
  ImageTypeDefinitions::ImageType::InternalPixelType level;
  ImageTypeDefinitions::ImageType::InternalPixelType window;
  ImageTypeDefinitions::ImageType::InternalPixelType min;
  ImageTypeDefinitions::ImageType::InternalPixelType max;
  ImageTypeDefinitions::ImageType::InternalPixelType mean;
  ImageTypeDefinitions::ImageType::InternalPixelType stddev;

public:

  inline bool Initialized();

  void SetViewerWindow(VizViewerBase *_view);
  bool InitializeRendering(QVTKWidget *renderonthis);

  bool SetScalarInput(int number, SceneData *si);
  bool SetTensorInput(SceneTensors *st);
  bool AddModel(SceneModel *sm);

  void Render(void) const;

  void ResetCamera();
  void GetCameraParameters(double *position, double *focalPoint, double *viewUp);
  void SetCameraParameters(double *position, double *focalPoint, double *viewUp);

  void SetSliceOrientationToXY();
  void SetSliceOrientationToYZ();
  void SetSliceOrientationToXZ();
  std::string GetSliceOrientation();
  void GetSliceRange(int &min, int &max);
  void SetSlice(int sliceindex);
  int GetSlice();
  void SetVisibleSlices(bool xy, bool yz, bool xz);

  void SetCrosshairPosition(double *worldPos);
  void UpdateCrosshairActor();
  void SetBlendValue(double b);

  void Set2DInteractor();
  void Set3DInteractor();
  bool Using3DInteractor();

  void SetDataType(int number,ImageTypeDefinitions::DataType t);

  void SetCursor(bool status);
  void SetCursorPosition(double pos[3]);

  void ShowMarker(bool status);
  void ShowBoundingBox(bool status);
  void ShowOrientationCube(double size);
  void ShowMousePointer(bool status);
  void SetSpecialLabelBlending(bool status);
  void ShowIGTTool(bool status);

  void SavePoints(std::string filename);
  void LoadPoints(std::string filename);
  void ClearPoints();

  void GetWindowLevel(double wl[4]);
  void SetWindowLevel(double wl[4]);

  void EnableROIBox();
  void DisableROIBox();
  void SaveROI(std::string filename);

  void SaveScreenshot(std::string filename);
  void SaveRIB(std::string filename);
  void SaveOBJ(std::string filename);
  void SaveMovie(std::string filename);

  void SetBackground(double r,double g,double b);

  void ToggleFullScreen();

  unsigned int GetNumberOfRenderersForObject(SceneData *sd);
  void SaveModel(SceneData *model,std::string fileName);
  std::set<SceneModel*> GetModels();

  void SetToolPosition(vtkMatrix4x4 *position, bool isValid);


protected:

  Render2D3D();
  ~Render2D3D();

private:

  static void On2DRightButtonPressEvent(vtkObject *caller, unsigned long event, void *data, void *calldata);
  static void On3DRightButtonPressEvent(vtkObject *caller, unsigned long event, void *data, void *calldata);
  static void On2DLeftButtonPressEvent(vtkObject *caller, unsigned long event, void *data, void *calldata);
  static void On3DMiddleButtonPressEvent(vtkObject *caller, unsigned long event, void *data, void *calldata);
  static void OnRenderStartEvent(vtkObject* caller,unsigned long event,void* data,void* calldata);
  static void WindowLevelCallback(vtkObject *caller, unsigned long event, void *data, void *calldata);
  static void On2DMouseMoveEvent(vtkObject *caller, unsigned long event, void *data, void *calldata);

  bool GetVoxelValueAtPosition(double *worldPos, float *imageValue, float *overlayValue, int *ijkPos=0);
  void DisplayVoxelValueAtPosition(double *world);
  void AddOrDeletePoint(double *world);
  void AdjustTensorActorPosition();
  void ComputeTrueWorldPosition(double *worldPos, double *trueWorldPos);
  void FixDepthSort();

  bool crosshairVisible;
  double wlInitialWindow, wlInitialLevel, wlInitialWindowOverlay, wlInitialLevelOverlay;

  QVTKInteractor *interactor;
  ImageTypeDefinitions::OrientationType orientation;

  SceneImageView *siv;
  SceneImageView *overlay_siv;

  vtkImageActor *vtkactor[3];
  vtkImageBlend *vtkblend;
  vtkImageMapToColors *vtkmapper, *vtkoverlaymapper;
  int lastSlice[3];
  vtkImageData *imageData, *overlayImageData;
  vtkActor *outlineActor;
  vtkAnnotatedCubeActor *cube;
  vtkOrientationMarkerWidget *orientationMarkers;
  vtkBoxWidget *roiBox;
  vtkActor *igtTool;

  SceneTensors *tensors;
  vtkActor *tensorActor;

  vtkCursor3D *cursor;
  vtkActor *cursorActor;
  bool cursorVisible;
  bool mousePointerVisible;
  unsigned int currentPointClass;
  float currentColor[3];
  bool labelBlending;

  bool in3D;
  std::map<SceneModel*,std::vector<vtkActor*> > modelActors;

  ImageTypeDefinitions::DataType scalarDataType[2];

  std::vector<PointType> points;
};

#endif
