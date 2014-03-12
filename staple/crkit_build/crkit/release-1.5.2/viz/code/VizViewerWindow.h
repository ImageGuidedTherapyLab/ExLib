/* This maintains the state of a viewer window. */

#include <QtGui>
#include <QMainWindow>
#include <QSettings>
#include <QEvent>
#include <QDropEvent>
#include <QDrag>
#include <QString>

#include "ui_vizviewer.h"
#include <set>

class QAction;
class QActionGroup;
class QLabel;
class QMenu;

// Viz functionality
#include <VizViewerBase.h>
#include "ImageTypeDefinitions.h"
#include "Render2D3D.h"
#include "SceneImage.h"
#include "ImageOperations.h"
#include <SceneTensors.h>

class VizViewerWindow : public VizViewerBase , private Ui::VizViewer
{
  Q_OBJECT

public:
  VizViewerWindow(VizMainWindow* main, bool stereo);
  virtual ~VizViewerWindow();

  bool DisplayImage(SceneData* si, ImageTypeDefinitions::DataType t);
  bool DisplayOverlayImage(SceneData* si, ImageTypeDefinitions::DataType t);
  bool DisplayTensors(SceneTensors *st);
  bool DisplayModelOnly(SceneModel *sm);
  virtual void SetMarker(unsigned int _group, double* worldPos);
  virtual void SetMarkerInViewers(unsigned int _group, double *worldPos);
  virtual void CursorPositionHasChanged(double pos[3]);
  void UpdateTitle(unsigned int, QString);
  void SetGroup(unsigned int _group);
  virtual unsigned int GetGroup();
  virtual unsigned int GetNumberOfRenderersForObject(SceneData *sd);
  virtual void SaveModel(SceneData *model,std::string fileName);
  virtual void RotateModel(SceneData *model);
  virtual void Render();

private slots:
// needs to initialize the slots for the menus
void on_actionExit_triggered();
  void Display_slice_changed(int s);
  void actionSliceUp_triggered();
  void actionSliceDown_triggered();
  void Display_orientation_xy();
  void Display_orientation_yz();
  void Display_orientation_xz();
  void on_actionRemove_Overlay_Image_triggered();
  void on_actionRemove_Tensors_triggered();
  void Blending_value_changed(int b);
  void Use2DInteractor();
  void Use3DInteractor();
  void UseGreyscale();
  void UseProbabilityMap();
  void UseSegmentation();
  void on_actionGroup1_triggered();
  void on_actionGroup2_triggered();
  void on_actionGroup3_triggered();
  void on_actionNo_Group_triggered();
  void on_actionDisplay_Cursor_triggered();
  void on_actionSave_Points_triggered();
  void on_actionLoad_Points_triggered();
  void on_actionClear_Points_triggered();
  void on_actionSelect_ROI_triggered();
  void on_actionSave_ROI_triggered();
  void on_actionSave_Screenshot_triggered();
  void on_actionSave_Movie_triggered();
  void on_actionSave_OBJ_File_triggered();
  void on_actionSave_RIB_File_triggered();
  void on_actionShow_Marker_triggered();
  void on_actionBlack_triggered();
  void on_actionWhite_triggered();
  void on_actionShow_Bounding_Box_triggered();
  void on_actionShow_Orientation_Cube_triggered();
  void on_actionSmall_Cube_triggered();
  void on_actionShow_Mouse_Pointer_triggered();
  void FullScreen();
  void on_actionSpecial_Label_Blending_triggered();
  void on_actionShow_Tool_triggered();
  void on_actionDisplay_in_VectorVision_triggered();

private:
  QSettings *appsettings;
  QAction *actionSliceUp;
  QAction *actionSliceDown;
  QAction *actionDisplayOrientationXY;
  QAction *actionDisplayOrientationYZ;
  QAction *actionDisplayOrientationXZ;
  QAction *actionFullScreen;

  // A class for rendering images in 2D and 3D
  Render2D3D *r2d;

  QString name;
  QString overlay_name;

  unsigned int group;
  bool showMarker;
  bool showBoundingBox;
  bool showOrientationCube;
  bool smallCube;
  bool showMousePointer;

  virtual void closeEvent(QCloseEvent *event);

  // Support drops of file names to load images
  void dropEvent(QDropEvent *event);
  void dragEnterEvent(QDragEnterEvent *e);

  friend class VizMainWindow;
};
