/* This maintains the state of a viewer window. */

#include <QtGui>
#include <QMainWindow>
#include <QSettings>
#include <QEvent>
#include <QDropEvent>
#include <QDrag>
#include <QString>

#include "ui_vizfourviewers.h"
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

class VizFourViewers : public VizViewerBase , private Ui::VizFourViewers
{
  Q_OBJECT

public:
  VizFourViewers(VizMainWindow* main);
  virtual ~VizFourViewers();

  bool DisplayImage(SceneImage* si, ImageTypeDefinitions::DataType t);
  bool DisplayOverlayImage(SceneImage* si, ImageTypeDefinitions::DataType t);
  bool DisplayTensors(SceneTensors *st);
  virtual void SetMarker(unsigned int _group, double* worldPos);
  virtual void SetMarkerInViewers(unsigned int _group, double *worldPos);
  virtual void CursorPositionHasChanged(double pos[3]) {};
  void UpdateTitle(unsigned int, QString);
  void SetGroup(unsigned int _group);
  virtual unsigned int GetGroup();
  virtual unsigned int GetNumberOfRenderersForObject(SceneData *sd);
  virtual void SaveModel(SceneData*,std::string) {};
  virtual void RotateModel(SceneData*) {};
  virtual void Render();


private slots:
// needs to initialize the slots for the menus
void on_actionExit_triggered();
  void Display_slice_changed1(int s);
  void Display_slice_changed2(int s);
  void Display_slice_changed3(int s);
  void Display_slice_changed4(int s);
  void actionSliceUp_triggered();
  void actionSliceDown_triggered();
  void on_actionRemove_Overlay_Image_triggered();
  void on_actionRemove_Tensors_triggered();
  void Blending_value_changed(int b);
  void UseGreyscale();
  void UseProbabilityMap();
  void UseSegmentation();
  void on_actionGroup1_triggered();
  void on_actionGroup2_triggered();
  void on_actionGroup3_triggered();
  void on_actionNo_Group_triggered();
  void on_actionDisplay_Cursor_triggered();
  void on_actionSave_Points_triggered();
  void on_actionClear_Points_triggered();
  void on_actionShow_Marker_triggered();
  void on_actionBlack_triggered();
  void on_actionWhite_triggered();
  void on_actionShow_Bounding_Box_triggered();
  void on_actionShow_Orientation_Cube_triggered();
  void on_actionSmall_Cube_triggered();

private:
  QSettings *appsettings;
  QAction *actionSliceUp;
  QAction *actionSliceDown;

  // A class for rendering images in 2D and 3D
  Render2D3D *r2d[4];

  QString name;
  QString overlay_name;

  unsigned int group;
  bool showMarker;
  bool showBoundingBox;
  bool showOrientationCube;
  bool smallCube;

  void Display_slice_changed(int viewer,int s);

  virtual void closeEvent(QCloseEvent *event);

  static void OnWindowLevelChangedEvent(vtkObject *caller,unsigned long event,void *data, void *calldata);

  // Support drops of file names to load images
  void dropEvent(QDropEvent *event);
  void dragEnterEvent(QDragEnterEvent *e);

  friend class VizMainWindow;
};
