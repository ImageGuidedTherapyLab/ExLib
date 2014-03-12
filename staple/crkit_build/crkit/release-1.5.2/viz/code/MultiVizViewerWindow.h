/* This maintains the state of a viewer window. */

#include <QtGui>
#include <QMainWindow>
#include <QSettings>
#include <QEvent>
#include <QDropEvent>
#include <QDrag>
#include <QString>

#include "ui_multivizviewer.h"
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

class MultiVizViewerWindow : public VizViewerBase , private Ui::MultiVizViewer
{
  Q_OBJECT

public:
  MultiVizViewerWindow();
  virtual ~MultiVizViewerWindow();

  bool LoadImage(std::string name, ImageTypeDefinitions::DataType t);
  bool DisplayImage(SceneData* si, ImageTypeDefinitions::DataType t,int number);
  virtual void SetMarker(unsigned int _group, double* worldPos);
  virtual void SetMarkerInViewers(unsigned int _group, double *worldPos);
  virtual void CursorPositionHasChanged(double pos[3]) {};
  virtual unsigned int GetGroup() {return 1;};
  virtual unsigned int GetNumberOfRenderersForObject(SceneData *sd) {return 0;};
  virtual void SaveModel(SceneData*,std::string) {};
  virtual void RotateModel(SceneData*) {};
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
  void Blending_value_changed(int b);
  void UseGreyscale();
  void UseProbabilityMap();
  void UseSegmentation();
  void on_actionShow_Marker_triggered();
  void on_actionBlack_triggered();
  void on_actionWhite_triggered();

private:
  QAction *actionSliceUp;
  QAction *actionSliceDown;
  QAction *actionDisplayOrientationXY;
  QAction *actionDisplayOrientationYZ;
  QAction *actionDisplayOrientationXZ;
  QAction *actionFullScreen;

  Render2D3D* r2d[12];
  int imageCounter;

  bool showMarker;
  unsigned int group;

  virtual void closeEvent(QCloseEvent *event);
};
