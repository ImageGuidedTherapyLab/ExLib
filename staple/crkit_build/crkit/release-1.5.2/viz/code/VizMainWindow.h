
/* This maintains the state of the main window. */

#include <QtGui>
#include <QMainWindow>
#include <QSettings>
#include <QEvent>
#include <QDropEvent>
#include <QDrag>
#include <QString>

#include "ui_viz.h"

#include <SceneModel.h>
#include <SceneImage.h>
class QAction;
class QActionGroup;
class QLabel;
class QMenu;
class VizViewerBase;
class SceneData;
class Render2D3D;

// Viz functionality
#include "ImageTypeDefinitions.h"
#include "ImageOperations.h"

// VVLink functionality
#ifdef HAVE_VVLINK
class vtkWindowToImageFilter;
#include <VVLConnection.h>
#include <common/VVLCoordinateSystem.h>
#include <VVLNTrackedInstrument.h>
#include <VVLQSetView.h>
#include <common/VVLInstrument.h>
#endif

typedef std::vector<VizViewerBase*> ViewerVector;

enum eOrientation {xy,yz,xz};

struct ImageInformation
{
  std::string fileName;
  SceneData *sd;
  ImageTypeDefinitions::DataType dataType;
  QListWidgetItem* listItem;
  ViewerVector viewers;
};

typedef QHash<QString, ImageInformation> ImageInformationMap;

class VizMainWindow : public QMainWindow , private Ui::Viz
{
  Q_OBJECT

public:
  VizMainWindow();
  ~VizMainWindow();

  bool FileOpenWithName(QString s,ImageTypeDefinitions::DataType t, QString handle="");
  void OpenScene(QString s);

  void PropagateSliceChange(unsigned int group, int slice);
  void PropagateOrientationChange(unsigned int group, eOrientation o);
  void PropagateCursorPosition(unsigned int group, double pos[3]);
  void SetMarkerInViewers(unsigned int _group, double* worldPos);
  void RenderAllViewers();

  void RemoveViewerFromMap(VizViewerBase *viewer);
  ImageInformation GetImageInformationForHandle(QString handle);

  bool stereo, ignoreLastDirectory;

#ifdef HAVE_VVLINK
  bool DisplayInVectorVision(Render2D3D *image);
  void TrackToolOn(Render2D3D *viewer);
  void TrackToolOff(Render2D3D *viewer);
#endif

private slots:
// needs to initialize the slots for the menus
void on_actionExit_triggered();
  void on_actionAbout_Viz_triggered();
  void on_actionOpen_Scene_triggered();
  void on_actionSave_Scene_triggered();
  void on_actionAdd_Volume_triggered();
  void on_actionAdd_Probability_Map_triggered();
  void on_actionAdd_Segmented_Volume_triggered();
  void on_actionAdd_Tensor_Volume_triggered();
  void on_actionAdd_RGB_Volume_triggered();
  void on_actionAdd_Model_triggered();
  void on_actionAdd_Model_with_Scalars_triggered();
  void on_actionAdd_fMRI_Activation_Map_triggered();
  void on_actionConnect_to_VectorVision_triggered();
  void on_actionDisconnect_triggered();
  void on_actionGet_Image_Data_triggered();
  void on_btnDelete_clicked();
  void on_btnFourViewers_clicked();
  VizViewerBase* FileListDoubleClick(QListWidgetItem* item);
  void FileListRightClick(const QPoint &pos);
  void timerCallback();
  void toolsMenuCallback(QAction *action);

  // Wizards
  void on_actionSTAPLE_triggered();

private:
  QSettings *appsettings;

  // Keep track of all of the images that are in the scene.
  // It would be lighter weight to keep pointers to SceneImage objects,
  // rather than copies of the objects themselves.
  ImageInformationMap imagesInScene;

  virtual void closeEvent(QCloseEvent *event);

  // Support drops of file names to load images
  void dropEvent(QDropEvent *event);
  void dragEnterEvent(QDragEnterEvent *e);

  bool inSliceRecursion,inOrientationRecursion,inCursorRecursion;

  void ChooseModelColor(SceneModel *sm,SceneModel::ColorType whichOne);
  void ChooseImageColor(SceneImage *si,SceneImage::ColorType whichOne);

  void SaveScene(QString s);

#ifdef HAVE_VVLINK
  vvlink::VVLConnection *server;
  vvlink::VVLCoordinateSystem *coord;
  vvlink::VVLNTrackedInstrument *trackedToolNotifier;
  vvlink::VVLQSetView *view;
  vvlink::VVLInstrument *tool;
  Render2D3D *displayInVectorVision;
  std::vector<Render2D3D*> trackToolViewers;
  bool inTimer;
  vtkWindowToImageFilter *wtif;
  bool connected;
  std::vector<QAction*> tools;
  vtkMatrix4x4 *lastToolPosition;
#endif
};
