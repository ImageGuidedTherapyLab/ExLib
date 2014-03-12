#include <QtGui>
#include <QDialog>

#include "ui_vizviewerwindowrotationdialog.h"


class Render2D3D;
class SceneModel;
class vtkTransform;
class vtkTransformPolyDataFilter;
class vtkPolyData;
class VizMainWindow;
class QCloseEvent;

class VizViewerWindowRotationDialog : public QDialog, private Ui::VizViewerWindowRotationDialog
{
  Q_OBJECT

public:
  VizViewerWindowRotationDialog();
  ~VizViewerWindowRotationDialog();


  void SetModel(SceneModel *m);
  void SetRenderer(Render2D3D *r);

private slots:
  void X_value_changed(int x);
  void Y_value_changed(int y);
  void Z_value_changed(int z);
  void save_clicked();

private:
  Render2D3D *ren;
  SceneModel *model;
  vtkTransform *transform;
  vtkTransformPolyDataFilter *filter;
  vtkPolyData *original;
  double centroid[3];
};
