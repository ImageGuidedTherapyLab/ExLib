/* This maintains the state of a viewer window. */

#include <MultiVizViewerWindow.h>

#include <qsettings.h>
#include <QFileDialog>
#include <QStatusBar>
#include <QKeySequence>
#include <QAction>

// To support drag and drop
#include <QMimeData>
#include <QList>
#include <QUrl>

#include "version.h"
#include <VizMainWindow.h>
#include <SceneImage.h>

MultiVizViewerWindow::MultiVizViewerWindow()
{
  setupUi(this);

  this->mainWin = NULL;
  group = 1;
  imageCounter = 0;
  showMarker = false;

  // Configure the Render2D3D object
  for (int i = 0; i < 12; i++)
    {
      r2d[i] = Render2D3D::New();
      r2d[i]->SetViewerWindow(this);
      r2d[i]->stereo = false;
    }

  r2d[0]->InitializeRendering(qvtkWidget);
  r2d[1]->InitializeRendering(qvtkWidget_2);
  r2d[2]->InitializeRendering(qvtkWidget_3);
  r2d[3]->InitializeRendering(qvtkWidget_4);
  r2d[4]->InitializeRendering(qvtkWidget_5);
  r2d[5]->InitializeRendering(qvtkWidget_6);
  r2d[6]->InitializeRendering(qvtkWidget_7);
  r2d[7]->InitializeRendering(qvtkWidget_8);
  r2d[8]->InitializeRendering(qvtkWidget_9);
  r2d[9]->InitializeRendering(qvtkWidget_10);
  r2d[10]->InitializeRendering(qvtkWidget_11);
  r2d[11]->InitializeRendering(qvtkWidget_12);


  // Configure keyboard shortcut actions for the qvtkWidget
  actionSliceDown = new QAction(this);
  actionSliceDown->setShortcutContext(Qt::WindowShortcut);
  QList<QKeySequence> shortCutList;
  shortCutList.append( QKeySequence(Qt::Key_Comma) );
  shortCutList.append( QKeySequence(Qt::Key_Left)  );
  shortCutList.append( QKeySequence(Qt::Key_Down)  );
  actionSliceDown->setShortcuts( shortCutList );
  connect( actionSliceDown, SIGNAL(triggered()),
	   this, SLOT(actionSliceDown_triggered()));
  this->addAction( actionSliceDown );

  actionSliceUp = new QAction(this);
  actionSliceUp->setShortcutContext(Qt::WindowShortcut);
  shortCutList[0] = QKeySequence(Qt::Key_Period);
  shortCutList[1] = QKeySequence(Qt::Key_Right);
  shortCutList[2] = QKeySequence(Qt::Key_Up);
  actionSliceUp->setShortcuts( shortCutList );
  this->addAction( actionSliceUp );
  connect( actionSliceUp, SIGNAL(triggered()),
	   this, SLOT(actionSliceUp_triggered()));

  sliceScrollBar->setMinimum(0);
  sliceScrollBar->setMaximum(255);
  sliceScrollBar->setValue(255);
  connect( sliceScrollBar, SIGNAL(valueChanged(int)),
	   this, SLOT(Display_slice_changed(int))
	   );

  connect( sliceSpinBox, SIGNAL(valueChanged(int)),
	   this, SLOT(Display_slice_changed(int))
	   );

  actionDisplayOrientationXY = new QAction(this);
  actionDisplayOrientationXY->setShortcutContext(Qt::WindowShortcut);
  actionDisplayOrientationXY->setShortcut(QKeySequence(Qt::Key_Z));
  this->addAction( actionDisplayOrientationXY);
  connect( actionDisplayOrientationXY, SIGNAL(triggered()),
	   this, SLOT(Display_orientation_xy()));

  actionDisplayOrientationYZ = new QAction(this);
  actionDisplayOrientationYZ->setShortcutContext(Qt::WindowShortcut);
  actionDisplayOrientationYZ->setShortcut(QKeySequence(Qt::Key_X));
  this->addAction( actionDisplayOrientationYZ);
  connect( actionDisplayOrientationYZ, SIGNAL(triggered()),
	   this, SLOT(Display_orientation_yz()));

  actionDisplayOrientationXZ = new QAction(this);
  actionDisplayOrientationXZ->setShortcutContext(Qt::WindowShortcut);
  actionDisplayOrientationXZ->setShortcut(QKeySequence(Qt::Key_Y));
  this->addAction( actionDisplayOrientationXZ);
  connect( actionDisplayOrientationXZ, SIGNAL(triggered()),
	   this, SLOT(Display_orientation_xz()));

  // Connect up the radio buttons to change the orientation
  connect( orientationXZ, SIGNAL( clicked() ), this,
	   SLOT(Display_orientation_xz()) );
  connect( orientationXY, SIGNAL( clicked() ), this,
	   SLOT(Display_orientation_xy()) );
  connect( orientationYZ, SIGNAL( clicked() ), this,
	   SLOT(Display_orientation_yz()) );

  // Connect radio buttons to change LUT
  connect(lutGreyscale, SIGNAL(clicked()), this, SLOT(UseGreyscale()));
  connect(lutProbabilityMap, SIGNAL(clicked()), this, SLOT(UseProbabilityMap()));
  connect(lutSegmentation, SIGNAL(clicked()), this, SLOT(UseSegmentation()));

  connect(blendingSlider,SIGNAL(valueChanged(int)),this,SLOT(Blending_value_changed(int)));
}

void MultiVizViewerWindow::on_actionExit_triggered()
{
  // Call the QWidget close event handler.
  this->close();
}

void MultiVizViewerWindow::closeEvent(QCloseEvent *event)
{
  if (r2d[0])
    {
      for (int i=0; i < 12;i ++)
	{
	  r2d[i]->Delete();
	  r2d[i]=NULL;
	}
    }
}

bool MultiVizViewerWindow::DisplayImage(SceneData *si, ImageTypeDefinitions::DataType t, int number)
{
  if (!si) {
    std::cerr << "SceneImage for handle is a NULL ptr " <<
      si->GetName().toStdString() << std::endl;
    return false;
  }
  if ((t != ImageTypeDefinitions::rgb) && (!dynamic_cast<SceneImage*>(si)->GetImage())) {
    std::cout << "si->getImage() is NULL for handle " << si->GetName().toStdString() <<
      std::endl;
    return false;
  }

  r2d[number]->SetDataType(0,t);
  if (t == ImageTypeDefinitions::probmap)
    lutProbabilityMap->setChecked(true);
  else if (t == ImageTypeDefinitions::segmentation)
    lutSegmentation->setChecked(true);

  // UI settings - only update if they are different
  int oldmin = 0;
  int oldmax = 64;
  int currentSlice = 0;
  currentSlice = r2d[number]->GetSlice();
  r2d[number]->GetSliceRange(oldmin, oldmax);

  if (!r2d[number]->SetScalarInput ( 0, si ) ) return false;

  // Update the UI controls to reflect the new image that has been loaded.
  // If it has a different extent, then reset the slice position if it
  // is not compatible with the existing extent.
  int min; int max; r2d[number]->GetSliceRange(min, max);
  if ((min != oldmin) || (max != oldmax)) {
    sliceScrollBar->setMinimum(min);
    sliceScrollBar->setMaximum(max);
    sliceSpinBox->setMinimum(min);
    sliceSpinBox->setMaximum(max);
    sliceScrollBar->setValue(r2d[number]->GetSlice());
    sliceSpinBox->setValue(r2d[number]->GetSlice());
  }

  return true;
}

void MultiVizViewerWindow::actionSliceDown_triggered()
{
  // Manipulate the current slice display
  if (!r2d[0]) return;
  Display_slice_changed( r2d[0]->GetSlice() - 1 );
}

void MultiVizViewerWindow::actionSliceUp_triggered()
{
  // Manipulate the current slice display
  if (!r2d[0]) return;
  Display_slice_changed( r2d[0]->GetSlice() + 1 );
}

void MultiVizViewerWindow::Display_slice_changed(int s)
{
  // Manipulate the slice that is currently displayed.
  if (!r2d[0]) return;
  for (int i=0;i<12;i++)
    r2d[i]->SetSlice(s);
  sliceScrollBar->setValue(s); // Make sure the scroll bar is up to date.
  sliceSpinBox->setValue(s); // Make sure the scroll bar is up to date.
}

void MultiVizViewerWindow::Display_orientation_xy()
{
  if (!r2d[0]) return;
  for (int i=0;i<12;i++)
    r2d[i]->SetSliceOrientationToXY();
  int min; int max; r2d[0]->GetSliceRange(min, max);
  sliceScrollBar->setMinimum( min );
  sliceScrollBar->setMaximum( max );
  sliceSpinBox->setMinimum(min);
  sliceSpinBox->setMaximum(max);
  sliceScrollBar->setValue( r2d[0]->GetSlice() );
  sliceSpinBox->setValue( r2d[0]->GetSlice() );
}

void MultiVizViewerWindow::Display_orientation_yz()
{
  if (!r2d[0]) return;
  for (int i=0;i<12;i++)
    r2d[i]->SetSliceOrientationToYZ();
  int min; int max; r2d[0]->GetSliceRange(min, max);
  sliceScrollBar->setMinimum( min );
  sliceScrollBar->setMaximum( max );
  sliceSpinBox->setMinimum(min);
  sliceSpinBox->setMaximum(max);
  sliceScrollBar->setValue( r2d[0]->GetSlice() );
  sliceSpinBox->setValue( r2d[0]->GetSlice() );
}

void MultiVizViewerWindow::Display_orientation_xz()
{
  if (!r2d[0]) return;
  for (int i=0;i<12;i++)
    r2d[i]->SetSliceOrientationToXZ();
  int min; int max; r2d[0]->GetSliceRange(min, max);
  sliceScrollBar->setMinimum( min );
  sliceScrollBar->setMaximum( max );
  sliceSpinBox->setMinimum(min);
  sliceSpinBox->setMaximum(max);
  sliceScrollBar->setValue( r2d[0]->GetSlice() );
  sliceSpinBox->setValue( r2d[0]->GetSlice() );
}

MultiVizViewerWindow::~MultiVizViewerWindow()
{
}

void MultiVizViewerWindow::SetMarker(unsigned int _group, double *worldPos)
{
  if (r2d[0])
    for (int i=0;i<12;i++)
      r2d[i]->SetCrosshairPosition(worldPos);
}

void MultiVizViewerWindow::SetMarkerInViewers(unsigned int _group, double *worldPos)
{
  SetMarker(_group,worldPos);
}

void MultiVizViewerWindow::Blending_value_changed(int b)
{
  if (r2d[0])
    for (int i=0;i<12;i++)
      r2d[i]->SetBlendValue((double)(b)/100);
}

void MultiVizViewerWindow::UseGreyscale()
{
  if (r2d[0])
    for (int i=0;i<12;i++)
      r2d[i]->SetDataType(0,ImageTypeDefinitions::greyscale);
}

void MultiVizViewerWindow::UseProbabilityMap()
{
  if (r2d[0])
    for (int i=0;i<12;i++)
      r2d[i]->SetDataType(0,ImageTypeDefinitions::probmap);
}

void MultiVizViewerWindow::UseSegmentation()
{
  if (r2d[0])
    for (int i=0;i<12;i++)
      r2d[i]->SetDataType(0,ImageTypeDefinitions::segmentation);
}

void MultiVizViewerWindow::on_actionShow_Marker_triggered()
{
  showMarker = !showMarker;
  actionShow_Marker->setChecked(showMarker);

  if (r2d[0])
    for (int i=0;i<12;i++)
      r2d[i]->ShowMarker(showMarker);
}

void MultiVizViewerWindow::on_actionBlack_triggered()
{
  if (r2d[0])
    for (int i=0;i<12;i++)
      r2d[i]->SetBackground(0,0,0);
}

void MultiVizViewerWindow::on_actionWhite_triggered()
{
  if (r2d[0])
    for (int i=0;i<12;i++)
      r2d[i]->SetBackground(255,255,255);
}

bool MultiVizViewerWindow::LoadImage(std::string name,ImageTypeDefinitions::DataType t)
{
  if (imageCounter > 11)
    return false;

  SceneImage *si = new SceneImage;
  si->LoadImage(QString(name.c_str()));

  this->DisplayImage(si,t,imageCounter);
  imageCounter++;

  return true;
}

void MultiVizViewerWindow::Render()
{
  if (r2d[0])
    for (int i=0;i<12;i++)
      r2d[i]->Render();
}
