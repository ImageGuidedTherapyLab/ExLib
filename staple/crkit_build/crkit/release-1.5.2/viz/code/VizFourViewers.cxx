#include <VizFourViewers.h>

#include <qsettings.h>
#include <QFileDialog>
#include <QStatusBar>
#include <QKeySequence>
#include <QAction>

// To support drag and drop
#include <QMimeData>
#include <QList>
#include <QUrl>

#include <vtkCallbackCommand.h>

#include "version.h"
#include <VizMainWindow.h>
#include <SceneModel.h>


VizFourViewers::VizFourViewers(VizMainWindow* main)
{
  setupUi(this);

  this->mainWin=main;
  group = 0;
  showMarker = false;
  showBoundingBox = true;
  showOrientationCube = true;
  smallCube = false;

  setAcceptDrops(TRUE);

  // Configure the Render2D3D objects
  r2d[0] = Render2D3D::New();
  r2d[0]->InitializeRendering(qvtkWidget1);
  r2d[1] = Render2D3D::New();
  r2d[1]->InitializeRendering(qvtkWidget2);
  r2d[2] = Render2D3D::New();
  r2d[2]->InitializeRendering(qvtkWidget3);
  r2d[3] = Render2D3D::New();
  r2d[3]->InitializeRendering(qvtkWidget4);
  for (unsigned int n=0; n < 4; n++)
    r2d[n]->SetViewerWindow(this);

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

  sliceScrollBar1->setMinimum(0);
  sliceScrollBar1->setMaximum(255);
  sliceScrollBar1->setValue(255);
  connect( sliceScrollBar1, SIGNAL(valueChanged(int)),
	   this, SLOT(Display_slice_changed1(int))
	   );

  connect( sliceSpinBox1, SIGNAL(valueChanged(int)),
	   this, SLOT(Display_slice_changed1(int))
	   );

  sliceScrollBar2->setMinimum(0);
  sliceScrollBar2->setMaximum(255);
  sliceScrollBar2->setValue(255);
  connect( sliceScrollBar2, SIGNAL(valueChanged(int)),
	   this, SLOT(Display_slice_changed2(int))
	   );

  connect( sliceSpinBox2, SIGNAL(valueChanged(int)),
	   this, SLOT(Display_slice_changed2(int))
	   );

  sliceScrollBar3->setMinimum(0);
  sliceScrollBar3->setMaximum(255);
  sliceScrollBar3->setValue(255);
  connect( sliceScrollBar3, SIGNAL(valueChanged(int)),
	   this, SLOT(Display_slice_changed3(int))
	   );

  connect( sliceSpinBox3, SIGNAL(valueChanged(int)),
	   this, SLOT(Display_slice_changed3(int))
	   );

  sliceScrollBar4->setMinimum(0);
  sliceScrollBar4->setMaximum(255);
  sliceScrollBar4->setValue(255);
  connect( sliceScrollBar4, SIGNAL(valueChanged(int)),
	   this, SLOT(Display_slice_changed4(int))
	   );

  connect( sliceSpinBox4, SIGNAL(valueChanged(int)),
	   this, SLOT(Display_slice_changed4(int))
	   );

  // Connect radio buttons to change LUT
  connect(lutGreyscale, SIGNAL(clicked()), this, SLOT(UseGreyscale()));
  connect(lutProbabilityMap, SIGNAL(clicked()), this, SLOT(UseProbabilityMap()));
  connect(lutSegmentation, SIGNAL(clicked()), this, SLOT(UseSegmentation()));

  connect(blendingSlider,SIGNAL(valueChanged(int)),this,SLOT(Blending_value_changed(int)));
}


VizFourViewers::~VizFourViewers()
{
}


void VizFourViewers::on_actionExit_triggered()
{
  // Call the QWidget close event handler.
  this->close();
}


void VizFourViewers::closeEvent(QCloseEvent *event)
{
  if (r2d[0])
    for (unsigned int v=0;v<4;v++)
      {
	r2d[v]->Delete();
	r2d[v]=NULL;
      }

  mainWin->RemoveViewerFromMap(this);
}


bool VizFourViewers::DisplayImage(SceneImage *si, ImageTypeDefinitions::DataType t)
{
  if (!si) {
    std::cerr << "SceneImage for handle is a NULL ptr " <<
      si->GetName().toStdString() << std::endl;
    return false;
  }
  if (!si->GetImage()) {
    std::cout << "si->getImage() is NULL for handle " << si->GetName().toStdString() <<
      std::endl;
    return false;
  }

  for (unsigned int v=0;v<4;v++)
    r2d[v]->SetDataType(0,t);

  if (t == ImageTypeDefinitions::probmap)
    lutProbabilityMap->setChecked(true);
  else if (t == ImageTypeDefinitions::segmentation)
    lutSegmentation->setChecked(true);

  // UI settings - only update if they are different
  int oldmin = 0;
  int oldmax = 64;
  int currentSlice = 0;
  currentSlice = r2d[0]->GetSlice();
  r2d[0]->GetSliceRange(oldmin, oldmax);

  if (!r2d[0]->SetScalarInput ( 0, si ) ) return false;
  r2d[1]->SetScalarInput(0,si);
  r2d[2]->SetScalarInput(0,si);
  r2d[3]->SetScalarInput(0,si);
  r2d[0]->SetSliceOrientationToXY();
  r2d[1]->SetSliceOrientationToYZ();
  r2d[2]->SetSliceOrientationToXZ();
  r2d[3]->Set3DInteractor();

  UpdateTitle(0,si->GetName());		// should we clear the title if setScalarInput fails?

  // Update the UI controls to reflect the new image that has been loaded.
  // If it has a different extent, then reset the slice position if it
  // is not compatible with the existing extent.
  int min; int max;
  r2d[0]->GetSliceRange(min, max);
  if ((min != oldmin) || (max != oldmax)) {
    sliceScrollBar1->setMinimum(min);
    sliceScrollBar1->setMaximum(max);
    sliceSpinBox1->setMinimum(min);
    sliceSpinBox1->setMaximum(max);
    sliceScrollBar1->setValue(r2d[0]->GetSlice());
    sliceSpinBox1->setValue(r2d[0]->GetSlice());
    r2d[1]->GetSliceRange(min, max);
    sliceScrollBar2->setMinimum(min);
    sliceScrollBar2->setMaximum(max);
    sliceSpinBox2->setMinimum(min);
    sliceSpinBox2->setMaximum(max);
    sliceScrollBar2->setValue(r2d[1]->GetSlice());
    sliceSpinBox2->setValue(r2d[1]->GetSlice());
    r2d[2]->GetSliceRange(min, max);
    sliceScrollBar3->setMinimum(min);
    sliceScrollBar3->setMaximum(max);
    sliceSpinBox3->setMinimum(min);
    sliceSpinBox3->setMaximum(max);
    sliceScrollBar3->setValue(r2d[2]->GetSlice());
    sliceSpinBox3->setValue(r2d[2]->GetSlice());
    r2d[3]->GetSliceRange(min, max);
    sliceScrollBar4->setMinimum(min);
    sliceScrollBar4->setMaximum(max);
    sliceSpinBox4->setMinimum(min);
    sliceSpinBox4->setMaximum(max);
    sliceScrollBar4->setValue(r2d[3]->GetSlice());
    sliceSpinBox4->setValue(r2d[3]->GetSlice());
  }


  // register callbacks
  vtkCallbackCommand *cbc = vtkCallbackCommand::New();
  cbc->SetCallback(this->OnWindowLevelChangedEvent);
  cbc->SetClientData(this);
  for (unsigned int n=0;n<4;n++)
    r2d[n]->AddObserver(Render2D3DWindowLevelEvent,cbc);
  cbc->Delete();

  return true;
}


bool VizFourViewers::DisplayOverlayImage(SceneImage *si, ImageTypeDefinitions::DataType t)
{
  UpdateTitle(1,si->GetName());

  for (unsigned int v=0;v<4;v++)
    {
      r2d[v]->SetDataType(1,t);
      r2d[v]->SetScalarInput(1, si);
    }
  return true;
}


void VizFourViewers::UpdateTitle(unsigned int img, QString s)
{
  // we set this here for convenience
  // img == 0 for main image and img == 1 for overlay
  if ( img == 0 ) {
    this->name = s;
  } else if ( img == 1 ) {
    this->overlay_name = s;
  }

  // either way, we update the title
  if ( this->overlay_name != "" ) {
    this->setWindowTitle("crlViz Viewer - " + this->name + " (" + this->overlay_name + ")" );
  } else {
    this->setWindowTitle("crlViz Viewer - " + this->name );
  }
}


bool VizFourViewers::DisplayTensors(SceneTensors *st)
{
  if (!r2d[0]->SetTensorInput(st)) return false;

  int min,max;
  r2d[0]->GetSliceRange(min,max);
  sliceScrollBar1->setMinimum(min);
  sliceScrollBar1->setMaximum(max);
  sliceSpinBox1->setMinimum(min);
  sliceSpinBox1->setMaximum(max);
  sliceScrollBar1->setValue(0);
  sliceSpinBox1->setValue(0);
  return true;
}


void VizFourViewers::actionSliceDown_triggered()
{
  // Manipulate the current slice display
  if (!r2d[0]) return;
  Display_slice_changed(0, r2d[0]->GetSlice() - 1 );
}


void VizFourViewers::actionSliceUp_triggered()
{
  // Manipulate the current slice display
  if (!r2d[0]) return;
  Display_slice_changed(0, r2d[0]->GetSlice() + 1 );
}


void VizFourViewers::Display_slice_changed(int viewer, int s)
{
  // Manipulate the slice that is currently displayed.
  if (!r2d[viewer]) return;
  r2d[viewer]->SetSlice(s);
  switch (viewer)
    {
    case 0:
      sliceScrollBar1->setValue(s); // Make sure the scroll bar is up to date.
      sliceSpinBox1->setValue(s); // Make sure the scroll bar is up to date.
      break;
    case 1:
      sliceScrollBar2->setValue(s); // Make sure the scroll bar is up to date.
      sliceSpinBox2->setValue(s); // Make sure the scroll bar is up to date.
      break;
    case 2:
      sliceScrollBar3->setValue(s); // Make sure the scroll bar is up to date.
      sliceSpinBox3->setValue(s); // Make sure the scroll bar is up to date.
      break;
    case 3:
      sliceScrollBar4->setValue(s); // Make sure the scroll bar is up to date.
      sliceSpinBox4->setValue(s); // Make sure the scroll bar is up to date.
    }
  //mainWin->PropagateSliceChange(s);
}


void VizFourViewers::Display_slice_changed1(int s)
{
  Display_slice_changed(0,s);
}


void VizFourViewers::Display_slice_changed2(int s)
{
  Display_slice_changed(1,s);
}


void VizFourViewers::Display_slice_changed3(int s)
{
  Display_slice_changed(2,s);
}


void VizFourViewers::Display_slice_changed4(int s)
{
  Display_slice_changed(3,s);
}


void VizFourViewers::dragEnterEvent(QDragEnterEvent *e)
{
  if (!e->mimeData()->hasFormat("application/x-qabstractitemmodeldatalist")) {
    e->ignore();
    return;
  }
  e->acceptProposedAction();
}


void VizFourViewers::dropEvent(QDropEvent *event)
{
  if (!event->mimeData()->hasFormat("application/x-qabstractitemmodeldatalist")) {
    event->ignore();
    return;
  }

  QByteArray data = event->mimeData()->data("application/x-qabstractitemmodeldatalist");
  QDataStream stream(&data, QIODevice::ReadOnly);
  QMap<int, QVariant> map;

  if (!stream.atEnd())
    {
      int r,c;
      stream >> r >> c >> map;
      QString handle = map.begin().value().toString();
      ImageInformation img = mainWin->GetImageInformationForHandle(handle);
      // add greyscale images, probability maps or segmentations as overlays
      if ((img.dataType == ImageTypeDefinitions::greyscale) || (img.dataType == ImageTypeDefinitions::probmap) || (img.dataType == ImageTypeDefinitions::segmentation))
	{
	  std::cout << "Adding overlay " << handle.toStdString() << std::endl;
	  DisplayOverlayImage(dynamic_cast<SceneImage*>(img.sd),img.dataType);
	}
      else if (img.dataType == ImageTypeDefinitions::polydata && r2d[3])
	{
	  r2d[3]->AddModel(dynamic_cast<SceneModel*>(img.sd));
	}
      else if (img.dataType == ImageTypeDefinitions::tensor && r2d[0])
	{
	  std::cout << "Adding tensors " << handle.toStdString() << std::endl;
	  for (unsigned int v=0;v<4;v++)
	    r2d[v]->SetTensorInput(dynamic_cast<SceneTensors*>(img.sd));
	}
    }

  event->acceptProposedAction();
}


void VizFourViewers::SetMarker(unsigned int _group, double *worldPos)
{
  if (r2d[0] && (_group == group))
    {
      for (unsigned int v=0;v<3;v++)
	r2d[v]->SetCrosshairPosition(worldPos);
    }
}

void VizFourViewers::SetMarkerInViewers(unsigned int _group, double *worldPos)
{
  if (mainWin)
    mainWin->SetMarkerInViewers(_group,worldPos);
}


void VizFourViewers::on_actionRemove_Overlay_Image_triggered()
{
  if (r2d[0])
    {
      for (unsigned int v=0;v<4;v++)
	r2d[v]->SetScalarInput(1,NULL);
    }
  this->UpdateTitle(1,"");
}


void VizFourViewers::on_actionRemove_Tensors_triggered()
{
  if (r2d[0])
    for (unsigned int v=0;v<4;v++)
      r2d[v]->SetTensorInput(NULL);
}


void VizFourViewers::Blending_value_changed(int b)
{
  if (r2d[0])
    for (unsigned int v=0;v<4;v++)
      r2d[v]->SetBlendValue((double)(b)/100);
}


void VizFourViewers::UseGreyscale()
{
  if (r2d[0])
    {
      for (unsigned int v=0;v<4;v++)
	r2d[v]->SetDataType(0,ImageTypeDefinitions::greyscale);
    }
}


void VizFourViewers::UseProbabilityMap()
{
  if (r2d[0])
    {
      for (unsigned int v=0;v<4;v++)
	r2d[v]->SetDataType(0,ImageTypeDefinitions::probmap);

    }
}


void VizFourViewers::UseSegmentation()
{
  if (r2d[0])
    {
      for (unsigned int v=0;v<4;v++)
	r2d[v]->SetDataType(0,ImageTypeDefinitions::segmentation);
    }
}


void VizFourViewers::SetGroup(unsigned int _group)
{
  this->group = _group;
}


unsigned int VizFourViewers::GetGroup()
{
  return this->group;
}


void VizFourViewers::on_actionGroup1_triggered()
{
  actionGroup1->setChecked(true);
  actionGroup2->setChecked(false);
  actionGroup3->setChecked(false);
  actionNo_Group->setChecked(false);

  SetGroup(1);
}


void VizFourViewers::on_actionGroup2_triggered()
{
  actionGroup1->setChecked(false);
  actionGroup2->setChecked(true);
  actionGroup3->setChecked(false);
  actionNo_Group->setChecked(false);

  SetGroup(2);
}


void VizFourViewers::on_actionGroup3_triggered()
{
  actionGroup1->setChecked(false);
  actionGroup2->setChecked(false);
  actionGroup3->setChecked(true);
  actionNo_Group->setChecked(false);

  SetGroup(3);
}


void VizFourViewers::on_actionNo_Group_triggered()
{
  actionGroup1->setChecked(false);
  actionGroup2->setChecked(false);
  actionGroup3->setChecked(false);
  actionNo_Group->setChecked(true);

  SetGroup(0);
}


void VizFourViewers::on_actionDisplay_Cursor_triggered()
{
  if (actionDisplay_Cursor->isChecked())
    {
      if (r2d[0]) r2d[0]->SetCursor(true);
    }
  else
    {
      if (r2d[0]) r2d[0]->SetCursor(false);
    }
}


void VizFourViewers::on_actionSave_Points_triggered()
{
  if (r2d[0])
    {
      QString s = QFileDialog::getSaveFileName(this,"Save Points",".","Point Files (*.map)");

      if (s.isEmpty())
	return;
      r2d[0]->SavePoints(s.toStdString());
    }
}


void VizFourViewers::on_actionClear_Points_triggered()
{
  if (r2d[0])
    {
      if (QMessageBox::question(this,"Clear Points","Do you really want to delete all points?",QMessageBox::Yes | QMessageBox::No,QMessageBox::No) == QMessageBox::Yes)
	r2d[0]->ClearPoints();
    }
}


void VizFourViewers::OnWindowLevelChangedEvent(vtkObject *caller,unsigned long event,void *data, void *calldata)
{
  Render2D3D *source = dynamic_cast<Render2D3D*>(caller);
  if (!source)
    // the event does not come from a Render2D3D object, give up
    return;

  VizFourViewers *me = reinterpret_cast<VizFourViewers*>(data);

  double wl[4];
  source->GetWindowLevel(wl);

  for (unsigned int n=0;n<4;n++)
    {
      me->r2d[n]->SetWindowLevel(wl);
      me->r2d[n]->Render();
    }
}


void VizFourViewers::on_actionShow_Marker_triggered()
{
  showMarker = !showMarker;
  actionShow_Marker->setChecked(showMarker);

  if (r2d[0])
    {
      for (unsigned int v=0;v<4;v++)
	r2d[v]->ShowMarker(showMarker);
    }
}


void VizFourViewers::on_actionShow_Bounding_Box_triggered()
{
  showBoundingBox = !showBoundingBox;
  actionShow_Bounding_Box->setChecked(showBoundingBox);

  if (r2d[3])
    r2d[3]->ShowBoundingBox(showBoundingBox);
}


void VizFourViewers::on_actionShow_Orientation_Cube_triggered()
{
  showOrientationCube = !showOrientationCube;
  actionShow_Orientation_Cube->setChecked(showOrientationCube);

  if (r2d[3] && showOrientationCube)
    r2d[3]->ShowOrientationCube(1.0);
  else if (r2d[3] && !showOrientationCube)
    r2d[3]->ShowOrientationCube(0.0);
}


void VizFourViewers::on_actionSmall_Cube_triggered()
{
  smallCube = !smallCube;
  actionSmall_Cube->setChecked(smallCube);

  if (r2d[3] && smallCube)
    r2d[3]->ShowOrientationCube(0.5);
  else if (r2d[3] && !smallCube)
    r2d[3]->ShowOrientationCube(1.0);
}


void VizFourViewers::on_actionBlack_triggered()
{
  if (r2d[0])
    {
      for (unsigned int v=0;v<4;v++)
	r2d[v]->SetBackground(0,0,0);
    }
}


void VizFourViewers::on_actionWhite_triggered()
{
  if (r2d[0])
    {
      for (unsigned int v=0;v<4;v++)
	r2d[v]->SetBackground(255,255,255);
    }
}


unsigned int VizFourViewers::GetNumberOfRenderersForObject(SceneData *sd)
{
  if (r2d[0])
    return r2d[0]->GetNumberOfRenderersForObject(sd)+r2d[1]->GetNumberOfRenderersForObject(sd)+r2d[2]->GetNumberOfRenderersForObject(sd)+r2d[3]->GetNumberOfRenderersForObject(sd);
  else
    return 0;
}


void VizFourViewers::Render()
{
  if (r2d[0])
    for (unsigned int v=0;v<4;v++)
      r2d[v]->Render();
}
