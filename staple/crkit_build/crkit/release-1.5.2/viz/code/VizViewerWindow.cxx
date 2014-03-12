/* This maintains the state of a viewer window. */

#include <VizViewerWindow.h>

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
#include <SceneModel.h>
#include <VizViewerWindowRotationDialog.h>

VizViewerWindow::VizViewerWindow(VizMainWindow* main, bool stereo)
{
  setupUi(this);

  this->mainWin=main;
  group = 1;
  showMarker = false;
  showBoundingBox = true;
  showOrientationCube = true;
  smallCube = false;
  showMousePointer = true;
  name = "";
  overlay_name = "";

  setAcceptDrops(TRUE);

  // Configure the Render2D3D object
  r2d = Render2D3D::New();
  r2d->SetViewerWindow(this);
  r2d->stereo = stereo;
  r2d->InitializeRendering(qvtkWidget);

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

  actionFullScreen = new QAction(this);
  actionFullScreen->setShortcutContext(Qt::WindowShortcut);
  actionFullScreen->setShortcut(QKeySequence(Qt::Key_F11));
  this->addAction(actionFullScreen);
  connect(actionFullScreen,SIGNAL(triggered()),this,SLOT(FullScreen()));

  // Connect up the radio buttons to change the orientation
  connect( orientationXZ, SIGNAL( clicked() ), this,
	   SLOT(Display_orientation_xz()) );
  connect( orientationXY, SIGNAL( clicked() ), this,
	   SLOT(Display_orientation_xy()) );
  connect( orientationYZ, SIGNAL( clicked() ), this,
	   SLOT(Display_orientation_yz()) );

  // Connect the radio buttons to change the render window interactor
  connect(interactor2D, SIGNAL(clicked()), this, SLOT(Use2DInteractor()));
  connect(interactor3D, SIGNAL(clicked()), this, SLOT(Use3DInteractor()));

  // Connect radio buttons to change LUT
  connect(lutGreyscale, SIGNAL(clicked()), this, SLOT(UseGreyscale()));
  connect(lutProbabilityMap, SIGNAL(clicked()), this, SLOT(UseProbabilityMap()));
  connect(lutSegmentation, SIGNAL(clicked()), this, SLOT(UseSegmentation()));

  connect(blendingSlider,SIGNAL(valueChanged(int)),this,SLOT(Blending_value_changed(int)));

#ifndef HAVE_VVLINK
  menuIGT->setEnabled(false);
#endif
}

void VizViewerWindow::on_actionExit_triggered()
{
  // Call the QWidget close event handler.
  this->close();
}

void VizViewerWindow::closeEvent(QCloseEvent *event)
{
  if (r2d)
    {
      r2d->Delete();
      r2d=NULL;
    }
  mainWin->RemoveViewerFromMap(this);
}

bool VizViewerWindow::DisplayImage(SceneData *si, ImageTypeDefinitions::DataType t)
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

  r2d->SetDataType(0,t);
  if (t == ImageTypeDefinitions::probmap)
    lutProbabilityMap->setChecked(true);
  else if (t == ImageTypeDefinitions::segmentation)
    lutSegmentation->setChecked(true);

  // UI settings - only update if they are different
  int oldmin = 0;
  int oldmax = 64;
  int currentSlice = 0;
  currentSlice = r2d->GetSlice();
  r2d->GetSliceRange(oldmin, oldmax);

  if (!r2d->SetScalarInput ( 0, si ) ) return false;

  UpdateTitle(0,si->GetName());		// should we clear the title if setScalarInput fails?

  // Update the UI controls to reflect the new image that has been loaded.
  // If it has a different extent, then reset the slice position if it
  // is not compatible with the existing extent.
  int min; int max; r2d->GetSliceRange(min, max);
  //if ((min != oldmin) || (max != oldmax)) {
  sliceScrollBar->setMinimum(min);
  sliceScrollBar->setMaximum(max);
  sliceSpinBox->setMinimum(min);
  sliceSpinBox->setMaximum(max);
  sliceScrollBar->setValue(r2d->GetSlice());
  sliceSpinBox->setValue(r2d->GetSlice());
  //}

  return true;
}


bool VizViewerWindow::DisplayOverlayImage(SceneData *si, ImageTypeDefinitions::DataType t)
{
  r2d->SetDataType(1,t);
  bool result = r2d->SetScalarInput(1, si);
  if (result)
    {
      Display_slice_changed(r2d->GetSlice());
      UpdateTitle(1,si->GetName());
    }

  return result;
}

void VizViewerWindow::UpdateTitle(unsigned int img, QString s)
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


bool VizViewerWindow::DisplayTensors(SceneTensors *st)
{
  if (!r2d->SetTensorInput(st)) return false;

  int min,max;
  r2d->GetSliceRange(min,max);
  sliceScrollBar->setMinimum(min);
  sliceScrollBar->setMaximum(max);
  sliceSpinBox->setMinimum(min);
  sliceSpinBox->setMaximum(max);
  sliceScrollBar->setValue(0);
  sliceSpinBox->setValue(0);
  return true;
}


bool VizViewerWindow::DisplayModelOnly(SceneModel *sm)
{
  if (!r2d->AddModel(sm))
    return false;

  r2d->ResetCamera();
  return true;
}


void VizViewerWindow::actionSliceDown_triggered()
{
  // Manipulate the current slice display
  if (!r2d) return;
  Display_slice_changed( r2d->GetSlice() - 1 );
}

void VizViewerWindow::actionSliceUp_triggered()
{
  // Manipulate the current slice display
  if (!r2d) return;
  Display_slice_changed( r2d->GetSlice() + 1 );
}

void VizViewerWindow::Display_slice_changed(int s)
{
  // Manipulate the slice that is currently displayed.
  if (!r2d) return;
  r2d->SetSlice(s);
  sliceScrollBar->setValue(s); // Make sure the scroll bar is up to date.
  sliceSpinBox->setValue(s); // Make sure the scroll bar is up to date.
  mainWin->PropagateSliceChange(this->group,s);
}

void VizViewerWindow::Display_orientation_xy()
{
  if (!r2d) return;
  r2d->SetSliceOrientationToXY();
  int min; int max; r2d->GetSliceRange(min, max);
  sliceScrollBar->setMinimum( min );
  sliceScrollBar->setMaximum( max );
  sliceSpinBox->setMinimum(min);
  sliceSpinBox->setMaximum(max);
  sliceScrollBar->setValue( r2d->GetSlice() );
  sliceSpinBox->setValue( r2d->GetSlice() );
  mainWin->PropagateOrientationChange(this->group,xy);
}

void VizViewerWindow::Display_orientation_yz()
{
  if (!r2d) return;
  r2d->SetSliceOrientationToYZ();
  int min; int max; r2d->GetSliceRange(min, max);
  sliceScrollBar->setMinimum( min );
  sliceScrollBar->setMaximum( max );
  sliceSpinBox->setMinimum(min);
  sliceSpinBox->setMaximum(max);
  sliceScrollBar->setValue( r2d->GetSlice() );
  sliceSpinBox->setValue( r2d->GetSlice() );
  mainWin->PropagateOrientationChange(this->group,yz);
}

void VizViewerWindow::Display_orientation_xz()
{
  if (!r2d) return;
  r2d->SetSliceOrientationToXZ();
  int min; int max; r2d->GetSliceRange(min, max);
  sliceScrollBar->setMinimum( min );
  sliceScrollBar->setMaximum( max );
  sliceSpinBox->setMinimum(min);
  sliceSpinBox->setMaximum(max);
  sliceScrollBar->setValue( r2d->GetSlice() );
  sliceSpinBox->setValue( r2d->GetSlice() );
  mainWin->PropagateOrientationChange(this->group,xz);
}

VizViewerWindow::~VizViewerWindow()
{
}

void VizViewerWindow::dragEnterEvent(QDragEnterEvent *e)
{
  if (!e->mimeData()->hasFormat("application/x-qabstractitemmodeldatalist")) {
    e->ignore();
    return;
  }
  e->acceptProposedAction();
}

void VizViewerWindow::dropEvent(QDropEvent *event)
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
      if ((img.dataType == ImageTypeDefinitions::greyscale) || (img.dataType == ImageTypeDefinitions::probmap) || (img.dataType == ImageTypeDefinitions::segmentation) || (img.dataType == ImageTypeDefinitions::fmri) || (img.dataType == ImageTypeDefinitions::rgb))
	{
	  //std::cout << "Adding overlay " << handle.toStdString() << std::endl;
	  DisplayOverlayImage(img.sd,img.dataType);
	}
      else if ((img.dataType == ImageTypeDefinitions::polydata || img.dataType == ImageTypeDefinitions::scalarpolydata) && r2d)
	{
	  r2d->AddModel(dynamic_cast<SceneModel*>(img.sd));
	}
      else if (img.dataType == ImageTypeDefinitions::tensor && r2d)
	{
	  //std::cout << "Adding tensors " << handle.toStdString() << std::endl;
	  r2d->SetTensorInput(dynamic_cast<SceneTensors*>(img.sd));
	}
    }

  event->acceptProposedAction();
}

void VizViewerWindow::SetMarker(unsigned int _group, double *worldPos)
{
  if (r2d && (_group == group))
    r2d->SetCrosshairPosition(worldPos);
}


void VizViewerWindow::SetMarkerInViewers(unsigned int _group, double *worldPos)
{
  if (mainWin)
    mainWin->SetMarkerInViewers(_group,worldPos);
}

void VizViewerWindow::CursorPositionHasChanged(double pos[3])
{
  mainWin->PropagateCursorPosition(group,pos);
}

void VizViewerWindow::on_actionRemove_Overlay_Image_triggered()
{
  if (r2d) {
    r2d->SetScalarInput(1,NULL);
  }
  this->UpdateTitle(1,"");
}

void VizViewerWindow::on_actionRemove_Tensors_triggered()
{
  if (r2d)
    r2d->SetTensorInput(NULL);
}

void VizViewerWindow::Blending_value_changed(int b)
{
  if (r2d)
    r2d->SetBlendValue((double)(b)/100);
}

void VizViewerWindow::Use2DInteractor()
{
  if (r2d)
    r2d->Set2DInteractor();

  actionShow_Bounding_Box->setEnabled(false);
  actionShow_Orientation_Cube->setEnabled(false);
  actionSmall_Cube->setEnabled(false);
}

void VizViewerWindow::Use3DInteractor()
{
  if (r2d)
    r2d->Set3DInteractor();

  actionShow_Bounding_Box->setEnabled(true);
  actionShow_Bounding_Box->setChecked(true);
  actionShow_Orientation_Cube->setEnabled(true);
  actionShow_Orientation_Cube->setChecked(true);
  actionSmall_Cube->setEnabled(true);
  actionSmall_Cube->setChecked(false);
}

void VizViewerWindow::UseGreyscale()
{
  if (r2d)
    r2d->SetDataType(0,ImageTypeDefinitions::greyscale);
}

void VizViewerWindow::UseProbabilityMap()
{
  if (r2d)
    r2d->SetDataType(0,ImageTypeDefinitions::probmap);
}

void VizViewerWindow::UseSegmentation()
{
  if (r2d)
    r2d->SetDataType(0,ImageTypeDefinitions::segmentation);
}

void VizViewerWindow::SetGroup(unsigned int _group)
{
  this->group = _group;
}

unsigned int VizViewerWindow::GetGroup()
{
  return this->group;
}

void VizViewerWindow::on_actionGroup1_triggered()
{
  actionGroup1->setChecked(true);
  actionGroup2->setChecked(false);
  actionGroup3->setChecked(false);
  actionNo_Group->setChecked(false);

  SetGroup(1);
}

void VizViewerWindow::on_actionGroup2_triggered()
{
  actionGroup1->setChecked(false);
  actionGroup2->setChecked(true);
  actionGroup3->setChecked(false);
  actionNo_Group->setChecked(false);

  SetGroup(2);
}

void VizViewerWindow::on_actionGroup3_triggered()
{
  actionGroup1->setChecked(false);
  actionGroup2->setChecked(false);
  actionGroup3->setChecked(true);
  actionNo_Group->setChecked(false);

  SetGroup(3);
}

void VizViewerWindow::on_actionNo_Group_triggered()
{
  actionGroup1->setChecked(false);
  actionGroup2->setChecked(false);
  actionGroup3->setChecked(false);
  actionNo_Group->setChecked(true);

  SetGroup(0);
}

void VizViewerWindow::on_actionDisplay_Cursor_triggered()
{
  if (actionDisplay_Cursor->isChecked())
    {
      if (r2d) r2d->SetCursor(true);
    }
  else
    {
      if (r2d) r2d->SetCursor(false);
    }
}

void VizViewerWindow::on_actionSave_Points_triggered()
{
  if (r2d)
    {
      QString s = QFileDialog::getSaveFileName(this,"Save Points",".","Point Files (*.map)");

      if (s.isEmpty())
	return;
      r2d->SavePoints(s.toStdString());
    }
}

void VizViewerWindow::on_actionLoad_Points_triggered()
{
  if (r2d)
    {
      QString s = QFileDialog::getOpenFileName(this,"Load Points",".","Point Files (*.map)");
      if (s.isEmpty())
	return;
      r2d->LoadPoints(s.toStdString());
    }
}

void VizViewerWindow::on_actionClear_Points_triggered()
{
  if (r2d)
    {
      if (QMessageBox::question(this,"Clear Points","Do you really want to delete all points?",QMessageBox::Yes | QMessageBox::No,QMessageBox::No) == QMessageBox::Yes)
	r2d->ClearPoints();
      r2d->Render();
    }
}

void VizViewerWindow::on_actionSelect_ROI_triggered()
{
  if (actionSelect_ROI->isChecked())
    {
      if (r2d) r2d->EnableROIBox();
    }
  else
    {
      if (r2d) r2d->DisableROIBox();
    }
}

void VizViewerWindow::on_actionSave_ROI_triggered()
{
  if (r2d && actionSelect_ROI->isChecked())
    {
      QString s = QFileDialog::getSaveFileName(this,"Save ROI",".","Volumes (*.nhdr *.nrrd *.mhd *.mha *.vtk *.hdr *.crl)");

      if (s.isEmpty())
	return;
      r2d->SaveROI(s.toStdString());
    }
}

void VizViewerWindow::on_actionSave_Screenshot_triggered()
{
  if (r2d)
    {
      QString s = QFileDialog::getSaveFileName(this,"Save Screenshot",".","PNG Images (*.png)");

      if (s.isEmpty())
	return;
      if (!(s.contains(".png",Qt::CaseInsensitive)))
	s = s+".png";
      // sleep for a second to make dialog box disappear before saving screenshot
#ifndef WIN32
      sleep(1);
#else
      Sleep(1000);
#endif
      r2d->SaveScreenshot(s.toStdString());
    }
}

void VizViewerWindow::on_actionSave_Movie_triggered()
{
  if (r2d)
    {
      QString s = QFileDialog::getSaveFileName(this,"Save Movie (file prefix)",".","Prefix (*.*)");

      if (s.isEmpty())
	return;
      r2d->SaveMovie(s.toStdString());
    }
}

void VizViewerWindow::on_actionSave_OBJ_File_triggered()
{
  if (r2d)
    {
      QString s = QFileDialog::getSaveFileName(this,"Save OBJ File",".","OBJ Files (*.obj)");

      if (s.isEmpty())
	return;
      r2d->SaveOBJ(s.toStdString());
    }
}

void VizViewerWindow::on_actionSave_RIB_File_triggered()
{
  if (r2d)
    {
      QString s = QFileDialog::getSaveFileName(this,"Save RIB File",".","RIB Files (*.rib)");

      if (s.isEmpty())
	return;
      r2d->SaveRIB(s.toStdString());
    }
}

void VizViewerWindow::on_actionShow_Marker_triggered()
{
  showMarker = !showMarker;
  actionShow_Marker->setChecked(showMarker);

  if (r2d)
    r2d->ShowMarker(showMarker);
}

void VizViewerWindow::on_actionShow_Bounding_Box_triggered()
{
  showBoundingBox = !showBoundingBox;
  actionShow_Bounding_Box->setChecked(showBoundingBox);

  if (r2d)
    r2d->ShowBoundingBox(showBoundingBox);
}

void VizViewerWindow::on_actionShow_Orientation_Cube_triggered()
{
  showOrientationCube = !showOrientationCube;
  actionShow_Orientation_Cube->setChecked(showOrientationCube);

  if (r2d && showOrientationCube)
    r2d->ShowOrientationCube(1.0);
  else if (r2d && !showOrientationCube)
    r2d->ShowOrientationCube(0.0);
}

void VizViewerWindow::on_actionSmall_Cube_triggered()
{
  smallCube = !smallCube;
  actionSmall_Cube->setChecked(smallCube);

  if (r2d && smallCube)
    r2d->ShowOrientationCube(0.5);
  else if (r2d && !smallCube)
    r2d->ShowOrientationCube(1.0);
}

void VizViewerWindow::on_actionBlack_triggered()
{
  if (r2d)
    r2d->SetBackground(0,0,0);
}

void VizViewerWindow::on_actionWhite_triggered()
{
  if (r2d)
    r2d->SetBackground(255,255,255);
}

void VizViewerWindow::FullScreen()
{
  if (r2d)
    r2d->ToggleFullScreen();
}

unsigned int VizViewerWindow::GetNumberOfRenderersForObject(SceneData *sd)
{
  if (r2d)
    return r2d->GetNumberOfRenderersForObject(sd);
  else
    return 0;
}

void VizViewerWindow::SaveModel(SceneData *model,std::string fileName)
{
  if (r2d)
    r2d->SaveModel(model,fileName);
}

void VizViewerWindow::on_actionShow_Mouse_Pointer_triggered()
{
  showMousePointer = !showMousePointer;
  actionShow_Mouse_Pointer->setChecked(showMousePointer);

  if (r2d)
    r2d->ShowMousePointer(showMousePointer);
}

void VizViewerWindow::RotateModel(SceneData *model)
{
  if (r2d)
    {
      VizViewerWindowRotationDialog *rdlg = new VizViewerWindowRotationDialog();
      rdlg->SetRenderer(r2d);
      rdlg->SetModel(dynamic_cast<SceneModel*>(model));
      rdlg->show();
    }
}

void VizViewerWindow::on_actionSpecial_Label_Blending_triggered()
{
  if (r2d)
    r2d->SetSpecialLabelBlending(actionSpecial_Label_Blending->isChecked());
}


void VizViewerWindow::Render()
{
  if (r2d)
    r2d->Render();
}

void VizViewerWindow::on_actionShow_Tool_triggered()
{
#ifdef HAVE_VVLINK
  if (r2d)
    {
      r2d->ShowIGTTool(actionShow_Tool->isChecked());
      if (actionShow_Tool->isChecked())
	mainWin->TrackToolOn(r2d);
      else
	mainWin->TrackToolOff(r2d);
    }
#endif
}

void VizViewerWindow::on_actionDisplay_in_VectorVision_triggered()
{
#ifdef HAVE_VVLINK
  if (!r2d) return;

  if (actionDisplay_in_VectorVision->isChecked())
    {
      if (!mainWin->DisplayInVectorVision(r2d))
	// this fails if another render window is already displayed in VectorVision, so remove the check mark here
	actionDisplay_in_VectorVision->setChecked(false);
    }
  else
    mainWin->DisplayInVectorVision(NULL);
#endif
}
