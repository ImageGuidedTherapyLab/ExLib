
/* This maintains the state of the main window. */

#include <VizMainWindow.h>
#include <QSettings>
#include <QFileDialog>
#include <QStatusBar>
#include <QKeySequence>
#include <QAction>
#include <QFileInfo>

// To support drag and drop
#include <QMimeData>
#include <QList>
#include <QUrl>

// To support XML reading and writing
#include <QXmlStreamReader>

#include "version.h"
#include <VizViewerWindow.h>
#include <SceneImage.h>
#include <SceneModel.h>
#include <SceneTensors.h>
#include <SceneRGBImage.h>
#include <SceneImageView.h>
#include <VizFourViewers.h>
#include <VizWizardSTAPLE.h>

#include <set>


#ifdef HAVE_VVLINK
#include <vtkWindowToImageFilter.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkMatrix4x4.h>
#include <vtkImageCast.h>
#include <vtkMath.h>
#include <VVLQGetAvailableInstruments.h>
#include <VVLQGetSliceSets.h>
#include <VVLSliceSet.h>
#endif

VizMainWindow::VizMainWindow() : stereo(false),inSliceRecursion(false),inOrientationRecursion(false),inCursorRecursion(false)
{
  setupUi(this);
  
  // Configure application specific settings.
  QCoreApplication::setOrganizationName("crl");
  QCoreApplication::setOrganizationDomain("crl.med.harvard.edu");
  QCoreApplication::setApplicationName("crlViz");
  
  appsettings = new QSettings();
  // Restore the main window position
  appsettings->beginGroup("MainWindow");
  resize(appsettings->value("size", QSize(330, 600)).toSize());
  move(appsettings->value("pos", QSize(200, 200)).toPoint());
  appsettings->endGroup();
  // Later we should define a data structure to hold user preferences
  // and a menu to allow them to be set. They could be stored in the
  // appsettings file each time they are modified.
  appsettings->beginGroup("UserPreferences");
  bool resetInitialPathOnStartup =
    appsettings->value("resetInitialPath", true).toBool();
  appsettings->endGroup();
  // Reset the initial path if the user desires.
  if (resetInitialPathOnStartup) {
    appsettings->setValue("RecentPath/dir", ".");
  }
  
  setAcceptDrops(TRUE);
  
  fileListWidget->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(fileListWidget,SIGNAL(itemDoubleClicked(QListWidgetItem*)),this,SLOT(FileListDoubleClick(QListWidgetItem*)));
  connect(fileListWidget,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(FileListRightClick(const QPoint&)));
  
#ifdef HAVE_VVLINK
  server = NULL;
  connected = false;
  coord = 0;
  trackedToolNotifier = NULL;
  view = vvlink::VVLQSetView::New();
  view->SetChannel("crlViz Research Data");
  tool = vvlink::VVLInstrument::New();
  displayInVectorVision = NULL;
  inTimer = false;
  wtif = vtkWindowToImageFilter::New();
  wtif->SetInputBufferTypeToRGB();
  lastToolPosition = vtkMatrix4x4::New();
  lastToolPosition->Identity();
  QTimer *timer = new QTimer(this);
  connect(timer,SIGNAL(timeout()),this,SLOT(timerCallback()));
  timer->start(1000);
  connect(menuTools,SIGNAL(triggered(QAction*)),this,SLOT(toolsMenuCallback(QAction*)));
  menuTools->addAction("not connected")->setEnabled(false);
#else
  menuIGT->setEnabled(false);
#endif
}

void VizMainWindow::on_actionExit_triggered()
{
  // Call the QWidget close event handler.
  this->close();
}

void VizMainWindow::on_actionAbout_Viz_triggered()
{
  QString messageText = QString("crlViz version ") +
    QString(Version::GetVersion()) +
    QString( "\n" ) +
    QString( "Visit http://www.crl.med.harvard.edu for more information.");
  QMessageBox::about(this, "crlViz", messageText);
}

void VizMainWindow::closeEvent(QCloseEvent *event)
{
  // Save the application settings when we close.
  appsettings->beginGroup("MainWindow");
  appsettings->setValue("size", size());
  appsettings->setValue("pos", pos());
  appsettings->endGroup();
  appsettings->sync();
  
  ImageInformationMap::iterator i = imagesInScene.begin();
  while (i != imagesInScene.end())
    {
      delete i.value().sd;
      
      while (!i.value().viewers.empty())
	{
	  (*i.value().viewers.begin())->close();
	}
      i++;
    }
}

void VizMainWindow::on_actionOpen_Scene_triggered()
{
  QSettings appsettings;
  appsettings.beginGroup("RecentPath");
  QString recentDirName = appsettings.value("dir",QString(".")).toString();
  if (ignoreLastDirectory) recentDirName = ".";
  
  QString s = QFileDialog::getOpenFileName(this,"Open Scene",recentDirName,"crlViz Scenes (*.xml)");
  
  if (s.isEmpty())
    return;
  
  if (!QFile::exists(s))
    {
      this->statusBar()->showMessage("File " + s + " doesn't exist.",9000);
      return;
    }
  
  QFileInfo info = QFileInfo(s);
  appsettings.setValue("dir",info.absolutePath());
  appsettings.sync();
  
  OpenScene(s);
}

void VizMainWindow::on_actionSave_Scene_triggered()
{
  QString recentDirName = ".";
  
  QString s = QFileDialog::getSaveFileName(this,"Save Scene",recentDirName,"crlViz Scenes (*.xml)");
  
  if (s.isEmpty())
    return;
  
  SaveScene(s);
}

void VizMainWindow::on_actionAdd_Volume_triggered()
{
  //std::cout << "Add Volume triggered." << std::endl;
  
  QSettings appsettings;
  appsettings.beginGroup("RecentPath");
  QString recentDirName = appsettings.value("dir",QString(".")).toString();
  if (ignoreLastDirectory) recentDirName = ".";
  
  QString s = QFileDialog::getOpenFileName(
					   this,
					   "Add Volume",
					   recentDirName,
					   "Volumes (*.nhdr *.nrrd *.nii *.mhd *.mha *.vtk *.hdr *.crl)"
					   );
  if ( s.isEmpty() ) {
    return;
  }
  if (!QFile::exists(s)) {
    this->statusBar()->showMessage( "File " + s + " doesn't exist.", 9000 );
    return;
  }
  
  QFileInfo info = QFileInfo(s);
  appsettings.setValue("dir",info.absolutePath());
  appsettings.sync();
  
  FileOpenWithName(s,ImageTypeDefinitions::greyscale);
}

void VizMainWindow::on_actionAdd_Probability_Map_triggered()
{
  QSettings appsettings;
  appsettings.beginGroup("RecentPath");
  QString recentDirName = appsettings.value("dir",QString(".")).toString();
  if (ignoreLastDirectory) recentDirName = ".";
  
  QString s = QFileDialog::getOpenFileName(
					   this,
					   "Add Probability Map",
					   recentDirName,
					   "Volumes (*.nhdr *.nrrd *.nii *.mhd *.mha *.vtk *.hdr *.crl)"
					   );
  if ( s.isEmpty() ) {
    return;
  }
  if (!QFile::exists(s)) {
    this->statusBar()->showMessage( "File " + s + " doesn't exist.", 9000 );
    return;
  }
  
  QFileInfo info = QFileInfo(s);
  appsettings.setValue("dir",info.absolutePath());
  appsettings.sync();
  
  FileOpenWithName(s,ImageTypeDefinitions::probmap);
}

void VizMainWindow::on_actionAdd_Segmented_Volume_triggered()
{
  QSettings appsettings;
  appsettings.beginGroup("RecentPath");
  QString recentDirName = appsettings.value("dir",QString(".")).toString();
  if (ignoreLastDirectory) recentDirName = ".";
  
  QString s = QFileDialog::getOpenFileName(
					   this,
					   "Add Segmented Volume",
					   recentDirName,
					   "Volumes (*.nhdr *.nrrd *.nii *.mhd *.mha *.vtk *.hdr *.crl)"
					   );
  if ( s.isEmpty() ) {
    return;
  }
  if (!QFile::exists(s)) {
    this->statusBar()->showMessage( "File " + s + " doesn't exist.", 9000 );
    return;
  }
  
  QFileInfo info = QFileInfo(s);
  appsettings.setValue("dir",info.absolutePath());
  appsettings.sync();
  
  FileOpenWithName(s,ImageTypeDefinitions::segmentation);
}

void VizMainWindow::on_actionAdd_fMRI_Activation_Map_triggered()
{
  QSettings appsettings;
  appsettings.beginGroup("RecentPath");
  QString recentDirName = appsettings.value("dir",QString(".")).toString();
  if (ignoreLastDirectory) recentDirName = ".";
  
  QString s = QFileDialog::getOpenFileName(
					   this,
					   "Add fMRI Activation Map",
					   recentDirName,
					   "Volumes (*.nhdr *.nrrd *.nii *.mhd *.mha *.vtk *.hdr *.crl)"
					   );
  if ( s.isEmpty() ) {
    return;
  }
  if (!QFile::exists(s)) {
    this->statusBar()->showMessage( "File " + s + " doesn't exist.", 9000 );
    return;
  }
  
  QFileInfo info = QFileInfo(s);
  appsettings.setValue("dir",info.absolutePath());
  appsettings.sync();
  
  FileOpenWithName(s,ImageTypeDefinitions::fmri);
}

void VizMainWindow::on_actionAdd_Tensor_Volume_triggered()
{
  QSettings appsettings;
  appsettings.beginGroup("RecentPath");
  QString recentDirName = appsettings.value("dir",QString(".")).toString();
  if (ignoreLastDirectory) recentDirName = ".";
  
  QString s = QFileDialog::getOpenFileName(
					   this,
					   "Add Tensor Volume",
					   recentDirName,
					   "Volumes (*.nhdr *.nrrd *.nii *.mhd *.mha *.vtk *.hdr)"
					   );
  if ( s.isEmpty() ) {
    return;
  }
  if (!QFile::exists(s)) {
    this->statusBar()->showMessage( "File " + s + " doesn't exist.", 9000 );
    return;
  }
  
  QFileInfo info = QFileInfo(s);
  appsettings.setValue("dir",info.absolutePath());
  appsettings.sync();
  
  FileOpenWithName(s,ImageTypeDefinitions::tensor);
}

void VizMainWindow::on_actionAdd_RGB_Volume_triggered()
{
  QSettings appsettings;
  appsettings.beginGroup("RecentPath");
  QString recentDirName = appsettings.value("dir",QString(".")).toString();
  if (ignoreLastDirectory) recentDirName = ".";
  
  QString s = QFileDialog::getOpenFileName(this,"Add RGB Tensor Volume",recentDirName,"Volumes (*.nhdr *.nrrd *.nii *.mhd *.mha *.vtk *.hdr)");
  
  if (s.isEmpty())
    {
      return;
    }
  if (!QFile::exists(s))
    {
      this->statusBar()->showMessage("File "+s+" doesn't exist.",9000);
      return;
    }
  
  QFileInfo info = QFileInfo(s);
  appsettings.setValue("dir",info.absolutePath());
  appsettings.sync();
  
  FileOpenWithName(s,ImageTypeDefinitions::rgb);
}

void VizMainWindow::on_actionAdd_Model_triggered()
{
  QSettings appsettings;
  appsettings.beginGroup("RecentPath");
  QString recentDirName = appsettings.value("dir",QString(".")).toString();
  if (ignoreLastDirectory) recentDirName = ".";
  
  QString s = QFileDialog::getOpenFileName(
					   this,
					   "Add Model",
					   recentDirName,
					   "Models (*.vtk *.vtp *.gii)"
					   );
  if ( s.isEmpty() ) {
    return;
  }
  if (!QFile::exists(s)) {
    this->statusBar()->showMessage( "File " + s + " doesn't exist.", 9000 );
    return;
  }
  
  QFileInfo info = QFileInfo(s);
  appsettings.setValue("dir",info.absolutePath());
  appsettings.sync();
  
  FileOpenWithName(s,ImageTypeDefinitions::polydata);
}

void VizMainWindow::on_actionAdd_Model_with_Scalars_triggered()
{
  QSettings appsettings;
  appsettings.beginGroup("RecentPath");
  QString recentDirName = appsettings.value("dir",QString(".")).toString();
  if (ignoreLastDirectory) recentDirName = ".";
  
  QString s = QFileDialog::getOpenFileName(
					   this,
					   "Add Model",
					   recentDirName,
					   "Models (*.vtk *.vtp *.gii)"
					   );
  if ( s.isEmpty() ) {
    return;
  }
  if (!QFile::exists(s)) {
    this->statusBar()->showMessage( "File " + s + " doesn't exist.", 9000 );
    return;
  }
  
  QFileInfo info = QFileInfo(s);
  appsettings.setValue("dir",info.absolutePath());
  appsettings.sync();
  
  FileOpenWithName(s,ImageTypeDefinitions::scalarpolydata);
}

VizMainWindow::~VizMainWindow()
{
}

void VizMainWindow::dragEnterEvent(QDragEnterEvent *e)
{
  if (e->mimeData()->hasUrls() || e->mimeData()->hasFormat("application/x-qabstractitemmodeldatalist")) {
    e->acceptProposedAction();
    return;
  }
  e->ignore();
}

void VizMainWindow::dropEvent(QDropEvent *event)
{
  if (event->mimeData()->hasUrls())
    {
      QList<QUrl> urls = event->mimeData()->urls();
      for (int i = 0; i < urls.size(); i++) {
	QUrl url = urls.at(i);
	if (url.toString().isEmpty()) continue;
	QString fileName = QDir::convertSeparators( url.toLocalFile() );
	
	if (!fileName.isEmpty()) {
	  if (!FileOpenWithName(fileName,ImageTypeDefinitions::greyscale)) {
	    QMessageBox::critical(this,"crlViz","File could not be opened.");
	  }
	}
      }
    }
  else if (event->mimeData()->hasFormat("application/x-qabstractitemmodeldatalist"))
    {
      QByteArray data = event->mimeData()->data("application/x-qabstractitemmodeldatalist");
      QDataStream stream(&data, QIODevice::ReadOnly);
      QMap<int, QVariant> map;
      
      if (!stream.atEnd())
	{
	  int r,c;
	  stream >> r >> c >> map;
	  QString handle = map.begin().value().toString();
	  ImageInformation img = imagesInScene[handle];
	  if ((img.dataType == ImageTypeDefinitions::greyscale) || (img.dataType == ImageTypeDefinitions::probmap) || (img.dataType == ImageTypeDefinitions::segmentation) || (img.dataType == ImageTypeDefinitions::fmri) || (img.dataType == ImageTypeDefinitions::rgb))
	    {
	      if (fileListWidget->itemAt(fileListWidget->mapFrom(this,event->pos())) == NULL)
		return;
	      
	      QString target = fileListWidget->itemAt(fileListWidget->mapFrom(this,event->pos()))->text();
	      // do nothing if the target is not a model
	      if (imagesInScene[target].dataType != ImageTypeDefinitions::polydata)
		return;
	      
	      //std::cout << "Using volume " << handle.toStdString() << " as texture map for " << target.toStdString() << std::endl;
	      
	      dynamic_cast<SceneModel*>(imagesInScene[target].sd)->AddTextureMap(dynamic_cast<SceneImage*>(imagesInScene[handle].sd)->GetSceneImageView()->GetVTKImageData());
	    }
	}
    }
  else
    {
      event->ignore();
      return;
    }
  
  event->acceptProposedAction();
}

bool VizMainWindow::FileOpenWithName(QString s, ImageTypeDefinitions::DataType t, QString handle)
{
  if (!QFile::exists(s)) return false;
  
  QFileInfo fi(s);
  
  // Make a unique handle for this file if none is given
  if (handle == "")
    {
      handle = fi.baseName();
      QHash<QString, ImageInformation>::const_iterator i = imagesInScene.find(handle);
      unsigned int j = 0;
      while (i != imagesInScene.end()) {
	// It is already in the scene, so add a numeric suffix and repeat
	j++;
	handle = fi.baseName() + QString::number(j);
	i = imagesInScene.find(handle);
      }
    }
  //std::cout << "Setting the handle to be " << handle.toStdString() << std::endl;
  
  ImageInformation siInfo;
  siInfo.fileName = s.toStdString();
  
  if ((t == ImageTypeDefinitions::greyscale) || (t == ImageTypeDefinitions::probmap) || (t == ImageTypeDefinitions::segmentation) || (t == ImageTypeDefinitions::fmri))
    {
      SceneImage *si = new SceneImage(handle);
      assert(si != NULL);
      
      QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
      if (!si->LoadImage(s))
	{
	  std::cout << "Failed to load image called " << s.toStdString() << std::endl;
	  delete si;
	  return false;
	}
      QApplication::restoreOverrideCursor();
      siInfo.sd = si;
    }
  else if ((t == ImageTypeDefinitions::polydata) || (t == ImageTypeDefinitions::scalarpolydata))
    {
      SceneModel *sm = new SceneModel(handle);
      assert(sm != NULL);
      bool scalars;
      if (t == ImageTypeDefinitions::polydata)
	scalars = false;
      else
	scalars = true;
      
      if (!sm->LoadModel(s,scalars))
	{
	  std::cout << "Failed to load model called " << s.toStdString() << std::endl;
	  delete sm;
	  return false;
	}
      siInfo.sd = sm;
    }
  else if (t == ImageTypeDefinitions::tensor)
    {
      SceneTensors *st = new SceneTensors(handle);
      assert(st != NULL);
      if (!st->LoadTensors(s))
	{
	  std::cout << "Failed to load tensor image called " << s.toStdString() << std::endl;
	  delete st;
	  return false;
	}
      siInfo.sd = st;
    }
  else if (t == ImageTypeDefinitions::rgb)
    {
      SceneRGBImage *sr = new SceneRGBImage(handle);
      assert(sr != NULL);
      if (!sr->LoadImage(s))
	{
	  std::cout << "Failed to load RGB image called " << s.toStdString() << std::endl;
	  delete sr;
	  return false;
	}
      siInfo.sd = sr;
    }
  
  // store information about image
  siInfo.dataType = t;
  siInfo.listItem = new QListWidgetItem(handle, fileListWidget);
  siInfo.listItem->setToolTip(s);
  
  if (t == ImageTypeDefinitions::probmap)
    {
      // probability maps are displayed in red
      QBrush brush = siInfo.listItem->foreground();
      brush.setColor(Qt::darkRed);
      brush.setStyle(Qt::SolidPattern);
      siInfo.listItem->setForeground(brush);
    }
  else if (t == ImageTypeDefinitions::segmentation)
    {
      // segmentations are displayed in green
      QBrush brush = siInfo.listItem->foreground();
      brush.setColor(Qt::darkGreen);
      brush.setStyle(Qt::SolidPattern);
      siInfo.listItem->setForeground(brush);
    }
  else if (t == ImageTypeDefinitions::tensor)
    {
      // tensor images are displayed in blue
      QBrush brush = siInfo.listItem->foreground();
      brush.setColor(Qt::darkBlue);
      brush.setStyle(Qt::SolidPattern);
      siInfo.listItem->setForeground(brush);
    }
  else if (t == ImageTypeDefinitions::rgb)
    {
      // RGB images are displayed in gray
      QBrush brush = siInfo.listItem->foreground();
      brush.setColor(Qt::gray);
      brush.setStyle(Qt::SolidPattern);
      siInfo.listItem->setForeground(brush);
    }
  else if (t == ImageTypeDefinitions::polydata || t == ImageTypeDefinitions::scalarpolydata)
    {
      // models are displayed in yellow
      QBrush brush = siInfo.listItem->foreground();
      brush.setColor(Qt::darkYellow);
      brush.setStyle(Qt::SolidPattern);
      siInfo.listItem->setForeground(brush);
    }
  else if (t == ImageTypeDefinitions::fmri)
    {
      // fMRI maps are displayed in magenta
      QBrush brush = siInfo.listItem->foreground();
      brush.setColor(Qt::darkMagenta);
      brush.setStyle(Qt::SolidPattern);
      siInfo.listItem->setForeground(brush);
    }
  
  imagesInScene.insert(handle,siInfo);
  return true;
}

void VizMainWindow::PropagateSliceChange(unsigned int group, int s)
{
  if (inSliceRecursion)
    return;
  
  if (group == 0)
    return;
  
  inSliceRecursion=true;
  ImageInformationMap::iterator it;
  for (it = imagesInScene.begin(); it != imagesInScene.end(); it++)
    {
      ViewerVector::iterator it2;
      for (it2=it.value().viewers.begin();it2!=it.value().viewers.end();it2++)
	{
	  // only process single viewers
	  if (dynamic_cast<VizViewerWindow*>(*it2) == 0)
	    continue;
	  
	  if ((*it2)->GetGroup() == group)
	    dynamic_cast<VizViewerWindow*>(*it2)->Display_slice_changed(s);
	}
    }
  inSliceRecursion=false;
}

void VizMainWindow::PropagateOrientationChange(unsigned int group, eOrientation o)
{
  if (inOrientationRecursion)
    return;
  
  if (group == 0)
    return;
  
  
  inOrientationRecursion=true;
  ImageInformationMap::iterator it;
  for (it = imagesInScene.begin(); it != imagesInScene.end(); it++)
    {
      ViewerVector::iterator it2;
      for (it2=it.value().viewers.begin();it2!=it.value().viewers.end();it2++)
	{
	  // only process single viewers
	  if (dynamic_cast<VizViewerWindow*>(*it2) == 0)
	    continue;
	  
	  if ((*it2)->GetGroup() == group)
	    {
	      switch (o)
		{
		case xy:
		  dynamic_cast<VizViewerWindow*>(*it2)->Display_orientation_xy();
		  break;
		case yz:
		  dynamic_cast<VizViewerWindow*>(*it2)->Display_orientation_yz();
		  break;
		case xz:
		  dynamic_cast<VizViewerWindow*>(*it2)->Display_orientation_xz();
		  break;
		}
	    }
	}
    }
  inOrientationRecursion=false;
}

void VizMainWindow::PropagateCursorPosition(unsigned int group, double pos[3])
{
  if (inCursorRecursion)
    return;
  if (group == 0)
    return;
  
  inCursorRecursion = true;
  ImageInformationMap::iterator it;
  for (it = imagesInScene.begin(); it != imagesInScene.end(); it++)
    {
      ViewerVector::iterator it2;
      for (it2 = it.value().viewers.begin(); it2 != it.value().viewers.end(); it2++)
	{
	  //only process single viewers
	  if (dynamic_cast<VizViewerWindow*>(*it2) == 0)
	    continue;
	  
	  if (((*it2)->GetGroup() == group) && (dynamic_cast<VizViewerWindow*>(*it2)->r2d))
	    dynamic_cast<VizViewerWindow*>(*it2)->r2d->SetCursorPosition(pos);
	}
    }
  
  inCursorRecursion = false;  
}

void VizMainWindow::SetMarkerInViewers(unsigned int _group, double *worldPos)
{
  ImageInformationMap::iterator it;
  for (it = imagesInScene.begin(); it!= imagesInScene.end(); it++)
    {
      
      ViewerVector::iterator it2;
      for (it2 = it.value().viewers.begin();it2 != it.value().viewers.end();it2++)
	{
	  (*it2)->SetMarker(_group,worldPos);
	}
    }
}

void VizMainWindow::RenderAllViewers()
{
  ImageInformationMap::iterator it;
  for (it = imagesInScene.begin(); it!= imagesInScene.end(); it++)
    {
      
      ViewerVector::iterator it2;
      for (it2 = it.value().viewers.begin();it2 != it.value().viewers.end();it2++)
	{
	  (*it2)->Render();
	}
    }
}

VizViewerBase* VizMainWindow::FileListDoubleClick(QListWidgetItem* item)
{
  // after double click, open a new viewer
  VizViewerWindow* w = new VizViewerWindow(this,stereo);
  
  w->show();
  
  bool result;
  if (imagesInScene[item->text()].dataType == ImageTypeDefinitions::tensor)
    result = w->DisplayTensors(dynamic_cast<SceneTensors*>(imagesInScene[item->text()].sd));
  else if (imagesInScene[item->text()].dataType == ImageTypeDefinitions::rgb)
    result = w->DisplayImage(dynamic_cast<SceneRGBImage*>(imagesInScene[item->text()].sd),ImageTypeDefinitions::rgb);
  else if ((imagesInScene[item->text()].dataType == ImageTypeDefinitions::polydata) || (imagesInScene[item->text()].dataType == ImageTypeDefinitions::scalarpolydata))
    result = w->DisplayModelOnly(dynamic_cast<SceneModel*>(imagesInScene[item->text()].sd));
  else
    result = w->DisplayImage(dynamic_cast<SceneImage*>(imagesInScene[item->text()].sd),imagesInScene[item->text()].dataType);
  if (!result)
    {
      delete w;
      return NULL;
    }
  imagesInScene[item->text()].viewers.push_back(w);
  return w;
}

void VizMainWindow::RemoveViewerFromMap(VizViewerBase *viewer)
{
  ImageInformationMap::iterator it;
  for (it=imagesInScene.begin();it!=imagesInScene.end();it++)
    {
      ViewerVector::iterator it2;
      for (it2 = it.value().viewers.begin();it2 != it.value().viewers.end();it2++)
	{
	  if ((*it2) == viewer)
	    {
	      // do not delete viewer object here
	      //delete *it2;
	      it.value().viewers.erase(it2);
	      return;
	    }
	}
    }
}


void VizMainWindow::FileListRightClick(const QPoint &pos)
{
  QListWidgetItem *currentItem = fileListWidget->currentItem();
  if (currentItem)
    {
      if (imagesInScene[currentItem->text()].dataType == ImageTypeDefinitions::polydata)
	{
	  QMenu *m = new QMenu(fileListWidget);
	  QAction *color = m->addAction("Change color...");
	  QAction *opacity = m->addAction("Change opacity...");
	  QAction *save = m->addAction("Save model...");
	  QAction *rotate = m->addAction("Rotate model...");
	  unsigned int count = 0;
	  VizViewerBase *viewer = 0;
	  
	  for (ViewerVector::iterator it = imagesInScene[currentItem->text()].viewers.begin(); it != imagesInScene[currentItem->text()].viewers.end(); it++)
	    {
	      count += (*it)->GetNumberOfRenderersForObject(imagesInScene[currentItem->text()].sd);
	      viewer = (*it);
	    }
	  if (count != 1)
	    {
	      save->setEnabled(false);
	      rotate->setEnabled(false);
	    }
	  m->move(fileListWidget->mapToGlobal(pos));
	  QAction *result = m->exec();
	  if (result == color)
	    {
	      ChooseModelColor(dynamic_cast<SceneModel*>(imagesInScene[currentItem->text()].sd),SceneModel::noScalars);
	    }
	  else if (result == opacity)
	    {
	      double o = QInputDialog::getDouble(this,"Set Opacity","Please enter the desired opacity (0.0 - 1.0)",1.0);
	      if ((o>=0) && (o<=1))
		{
		  dynamic_cast<SceneModel*>(imagesInScene[currentItem->text()].sd)->SetOpacity(o);
		  RenderAllViewers();
		}
	    }
	  else if (result == save)
	    {
	      QString recentDirName = ".";
	      QString s = QFileDialog::getSaveFileName(this,"Save Model",recentDirName,"Models (*.gii *.vtk)");
	      if (s.isEmpty())
		return;
	      viewer->SaveModel(imagesInScene[currentItem->text()].sd,s.toStdString());
	    }
	  else if (result == rotate)
	    {
	      viewer->RotateModel(imagesInScene[currentItem->text()].sd);
	    }
	  delete m;
	}
      else if ((imagesInScene[currentItem->text()].dataType == ImageTypeDefinitions::scalarpolydata) || (imagesInScene[currentItem->text()].dataType == ImageTypeDefinitions::probmap))
	{
	  QMenu *m = new QMenu(fileListWidget);
	  QAction *toggle = m->addAction("Logarithmic lookup table");
	  QAction *color1 = m->addAction("Change lower color...");
	  QAction *color2 = m->addAction("Change upper color...");
	  QAction *value1 = m->addAction("Change lower value...");
	  QAction *value2 = m->addAction("Change upper value...");
	  QAction *save = NULL;
	  if (imagesInScene[currentItem->text()].dataType == ImageTypeDefinitions::scalarpolydata)
	    save = m->addAction("Save model...");
	  m->move(fileListWidget->mapToGlobal(pos));
	  if (imagesInScene[currentItem->text()].dataType == ImageTypeDefinitions::probmap)
	    {
	      if (dynamic_cast<SceneImage*>(imagesInScene[currentItem->text()].sd)->GetLogLookupTable())
		{
		  toggle->setCheckable(true);
		  toggle->setChecked(true);
		}
	    }
	  else
	    if (dynamic_cast<SceneModel*>(imagesInScene[currentItem->text()].sd)->GetLogLookupTable())
	      {
		toggle->setCheckable(true);
		toggle->setChecked(true);
	      }
	  
	  QAction *result = m->exec();
	  if (result == toggle)
	    {
	      if (imagesInScene[currentItem->text()].dataType == ImageTypeDefinitions::scalarpolydata)
		dynamic_cast<SceneModel*>(imagesInScene[currentItem->text()].sd)->ToggleLookupTable();
	      else
		dynamic_cast<SceneImage*>(imagesInScene[currentItem->text()].sd)->ToggleLookupTable();
	      RenderAllViewers();
	    }
	  else if (result == color1)
	    {
	      if (imagesInScene[currentItem->text()].dataType == ImageTypeDefinitions::scalarpolydata)
		ChooseModelColor(dynamic_cast<SceneModel*>(imagesInScene[currentItem->text()].sd),SceneModel::lower);
	      else
		ChooseImageColor(dynamic_cast<SceneImage*>(imagesInScene[currentItem->text()].sd),SceneImage::lower);
	    }
	  else if (result == color2)
	    {
	      if (imagesInScene[currentItem->text()].dataType == ImageTypeDefinitions::scalarpolydata)
		ChooseModelColor(dynamic_cast<SceneModel*>(imagesInScene[currentItem->text()].sd),SceneModel::upper);
	      else
		ChooseImageColor(dynamic_cast<SceneImage*>(imagesInScene[currentItem->text()].sd),SceneImage::upper);
	    }
	  else if (result == value1)
	    {
	      if (imagesInScene[currentItem->text()].dataType == ImageTypeDefinitions::scalarpolydata)
		{
		  SceneModel* model = dynamic_cast<SceneModel*>(imagesInScene[currentItem->text()].sd);
		  double value = QInputDialog::getDouble(this,"Set Lower Value","Please enter the lower value",model->GetLowerValue());
		  model->SetLowerValue(value);
		  RenderAllViewers();
		}
	      else
		{
		  SceneImage* image = dynamic_cast<SceneImage*>(imagesInScene[currentItem->text()].sd);
		  double value = QInputDialog::getDouble(this,"Set Lower Value","Please enter the lower value",image->GetLowerValue());
		  image->SetLowerValue(value);
		  RenderAllViewers();
		}
	    }
	  else if (result == value2)
	    {
	      if (imagesInScene[currentItem->text()].dataType == ImageTypeDefinitions::scalarpolydata)
		{
		  SceneModel* model = dynamic_cast<SceneModel*>(imagesInScene[currentItem->text()].sd);
		  double value = QInputDialog::getDouble(this,"Set Upper Value","Please enter the upper value",model->GetUpperValue());
		  model->SetUpperValue(value);
		  RenderAllViewers();
		}
	      else
		{
		  SceneImage* image = dynamic_cast<SceneImage*>(imagesInScene[currentItem->text()].sd);
		  double value = QInputDialog::getDouble(this,"Set Upper Value","Please enter the upper value",image->GetUpperValue());
		  image->SetUpperValue(value);
		  RenderAllViewers();
		}
	    }
	  else if (result == save)
	    {
	      unsigned int count = 0;
	      VizViewerBase *viewer = 0;
	      
	      for (ViewerVector::iterator it = imagesInScene[currentItem->text()].viewers.begin(); it != imagesInScene[currentItem->text()].viewers.end(); it++)
		{
		  count += (*it)->GetNumberOfRenderersForObject(imagesInScene[currentItem->text()].sd);
		  viewer = (*it);
		}
	      if (count == 1)
		{
		  QString recentDirName = ".";
		  QString s = QFileDialog::getSaveFileName(this,"Save Model",recentDirName,"Models (*.gii *.vtk)");
		  if (s.isEmpty())
		    return;
		  viewer->SaveModel(imagesInScene[currentItem->text()].sd,s.toStdString());
		}
	    }
	  delete m;
	}
      else if (imagesInScene[currentItem->text()].dataType == ImageTypeDefinitions::fmri)
	{
	  QMenu *m = new QMenu(fileListWidget);
	  QAction *value1 = m->addAction("Change lower value...");
	  QAction *value2 = m->addAction("Change upper value...");
	  m->move(fileListWidget->mapToGlobal(pos));
	  QAction *result = m->exec();
	  if (result == value1)
	    {
	      SceneImage *image = dynamic_cast<SceneImage*>(imagesInScene[currentItem->text()].sd);
	      double value = QInputDialog::getDouble(this,"Set Lower Value","Please enter the lower value",image->GetLowerValue());
	      image->SetLowerValue(value);
	      RenderAllViewers();
	    }
	  else if (result == value2)
	    {
	      SceneImage *image = dynamic_cast<SceneImage*>(imagesInScene[currentItem->text()].sd);
	      double value = QInputDialog::getDouble(this,"Set Upper Value","Please enter the upper value",image->GetUpperValue());
	      image->SetUpperValue(value);
	      RenderAllViewers();
	    }
	  delete m;
	}
      else if (imagesInScene[currentItem->text()].dataType == ImageTypeDefinitions::segmentation)
	{
	  QMenu *m = new QMenu(fileListWidget);
	  QAction *model = m->addAction("Create model");
	  m->move(fileListWidget->mapToGlobal(pos));
	  QAction *result = m->exec();
	  if (result == model)
	    {
	    }
	  delete m;
	}
      else if (imagesInScene[currentItem->text()].dataType == ImageTypeDefinitions::tensor)
	{
	  QMenu *m = new QMenu(fileListWidget);
	  QAction *scale = m->addAction("Set scale factor...");
	  m->move(fileListWidget->mapToGlobal(pos));
	  QAction *result = m->exec();
	  if (result == scale)
	    {
	      SceneTensors *tensors = dynamic_cast<SceneTensors*>(imagesInScene[currentItem->text()].sd);
	      double value = QInputDialog::getDouble(this,"Set Scale Factor","Please enter the scale factor",tensors->GetScale());
	      tensors->SetScale(value);
	      RenderAllViewers();
	    }
	  delete m;
	}
    }
}

void VizMainWindow::on_btnDelete_clicked()
{
  QListWidgetItem *currentItem = fileListWidget->currentItem();
  
  if (currentItem)
    {
      fileListWidget->takeItem(fileListWidget->row(currentItem));
      // look for viewers that still display this image and delete them
      while (!imagesInScene[currentItem->text()].viewers.empty())
	{
	  (*imagesInScene[currentItem->text()].viewers.begin())->close();
	}
      delete imagesInScene[currentItem->text()].sd;
      
      imagesInScene.remove(currentItem->text());
      delete currentItem;
    }
}

void VizMainWindow::on_btnFourViewers_clicked()
{
  QListWidgetItem *currentItem = fileListWidget->currentItem();
  if (currentItem)
    {
      VizFourViewers *fv = new VizFourViewers(this);
      fv->show();
      fv->SetGroup(1);
      if (imagesInScene[currentItem->text()].dataType == ImageTypeDefinitions::greyscale)
	if (!fv->DisplayImage(dynamic_cast<SceneImage*>(imagesInScene[currentItem->text()].sd),imagesInScene[currentItem->text()].dataType))
	  {
	    delete fv;
	    return;
	  }
      
      imagesInScene[currentItem->text()].viewers.push_back(fv);
    }
}


ImageInformation VizMainWindow::GetImageInformationForHandle(QString handle)
{
  return imagesInScene[handle];
}


void VizMainWindow::ChooseModelColor(SceneModel *sm, SceneModel::ColorType whichOne)
{
  QColor newColor = QColorDialog::getColor(sm->GetColor(whichOne),this);
  
  if (newColor.isValid())
    sm->SetColor(newColor,whichOne);
}


void VizMainWindow::ChooseImageColor(SceneImage *si, SceneImage::ColorType whichOne)
{
  QColor newColor = QColorDialog::getColor(si->GetColor(whichOne),this);
  
  if (newColor.isValid())
    {
      si->SetColor(newColor,whichOne);
      RenderAllViewers();
    }
}


void VizMainWindow::OpenScene(QString s)
{
  bool inViz = false;
  bool inData = false;
  bool inScene = false;
  int inViewer = -1;
  bool in3D = false;

  VizViewerWindow *currentViewer = 0;
  
  if (!QFile::exists(s)) return;
  
  QFile *file = new QFile(s);
  file->open(QIODevice::ReadOnly);
  QXmlStreamReader *reader = new QXmlStreamReader(file);
  
  reader->setNamespaceProcessing(false);
  
  while (!reader->atEnd())
    {
      reader->readNext();
      if (reader->isStartElement())
	{
	  if (reader->name() == "crlViz")
	    inViz = true;
	  
	  if ((reader->name() == "data") && inViz)
	    inData = true;
	  
	  if ((reader->name() == "scene") && inViz)
	    inScene = true;
	  
	  if ((reader->name() == "volume") && inData)
	    {
	      QString handle = reader->attributes().value("handle").toString();
	      QString file = reader->attributes().value("file").toString();
	      if (reader->attributes().value("type") == "greyscale")
		FileOpenWithName(file,ImageTypeDefinitions::greyscale,handle);
	      else if (reader->attributes().value("type") == "probmap")
		FileOpenWithName(file,ImageTypeDefinitions::probmap,handle);
	      else if (reader->attributes().value("type") == "segmentation")
		FileOpenWithName(file,ImageTypeDefinitions::segmentation,handle);
	      else if (reader->attributes().value("type") == "fmri")
		FileOpenWithName(file,ImageTypeDefinitions::fmri,handle);
	      else if (reader->attributes().value("type") == "tensor")
		FileOpenWithName(file,ImageTypeDefinitions::tensor,handle);
	      else if (reader->attributes().value("type") == "rgb")
		FileOpenWithName(file,ImageTypeDefinitions::rgb,handle);
	    };
	  
	  if ((reader->name() == "model") && inData)
	    {
	      QString handle = reader->attributes().value("handle").toString();
	      QString file = reader->attributes().value("file").toString();
	      if (reader->attributes().value("scalars") == "true")
		{
		  FileOpenWithName(file,ImageTypeDefinitions::scalarpolydata,handle);
		  QColor color;
		  color.setRed(reader->attributes().value("colorLowerR").toString().toUInt());
		  color.setGreen(reader->attributes().value("colorLowerG").toString().toUInt());
		  color.setBlue(reader->attributes().value("colorLowerB").toString().toUInt());
		  dynamic_cast<SceneModel*>(imagesInScene[handle].sd)->SetColor(color,SceneModel::lower);
		  color.setRed(reader->attributes().value("colorUpperR").toString().toUInt());
		  color.setGreen(reader->attributes().value("colorUpperG").toString().toUInt());
		  color.setBlue(reader->attributes().value("colorUpperB").toString().toUInt());
		  dynamic_cast<SceneModel*>(imagesInScene[handle].sd)->SetColor(color,SceneModel::upper);
		  dynamic_cast<SceneModel*>(imagesInScene[handle].sd)->SetLowerValue(reader->attributes().value("lutLowerValue").toString().toFloat());
		  dynamic_cast<SceneModel*>(imagesInScene[handle].sd)->SetUpperValue(reader->attributes().value("lutUpperValue").toString().toFloat());
		  if (reader->attributes().hasAttribute("opacity"))
		    dynamic_cast<SceneModel*>(imagesInScene[handle].sd)->SetOpacity(reader->attributes().value("opacity").toString().toFloat());
		  if (reader->attributes().value("logLUT").toString() == "1")
		    dynamic_cast<SceneModel*>(imagesInScene[handle].sd)->ToggleLookupTable();
		  
		}
	      else
		{
		  FileOpenWithName(file,ImageTypeDefinitions::polydata,handle);
		  QColor color;
		  color.setRed(reader->attributes().value("colorR").toString().toUInt());
		  color.setGreen(reader->attributes().value("colorG").toString().toUInt());
		  color.setBlue(reader->attributes().value("colorB").toString().toUInt());
		  dynamic_cast<SceneModel*>(imagesInScene[handle].sd)->SetColor(color,SceneModel::noScalars);
		  if (reader->attributes().hasAttribute("opacity"))
		    dynamic_cast<SceneModel*>(imagesInScene[handle].sd)->SetOpacity(reader->attributes().value("opacity").toString().toFloat());
		}
	    }
	  
	  if ((reader->name() == "viewer") && inScene)
	    {
	      inViewer = reader->attributes().value("group").toString().toUInt();
	      if (reader->attributes().hasAttribute("in3D"))
		in3D = reader->attributes().value("in3D").toString().toUInt();
	      else
		in3D = false;
	    }
	  
	  if ((reader->name() == "image") && (inViewer != -1))
	    {
	      currentViewer = dynamic_cast<VizViewerWindow*>(this->FileListDoubleClick(imagesInScene[reader->attributes().value("handle").toString()].listItem));
	  currentViewer->SetGroup(inViewer);
	  if (reader->attributes().value("sliceOrientation").toString() == "XY")
	    currentViewer->Display_orientation_xy();
	  else if (reader->attributes().value("sliceOrientation").toString() == "YZ")
	    currentViewer->Display_orientation_yz();
	  else
	    currentViewer->Display_orientation_xz();
	  currentViewer->Display_slice_changed(reader->attributes().value("slicePosition").toString().toUInt());
	  if (in3D)
	    {
	      currentViewer->interactor3D->setChecked(true);
	      currentViewer->Use3DInteractor();
	    }
	}
      
      if ((reader->name() == "overlay") && (currentViewer != 0))
	{
	  currentViewer->DisplayOverlayImage(imagesInScene[reader->attributes().value("handle").toString()].sd,imagesInScene[reader->attributes().value("handle").toString()].dataType);
	  if (reader->attributes().value("sliceOrientation").toString() == "XY")
	    currentViewer->Display_orientation_xy();
	  else if (reader->attributes().value("sliceOrientation").toString() == "YZ")
	    currentViewer->Display_orientation_yz();
	  else
	    currentViewer->Display_orientation_xz();
	  currentViewer->Display_slice_changed(reader->attributes().value("slicePosition").toString().toUInt());
	  currentViewer->blendingSlider->setSliderPosition(reader->attributes().value("opacity").toString().toUInt());
	}
      if ((reader->name() == "model") && (inViewer != -1) && (currentViewer == 0))
	{
	  currentViewer = dynamic_cast<VizViewerWindow*>(this->FileListDoubleClick(imagesInScene[reader->attributes().value("handle").toString()].listItem));
	  if (in3D)
	    currentViewer->Use3DInteractor();
	}
      if ((reader->name() == "model") && (currentViewer != 0))
	{
	  currentViewer->r2d->AddModel(dynamic_cast<SceneModel*>(imagesInScene[reader->attributes().value("handle").toString()].sd));
	}
      if ((reader->name() == "camera") && (currentViewer != 0))
	{
	  double pos[3],fp[3],vu[3];
	  pos[0] = reader->attributes().value("positionX").toString().toFloat();
	  pos[1] = reader->attributes().value("positionY").toString().toFloat();
	  pos[2] = reader->attributes().value("positionZ").toString().toFloat();
	  fp[0] = reader->attributes().value("focalPointX").toString().toFloat();
	  fp[1] = reader->attributes().value("focalPointY").toString().toFloat();
	  fp[2] = reader->attributes().value("focalPointZ").toString().toFloat();
	  vu[0] = reader->attributes().value("viewUpX").toString().toFloat();
	  vu[1] = reader->attributes().value("viewUpY").toString().toFloat();
	  vu[2] = reader->attributes().value("viewUpZ").toString().toFloat();
	  currentViewer->r2d->SetCameraParameters(pos,fp,vu);
	}
      
    }
  else if (reader->isEndElement())
    {
      if ((reader->name() == "scene") && inViz)
	inScene = false;
      
      if ((reader->name() == "data") && inViz)
	inData = false;
      
      if ((reader->name() == "viewer") && inScene)
	{
	  currentViewer = 0;
	  inViewer = -1;
	  in3D = false;
	}
      
      if (reader->name() == "crlViz")
	inViz = false;
    };
  // otherwise do nothing, we ignore everything that is not a start or end element
}
if (reader->hasError())
  {
    std::cout << reader->errorString().toStdString() << std::endl;
  }

delete reader;
file->close();
delete file;
}

void VizMainWindow::SaveScene(QString s)
{
  if (!s.contains('.'))
    s += ".xml";
  
  QFile *file = new QFile(s);
  file->open(QIODevice::WriteOnly);
  QXmlStreamWriter *writer = new QXmlStreamWriter(file);
  
  writer->setAutoFormatting(true);
  writer->writeStartDocument();
  writer->writeStartElement("crlViz");
  writer->writeAttribute("version",Version::GetVersion());
  writer->writeStartElement("data");
  
  ImageInformationMap::iterator it;
  for (it = imagesInScene.begin(); it != imagesInScene.end(); it++)
    {
      if ((it->dataType != ImageTypeDefinitions::polydata) && (it->dataType != ImageTypeDefinitions::scalarpolydata))
	{
	  writer->writeStartElement("volume");
	  writer->writeAttribute("file",QString(it->fileName.c_str()));
	  switch(it->dataType)
	    {
	    case ImageTypeDefinitions::greyscale:
	      writer->writeAttribute("type","greyscale");
	      break;
	    case ImageTypeDefinitions::probmap:
	      writer->writeAttribute("type","probmap");
	      break;
	    case ImageTypeDefinitions::segmentation:
	      writer->writeAttribute("type","segmentation");
	      break;
	    case ImageTypeDefinitions::fmri:
	      writer->writeAttribute("type","fmri");
	      break;
	    case ImageTypeDefinitions::tensor:
	      writer->writeAttribute("type","tensor");
	      break;
	    case ImageTypeDefinitions::rgb:
	      writer->writeAttribute("type","rgb");
	    default:
	      break;
	    }
	  writer->writeAttribute("handle",it.key());
	  writer->writeEndElement();
	}
      else
	{
	  writer->writeStartElement("model");
	  writer->writeAttribute("file",QString(it->fileName.c_str()));
	  if (it->dataType == ImageTypeDefinitions::polydata)
	    {
	      writer->writeAttribute("scalars","false");
	      QString s;
	      writer->writeAttribute("colorR",s.setNum(dynamic_cast<SceneModel*>(it->sd)->GetColor(SceneModel::noScalars).red()));
	      writer->writeAttribute("colorG",s.setNum(dynamic_cast<SceneModel*>(it->sd)->GetColor(SceneModel::noScalars).green()));
	      writer->writeAttribute("colorB",s.setNum(dynamic_cast<SceneModel*>(it->sd)->GetColor(SceneModel::noScalars).blue()));
	      writer->writeAttribute("opacity",s.setNum(dynamic_cast<SceneModel*>(it->sd)->GetOpacity()));
	    }
	  else
	    {
	      writer->writeAttribute("scalars","true");
	      QString s;
	      writer->writeAttribute("colorLowerR",s.setNum(dynamic_cast<SceneModel*>(it->sd)->GetColor(SceneModel::lower).red()));
	      writer->writeAttribute("colorLowerG",s.setNum(dynamic_cast<SceneModel*>(it->sd)->GetColor(SceneModel::lower).green()));
	      writer->writeAttribute("colorLowerB",s.setNum(dynamic_cast<SceneModel*>(it->sd)->GetColor(SceneModel::lower).blue()));
	      writer->writeAttribute("colorUpperR",s.setNum(dynamic_cast<SceneModel*>(it->sd)->GetColor(SceneModel::upper).red()));
	      writer->writeAttribute("colorUpperG",s.setNum(dynamic_cast<SceneModel*>(it->sd)->GetColor(SceneModel::upper).green()));
	      writer->writeAttribute("colorUpperB",s.setNum(dynamic_cast<SceneModel*>(it->sd)->GetColor(SceneModel::upper).blue()));
	      writer->writeAttribute("lutLowerValue",s.setNum(dynamic_cast<SceneModel*>(it->sd)->GetLowerValue()));
	      writer->writeAttribute("lutUpperValue",s.setNum(dynamic_cast<SceneModel*>(it->sd)->GetUpperValue()));
	      writer->writeAttribute("opacity",s.setNum(dynamic_cast<SceneModel*>(it->sd)->GetOpacity()));
	      writer->writeAttribute("logLUT",s.setNum(dynamic_cast<SceneModel*>(it->sd)->GetLogLookupTable()));
	    }
	  writer->writeAttribute("handle",it.key());
	  writer->writeEndElement();
	}
    }
  
  writer->writeEndElement(); // </data>
  
  writer->writeStartElement("scene");
  // Information about open viewers, slice positions, etc. goes here
  
  std::set<VizViewerBase*> visitedViewers;
  unsigned int viewerId = 0;
  
  for (it = imagesInScene.begin(); it != imagesInScene.end(); it++)
    {
      ViewerVector::iterator it2;
      for (it2 = it->viewers.begin(); it2 != it->viewers.end(); it2++)
	{
	  if (visitedViewers.find(*it2) != visitedViewers.end()) continue;
	  visitedViewers.insert(*it2);
	  VizViewerWindow* viewer = dynamic_cast<VizViewerWindow*>(*it2);
	  if (!viewer) continue;
	  viewerId++;
	  writer->writeStartElement("viewer");
	  writer->writeAttribute("id",QString("%1").arg(viewerId));
	  writer->writeAttribute("group",QString("%1").arg(viewer->group));
	  writer->writeAttribute("in3D",QString("%1").arg(viewer->r2d->Using3DInteractor()));
	  if (viewer->name != "")
	    {
	      writer->writeStartElement("image");
	      writer->writeAttribute("handle",viewer->name);
	      writer->writeAttribute("sliceOrientation",viewer->r2d->GetSliceOrientation().c_str());
	      writer->writeAttribute("slicePosition",QString("%1").arg(viewer->r2d->GetSlice()));
	      writer->writeEndElement(); // </volume>
	    }
	  if (viewer->overlay_name != "")
	    {
	      writer->writeStartElement("overlay");
	      writer->writeAttribute("handle",viewer->overlay_name);
	      writer->writeAttribute("sliceOrientation",viewer->r2d->GetSliceOrientation().c_str());
	      writer->writeAttribute("slicePosition",QString("%1").arg(viewer->r2d->GetSlice()));
	      writer->writeAttribute("opacity",QString("%1").arg(viewer->blendingSlider->value()));
	      writer->writeEndElement(); // </overlay>
	    }
	  
	  std::set<SceneModel*> models = viewer->r2d->GetModels();
	  for (std::set<SceneModel*>::iterator it = models.begin(); it != models.end(); it++)
	    {
	      writer->writeStartElement("model");
	      writer->writeAttribute("handle",(*it)->GetName());
	      writer->writeEndElement(); // </model>
	    }
	  writer->writeStartElement("camera");
	  double pos[3],fp[3],vu[3];
	  viewer->r2d->GetCameraParameters(pos,fp,vu);
	  writer->writeAttribute("positionX",QString("%1").arg(vu[0]));
	  writer->writeAttribute("positionY",QString("%1").arg(vu[1]));
	  writer->writeAttribute("positionZ",QString("%1").arg(vu[2]));
	  writer->writeAttribute("focalPointX",QString("%1").arg(fp[0]));
	  writer->writeAttribute("focalPointY",QString("%1").arg(fp[1]));
	  writer->writeAttribute("focalPointZ",QString("%1").arg(fp[2]));
	  writer->writeAttribute("viewUpX",QString("%1").arg(vu[0]));
	  writer->writeAttribute("viewUpY",QString("%1").arg(vu[1]));
	  writer->writeAttribute("viewUpZ",QString("%1").arg(vu[2]));
	  writer->writeEndElement(); // </camera>
	  
	  writer->writeEndElement(); // </viewer>
	}
    }
  
  writer->writeEndElement(); // </scene>
  writer->writeEndElement(); // </crlViz>
  writer->writeEndDocument();
  
  delete writer;
  file->close();
  delete file;
}


// Wizards

void VizMainWindow::on_actionSTAPLE_triggered()
{
  VizWizardSTAPLE *wizard = new VizWizardSTAPLE(this);
  
  wizard->show();
}

void VizMainWindow::on_actionConnect_to_VectorVision_triggered()
{
#ifdef HAVE_VVLINK
  QString host = QInputDialog::getText(this,"crlViz","Host name or IP address:",QLineEdit::Normal,"crlnav.tch.harvard.edu");
  QString password = QInputDialog::getText(this,"crlViz","Password:",QLineEdit::Normal,"711DE3");
  if (host.isEmpty() || password.isEmpty())
    return;
  
  if (!server) server = vvlink::VVLConnection::New();
  // connect to Brainlab in OR 25
  try
    {
      server->ConnectTo(host.toStdString(),"crlNav",password.toStdString(),22222,1500);
    }
  catch(vvlink::VVLEx &ex)
    {
      QMessageBox::critical(this,"crlNav","Connection could not be established: "+tr(ex.What().c_str()));
      return;
    }
  
  if (!coord)
    coord = vvlink::VVLCoordinateSystem::New(vvlink::VVLCoordinateSystem::PATIENT);
  
  actionDisconnect->setEnabled(true);
  actionConnect_to_VectorVision->setEnabled(false);
  
  // retrieve instrument list
  vvlink::VVLQGetAvailableInstruments *instruments = vvlink::VVLQGetAvailableInstruments::New();
  instruments->Execute(server);
  
  menuTools->clear();
  for (int i =0; i < instruments->GetNumberOfInstruments(); i++)
    {
      QAction *action = menuTools->addAction(tr(instruments->GetInstrumentID(i).c_str()));
      action->setCheckable(true);
      tools.push_back(action);
    }
  
  connected = true;
#endif
}

void VizMainWindow::on_actionDisconnect_triggered()
{
#ifdef HAVE_VVLINK
  if (connected)
    {
      while (inTimer) {}; // wait until we are out of the timer loop
      inTimer = true; // prevent the timer loop from executing while we are disconnecting
      connected = false;
      displayInVectorVision = NULL;
      if (trackedToolNotifier)
	{
	  trackedToolNotifier->UnregisterFromServer();
	  trackedToolNotifier->Delete();
	  trackedToolNotifier = NULL;
	  menuTools->clear();
	  menuTools->addAction("not connected")->setEnabled(false);
	}
      if (view)
	view->SetImage(NULL);
      coord->Delete();
      inTimer = false; // now it's safe to run again
      actionDisconnect->setEnabled(false);
      actionConnect_to_VectorVision->setEnabled(true);
    }
#endif
}

void VizMainWindow::on_actionGet_Image_Data_triggered()
{
#ifdef HAVE_VVLINK
  if (!connected)
    return;
  
  vvlink::VVLQGetSliceSets *command = vvlink::VVLQGetSliceSets::New();
  command->SetCoordinateSystem(coord);
  command->Execute(server);
  
  const vvlink::VVLSliceSet *sliceSet = command->GetSliceSet(0);
  vtkImageData *image = sliceSet->GetOrthogonalVolume();
  
  vtkImageCast *cast = vtkImageCast::New();
  cast->SetInput(sliceSet->GetOrthogonalVolume());
  cast->SetOutputScalarTypeToFloat();
  cast->Update();
  image = cast->GetOutput();
  
  QString handle = tr(sliceSet->GetName().c_str());
  QHash<QString, ImageInformation>::const_iterator i = imagesInScene.find(handle);
  unsigned int j = 0;
  while (i != imagesInScene.end())
    {
      j++;
      handle = tr(sliceSet->GetName().c_str()) + QString::number(j);
      i = imagesInScene.find(handle);
    }
  
  ImageInformation vvInfo;
  SceneImage *si = new SceneImage(handle);
  vvInfo.fileName = "";
  vvInfo.dataType = ImageTypeDefinitions::greyscale;
  vvInfo.sd = si;
  
  ImageTypeDefinitions::VTK2ITKConnectorFilterType::Pointer VTK2ITKConnector = ImageTypeDefinitions::VTK2ITKConnectorFilterType::New();
  VTK2ITKConnector->SetInput(image);
  VTK2ITKConnector->Update();
  
  // make a copy of the image to store in the SceneImage object
  ImageTypeDefinitions::ImageType::Pointer itkImage = ImageTypeDefinitions::ImageType::New();
  itkImage->SetRegions(VTK2ITKConnector->GetOutput()->GetLargestPossibleRegion());
  itkImage->SetOrigin(VTK2ITKConnector->GetOutput()->GetOrigin());
  itkImage->SetSpacing(VTK2ITKConnector->GetOutput()->GetSpacing());
  ImageTypeDefinitions::ImageType::DirectionType direction;
  vtkMatrix4x4* vvDir = sliceSet->GetOrthogonalTransformation();
  direction[0][0] = -vvDir->Element[0][0]; // flip x
  direction[0][1] = vvDir->Element[0][1];
  direction[0][2] = vvDir->Element[0][2];
  direction[1][0] = vvDir->Element[1][0];
  direction[1][1] = -vvDir->Element[1][1]; // flip y
  direction[1][2] = vvDir->Element[1][2];
  direction[2][0] = vvDir->Element[2][0];
  direction[2][1] = vvDir->Element[2][1];
  direction[2][2] = vvDir->Element[2][2];
  itkImage->SetDirection(direction);
  itkImage->Allocate();
  
  itk::ImageRegionConstIterator<ImageTypeDefinitions::ImageType> iit(VTK2ITKConnector->GetOutput(),VTK2ITKConnector->GetOutput()->GetLargestPossibleRegion());
  itk::ImageRegionIterator<ImageTypeDefinitions::ImageType> oit(itkImage,itkImage->GetLargestPossibleRegion());
  while (!iit.IsAtEnd())
    {
      oit.Set(iit.Get());
      ++iit;
      ++oit;
    }
  
  itk::OrientImageFilter<ImageTypeDefinitions::ImageType,ImageTypeDefinitions::ImageType>::Pointer orienter = itk::OrientImageFilter<ImageTypeDefinitions::ImageType,ImageTypeDefinitions::ImageType>::New();
  orienter->SetInput(itkImage);
  orienter->UseImageDirectionOn();
  orienter->SetDesiredCoordinateOrientationToAxial();
  orienter->Update();
  si->SetImage(orienter->GetOutput());
  
  vvInfo.listItem = new QListWidgetItem(handle, fileListWidget);
  vvInfo.listItem->setToolTip("VectorVision image");
  imagesInScene.insert(handle,vvInfo);
#endif
}

void VizMainWindow::timerCallback()
{
#ifdef HAVE_VVLINK
  if (inTimer)
    return;
  
  inTimer = true;
  
  if (displayInVectorVision && connected)
    {
      // send a viewer window to VectorVision
      wtif->Modified();
      wtif->Update();
      view->Execute(server);
    }
  
  if (trackedToolNotifier && connected)
    {
      // track tool position
      trackedToolNotifier->GetInstrument(tool);
      std::vector<Render2D3D*>::iterator it;
      double tip[4];
      double direction[4], perp1[4], perp2[4];
      for (it = trackToolViewers.begin(); it != trackToolViewers.end(); it++)
	if (tool->IsValid())
	  {
	    tool->GetTip(tip);
	    tool->GetDirection(direction);
	    // translation
	    lastToolPosition->SetElement(0,3,-tip[0]);
	    lastToolPosition->SetElement(1,3,-tip[1]);
	    lastToolPosition->SetElement(2,3,tip[2]);
	    // rotation
	    direction[0] = direction[0];
	    direction[1] = direction[1];
	    direction[2] = -direction[2];
	    vtkMath::Perpendiculars(direction,perp1,perp2,0);
	    lastToolPosition->SetElement(0,0,direction[0]);
	    lastToolPosition->SetElement(1,0,direction[1]);
	    lastToolPosition->SetElement(2,0,direction[2]);
	    lastToolPosition->SetElement(0,1,perp1[0]);
	    lastToolPosition->SetElement(1,1,perp1[1]);
	    lastToolPosition->SetElement(2,1,perp1[2]);
	    lastToolPosition->SetElement(0,2,perp2[0]);
	    lastToolPosition->SetElement(1,2,perp2[1]);
	    lastToolPosition->SetElement(2,2,perp2[2]);
	    (*it)->SetToolPosition(lastToolPosition,true);
	  }
	else
	  (*it)->SetToolPosition(lastToolPosition,false);
    }
  
  inTimer = false;
#endif
}

#ifdef HAVE_VVLINK
bool VizMainWindow::DisplayInVectorVision(Render2D3D *image)
{
  if (!image && displayInVectorVision)
    {
      view->SetImage(NULL);
      displayInVectorVision = NULL;
      return true;
    }
  
  if (displayInVectorVision)
    return false;
  
  displayInVectorVision = image;
  wtif->SetInput(image->ren->GetRenderWindow());
  view->SetImage(wtif->GetOutput());
  return true;
}

void VizMainWindow::TrackToolOn(Render2D3D *viewer)
{
  trackToolViewers.push_back(viewer);
}

void VizMainWindow::TrackToolOff(Render2D3D * viewer)
{
  std::vector<Render2D3D*>::iterator it;
  for (it = trackToolViewers.begin(); it != trackToolViewers.end(); it++)
    if ((*it) == viewer)
      {
	trackToolViewers.erase(it);
	return;
      }
}
#endif

void VizMainWindow::toolsMenuCallback(QAction *action)
{
#ifdef HAVE_VVLINK
  if (!connected)
    return;
  
  if (trackedToolNotifier)
    {
      trackedToolNotifier->UnregisterFromServer();
      trackedToolNotifier->Delete();
    }
  
  trackedToolNotifier = vvlink::VVLNTrackedInstrument::New();
  trackedToolNotifier->SetCoordinateSystem(coord);
  trackedToolNotifier->SetInstrumentID(action->text().toStdString());
  trackedToolNotifier->RegisterOnServer(server);
  
  std::vector<QAction*>::iterator it;
  for (it = tools.begin(); it != tools.end(); it++)
    {
      (*it)->setChecked(false);
    }
  action->setChecked(true);
#endif
}
