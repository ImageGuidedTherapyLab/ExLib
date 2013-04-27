#include <Render2D3D.h>

#include <QtGui>
#include <QMainWindow>
#include <VizViewerWindow.h>
#include <VizMainWindow.h>
#include <CrosshairActor.h>

#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageActor.h>
#include <vtkImageReslice.h>
#include <vtkCamera.h>
#include <vtkColorTransferFunction.h>
#include <vtkCallbackCommand.h>
#include <vtkInteractorObserver.h>
#include <vtkImageMapToColors.h>
#include <vtkImageMapToWindowLevelColors.h>
#include <vtkLookupTable.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkPolyDataMapper.h>
#include <vtkLogLookupTable.h>
#include <vtkProperty.h>
#include <vtkOutlineFilter.h>
#include <vtkCursor3D.h>
#include <vtkAnnotatedCubeActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPropPicker.h>
#include <vtkCellPicker.h>
#include <vtkBoxWidget.h>
#include <vtkProperty2D.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkMatrix4x4.h>
#include <vtkDepthSortPolyData.h>
#include <vtkWorldPointPicker.h>
#include <vtkOBJExporter.h>
#include <vtkRIBExporter.h>
#include <vtkCell.h>
#include <vtkMath.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPointData.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkImageBlend.h>
#include <vtkConeSource.h>
#include <vtkCylinderSource.h>
#include <vtkArrowSource.h>

// Relevant Viz classes
#include <ImageOperations.h>
#include <SceneData.h>
#include <SceneImage.h>
#include <SceneImageView.h>
#include <SceneRGBImage.h>
#include <SceneModel.h>
#include <SceneTensors.h>
#include <vtkGIFTIWriter.h>

// QVTK Widget on the Mac brings in Carbon, which can conflict with some of
// The Qt header files. Therefore we bring this in last.
#include <QVTKWidget.h>


// round is used in this class, but does not exist under Windows
#ifdef WIN32
double round(double x)
{
  return floor(x + .5);
}
#endif


vtkStandardNewMacro(Render2D3D);
vtkCxxRevisionMacro(Render2D3D,"$Revision: 1.0 $");

Render2D3D::Render2D3D()
{
  ren = 0;
  interactor = 0;
  qvtk = 0;
  siv = 0;
  overlay_siv = 0;
  scalarDataType[0] = ImageTypeDefinitions::greyscale;
  scalarDataType[1] = ImageTypeDefinitions::greyscale;
  viewWin = 0;
  orientation = ImageTypeDefinitions::axial;
  vtkactor[ImageTypeDefinitions::axial] = 0;
  vtkactor[ImageTypeDefinitions::sagittal] = 0;
  vtkactor[ImageTypeDefinitions::coronal] = 0;
  lastSlice[ImageTypeDefinitions::axial] = 0;
  lastSlice[ImageTypeDefinitions::sagittal] = 0;
  lastSlice[ImageTypeDefinitions::coronal] = 0;
  vtkblend = 0;
  vtkmapper = 0;
  vtkoverlaymapper = 0;
  imageData = 0;
  overlayImageData = 0;
  in3D = false;
  tensors = 0;
  tensorActor = 0;
  outlineActor = 0;
  cursor = 0;
  cursorActor = 0;
  cursorVisible = false;
  mousePointerVisible = true;
  cube = 0;
  orientationMarkers = 0;
  roiBox = 0;
  currentPointClass = 0;
  currentColor[0] = 1.0; currentColor[1] = 0.0; currentColor[2] = 0.0; // red
  stereo = false;
  labelBlending = true;
  
  // Create IGT tool
  igtTool = vtkActor::New();
  /*vtkConeSource *cone = vtkConeSource::New();
    cone->SetHeight(10);
    cone->SetCenter(5,0,0);
    cone->SetDirection(1,0,0);
    cone->SetRadius(5);
    cone->SetResolution(20);
    
    vtkCylinderSource *cyl = vtkCylinderSource::New();
    cyl->CappingOn();
    cyl->SetHeight(60);
    cyl->SetCenter(40,0,0);
    cyl->SetRadius(5);
    cyl->SetResolution(20);
    
    vtkAppendPolyData *append = vtkAppendPolyData::New();
    append->AddInput(cone->GetOutput());
    append->AddInput(cyl->GetOutput());
    append->Update();
  */
  vtkArrowSource *arrow = vtkArrowSource::New();
  arrow->InvertOn();
  arrow->SetTipRadius(3);
  arrow->SetShaftRadius(3);
  vtkPolyDataMapper *pdm = vtkPolyDataMapper::New();
  pdm->SetInput(arrow->GetOutput());
  igtTool->SetMapper(pdm);
  igtTool->GetProperty()->SetColor(1,0,0);
  pdm->Delete();
  //append->Delete();
  //cyl->Delete();
  //cone->Delete();
  arrow->Delete();
}


Render2D3D::~Render2D3D()
{
  // We don't delete the QVTKWidget - this class simply uses it as a place
  // to display, it doesn't manage allocation and deallocation of QVTKWidget
  if (qvtk)
    qvtk->GetRenderWindow()->SetInteractor(NULL);
  if (ren)
    qvtk->GetRenderWindow()->RemoveRenderer(ren);
  if (interactor) interactor->Delete();
  if (vtkactor[ImageTypeDefinitions::axial]) vtkactor[ImageTypeDefinitions::axial]->Delete();
  if (vtkactor[ImageTypeDefinitions::sagittal]) vtkactor[ImageTypeDefinitions::sagittal]->Delete();
  if (vtkactor[ImageTypeDefinitions::coronal]) vtkactor[ImageTypeDefinitions::coronal]->Delete();
  if (vtkblend) vtkblend->Delete();
  if (vtkmapper) vtkmapper->Delete();
  if (vtkoverlaymapper) vtkoverlaymapper->Delete();
  if (tensorActor) tensorActor->Delete();
  if (outlineActor) outlineActor->Delete();
  if (cursorActor) cursorActor->Delete();
  if (cursor) cursor->Delete();
  if (cube) cube->Delete();
  if (orientationMarkers) orientationMarkers->Delete();
  if (igtTool) igtTool->Delete();
  if (roiBox) DisableROIBox();
  modelActors.clear();
  ClearPoints();
  if (ren) ren->Delete();
}


void Render2D3D::SetViewerWindow(VizViewerBase *_view)
{
  viewWin = _view;
}


inline bool Render2D3D::Initialized()
{
  return (ren != 0);
}


bool Render2D3D::InitializeRendering(QVTKWidget *renderonthis)
{
  if (!renderonthis) return false;
  if (!ren) {
    ren = vtkRenderer::New();
    ren->SetUseDepthPeeling(1);
    ren->SetMaximumNumberOfPeels(400);
    ren->SetOcclusionRatio(0.1);
    if (stereo)
      {
	renderonthis->GetRenderWindow()->SetStereoCapableWindow(1);
	renderonthis->GetRenderWindow()->SetStereoRender(1);
	//renderonthis->GetRenderWindow()->SetStereoTypeToDresden();
	renderonthis->GetRenderWindow()->SetStereoTypeToCrystalEyes();
      }
    renderonthis->GetRenderWindow()->SetAlphaBitPlanes(1);
    renderonthis->GetRenderWindow()->SetMultiSamples(0);
    renderonthis->GetRenderWindow()->AddRenderer(ren);
    interactor = QVTKInteractor::New();
    interactor->SetRenderWindow(renderonthis->GetRenderWindow());
    PointType crosshair;
    crosshair.actor = new CrosshairActor;
    crosshair.x = 0; crosshair.y = 0; crosshair.z = 0;
    ren->AddViewProp(crosshair.actor->GetActor());
    crosshair.actor->GetActor()->SetVisibility(false);
    crosshairVisible = false;
    // the crosshair actor is always the first element in the points vector
    points.push_back(crosshair);
  } else {
    std::cout << "view already constructed." << std::endl;
  }
  
  ren->GetActiveCamera()->ParallelProjectionOn();
  ren->SetBackground(255,255,255);
  
  // Now we have qvtk being allocated by the UI constructor, so this is a
  // tautology.
  qvtk = renderonthis;
  Set2DInteractor();
  
  cursor = vtkCursor3D::New();
  cursor->OutlineOff();
  //cursor->XShadowsOff();
  //cursor->YShadowsOff();
  //cursor->ZShadowsOff();
  cursorActor = vtkActor::New();
  vtkPolyDataMapper *cursorMapper = vtkPolyDataMapper::New();
  cursorMapper->SetInput(cursor->GetOutput());
  cursorActor->SetMapper(cursorMapper);
  cursorMapper->Delete();
  
  ren->GetRenderWindow()->Render();
  
  return true;
}


bool Render2D3D::SetScalarInput(int number, SceneData *si)
{
  if (number == 1 && si == NULL)
    {
      // if this function is called with a NULL pointer for the overlay image,
      // just remove the current overlay from the render window, do nothing else
      if (overlay_siv)
	{
	  vtkblend->RemoveInputConnection(0,vtkoverlaymapper->GetOutputPort(0));
	}
      overlay_siv = 0; // Reset the pointer to have no overlay.
      // There may be other state that needs to be reset to ensure
      // adding a second overlay does the correct thing here.
      return true;
    }
  
  if (!si) return false;
  if (number < 0 || number > 1) return false;
  if (!siv && number == 1) return false; // don't allow overlays if there no actual scalar image already (e.g. in case of tensors)
  
  if (!Initialized()) {
    std::cerr << "rendering widget not initialized." << std::endl;
    return false;
  }
  
  SceneImageView *currentSiv;
  vtkImageMapToColors *currentMapper;
  vtkImageData *currentImageData;
  bool resampleOverlay = false;
  
  if (number == 0)
    {
      // we want to set a new primary image
      if (scalarDataType[number] == ImageTypeDefinitions::rgb)
	siv =dynamic_cast<SceneRGBImage*>(si)->GetSceneImageView();
      else
	siv = dynamic_cast<SceneImage*>(si)->GetSceneImageView();
      currentSiv = siv;
      // delete any already existing actors, blends and mappers
      if (vtkactor[ImageTypeDefinitions::axial]) vtkactor[ImageTypeDefinitions::axial]->Delete();
      if (vtkactor[ImageTypeDefinitions::sagittal]) vtkactor[ImageTypeDefinitions::sagittal]->Delete();
      if (vtkactor[ImageTypeDefinitions::coronal]) vtkactor[ImageTypeDefinitions::coronal]->Delete();
      if (vtkblend) vtkblend->Delete();
      if (vtkmapper) vtkmapper->Delete();
      
      vtkactor[ImageTypeDefinitions::axial] = vtkImageActor::New();
      vtkactor[ImageTypeDefinitions::sagittal] = vtkImageActor::New();
      vtkactor[ImageTypeDefinitions::coronal] = vtkImageActor::New();
      vtkblend = vtkImageBlend::New();
      vtkblend->SetOpacity(0,1.0);
      vtkblend->SetOpacity(1,0.0);
      if (scalarDataType[number] == ImageTypeDefinitions::greyscale)
	vtkmapper = vtkImageMapToWindowLevelColors::New();
      else
	vtkmapper = vtkImageMapToColors::New();
      
      currentMapper = vtkmapper;
      imageData = siv->GetVTKImageData();
      if (!imageData) {
	std::cout << "Failed to find any image data in the SceneImageView." << std::endl;
	return false;
      }
      
      currentImageData = imageData;
      orientation = ImageTypeDefinitions::axial;
    }
  else
    {
      // Check to see if we already have an overlay.
      if (overlay_siv != 0) {
        // Since we have one already, remove it, so the new overlay can be
        // properly displayed.
	vtkblend->RemoveInputConnection(0,vtkoverlaymapper->GetOutputPort(0));
        overlay_siv = 0;
      }
      // we want to set a new overlay image
      if (scalarDataType[number] == ImageTypeDefinitions::rgb)
	overlay_siv = dynamic_cast<SceneRGBImage*>(si)->GetSceneImageView();
      else
	overlay_siv = dynamic_cast<SceneImage*>(si)->GetSceneImageView();
      overlayImageData = overlay_siv->GetVTKImageData();
      currentImageData = overlayImageData;
      if (!overlayImageData) {
	std::cout << "Failed to find any image data in the SceneImageView." << std::endl;
	return false;
      }
      
      currentSiv = overlay_siv;
      // delete any already existing mappers
      if (vtkoverlaymapper) vtkoverlaymapper->Delete();
      
      //check if the extent is the same for both images
      //if not, set resample flag
      int *e0, *e1;
      e0 = imageData->GetExtent();
      // we need to update the overlay image data object in order to get the extent
      overlayImageData->Update();
      e1 = overlayImageData->GetExtent();
      if (e0[0]!=e1[0] || e0[1]!=e1[1] || e0[2]!=e1[2] || e0[3]!=e1[3] || e0[4]!=e1[4] || e0[5]!=e1[5])
	{
	  std::cout << "Overlay image will be automatically resampled." << std::endl;
          std::cout << "If the voxel dimensions are different, this is a problem." << std::endl;
	  resampleOverlay = true;
	}
      
      if (scalarDataType[number] == ImageTypeDefinitions::greyscale)
	vtkoverlaymapper = vtkImageMapToWindowLevelColors::New();
      else
	vtkoverlaymapper = vtkImageMapToColors::New();
      
      currentMapper = vtkoverlaymapper;
    }
  
  currentSiv->SetDataType(scalarDataType[number]);
  
  if (!resampleOverlay)
    currentMapper->SetInput(currentImageData);
  else
    {
      // This should probably be an ITK resample.
      //this is an overlay image, and it needs to be resampled
      vtkImageReslice *resampler = vtkImageReslice::New();
      resampler->SetInput(currentImageData);
      resampler->SetInformationInput(imageData);
      resampler->Update();
      overlayImageData = resampler->GetOutput();
      currentImageData = overlayImageData;
      currentMapper->SetInput(overlayImageData);
      resampler->Delete();
    }
  if (scalarDataType[number] == ImageTypeDefinitions::greyscale)
    {
      dynamic_cast<vtkImageMapToWindowLevelColors*>(currentMapper)->SetWindow(currentSiv->GetInitialWindow());
      dynamic_cast<vtkImageMapToWindowLevelColors*>(currentMapper)->SetLevel(currentSiv->GetInitialLevel());
    }
  
  vtkblend->RemoveAllInputs();
  vtkblend->AddInputConnection(0,vtkmapper->GetOutputPort(0));
  if (overlay_siv)
    vtkblend->AddInputConnection(0,vtkoverlaymapper->GetOutputPort(0));
  
  vtkactor[ImageTypeDefinitions::axial]->SetInput(vtkblend->GetOutput());
  vtkactor[ImageTypeDefinitions::sagittal]->SetInput(vtkblend->GetOutput());
  vtkactor[ImageTypeDefinitions::coronal]->SetInput(vtkblend->GetOutput());
  
  // rotate the image actor according to the original image directions (vtkImageData does not support arbitrary directions)
  /*vtkactor[ImageTypeDefinitions::axial]->SetOrigin(currentSiv->GetTranslation());
    vtkactor[ImageTypeDefinitions::sagittal]->SetOrigin(currentSiv->GetTranslation());
    vtkactor[ImageTypeDefinitions::coronal]->SetOrigin(currentSiv->GetTranslation());
    vtkactor[ImageTypeDefinitions::axial]->SetUserMatrix(currentSiv->GetRotation());
    vtkactor[ImageTypeDefinitions::sagittal]->SetUserMatrix(currentSiv->GetRotation());
    vtkactor[ImageTypeDefinitions::coronal]->SetUserMatrix(currentSiv->GetRotation());*/
  
  if (number != 1)
    {
      int *extent = imageData->GetWholeExtent();
      vtkactor[ImageTypeDefinitions::axial]->SetDisplayExtent(extent[0],extent[1],extent[2],extent[3],extent[4],extent[4]);
      vtkactor[ImageTypeDefinitions::sagittal]->SetDisplayExtent(extent[0],extent[0],extent[2],extent[3],extent[4],extent[5]);
      vtkactor[ImageTypeDefinitions::coronal]->SetDisplayExtent(extent[0],extent[1],extent[2],extent[2],extent[4],extent[5]);
    }
  
  vtkLookupTable *vtklut;
  if ((number > 0) && (siv->GetDataType() == ImageTypeDefinitions::greyscale))
    vtklut = currentSiv->GetLookupTable();
  else
    vtklut = currentSiv->GetLookupTable(1.0);
  currentMapper->SetLookupTable(vtklut);
  currentMapper->SetOutputFormatToRGBA();
  currentMapper->Update();
  
  cursor->SetModelBounds(imageData->GetBounds());
  
  vtkactor[ImageTypeDefinitions::axial]->SetVisibility(0);
  ren->AddViewProp(vtkactor[ImageTypeDefinitions::axial]);
  vtkactor[ImageTypeDefinitions::sagittal]->SetVisibility(0);
  ren->AddViewProp(vtkactor[ImageTypeDefinitions::sagittal]);
  vtkactor[ImageTypeDefinitions::coronal]->SetVisibility(0);
  ren->AddViewProp(vtkactor[ImageTypeDefinitions::coronal]);
  
  vtkactor[orientation]->SetVisibility(1);
  
  if (number == 1)
    {
      vtkblend->SetOpacity(1,1.0-vtkblend->GetOpacity(0));
    }
  
  if (number != 1)
    {
      //ResetCamera();
      // VTK 5.8.0 ignored ResetCamera
      SetSliceOrientationToXY();
    }
  
  ren->GetRenderWindow()->Render();
  return true; // success
}


void Render2D3D::Render(void) const {
  if (!ren) return;
  ren->GetRenderWindow()->Render();
}


void Render2D3D::SetSliceOrientationToXY() {
  if (!siv && !tensors) return;
  
  lastSlice[orientation] = GetSlice();
  
  if (siv)
    {
      vtkactor[ImageTypeDefinitions::axial]->SetVisibility(1);
      vtkactor[ImageTypeDefinitions::sagittal]->SetVisibility(0);
      vtkactor[ImageTypeDefinitions::coronal]->SetVisibility(0);
    }
  
  orientation = ImageTypeDefinitions::axial;
  SetSlice(lastSlice[ImageTypeDefinitions::axial]);
  
  ResetCamera();
}


void Render2D3D::SetSliceOrientationToYZ() {
  if (!siv && !tensors) return;
  
  lastSlice[orientation] = GetSlice();
  
  if (siv)
    {
      vtkactor[ImageTypeDefinitions::axial]->SetVisibility(0);
      vtkactor[ImageTypeDefinitions::sagittal]->SetVisibility(1);
      vtkactor[ImageTypeDefinitions::coronal]->SetVisibility(0);
    }
  
  orientation = ImageTypeDefinitions::sagittal;
  SetSlice(lastSlice[ImageTypeDefinitions::sagittal]);
  
  ResetCamera();
}


void Render2D3D::SetSliceOrientationToXZ() {
  if (!siv && !tensors) return;
  
  lastSlice[orientation] = GetSlice();
  
  if (siv)
    {
      vtkactor[ImageTypeDefinitions::axial]->SetVisibility(0);
      vtkactor[ImageTypeDefinitions::sagittal]->SetVisibility(0);
      vtkactor[ImageTypeDefinitions::coronal]->SetVisibility(1);
    }
  
  orientation = ImageTypeDefinitions::coronal;
  SetSlice(lastSlice[ImageTypeDefinitions::coronal]);
  
  ResetCamera();
}


void Render2D3D::SetVisibleSlices(bool xy, bool yz, bool xz)
{
  if (!siv && !tensors) return;
  
  if (siv)
    {
      vtkactor[ImageTypeDefinitions::axial]->SetVisibility(xy);
      vtkactor[ImageTypeDefinitions::sagittal]->SetVisibility(yz);
      vtkactor[ImageTypeDefinitions::coronal]->SetVisibility(xz);
    }
}


std::string Render2D3D::GetSliceOrientation()
{
  switch(orientation)
    {
    case ImageTypeDefinitions::axial:
      return "XY";
    case ImageTypeDefinitions::sagittal:
      return "YZ";
    case ImageTypeDefinitions::coronal:
      return "XZ";
    }
  // we will never get here, but the compiler complains otherwise
  return "XY";
}


void Render2D3D::GetSliceRange(int &min, int &max) {
  if (siv)
    {
      int *extent = imageData->GetWholeExtent();
      switch (orientation)
	{
	case ImageTypeDefinitions::axial:
	  min = extent[4];
	  max = extent[5];
	  break;
	case ImageTypeDefinitions::sagittal:
	  min = extent[0];
	  max = extent[1];
	  break;
	case ImageTypeDefinitions::coronal:
	  min = extent[2];
	  max = extent[3];
	}
    }
  else if (tensors)
    {
      min = tensors->GetMinSlice();
      max = tensors->GetMaxSlice();
    }
}


void Render2D3D::SetSlice(int sliceindex) {
  if (siv)
    {
      int min=0,max=0;
      GetSliceRange(min,max);
      if (sliceindex < min || sliceindex > max)
	return;
      
      int *extent = imageData->GetWholeExtent();
      switch (orientation)
	{
	case ImageTypeDefinitions::axial:
	  vtkactor[ImageTypeDefinitions::axial]->SetDisplayExtent(extent[0],extent[1],extent[2],extent[3],sliceindex,sliceindex);
	  break;
	case ImageTypeDefinitions::sagittal:
	  vtkactor[ImageTypeDefinitions::sagittal]->SetDisplayExtent(sliceindex,sliceindex,extent[2],extent[3],extent[4],extent[5]);
	  break;
	case ImageTypeDefinitions::coronal:
	  vtkactor[ImageTypeDefinitions::coronal]->SetDisplayExtent(extent[0],extent[1],sliceindex,sliceindex,extent[4],extent[5]);
	}
    }
  
  if (tensors)
    {
      tensors->SetSlice(sliceindex,orientation);
      AdjustTensorActorPosition();
    }
  
  ren->ResetCameraClippingRange();
  ren->GetRenderWindow()->Render();
}


int Render2D3D::GetSlice() {
  if (!siv) return 0;
  return vtkactor[orientation]->GetSliceNumber();
}


void Render2D3D::On2DRightButtonPressEvent(vtkObject *caller, unsigned long event, void *data, void *calldata)
{
  Render2D3D *me = reinterpret_cast<Render2D3D*>(data);
  
  if (me->interactor->GetShiftKey())
    {
      int *pos = me->interactor->GetEventPosition();
      vtkWorldPointPicker *pick = vtkWorldPointPicker::New();
      pick->Pick(pos[0],pos[1],0,me->ren);
      double world[3],world2[3];
      pick->GetPickPosition(world);
      pick->Delete();
      me->ComputeTrueWorldPosition(world,world2);
      if (me->viewWin)
	me->viewWin->SetMarkerInViewers(me->viewWin->GetGroup(),world2);
      me->Render();
    }
  else
    if (me->cursorVisible && me->interactor->GetControlKey())
      {
	if (me->viewWin)
	  {
	    int newClass = QInputDialog::getInteger(me->viewWin,"Set Label","Please enter the desired label (0 - 255)",me->currentPointClass,0,255);
	    me->currentPointClass = newClass;
	  }
      }
    else
      dynamic_cast<vtkInteractorStyleImage*>(me->interactor->GetInteractorStyle())->OnRightButtonDown();
}


void Render2D3D::On2DLeftButtonPressEvent(vtkObject *caller, unsigned long event, void *data, void *calldata)
{
  Render2D3D *me = reinterpret_cast<Render2D3D*>(data);
  
  if (!me->cursorVisible)
    dynamic_cast<vtkInteractorStyleImage*>(me->interactor->GetInteractorStyle())->OnLeftButtonDown();
  else
    {
      int *pos = me->interactor->GetEventPosition();
      double world[3];
      vtkWorldPointPicker *pick = vtkWorldPointPicker::New();
      pick->Pick(pos[0],pos[1],0,me->ren);
      pick->GetPickPosition(world);
      pick->Delete();
      me->AddOrDeletePoint(world);
    }
}


void Render2D3D::On3DMiddleButtonPressEvent(vtkObject *caller, unsigned long event, void *data, void *calldata)
{
  Render2D3D *me = reinterpret_cast<Render2D3D*>(data);
  if (me->interactor->GetControlKey())
    {
      int *pos = me->interactor->GetEventPosition();
      double x = pos[0];
      double y = pos[1];
      vtkWorldPointPicker *picker = vtkWorldPointPicker::New();
      picker->Pick(x,y,0,me->ren);
      double world[3];
      picker->GetPickPosition(world);
      picker->Delete();
      me->AddOrDeletePoint(world);
    }
  else
    dynamic_cast<vtkInteractorStyleTrackballCamera*>(dynamic_cast<vtkInteractorStyleSwitch*>(me->interactor->GetInteractorStyle())->GetCurrentStyle())->OnRightButtonDown();
}


void Render2D3D::AddOrDeletePoint(double *world)
{
  // this method adds a point at specific world coordinates if no point there exists yet
  // if a point exists already at this position, that point will get deleted
  
  if (!imageData)
    return;
  
  float value,overlay;
  int ijk[3];
  if (!GetVoxelValueAtPosition(world,&value,&overlay,ijk))
    return;
  
  // check if there is already a point with the same label at this position
  // if so, do not add a new point but delete the existing one
  bool deleteIt = false;
  std::vector<PointType>::iterator it,it2;
  it = points.begin();
  it2 = points.end();
  it++; // skip marker
  
  while (it != points.end())
    {
      if ((it->x == ijk[0]) && (it->y == ijk[1]) && (it->z == ijk[2]) && (it->pointClass == currentPointClass))
	{
	  it2 = it;
	  deleteIt = true;
	}
      it++;
    }
  if (!deleteIt)
    {
      PointType point;
      point.x = ijk[0];
      point.y = ijk[1];
      point.z = ijk[2];
      
      point.pixelValue = value;
      point.overlayValue = overlay;
      point.pointClass = currentPointClass;
      point.actor = new CrosshairActor;
      double rgb[3];
      vtkMath::HSVToRGB(((currentPointClass*30)%360)/360.0,1.0,1.0,&rgb[0],&rgb[1],&rgb[2]);
      point.actor->GetActor()->GetProperty()->SetColor(rgb);
      ren->AddViewProp(point.actor->GetActor());
      ComputeTrueWorldPosition(world,point.world);
      
      point.actor->SetWorldPosition(point.world);
      points.push_back(point);
      Render();
      
      std::cout << "Added point to class " << point.pointClass << " at " << point.x << " " << point.y << " " << point.z << " with value " << point.pixelValue << std::endl;
    }
  else
    {
      ren->RemoveViewProp(it2->actor->GetActor());
      delete it2->actor;
      std::cout << "Removed point from class " << it2->pointClass << " at " << it2->x << " " << it2->y << " " << it2->z << " with value " << it2->pixelValue << std::endl;
      points.erase(it2);
      Render();
    }  
}


void Render2D3D::SetCrosshairPosition(double *worldPos)
{
  if (siv)
    {
      points[0].actor->SetWorldPosition(worldPos);
      this->Render();
    }
}


void Render2D3D::ComputeTrueWorldPosition(double *worldPos, double *trueWorldPos)
{
  if (!imageData)
    return;
  
  double bounds[6];
  
  vtkactor[orientation]->GetDisplayBounds(bounds);
  
  switch(orientation)
    {
    case ImageTypeDefinitions::axial:
      trueWorldPos[0] = worldPos[0];
      trueWorldPos[1] = worldPos[1];
      trueWorldPos[2] = bounds[4];
      break;
      
    case ImageTypeDefinitions::sagittal:
      trueWorldPos[0] = bounds[0];
      trueWorldPos[1] = worldPos[1];
      trueWorldPos[2] = worldPos[2];
      break;
      
    case ImageTypeDefinitions::coronal:
      trueWorldPos[0] = worldPos[0];
      trueWorldPos[1] = bounds[2];
      trueWorldPos[2] = worldPos[2];
    }
}


bool Render2D3D::GetVoxelValueAtPosition(double *worldPos, float *imageValue, float *overlayValue, int *ijkPos)
{
  if (!imageData)
    return false;
  
  double world[3];
  int ijk[3];
  double pcoords[3];
  
  pcoords[0] = 0;
  pcoords[1] = 0;
  pcoords[2] = 0;
  
  ComputeTrueWorldPosition(worldPos,world);
  
  if (!imageData->ComputeStructuredCoordinates(world,ijk,pcoords))
    {
      return false;
    }
  else
    {
      *imageValue = imageData->GetScalarComponentAsFloat(ijk[0],ijk[1],ijk[2],0);
      if (overlayImageData)
	*overlayValue = overlayImageData->GetScalarComponentAsFloat(ijk[0],ijk[1],ijk[2],0);
      else
	*overlayValue = 0.0;
    }
  if (ijkPos)
    {
      ijkPos[0] = ijk[0];
      ijkPos[1] = ijk[1];
      ijkPos[2] = ijk[2];
    }
  return true;
}


void Render2D3D::DisplayVoxelValueAtPosition(double *world)
{
  float imageValue, overlayValue;
  int ijk[3];
  if (GetVoxelValueAtPosition(world,&imageValue,&overlayValue,ijk))
    {
      std::stringstream ss;
      ss << "Position: " << ijk[0] << ", " << ijk[1] << ", " << ijk[2] << "        ";
      ss << "Image: " << imageValue;
      ss << "     Overlay: ";
      if (!overlayImageData)
	ss << "N/A";
      else
	ss << overlayValue;
      if (viewWin)
	viewWin->statusBar()->showMessage(QString(ss.str().c_str()));
    }
  else
    if (viewWin)
      viewWin->statusBar()->showMessage("");
}


void Render2D3D::UpdateCrosshairActor()
{
  if (!siv)
    return;
  
  if (!in3D)
    {
      std::vector<PointType>::iterator it;
      for (it = points.begin(); it != points.end(); it++)
	{
	  it->actor->GetActor()->VisibilityOff();
	  double w[3],w2[3];
	  w[0]=0; w[1]=0; w[2]=0;w2[0]=0;w2[1]=0;w2[2]=0;
	  ComputeTrueWorldPosition(w,w2);
	  bool visible = false;
	  switch(orientation)
	    {
	    case ImageTypeDefinitions::axial:
	      if (it->actor->GetWorldPosition()[2] == w2[2]) visible = true;
	      break;
	    case ImageTypeDefinitions::sagittal:
	      if (it->actor->GetWorldPosition()[0] == w2[0]) visible = true;
	      break;
	    case ImageTypeDefinitions::coronal:
	      if (it->actor->GetWorldPosition()[1] == w2[1]) visible = true;
	      break;
	    }
	  
	  if (points[0].actor == it->actor)
	    {
	      // crosshair
	      if (crosshairVisible && visible)
		{
		  points[0].actor->GetActor()->VisibilityOn();
		  
		  if (imageData && !cursorVisible)
		    // if image data is loaded, display values of voxels under crosshair in status bar
		    DisplayVoxelValueAtPosition(points[0].actor->GetWorldPosition());
		}
	      else
		points[0].actor->GetActor()->VisibilityOff();
	    }
	  else
	    {
	      // other points
	      if (visible)
		{
		  it->actor->GetActor()->VisibilityOn();
		}
	      else
		it->actor->GetActor()->VisibilityOff();
	    }
	}
    }
  else
    {
      std::vector<PointType>::iterator it;
      
      if (crosshairVisible)
	points[0].actor->GetActor()->VisibilityOn();
      else
	points[0].actor->GetActor()->VisibilityOff();
      
      // update points
      it = points.begin();
      it++;
      while (it != points.end())
	{
	  it->actor->GetActor()->VisibilityOn();
	  it++;
	}
    }
}


void Render2D3D::OnRenderStartEvent(vtkObject *caller, unsigned long event, void *data, void *calldata)
{
  Render2D3D *me = reinterpret_cast<Render2D3D*>(data);
  
  me->UpdateCrosshairActor();
}


void Render2D3D::On2DMouseMoveEvent(vtkObject *caller, unsigned long event, void *data, void *calldata)
{
  Render2D3D *me = reinterpret_cast<Render2D3D*>(data);
  
  if (me->cursorVisible)
    {
      int *pos = me->interactor->GetEventPosition();
      me->ren->SetDisplayPoint(pos[0],pos[1],me->GetSlice());
      me->ren->DisplayToWorld();
      double *world = me->ren->GetWorldPoint();
      me->cursor->SetFocalPoint(world);
      me->viewWin->CursorPositionHasChanged(world);
      if (me->imageData)
	me->DisplayVoxelValueAtPosition(world);
      
      me->Render();
    }
  
  dynamic_cast<vtkInteractorStyleImage*>(me->interactor->GetInteractorStyle())->OnMouseMove();
}


void Render2D3D::SetBlendValue(double b)
{
  if (!siv) return;
  if ((!overlay_siv) || (overlay_siv->GetDataType() == ImageTypeDefinitions::greyscale) || (siv->GetDataType() != ImageTypeDefinitions::greyscale && overlay_siv->GetDataType() != ImageTypeDefinitions::greyscale) || (!labelBlending))
    {
      vtkblend->SetOpacity(0,1.0-b);
    }
  if (overlay_siv)
    {
      vtkblend->SetOpacity(1,b);
    }
  
  ren->GetRenderWindow()->Render();
}


void Render2D3D::ResetCamera()
{
  vtkCamera *cam = ren->GetActiveCamera();
  
  if (!in3D)
    {
      switch(orientation)
	{
	case ImageTypeDefinitions::axial:
	  cam->SetViewUp(0,-1,0);
	  cam->SetFocalPoint(0,0,0);
	  cam->SetPosition(0,0,-1);
	  break;
	case ImageTypeDefinitions::sagittal:
	  cam->SetViewUp(0,0,1);
	  cam->SetFocalPoint(0,0,0);
	  cam->SetPosition(1,0,0);
	  break;
	case ImageTypeDefinitions::coronal:
	  cam->SetViewUp(0,0,1);
	  cam->SetFocalPoint(0,0,0);
	  cam->SetPosition(0,-1,0);
	}
    }
  
  // reset the camera
  ren->ResetCamera();
  
  double scale = cam->GetParallelScale();
  cam->SetParallelScale( scale * 0.7 );
}


void Render2D3D::GetCameraParameters(double *position, double *focalPoint, double *viewUp)
{
  ren->GetActiveCamera()->GetPosition(position);
  ren->GetActiveCamera()->GetFocalPoint(focalPoint);
  ren->GetActiveCamera()->GetViewUp(viewUp);
}


void Render2D3D::SetCameraParameters(double *position, double *focalPoint, double *viewUp)
{
  vtkCamera *cam = ren->GetActiveCamera();
  cam->SetPosition(position);
  cam->SetFocalPoint(focalPoint);
  cam->SetViewUp(viewUp);
  ren->ResetCamera();
  ren->GetRenderWindow()->Render();
}


void Render2D3D::Set2DInteractor()
{
  if (!ren)
    return;
  
  // set callbacks for 2D mode
  vtkInteractorStyleImage *is = vtkInteractorStyleImage::New();
  qvtk->GetRenderWindow()->GetInteractor()->SetInteractorStyle(is);
  
  vtkCallbackCommand *cbc = vtkCallbackCommand::New();
  cbc->SetCallback(this->On2DRightButtonPressEvent);
  cbc->SetClientData(this);
  is->AddObserver(vtkCommand::RightButtonPressEvent,cbc);
  cbc->Delete();
  cbc = vtkCallbackCommand::New();
  cbc->SetCallback(this->On2DLeftButtonPressEvent);
  cbc->SetClientData(this);
  is->AddObserver(vtkCommand::LeftButtonPressEvent,cbc);
  cbc->Delete();
  cbc = vtkCallbackCommand::New();
  cbc->SetCallback(this->OnRenderStartEvent);
  cbc->SetClientData(this);
  ren->AddObserver(vtkCommand::StartEvent,cbc);
  cbc->Delete();
  cbc = vtkCallbackCommand::New();
  cbc->SetCallback(this->On2DMouseMoveEvent);
  cbc->SetClientData(this);
  is->AddObserver(vtkCommand::MouseMoveEvent,cbc);
  cbc->Delete();
  cbc = vtkCallbackCommand::New();
  cbc->SetCallback(this->WindowLevelCallback);
  cbc->SetClientData(this);
  is->AddObserver(vtkCommand::WindowLevelEvent,cbc);
  is->AddObserver(vtkCommand::StartWindowLevelEvent,cbc);
  is->AddObserver(vtkCommand::ResetWindowLevelEvent,cbc);
  cbc->Delete();
  is->Delete();
  
  in3D = false;
  if (outlineActor)
    {
      ren->RemoveViewProp(outlineActor);
      orientationMarkers->SetEnabled(0);
    }
  
  ResetCamera();
  ren->GetRenderWindow()->Render();
}


void Render2D3D::Set3DInteractor()
{
  if (!ren)
    return;
  
  // turn off cursor when switching to 3D mode
  //if (cursorVisible)
  //SetCursor(false);
  
  vtkInteractorStyleSwitch *is = vtkInteractorStyleSwitch::New();
  is->SetCurrentStyleToTrackballCamera();
  //is->SetCurrentStyleToTrackballActor();
  qvtk->GetRenderWindow()->GetInteractor()->SetInteractorStyle(is);
  is->Delete();
  in3D = true;
  
  vtkCallbackCommand *cbc = vtkCallbackCommand::New();
  cbc->SetCallback(this->On3DRightButtonPressEvent);
  cbc->SetClientData(this);
  is->GetCurrentStyle()->AddObserver(vtkCommand::RightButtonPressEvent,cbc);
  cbc->Delete();
  cbc = vtkCallbackCommand::New();
  cbc->SetCallback(this->On3DMiddleButtonPressEvent);
  cbc->SetClientData(this);
  is->GetCurrentStyle()->AddObserver(vtkCommand::MiddleButtonPressEvent,cbc);
  cbc->Delete();
  
  if (!imageData)
    return;
  // create outline and annotation cube if necessary
  if (!outlineActor)
    {
      vtkOutlineFilter *outline = vtkOutlineFilter::New();
      outline->SetInput(imageData);
      vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
      mapper->SetInput(outline->GetOutput());
      outlineActor = vtkActor::New();
      outlineActor->SetMapper(mapper);
      outlineActor->GetProperty()->SetColor(1,0,0);
      mapper->Delete();
      outline->Delete();
      cube = vtkAnnotatedCubeActor::New();
      orientationMarkers = vtkOrientationMarkerWidget::New();
      cube->SetXPlusFaceText("L");
      cube->SetXMinusFaceText("R");
      cube->SetYPlusFaceText("P");
      cube->SetYMinusFaceText("A");
      cube->SetZPlusFaceText("S");
      cube->SetZMinusFaceText("I");
      orientationMarkers->SetOrientationMarker(cube);
      orientationMarkers->SetInteractor(qvtk->GetRenderWindow()->GetInteractor());
    }
  ren->AddViewProp(outlineActor);
  orientationMarkers->SetEnabled(1);
  orientationMarkers->SetInteractive(0);
  points[0].actor->GetActor()->VisibilityOff();
  
  ResetCamera();
  ren->GetRenderWindow()->Render();
}


bool Render2D3D::Using3DInteractor()
{
  return in3D;
}


bool Render2D3D::SetTensorInput(SceneTensors *st)
{
  if (tensorActor)
    {
      ren->RemoveViewProp(tensorActor);
      tensorActor->Delete();
      tensorActor = 0;
    }
  
  tensors = st;
  if (tensors)
    {
      tensorActor = vtkActor::New();
      vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
      if (siv)
	mapper->SetInput(tensors->GetVTKTensorGlyphs(siv));
      else
	mapper->SetInput(tensors->GetVTKTensorGlyphs());
      mapper->Update();
      tensorActor->SetMapper(mapper);
      tensorActor->GetProperty()->SetAmbient(1.0);
      tensorActor->GetProperty()->SetDiffuse(0.0);
      ren->AddViewProp(tensorActor);
      mapper->Delete();
    }
  
  if (vtkactor[orientation])
    tensors->SetSlice(this->GetSlice(),this->orientation);
  
  //AdjustTensorActorPosition();
  ResetCamera();
  
  return true;
}


bool Render2D3D::AddModel(SceneModel *sm)
{
  if (modelActors.find(sm) == modelActors.end())
    {
      // std::cout << "Adding model " << sm->GetName().toStdString() << std::endl;
      vtkActor *actor;
      if (siv)
	actor = sm->CreateNewActor(siv);
      else
	actor = sm->CreateNewActor();
      
      ren->AddViewProp(actor);
      modelActors[sm].push_back(actor);
    }
  else
    {
      // std::cout << "Removing model " << sm->GetName().toStdString() << std::endl;
      std::vector<vtkActor*>::iterator it;
      for (it = modelActors[sm].begin(); it != modelActors[sm].end(); it++)
	ren->RemoveViewProp(*it);
      modelActors.erase(sm);
    }
  return true;
}


void Render2D3D::SetDataType(int number,ImageTypeDefinitions::DataType t)
{
  if (number == 0 || number == 1)
    {
      scalarDataType[number] = t;
      
      if ((number == 1) || (!vtkmapper))
	return;
      
      siv->SetDataType(t);
      vtkLookupTable *vtklut = siv->GetLookupTable();
      vtkmapper->SetLookupTable(vtklut);
      Render();
    }
}


void Render2D3D::SetCursor(bool status)
{
  if (!siv)
    return;
  
  cursorVisible = status;
  
  if (status)
    {
      ren->AddViewProp(cursorActor);
    }
  else
    {
      ren->RemoveViewProp(cursorActor);
      // reset status bar
      if (viewWin)
	viewWin->statusBar()->showMessage("");
      UpdateCrosshairActor();
    }
  ren->ResetCameraClippingRange();
  Render();
}


void Render2D3D::SetCursorPosition(double pos[3])
{
  if (siv && cursorVisible)
    {
      cursor->SetFocalPoint(pos);
      if (imageData)
	DisplayVoxelValueAtPosition(pos);
      Render();
    }
}


void Render2D3D::AdjustTensorActorPosition()
{
  // set correct actor position
  if (vtkactor[orientation])
    {
      double bounds[6];
      vtkactor[orientation]->GetDisplayBounds(bounds);
      tensorActor->SetPosition(bounds[0],bounds[2],bounds[4]);
    }
}


void Render2D3D::WindowLevelCallback(vtkObject *caller,unsigned long event,void *data,void *calldata)
{
  Render2D3D *me = reinterpret_cast<Render2D3D*>(data);
  
  if (!me->siv)
    return;
  
  // Reset
  
  if (event == vtkCommand::ResetWindowLevelEvent)
    {
      if (me->overlay_siv && (me->vtkblend->GetOpacity(1) == 1.0))
	{
	  if (me->scalarDataType[1] == ImageTypeDefinitions::greyscale)
	    {
	      dynamic_cast<vtkImageMapToWindowLevelColors*>(me->vtkoverlaymapper)->SetWindow(me->overlay_siv->GetInitialWindow());
	      dynamic_cast<vtkImageMapToWindowLevelColors*>(me->vtkoverlaymapper)->SetLevel(me->overlay_siv->GetInitialLevel());
	    }
	}
      else
	{
	  if (me->scalarDataType[0] == ImageTypeDefinitions::greyscale)
	    {
	      dynamic_cast<vtkImageMapToWindowLevelColors*>(me->vtkmapper)->SetWindow(me->siv->GetInitialWindow());
	      dynamic_cast<vtkImageMapToWindowLevelColors*>(me->vtkmapper)->SetLevel(me->siv->GetInitialLevel());
	    }
	}
      me->ren->GetRenderWindow()->Render();
      return;
    }
  
  // Start
  
  if (event == vtkCommand::StartWindowLevelEvent)
    {
      if (me->overlay_siv && (me->vtkblend->GetOpacity(1) == 1.0))
	{
	  if (me->scalarDataType[1] == ImageTypeDefinitions::greyscale)
	    {
	      me->wlInitialWindowOverlay = dynamic_cast<vtkImageMapToWindowLevelColors*>(me->vtkoverlaymapper)->GetWindow();
	      me->wlInitialLevelOverlay = dynamic_cast<vtkImageMapToWindowLevelColors*>(me->vtkoverlaymapper)->GetLevel();
	    }
	}
      else
	{
	  if (me->scalarDataType[0] == ImageTypeDefinitions::greyscale)
	    {
	      me->wlInitialWindow = dynamic_cast<vtkImageMapToWindowLevelColors*>(me->vtkmapper)->GetWindow();
	      me->wlInitialLevel = dynamic_cast<vtkImageMapToWindowLevelColors*>(me->vtkmapper)->GetLevel();
	    }
	}
      return;
    }
  
  // Adjust the window level here
  
  vtkInteractorStyleImage *isi = dynamic_cast<vtkInteractorStyleImage*>(caller);
  
  int *size = me->ren->GetRenderWindow()->GetSize();
  
  // Compute normalized delta
  
  double dx = 4.0 *
    (isi->GetWindowLevelCurrentPosition()[0] -
     isi->GetWindowLevelStartPosition()[0]) / size[0];
  double dy = 4.0 *
    (isi->GetWindowLevelStartPosition()[1] -
     isi->GetWindowLevelCurrentPosition()[1]) / size[1];
  
  double odx = dx;
  double ody = dy;
  
  if (me->overlay_siv && (me->vtkblend->GetOpacity(1) == 1.0))
    {
      if (fabs(me->wlInitialWindowOverlay) > 0.01)
	{
	  odx = odx * me->wlInitialWindowOverlay;
	}
      else
	{
	  odx = odx * (me->wlInitialWindowOverlay < 0 ? -0.01 : 0.01);
	}
      if (fabs(me->wlInitialLevelOverlay) > 0.01)
	{
	  ody = ody * me->wlInitialLevelOverlay;
	}
      else
	{
	  ody = ody * (me->wlInitialLevelOverlay < 0 ? -0.01 : 0.01);
	}
      
      // Abs so that direction does not flip
      
      if (me->wlInitialWindowOverlay < 0.0)
	{
	  odx = -1*odx;
	}
      if (me->wlInitialLevelOverlay < 0.0)
	{
	  ody = -1*ody;
	}
    }
  else
    {
      // Scale by current values
      
      if (fabs(me->wlInitialWindow) > 0.01)
	{
	  dx = dx * me->wlInitialWindow;
	}
      else
	{
	  dx = dx * (me->wlInitialWindow < 0 ? -0.01 : 0.01);
	}
      if (fabs(me->wlInitialLevel) > 0.01)
	{
	  dy = dy * me->wlInitialLevel;
	}
      else
	{
	  dy = dy * (me->wlInitialLevel < 0 ? -0.01 : 0.01);
	}
      
      // Abs so that direction does not flip
      
      if (me->wlInitialWindow < 0.0)
	{
	  dx = -1*dx;
	}
      if (me->wlInitialLevel < 0.0)
	{
	  dy = -1*dy;
	}
    }
  
  if (me->overlay_siv && (me->vtkblend->GetOpacity(1) == 1.0))
    {
      if (me->scalarDataType[1] == ImageTypeDefinitions::greyscale)
	{
	  dynamic_cast<vtkImageMapToWindowLevelColors*>(me->vtkoverlaymapper)->SetWindow(odx + me->wlInitialWindowOverlay);
	  dynamic_cast<vtkImageMapToWindowLevelColors*>(me->vtkoverlaymapper)->SetLevel(me->wlInitialLevelOverlay - ody);
	}
    }
  else
    {
      if (me->scalarDataType[0] == ImageTypeDefinitions::greyscale)
	{
	  dynamic_cast<vtkImageMapToWindowLevelColors*>(me->vtkmapper)->SetWindow(dx + me->wlInitialWindow);
	  dynamic_cast<vtkImageMapToWindowLevelColors*>(me->vtkmapper)->SetLevel(me->wlInitialLevel - dy);
	}
    }
  me->ren->GetRenderWindow()->Render();
  me->InvokeEvent(Render2D3DWindowLevelEvent);
}


void Render2D3D::SavePoints(std::string filename)
{
  if (!imageData)
    return;
  
  std::vector<PointType>::iterator it;
  
  std::ofstream pointsFile(filename.c_str());
  
  it = points.begin();
  it++; // don't save the crosshair marker
  
  while (it != points.end())
    {
      pointsFile << it->pointClass << std::endl;
      pointsFile << it->pixelValue << " " << it->overlayValue << " ";
      pointsFile << ":row-" << it->y << "-col-" << it->x << "-sliceindex-" << it->z << std::endl;
      it++;
    }
  pointsFile.close();
}


void Render2D3D::LoadPoints(std::string filename)
{
  if (!imageData)
    return;
  
  std::ifstream pointsFile(filename.c_str());
  std::string line;
  
  PointType point;
  unsigned int count = 0;
  
  while (!pointsFile.eof())
    {
      getline(pointsFile,line);
      if (line.length() < 10)
	sscanf(line.c_str(),"%d",&point.pointClass);
      else
	{
	  float pixelValue, overlayValue;
	  sscanf(line.c_str(),"%f %f :row-%d-col-%d-sliceindex-%d",&pixelValue,&overlayValue,&point.y,&point.x,&point.z);
	  point.pixelValue = pixelValue;
	  point.overlayValue = overlayValue;
	  point.actor = new CrosshairActor;
	  double rgb[3];
	  vtkMath::HSVToRGB(((point.pointClass*30)%360)/360.0,1.0,1.0,&rgb[0],&rgb[1],&rgb[2]);
	  point.actor->GetActor()->GetProperty()->SetColor(rgb);
	  ren->AddViewProp(point.actor->GetActor());
	  
	  point.world[0] = imageData->GetOrigin()[0]+imageData->GetSpacing()[0]*point.x;
	  point.world[1] = imageData->GetOrigin()[1]+imageData->GetSpacing()[1]*point.y;
	  point.world[2] = imageData->GetOrigin()[2]+imageData->GetSpacing()[2]*point.z;
	  point.actor->SetWorldPosition(point.world);
	  points.push_back(point);
	  count++;
	}
    }
  pointsFile.close();
  std::cout << count << " points added." << std::endl;
}


void Render2D3D::ClearPoints()
{
  if (!imageData)
    return;
  
  std::vector<PointType>::iterator it;
  it = points.begin();
  it++; // don't delete the crosshair
  while (it != points.end())
    {
      ren->RemoveViewProp(it->actor->GetActor());
      delete it->actor;
      it++;
    }
  PointType crosshair = points[0];
  points.clear();
  points.push_back(crosshair);
}


void Render2D3D::On3DRightButtonPressEvent(vtkObject *caller, unsigned long event, void *data, void *calldata)
{
  Render2D3D *me = reinterpret_cast<Render2D3D*>(data);
  
  
  if (me->interactor->GetShiftKey() && (!me->interactor->GetControlKey()))
    {
      bool picked = false;
      
      int *pos = me->interactor->GetEventPosition();
      double x = pos[0];
      double y = pos[1];
      vtkPropPicker *picker = vtkPropPicker::New();
      if (picker->PickProp(x,y,me->ren))
	{
	  // an object was picked, let's find out if it is a model
	  vtkProp3D *prop = dynamic_cast<vtkProp3D*>(picker->GetViewProp());
	  vtkActor *actor = dynamic_cast<vtkActor*>(prop);
	  if (actor)
	    {
	      std::map<SceneModel*,std::vector<vtkActor*> >::iterator it;
	      std::vector<vtkActor*>::iterator it2,it3;
	      for (it = me->modelActors.begin(); it != me->modelActors.end(); it++)
		for (it2 = it->second.begin(); it2 != it->second.end(); it2++)
		  if (*it2 == actor)
		    {
		      if (me->viewWin)
			{
			  // change the opacity for this model
			  double opacity = QInputDialog::getDouble(me->viewWin,"Set Opacity","Please enter the desired opacity (0.0 - 1.0)",actor->GetProperty()->GetOpacity(),0.0,1.0,1);
			  //actor->GetProperty()->SetOpacity(opacity);
			  for (it3 = it->second.begin(); it3 != it->second.end(); it3++)
			    (*it3)->GetProperty()->SetOpacity(opacity);
			  picked = true;
			}
		    }
	    }
	  else
	    {
	      if (prop == me->vtkactor[ImageTypeDefinitions::axial])
		{
		  picked = true;
		}
	      if (prop == me->vtkactor[ImageTypeDefinitions::sagittal])
		{
		  picked = true;
		}
	      if (prop == me->vtkactor[ImageTypeDefinitions::coronal])
		{
		  picked = true;
		}
	      if (picked)
		{
		  double opacity = QInputDialog::getDouble(me->viewWin,"Set Opacity","Please enter the desired opacity (0.0 - 1.0)",dynamic_cast<vtkImageActor*>(prop)->GetOpacity(),0.0,1.0,1);
		  dynamic_cast<vtkImageActor*>(prop)->SetOpacity(opacity);
		}
	    }
	}
      
      picker->Delete();
      
      if (!picked)
	// if nothing was actually picked, forward the event to the interactor style for zooming
	dynamic_cast<vtkInteractorStyleTrackballCamera*>(dynamic_cast<vtkInteractorStyleSwitch*>(me->interactor->GetInteractorStyle())->GetCurrentStyle())->OnRightButtonDown();
    }
  else if (me->interactor->GetControlKey() && (!me->interactor->GetShiftKey()))
    {
      bool picked = false;
      int *pos = me->interactor->GetEventPosition();
      double x = pos[0];
      double y = pos[1];
      vtkCellPicker *picker = vtkCellPicker::New();
      if (picker->Pick(x,y,0,me->ren))
	{
	  if (picker->GetCellId() != -1)
	    {
	      vtkActor *actor = dynamic_cast<vtkActor*>(picker->GetViewProp());
	      for (std::map<SceneModel*,std::vector<vtkActor*> >::iterator it = me->modelActors.begin(); it != me->modelActors.end(); it++)
		for (std::vector<vtkActor*>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
		  if (*it2 == actor)
		    {
		      vtkPolyDataConnectivityFilter *pdcf = vtkPolyDataConnectivityFilter::New();
		      pdcf->SetInput(dynamic_cast<vtkPolyData*>(actor->GetMapper()->GetInput()));
		      pdcf->SetExtractionModeToCellSeededRegions();
		      pdcf->AddSeed(picker->GetCellId());
		      pdcf->ColorRegionsOn();
		      pdcf->Update();
		      vtkPolyDataMapper *pdcfm = vtkPolyDataMapper::New();
		      pdcfm->SetInput(pdcf->GetOutput());
		      vtkLookupTable *pdcflut = vtkLookupTable::New();
		      pdcflut->SetTableValue(0,me->currentColor[0],me->currentColor[1],me->currentColor[2]);
		      pdcflut->Build();
		      pdcfm->SetLookupTable(pdcflut);
		      pdcflut->Delete();
		      vtkActor *pdcfa = vtkActor::New();
		      pdcfa->SetMapper(pdcfm);
		      me->ren->AddViewProp(pdcfa);
		      it->second.push_back(pdcfa);
		      pdcfm->Delete();
		      
		      picked = true;
		      me->Render();
		      break;
		    }
	    }
	}
      if (!picked)
	dynamic_cast<vtkInteractorStyleTrackballCamera*>(dynamic_cast<vtkInteractorStyleSwitch*>(me->interactor->GetInteractorStyle())->GetCurrentStyle())->OnRightButtonDown();
    }
  else if (me->interactor->GetControlKey() && me->interactor->GetShiftKey())
    {
      QColor newColor = QColorDialog::getColor(QColor((int)floor(me->currentColor[0]*255),(int)floor(me->currentColor[1]*255),(int)floor(me->currentColor[2]*255)),me->viewWin);
      if (newColor.isValid())
	{
	  me->currentColor[0] = newColor.redF();
	  me->currentColor[1] = newColor.greenF();
	  me->currentColor[2] = newColor.blueF();
	}
    }
  else
    dynamic_cast<vtkInteractorStyleTrackballCamera*>(dynamic_cast<vtkInteractorStyleSwitch*>(me->interactor->GetInteractorStyle())->GetCurrentStyle())->OnRightButtonDown();
}


void Render2D3D::GetWindowLevel(double wl[4])
{
  if (!this->siv)
    return;
  
  if (this->scalarDataType[0] == ImageTypeDefinitions::greyscale)
    {
      wl[0] = dynamic_cast<vtkImageMapToWindowLevelColors*>(this->vtkmapper)->GetWindow();
      wl[1] = dynamic_cast<vtkImageMapToWindowLevelColors*>(this->vtkmapper)->GetLevel();
    }
  else
    {
      wl[0] = 0.0;
      wl[1] = 0.0;
    }
  
  if (this->overlay_siv && (this->scalarDataType[1] == ImageTypeDefinitions::greyscale))
    {
      wl[2] = dynamic_cast<vtkImageMapToWindowLevelColors*>(this->vtkoverlaymapper)->GetWindow();
      wl[3] = dynamic_cast<vtkImageMapToWindowLevelColors*>(this->vtkoverlaymapper)->GetLevel();
    }
  else
    {
      wl[2] = 0.0;
      wl[3] = 0.0;
    }
}


void Render2D3D::SetWindowLevel(double wl[4])
{
  if (!this->siv)
    return;
  
  if (this->scalarDataType[0] == ImageTypeDefinitions::greyscale)
    {
      dynamic_cast<vtkImageMapToWindowLevelColors*>(this->vtkmapper)->SetWindow(wl[0]);
      dynamic_cast<vtkImageMapToWindowLevelColors*>(this->vtkmapper)->SetLevel(wl[1]);
    }
  
  if (this->overlay_siv && (this->scalarDataType[1] == ImageTypeDefinitions::greyscale))
    {
      dynamic_cast<vtkImageMapToWindowLevelColors*>(this->vtkoverlaymapper)->SetWindow(wl[2]);
      dynamic_cast<vtkImageMapToWindowLevelColors*>(this->vtkoverlaymapper)->SetLevel(wl[3]);
    }
}


void Render2D3D::EnableROIBox()
{
  if (!imageData)
    return;
  
  roiBox = vtkBoxWidget::New();
  roiBox->PlaceWidget(imageData->GetBounds());
  roiBox->RotationEnabledOff();
  roiBox->TranslationEnabledOn();
  roiBox->ScalingEnabledOn();
  roiBox->OutlineFaceWiresOn();
  roiBox->SetInteractor(qvtk->GetRenderWindow()->GetInteractor());
  roiBox->SetEnabled(1);
}


void Render2D3D::DisableROIBox()
{
  if (roiBox)
    {
      roiBox->SetEnabled(0);
      roiBox->SetInteractor(NULL);
      roiBox->Delete();
      roiBox = 0;
    }
}


void Render2D3D::SaveROI(std::string filename)
{
  if (!roiBox)
    return;
  
  vtkPolyData *points = vtkPolyData::New();
  roiBox->GetPolyData(points);
  int ijk[6][3];
  double pc[3] = {0,0,0};
  int bounds[6];
  
  imageData->ComputeStructuredCoordinates(points->GetPoint(8),ijk[0],pc);
  imageData->ComputeStructuredCoordinates(points->GetPoint(9),ijk[1],pc);
  imageData->ComputeStructuredCoordinates(points->GetPoint(10),ijk[2],pc);
  imageData->ComputeStructuredCoordinates(points->GetPoint(11),ijk[3],pc);
  imageData->ComputeStructuredCoordinates(points->GetPoint(12),ijk[4],pc);
  imageData->ComputeStructuredCoordinates(points->GetPoint(13),ijk[5],pc);
  
  bounds[0] = ijk[0][0];
  bounds[1] = ijk[1][0]+1;
  //int *dimensions = imageData->GetDimensions();
  // x axis is flipped
  //bounds[0] = dimensions[0]-ijk[1][0];
  //bounds[1] = dimensions[0]-ijk[0][0];
  bounds[2] = ijk[2][1];
  bounds[3] = ijk[3][1]+1;
  // y axis is flipped
  //bounds[2] = dimensions[1]-ijk[3][1];
  //bounds[3] = dimensions[1]-ijk[2][1];
  bounds[4] = ijk[4][2];
  bounds[5] = ijk[5][2]+1;
  
  ImageTypeDefinitions::WriterType::Pointer writer = ImageTypeDefinitions::WriterType::New();
  ImageTypeDefinitions::ImageType::Pointer roi = siv->ExtractROI(bounds);
  
  std::stringstream output;
  output << bounds[0] << " " << bounds[2] << " " << bounds[4] << " " << bounds[1]-bounds[0] << " " << bounds[3]-bounds[2] << " " << bounds[5]-bounds[4];
  
  std::cout << output.str() << std::endl;
  if (!roi)
    {
      points->Delete();
      return;
    }
  
  writer->SetInput(roi);
  writer->SetFileName(filename.c_str());
  writer->UseCompressionOn();
  
  itk::MetaDataDictionary &mdd = roi->GetMetaDataDictionary();
  itk::EncapsulateMetaData<std::string>(mdd,"CRKIT_ROI",output.str());
  try
    {
      writer->Write();
    }
  catch (itk::ExceptionObject &e)
    {
      std::cout << "Writing of ROI failed." << std::endl;
      std::cout << e << std::endl;
    }
  
  points->Delete();
}


void Render2D3D::SaveScreenshot(std::string filename)
{
  vtkWindowToImageFilter *wtif = vtkWindowToImageFilter::New();
  wtif->SetInputBufferTypeToRGBA();
  wtif->SetInput(ren->GetRenderWindow());
  wtif->Modified();
  
  vtkPNGWriter *writer = vtkPNGWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInput(wtif->GetOutput());
  writer->Write();
  
  writer->Delete();
  wtif->Delete();
}


void Render2D3D::SaveOBJ(std::string filename)
{
  if (ren)
    {
      if (orientationMarkers) orientationMarkers->SetEnabled(0);
      vtkOBJExporter *exporter = vtkOBJExporter::New();
      exporter->SetInput(ren->GetRenderWindow());
      exporter->SetFilePrefix(filename.c_str());
      exporter->Write();
      exporter->Delete();
      if (orientationMarkers) orientationMarkers->SetEnabled(1);
    }
}


void Render2D3D::SaveRIB(std::string filename)
{
  if (ren)
    {
      if (orientationMarkers) orientationMarkers->SetEnabled(0);
      vtkRIBExporter *exporter = vtkRIBExporter::New();
      exporter->SetInput(ren->GetRenderWindow());
      exporter->SetFilePrefix(filename.c_str());
      exporter->Write();
      exporter->Delete();
      if (orientationMarkers) orientationMarkers->SetEnabled(1);
    }
}


void Render2D3D::ShowMarker(bool status)
{
  if (!siv)
    return;
  
  crosshairVisible = status;
  points[0].actor->GetActor()->SetVisibility(status);
}


void Render2D3D::ShowBoundingBox(bool status)
{
  if (in3D && outlineActor)
    {
      outlineActor->SetVisibility(status);
      Render();
    }
}


void Render2D3D::ShowOrientationCube(double size)
{
  if (in3D && cube)
    {
      if (size == 0.0)
	cube->SetVisibility(false);
      else if (size == 1.0)
	cube->SetVisibility(true);
      
      if (size > 0.0)
	cube->SetScale(size);
      
      Render();
    }
}


void Render2D3D::SetBackground(double r,double g,double b)
{
  if (ren)
    ren->SetBackground(r,g,b);
  Render();
}


void Render2D3D::FixDepthSort()
{
  /*  std::map<SceneModel*,vtkActor*>::iterator it;
      for (it = modelActors.begin(); it != modelActors.end(); it++)
      {
      // check if there is already a sorter
      vtkDepthSortPolyData *sorter = dynamic_cast<vtkDepthSortPolyData*>(it->second->GetMapper()->GetInput());
      if (sorter)
      {
      dynamic_cast<vtkPolyDataMapper*>(it->second->GetMapper())->SetInput(it->first->GetPolyData());
      sorter->Delete();
      }
      
      sorter = vtkDepthSortPolyData::New();
      sorter->SetInput(it->first->GetPolyData());
      sorter->SetCamera(ren->GetActiveCamera());
      sorter->Update();
      
      dynamic_cast<vtkPolyDataMapper*>(it->second->GetMapper())->SetInput(sorter->GetOutput());
      }*/
}

void Render2D3D::ToggleFullScreen()
{
  if (ren)
    {
      ren->GetRenderWindow()->SetFullScreen(!(ren->GetRenderWindow()->GetFullScreen()));
      Render();
    }
}


unsigned int Render2D3D::GetNumberOfRenderersForObject(SceneData *sd)
{
  unsigned int count = 0;
  
  if (dynamic_cast<SceneImage*>(sd) || dynamic_cast<SceneRGBImage*>(sd))
    {
      if (siv->GetData() == sd) count++;
      if (overlay_siv->GetData() == sd) count++;
    }
  else if (dynamic_cast<SceneModel*>(sd))
    {
      if (modelActors.find(dynamic_cast<SceneModel*>(sd)) != modelActors.end())
	count++;
    }
  else if (dynamic_cast<SceneTensors*>(sd))
    {
      if (tensors == sd) count++;
    }
  return count;
}


void Render2D3D::SaveModel(SceneData *model,std::string fileName)
{
  vtkAppendPolyData *sum = vtkAppendPolyData::New();
  bool hasScalars = false;
  
  for (std::vector<vtkActor*>::iterator it = modelActors[dynamic_cast<SceneModel*>(model)].begin(); it != modelActors[dynamic_cast<SceneModel*>(model)].end(); it++)
    {
      vtkPolyData *pd = dynamic_cast<vtkPolyDataMapper*>((*it)->GetMapper())->GetInput();
      if (!pd->GetPointData()->GetScalars())
	{
	  // add gray scalars
	  vtkUnsignedCharArray *scalars = vtkUnsignedCharArray::New();
	  scalars->SetNumberOfComponents(3);
	  
	  unsigned char grey[3] = {192,192,192};
	  for (vtkIdType id = 0; id < pd->GetNumberOfPoints(); id++)
	    scalars->InsertNextTupleValue(grey);
	  pd->GetPointData()->SetScalars(scalars);
	  scalars->Delete();
	}
      else
	{
	  hasScalars = true;
	  if (pd->GetPointData()->GetScalars()->GetNumberOfComponents() == 1)
	    {
	      vtkIdTypeArray *scalars = dynamic_cast<vtkIdTypeArray*>(pd->GetPointData()->GetScalars());
	      if (scalars)
		{
		  vtkUnsignedCharArray *newScalars = vtkUnsignedCharArray::New();
		  newScalars->SetNumberOfComponents(3);
		  unsigned char color[3];
		  double dcolor[4];
		  for (vtkIdType id = 0; id < pd->GetNumberOfPoints(); id++)
		    {
		      dynamic_cast<vtkLookupTable*>((*it)->GetMapper()->GetLookupTable())->GetTableValue(scalars->GetValue(id),dcolor);
		      color[0]=(unsigned char)(dcolor[0]*255);
		      color[1]=(unsigned char)(dcolor[1]*255);
		      color[2]=(unsigned char)(dcolor[2]*255);
		      newScalars->InsertNextTupleValue(color);
		    }
		  pd->GetPointData()->SetScalars(newScalars);
		  newScalars->Delete();
		}
	    }
	}
      sum->AddInput(pd);
    }
  
  sum->Update();
  //vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  vtkGIFTIWriter *writer = vtkGIFTIWriter::New();
  writer->SetFileName(fileName.c_str());
  //writer->SetFileTypeToBinary();
  writer->SetInput(sum->GetOutput());
  writer->Write();
}


void Render2D3D::SaveMovie(std::string fileName)
{
#ifndef WIN32
  timespec delay = {0,100000000};
#endif
  std::stringstream name;
  for (int i = 0; i < 36 ; i++)
    {
      name.str("");
      name << fileName << std::setfill('0') << std::setw(3) << i+1 << ".png";

      vtkWindowToImageFilter *wtif = vtkWindowToImageFilter::New();
      wtif->SetInputBufferTypeToRGBA();
      wtif->SetInput(ren->GetRenderWindow());
      wtif->Modified();
      vtkPNGWriter *writer = vtkPNGWriter::New();
      writer->SetFileName(name.str().c_str());
      writer->SetInput(wtif->GetOutput());
      writer->Write();
      wtif->Delete(); 
      // Rotate camera, render and sleep for 100 ms
      ren->GetActiveCamera()->Azimuth(10.0);
      Render();
#ifndef WIN32
      nanosleep(&delay,NULL);
#else
      Sleep(100);
#endif
    }
}

void Render2D3D::ShowMousePointer(bool status)
{
  mousePointerVisible = status;
  
  if (mousePointerVisible)
    ren->GetRenderWindow()->ShowCursor();
  else
    ren->GetRenderWindow()->HideCursor();
}

void Render2D3D::SetSpecialLabelBlending(bool status)
{
  labelBlending = status;
}


void Render2D3D::SetToolPosition(vtkMatrix4x4 *position, bool isValid)
{
  position->SetElement(0,0,(position->Element[0][0])*50);
  position->SetElement(1,0,(position->Element[1][0])*50);
  position->SetElement(2,0,(position->Element[2][0])*50);
  igtTool->SetUserMatrix(position);
  
  if (isValid)
    igtTool->GetProperty()->SetColor(0,1,0);
  else
    igtTool->GetProperty()->SetColor(1,0,0);
  Render();
}

void Render2D3D::ShowIGTTool(bool status)
{
  if (status)
    ren->AddViewProp(igtTool);
  else
    ren->RemoveViewProp(igtTool);
  
  Render();
}

std::set<SceneModel*> Render2D3D::GetModels()
{
  std::set<SceneModel*> models;
  
  std::map<SceneModel*,std::vector<vtkActor*> >::iterator it;
  for (it = modelActors.begin(); it != modelActors.end(); it++)
    {
      models.insert(it->first);
    }
  
  return models;
}
