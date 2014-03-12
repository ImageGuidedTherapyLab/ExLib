#include "VizViewerWindowRotationDialog.h"

#include <SceneModel.h>
#include <Render2D3D.h>

#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyData.h>

#include <itkAffineTransform.h>
#include <itkTransformFileWriter.h>

VizViewerWindowRotationDialog::VizViewerWindowRotationDialog()
{
  setupUi(this);
  
  ren = 0;
  model = 0;
  filter = vtkTransformPolyDataFilter::New();
  transform = vtkTransform::New();
  transform->PostMultiply();
  original = vtkPolyData::New();
  filter->SetTransform(transform);
  
  connect(sliderX,SIGNAL(valueChanged(int)),this,SLOT(X_value_changed(int)));
  connect(sliderY,SIGNAL(valueChanged(int)),this,SLOT(Y_value_changed(int)));
  connect(sliderZ,SIGNAL(valueChanged(int)),this,SLOT(Z_value_changed(int)));
  connect(save,SIGNAL(clicked()),this,SLOT(save_clicked()));
  
  this->setAttribute(Qt::WA_DeleteOnClose);
}


VizViewerWindowRotationDialog::~VizViewerWindowRotationDialog()
{
  filter->Delete();
  transform->Delete();
  original->Delete();
}


void VizViewerWindowRotationDialog::SetModel(SceneModel *m)
{
  model = m;
  original->DeepCopy(model->GetPolyData());
  filter->SetInput(original);
  filter->Update();
  
  // compute centroid
  centroid[0] = 0.0;
  centroid[1] = 0.0;
  centroid[2] = 0.0;
  
  vtkIdType numpoints = original->GetPoints()->GetNumberOfPoints();
  for (vtkIdType i = 0; i < numpoints; i++)
    {
      double point[3];
      original->GetPoints()->GetPoint(i,point);
      centroid[0] += point[0];
      centroid[1] += point[1];
      centroid[2] += point[2];
    }
  
  centroid[0] /= numpoints;
  centroid[1] /= numpoints;
  centroid[2] /= numpoints;
}


void VizViewerWindowRotationDialog::SetRenderer(Render2D3D *r)
{
  ren = r;
}


void VizViewerWindowRotationDialog::X_value_changed(int x)
{
  if (model && ren)
    {
      labelX->setText(QString().setNum(x));
      transform->Identity();
      transform->Translate(-centroid[0],-centroid[1],-centroid[2]);
      transform->RotateX(x);
      transform->Translate(centroid);
      filter->Update();
      model->SetPolyData(filter->GetOutput());
      ren->Render();
    }
}


void VizViewerWindowRotationDialog::Y_value_changed(int y)
{
  if (model && ren)
    {
      labelY->setText(QString().setNum(y));
      transform->Identity();
      transform->Translate(-centroid[0],-centroid[1],-centroid[2]);
      transform->RotateY(y);
      transform->Translate(centroid);
      filter->Update();
      model->SetPolyData(filter->GetOutput());
      ren->Render();
    }
}


void VizViewerWindowRotationDialog::Z_value_changed(int z)
{
  if (model && ren)
    {
      labelZ->setText(QString().setNum(z));
      transform->Identity();
      transform->Translate(-centroid[0],-centroid[1],-centroid[2]);
      transform->RotateZ(z);
      transform->Translate(centroid);
      filter->Update();
      model->SetPolyData(filter->GetOutput());
      ren->Render();
    }
}


void VizViewerWindowRotationDialog::save_clicked()
{
  if (model && ren)
    {
      QString s = QFileDialog::getSaveFileName(this,"Save Transform",".","Transform files (*.trsf)");
      if (s.isEmpty())
	return;
      if (!(s.contains(".trsf")))
	s.append(".trsf");
      
      vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
      matrix->DeepCopy(transform->GetMatrix());
      matrix->Invert();
      itk::AffineTransform<double,3>::Pointer itkTransform = itk::AffineTransform<double,3>::New();
      itk::AffineTransform<double,3>::ParametersType params(itkTransform->GetNumberOfParameters());
      itk::TransformFileWriter::Pointer writer = itk::TransformFileWriter::New();
      
      itk::AffineTransform<double,3>::MatrixType itkMatrix;
      itk::AffineTransform<double,3>::OutputVectorType itkOffset;
      
      for (unsigned int i=0; i<3; i++)
	{
	  for (unsigned int j=0; j<3; j++)
	    itkMatrix.GetVnlMatrix().put(i,j,matrix->GetElement(i,j));

	  itkOffset[i] = matrix->GetElement(i,3);
	}

      itkTransform->SetMatrix(itkMatrix);
      itkTransform->SetOffset(itkOffset);

      writer->SetInput(itkTransform);
      writer->SetFileName(s.toStdString().c_str());
      writer->Update();
    }
}
