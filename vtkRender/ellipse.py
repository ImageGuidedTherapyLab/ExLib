import vtk
 
# FIXME - input your stl file here
stlReader = vtk.vtkSTLReader()
stlReader.SetFileName("ellipse.stl")
#stlReader.SetFileTypeToBinary()
stlReader.Update()

#coneSource = vtk.vtkConeSource()
#coneSource.SetResolution(60)
#coneSource.SetCenter(-2,0,0)
 
# Create a mapper and actor
mapper = vtk.vtkPolyDataMapper()
#mapper.SetInput(coneSource.GetOutput())
mapper.SetInput(stlReader.GetOutput())
actor = vtk.vtkActor()
actor.SetMapper(mapper)
 
# Visualize
renderer = vtk.vtkRenderer()
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)
 
renderer.AddActor(actor)
renderer.SetBackground(.1, .2, .3) # Background color dark blue
renderer.SetBackground(.3, .2, .1) # Background color dark red
renderWindow.Render()

# write render window to file
windowToImage = vtk.vtkWindowToImageFilter() 
windowToImage.SetInput(renderWindow )
windowToImage.Update()
jpgWriter     = vtk.vtkJPEGWriter() 
jpgWriter.SetFileName( 'testtwo.jpg' )
jpgWriter.SetInput(windowToImage.GetOutput())
jpgWriter.Write()


# interactive render
#renderWindowInteractor.Start()
