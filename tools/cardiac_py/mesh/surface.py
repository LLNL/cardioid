from __future__ import division
import vtk


class Surface_Smoother():
    """ Adjusting using a surface smoother on the vtk mesh provided by Slava
    """
    def __init__(self, my_file):
        self.vtk_file = my_file
        self.smooth_algorithm = 'Laplacian'
        self.subdivide = 2
        self.divide_algorith = "Butterfly"

    def read_file(self):
        self.vtk_reader = vtk.vtkUnstructuredGridReader()
        self.vtk_reader.SetFileName(self.vtk_file)
        self.vtk_reader.Update()

    def extract_surface(self):
        self.vtk_surface = vtk.vtkDataSetSurfaceFilter()
        self.vtk_surface.SetInputConnection(self.vtk_reader.GetOutputPort())
        self.vtk_surface.Update()

    def subdivide_surface(self):
        if self.divide_algorith == "Butterfly":
            self.vtk_sd = vtk.vtkButterflySubdivisionFilter()
        else:
            self.vtk_sd = vtk.vtkInterpolatingSubdivisionFilter()
        self.vtk_sd.SetNumberOfSubdivisions(self.subdivide)
        self.vtk_sd.SetInputConnection(self.vtk_surface.GetOutputPort())

    def write_vtk(self, file_name):
        """ Write vtk data.
        """
        vtk_grid_writer = vtk.vtkXMLPolyDataWriter()
        vtk_grid_writer.SetInputConnnection(self.vtk_surface.GetOutputPort())
        vtk_grid_writer.SetFileName(file_name)
        vtk_grid_writer.Write()

    def vtk_smooth(self, iterations):
        """ Run smooth on contour object and create a surface object"""
        if self.smooth_algorithm == 'Laplacian':
            self.smoother = vtk.vtkSmoothPolyDataFilter()
            self.smoother.SetInputConnection(self.vtk_sd.GetOutputPort())
            self.smoother.SetConvergence(0.0)
            self.smoother.SetNumberOfIterations(iterations)
        elif self.smooth_algorithm == 'WindowedSinc':
            print 'Smoother is: ', self.smooth_algorithm
            self.smoother = vtk.vtkWindowedSincPolyDataFilter()
            self.smoother.SetInputConnection(self.vtk_sd.GetOutputPort())
            self.smoother.SetNumberOfIterations(iterations)
        else:
            print 'Incorrect smoothing algorithm name'
            return

    def view_vtk(self):
        """ Generate winder and render the smoothed contour """
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(self.smoother.GetOutput())
        mapper.ScalarVisibilityOff()
        # Create an actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetInterpolationToFlat()
        actor.GetProperty().EdgeVisibilityOn()
        actor.GetProperty().SetColor(0.7, 0.7, 0.7)
        # Now create the rendering windows etc
        renderer = vtk.vtkRenderer()
        window = vtk.vtkRenderWindow()
        window.SetSize(400, 400)
        window.AddRenderer(renderer)
        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(window)
        renderer.SetBackground(1, 1, 1)
        renderer.AddActor(actor)
        renderer.ResetCamera()
        renderer.GetActiveCamera().Azimuth(30)
        renderer.GetActiveCamera().Elevation(20)
        renderer.GetActiveCamera().Dolly(0.9)
        renderer.ResetCameraClippingRange()
        interactor.Initialize()
        renderer.Render()
        interactor.Start()

    def write_vts(self, file_name):
        pass

    def write_stl(self, my_file):
        """ write the smoothed object
        """
        # enforce triangle structure
        triangled = vtk.vtkTriangleFilter()
        triangled.SetInputConnection(self.smoother.GetOutputPort())
        stl = vtk.vtkSTLWriter()
        stl.SetFileName(my_file)
        stl.SetInput(triangled.GetOutput())
        stl.Write()
