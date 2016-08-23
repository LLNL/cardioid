
import double_res
import read_anat
import scale_conductivity
__all__ = ['double_res', 'read_anat', 'scale_conductivity']
vtk_available = False
try:
    import vtk
    vtk_available = True
except ImportError:
    pass

if vtk_available:
    import crop
    import surface
    __all__.append("crop")
    __all__.append("surface")
