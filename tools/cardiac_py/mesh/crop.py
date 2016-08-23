'''
Created on 16/01/2013

@author: butler
'''
import vtk
from . import read_anat
import os


class Crop(object):
    '''
    Simple class to crop out snapshot or anat files.
    Note that it does not
    '''

    def __init__(self, old_dir, new_dir):
        """
        """
        self.old_dir = old_dir
        self.debug = False
        self.new_dir = new_dir
        self.scale = int(2)
        self.anat_match = 'anatomy#*'
        self.mapping_file = 'mapping.txt'
        self.sensor_file = 'sensor.txt'
        self.n_files = 1
        self.new_dir
        old_g_matcher = self.old_dir + os.sep + self.anat_match
        self.old_file = read_anat.AnatReader(old_g_matcher)
        self.new_sshot_dir = self.new_dir
        self.o_num_records = int(self.old_file.header_vars['nrecords'])
        self.o_nx = 480  # int(self.old_file.header_vars["nx"])
        self.o_ny = 520  # int(self.old_file.header_vars["ny"])
        self.o_nz = 460  # int(self.old_file.header_vars["nz"])
        self.x_min = 150
        self.x_max = 222
        self.y_min = 130
        self.y_max = 208
        self.z_min = 300
        self.z_max = 392

    def __write_header(self, fd, ftype):
        fd.write("anatomy FILEHEADER {\n")
        fd.write("  datatype = FIXRECORDASCII;\n")
        fd.write(("  nfiles = " + repr(self.n_files) + ";\n"))
        fd.write(("  nrecord = " + repr(self.n_new_records) + ";\n"))
        fd.write(("  nfields = " + self.old_file.header_vars["nfields"] + ";\n"
                  ))
        fd.write(("  lrec = 84;\n"))
        fd.write(("  endian_key = " + self.old_file.header_vars["endian_key"] +
                  ";\n"))
        fd.write(("  field_names = " + self.old_file.header_vars["field_names"]
                  + ";\n"))
        fd.write(("  field_types = " + self.old_file.header_vars["field_types"]
                  + ";\n"))
        fd.write(("  field_units = " + self.old_file.header_vars["field_units"]
                  + ";\n"))
        fd.write(("  nx = " + repr(self.n_nx) + ";\n"))
        fd.write(("  ny = " + repr(self.n_ny) + ";\n"))
        fd.write(("  nz = " + repr(self.n_nz) + ";\n"))
        fd.write(("  dx = " + repr(self.n_dx) + ";\n"))
        fd.write(("  dy = " + repr(self.n_dy) + ";\n"))
        fd.write(("  dz = " + repr(self.n_dz) + ";\n"))
        fd.write(("}\n\n"))

    def _WriteSampledSnapshot(self):
        print "Writing sampled snapshot. n_records is to be changed"
        pass

    def _WriteSampledAnatomy(self):
        pass

    def WriteSampledFile(self):
        """ Uses whether gid exists to determine whether a snapshot.initial
        anatomy or if the file is a snapshot.00* cellVis
        """
        pass
        '''
        fields = self.old_file.header_vars["field_names"]
        fields = fields.strip(" ").split()
        if "gid" in fields:
            self._WriteSampledAnatomy()
        else:
            self.
        '''

    def CreateSampledvtkFile(self):
        # 1: presume we are reading a cellVis snapshot
        # 2: use x-axis limits to create a smaller matrix
        # 3: be lazy and use a full matrix - preallocation may make it faster
        # 4: This shirt is too tight.
        # 5: Use vtk python bindings.
        n_x = self.x_max - self.x_min + 1
        n_y = self.y_max - self.y_min + 1
        n_z = self.z_max - self.z_min + 1
        # this offset is designed to put x_min, y_min, z_min at 0,0,0
        offset_x = self.o_nx / 2 - self.x_min
        offset_y = self.o_ny / 2 - self.y_min
        offset_z = self.o_nz / 2 - self.z_min
        # determine var locations
        field_names = self.old_file.header_vars["field_names"].split()
        i_rx = field_names.index("rx")
        i_ry = field_names.index("ry")
        i_rz = field_names.index("rz")
        i_vm = field_names.index("Vm")
        grid = vtk.vtkImageData()
        grid.SetOrigin(0, 0, 0)
        dx = 0.2  # float(self.old_file.header_vars["dx"])
        dy = 0.2  # float(self.old_file.header_vars["dy"])
        dz = 0.2  # float(self.old_file.header_vars["dz"])
        grid.SetSpacing(dx, dy, dz)
        grid.SetDimensions((n_x + 1), (n_y + 1), (n_z + 1))
        array = vtk.vtkDoubleArray()
        array.SetNumberOfComponents(1)
        array.SetNumberOfTuples(grid.GetNumberOfCells())
        print "N tuples: ", array.GetNumberOfTuples()
        for ii in range(array.GetNumberOfTuples()):
            array.SetValue(ii, -1000.0)
        array.SetName("Vm")
        print 'parsing records'
        for record in self.old_file:
            fields = record.split()
            x = int(fields[i_rx]) + offset_x
            y = int(fields[i_ry]) + offset_y
            z = int(fields[i_rz]) + offset_z
            vm = float(fields[i_vm])
            if ((x < n_x) and (x >= 0) and (y < n_y) and
                (y >= 0) and (z < n_z) and (z >= 0)):
                # Now lets calculate the effective GID
                vtkGID = x + n_x * (y + n_y * z)
                if vtkGID < grid.GetNumberOfCells():
                    array.SetValue(vtkGID, vm)
                else:
                    print 'to big'
        print 'Array constructed'
        grid.GetCellData().AddArray(array)
        print 'Added to grid'
        writer = vtk.vtkXMLImageDataWriter()
        writer.SetFileName((self.new_dir + os.sep + "subset.vti"))
        writer.SetInputConnection(grid.GetProducerPort())
        print 'writing'
        writer.Write()
        print 'done'
