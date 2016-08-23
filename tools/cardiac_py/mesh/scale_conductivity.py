'''
Created on 14/01/2013

@author: butler
'''
import pio
import os
import math


class Scaler():
    '''
    classdocs
    '''

    def __init__(self, new_dir, old_dir, scale_factor):
        '''
        Constructor
        '''
        self.scale_factor = scale_factor
        self.old_dir = old_dir
        self.debug = False
        self.new_dir = new_dir
        self.sshot_dir_name = 'snapshot.initial'
        self.scale = int(1)
        self.anat_match = 'anatomy#*'
        self.mapping_file = 'mapping.txt'
        self.sensor_file = 'sensor.txt'
        self.n_files = 50
        old_g_matcher = self.old_dir + os.sep + self.sshot_dir_name + os.sep \
        + self.anat_match
        self.old_file = pio.iter_read.Reader(old_g_matcher)
        self.new_sshot_dir = self.new_dir + os.sep + self.sshot_dir_name
        self.o_num_records = long(self.old_file.header_vars['nrecord'])
        self.n_new_records = (self.scale ** 3) * self.o_num_records
        self.o_nx = int(self.old_file.header_vars["nx"])
        self.o_ny = int(self.old_file.header_vars["ny"])
        self.o_nz = int(self.old_file.header_vars["nz"])
        self.n_nx = self.scale * self.o_nx
        self.n_ny = self.scale * self.o_ny
        self.n_nz = self.scale * self.o_nz
        self.n_dx = float(self.old_file.header_vars["dx"])
        self.n_dy = float(self.old_file.header_vars["dy"])
        self.n_dz = float(self.old_file.header_vars["dz"])

    def __write_anat_header(self, fd):
        fd.write("anatomy FILEHEADER {\n")
        fd.write("  datatype = FIXRECORDASCII;\n")
        fd.write(("  nfiles = " + repr(self.n_files) + ";\n"))
        fd.write(("  nrecord = " + '{0:n}'.format(self.n_new_records) + ";\n"))
        fd.write(("  nfields = " + self.old_file.header_vars["nfields"]
                  + ";\n"))
        fd.write(("  lrec = 84;\n"))
        fd.write(("  endian_key = " + self.old_file.header_vars["endian_key"]
                  + ";\n"))
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

    def WriteAnat(self):
        try:
            if not os.path.exists(self.new_dir):
                os.mkdir(self.new_dir)
            if not os.path.exists(self.new_sshot_dir):
                os.mkdir(self.new_sshot_dir)
            anat_f_index = 0
            rec_count = 0
            n_rec_split = math.ceil(self.o_num_records / self.n_files)
            anat_file = open(self.new_sshot_dir + os.sep +
                             "anatomy#{0:06d}".format(anat_f_index), 'w')
            self.__write_anat_header(anat_file)
            for record in self.old_file:
                if len(record) < 60:
                    print "warning failed line"
                    continue
                split_strings = record.strip("\n").strip().split()
                record_out = split_strings[0] + " " + split_strings[1] + " "
                t1 = self.scale_factor * float(split_strings[2])
                t2 = self.scale_factor * float(split_strings[3])
                t3 = self.scale_factor * float(split_strings[4])
                t4 = self.scale_factor * float(split_strings[5])
                t5 = self.scale_factor * float(split_strings[6])
                t6 = self.scale_factor * float(split_strings[7])
                record_out = record_out + '{0:10f}'.format(t1) + " " + \
                '{0:10f}'.format(t2) + " " + '{0:10f}'.format(t3) + " " + \
                '{0:10f}'.format(t4) + " " + '{0:10f}'.format(t5) + " " + \
                '{0:10f}'.format(t6) + "\n"
                record_out = record_out.rjust(84)
                anat_file.write(record_out)
                # Use old record count to split files
                rec_count = rec_count + 1
                if (rec_count % n_rec_split == 0) and ((anat_f_index + 1) <
                                                       self.n_files):
                    anat_file.close()
                    anat_f_index = anat_f_index + 1
                    anat_file = open(self.new_sshot_dir + os.sep +
                                "anatomy#{0:06d}".format(anat_f_index), 'w')
            anat_file.close()
        except ValueError:
            print self.old_file.f_list[self.old_file.file_index]
            print 'record: ', record
            print 'rec_count: ', rec_count
            print 'rec_out was: ', record_out
            raise
