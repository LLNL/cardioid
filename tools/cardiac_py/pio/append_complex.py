'''
Created on 14/05/2013

@author: butler
'''

import iter_read
import os
import math


class Piofy():
    '''
    pio.append_complex.Piofy
    Take a vaild pio file, take a complex numpy vector where the data has been
    derived from the same set of pio files, and write a new pio file where the
    complex vector is appended as 'Re' and 'Iz' fields to the dataset.
    '''

    def __init__(self, old_pio_match, new_stem, data):
        '''
        Constructor
        '''
        self.data = data[0, :]
        self.debug = False
        self.old_pio = old_pio_match
        self.new_stem = new_stem
        self.sshot_dir = self.old_pio.rpartition(os.sep)[0]
        self.n_files = 5
        self.old_file = iter_read.Reader(self.old_pio)
        self.new_vars = 2
        self.new_rec_length = int(self.old_file.header_vars["lrec"]) + \
                              self.new_vars * 24

    def __write_anat_header(self, fd):
        fd.write("cellViz FILEHEADER {\n")
        fd.write("  datatype = FIXRECORDASCII;\n")
        fd.write(("  nfiles = " + repr(self.n_files) + ";\n"))
        fd.write(("  nrecords = " + self.old_file.header_vars["nrecords"] +
                  ";\n"))
        fd.write(("  nfields = " +
                  repr(int(self.old_file.header_vars["nfields"]) +
                       self.new_vars) + ";\n"))
        fd.write(("  lrec = " + repr(self.new_rec_length) + ";\n"))
        fd.write(("  field_names = " + self.old_file.header_vars["field_names"]
                  + " Re Iz" + ";\n"))
        fd.write(("  field_types = " + self.old_file.header_vars["field_types"]
                  + " f f" + ";\n"))
        fd.write(("  h = " + self.old_file.h[0] + "\n"))
        fd.write(("      " + self.old_file.h[1] + "\n"))
        fd.write(("      " + self.old_file.h[2] + ";\n"))
        fd.write(("}\n\n"))

    def Write(self):
        try:
            anat_f_index = 0
            rec_count = 0
            n_rec_split = math.ceil(float(self.old_file.header_vars["nrecords"]
                                          ) / self.n_files)
            anat_file = open(self.sshot_dir + os.sep + self.new_stem
                             + "#{0:06d}".format(anat_f_index), 'w')
            self.__write_anat_header(anat_file)
            for record in self.old_file:
                if len(record) < 20:
                    print "warning failed line"
                    continue
                else:
                    new_record = record
                    new_record = (record.rstrip("\n") + " " +
                                  repr(self.data[rec_count].real) + " " +
                                  repr(self.data[rec_count].imag) + '\n')
                    new_record = new_record.rjust(self.new_rec_length)
                    anat_file.write(new_record)
                # Use old record count to split files
                rec_count = rec_count + 1
                if (rec_count % n_rec_split == 0) and ((anat_f_index + 1) <
                                                       self.n_files):
                    anat_file.close()
                    anat_f_index = anat_f_index + 1
                    anat_file = open(self.sshot_dir + os.sep + self.new_stem
                             + "#{0:06d}".format(anat_f_index), 'w')
            anat_file.close()
        except ValueError:
            print self.old_file.f_list[self.old_file.file_index]
            print 'record: ', record
            print 'rec_count: ', rec_count
