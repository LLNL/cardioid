'''
Created on 04/03/2013

@author: butler
'''
import glob
import os

class Reader():
    '''
    seeker_read.Reader: Allows reading of a PIO fileset as an indexable object.
    This is nice for random access across files.
    '''

    def __init__(self, file_stem):
        g_matcher = file_stem
        self.f_list = glob.glob(g_matcher)
        self.f_list.sort()
        if len(self.f_list) < 0:
            print "No files matching"
            raise
        self.header_vars = {}
        self.__parse_header()
        self.__getfile_sizes()

    def __iter__(self):
        return self

    def __getfile_sizes(self):
        self.file_sizes = []
        for file_name in self.f_list:
            self.file_sizes.append(os.stat(file_name).st_size)

    def __parse_header(self):
        """ Very simple parser"""
        self.header_open = False
        self.file_index = 0
        print self.f_list
        self.fd = open(self.f_list[self.file_index])
        while True:
            line = self.fd.readline()
            if not self.header_open:
                if "{" in line:
                    self.header_open = True
            else:
                if "//" in line:
                    continue
                if "}" in line:
                    break
                fields = line.rstrip("\n").split(";")
                for field in fields:
                    print field
                    if len(field) == 0:
                        continue
                    keypair = field.split("=")
                    if len(keypair) == 2:
                        self.header_vars[keypair[0].strip(" ")] = \
                        keypair[1].strip(" ")
        line = self.fd.readline()
        print "should be nothing: ", line
        self.header_end = self.fd.tell()
        self.field_width = int(self.header_vars['lrec'].rstrip(";").strip(' '))
        self.records = int(self.header_vars['nrecords'].rstrip(";").strip(' '))
        self.fd.close()

    def __getitem__(self, index):
        assert index < self.records
        g_offset = self.header_end + self.field_width * index
        accumulated = 0
        l_offset = 0
        ii = 0
        for l_size in self.file_sizes:
            l_offset = g_offset - accumulated
            if l_offset < l_size:
                break
            else:
                accumulated = accumulated + l_size
                ii = ii + 1
        fd = open(self.f_list[ii], 'rb')
        fd.seek(l_offset)
        return fd.read(self.field_width)
