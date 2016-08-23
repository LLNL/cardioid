import glob
import numpy as np


class Reader():
    '''
    Reader: Reads PIO data records via a iterator.
    It is designed to maintain some backwards compatibility with existing tools
    for text data. Mixed binary / ascii files makes this look a little messy

    Note that this is not the fastest way to read a file. However, it is low
    memory utilization.

    This is a nicer implementation when you want to read through the total
    content of a pio file collection as it does it seamlessly w/o bothering the
    user.
    Header data is parsed and held in a dictionary.
    '''

    def __init__(self, file_stem, debug=False):
        g_matcher = file_stem
        self.debug = debug
        self.f_list = glob.glob(g_matcher)
        self.f_list.sort()
        if len(self.f_list) <= 0:
            print "No files matching"
            raise
        self.header_vars = {}
        self.__parse_header()

    def __iter__(self):
        return self

    def __parse_header(self):
        """ Very simple parser"""
        self.header_open = False
        self.file_index = 0
        if self.debug:
            print self.f_list
        self.fd = open(self.f_list[self.file_index], 'rb')
        while True:
            line = self.read_a_line()
            if not self.header_open:
                if "{" in line:
                    self.header_open = True
            else:
                if "//" in line:
                    continue
                if "}" in line:
                    break
                fields = line.rstrip("\n").strip(" ").split(";")
                # check for split b
                if fields[-1] == "":
                    # we have a standard case and we can treat it how we want
                    for field in fields:
                        if self.debug:
                            print field
                        if len(field) == 0:
                            continue
                        keypair = field.split("=")
                        if len(keypair) == 2:
                            self.header_vars[keypair[0].strip(" ")] = \
                            keypair[1].strip(" ")
                elif len(fields) == 1:
                    # we have something else check if h. .. will be first field
                    field1 = fields[0]
                    prospective_h = field1.split("=")[0].strip(" ")
                    prospective_rhs = field1.split("=")[1].strip(" ")
                    if prospective_h == "h":
                        self.h = [prospective_rhs]
                        line = self.read_a_line()
                        l2 = line.rstrip("\n").strip(" ")
                        self.h.append(l2)
                        line = self.read_a_line()
                        l3 = line.rstrip("\n").rstrip(";").strip(" ")
                        self.h.append(l3)
                        # Done with H
        line = self.read_a_line()
        if self.debug:
            print "should be nothing: ", line
        # Done the actual reading now setup parameters
        self.lrec = int(self.header_vars['lrec'])
        if self.header_vars['datatype'] == 'FIXRECORDASCII':
            self.type = 1
        elif self.header_vars['datatype'] == 'FIXRECORDBINARY':
            self.type = 2
        elif self.header_vars['datatype'] == 'VARRECORDASCII':
            self.type = 3
        else:
            print 'Warning!!!! unpredictable / unknown result as datatype'
            print 'is unknown.'

    def next(self):
        line = self.read_record()
        if len(line) == 0:
            if self.file_index < len(self.f_list) - 1:
                self.file_index = self.file_index + 1
                self.fd.close()
                self.fd = open(self.f_list[self.file_index])
                line = self.read_record()
            else:
                raise StopIteration
        return line

    def read_a_line(self):
        # hello
        buf = ""
        while True:
            buf += self.fd.read(1)
            if buf[-1] == "\n":
                return buf

    def read_record(self):
        if self.type == 1:
            return self.fd.read(self.lrec)
        elif self.type == 2:
            return np.fromfile(self.fd, dtype='f4,f4,f4', count=1)
        elif self.type == 3:
            return self.fd.read_a_line()
        else:
            pass
