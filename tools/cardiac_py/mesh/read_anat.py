'''
Created on 14/01/2013

@author: butler
'''
import glob


class AnatReader():
    '''
    classdocs
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

    def __iter__(self):
        return self

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
                fields = line.rstrip("\n").strip(" ").split(";")
                # check for split b
                if fields[-1] == "":
                    # we have a standard case and we can treat it how we want
                    for field in fields:
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
                        line = self.fd.readline()
                        l2 = line.rstrip("\n").strip(" ")
                        self.h.append(l2)
                        line = self.fd.readline()
                        l3 = line.rstrip("\n").rstrip(";").strip(" ")
                        self.h.append(l3)
                        # Done with H
        line = self.fd.readline()
        print "should be nothing: ", line

    def next(self):
        line = self.fd.readline()
        if not line:
            if self.file_index < len(self.f_list) - 1:
                self.file_index = self.file_index + 1
                self.fd.close()
                self.fd = open(self.f_list[self.file_index])
                line = self.fd.readline()
            else:
                raise StopIteration
        return line
