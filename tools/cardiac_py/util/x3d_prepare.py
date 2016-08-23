'''
Created on Sep 23, 2012

@author: butler
'''
import glob
import os
import shutil


class X3DTools():
    '''
    classdocs
    '''
    def __init__(self, directory, substitute_file, sequence_stub):
        '''
        Constructor
        '''
        self.directory = directory
        self.substitute_file = substitute_file
        self.sequence_stub = sequence_stub
        self.threshold = 1500

    def build_valid_sequence(self):
        """ Here we will build a Valid x3d sequence from a series of x3d files
        that we have got from paraview. Here it is recognised that 'zero
        object' x3d files are problematic to blender. Therefore we substitute
        any object that is too small to contain a mesh with a valid object that
        will be outside of the expected field of view
        """
        sub_path = self.directory + os.sep + self.substitute_file
        x3d_glob = self.directory + os.sep + self.sequence_stub + "*"
        x3ds = glob.glob(x3d_glob)
        for x3d in x3ds:
            if os.path.getsize(x3d) < self.threshold:
                shutil.copy(sub_path, x3d)
        print "done"
