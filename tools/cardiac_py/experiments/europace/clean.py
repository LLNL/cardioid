'''
Created on 11/07/2013

@author: butler
'''
import glob
import util
import os


class Cleaner(object):
    '''
    Cleaner: Utility class to remove unnecessary checkpoints in the europace
    parameter space
    '''

    def __init__(self, directory, mode=0):
        '''
        Mode defaults to safe mode.
        '''
        self.directory = directory
        self.mode = mode

    def set_s2_location(self, loc_num):
        self.s2_path = self.__add_sep(self.directory) + \
        "s_2_loc_{0:1}".format(loc_num)

    def clean_BCL(self, BCL):
        BCL_path = self.__add_sep(self.s2_path) + "BCL_{0:04}ms".format(BCL)
        IKr_glob = self.__add_sep(BCL_path) + "IKr_block*"
        IKr_dirs = glob.glob(IKr_glob)
        for IKr_dir in IKr_dirs:
            priming_dir = self.__add_sep(IKr_dir) + "priming"
            util.file_manipulation.remove_checkpoints(priming_dir, self.mode)
            s2_glob = self.__add_sep(IKr_dir) + "s2_offset_*"
            s2_dirs = glob.glob(s2_glob)
            for s2_dir in s2_dirs:
                util.file_manipulation.remove_checkpoints(s2_dir, self.mode)

    def __add_sep(self, input_path):
        '''
        Apply this function to a path to add a os.sep if no os.sep is trailing
        '''
        if input_path[-1] == os.sep:
            return input_path
        else:
            return (input_path + os.sep)
