'''
Created on 12/06/2013

@author: butler
'''
import sys
import pio
import glob
import os
import numpy as np


class Cleaner():
    ''' This class is designed to strip vm out of snapshot files into a regular
    array that can be handed to a ECG calculator
    '''

    def __init__(self, directory):
        self.directory = directory
        self.pio_stem = 'anatomy#'
        self.sshot_match = 'snapshot.0*'

    def collect_directory_list(self):
        sshot_glob = self.directory + os.sep + self.sshot_match
        sshot_list = glob.glob(sshot_glob)
        self.sshot_list = sshot_list

    def read_and_matrix_VM(self):
        first = True
        n_times = len(self.sshot_list)
        for jj in range(n_times):
            pio_path = self.sshot_list[jj] + os.sep + self.pio_stem + "*"
            reader = pio.iter_read.Reader(pio_path)
            if first:
                n_points = int(reader.header_vars["nrecords"])
                self.data = np.zeros((n_times, n_points))
                first = False
            ii = 0
            for record in reader:
                vm = float(record.strip('\n').strip(' ').rsplit(' ')[-1])
                self.data[jj, ii] = vm
                ii = ii + 1
    def return_data(self):
        return self.data
