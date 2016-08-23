'''
Created on 18/06/2013

@author: butler
'''
'''
Created on 14/06/2013

@author: butler
'''
import strip_Vm
import numpy as np
import matplotlib.pyplot as plt
import os


class APD():
    '''
    classdocs:
    Calculate the APD an associated properties
    '''
    def __init__(self, directory, exclude = 0):
        '''
        Units:
        All displacements are in mm
        All voltages are in mV
        All times are in seconds
        '''
        self.directory = directory
        self.conductivity_ratio = 0.03 / 0.12
        self.fiber_r = 10
        ''' Voxel size mm '''
        self.dx = 0.1
        ''' How many points are not used in the calculation. Not used includes
        Not using the distance '''
        self.exclude = exclude

        '''Time between each time sample'''
        self.dt = 0.001
        self.period_s = 1.0
        ''' Membrane voltage threshold in mV '''
        self.threshold = -75

    def get_data(self):
        grabber = strip_Vm.Cleaner(self.directory)
        grabber.collect_directory_list()
        grabber.read_and_matrix_VM()
        self.VM = grabber.return_data()
        
    def calc_apd(self):
        (n_times, n_cells) = self.VM.shape
        n_t_req = int(round(self.period_s / self.dt))
        self.last_beat = self.VM[-n_t_req:, :]
        self.activation = np.zeros((n_cells))
        self.deactivation = np.zeros((n_cells))
        self.apd = np.zeros((n_cells))
        # set neg vals as status not changed
        self.activation[:] = -1.0
        self.deactivation[:] = -1.0
        self.apd[:] = -1.0
        for ii in range(n_t_req):
            time = ii * self.dt
            for jj in range(n_cells):
                if (self.activation[jj] < 0 and
                    self.last_beat[ii, jj] > self.threshold):
                    self.activation[jj] = time
                elif (self.activation[jj] > 0 and
                      self.last_beat[ii, jj] < self.threshold and
                      self.deactivation[jj] < 0):
                    self.deactivation[jj] = time

    def calc_conduct_vel(self):
        ''' this presumes that there is a roughly linear relationship between
        location and activation time
        '''
        (n_times, n_cells) = self.VM.shape
        self.location = np.zeros((n_cells))
        for ii in range(self.location.size):
            self.location[ii] = ii * self.dx
        A = np.array([self.activation[:], np.ones(self.location.size)])
        w = np.linalg.lstsq(A.T, self.location)[0]
        print w
