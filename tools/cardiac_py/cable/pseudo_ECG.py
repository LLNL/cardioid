'''
Created on 14/06/2013

@author: butler
'''
import strip_Vm
import numpy as np
import matplotlib.pyplot as plt
import os


class PseudoECG():
    '''
    classdocs:
    Calculate the pseudo ECG according to the method of Brennan et al.
    '''
    def __init__(self, directory, exclude):
        '''
        Units:
        All displacements are in mm
        All voltages are in mV
        All times are in seconds
        '''
        self.directory = directory
        self.conductivity_ratio = 0.03 / 0.12
        self.fiber_r = 10
        ''' Voxel size '''
        self.dx = 0.1
        ''' How many points are not used in the calculation. Not used includes
        Not using the distance '''
        self.exclude = exclude
        '''Offset in mm which the '''
        self.offset = 10
        '''Time between each time sample'''
        self.dt = 0.001

    def get_data(self):
        grabber = strip_Vm.Cleaner(self.directory)
        grabber.collect_directory_list()
        grabber.read_and_matrix_VM()
        self.VM = grabber.return_data()

    def calc(self):
        scale_factor = (self.fiber_r ** 2) / (4)
        (n_times, n_cells) = self.VM.shape
        self.endo_electrode = np.zeros((n_times))
        self.epi_electrode = np.zeros((n_times))
        for tt in range(n_times):
            # calc gradV
            max_i = n_cells - 1
            dv_dx = np.zeros((n_cells))
            dv_dx[0] = (self.VM[tt, 1] - self.VM[tt, 0]) / self.dx
            dv_dx[max_i] = ((self.VM[tt, max_i] - self.VM[tt, (max_i - 1)]) /
                            self.dx)
            for ii in range(1, max_i):
                dv_dx[ii] = ((self.VM[tt, (ii + 1)] - self.VM[tt, (ii - 1)]) /
                             (2 * self.dx))
            # Now calculate the endo ECG
            endo_V = 0.0
            for ii in range(self.exclude, (n_cells - self.exclude)):
                delta_d = ((1 / (self.offset + (ii + 1) * self.dx)) -
                           (1 / (self.offset + ii * self.dx)))
                endo_V = endo_V - dv_dx[ii] * delta_d * self.dx
            self.endo_electrode[tt] = scale_factor * endo_V
            epi_V = 0.0
            epi_loc = n_cells * self.dx + self.offset
            for ii in range(self.exclude, (n_cells - self.exclude)):
                delta_d = ((1 / (epi_loc - (ii + 1) * self.dx)) -
                           (1 / (epi_loc - ii * self.dx)))
                epi_V = epi_V - dv_dx[ii] * delta_d * self.dx
            self.epi_electrode[tt] = scale_factor * epi_V

    def plotECG(self):
        n_times = self.endo_electrode.size
        time = np.zeros((n_times))
        for ii in range(n_times):
            time[ii] = self.dt * ii
        transmural = self.epi_electrode[:] - self.endo_electrode[:]
        plt.plot(time, transmural)
        plt.show()
        plt.savefig((self.directory + os.sep + 'test.png'))

    def plotbeatECG(self):
        n_times = self.endo_electrode.size
        time = np.zeros((n_times))
        for ii in range(n_times):
            time[ii] = self.dt * ii
        transmural = self.epi_electrode[:] - self.endo_electrode[:]
        time_beat = time[-1000:]
        transmural_beat = transmural[-1000:]
        plt.figure()
        plt.plot(time_beat, transmural_beat)
        plt.show()
        plt.savefig((self.directory + os.sep + 'test_beat.png'))

    def save(self, save_csv):
        
        pass
    
