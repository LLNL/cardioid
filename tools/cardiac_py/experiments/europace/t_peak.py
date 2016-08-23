'''
Created on 01/07/2013

@author: butler
'''
import numpy as np
import analysis
import os


class TPeak():
    '''
    Post processing priming information for experiment 18 parameter space with
    an eye cast forward to a 'generalizable' case in the future.
    Experiment directory structure:
    priming/BCL_BCLVAL/l2_param_change/priming
    '''
    # Notes: This code base was adapted from experiments_18/priming.py
    # Some of the naming conventions are a litle strange because of it.
    #

    def __init__(self,
                 directory="/Users/butler/Desktop/CARDIOID/europace",
                 wedge_e_ver=1):
        '''
        Constructor: We will try to set all parameters in here.
        '''
        self.directory = directory
        self.ex_dir = "priming"
        self.param_dirs = ['IKr_block_000', 'IKr_block_025',
                           'IKr_block_033', 'IKr_block_050',
                           'IKr_block_066', 'IKr_block_075',
                           'IKr_block_100']
        self.case_suffix = 'priming'
        if wedge_e_ver == 1:
            self.electrodes = ['electrode#000094150', 'electrode#000092294']
        elif wedge_e_ver == 2:
            self.electrodes = ['electrode#000094790', 'electrode#000105838']
        elif wedge_e_ver == 3:
            self.electrodes = ['electrode#000093815', 'electrode#000095757']
        print analysis.__file__

    def create_list(self, BCL):
        ''' for each l1 parameter plot the set of l2 parameters as subplots
        '''
        self.BCL = BCL
        # This outer level will drive all of the
        self.case_paths = []
        self.plot_legend = []
        base_path = self.__add_sep(self.directory) + self.ex_dir  + os.sep + "BCL_{:04}".format(BCL) + "ms" + os.sep
        # Now we create the plot directories on l2
        for jj in range(len(self.param_dirs)):
            case_dir = self.__add_sep(base_path) + self.param_dirs[jj] + os.sep + self.case_suffix
            self.case_paths.append(case_dir)
        # Construct graph here

    def calc_extrema(self, beat):
        n_cases = len(self.param_dirs)
        self.t_peak_t = np.zeros((n_cases))
        self.t_peak_v = np.zeros((n_cases))
        for ii in range(n_cases):
            print 'Finding t_peak for case', self.param_dirs[ii]
            identifier = analysis.l1_analysis.Analyse(self.case_paths[ii])
            identifier.electrodes = self.electrodes
            identifier.period_ms = self.BCL
            identifier.load_data()
            identifier.construct_lead_data(flipper = -1)
            identifier.pick_beat(beat)
            t_peak, v_peak = identifier.calculate_extema_and_estimate_T_wave_max()
            print 'T wave peak time (offset in [s]): ', t_peak
            print 'T wave peak time [mV]: ', v_peak

    def __add_sep(self, input_path):
        '''
        Apply this function to a path to add a os.sep if no os.sep is trailing
        '''
        if input_path[-1] == os.sep:
            return input_path
        else:
            return (input_path + os.sep)
