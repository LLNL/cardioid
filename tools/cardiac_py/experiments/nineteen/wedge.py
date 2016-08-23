'''
Created on 29/04/2013

@author: butler
'''
import gecg
import os
import numpy as np


class Priming():
    '''
    Post processing priming information for experiment 18 parameter space with
    an eye cast forward to a 'generalizable' case in the future.
    Experiment directory structure:
    root_directory/BCL_BCLVAL/l1_param_change/l2_param_change/priming
    ""             ""             ""              ""         /s2_run_*
    '''

    def __init__(self,
                 directory="/Users/butler/Desktop/CARDIOID/ECG_data"):
        '''
        Constructor: We will try to set all parameters in here.
        '''
        self.directory = directory
        self.ex_dir = "wedge_ex_19_a"
        self.l2_param_dirs = ['ikr_block_0/gNaL015', 'ikr_block_66/c_3',
                              'ikr_block_66/ranolazine']
        self.l2_param_labels = ['WT', 'E4031', 'Ranolazine']
        self.plotter = gecg.GECGs(self.directory, full=False)
        self.plotter.template.figfractionalwidth = 0.85

    def set_analysis_type(self, BCL, sweep):
        ''' set_analysis_type: Driver routine for the selection of graphing
        analysis. Specific types are listed below and selected via option.
        Given these are priming runs we effectively have two choices:
        1) alter the l1 between files and allow the l2 parameter vary
        for a graph
        2) alter the l2 between files and allow the l1 parameter vary
        for a grapjh
        '''
        pass

    def subplot_l2_comparison(self, BCL):
        ''' for each l1 parameter plot the set of l2 parameters as subplots
        '''
        self.plotter.period_ms = BCL

        # This outer level will drive all of the
        for ii in range(len(self.l1_param_dirs)):
            self.plot_directories = []
            self.plot_legend = []
            l1_path = self.ex_dir + os.sep + "BCL_" \
            + repr(BCL) + os.sep + self.l1_param_dirs[ii] + os.sep
            plot_name = ('BCL_' + repr(BCL) + "_" + self.l1_param_labels[ii] +
                         '.pdf')
            # Now we create the plot directories on l2
            for jj in range(len(self.l2_param_dirs)):
                plot_dir = l1_path + os.sep + self.l2_param_dirs[jj]
                self.plot_directories.append(plot_dir)
                self.plot_legend.append(self.l2_param_labels[jj])
            # Construct graph here
            self.plotter.case_directories = self.plot_directories
            self.plotter.case_legends = self.plot_legend
            self.plotter.load_data()
            self.plotter.set_ECG_type(0, -1)
            self.plotter.full_set_as_subplots(plot_name)

    def overlay_l2_comparison(self, BCL):
        ''' for each l1 parameter plot the set of l2 parameters as subplots
        '''
        self.plotter.period_ms = BCL

        # This outer level will drive all of the
        for ii in range(len(self.l1_param_dirs)):
            self.plot_directories = []
            self.plot_legend = []
            l1_path = self.ex_dir + os.sep + "BCL_" \
            + repr(BCL) + os.sep + self.l1_param_dirs[ii] + os.sep
            plot_name = ('BCL_' + repr(BCL) + "_" + self.l1_param_labels[ii] +
                         '.pdf')
            # Now we create the plot directories on l2
            for jj in range(len(self.l2_param_dirs)):
                plot_dir = l1_path + os.sep + self.l2_param_dirs[jj]
                self.plot_directories.append(plot_dir)
                self.plot_legend.append(self.l2_param_labels[jj])
            # Construct graph here
            self.plotter.case_directories = self.plot_directories
            self.plotter.case_legends = self.plot_legend
            self.plotter.load_data()
            self.plotter.set_ECG_type(0, -1)
            self.plotter.full_ECG_set(plot_name)

    def beat_comparison(self, BCL):
        ''' for each l1 parameter plot the set of l2 parameters as subplots
        '''
        self.plotter.period_ms = BCL
        self.plot_directories = []
        self.plot_legend = []
        l1_path = self.ex_dir + os.sep + "BCL_" \
        + repr(BCL) + os.sep
        plot_name = 'BCL_' + repr(BCL)
        # Now we create the plot directories on l2
        for jj in range(len(self.l2_param_dirs)):
            plot_dir = l1_path + os.sep + self.l2_param_dirs[jj]
            print plot_dir
            self.plot_directories.append(plot_dir)
            self.plot_legend.append(self.l2_param_labels[jj])
        # Construct graph here
        self.plotter.case_directories = self.plot_directories
        self.plotter.case_legends = self.plot_legend
        self.plotter.load_data()
        self.plotter.set_ECG_type(0, 1)
        np.savetxt('wedge_data.csv',self.plotter.ECG_data,delimiter=',')
        #cd self.plotter.beat_ECG_Set(9, plot_name + "_beat_9.pdf")
        #self.plotter.beat_ECG_Set(8, plot_name + "_beat_8.pdf")
        #self.plotter.beat_ECG_Set(7, plot_name + "_beat_7.pdf")
        #self.plotter.beat_ECG_Set(6, plot_name + "_beat_6.pdf")
        self.plotter.beat_ECG_Set(5, plot_name + "_beat_5.pdf")

    def beat_l1_comparison(self, BCL):
        ''' for each l2 parameter plot the set of l1 parameters as subplots
        '''
        self.plotter.period_ms = BCL
        for ii in range(len(self.l2_param_dirs)):
            self.plot_directories = []
            self.plot_legend = []
            l2_dir = self.l2_param_dirs[ii]
            plot_name = 'BCL_' + repr(BCL) + "_" + self.l2_param_labels[ii]
            # Now we create the plot directories on l2
            for jj in range(len(self.l1_param_dirs)):
                plot_dir = self.ex_dir + os.sep + "BCL_" + repr(BCL) + os.sep \
                           + self.l1_param_dirs[jj] + os.sep + l2_dir
                self.plot_directories.append(plot_dir)
                self.plot_legend.append(self.l1_param_labels[jj])
            # Construct graph here
            self.plotter.case_directories = self.plot_directories
            self.plotter.case_legends = self.plot_legend
            self.plotter.load_data()
            self.plotter.set_ECG_type(0, -1)
            self.plotter.beat_ECG_Set(9, plot_name + "_beat_9.pdf")
            self.plotter.beat_ECG_Set(8, plot_name + "_beat_8.pdf")
            self.plotter.beat_ECG_Set(7, plot_name + "_beat_7.pdf")
            self.plotter.beat_ECG_Set(6, plot_name + "_beat_6.pdf")
            self.plotter.beat_ECG_Set(5, plot_name + "_beat_5.pdf")
