'''
Created on 29/04/2013

@author: butler
'''
import gecg
import os


class Priming():
    '''
    Post processing priming information for experiment 18 parameter space with
    an eye cast forward to a 'generalizable' case in the future.
    Experiment directory structure:
    priming/BCL_BCLVAL/l2_param_change/priming
    '''
    # Notes: This code base was adapted from experiments_18/priming.py
    # Some of the naming conventions are a litle strange because of it.

    def __init__(self,
                 directory="/Users/butler/Desktop/CARDIOID/europace"):
        '''
        Constructor: We will try to set all parameters in here.
        '''
        self.directory = directory
        self.ex_dir = "priming"
        self.l2_param_dirs = ['IKr_block_000', 'IKr_block_025',
                              'IKr_block_033', 'IKr_block_050',
                              'IKr_block_066', 'IKr_block_075',
                              'IKr_block_100']
        self.l2_param_labels = [r'$I_{Kr}Block=0\%$', r'$I_{Kr}Block=25\%$',
                                r'$I_{Kr}Block=33\%$', r'$I_{Kr}Block=50\%$',
                                r'$I_{Kr}Block=66\%$', r'$I_{Kr}Block=75\%$',
                                r'$I_{Kr}Block=100\%$']
        self.l2_suffix = 'priming'
        self.plotter = gecg.GECGs(self.directory,full=False,)
        self.plotter.template.figfractionalwidth = 0.85

    def complete_as_subplots(self, BCL):
        ''' for each l1 parameter plot the set of l2 parameters as subplots
        '''
        self.plotter.period_ms = BCL

        # This outer level will drive all of the
        self.plot_directories = []
        self.plot_legend = []
        base_path = self.__add_sep(self.ex_dir) + "BCL_{:04}".format(BCL) + "ms" + os.sep
        plot_name = ('BCL_' + repr(BCL) + "_" + self.l2_suffix + "_" +
                     "subplots" + '.pdf')
        # Now we create the plot directories on l2
        for jj in range(len(self.l2_param_dirs)):
            plot_dir = self.__add_sep(base_path) + self.l2_param_dirs[jj] + os.sep + self.l2_suffix
            self.plot_directories.append(plot_dir)
            self.plot_legend.append(self.l2_param_labels[jj])
        # Construct graph here
        self.plotter.case_directories = self.plot_directories
        self.plotter.case_legends = self.plot_legend
        self.plotter.load_data()
        self.plotter.set_ECG_type(0, 1)
        self.plotter.full_set_as_subplots(plot_name)

    def complete_as_overlay(self, BCL):
        self.plotter.period_ms = BCL
        # This outer level will drive all of the
        self.plot_directories = []
        self.plot_legend = []
        base_path = self.__add_sep(self.ex_dir) + "BCL_{:04}".format(BCL) + "ms" + os.sep
        plot_name = ('BCL_' + repr(BCL) + "_" + self.l2_suffix + "_" +
                     "overlay" + '.pdf')
        # Now we create the plot directories on l2
        for jj in range(len(self.l2_param_dirs)):
            plot_dir = self.__add_sep(base_path) + self.l2_param_dirs[jj] + os.sep + self.l2_suffix
            self.plot_directories.append(plot_dir)
            self.plot_legend.append(self.l2_param_labels[jj])
        # Construct graph here
        self.plotter.case_directories = self.plot_directories
        self.plotter.case_legends = self.plot_legend
        self.plotter.load_data()
        self.plotter.set_ECG_type(0, 1)
        self.plotter.full_ECG_set(plot_name)

    def complete_as_single(self, BCL):
        self.plotter.period_ms = BCL
        # This outer level will drive all of the
        self.plot_directories = []
        self.plot_legend = []
        base_path = self.__add_sep(self.ex_dir) + "BCL_{:04}".format(BCL) + "ms" + os.sep
        # Now we create the plot directories on l2
        for jj in range(len(self.l2_param_dirs)):
            plot_dir = self.__add_sep(base_path) + self.l2_param_dirs[jj] + os.sep + self.l2_suffix
            self.plot_directories.append(plot_dir)
            self.plot_legend.append(self.l2_param_labels[jj])
        self.plotter.case_directories = self.plot_directories
        self.plotter.case_legends = self.plot_legend
        self.plotter.load_data()
        self.plotter.set_ECG_type(0, 1)
        for ii in range(len(self.l2_param_dirs)):
            plot_name = ('BCL_' + repr(BCL) + "_" + self.l2_suffix + "_" +
                     self.l2_param_dirs[ii] + '.pdf')
            self.plotter.full_ECG_case(ii, plot_name)


    def beat_as_subplots(self, BCL, beat):
        ''' for each l1 parameter plot the set of l2 parameters as subplots
        '''
        self.plotter.period_ms = BCL
        self.plot_directories = []
        self.plot_legend = []
        base_path = self.__add_sep(self.ex_dir) + "BCL_{:04}".format(BCL) + "ms" + os.sep
        plot_name = ('BCL_' + repr(BCL) + "_" + self.l2_suffix + "_" +
                     "_beat_" + repr(beat) + ".pdf")
        # Now we create the plot directories on l2
        for jj in range(len(self.l2_param_dirs)):
            plot_dir = self.__add_sep(base_path) + self.l2_param_dirs[jj] + os.sep + self.l2_suffix
            self.plot_directories.append(plot_dir)
            self.plot_legend.append(self.l2_param_labels[jj])
        # Construct graph here
        self.plotter.case_directories = self.plot_directories
        self.plotter.case_legends = self.plot_legend
        self.plotter.load_data()
        self.plotter.set_ECG_type(0, 1)
        self.plotter.beat_set_as_subplots(beat, plot_name)

    def beat_as_overlay(self, BCL, beat):
        ''' for each l1 parameter plot the set of l2 parameters as subplots
        '''
        self.plotter.period_ms = BCL
        self.plot_directories = []
        self.plot_legend = []
        base_path = self.__add_sep(self.ex_dir) + "BCL_{:04}".format(BCL) + "ms" + os.sep
        plot_name = ('BCL_' + repr(BCL) + "_" + self.l2_suffix + "_" +
                     "_beat_" + repr(beat) + ".pdf")
        # Now we create the plot directories on l2
        for jj in range(len(self.l2_param_dirs)):
            plot_dir = self.__add_sep(base_path) + self.l2_param_dirs[jj] + os.sep + self.l2_suffix
            self.plot_directories.append(plot_dir)
            self.plot_legend.append(self.l2_param_labels[jj])
        # Construct graph here
        self.plotter.case_directories = self.plot_directories
        self.plotter.case_legends = self.plot_legend
        self.plotter.load_data()
        self.plotter.set_ECG_type(0, 1)
        self.plotter.beat_ECG_Set(beat, plot_name)

    def beat_as_single(self, BCL, beat):
        ''' for each l1 parameter plot the set of l2 parameters as subplots
        '''
        self.plotter.period_ms = BCL
        self.plot_directories = []
        self.plot_legend = []
        base_path = self.__add_sep(self.ex_dir) + "BCL_{:04}".format(BCL) + "ms" + os.sep
        plot_name = ('BCL_' + repr(BCL) + "_" + self.l2_suffix + "_" +
                     self.l2_param_dirs[ii] + "_beat_" + repr(beat))
        # Now we create the plot directories on l2
        for jj in range(len(self.l2_param_dirs)):
            plot_dir = self.__add_sep(base_path) + self.l2_param_dirs[jj] + os.sep + self.l2_suffix
            self.plot_directories.append(plot_dir)
            self.plot_legend.append(self.l2_param_labels[jj])
        # Construct graph here
        self.plotter.case_directories = self.plot_directories
        self.plotter.case_legends = self.plot_legend
        self.plotter.load_data()
        self.plotter.set_ECG_type(0, 1)
        for jj in range(len(self.l2_param_dirs)):
            plot_name_r = plot_name + "_" + self.l2_param_dirs[jj] + ".pdf"
            self.plotter.beat_ECG_case(beat, jj, plot_name_r)

    def __add_sep(self, input_path):
        '''
        Apply this function to a path to add a os.sep if no os.sep is trailing
        '''
        if input_path[-1] == os.sep:
            return input_path
        else:
            return (input_path + os.sep)
