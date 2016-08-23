'''
Created on 29/04/2013

@author: butler
'''
import gecg
import os


class s2_plotter():
    '''
    Post processing priming information for experiment 18 parameter space with
    an eye cast forward to a 'generalizable' case in the future.
    Experiment directory structure:
    root_directory/BCL_BCLVAL/l1_param_change/l2_param_change/priming
    ""             ""             ""              ""         /s2_run_*
    '''

    def __init__(self,
                 directory="/Users/butler/Desktop/CARDIOID/europace"):
        '''
        Constructor: We will try to set all parameters in here.
        '''
        self.directory = directory
        self.ex_dir = "s_2_loc_2"
        self.l1_param_dirs = ['IKr_block_000', 'IKr_block_025',
                              'IKr_block_033', 'IKr_block_050',
                              'IKr_block_066', 'IKr_block_075',
                              'IKr_block_100']
        self.l1_param_labels = [r'$I_{Kr}Block=0\%$', r'$I_{Kr}Block=25\%$',
                                r'$I_{Kr}Block=33\%$', r'$I_{Kr}Block=50\%$',
                                r'$I_{Kr}Block=66\%$', r'$I_{Kr}Block=75\%$',
                                r'$I_{Kr}Block=100\%$']

        self.plotter = gecg.GECGs(self.directory, full=False)
        self.plotter.template.figfractionalwidth = 0.85

    def construct_s2_stimulus_directories(self, n_cases=22, delta=10):
        self.l2_param_dirs = []
        self.l2_param_labels = []
        for ii in range(n_cases):
            s2_dir = "s2_offset_{0:03}".format((delta * ii))
            s2_legend  = 'Tpeak offset = {0:03}ms'.format((delta * ii))
            self.l2_param_dirs.append(s2_dir)
            self.l2_param_labels.append(s2_legend)

    def subplot_l2_comparison(self, BCL):
        ''' for each l1 parameter plot the set of l2 parameters as subplots
        '''
        self.plotter.period_ms = BCL

        # This outer level will drive all of the
        for ii in range(len(self.l1_param_dirs)):
            self.plot_directories = []
            self.plot_legend = []
            l1_path = self.ex_dir + os.sep + "BCL_{0:04}".format(BCL) + "ms" \
            + os.sep + self.l1_param_dirs[ii] + os.sep
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
            self.plotter.set_ECG_type(0, 1)
            self.plotter.full_set_as_subplots(plot_name)

    def overlay_l2_comparison(self, BCL):
        ''' for each l1 parameter plot the set of l2 parameters as subplots
        '''
        self.plotter.period_ms = BCL

        # This outer level will drive all of the
        for ii in range(len(self.l1_param_dirs)):
            self.plot_directories = []
            self.plot_legend = []
            l1_path = self.ex_dir + os.sep + "BCL_{0:04}".format(BCL) + "ms" \
            + os.sep + self.l1_param_dirs[ii] + os.sep
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
            self.plotter.set_ECG_type(0, 1)
            self.plotter.full_ECG_set(plot_name)

    def beat_comparison(self, BCL):
        ''' for each l1 parameter plot the set of l2 parameters as subplots
        '''
        self.plotter.period_ms = BCL
        for ii in range(len(self.l1_param_dirs)):
            self.plot_directories = []
            self.plot_legend = []
            l1_path = self.ex_dir + os.sep + "BCL_{0:04}".format(BCL) + "ms" \
            + os.sep + self.l1_param_dirs[ii] + os.sep
            plot_name = 'BCL_' + repr(BCL) + "_" + self.l1_param_labels[ii]
            # Now we create the plot directories on l2
            for jj in range(len(self.l2_param_dirs)):
                plot_dir = l1_path + os.sep + self.l2_param_dirs[jj]
                self.plot_directories.append(plot_dir)
                self.plot_legend.append(self.l2_param_labels[jj])
            # Construct graph here
            self.plotter.case_directories = self.plot_directories
            self.plotter.case_legends = self.plot_legend
            self.plotter.load_data()
            self.plotter.set_ECG_type(0, 1)
            self.plotter.beat_ECG_Set(1, plot_name + "_beat_0.pdf")
