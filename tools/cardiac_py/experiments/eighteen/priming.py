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
    root_directory/BCL_BCLVAL/l1_param_change/l2_param_change/priming
    ""             ""             ""              ""         /s2_run_*
    '''

    def __init__(self,
                 directory="/Users/butler/Desktop/CARDIOID/ECG_data"):
        '''
        Constructor: We will try to set all parameters in here.
        '''
        self.directory = directory
        self.ex_dir = "experiment_18"
        self.l1_param_dirs = ['ikr_block_0', 'ikr_block_100']
        self.l1_param_labels = ["Control", "E-4031"]
        self.l2_param_dirs = ['g_NaL_M_0_25',
                              'g_NaL_M_0_275',
                              'multi_variable_1',
                              'multi_variable_2',
                              'multi_variable_3',
                              'multi_variable_4']
        '''
                              'multi_variable_5',
                              'multi_variable_6',
                              'multi_variable_7']
        '''
        self.l2_param_labels = [r'\[ \begin{split}g_{NaL, endo} &= 0.20\\ g_{NaL, mid} &= 0.25\\g_{NaL, epi} &= 0.15\\g_{Ks, endo} &= 0.392 \end{split}\]',
                                r'\[ \begin{split}g_{NaL, endo} &= 0.20\\ g_{NaL, mid} &= 0.275\\g_{NaL, epi} &= 0.15\\g_{Ks, endo} &= 0.392 \end{split}\]',
                                r'\[ \begin{split}g_{NaL, endo} &= 0.20\\ g_{NaL, mid} &= 0.325\\g_{NaL, epi} &= 0.15\\g_{Ks, endo} &= 0.392 \end{split}\]',
                                r'\[ \begin{split}g_{NaL, endo} &= 0.20\\ g_{NaL, mid} &= 0.35\\g_{NaL, epi} &= 0.15\\g_{Ks, endo} &= 0.392 \end{split}\]',
                                r'\[ \begin{split}g_{NaL, endo} &= 0.175\\ g_{NaL, mid} &= 0.3\\g_{NaL, epi} &= 0.15\\g_{Ks, endo} &= 0.392 \end{split}\]',
                                r'\[ \begin{split}g_{NaL, endo} &= 0.2\\ g_{NaL, mid} &= 0.3\\g_{NaL, epi} &= 0.15\\g_{Ks, endo} &= 0.196 \end{split}\]']
        '''
                                r'\[ \begin{split}g_{NaL, endo} &= 0.2\\ g_{NaL, mid} &= 0.3\\g_{NaL, epi} &= 0.15\\g_{Ks, endo} &= 0.098 \end{split}\]',
                                r'\[ \begin{split}g_{NaL, endo} &= 0.25\\ g_{NaL, mid} &= 0.3\\g_{NaL, epi} &= 0.15\\g_{Ks, endo} &= 0.098 \end{split}\]',
                                r'\[ \begin{split}g_{NaL, endo} &= 0.3\\ g_{NaL, mid} &= 0.3\\g_{NaL, epi} &= 0.15\\g_{Ks, endo} &= 0.098 \end{split}\]']
        '''
        self.plotter = gecg.GECGs(self.directory)
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
                plot_dir = l1_path + os.sep + self.l2_param_dirs[jj] + os.sep\
                + "priming"
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
                plot_dir = l1_path + os.sep + self.l2_param_dirs[jj] + os.sep\
                + "priming"
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
        for ii in range(len(self.l1_param_dirs)):
            self.plot_directories = []
            self.plot_legend = []
            l1_path = self.ex_dir + os.sep + "BCL_" \
            + repr(BCL) + os.sep + self.l1_param_dirs[ii] + os.sep
            plot_name = 'BCL_' + repr(BCL) + "_" + self.l1_param_labels[ii]
            # Now we create the plot directories on l2
            for jj in range(len(self.l2_param_dirs)):
                plot_dir = l1_path + os.sep + self.l2_param_dirs[jj] + os.sep\
                + "priming"
                self.plot_directories.append(plot_dir)
                self.plot_legend.append(self.l2_param_labels[jj])
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
