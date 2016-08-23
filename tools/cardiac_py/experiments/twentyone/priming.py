'''
Created on 29/04/2013

@author: butler
'''
import gecg
import os

# FIXME: Some of the paths need to be corrected.
class Priming():
    '''
    Post processing priming information for experiment 18 parameter space with
    an eye cast forward to a 'generalizable' case in the future.
    Experiment directory structure:
    root_directory//l1_param_change/BCL_BCLVAL/l2_param_change/priming
    ""             ""             ""              ""         /s2_run_*
    '''
    

    def __init__(self,
                 directory="/Users/butler/Desktop/CARDIOID_SEQ/ECG_data"):
        '''
        Constructor: We will try to set all parameters in here.
        '''
        self.directory = directory
        self.ex_dir = "experiment_21"
        self.l1_param_dirs = ['thick_apex']
        self.l1_param_labels = ["thick_apex"]
        self.l2_param_dirs = ['gKr_0.102',
                              'gKr_0.153']

        self.l2_param_labels = [r'33 \% gKr block',
                                r'Control']

        self.plotter = gecg.GECGs(self.directory)
        self.plotter.template.figfractionalwidth = 0.85
        
    def set_l1_l2(self, BCL):
        '''
        Optional function to over-ride behaviour based on what the BCL is. In
        the future it would be nice to do something different. Like specify what 
        BCL's are acceptable and update.
        TODO: rework to the afformentioned idea
        TODO: Use inheretance more. Perhaps we need a standard prototype class
        which would work with common viz patterns to minimize the codebase further.
        '''
        if BCL == 500:
            self.l2_param_dirs = ["gKr_0.153"]
            self.l2_param_labels = ["Control"]
        elif BCL == 1000:
            self.l2_param_dirs = ["gKr_0.153","gKr_0.0","gKr_0.153_gNaL_0.0"]
            self.l2_param_labels = ["Control","100\% gKr block","TT06 standard"]
        elif BCL == 2000:
            self.l2_param_dirs = ["gKr_0.153", "gKr_0.153_r1685"]
            self.l2_param_labels = ["Control","Control updated"]
        else:
            pass
        
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
        self.set_l1_l2(BCL)
        for ii in range(len(self.l1_param_dirs)):
            self.plot_directories = []
            self.plot_legend = []
            l1_path = self.ex_dir + os.sep + self.l1_param_dirs[ii] + os.sep + \
              "BCL_" + repr(BCL) + os.sep
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
            self.plotter.beat_ECG_Set(9, plot_name + "_beat_lead_1.pdf")
            self.plotter.set_ECG_type(1, 1)
            self.plotter.beat_ECG_Set(9, plot_name + "_beat_lead_v1.pdf")
            self.plotter.set_ECG_type(2, 1)
            self.plotter.beat_ECG_Set(9, plot_name + "_beat_lead_v2.pdf")
            self.plotter.set_ECG_type(3, 1)
            self.plotter.beat_ECG_Set(9, plot_name + "_beat_lead_v3.pdf")
            self.plotter.set_ECG_type(4, 1)
            self.plotter.beat_ECG_Set(9, plot_name + "_beat_lead_v4.pdf")
            self.plotter.set_ECG_type(5, 1)
            self.plotter.beat_ECG_Set(9, plot_name + "_beat_lead_v5.pdf")
            self.plotter.set_ECG_type(6, 1)
            self.plotter.beat_ECG_Set(9, plot_name + "_beat_lead_v6.pdf")
