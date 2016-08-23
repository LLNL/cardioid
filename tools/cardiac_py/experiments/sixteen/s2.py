'''
Created on 25/02/2013

@author: butler
'''
import gecg
import os


class Ikr():
    '''
    This class produces plots which compare priming beats across multiple ikr
    for a single BCL
    '''

    def __init__(self,
                 directory="/Users/butler/Desktop/CARDIOID_SEQ/ecg_data"):
        '''
        Constructor
        '''
        self.directory = directory
        self.ex_16 = "experiment_16"
        #self.ex_16_ikr_dirs = ['ikr_block_0','ikr_block_33','ikr_block_66', \
        #'ikr_block_100']
        #self.ex_16_legends = [r'$I_{kr}Block = 0 $',r'$I_{kr}Block = 33 $', \
        #r'$I_{kr}Block = 66 [E-4031]$ ',r'$I_{kr}Block = 100 $']
        self.ex_16_ikr_dirs = ['ikr_block_0', 'ikr_block_66']
        self.ex_16_legends = ["Control", "E-4031"]
        self.ex_17 = "experiment_17"
        self.ex_17_ikr_dirs = ['ikr_block_66']
        self.ex_17_legends = ['Ranolazine']
        self.plotter = gecg.GECGs(self.directory)
        self.plotter.template.figfractionalwidth = 1.15

    def pick_BCL_req_paths(self, BCL, s2_run):
        self.plotter.period_ms = BCL
        BCL_16_path = self.ex_16 + os.sep + "BCL_" + repr(BCL)
        BCL_17_path = self.ex_17 + os.sep + "BCL_" + repr(BCL)
        self.plot_dirs = []
        self.plot_legend = []
        ii = 0
        for ikr_dir in self.ex_16_ikr_dirs:
            ikr_path = BCL_16_path + os.sep + ikr_dir + os.sep + "s2_run_" + \
            repr(s2_run)
            self.plot_dirs.append(ikr_path)
            self.plot_legend.append(self.ex_16_legends[ii])
            ii = ii + 1

        for ikr_dir in self.ex_17_ikr_dirs:
            ikr_path = BCL_17_path + os.sep + ikr_dir + os.sep + "s2_run_" + \
            repr(s2_run)
            self.plot_dirs.append(ikr_path)
            self.plot_legend.append(self.ex_17_legends[0])

    def plot_full_as_subplots(self, save_name):
        self.plotter.case_directories = self.plot_dirs
        self.plotter.case_legends = self.plot_legend
        self.plotter.load_data()
        self.plotter.set_ECG_type(0, -1)
        self.plotter.full_set_as_subplots(save_name)

    def plot_beat_overlay(self):
        self.plotter.case_directories = self.plot_dirs
        self.plotter.case_legends = self.plot_legend
        self.plotter.load_data()
        self.plotter.set_ECG_type(0, -1)
        self.plotter.beat_ECG_Set(10, 'One_pacing_beat_10.pdf')
        self.plotter.beat_ECG_Set(9, 'One_pacing_beat_9.pdf')

    def plot_frac_ECG_set(self, frac, save_file):
        try:
            self.plotter.case_directories = self.plot_dirs
            self.plotter.case_legends = self.plot_legend
            self.plotter.load_data()
            self.plotter.set_ECG_type(0, -1)
            self.plotter.frac_ECG_set(frac, save_file)
        except IOError:
            pass

    def plot_all_s2_For_BCL(self, BCL, frac):
        for s2_run in range(7, 13):
            s2_run = s2_run + 1
            self.pick_BCL_req_paths(BCL, s2_run)
            save_name = 'BCL_' + repr(BCL) + '_s2_run_' + repr(s2_run) + '.pdf'
            self.plot_frac_ECG_set(frac, save_name)

    def plot_all_s2_subplots_BCL(self, BCL, frac):
        for s2_run in range(13):
            s2_run = s2_run + 1
            self.pick_BCL_req_paths(BCL, s2_run)
            save_name = 'SPlot_BCL_' + repr(BCL) + '_s2_run_' + repr(s2_run) +\
            '.pdf'
            self.plot_full_as_subplots(save_name)

    def plot_all_leads(self, ikr_case, save_file):
        self.plotter.case_directories = self.plot_dirs
        self.plotter.case_legends = self.plot_legend
        self.plotter.load_data()
        self.plotter.set_ECG_type(0, -1)
        self.plotter.all_leads_subplot(0, save_file)

    def get_all_lead_csv(self, ikr_case, save_file):
        self.plotter.case_directories = self.plot_dirs
        self.plotter.case_legends = self.plot_legend
        self.plotter.load_data()
        self.plotter.set_ECG_type(0, -1)
        self.plotter.write_csv_all_leads(0, save_file)


class s1_variation_tests():
    """
    This documents a small series of sensitivity tests:
    1) Testing 5ms offsets for
    ./experiment_16/BCL_2000/ikr_block_0/s2_run_3
    2) Testing effect of only one (of two) reentry stimuli on
    ./experiment_16/BCL_2000/ikr_block_0/s2_run_3
    3) Testing effect of s1 continuation on TdP produced for
    ./experiment_16/BCL_2000/ikr_block_0/s2_run_3
    ./experiment_16/BCL_2000/ikr_block_0/s2_run_5
    ./experiment_16/BCL_2000/ikr_block_0/s2_run_6
    ./experiment_16/BCL_2000/ikr_block_0/s2_run_7
    ./experiment_16/BCL_2000/ikr_block_66/s2_run_13
    4) Test highly likely E-4031 arrhythmia
    ./experiment_16/BCL_2000/ikr_block_66/s2_run_14
    """

    def __init__(self,
                 directory="/Users/butler/Desktop/CARDIOID_SEQ/ecg_data"):
        '''
        Constructor
        '''
        self.directory = directory
        self.plotter = gecg.GECGs(self.directory)
        self.plotter.template.figfractionalwidth = 1.15

    def TdP_case_s1_r1_only(self):
        pass

    def TdP_case_s1_r2_only(self):
        pass

    def TdP_case_minus_5ms(self):
        pass

    def TdP_case_plus_5ms(self):
        pass

    def potential_s1_interference_1(self):
        normal_dir = "experiment_16/BCL_2000/ikr_block_0/s2_run_3"
        altered_dir = "experiment_16/BCL_2000/ikr_block_0/s2_run_3_no_s1"
        normal_legend = "With s1"
        altered_legend = "Without s1"
        self.plotter.case_directories = [normal_dir, altered_dir]
        self.plotter.case_legends = [normal_legend, altered_legend]
        self.plotter.load_data()
        self.plotter.set_ECG_type(0, -1)
        self.plotter.frac_ECG_set(1.0,
                                  "BCL_2000_ikr_block_0_s2_run_3_no_s1.pdf")

    def potential_s1_interference_2(self):
        pass

    def potential_s1_interference_3(self):
        pass

    def potential_s1_interference_4(self):
        pass

    def potential_s1_interference_5(self):
        pass

    def extra_case(self):
        pass
