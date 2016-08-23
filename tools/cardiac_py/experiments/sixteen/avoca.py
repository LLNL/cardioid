'''
Created on 21/04/2013

@author: butler
'''
import gecg
import os


class S1Cont():
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
        self.ex_16 = "s1_comparison"
        #self.ex_16_ikr_dirs = ['ikr_block_0','ikr_block_33','ikr_block_66', \
        #'ikr_block_100']
        #self.ex_16_legends = [r'$I_{kr}Block = 0 $',r'$I_{kr}Block = 33 $', \
        #r'$I_{kr}Block = 66 [E-4031]$ ',r'$I_{kr}Block = 100 $']
        self.ex_16_ikr_dirs = ['s1_cont', 's1_stop']
        self.ex_16_legends = ["Continued S1", "Stopped S1"]
        self.plotter = gecg.GECGs(self.directory)
        self.plotter.template.figfractionalwidth = 1.15

    def pick_BCL_req_paths(self, BCL):
        self.plotter.period_ms = BCL
        BCL_16_path = self.ex_16
        self.plot_dirs = []
        self.plot_legend = []
        ii = 0
        for ikr_dir in self.ex_16_ikr_dirs:
            ikr_path = BCL_16_path + os.sep + ikr_dir + os.sep
            self.plot_dirs.append(ikr_path)
            self.plot_legend.append(self.ex_16_legends[ii])
            ii = ii + 1

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

    def plot_as_subplots_BCL(self, BCL):
        self.pick_BCL_req_paths(BCL)
        save_name = 's1_s2_comparsion_at_BCL_' + repr(BCL)
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
