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

    def __init__(self, direct="/Users/butler/Desktop/CARDIOID_SEQ/ecg_data"):
        '''
        Constructor
        '''
        self.directory = direct
        self.ex_16 = "experiment_16"
        #self.ex_16_ikr_dirs = ['ikr_block_0','ikr_block_33','ikr_block_66',
        # 'ikr_block_100']
        #self.ex_16_legends = [r'$I_{kr}Block = 0 $',r'$I_{kr}Block = 33 $', \
        #r'$I_{kr}Block = 66 [E-4031]$ ',r'$I_{kr}Block = 100 $']
        self.ex_16_ikr_dirs = ['ikr_block_0', 'ikr_block_66']
        self.ex_16_legends = ["Control", "E-4031"]
        self.ex_17 = "experiment_17"
        self.ex_17_ikr_dirs = ['ikr_block_66']
        self.ex_17_legends = ['Ranolazine']
        self.plotter = gecg.GECGs(self.directory)
        self.plotter.template.figfractionalwidth = 0.85

    def pick_BCL_req_paths(self, BCL):
        self.plotter.period_ms = BCL
        BCL_16_path = self.ex_16
        BCL_17_path = self.ex_17
        self.plot_dirs = []
        self.plot_legend = []
        ii = 0
        for ikr_dir in self.ex_16_ikr_dirs:
            ikr_path = BCL_16_path + os.sep + ikr_dir + os.sep + "priming"
            self.plot_dirs.append(ikr_path)
            self.plot_legend.append(self.ex_16_legends[ii])
            ii = ii + 1

        for ikr_dir in self.ex_17_ikr_dirs:
            ikr_path = BCL_17_path + os.sep + ikr_dir + os.sep + "priming"
            self.plot_dirs.append(ikr_path)
            self.plot_legend.append(self.ex_17_legends[0])

    def plot_full_as_subplots(self):
        self.plotter.case_directories = self.plot_dirs
        self.plotter.case_legends = self.plot_legend
        self.plotter.load_data()
        self.plotter.set_ECG_type(0, -1)
        self.plotter.full_set_as_subplots('test_subplots.png')

    def plot_beat_overlay(self):
        self.plotter.case_directories = self.plot_dirs
        self.plotter.case_legends = self.plot_legend
        self.plotter.load_data()
        self.plotter.set_ECG_type(0, -1)
        self.plotter.beat_ECG_Set(10, 'One_pacing_beat_10.pdf')
        self.plotter.beat_ECG_Set(9, 'One_pacing_beat_9.pdf')
