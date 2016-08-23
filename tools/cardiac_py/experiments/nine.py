'''
Created on 08/01/2013

@author: butler
'''
import numpy as np
import plot_schema
import matplotlib.pyplot as plt


class PlotsNine():
    '''
    classdocs
    '''
    def __init__(self, directory, full=True):
        '''
        Constructor
        '''
        self.t_gap_ms = 5.0
        self.directory = directory
        self.full = full
        if self.full:
            self.electrodes = ['electrode#000448302', 'electrode#000451300',
                               'electrode#000452730', 'electrode#000453393',
                               'electrode#000457525', 'electrode#000458894',
                               'electrode#000438028', 'electrode#000460291']
        else:
            self.electrodes = ['electrode#000094150', 'electrode#000092294']
        self.period_ms = 1000
        self.template = plot_schema.PlotSchema()

    def load_data(self):
        n_data_path = self.directory
        self.n_data = self.load_case_data(n_data_path)

    def load_case_data(self, case_path):
        n_electrodes = len(self.electrodes)
        for ii in range(n_electrodes):
            electrode_path = case_path + '/' + self.electrodes[ii]
            temp_data = np.loadtxt(electrode_path)
            if ii == 0:
                n_times = temp_data.size
                data = np.zeros((n_times, n_electrodes))
            data[:, ii] = temp_data[:]
        return data

    def set_ECG_type(self, ECG_lead, flipper=1):
        ''' SETUP plot vectors for each of the different ECG types
        '''
        if self.full:
            lead_1 = 7
            lead_2 = 8
        else:
            lead_1 = 1
            lead_2 = 2
        ECG_type = ECG_lead - 1
        assert((len(self.electrodes) >= ECG_type))
        n_points = self.n_data[:, 0].size
        # 1 setup time series
        max_time = self.t_gap_ms / 1000 * n_points
        self.time = np.linspace(0, max_time, n_points)
        # first normalise
        if ECG_type == -1:
            col_1 = lead_1 - 1
            col_2 = lead_2 - 1
            self.n_ECG_data = flipper * (self.n_data[:, col_1] -
                                         self.n_data[:, col_2])
            self.y_label = r'$ \Delta V$'
        else:
            self.n_ECG_data = flipper * self.n_data[:, ECG_type]
            self.y_label = r'$ V$'

    def plot_normal_ECG_final(self, save_file):
        """ """
        index_end = self.time.size / 2
        index_start = 0
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.grid()
        plt.plot(self.time[index_start:index_end],
                 self.n_ECG_data[index_start:index_end])
        plt.xlabel('$t (s)$', fontsize=14)
        plt.ylabel(self.y_label, fontsize=14, rotation='horizontal')
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(14)
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(14)
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=100, bbox_inches='tight')

    def plot_normal_ECG_full(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.grid()
        plt.plot(self.time, self.n_ECG_data)
        plt.xlabel('$t (s)$', fontsize=14)
        plt.ylabel(self.y_label, fontsize=14, rotation='horizontal')
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(14)
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(14)
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=100, bbox_inches='tight')
        plt.show()
