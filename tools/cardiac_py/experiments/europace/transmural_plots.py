'''
Created on Sep 11, 2012

@author: butler
'''
import numpy as np
import plot_schema
import matplotlib.pyplot as plt


class TransmuralPlots():
    '''
    classdocs
    '''
    def __init__(self, directory):
        '''
        Constructor
        '''
        self.t_gap_ms = 5.0
        self.directory = directory
        self.normal_data_name = "normal_ikr.dat"
        self.modified_data_name = 'zeroed_ikr.dat'
        self.period_ms = 1000
        self.template = plot_schema.PlotSchema()

    def load_data(self):
        n_data_path = self.directory + '/' + self.normal_data_name
        self.n_data = np.loadtxt(n_data_path)
        m_data_path = self.directory + '/' + self.modified_data_name
        self.m_data = np.loadtxt(m_data_path)

    def set_ECG_type(self, ECG_type, flipper):
        ''' SETUP plot vectors for each of the different ECG types
        '''
        n_points = self.n_data[:, 0].size
        # 1 setup time series
        max_time = self.t_gap_ms / 1000 * n_points
        self.time = np.linspace(0, max_time, n_points)
        # first normalise
        if ECG_type == 0:
            pass
            self.n_ECG_data = flipper * (self.n_data[:, 0] - self.n_data[:, 1])
            self.m_ECG_data = flipper * (self.m_data[:, 0] - self.m_data[:, 1])
            self.y_label = r'$ \Delta V$'
        elif ECG_type == 1:
            self.n_ECG_data = flipper * self.n_data[:, 0]
            self.m_ECG_data = flipper * self.m_data[:, 0]
            self.y_label = r'$ V$'
        elif ECG_type == 2:
            self.n_ECG_data = flipper * self.n_data[:, 1]
            self.m_ECG_data = flipper * self.m_data[:, 1]
            self.y_label = r'$V$'
        else:
            print "ECG type is incorrect"
            return

    def plot_normal_ECG_final(self, save_file):

        """ """
        index_end = self.time.size - 1
        index_start = self.time.size - 1 - (self.period_ms / self.t_gap_ms)
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.time[index_start:index_end],
                 self.n_ECG_data[index_start:index_end])
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=100, bbox_inches='tight')
        plt.show()

    def plot_normal_ECG_full(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.time, self.n_ECG_data)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=100, bbox_inches='tight')
        plt.show()

    def plot_modified_ECG_full(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.time, self.m_ECG_data)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=100, bbox_inches='tight')
        plt.show()

    def plot_modified_ECG_final(self, save_file):
        index_end = self.time.size - 1
        index_start = self.time.size - 1 - (self.period_ms / self.t_gap_ms)
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.time[index_start:index_end],
                 self.m_ECG_data[index_start:index_end])
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=100, bbox_inches='tight')
        plt.show()

    def overlay_ECG_full(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.time, self.n_ECG_data)
        plt.plot(self.time, self.m_ECG_data)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=100, bbox_inches='tight')
        plt.show()

    def overlay_ECG_final(self, save_file):
        index_end = self.time.size - 1
        index_start = self.time.size - 1 - (self.period_ms / self.t_gap_ms)
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.time[index_start:index_end],
                 self.n_ECG_data[index_start:index_end])
        plt.plot(self.time[index_start:index_end],
                 self.m_ECG_data[index_start:index_end])
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=100, bbox_inches='tight')
        plt.show()
