'''
Created on 08/01/2013

@author: butler
'''
import numpy as np
import plot_schema
import matplotlib.pyplot as plt


class PlotECGs():
    '''
    classdocs
    '''
    def __init__(self, directory, full=True):
        '''
        Constructor
        '''
        self.t_gap_ms = 5.0
        self.directory = directory
        self.normal_dir_name = "/ECG_OLD"
        self.modified_dir_name = '/ECG_NEW'
        self.normal_legend = 'OLD ECG'
        self.modified_legend = 'NEW ECG'
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
        self.lfont = matplotlib.font_manager.FontProperties(size=(0.75 * self.template.fontsize))

    def load_data(self):
        n_data_path = self.directory + '/' + self.normal_dir_name
        self.n_data = self.load_case_data(n_data_path)
        m_data_path = self.directory + '/' + self.modified_dir_name
        self.m_data = self.load_case_data(m_data_path)

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
        n_points_n = self.n_data[:, 0].size
        n_points_m = self.m_data[:, 0].size
        # 1 setup time series
        max_time_n = self.t_gap_ms / 1000 * n_points_n
        max_time_m = self.t_gap_ms / 1000 * n_points_m
        self.time_n = np.linspace(0, max_time_n, n_points_n)
        self.time_m = np.linspace(0, max_time_m, n_points_m)
        # first normalise
        if ECG_type == -1:
            col_1 = lead_1 - 1
            col_2 = lead_2 - 1
            self.n_ECG_data = flipper * (self.n_data[:, col_1] -
                                         self.n_data[:, col_2])
            self.m_ECG_data = flipper * (self.m_data[:, col_1] -
                                         self.m_data[:, col_2])
            self.y_label = r'$ \Delta V$'
        else:
            self.n_ECG_data = flipper * self.n_data[:, ECG_type]
            self.m_ECG_data = flipper * self.m_data[:, ECG_type]
            self.y_label = r'$ V$'

    def plot_normal_ECG_final(self, save_file):
        """ Hello world"""
        index_end = self.time.size - 1
        index_start = self.time.size - 1 - (self.period_ms / self.t_gap_ms)
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.time_n[index_start:index_end],
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
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')
        plt.show()

    def plot_normal_ECG_full(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.time_n, self.n_ECG_data)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')
        plt.show()

    def plot_modified_ECG_full(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.time_m, self.m_ECG_data)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')
        plt.show()

    def plot_modified_ECG_final(self, save_file):
        index_end = self.time.size - 1
        index_start = self.time.size - 1 - (self.period_ms / self.t_gap_ms)
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.time_m[index_start:index_end],
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
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')
        plt.show()

    def overlay_ECG_full(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.time_n, self.n_ECG_data, label=self.normal_legend)
        plt.plot(self.time_m, self.m_ECG_data, label=self.modified_legend)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        plt.legend(prop=self.lfont)
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')
        plt.show()

    def overlay_ECG_beat(self, save_file, beat):
        n_beats = self.time_n.max() / (self.period_ms / 1000)
        assert(beat < n_beats)
        index_start = (self.period_ms / self.t_gap_ms) * (beat - 1)
        index_end = (self.period_ms / self.t_gap_ms) * (beat)
        if beat > 1:
            index_start = index_start - 1
            index_end = index_end - 1
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.grid()
        plt.plot(self.time_n[index_start:index_end],
                 self.n_ECG_data[index_start:index_end],
                 label=self.normal_legend)
        plt.plot(self.time_m[index_start:index_end],
                 self.m_ECG_data[index_start:index_end],
                 label=self.modified_legend)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        plt.legend(prop=self.lfont)
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')

    def overlay_ECG_final(self, save_file):
        index_end = self.time_n.size - 1
        index_start = self.time_n.size - 1 - (self.period_ms / self.t_gap_ms)
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.time_n[index_start:index_end],
                 self.n_ECG_data[index_start:index_end],
                 label=self.normal_legend)
        plt.plot(self.time_m[index_start:index_end],
                 self.m_ECG_data[index_start:index_end],
                 label=self.modified_legend)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        plt.legend(prop=self.lfont)
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')

    def overlay_ECG_rapid(self, save_file, colour='black', l1=False):
        index_end = self.time.size - 1
        index_start = self.time.size - 1 - (self.period_ms / self.t_gap_ms)
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        if self.fix_axes:
            axes.set_ylim((self.e_min, self.e_max))
            axes.set_xlim((self.time[index_start]), self.time[index_end])
        axes.set_frame_on(False)
        plt.plot(self.time[index_start:index_end],
                 self.m_ECG_data[index_start:index_end],
                 label=self.modified_legend, color=colour, linestyle="--",
                 linewidth=2)
        plt.plot(self.time[index_start:index_end],
                 self.n_ECG_data[index_start:index_end],
                 label=self.normal_legend, color=colour, linestyle="-",
                 linewidth=2)
        axes.get_xaxis().set_visible(l1)
        axes.get_yaxis().set_visible(l1)

        self.template.apply_figuresize_settings(f)
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight', transparent=True)

    def overlay_plus_difference(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.subplot(211)
        plt.plot(self.time_n, self.n_ECG_data, label=self.normal_legend)
        plt.plot(self.time_m, self.m_ECG_data, label=self.modified_legend)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        plt.legend(prop=self.lfont)
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        plt.subplot(212)
        # First we must check that vectors are the same size.
        if self.n_ECG_data.size == self.m_ECG_data.size:
            # simple case
            ECG_diff_mag = np.abs(self.n_ECG_data[:] - self.m_ECG_data)
            time = self.time_n
        elif self.n_ECG_data.size < self.m_ECG_data.size:
            dat_size = self.n_ECG_data.size
            ECG_diff_mag = np.abs(self.n_ECG_data[:] -
                                  self.m_ECG_data[:dat_size])
            time = self.time_n
        else:
            dat_size = self.m_ECG_data.size
            ECG_diff_mag = np.abs(self.n_ECG_data[:dat_size] -
                                  self.m_ECG_data[:])
            time = self.time_m
        plt.plot(time, ECG_diff_mag)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel('$\Delta$', fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')
        plt.show()

    def overlay_plus_perc_error(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.subplot(211)
        plt.plot(self.time_n, self.n_ECG_data, label=self.normal_legend)
        plt.plot(self.time_m, self.m_ECG_data, label=self.modified_legend)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        plt.legend(prop=self.lfont)
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        plt.subplot(212)
        # First we must check that vectors are the same size.
        if self.n_ECG_data.size == self.m_ECG_data.size:
            # simple case
            ECG_diff_mag = np.abs((self.n_ECG_data[:] - self.m_ECG_data) /
                                  self.n_ECG_data.max() * 100)
            time = self.time_n
        elif self.n_ECG_data.size < self.m_ECG_data.size:
            dat_size = self.n_ECG_data.size
            ECG_diff_mag = np.abs((self.n_ECG_data[:] -
                                   self.m_ECG_data[:dat_size]) /
                                  self.n_ECG_data.max() * 100)
            time = self.time
        else:
            dat_size = self.m_ECG_data.size
            ECG_diff_mag = np.abs((self.n_ECG_data[:dat_size] -
                                   self.m_ECG_data[:]) /
                                  self.n_ECG_data.max() * 100)
            time = self.time_m
        plt.semilogy(time, ECG_diff_mag)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel('$\Delta$', fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')

    def overlay_ECG_reent(self, save_file, beat):
        n_beats = self.time_n.max() / (self.period_ms / 1000)
        assert(beat < n_beats)
        index_start = (self.period_ms / self.t_gap_ms) * (9.22 - 1)
        index_end = (self.period_ms / self.t_gap_ms) * (9.28 - 1)
        if beat > 1:
            index_start = index_start - 1
            index_end = index_end - 1
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.grid()
        plt.plot(self.time_n[index_start:index_end],
                 self.n_ECG_data[index_start:index_end],
                 label=self.normal_legend)
        plt.plot(self.time_m[index_start:index_end],
                 self.m_ECG_data[index_start:index_end],
                 label=self.modified_legend)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        plt.legend(prop=self.lfont)
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')

    def write_normal_lead_csv(self, save_file_name):
        #self.time_n, self.n_ECG_data
        a = np.zeros((self.time_n.size, 2))
        a[:, 0] = self.time_n[:]
        a[:, 1] = self.n_ECG_data[:]
        np.savetxt(save_file_name, a, delimiter=",")

    def overlay_plus_perc_error_beat(self, save_file, beat):
        n_beats = self.time_n.max() / (self.period_ms / 1000)
        assert(beat < n_beats)
        index_start = (self.period_ms / self.t_gap_ms) * (beat - 1)
        index_end = (self.period_ms / self.t_gap_ms) * (beat)
        if beat > 1:
            index_start = index_start - 1
            index_end = index_end - 1
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.subplot(211)
        plt.plot(self.time_n[index_start:index_end],
                 self.n_ECG_data[index_start:index_end],
                 label=self.normal_legend)
        plt.plot(self.time_m[index_start:index_end],
                 self.m_ECG_data[index_start:index_end],
                 label=self.modified_legend)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        plt.legend(prop=self.lfont)
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        plt.subplot(212)
        ECG_diff_mag = np.abs((self.n_ECG_data[index_start:index_end] -
                               self.m_ECG_data[index_start:index_end]) /
                              np.abs(self.n_ECG_data).max())
        time = self.time_m[index_start:index_end]
        plt.semilogy(time, ECG_diff_mag)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel('$\Delta$', fontsize=self.template.get_fontsize(),
                   rotation='horizontal')

        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')
