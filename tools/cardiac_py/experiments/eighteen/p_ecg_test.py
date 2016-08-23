'''
Created on 05/05/2013

@author: butler
'''
import numpy as np
import plot_schema
import matplotlib.pyplot as plt
import os
import glob
from matplotlib.ticker import MultipleLocator


class Compare_parallel_v_serial():
    '''
    This class is a reduction of ../gecgs module. Here the idea is to compare a
    case were we have both
    '''
    def __init__(self, directory, full=True):
        '''
        Constructor
        '''
        self.t_gap_ms = 5.0
        self.directory = directory
        self.full = full
        self.case_legends = ["serial", "parallel"]
        if self.full:
            self.electrodes = ['electrode#000448302', 'electrode#000451300',
                               'electrode#000452730', 'electrode#000453393',
                               'electrode#000457525', 'electrode#000458894',
                               "electrode#000438028", "electrode#000460291"]
        else:
            self.electrodes = ['electrode#000094150', 'electrode#000092294']
        self.period_ms = 1000
        self.template = plot_schema.PlotSchema()
        self.lfont = matplotlib.font_manager.FontProperties(size=(0.75 * self.template.fontsize))
    def load_data(self):
        n_cases = 2
        self.stat_flag = np.zeros((n_cases))
        print 'no of cases: ', repr(n_cases)
        data_set_index = 0
        for ii in range(n_cases):
            if ii == 0:
                case_data = self.load_case_data_serial(self.directory)
            else:
                case_data = self.load_case_data_serial(self.directory)
            if ii == data_set_index:
                n_steps, n_electrodes = case_data.shape
                if n_steps == 1:
                    print 'Ohoh'
                    self.stat_flag[ii] = 0
                    data_set_index = data_set_index + 1
                else:
                    self.data = np.zeros((n_steps, n_electrodes, n_cases))
                    print 'n tsteps: ', n_steps
                    self.stat_flag[ii] = 1
                self.data[:, :, ii] = case_data
            else:
                n_steps_c, n_electrodes = case_data.shape
                n_steps, n_electrodes, n_cases = self.data.shape
                print 'n tsteps: ', n_steps_c
                if n_steps_c == n_steps:
                    self.data[:, :, ii] = case_data
                    self.stat_flag[ii] = 1
                else:
                    self.data[:n_steps_c, :, ii] = case_data[:, :]
                    self.stat_flag[ii] = 1
        print "data loaded"

    def load_case_data_parallel(self, case_path):
        """ Work with parallel data sets """
        n_electrodes = len(self.electrodes)
        try:
            for ii in range(n_electrodes):
                #
                electrode_path_glob = case_path + '/' + self.electrodes[ii] + \
                ".*"
                electrode_file_list = glob.glob(electrode_path_glob)
                electrode_file_list.sort()
                convertfunc = lambda x: float(x.strip("snapshot."))
                jj = 0
                for part_electode in electrode_file_list:
                    my_temp_data = np.genfromtxt(part_electode, delimiter="",
                                            converters={0: convertfunc})
                    if jj == 0:
                        temp_data = my_temp_data
                    else:
                        temp_data = np.concatenate((temp_data, my_temp_data),
                                                   axis=0)
                    jj = jj + 1
                # In place sorts are bad...bad bad. make sure you don't
                #re-assign while sorting.
                temp_data_2 = temp_data[temp_data[:, 0].argsort(), :]
                print 'temp_data_2'
                print temp_data_2
                temp_data_3 = temp_data_2[:, 1]
                print 'temp_data_3'
                print temp_data_3
                if ii == 0:
                    n_times = temp_data_3.size
                    data = np.zeros((n_times, n_electrodes))
                data[:, ii] = temp_data_3[:]
            return data
        except IOError:
            print 'IOError case: ', case_path
            return np.zeros((1, 1))
        except ValueError:
            print 'ValueError case: ', case_path
            print 'Electrode: ', self.electrodes[ii]
            return np.zeros((1, 1))

    def load_case_data_serial(self, case_path):
        n_electrodes = len(self.electrodes)
        try:
            for ii in range(n_electrodes):
                electrode_path = case_path + '/' + self.electrodes[ii]
                temp_data = np.loadtxt(electrode_path)
                if ii == 0:
                    n_times = temp_data.size
                    data = np.zeros((n_times, n_electrodes))
                    data[:, ii] = temp_data[:]
                else:
                    if temp_data.size == n_times:
                        data[:, ii] = temp_data[:]
                    elif temp_data.size < n_times:
                        data[:temp_data.size, ii] = temp_data[:]
                    elif temp_data.size > n_times:
                        data[:, ii] = temp_data[:n_times]
            return data
        except IOError:
            print 'IOError case: ', case_path
            return np.zeros((1, 1))
        except ValueError:
            print 'ValueError case: ', case_path
            return np.zeros((1, 1))

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

        self.t_gap_s = self.t_gap_ms / 1000
        n_points = self.data[:, 0, 0].size
        self.max_time = self.t_gap_s * (n_points)
        self.time = np.linspace(self.t_gap_s, self.max_time, n_points)
        # first normalise
        if ECG_type == -1:
            col_1 = lead_1 - 1
            col_2 = lead_2 - 1
            self.ECG_data = flipper * (self.data[:, col_1, :] - \
                                       self.data[:, col_2, :])
            self.y_label = r'$ \Delta V [mV]$'
        else:
            self.ECG_data = flipper * self.data[:, ECG_type, :]
            self.y_label = r'$ V [mV]$'
        (rows, cols) = self.ECG_data.shape
        for jj in range(rows):
            for ii in range(cols):
                if np.abs(self.ECG_data[jj, ii]) > 10:
                    self.ECG_data[jj, ii] = 0.0

    def full_ECG_set(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        for case in range(2):
            plt.plot(self.time, self.ECG_data[:, case],
                    label=self.case_legends[case])
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        plt.legend(prop=self.lfont)
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        axes.set_aspect(2.5)
        #axes.xaxis.set_major_locator(MultipleLocator(0.2))
        #axes.xaxis.set_minor_locator(MultipleLocator(0.04))
        #axes.yaxis.set_major_locator(MultipleLocator(0.5))
        #axes.yaxis.set_minor_locator(MultipleLocator(0.1))
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')
        plt.show()

    def beat_ECG_Set(self, beat, save_file):
        n_beats = self.time.max() / (self.period_ms / 1000)
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
        plt.grid(which='both')
        for case in range(2):
            plt.plot(self.time[index_start:index_end],
                    self.ECG_data[index_start:index_end, case],
                    label=self.case_legends[case])
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        axes.set_aspect(2.5)
        axes.xaxis.set_major_locator(MultipleLocator(0.2))
        axes.xaxis.set_minor_locator(MultipleLocator(0.04))
        axes.yaxis.set_major_locator(MultipleLocator(0.5))
        axes.yaxis.set_minor_locator(MultipleLocator(0.1))
        plt.legend(prop=self.lfont)
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')

    def full_set_as_subplots(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        f.subplots_adjust(wspace=0.3)
        self.template.apply_figuresize_settings(f)
        n_plots = 2
        for case in range(n_plots):
            plt.subplot(n_plots, 1, (case + 1))
            plt.plot(self.time, self.ECG_data[:, case],
                     label=self.case_legends[case])
            plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                       rotation='horizontal')
            plt.legend(prop=self.lfont)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        self.template.apply_figuresize_settings(f)
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')

    def beat_set_as_subplots(self, beat, save_file):
        n_beats = self.time.max() / (self.period_ms / 1000)
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
        n_plots = len(self.case_diretories)
        for case in range(n_plots):
            axes = f.add_subplot(n_plots, 1, (case + 1))
            axes.plot(self.time[index_start:index_end],
                      self.ECG_data[index_start:index_end, case],
                      label=self.case_legends[case])
            plt.legend(prop=self.lfont)
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


class Compare_serial_versions():
    '''
    This class is a reduction of ../gecgs module. Here the idea is to compare a
    case were we have both
    '''
    def __init__(self, directory, full=True):
        '''
        Constructor
        '''
        self.t_gap_ms = 5.0
        self.directory = directory
        self.full = full
        self.case_legends = ["old gcc", "XL"]
        if self.full:
            self.electrodes = ['electrode#000448302', 'electrode#000451300',
                               'electrode#000452730', 'electrode#000453393',
                               'electrode#000457525', 'electrode#000458894',
                               "electrode#000438028", "electrode#000460291"]
        else:
            self.electrodes = ['electrode#000094150', 'electrode#000092294']
        self.period_ms = 1000
        self.template = plot_schema.PlotSchema()

    def load_data(self):
        n_cases = 2
        self.stat_flag = np.zeros((n_cases))
        print 'no of cases: ', repr(n_cases)
        data_set_index = 0
        for ii in range(n_cases):
            if ii == 0:
                case_data = self.load_case_data_serial(self.directory + os.sep
                                                       + 'ecg_serial_data')
            else:
                case_data = self.load_case_data_serial(self.directory)
            if ii == data_set_index:
                n_steps, n_electrodes = case_data.shape
                if n_steps == 1:
                    self.stat_flag[ii] = 0
                    data_set_index = data_set_index + 1
                else:
                    self.data = np.zeros((n_steps, n_electrodes, n_cases))
                    print 'n tsteps: ', n_steps
                    self.stat_flag[ii] = 1
                self.data[:, :, ii] = case_data
            else:
                n_steps_c, n_electrodes = case_data.shape
                n_steps, n_electrodes, n_cases = self.data.shape
                print 'n tsteps: ', n_steps_c
                if n_steps_c == n_steps:
                    self.data[:, :, ii] = case_data
                    self.stat_flag[ii] = 1
                else:
                    self.data[:n_steps_c, :, ii] = case_data[:, :]
                    self.stat_flag[ii] = 1
        print "data loaded"

    def load_case_data_serial(self, case_path):
        n_electrodes = len(self.electrodes)
        try:
            for ii in range(n_electrodes):
                electrode_path = case_path + '/' + self.electrodes[ii]
                temp_data = np.loadtxt(electrode_path)
                if ii == 0:
                    n_times = temp_data.size
                    data = np.zeros((n_times, n_electrodes))
                data[:, ii] = temp_data[:]
            return data
        except IOError:
            print 'IOError case: ', case_path
            return np.zeros((1, 1))
        except ValueError:
            print 'ValueError case: ', case_path
            return np.zeros((1, 1))

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

        self.t_gap_s = self.t_gap_ms / 1000
        n_points = self.data[:, 0, 0].size
        self.max_time = self.t_gap_s * (n_points)
        self.time = np.linspace(self.t_gap_s, self.max_time, n_points)
        # first normalise
        if ECG_type == -1:
            col_1 = lead_1 - 1
            col_2 = lead_2 - 1
            self.ECG_data = flipper * (self.data[:, col_1, :] - \
                                       self.data[:, col_2, :])
            self.y_label = r'$ \Delta V [mV]$'
        else:
            self.ECG_data = flipper * self.data[:, ECG_type, :]
            self.y_label = r'$ V [mV]$'

    def full_ECG_set(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        for case in range(2):
            plt.plot(self.time, self.ECG_data[:, case],
                    label=self.case_legends[case])
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

    def beat_ECG_Set(self, beat, save_file):
        n_beats = self.time.max() / (self.period_ms / 1000)
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
        for case in range(2):
            plt.plot(self.time[index_start:index_end],
                    self.ECG_data[index_start:index_end, case],
                    label=self.case_legends[case])
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

    def full_set_as_subplots(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        f.subplots_adjust(wspace=0.3)
        self.template.apply_figuresize_settings(f)
        n_plots = 2
        for case in range(n_plots):
            plt.subplot(n_plots, 1, (case + 1))
            plt.plot(self.time, self.ECG_data[:, case],
                     label=self.case_legends[case])
            plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                       rotation='horizontal')
            plt.legend(prop=self.lfont)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        self.template.apply_figuresize_settings(f)
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')

    def beat_set_as_subplots(self, beat, save_file):
        n_beats = self.time.max() / (self.period_ms / 1000)
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
        n_plots = len(self.case_diretories)
        for case in range(n_plots):
            axes = f.add_subplot(n_plots, 1, (case + 1))
            axes.plot(self.time[index_start:index_end],
                      self.ECG_data[index_start:index_end, case],
                      label=self.case_legends[case])
            plt.legend(prop=self.lfont)
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
