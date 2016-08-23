'''
Created on 08/01/2013

@author: butler
'''
import numpy as np
import plot_schema
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager
import os
import glob

''' MATLPLOTLib preamble code'''


class GECGs():
    '''
    This class is designed to be more generic than the previous ECG
    incarnation. Rather than attemping to understand the nature of the data
    plots, this class gets passed the information into more generic interfaces.
    This still presumes that data generated for the ECG runs from the start of
    the simulation.
    '''
    def __init__(self, directory, full=True, parallel=False, wedge_e_ver=1):
        '''
        Constructor
        '''
        self.parallel = parallel
        self.t_gap_ms = 5.0
        self.directory = directory
        self.case_directories = []
        self.case_legends = []
        self.full = full
        if self.full:
            self.electrodes = ['electrode#000448302', 'electrode#000451300',
                               'electrode#000452730', 'electrode#000453393',
                               'electrode#000457525', 'electrode#000458894',
                               "electrode#000438028", "electrode#000460291"]
        else:
            '''
            wedge_e_ver
            There has been an unfortunate need to introduce different wedge
            electrodes.
            The first generation (wedge_e_ver = 1) will hold up well.
            Versions 2 and 3 are part of a set:
            Version 2 is considered a direct replacement for version 1 the
            location may be slightly different.
            Version 3 is in the same alignment within a intramural plane,
            however, the electrodes are close to the wedge than version 1 or 2
            which are at the outer edge of the water box.
            '''
            if wedge_e_ver == 1:
                self.electrodes = ['electrode#000094150',
                                   'electrode#000092294']
            elif wedge_e_ver == 2:
                self.electrodes = ['electrode#000094790',
                                   'electrode#000105838']
            elif wedge_e_ver == 3:
                self.electrodes = ['electrode#000093815',
                                   'electrode#000095757']
        self.period_ms = 1000
        self.template = plot_schema.PlotSchema()
        self.template.apply_fontsettings(plt)
        self.lfont = matplotlib.font_manager.FontProperties(size=(0.75 * self.template.fontsize))

    def load_data(self):
        self.stat_flag = np.zeros((len(self.case_directories)))
        n_cases = len(self.case_directories)
        print 'no of cases: ', repr(n_cases)
        data_set_index = 0
        for ii in range(n_cases):
            case_path = self.__add_sep(self.directory) + self.case_directories[ii]
            case_data = self.load_case_data(case_path)
            if ii == data_set_index:
                n_steps, n_electrodes = case_data.shape
                if n_steps == 1:
                    self.stat_flag[ii] = 0
                    data_set_index = data_set_index + 1
                    print 'this is bad'
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
                elif n_steps_c > n_steps:
                    temp = self.data
                    self.data = np.zeros((n_steps_c, n_electrodes, n_cases))
                    self.data[:n_steps, :, :] = temp[:, :, :]
                    self.stat_flag[ii] = 1
                    self.data[:, :, ii] = case_data[:, :]
                else:
                    self.data[:n_steps_c, :, ii] = case_data[:, :]
                    self.stat_flag[ii] = 1
        print "data loaded"

    def load_case_data(self, case_path):
        if self.parallel:
            return self.load_case_data_parallel(case_path)
        else:
            return self.load_case_data_serial(case_path)

    def load_case_data_parallel(self, case_path):
        """ Work with parallel data sets """
        n_electrodes = len(self.electrodes)
        try:
            for ii in range(n_electrodes):
                #
                electrode_path_glob = case_path + self.electrodes[ii] + \
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
                temp_data = temp_data[temp_data[:, 0].argsort()]
                temp_data = temp_data[:, 1]
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

    def load_case_data_serial(self, case_path):
        n_electrodes = len(self.electrodes)
        try:
            for ii in range(n_electrodes):
                electrode_path = self.__add_sep(case_path) + self.electrodes[ii]
                print 'Electrode path: ', electrode_path
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
    # desired interfaces

    def all_leads_subplot(self, case, save_file):
        frac = 1
        index_start = 0
        index_end = int(self.time.size * frac)
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        n_plots = 7
        for ii in range(n_plots):
            axes = f.add_subplot(n_plots, 1, (ii + 1))
            if ii == 0:
                flipper = -1
            else:
                flipper = 1
            self.set_ECG_type(ii, flipper)
            axes.plot(self.time[index_start:index_end],
                      self.ECG_data[index_start:index_end, case])
            plt.legend(prop=self.lfont)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.__add_sep(self.directory) + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')

    def full_ECG_case(self, case, save_file):
        if self.stat_flag[case] == 1:
            self.template.apply_fontsettings(plt)
            f = plt.figure()
            self.template.apply_figuresize_settings(f)
            axes = plt.axes()
            plt.plot(self.time, self.ECG_data[:, case])
            plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
            plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                       rotation='horizontal')
            self.template.apply_figuresize_settings(f)
            for x_ticl_i in axes.get_xticklabels():
                x_ticl_i.set_fontsize(self.template.get_fontsize())
            for y_ticl_i in axes.get_yticklabels():
                y_ticl_i.set_fontsize(self.template.get_fontsize())
            save_loc = self.__add_sep(self.directory) + save_file
            plt.savefig(save_loc, dpi=350, bbox_inches='tight')

    def full_ECG_set(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        for case in range(len(self.case_directories)):
            if self.stat_flag[case] == 1:
                plt.plot(self.time, self.ECG_data[:, case],
                         label=self.case_legends[case])
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        plt.legend(prop=self.lfont, bbox_to_anchor=(1.40, 0.40), loc=5)
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.__add_sep(self.directory) + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')

    def beat_ECG_case(self, beat, case, save_file):
        if self.stat_flag[case] == 1:
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
            plt.plot(self.time[index_start:index_end],
                     self.ECG_data[index_start:index_end, case])
            plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
            plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                       rotation='horizontal')
            plt.legend(props=self.lfont)
            self.template.apply_figuresize_settings(f)
            for x_ticl_i in axes.get_xticklabels():
                x_ticl_i.set_fontsize(self.template.get_fontsize())
            for y_ticl_i in axes.get_yticklabels():
                y_ticl_i.set_fontsize(self.template.get_fontsize())
            save_loc = self.__add_sep(self.directory) + save_file
            plt.savefig(save_loc, dpi=350, bbox_inches='tight')

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
        for case in range(len(self.case_directories)):
            if self.stat_flag[case] == 1:
                plt.plot(self.time[index_start:index_end],
                         self.ECG_data[index_start:index_end, case],
                         label=self.case_legends[case])
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                   rotation='horizontal')
        plt.legend(props=self.lfont,
                   bbox_to_anchor=(1.40, 0.40), loc=5)
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.__add_sep(self.directory) + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')

    def arbitary_t_ECG_case(self, min_t_s, max_t_s, case, save_file):
        print 'Not yet implemented'
        pass

    def full_set_as_subplots(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        f.subplots_adjust(wspace=0.3)
        self.template.apply_figuresize_settings(f)
        n_plots = len(self.case_directories)
        for case in range(n_plots):
            plt.subplot(n_plots, 1, (case + 1))
            plt.plot(self.time, self.ECG_data[:, case],
                     label=self.case_legends[case])
            plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                       rotation='horizontal')
            plt.legend(prop=self.lfont)
        plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
        self.template.apply_figuresize_settings(f)
        save_loc = self.__add_sep(self.directory) + save_file
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
        save_loc = self.__add_sep(self.directory) + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')

    def frac_ECG_subplots(self, frac, save_file):
        index_start = 0
        index_end = int(self.time.size * frac)
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
        save_loc = self.__add_sep(self.directory) + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')

    def frac_ECG_set(self, frac, save_file):
        index_start = 0
        index_end = int(self.time.size * frac)
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.grid()
        for case in range(len(self.case_directories)):
            if self.stat_flag[case] == 1:
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
        save_loc = self.__add_sep(self.directory) + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')

    def frac_ECG_case(self, frac, case, save_file):
        if self.stat_flag[case] == 1:
            index_start = 0
            index_end = int(self.time.size * frac)
            self.template.apply_fontsettings(plt)
            f = plt.figure()
            self.template.apply_figuresize_settings(f)
            axes = plt.axes()
            plt.grid()
            plt.plot(self.time[index_start:index_end],
                     self.ECG_data[index_start:index_end, case])
            plt.xlabel('$t (s)$', fontsize=self.template.get_fontsize())
            plt.ylabel(self.y_label, fontsize=self.template.get_fontsize(),
                       rotation='horizontal')
            plt.legend(prop=self.lfont)
            self.template.apply_figuresize_settings(f)
            for x_ticl_i in axes.get_xticklabels():
                x_ticl_i.set_fontsize(self.template.get_fontsize())
            for y_ticl_i in axes.get_yticklabels():
                y_ticl_i.set_fontsize(self.template.get_fontsize())
            save_loc = self.__add_sep(self.directory) + save_file
        plt.savefig(save_loc, dpi=350, bbox_inches='tight')

    def write_csv_case(self, case, save_file):
        a = np.zeros((self.time.size, 2))
        a[:, 0] = self.time[:]
        a[:, 1] = self.ECG_data[:, case]
        save_file_path = self.__add_sep(self.directory) + save_file
        np.savetxt(save_file_path, a, delimiter=",")

    def write_csv_all_leads(self, case, save_file):
        # Force 'full' ECG parameters
        assert(len(self.electrodes) == 8)
        # get time from whatever data we have
        csv_dat = np.zeros((self.time.size, 8))
        csv_dat[:, 0] = self.time[:]
        for ii in range(7):
            if ii == 0:
                flipper = -1
            else:
                flipper = 1
            self.set_ECG_type(ii, flipper)
            csv_dat[:, (ii + 1)] = self.ECG_data[:, case]
        save_file_path = self.__add_sep(self.directory) + save_file
        np.savetxt(save_file_path, csv_dat, delimiter=",")

    def __add_sep(self, input_path):
        '''
        Apply this function to a path to add a os.sep if no os.sep is trailing
        '''
        if input_path[-1] == os.sep:
            return input_path
        else:
            return (input_path + os.sep)
