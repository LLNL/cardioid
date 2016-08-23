'''
Created on 17/02/2013

@author: butler
'''
import numpy as np
import os
from . import akima


class Analyse():
    '''
    Tool for analysing the lead 1 wave structure. Primary goal is for finding
    feature points on the wave to be used as offsets. Works on an individual
    cardioid simulation basis. Has graphical and programattic output.
    '''

    def __init__(self, directory, debug=False):
        '''
        Constructor
        '''
        self.directory = directory
        self.electrodes = ["electrode#000438028", "electrode#000460291"]
        self.period_ms = 1000
        self.t_gap_ms = 5.0
        self.debug = True

    def load_data(self):
        n_electrodes = len(self.electrodes)
        for ii in range(n_electrodes):
            electrode_path = self.directory + os.sep + self.electrodes[ii]
            temp_data = np.loadtxt(electrode_path)
            if ii == 0:
                n_times = temp_data.size
                self.electrode_data = np.zeros((n_times, n_electrodes))
            self.electrode_data[:, ii] = temp_data[:]

    def construct_lead_data(self, flipper=1):
        '''
        Construct ECG and put in correct timeframe. Note that is presumes the
        ECG sensor is running from the first interval into the simulation
        '''
        self.t_gap_s = self.t_gap_ms / 1000
        n_points = self.electrode_data[:, 0].size
        self.max_time = self.t_gap_s * (n_points)
        self.time = np.linspace(self.t_gap_s, self.max_time, n_points)
        self.L1 = flipper * (self.electrode_data[:, 1] -
                             self.electrode_data[:, 0])

    def pick_beat(self, beat):
        """
        Cannot be first beat as we are missing the data we need
        """
        n_beats = self.time.max() / (self.period_ms / 1000)
        assert(beat < n_beats)
        index_start = (self.period_ms / self.t_gap_ms) * (beat - 1)
        index_end = (self.period_ms / self.t_gap_ms) * (beat)
        if beat > 1:
            index_start = index_start - 1
            index_end = index_end - 1
        self.beat_t = self.time[index_start:index_end]
        self.beat_L1 = self.L1[index_start:index_end]
        # Do some simple plot checks here

    def calculate_extema_and_estimate_T_wave_max(self, frac_of_peak_filter=0.05
                                                 , QRS_exclude=0.1):
        """
        """
        frac_of_peak_filter = 0.05
        self.spline_coeffs = akima.create(self.beat_t, self.beat_L1)
        self.extrema = akima.extremes(self.beat_t, self.spline_coeffs)
        self.voltage_maxima = self.extrema[0]
        self.voltage_maxima_t = self.extrema[1]
        # first let us do some checking
        print 'beat start: ', repr(self.beat_t[0])
        print 'beat end: ', repr(self.beat_t[-1])
        # assume that maxima we are looking for are +ve and at least 5% of peak
        v_maximum = self.voltage_maxima.max()
        vfilter = v_maximum * frac_of_peak_filter
        indices = (self.voltage_maxima > vfilter)
        self.voltage_maxima = self.voltage_maxima[indices]
        self.voltage_maxima_t = self.voltage_maxima_t[indices]
        self.offset_t = self.voltage_maxima_t - self.beat_t[0]
        # Exclude A fraction of the beat assuming it is during the QRS complex.
        # By doing this we should be able to say t-peak is max(maxima)
        if self.debug:
            print 'Offset times', self.offset_t
        indices = (self.offset_t > QRS_exclude)
        self.voltage_maxima = self.voltage_maxima[indices]
        self.offset_t = self.offset_t[indices]
        if self.debug:
            print 'Offset times', self.offset_t
            print 'Voltages', self.voltage_maxima
            print 'Delta t for maxima ', self.offset_t
            print 'Corresponding voltage ', self.voltage_maxima
        # Assume t-wave peak is the largest
        indices = self.voltage_maxima.argsort()
        self.offset_t = self.offset_t[indices]
        self.voltage_maxima = self.voltage_maxima[indices]
        return self.offset_t[-1], self.voltage_maxima[-1]

