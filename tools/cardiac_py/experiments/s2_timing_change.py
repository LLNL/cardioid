'''
Created on 04/02/2013

@author: butler
'''
import control


class Experiment():
    '''
    classdocs
    '''
    def __init__(self):
        self.bgq_sim_dir = '/hsm/VR0236/shared/cardioid_sims/'
        self.local_sim_dir = '/Users/butler/Desktop/CARDIOID/'
        self.ecg_relative_dir = 'ECG_data'
        self.experiment_no = 15
        self.gkr_normal_bin_name = 'cardioid_avoca_spi_normal'
        self.gkr_blocked_bin_name = 'cardioid_avoca_spi_zeroed'
        self.mesh_data_dir = '/vlsci/VR0236/butler/geometries'
        self.pacing_frequency_ms = 2000
        self.sim_duration_beats = 15
        self.s2_beat = 10
        self.s2_min_offset_ms = 450
        self.s2_max_offset_ms = 750
        self.s2_interval_ms = 5

    def build_object_writer(self):
        self.o_writer = control.object_creator.Creator()
        dt = float(self.o_writer.dt_ms.split(' ')[0])
        self.o_writer.maxLoop = int(self.sim_duration_beats *
                                    self.pacing_frequency_ms / dt)
        #Checkpoint
        self.o_writer.checkpointRate = int(self.s2_beat *
                                           self.pacing_frequency_ms / dt)

    def construct_object_file_with_s2(self, case_directory, s2_timing):
        pass

    def initialize_bgq(self):
        pass

    def pull_ecg_data(self):
        pass
