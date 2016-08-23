'''
Created on 03/03/2013

@author: butler
'''
import JParser
import os
import modred
import cPickle
import vec_containers


class Controller():
    '''
    DMD controller.
    This class forms the primary interface to the DMD calculation, wrapping
    modred and the io classes into convenient buckets. Currently the unknown
    factor is determining which modes we need to write out.
    '''

    def __init__(self, directory):
        '''
        Constructor:
        Instanciate
        '''
        self.directory = directory
        self.jfile = 'DMD.json'
        self.sshot_stem = 'snapshot.'
        self.pio_stem = 'anatomy#'

    def setup_problem(self):
        """ Construct file list required for dmd analysis and initialise the
        modred classes
        """
        config_file = self.directory + os.sep + self.jfile
        self.config = JParser.JParser(config_file)
        t_min_steps = self.config.t_min_s / self.config.delta_t_ms * 1000
        t_max_steps = self.config.t_max_s / self.config.delta_t_ms * 1000
        sample_step_min = int(t_min_steps / self.config.sample_rate)
        sample_step_max = int(t_max_steps / self.config.sample_rate)
        self.sshot_list = []
        for ii in range(sample_step_min, (sample_step_max + 1)):
            sshot_dir = self.sshot_stem + "{0:012}".format((ii *
                self.config.sample_rate))
            sshot_path = self.directory + os.sep + sshot_dir
            # we now
            self.sshot_list.append(sshot_path)
        self.direct_snapshots = []
        for sshot_dir in self.sshot_list:
            self.direct_snapshots.append(vec_containers.CustomVecHandle(
                (sshot_dir + os.sep + self.pio_stem + "*")))

    def solve_for_DMD(self):
        '''
        Solve for the ritz_values, mode_norms and build_coefficients.
        '''
        self.dmd = modred.DMDHandles(vec_containers.inner_product)
        ritz_vals, mode_norms, build_coeffs = self.dmd.compute_decomp(
            self.direct_snapshots)
        fd_ritz = open((self.directory + 'ritz.cPickle'), 'w')
        cPickle.dump(ritz_vals, fd_ritz)
        fd_ritz.close()
        fd_mode = open((self.directory + 'mode.cPickle'), 'w')
        cPickle.dump(mode_norms, fd_mode)
        fd_mode.close()
        fd_coeffs = open((self.directory + 'coeffs.cPickle'), 'w')
        cPickle.dump(build_coeffs, fd_coeffs)
        fd_coeffs.close()
        print 'ritz_vals', ritz_vals
        print 'mode_norms', mode_norms

    def solve_for_POD(self):
        self.dmd = modred.pod.PODHandles(vec_containers.inner_product)
        eigen_vec_handles, eigen_vals = self.dmd.compute_decomp(
                                            self.direct_snapshots)
        print 'eigen_vals', eigen_vals

    def compute_and_write_modes(self, mode_list):
        '''
        '''
        # First get reference case to work with
        modes = []
        for ii in range(len(mode_list)):
            sshot_spec = self.sshot_list[ii] + os.sep + self.pio_stem + "*"
            print 'Mode being saved to: ', sshot_spec
            modes.append(vec_containers.CustomVecHandle(sshot_spec))
        self.dmd.compute_modes(mode_list, modes)
