from __future__ import division
'''
Created on 03/03/2013

@author: butler
'''
import numpy as np
import jparser
import os
import pio
import numpy.fft as FFT
import mpi4py.MPI as mpi


class Controller():
    '''
    classdocs
    '''

    def __init__(self, comm):
        '''
        Constructor
        '''
        self.comm = comm
        self.directory = './'
        self.jfile = 'psd.json'
        self.sshot_stem = 'snapshot.'
        self.pio_glob = 'anatomy#*'

    def setup(self):
        # FIRST get my rank etc
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self.comm.Barrier()
        if self.rank == 0:
            print 'Hello world'
            self.setup_problem()
        else:
            self.sshot_list = None
            self.n_points_t = None
            self.n_points_x = None
        self.comm.Barrier()
        self.sshot_list = self.comm.bcast(self.sshot_list, root=0)
        self.n_points_t = self.comm.bcast(self.n_points_t, root=0)
        self.n_points_x = self.comm.bcast(self.n_points_x, root=0)

    def load_data(self):
        '''
        Load data in a distributed fashion.
        Each rank loads approximately the same size dataset. Where in this case
        the data loaded on each rank is all points in space and a sub-set of
        the data through time.

        Note that at extreme scales this will not work well
        as it is possible the rounding will screw up things for the last rank.
        This is an indication you are using too many processors for this code.
        '''
        self.default_n_files = int(np.round(len(self.sshot_list) / self.size))
        self.max_n_files = int(len(self.sshot_list) - (self.default_n_files *
                                                       (self.size - 1)))
        if self.rank < self.size - 1:
            self.my_n_files = self.default_n_files
        else:
            self.my_n_files = self.max_n_files
        self.my_time_steps = self.my_n_files
        print 'My n files', self.my_n_files
        # C-index
        self.my_data = np.zeros((self.n_points_x, self.my_n_files),
                                dtype=np.dtype('float64'))
        # loop over my data
        for ii in range(self.my_n_files):
            file_selector = self.rank * self.default_n_files + ii
            A = pio.iter_read.Reader(self.sshot_list[file_selector])
            jj = 0
            for record in A:
                val = float(record.rstrip("\n").strip('').split()[-1])
                self.my_data[jj, ii] = val
                jj = jj + 1
            print 'File: ', self.sshot_list[file_selector]
            print 'loaded'
        self.comm.Barrier()
        if self.rank == 0:
            print 'Loading Done'
        self.comm.Barrier()

    def transpose(self):
        '''
        Transpose the dataset such that instead of every point in space and
        a subset of points in time on a given MPI rank, every point in time and
        a partial set of points in space is on a given MPI rank. Required for
        the FFT to work.
        '''
        '''
        This should be re-written with scatter where data recieved is based on
        rank sent etc. Do this tonight.
        '''
        self.def_spat_points = int(np.round(self.n_points_x / self.size))
        self.max_spat_points = int(self.n_points_x - (self.def_spat_points *
                                                (self.size - 1)))
        if self.rank < self.size - 1:
            self.local_spat_points = self.def_spat_points
        else:
            self.local_spat_points = self.max_spat_points
        # create local array to insert data into
        self.rotated_data = np.zeros((self.local_spat_points, self.n_points_t),
                                     dtype=np.dtype('float64'))
        # my rank is
        for ii in range(self.n_points_x):
            if self.rank == 0:
                if (ii % 1000) == 0:
                    print "reordering point: ", ii
                    print "Of ", self.n_points_x
            # determine rank to send data too
            # Something here wrong
            recieve_rank = np.floor(ii / self.def_spat_points)
            if (recieve_rank == self.size):
                # correction for bad rounding at max. Occurs if n_points
                # on max is > normal_points
                recieve_rank = self.size - 1
            for send_rank in range(self.size):
                t_start = send_rank * self.default_n_files
                if send_rank < (self.size - 1):
                    send_size = self.default_n_files
                else:
                    send_size = self.max_n_files
                # this is probably not unique
                comm_id = ii + send_rank + recieve_rank
                if self.rank == recieve_rank:
                    recieve_buffer = np.zeros((1, send_size),
                                              dtype=np.dtype('float64'))
                    if send_rank != recieve_rank:
                        self.comm.Recv([recieve_buffer, mpi.DOUBLE],
                                    source=send_rank, tag=comm_id)
                if self.rank == send_rank:
                    send_buffer = self.my_data[ii, :]
                    if send_rank != recieve_rank:
                        self.comm.Send([send_buffer, mpi.DOUBLE],
                                    dest=recieve_rank, tag=comm_id)
                    else:
                        recieve_buffer = send_buffer
                self.comm.Barrier()
                # Should have data in buffer
                if self.rank == recieve_rank:
                    # where does the data go.
                    local_index = ii - (self.rank * self.def_spat_points)
                    self.rotated_data[local_index,
                                      t_start:(t_start + send_size)] = \
                                      recieve_buffer[:]
                self.comm.Barrier()
            self.comm.Barrier()
        self.comm.Barrier()
        del self.my_data

    def transpose_v2(self):
        '''
        Transpose the dataset such that instead of every point in space and
        a subset of points in time on a given MPI rank, every point in time and
        a partial set of points in space is on a given MPI rank. Required for
        the FFT to work.
        '''
        self.def_spat_points = int(np.round(self.n_points_x / self.size))
        self.max_spat_points = int(self.n_points_x - (self.def_spat_points *
                                                  (self.size - 1)))
        if self.rank < self.size - 1:
            self.local_spat_points = self.def_spat_points
        else:
            self.local_spat_points = self.max_spat_points
        # create local array to insert data into
        self.rotated_data = np.zeros((self.local_spat_points, self.n_points_t),
                                     dtype=np.dtype('float64'))
        # For each rank break my local datasets into
        send_buffer = []
        for ii in range(self.size):
            # slice out each ranks data
            min_spat_index = ii * self.def_spat_points
            if ((ii + 1) == self.size):
                max_spat_index = min_spat_index + self.max_spat_points
            else:
                max_spat_index = min_spat_index + self.def_spat_points
            sub_array = self.my_data[min_spat_index:max_spat_index, :]
            send_buffer.append(sub_array)
        for ii in range(self.size):
            self.comm.Barrier()
            local_data = self.comm.scatter(send_buffer, root=ii)
            min_t_index = self.default_n_files * ii
            if ((ii + 1) == self.size):
                max_t_index = min_t_index + self.max_n_files
            else:
                max_t_index = min_t_index + self.default_n_files
            self.rotated_data[:, min_t_index:max_t_index] = local_data[:, :]
        del self.my_data

    def do_DFTs(self):
        '''
        Calculate the DFT of the dataset provided through time. Store data such
        that the data may be sampled employed in multiple ways.
        '''
        self.comm.Barrier()
        # The DFT routine does some rounding which reduces the data to the
        # nearest power of two. Therefore we will do a dft first before storing
        # the data
        perc_complete = 0.0
        for ii in range(self.local_spat_points):
            # Progress block
            if self.rank == 0:
                perc_frac = float(ii) / self.local_spat_points * 100.0
                if ((perc_frac - 1) >= perc_complete):
                    print 'Rank 0 FFT completion: ', perc_frac
                    perc_complete = perc_frac
            dft_temp = FFT.rfft(self.rotated_data[ii, :])
            # create dataset based on first set
            if ii == 0:
                self.local_dft = np.zeros((self.local_spat_points,
                                        dft_temp.size), dtype=dft_temp.dtype)
            self.local_dft[ii, :] = dft_temp[:]
        self.comm.Barrier()

    def estimate_global_PSD(self, output_csv):
        '''
        Use the calculated DFTs to do a global PSD.
        Global PSD is calculated as the sum of pointwise PSD's
        '''
        self.comm.barrier()
        # Reduce data on each rank first
        (n_points, spectra_size) = self.local_dft.shape
        self.local_dft_energy = np.zeros((1, spectra_size),
                                         dtype=np.dtype('float64'))
        for ii in range(n_points):
            self.local_dft_energy[:] = self.local_dft_energy[:] + \
            (np.abs(self.local_dft[ii, :]) ** 2)
        # Reduce to rank zero before printing
        self.comm.Barrier()
        self.energy_estimate = np.zeros((1, self.local_dft_energy.size),
                                        dtype=np.dtype('float64'))
        self.comm.Reduce([self.local_dft_energy, mpi.DOUBLE],
                         [self.energy_estimate, mpi.DOUBLE], op=mpi.SUM,
                         root=0)
        # Shove data into an array to a csv can be written
        if self.rank == 0:
            timestep_s = (self.config.delta_t_ms * self.config.sample_rate /
                          1000.0)
            freq_t = FFT.fftfreq(((self.local_dft_energy.size - 1) * 2),
                                 d=timestep_s)
            freq = freq_t[0:self.local_dft_energy.size]
            write_array = np.zeros((freq.size, 2))
            write_array[:, 0] = freq[:]
            write_array[:, 1] = self.energy_estimate[:]
            np.savetxt(output_csv, write_array, delimiter=",")

    def setup_problem(self):
        """ Construct file list required for dmd analysis and initialise the
        modred classes
        """
        config_file = self.directory + os.sep + self.jfile
        self.config = jparser.JParser(config_file)
        t_min_steps = self.config.t_min_s / self.config.delta_t_ms * 1000
        t_max_steps = self.config.t_max_s / self.config.delta_t_ms * 1000
        sample_step_min = int(t_min_steps / self.config.sample_rate)
        sample_step_max = int(t_max_steps / self.config.sample_rate)
        self.sshot_list = []
        for ii in range(sample_step_min, (sample_step_max + 1)):
            # snapshot.000002394000 12 wide
            # snapshot.000000001000
            sshot_dir = self.sshot_stem + "{0:012}".format((ii *
                self.config.sample_rate))
            sshot_path = (self.directory + os.sep + sshot_dir + os.sep +
                          self.pio_glob)
            self.sshot_list.append(sshot_path)
        # Now we know how long the data is s
        self.n_points_t = len(self.sshot_list)
        # get number of points in space
        sample_pio = pio.seeker_read.Reader(self.sshot_list[0])
        self.n_points_x = sample_pio.records
        print 'If loading all data into memory'
        print 'Total number of points are: ', (self.n_points_x *
                                               self.n_points_t)

    def save_modal_data(self, mode_list):
        '''
        Save the dft information for a particular 'mode' (as in frequency).
        '''
        if self.rank == 0:
            timestep_s = (self.config.delta_t_ms * self.config.sample_rate /
                          1000.0)
            freq_t = FFT.fftfreq(((self.local_dft_energy.size - 1) * 2),
                                 d=timestep_s)
            self.freq = freq_t[0:self.local_dft_energy.size]
        mode_count = 0
        for mode in mode_list:
            if self.rank == 0:
                print "Dumping mode: ", mode
                print "At a frequency of [Hz]: ", self.freq[mode]
            self.comm.Barrier()
            # construct local mode vector
            local_mode = np.zeros((1, self.local_spat_points),
                                  dtype=self.local_dft.dtype)
            local_mode[:] = self.local_dft[:, mode]
            # create full size mode vector on all - just easier
            self.total_mode = np.zeros((1, self.n_points_x),
                                       dtype=local_mode.dtype)
            # reduce data
            for ii in range(self.size):
                start_index = ii * self.def_spat_points
                if (ii + 1) == self.size:
                    end_index = start_index + self.max_spat_points
                else:
                    end_index = start_index + self.def_spat_points
                if self.rank == 0:
                    if ii == 0:
                        print 'Mode dtype, ', self.total_mode.dtype
                        print 'start_index', start_index
                        print 'end_index', end_index
                        temp = self.total_mode[0, start_index:end_index]
                        print 'temp size', temp.size
                        print 'temp.shape',
                        print 'local_mode size', local_mode.size
                        self.total_mode[0, start_index:end_index] = \
                        local_mode[:]
                    else:
                        recv_buffer = np.zeros((1, (end_index - start_index)),
                                               dtype=local_mode.dtype)
                        recv_buffer = self.comm.recv(source=ii)
                        self.total_mode[0, start_index:end_index] = \
                        recv_buffer[:]
                elif self.rank == ii:
                    print 'Rank sending', ii
                    self.comm.send(local_mode, dest=0)
                self.comm.Barrier()
            self.comm.Barrier()
            # Data is ready to write
            if self.rank == 0:
                '''
                Calculated
                '''
                headers = ['Re_comp', 'Im_comp', 'Amplitude', 'phase']
                print 'Total mode size ', self.total_mode[0, :].size
                np.save('TOTAL_MODE', self.total_mode)
                data = np.zeros((self.total_mode[0, :].size, 4))
                for ii in range(self.total_mode[0, :].size):
                    data[ii, 0] = self.total_mode[0, ii].real
                    data[ii, 1] = self.total_mode[0, ii].imag
                    data[ii, 2] = np.sqrt((self.total_mode[0, ii].real ** 2) +
                                          (self.total_mode[0, ii].imag ** 2))
                    data[ii, 3] = np.arctan2(self.total_mode[0, ii].imag,
                                             self.total_mode[0, ii].real)
                mode_stem = 'mode'
                mode_stem_2 = 'mode_matrix'
                A = pio.append_complex.Piofy(self.sshot_list[mode_count],
                                         mode_stem, self.total_mode)
                A = pio.append_multiple.Piofy(self.sshot_list[mode_count],
                                           mode_stem_2, data, headers)
                A.Write()
            mode_count = mode_count + 1
        print 'Done writing modes'
        self.comm.Barrier()

    def save_trace_data(self):
        '''
        Save the time series of a given spatial point.
        Caveats:
        Must be completed after the rotation occurs.
        A snapshot file must be loaded to calculate metadata.
        '''
        pass
