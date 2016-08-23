from __future__ import division
'''
Created on Aug 29, 2012

@author: butler
'''
''' This is now obsolete. This is not required.
'''
import mpi4py.MPI as mpi
import numpy as np
import glob
import sys
import shlex


def average_nonzero(comm, rank_vector):
    rank = comm.Get_rank()
    procs = comm.Get_size()
    # prep on root
    if rank == 0:
        temp_v = np.zeros(rank_vector.size)
        vector = np.zeros(rank_vector.size)
    else:
        vector = None
    for ii in range(procs):
        comm.Barrier()
        if ii == 0:
            # assign straight to vector
            if rank == 0:
                vector = rank_vector

        else:
            # Perform communication
            comm.Barrier()
            if rank == ii:
                comm.Ssend([rank_vector, mpi.DOUBLE], dest=0)
            elif rank == 0:
                comm.Recv([temp_v, mpi.DOUBLE], source=ii)
                for jj in xrange(rank_vector.size):
                    if temp_v[jj] != 0:
                        if vector[jj] == 0.0:
                            vector[jj] = temp_v[jj]
                        else:
                            vector[jj] = 0.5 * (vector[jj] + temp_v[jj])
            comm.Barrier()
    comm.Barrier()
    return vector


def read_files(comm, electrode_list):
    """ Read all of the ecg files and return two numpy vectors """
    comm.Barrier()
    procs = comm.Get_size()
    rank = comm.Get_rank()
    file_list = glob.glob("./ecg*")
    file_list.sort()
    n_files = len(file_list)
    n_electrodes = len(electrode_list)
    comm.Barrier()
    if n_files < procs:
        print "To few files for number of ranks"
        comm.Abort()
    default_files_read = np.ceil(n_files / procs)
    # Rounding errors here should not be a problem unless
    # default size is very small
    end_size = n_files - (default_files_read * (procs - 1))
    assert end_size >= 1
    if (rank == (procs - 1)):
        n_files_2_read = end_size
    else:
        n_files_2_read = default_files_read
    first_file = rank * default_files_read
    last_file = rank * default_files_read + n_files_2_read
    # Create individual processor pointers
    potentials = np.zeros((n_files, n_electrodes))
    comm.Barrier()

    for ii in range(first_file, last_file + 1):
        file_name = file_list[ii]
        print file_name
        fd = open(file_name, 'r')
        jj = 0
        for line in fd:
            index, val = shlex.split(line)
            try:
                electrode_list.index(int(index))
                potentials[electrode_list.index(index), ii] = float(val)
                jj = jj + 1
            except ValueError:
                pass
            if jj == n_electrodes:
                break
        fd.close()
    print "this rank has read its files: ", rank
    comm.Barrier()
    print 'Reading done'
    # we now have data in all ranks
    print 'pool first electrode'
    potentials.shape = potentials.size
    n_potentials = average_nonzero(comm, potentials)
    print 'done'
    comm.Barrier()
    if rank == 0:
        n_potentials.shape = (n_electrodes, n_files)
    return n_potentials

if __name__ == "__main__":
    comm = mpi.COMM_WORLD
    procs = comm.Get_size()
    rank = comm.Get_rank()

    args = sys.argv
    if not len(args) >= 3:
        print "incorrrect arguments"
        print "correct arguments: save_file"
        print "Followed by list of electrode locations"
    comm.Barrier()
    electrode_list = map(int, args[2:])
    potentials = read_files(comm, electrode_list)
    comm.Barrier()
    if rank == 0:
        print 'saving and plotting'
        save_file = args[1]
        (n_electrodes, n_times) = potentials.shape
        np.savetxt(save_file, n_electrodes)
    comm.Barrier()
    print 'All done'
