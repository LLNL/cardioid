'''
Created on 09/05/2013

@author: butler
'''
import mpi4py.MPI as mpi
import controller

if __name__ == '__main__':
    print "Hello world"
    comm = mpi.COMM_WORLD
    comm.Barrier()
    rank = comm.Get_rank()
    if rank == 0:
        print "MPI initialized"
    A = controller.Controller(comm)
    if rank == 0:
        print "Controller instanciated"
    comm.Barrier()
    A.setup()
    if rank == 0:
        print "Setup Done"
    comm.Barrier()
    A.load_data()
    comm.Barrier()
    if rank == 0:
        print "data loaded"
    A.transpose_v2()
    comm.Barrier()
    if rank == 0:
        print "Data transposed"
    A.do_DFTs()
    if rank == 0:
        print "DFTs calculated"
    comm.Barrier()
    A.estimate_global_PSD('out_PSD.csv')
    comm.Barrier()
    if rank == 0:
        print "data gathered"
    A.save_modal_data([0, 35, 40])
    #A.save_modal_data([1, 2])
    comm.Barrier()
    mpi.Finalize()
