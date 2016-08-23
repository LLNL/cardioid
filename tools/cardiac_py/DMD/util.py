'''
Created on 03/07/2013

@author: butler
'''
import numpy as np
import cPickle
import os


def dump_Ritz_vals_as_csv(directory):
    """
    Take a directory with the ritz coefficients and dump the data out in a form
    Usable with excel etc.
    """
    fd_ritz = open((directory + os.sep + 'ritz.cPickle'), 'r')
    ritz_p_arr = cPickle.load(fd_ritz)
    ritz_np = np.array(ritz_p_arr)
    n_modes = ritz_np.size
    separated = np.zeros((n_modes, 2))
    separated[:, 0] = ritz_np.real
    separated[:, 1] = ritz_np.imag
    np.savetxt("ritz_val.csv", separated, delimiter=',')


def dump_Mode_amplitude_as_csv(directory):
    fd_mode = open((directory + os.sep + 'mode.cPickle'), 'r')
    mode_p_arr = cPickle.load(fd_mode)
    mode_np = np.array(mode_p_arr)
    np.savetxt("mode_val.csv", mode_np, delimiter=",")
