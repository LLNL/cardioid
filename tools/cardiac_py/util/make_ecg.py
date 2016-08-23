from __future__ import division
'''
Created on Aug 29, 2012

@author: butler
'''
''' This is now obsolete. This is not required.
'''
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import glob
import sys
import shlex


def read_files(electrode_1, electrode_2):
    """ Read all of the ecg files and return two numpy vectors """
    file_list = glob.glob("./ecg*")
    file_list.sort()
    n_files = len(file_list)
    potential_1 = np.zeros(n_files)
    potential_2 = np.zeros(n_files)
    ii = 0
    for file_name in file_list:
        print file_name
        fd = open(file_name, 'r')
        jj = 0
        for line in fd:
            index, val = shlex.split(line)
            if int(index) == electrode_1:
                potential_1[ii] = float(val)
                jj = jj + 1
            if int(index) == electrode_2:
                potential_2[ii] = float(val)
                jj = jj + 1
            if jj == 2:
                break
        ii = ii + 1
        fd.close()
    return potential_1, potential_2


def save(file_name_prefix, e_data_1, e_data_2):
    fname_1 = file_name_prefix + "_1"
    fname_2 = file_name_prefix + "_2"
    np.save(fname_1, e_data_1)
    np.save(fname_2, e_data_2)
    return


def plot(file_name, x, y):
    plt.plot(x, y)
    plt.xlabel('time (s)')
    plt.ylabel('Delta v')
    plt.savefig(file_name, dpi=400, bbox_inches='tight')


if __name__ == "__main__":
    args = sys.argv
    if not len(args) == 6:
        print "incorrrect arguments"
        print "correct arguments: First_electrode_location" +  \
        "second_electrode_Location"
        print "save_file_prefix plot_file_name time_spacing (ms)"
    potential_1, potential_2 = read_files(int(args[1]), int(args[2]))
    save_prefix = args[3]
    plot_fname = args[4]
    time_spacing = float(args[5]) / 1000.0
    save('test_fname', potential_1, potential_2)
    pot_diff = potential_2 - potential_1
    times = np.linspace(0, (pot_diff.size * time_spacing), pot_diff.size)
    plot(plot_fname, times, pot_diff)
    print "done"
