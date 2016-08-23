#!/bin/env python
'''
Created on 19/02/2013

@author: butler
'''
import os
import sys
import fileinput
import shutil

if __name__ == '__main__':
    args = sys.argv
    assert(len(args) == 4)
    s2_lower_ms = int(args[1])
    delta_ms = int(args[2])
    n_s2 = int(args[3])
    s2temp_dir = 's2temp'
    ikr_dirs = ["ikr_block_0", "ikr_block_33", "ikr_block_66", "ikr_block_100"]
    for ikr_dir in ikr_dirs:
        for ii in range(1, (n_s2 + 1)):
            s2run_dir = 's2_run_' + repr(ii)
            s2run_path = ikr_dir + os.sep + s2run_dir
            s2temp_path = ikr_dir + os.sep + s2temp_dir
            shutil.copytree(s2temp_path, s2run_path, symlinks=True)
            object_data = s2run_path + os.sep + "object.data"
            s2_time = s2_lower_ms + delta_ms * (ii - 1)
            for line in fileinput.input(object_data, inplace=1):
                print line.rstrip('\n').replace('20510', repr(s2_time))
