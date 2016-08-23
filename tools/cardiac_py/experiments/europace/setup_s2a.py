'''
Created on 05/07/2013

@author: butler
'''

import glob
import os
import sys
import shutil


def __add_sep(input_path):
    '''
    Apply this function to a path to add a os.sep if no os.sep is trailing
    '''
    if input_path[-1] == os.sep:
        return input_path
    else:
        return (input_path + os.sep)


def setup_s2_runs(BCL_dir, s2_zero_time, s2_gap_ms, n_s2):
    ''' For all IKr blocks of a given BCL create s2 stimulus runs and submit them
    '''
    object_proto = __add_sep(BCL_dir) + "object.data.s2proto_a"
    ikr_glob = __add_sep(BCL_dir) + "IKr_block_100a*"
    ikr_dirs = glob.glob(ikr_glob)
    for ikr_path in ikr_dirs:
        priming_path = __add_sep(ikr_path) + "priming"
        print 'Priming path', priming_path
        sshot_glob = __add_sep(priming_path) + 'snapshot.0*'
        p_sshot_dirs = glob.glob(sshot_glob)
        p_sshot_dirs.sort()
        p_final_shot_dir = p_sshot_dirs[-1]
        print 'p_final_shot_dir', p_final_shot_dir
        p_final_shot_name = p_final_shot_dir.split(os.sep)[-1]
        for s2_i in range(n_s2):
            s2_offset = s2_gap_ms * s2_i
            s2_dir_name = "s2_offset_{:03}".format(s2_offset)
            s2_path = __add_sep(ikr_path) + s2_dir_name
            # 1 copy and link the files we require in
            shutil.copytree(priming_path, s2_path, symlinks=True,
                            ignore=shutil.ignore_patterns("snapshot.0*",
                                                          "cardioid.out",
                                                          "ecg.out",
                                                          "electrode#*"))
            # now we need to cp and mod the prototype object.data
            s2_object = __add_sep(s2_path) + "object.data"
            s2_time = s2_offset + s2_zero_time
            transform_object_data(s2_time, object_proto, s2_object)
            # now we need to do a symlink of the final timestep's snapshot
            # directory to the correct location
            new_sshot_loc = __add_sep(s2_path) + p_final_shot_name
            os.symlink(p_final_shot_dir, new_sshot_loc)


def transform_object_data(s2_val, prototype_location, new_location):
    new_OD = open(new_location, 'w')
    for line in open(prototype_location):
        if line.find("99999") > 0:
            split_line = line.split("99999")
            assert (len(split_line) == 2)
            line = split_line[0] + repr(s2_val) + split_line[1]
        new_OD.write(line)

if __name__ == '__main__':
    assert(len(sys.argv) == 5)
    setup_s2_runs(sys.argv[1],int(sys.argv[2]),int(sys.argv[3]), int(sys.argv[4]))