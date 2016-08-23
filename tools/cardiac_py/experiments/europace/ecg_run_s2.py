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


def submit_ECG_as_req(BCL_dir, expect_sshots, check_electrode=False):
    ''' Loop over all of the s2_directories if n_snapshots not the expected
    number resubmit the EP job. If n_snapshots is as expected check whether a
    electrode file is the right length. If not submit ECG job.
    '''

    ikr_glob = __add_sep(BCL_dir) + "IKr_*"
    ikr_dirs = glob.glob(ikr_glob)
    for ikr_path in ikr_dirs:
        s2_glob = __add_sep(ikr_path) + "s2_offset_*"
        s2_dirs = glob.glob(s2_glob)
        for s2_path in s2_dirs:
            sshot_glob = __add_sep(s2_path) + 'snapshot.0*'
            s2_sshot_dirs = glob.glob(sshot_glob)
            n_sshots = len(s2_sshot_dirs)
            if n_sshots == expect_sshots:
                if check_electrode:
                    # Here we will optionally check whether
                    pass
                print 'submit ECG'
                os.chdir(s2_path)
                os.system("sbatch slurm_ECG_both")
            else:
                print 'CHECK EP'
                print 'for s2_dir: ', s2_path
                print 'n_sshots: ', n_sshots
                # Incorrect number of sshot dirs
                #os.system('sbatch slurm_EP')

if __name__ == '__main__':
    assert(len(sys.argv) == 3)
    submit_ECG_as_req(sys.argv[1], int(sys.argv[2]))
