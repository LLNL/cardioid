'''
Created on 10/07/2013

@author: butler
'''

import os
import glob


def remove_checkpoints(directory, mode):
    '''
    Provide a cardioid case directory and a mode which interprets what to do:
    mode == 0:
        Print what mode 1 and mode 2 will do
    mode == 1:
        Delete all state files except if they occur in the final directory.
        Will ignore sshot directories which are sym-links
    mode == 2:
        Delete all state files contained in the cardioid case directory.
        Will ignore sshot directories which are sym-links
    mode == 3:
        Delete all state files in a cardioid directory with the following exceptions
        1) Will ignore directories which are sim-links
        2) Will ignore a 
    '''
    sshot_glob = directory + os.sep + 'snapshot.0*'
    state_glob = 'state#*'
    profile_glob = 'profile#*'
    sshot_list = glob.glob(sshot_glob)
    sshot_list.sort()
    for sshot_dir in sshot_list:
        state_files = glob.glob((sshot_dir + os.sep + state_glob))
        profile_files = glob.glob((sshot_dir + os.sep + profile_glob))
        if len(profile_files) > 0:
            if mode > 0:
                for profile_file in profile_files:
                    os.remove(profile_file)
        if len(state_files) > 0:
            if os.path.islink(sshot_dir):
                print "This directory is sym-linked from somewhere else", \
                sshot_dir
                print "It will be skipped"
                continue
            if mode == 0:
                if sshot_dir == sshot_list[-1]:
                    print 'State files in this directory would be kept with'
                    print 'mode one', sshot_dir
                else:
                    print 'Mode > 0 will delete state files in this'
                    print 'directory', sshot_dir
            elif mode == 1:
                if sshot_dir != sshot_list[-1]:
                    for state_file in state_files:
                        os.remove(state_file)
            elif mode == 2:
                for state_file in state_files:
                    os.remove(state_file)
            elif mode == 3:
                if 'priming' in directory:
                    #behave like mode one:
                    pass
                else:
                    #behave like mode 2:
                    pass
