#!/bin/env python
'''
Created on 07/03/2013

@author: butler
'''
import os
import sys
import glob

if __name__ == "__main__":
    assert(len(sys.argv) == 2)
    mode = int(sys.argv[1])
    sshot_glob = './snapshot.0*'
    state_glob = 'state*'
    profile_glob = 'profile*'
    sshot_list = glob.glob(sshot_glob)
    sshot_list.sort()
    for sshot_dir in sshot_list:
        state_files = glob.glob(())
        profile_files = glob.glob(())
    for level_1 in list_1:
        list_2 = glob.glob((level_1 + os.sep + tier_2))
        for level_2 in list_2:
            list_3 = glob.glob((level_2 + os.sep + tier_3))
            for level_3 in list_3:
                list_4 = glob.glob((level_3 + os.sep + tier_4))
                for level_4 in list_4:
                    list_5 = glob.glob((level_4 + os.sep + tier_5))
                    list_5.sort()
                    final_sshot_dir = list_5[-1]
                    state_files = glob.glob((final_sshot_dir + os.sep +
                                             state_glob))
                    profile_files = glob.glob((final_sshot_dir + os.sep +
                                               profile_glob))
                    if execute_flag == 1:
                        for state_file in state_files:
                            #os.remove(state_file)
                            pass
                        for profile_file in profile_files:
                            #os.remove(profile_file)
                            pass
                    else:
                        if len(state_files) > 1:
                            print final_sshot_dir
