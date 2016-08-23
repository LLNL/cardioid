'''
Created on 10/06/2013

@author: butler
'''

import sys
import controller

if __name__ == '__main__':
    if (len(sys.argv) == 2):
        directory = sys.argv[1]
    elif (len(sys.argv) == 1):
        directory = './'
    else:
        sys.exit(1)
    controller = controller.Controller(directory)
    # mode agnostic
    controller.setup_problem()
    if controller.config.solver == 'DMD':
        controller.solve_for_DMD()
    elif controller.config.solver == 'POD':
        controller.solve_for_POD()
        controller.compute_and_write_modes([0, 1, 2, 3])
    else:
        print "Unknown solver type: ", controller.config.solver
        print "aborting"
        sys.exit(1)
    if controller.config.calc_modes:
        controller.compute_and_write_modes(controller.config.modes)
    sys.exit(0)
