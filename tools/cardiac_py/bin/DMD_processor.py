#!/bin/env python
'''
Created on 27/02/2013

@author: butler
'''
import DMD
import sys

if __name__ == '__main__':
    assert(len(sys.argv) == 2)
    directory = sys.argv[1]
    controller = DMD.controller.Controller(directory)
    controller.setup_problem()
    controller.solve_for_DMD()
