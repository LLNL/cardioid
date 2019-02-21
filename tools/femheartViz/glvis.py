#!/usr/bin/env python3

import os
import sys
import argparse
import glob
import re

# To create a GIF movie with GLVIS do:
# 1. python3 glvis.py niederer.vtk sol.*.gf
# 2. glvis -run GIFScript.glvis
# 3. convert -delay 20 image.????.png Movie.gif

def get_numbers_from_filename(filename):
    return re.search(r'\d+', filename).group(0)

class glvis:
    '''
    glivis solution reader
    '''
    def __init__(self, file_stem = 'sol.*.gf', mesh_file = 'mesh.vtk'):
        # file stem (we can pass something like *.0000)
        g_matcher = file_stem
        self.file_list = glob.glob(g_matcher)
        self.file_list.sort()

        self.mesh_file = mesh_file

        # Check if we found files
        if len(self.file_list) == 0:
            print('No files matching {}'.format(file_stem))
            self.foundFiles = False
        else:
            print('Files matching {}'.format(self.file_list))
            self.foundFiles = True

    def getFoundFiles(self):
        '''Check if files are detected or not'''
        return self.foundFiles

    def create_GIFScript(self):
        '''Create the script for animation in glvis'''

        print('Creating GIFScript.glvis...'.format(self.file_list), end='')
        out_file = open('GIFScript.glvis','w')

        textBlock = [
        '# Visualization window geometry',
        'window 0 0 300 300',
        '# Initial solution',
        'solution' + ' ' + self.mesh_file + ' ' + self.file_list[0],
        ' ',
        '# Setup the GLVis scene. Executed after pressing the space bar.',
        '{',
        ' perspective off',
        '# view 0 0',
        '# viewcenter 0 0',
        ' zoom 1.95',
        '#   keys fAmeIIiiiiiiiiiiibbvuuuuuuuuuuu',
        ' valuerange -83.0 30',
        '}',
        ' '
        ]
        for line in textBlock:
            out_file.write(line)
            out_file.write('\n')
        out_file.write('# Take multiple screenshots. Executed after pressing the space bar. \n')
        out_file.write('{ \n')
        for file in self.file_list:
            itime = get_numbers_from_filename(file)
            out_file.write('  solution' + ' ' + self.mesh_file + ' ' + file + ' ' + \
                           'valuerange -83.0 30' + ' ' + \
                           'screenshot' + ' ' + 'image.'+itime+'.png\n'
                           )
        out_file.write('} \n')

        out_file.close()
        print('done')

if __name__ == "__main__":

    nargs = len(sys.argv)

    if nargs == 1:
        mesh_file = 'niederer.vtk'
        sol_stem = 'sol.*.gf'
    elif nargs == 2:
        mesh_file = sys.argv[1]
        sol_stem = 'sol.*.gf'
    elif nargs == 3:
        mesh_file = sys.argv[1]
        sol_stem = sys.argv[2]
    else:
        raise Exception('usage: python3 ~/tools_py/{} snapshot_dir_stem anatomy_dir'.format(os.path.basename(sys.argv[0])))

    # Construct of the glvis interface
    glvis_sol = glvis(sol_stem, mesh_file)

    # Construct the script GIFScript.glvis to construct the animation
    glvis_sol.create_GIFScript()
