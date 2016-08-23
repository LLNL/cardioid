'''
Created on 14/01/2013

@author: butler
'''
from . import iter_read
import os
import math

class SensorCheck():
    '''
    classdocs
    '''


    def __init__(self, old_dir):
        '''
        Constructor
        '''
        self.old_dir = old_dir
        self.debug = False
        self.sshot_dir_name = 'snapshot.initial'
        self.anat_match = 'anatomy#*'
        self.mapping_file = 'mapping.txt'
        self.sensor_file = 'sensor.txt'
        old_g_matcher = self.old_dir + os.sep + self.sshot_dir_name + os.sep + self.anat_match
        self.old_file = iter_read.Reader(old_g_matcher)
        self.o_num_records = long(self.old_file.header_vars['nrecord'])
        self.o_nx = int(self.old_file.header_vars["nx"])
        self.o_ny = int(self.old_file.header_vars["ny"])
        self.o_nz = int(self.old_file.header_vars["nz"])
        self.list_of_sensors = []
        
    def readsensor(self):
        old_sensor_fname = self.old_dir + os.sep + self.sensor_file
        old_sensor_file = open(old_sensor_fname,'r')

        for line in old_sensor_file:
            oGIDs = line.strip("\n").split()
            for oGID in oGIDs:
                
                try:
                    if len(oGID) > 1:
                        self.list_of_sensors.append(long(oGID))
                except ValueError:
                    print 'oGID, ', oGID
                    raise
        
    def CheckSensor(self):
        try:
            rec_count = 0
            self.sensor_found = 0
            
            for record in self.old_file:
                if len(record) < 60:
                    print "warning failed line"
                    continue
                else:
                    GID, data = record.strip(" ").split(" ",1)
                # Use old record count to split files
                if GID in self.list_of_sensors:
                    self.sensor_found = self.sensor_found + 1
                rec_count = rec_count + 1
        except ValueError:
            print self.old_file.f_list[self.old_file.file_index]
            print 'record: ', record
            print 'rec_count: ', rec_count
        print 'N sensors: ', len(self.list_of_sensors)
        print 'Sensors found: ', self.sensor_found

    