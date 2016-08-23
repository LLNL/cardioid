'''
Created on 27/02/2013

@author: butler
'''
import json


class JParser(object):
    '''
    Parser for json config file. As this is designed to be '
    '''

    def __init__(self, file_name):
        '''
        '''
        self.file_name = file_name
        """
        self.delta_t_ms
        self.sample_rate
        self.t_min_s
        self.t_max_s
        """
        self.parse()

    def parse(self):
        try:
            self.fd = open(self.file_name, 'r')
            jdict = json.load(self.fd)
            jdict = jdict["PSD"]
            self.delta_t_ms = float(jdict['delta_t_ms'])
            self.sample_rate = int(jdict['sample_rate'])
            self.t_min_s = jdict['t_min_s']
            self.t_max_s = jdict['t_max_s']
        finally:
            self.fd.close()
