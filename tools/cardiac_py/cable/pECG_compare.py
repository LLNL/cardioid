'''
Created on 18/06/2013

@author: butler
'''
import os
import pseudo_ECG
import matplotlib.pyplot as plt
import numpy as np

class Comparer():
    def __init__(self, directory="/Users/butler/Desktop/cable"):
        self.directory = directory
        self.WT_dir = 'cable_het_ex3'
        self.WT_legend = 'Control'
        self.drug_dir = 'cable_het_ex3_e4031'
        self.drug_legend = 'E-4031'

    def calc_data(self):
        wt_path = self.directory + os.sep + self.WT_dir
        wt_calc = pseudo_ECG.PseudoECG(wt_path,10)
        wt_calc.offset = 20
        wt_calc.get_data()
        wt_calc.calc()
        self.wt_data = wt_calc.epi_electrode[:] - wt_calc.endo_electrode[:]

        drug_path = self.directory + os.sep + self.drug_dir
        drug_calc = pseudo_ECG.PseudoECG(drug_path,10)
        drug_calc.offset = 20
        drug_calc.get_data()
        drug_calc.calc()
        self.drug_data = drug_calc.epi_electrode[:] - drug_calc.endo_electrode[:]

        self.time = np.zeros((self.wt_data.size))
        for ii in range(self.time.size):
            self.time[ii] = drug_calc.dt * ii

    def plotbeatECG(self):
        # Please not ebeing lazy.. hardcoding values
        time_beat = self.time[-1000:]
        wt_beat = self.wt_data[-1000:]
        drug_beat = self.drug_data[-1000:]
        plt.figure()
        plt.grid()
        plt.plot(time_beat, wt_beat, label=self.WT_legend)
        plt.plot(time_beat, drug_beat, label=self.drug_legend)
        plt.show()
        plt.savefig((self.directory + os.sep + 'cable_drug_compare.pdf'))
