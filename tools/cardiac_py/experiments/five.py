'''
Created on Sep 24, 2012

@author: butler
'''
import numpy as np
import plot_schema
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib import rc
from matplotlib.mlab import griddata
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import LinearLocator
from matplotlib.transforms import  Bbox
import scipy.interpolate as interpolate
import matplotlib.colors as color
import os

class PlotECGs():
    '''
    classdocs
    '''
    def __init__(self, directory, full=True):
        '''
        Constructor
        '''
        self.t_gap_ms = 5.0
        self.directory = directory
        self.normal_dir_name = "run_14"
        self.modified_dir_name = 'run_15'
        self.normal_legend = '$I_{Kr} = 0.153$'
        self.modified_legend = '$I_{Kr} = 0.153$'
        self.full = full
        if self.full:
            self.electrodes = ['electrode#000448302','electrode#000451300','electrode#000452730','electrode#000453393','electrode#000457525','electrode#000458894',"electrode#000438028","electrode#000460291"]
        else:
            self.electrodes = ['electrode#000094150','electrode#000092294']
        self.period_ms = 1000
        self.template = plot_schema.PlotSchema()
        self.template.set_fontsize(16)
        self.fix_axes = True
        self.lfont = matplotlib.font_manager.FontProperties(size=(0.75 * self.template.fontsize))
        
        
    def load_data(self):
        n_data_path = self.directory + '/' + self.normal_dir_name
        self.n_data = self.load_case_data(n_data_path)
        m_data_path = self.directory + '/' + self.modified_dir_name 
        self.m_data = self.load_case_data(m_data_path)
        if self.fix_axes:
            #Presume normal is an ok max
            e_max = np.max(self.n_data[:,:])
            e_min = np.min(self.n_data[:,:])
            if self.full:
                l1_max = np.max(self.n_data[:,7] - self.n_data[:,6])
                l1_min = np.min(self.n_data[:,7] - self.n_data[:,6])
            else:
                l1_max = 0.0
                l2_min = 0.0
            self.e_max = max(e_max, l1_max)
            self.e_min = min(e_min, l1_min)
        
    def load_case_data(self, case_path):
        n_electrodes = len(self.electrodes)
        for ii in range(n_electrodes):
            electrode_path = case_path + '/' + self.electrodes[ii]
            temp_data = np.loadtxt(electrode_path)
            if ii == 0:
                n_times = temp_data.size
                data = np.zeros((n_times, n_electrodes))
            data[:,ii] = temp_data[:]
        return data
    
    
    def set_ECG_type(self,ECG_lead, flipper= -1):
        ''' SETUP plot vectors for each of the different ECG types
        '''
        if self.full:
            lead_1 = 7
            lead_2 = 8
        else:
            lead_1 = 1
            lead_2 = 2
        ECG_type = ECG_lead - 1
        
        assert((len(self.electrodes) >= ECG_type))
        n_points = self.n_data[:,0].size
        # 1 setup time series
        max_time = self.t_gap_ms / 1000 * n_points
        self.time = np.linspace(0, max_time, n_points)
        # first normalise
        if ECG_type == -1:
            col_1 = lead_1 - 1
            col_2 = lead_2 - 1
            self.n_ECG_data = flipper * (self.n_data[:,col_1] - self.n_data[:,col_2])
            self.m_ECG_data = flipper * (self.m_data[:,col_1] - self.m_data[:,col_2])
            self.y_label = r'$ \Delta V$'
        else:
            self.n_ECG_data = flipper * self.n_data[:,ECG_type]
            self.m_ECG_data = flipper * self.m_data[:,ECG_type]
            self.y_label = r'$ V$'
            
    def plot_normal_ECG_final(self, save_file):           

        """ """
        index_end = self.time.size -1
        index_start = self.time.size - 1 - (self.period_ms / self.t_gap_ms)
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.time[index_start:index_end], self.n_ECG_data[index_start:index_end])
        plt.xlabel('$t (s)$',fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label,fontsize=self.template.get_fontsize(), rotation='horizontal')
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=100, bbox_inches='tight')
        plt.show()
        
    def plot_normal_ECG_full(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.time, self.n_ECG_data)
        plt.xlabel('$t (s)$',fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label,fontsize=self.template.get_fontsize(), rotation='horizontal')
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=100, bbox_inches='tight')
        plt.show()
        
    def plot_modified_ECG_full(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.time, self.m_ECG_data)
        plt.xlabel('$t (s)$',fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label,fontsize=self.template.get_fontsize(), rotation='horizontal')
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=100, bbox_inches='tight')
        plt.show()
    def plot_modified_ECG_final(self, save_file):
        index_end = self.time.size -1
        index_start = self.time.size - 1 - (self.period_ms / self.t_gap_ms)
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.time[index_start:index_end], self.m_ECG_data[index_start:index_end])
        plt.xlabel('$t (s)$',fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label,fontsize=self.template.get_fontsize(), rotation='horizontal')
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=100, bbox_inches='tight')
        plt.show()
    def overlay_ECG_full(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        plt.plot(self.n_ECG_data,label=self.normal_legend)
        plt.plot(self.m_ECG_data,label=self.modified_legend)
        plt.xlabel('$t (s)$',fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label,fontsize=self.template.get_fontsize(), rotation='horizontal')
        plt.legend(prop=self.lfont)
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=100, bbox_inches='tight')
        plt.show()
        
    def overlay_ECG_rapid(self, save_file,colour='black',l1=False):
        index_end = self.time.size -1
        index_start = self.time.size - 1 - (self.period_ms / self.t_gap_ms)
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        if self.fix_axes:
            axes.set_ylim((self.e_min, self.e_max))
            axes.set_xlim((self.time[index_start]),self.time[index_end])
        axes.set_frame_on(False)
        plt.plot(self.time[index_start:index_end], self.m_ECG_data[index_start:index_end],label=self.modified_legend,color=colour,linestyle="--",linewidth=2)
        plt.plot(self.time[index_start:index_end], self.n_ECG_data[index_start:index_end],label=self.normal_legend,color=colour,linestyle="-",linewidth=2)
        axes.get_xaxis().set_visible(l1)
        axes.get_yaxis().set_visible(l1)
            
        self.template.apply_figuresize_settings(f)
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=100, bbox_inches='tight',transparent=True)
        
    def plot_for_rapid_fig_4(self):
        self.set_ECG_type(0, -1)
        self.overlay_ECG_rapid('lead_1.eps','black',False)
        self.set_ECG_type(1,1)
        self.overlay_ECG_rapid('V1.eps','red')
        self.set_ECG_type(2,1)
        self.overlay_ECG_rapid('V2.eps','green')
        self.set_ECG_type(3,1)
        self.overlay_ECG_rapid('V3.eps','brown')
        self.set_ECG_type(4,1)
        self.overlay_ECG_rapid('V4.eps','blue')
        self.set_ECG_type(5,1)
        self.overlay_ECG_rapid('V5.eps','orange')
        self.set_ECG_type(6,1)
        self.overlay_ECG_rapid('V6.eps','purple')
     
class Plot_G_NaL_Changes():
    '''
    classdocs
    '''
    def __init__(self, directory, full=True):
        '''
        Constructor
        '''
        self.t_gap_ms = 5.0
        self.directory = directory
        #self.g_NaL_dirs = ['g_NaL_0_15','g_NaL_0_20','g_NaL_0_25','g_NaL_0_26','g_NaL_0_27','g_NaL_0_28','g_NaL_0_29','g_NaL_0_3']
        self.g_NaL_dirs = ['g_NaL_0_15_a','g_NaL_0_20_a','g_NaL_0_25_a','g_NaL_0_275','g_NaL_0_29_a','g_NaL_0_2925','g_NaL_0_3_a']
        #self.legends = ['$g_{Nal}^{old} = 0.15$','$g_{Nal}^{old} = 0.20$','$g_{Nal}^{old} = 0.25$','$g_{Nal}^{old} = 0.26$','$g_{Nal}^{old} = 0.27$','$g_{Nal}^{old} = 0.28$','$g_{Nal}^{old} = 0.29$','$g_{Nal}^{old} = 0.30$']
        self.legends = ['$g_{Nal}^{new} = 0.15$','$g_{Nal}^{new} = 0.20$','$g_{Nal}^{new} = 0.25$','$g_{Nal}^{new} = 0.2725$','$g_{Nal}^{new} = 0.29$','$g_{Nal}^{new} = 0.2925$','$g_{Nal}^{new} = 0.30$']
        #self.g_NaL_dirs = ['g_NaL_0_28','g_NaL_0_29','g_NaL_0_3','g_NaL_0_29_a','g_NaL_0_2925','g_NaL_0_3_a']
        #self.legends = ['$g_{Nal}^{old} = 0.28$','$g_{Nal}^{old} = 0.29$','$g_{Nal}^{old} = 0.30$','$g_{Nal}^{new} = 0.29$','$g_{Nal}^{new} = 0.2925$','$g_{Nal}^{new} = 0.30$']

        self.case_type_dir_name = "normal_ikr/full"
        
        self.full = full
        if self.full:
            self.electrodes = ['electrode#000448302','electrode#000451300','electrode#000452730','electrode#000453393','electrode#000457525','electrode#000458894',"electrode#000438028","electrode#000460291"]
        else:
            self.electrodes = ['electrode#000094150','electrode#000092294']
        self.period_ms = 1000
        self.template = plot_schema.PlotSchema()
    def load_data(self):
        n_cases = len(self.g_NaL_dirs)
        self.ECG_data = np.zeros((0))
        ii = 0
        for g_NaL_dir in self.g_NaL_dirs:
            print g_NaL_dir
            data_path = self.directory + os.sep + g_NaL_dir + os.sep + self.case_type_dir_name
            case_data = self.load_case_data(data_path)
            (n_x, n_y) = case_data.shape
            print case_data.shape
            if self.ECG_data.size == 0:
                self.ECG_data = np.zeros((n_cases,n_x + 5,n_y))
                print 'shape changed'
                
            self.ECG_data[ii,0:n_x,:] = case_data[:,:]
            ii = ii + 1
    def load_case_data(self, case_path):
        n_electrodes = len(self.electrodes)
        for ii in range(n_electrodes):
            electrode_path = case_path + '/' + self.electrodes[ii]
            temp_data = np.loadtxt(electrode_path)
            if ii == 0:
                n_times = temp_data.size
                data = np.zeros((n_times, n_electrodes))
            data[:,ii] = temp_data[:]
        return data
    
    
    def set_ECG_type(self,ECG_lead, flipper= 1):
        ''' SETUP plot vectors for each of the different ECG types
        '''
        if self.full:
            lead_1 = 7
            lead_2 = 8
        else:
            lead_1 = 1
            lead_2 = 2
        ECG_type = ECG_lead - 1
        
        assert((len(self.electrodes) >= ECG_type))
        n_points = self.ECG_data[0,:,0].size
        # 1 setup time series
        max_time = self.t_gap_ms / 1000 * n_points
        self.time = np.linspace(0, max_time, n_points)
        # first normalise
        self.ECG_plot_data = np.zeros((len(self.g_NaL_dirs),n_points))
        if ECG_type == -1:
            for ii in range(len(self.g_NaL_dirs)):
                col_1 = lead_1 - 1
                col_2 = lead_2 - 1
                self.ECG_plot_data[ii,:] = flipper * (self.ECG_data[ii,:,col_1] - self.ECG_data[ii,:,col_2])
            self.y_label = r'$ \Delta V$'
        else:
            assert(1 == -1)
            self.n_ECG_data = flipper * self.n_data[:,ECG_type]
            self.m_ECG_data = flipper * self.m_data[:,ECG_type]
            self.y_label = r'$ V$'
            

        
    def overlay_ECG_full(self, save_file):
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        for ii in range(len(self.g_NaL_dirs)):
            legend_a = self.legends[ii]
            plt.plot(self.time, self.ECG_plot_data[ii,:],label=legend_a)
            
        plt.xlabel('$t (s)$',fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label,fontsize=self.template.get_fontsize(), rotation='horizontal')
        plt.legend(prop=self.lfont)
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=100, bbox_inches='tight')
        plt.show()
    def overlay_ECG_final(self, save_file):
        index_end = self.time.size -1
        index_start = self.time.size - 1 - (self.period_ms / self.t_gap_ms)
        self.template.apply_fontsettings(plt)
        f = plt.figure()
        self.template.apply_figuresize_settings(f)
        axes = plt.axes()
        for ii in range(len(self.g_NaL_dirs)):
            legend_a = self.legends[ii]
            plt.plot(self.time[index_start:index_end], self.ECG_plot_data[ii,index_start:index_end],label=legend_a)
        plt.xlabel('$t (s)$',fontsize=self.template.get_fontsize())
        plt.ylabel(self.y_label,fontsize=self.template.get_fontsize(), rotation='horizontal')
        plt.legend(prop=self.lfont)
        self.template.apply_figuresize_settings(f)
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.template.get_fontsize())
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.template.get_fontsize())
        save_loc = self.directory + '/' + save_file
        plt.savefig(save_loc, dpi=100, bbox_inches='tight')