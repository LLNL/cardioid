from .. import ECG_plots
import os


class Experiment():
    """
    hello world
    """

    def __init__(self, current,
          sim_dir="/Users/butler/Desktop/CARDIOID_SEQ/ecg_data/experiment_16"):
        self.root_dir = sim_dir
        self.june = 'june_BCL_2000'
        self.current = 'BCL_2000'
        if current:
            self.experiment_dir = self.current
        else:
            self.experiment_dir = self.june
        self.period_ms = 2000
        self.stimulus_time_ms = 20510

    def plot_A_s2_ECG(self, save_name, ikr_val, s2_run_num):
        """
        This plots two runs where we try to trigger reentry with & without I_kr
        block. Here the g_NaL parameters were set using the object.data clone
        properties
        """

        ikr_dir = self.root_dir + os.sep + self.experiment_dir + os.sep + \
        "ikr_block_" + repr(ikr_val)
        print 'Getting ECG data from: ' + ikr_dir
        plotter = ECG_plots.PlotECGs(ikr_dir)
        #set modified and normal ECG to same s2 dirs
        if s2_run_num > 0:
            s2_dir = 's2_run_' + repr(s2_run_num)
        else:
            s2_dir = 'priming'
        plotter.normal_dir_name = s2_dir
        plotter.normal_legend = "$ikr block " + repr(ikr_val) + "\%$"
        plotter.modified_dir_name = s2_dir
        plotter.modified_legend = "$ikr block " + repr(ikr_val) + "\%$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.write_normal_lead_csv("ecg.txt")
        plotter.plot_normal_ECG_full("lead_1_full.eps")
