'''
Created on 08/01/2013

@author: butler
'''
from . import ECG_plots
import os


class Experiment():

    def __init__(self, root_dir="/Users/butler/Desktop/CARDIOID/ECG_Data"):
        self.root_dir = root_dir
        self.experiment_dir = "experiment_14"
        self.period_ms = 2000
        self.stimulus_time_ms = 20510

    def plot_ECGs_for_object_data_set_g_NaL(self):
        """
        This plots two runs where we try to trigger reentry with & without I_k
         block.
        Here the g_NaL parameters were set using the object.data clone
        properties
        """
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.normal_dir_name = "run_5"
        plotter.normal_legend = "$I_{Kr} = 0.15$"
        plotter.modified_dir_name = "run_6"
        plotter.modified_legend = "$I_{Kr} = 0.0$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_ECG_full("Lead_1_object_data_set_g_NaL.eps")

    def plot_ECGs_for_binary_set_g_NaL(self):
        """
        This plots two runs where we try to trigger reentry with & without I_kr
        block.
        Here the g_NaL parameters were set using source code only
        """
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.normal_dir_name = "run_9"
        plotter.normal_legend = "$I_{Kr} = 0.15$"
        plotter.modified_dir_name = "run_10"
        plotter.modified_legend = "$I_{Kr} = 0.0$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_ECG_full("Lead_1_source_set_g_NaL.eps")

    def plot_ECG_overlay_set_beats(self):
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.normal_dir_name = "run_9"
        plotter.normal_legend = "$I_{Kr} = 0.15$"
        plotter.modified_dir_name = "run_10"
        plotter.modified_legend = "$I_{Kr} = 0.0$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.period_ms = 2000
        plotter.overlay_ECG_beat("Beat_9.eps", 9)
        plotter.overlay_ECG_beat("Beat_10.eps", 10)

    def plot_ECGs_for_r1310_comparison(self):
        """
        This plots two runs where we try to trigger reentry with & without I_kr
        block.
        Here the g_NaL parameters were set using source code only
        """
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.normal_dir_name = "run_10"
        plotter.normal_legend = "Precompiled"
        plotter.modified_dir_name = "run_13"
        plotter.modified_legend = "r1310"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_ECG_full("r1310_comparison.eps")

    def plot_fast_pacing(self):
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.normal_dir_name = "run_11"
        plotter.normal_legend = "$I_{Kr} = 0.15$"
        plotter.modified_dir_name = "run_12"
        plotter.modified_legend = "$I_{Kr} = 0.0$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_ECG_full("Lead_1_1hz.eps")

    def plot_fast_pacing_single_beat(self):
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.normal_dir_name = "run_11"
        plotter.normal_legend = "$I_{Kr} = 0.15$"
        plotter.modified_dir_name = "run_12"
        plotter.modified_legend = "$I_{Kr} = 0.0$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_ECG_beat("Beat_9_1_hz.eps", 9)

    def plot_ECG_0_1_mm(self):
        """
        This plots two runs where we try to trigger reentry with & without I_kr
        block.
        Here the g_NaL parameters were set using source code only
        """
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.normal_dir_name = "run_7"
        plotter.normal_legend = "$I_{Kr} = 0.00$"
        plotter.modified_dir_name = "run_7"
        plotter.modified_legend = "$I_{Kr} = 0.0$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_ECG_full("Lead_1_source_set_01mm.eps")

    def plot_run_15(self):
        """
        This plots two runs where we try to trigger reentry with & without I_kr
        block.
        Here the g_NaL parameters were set using source code only
        """
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.normal_dir_name = "run_15"
        plotter.normal_legend = "$I_{Kr} = 0.00$"
        plotter.modified_dir_name = "run_15"
        plotter.modified_legend = "$I_{Kr} = 0.0$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_ECG_full("Plot_run_15.eps")

    def plot_ECG_0_1_m_old_mesh(self):
        """
        This plots two runs where we try to trigger reentry with & without I_kr
        block.
        Here the g_NaL parameters were set using source code only
        """
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.normal_dir_name = "run_20"
        plotter.normal_legend = "$I_{Kr} = 0.153$"
        plotter.modified_dir_name = "run_21"
        plotter.modified_legend = "$I_{Kr} = 0.0$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_ECG_full("OLD_MESH.eps")

    def plot_ECG_0_1a_m_new(self):
        """
        This plots two runs where we try to trigger reentry with & without I_kr
        block.
        Here the g_NaL parameters were set using source code only
        """
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.normal_dir_name = "run_22"
        plotter.normal_legend = "$I_{Kr} = 0.153$"
        plotter.modified_dir_name = "run_23"
        plotter.modified_legend = "$I_{Kr} = 0.0$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_ECG_full("new_0.1mm_mesh.eps")
        print 'This was for runs: ', plotter.normal_dir_name
        print 'and: ', plotter.modified_dir_name

    def plot_ECG_0_2_0_1_compare(self):
        """
        This plots two runs where we try to trigger reentry with & without I_kr
        block.
        Here the g_NaL parameters were set using source code only
        """
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.normal_dir_name = "run_5"
        plotter.normal_legend = "$\Delta x = 0.2$"
        plotter.modified_dir_name = "run_22"
        plotter.modified_legend = "$\Delta x = 0.1$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_plus_difference("resolution_comparison.eps")
        print 'This was for runs: ', plotter.normal_dir_name
        print 'and: ', plotter.modified_dir_name

    def plot_ECG_0_2_0_1_compare_beat_all_leads(self):
        """
        This plots two runs where we try to trigger reentry with & without I_kr
        block.
        Here the g_NaL parameters were set using source code only
        """
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.normal_dir_name = "run_5"
        plotter.normal_legend = "$\Delta x = 0.2$"
        plotter.modified_dir_name = "run_22"
        plotter.modified_legend = "$\Delta x = 0.1$"
        plotter.period_ms = self.period_ms
        plotter.load_data()

        for ii in range(7):
            if ii == 0:
                plotter.set_ECG_type(0, -1)
            else:
                plotter.set_ECG_type(ii, 1)
            plotter.overlay_plus_perc_error_beat("frac_error_beat_9_lead_" +
                                                 repr(ii) + ".pdf", 9)
            plotter.overlay_plus_perc_error_beat("frac_error_beat_10_lead_" +
                                                 repr(ii) + ".pdf", 10)
            plotter.overlay_plus_perc_error_beat("frac_error_beat_11_lead_" +
                                                 repr(ii) + ".pdf", 11)
            plotter.overlay_plus_perc_error_beat("frac_error_beat_12_lead_" +
                                                 repr(ii) + ".pdf", 12)
        print 'This was for runs: ', plotter.normal_dir_name
        print 'and: ', plotter.modified_dir_name

    def plot_ECG_0_2_0_1_compare_error(self):
        """
        This plots two runs where we try to trigger reentry with & without I_kr
        block. Here the g_NaL parameters were set using source code only
        """
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.normal_dir_name = "run_5"
        plotter.normal_legend = "$\Delta x = 0.2$"
        plotter.modified_dir_name = "run_22"
        plotter.modified_legend = "$\Delta x = 0.1$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_plus_perc_error(
                                    "resolution_comparison_as_perc_error.eps")
        print 'This was for runs: ', plotter.normal_dir_name
        print 'and: ', plotter.modified_dir_name

    def plot_1hz_vs_2hz_ikr_normal(self):
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.period_ms = self.period_ms
        plotter.normal_dir_name = "run_5"
        plotter.normal_legend = "$0.5 Hz$"
        plotter.modified_dir_name = "run_11"
        plotter.modified_legend = "$1 hz$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_ECG_beat("freq_comparison_ikr_normal_a.eps", 2)
        plotter.overlay_ECG_beat("freq_comparison_ikr_normal_b.eps", 3)

    def plot_1hz_vs_2hz_ikr_zeroed(self):
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.period_ms = self.period_ms
        plotter.normal_dir_name = "run_6"
        plotter.normal_legend = "$0.5 Hz$"
        plotter.modified_dir_name = "run_12"
        plotter.modified_legend = "$1 hz$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_ECG_beat("freq_comparison_ikr_zeroed_a.eps", 2)
        plotter.overlay_ECG_beat("freq_comparison_ikr_zeroed_b.eps", 3)

    def plot_rentry_timing(self):
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.period_ms = self.period_ms
        plotter.normal_dir_name = "run_5"
        plotter.normal_legend = "$i_{kr} = 0.153$"
        plotter.modified_dir_name = "run_6"
        plotter.modified_legend = "$i_{kr} = 0.0$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_ECG_reent("Reentry_intiation_capture.eps", 11)

    def plot_g_NaL_0_25(self):
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.normal_dir_name = "run_24"
        plotter.normal_legend = "$I_{Kr} = 0.153$"
        plotter.modified_dir_name = "run_25"
        plotter.modified_legend = "$I_{Kr} = 0.0$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_ECG_full("Lead_g_NaL_0_25.eps")
        print 'This was for runs: ', plotter.normal_dir_name
        print 'and: ', plotter.modified_dir_name

    def plot_g_NaL_0_15(self):
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.normal_dir_name = "run_26"
        plotter.normal_legend = "$I_{Kr} = 0.153$"
        plotter.modified_dir_name = "run_27"
        plotter.modified_legend = "$I_{Kr} = 0.0$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_ECG_full("Lead_g_NaL_0_15.eps")
        print 'This was for runs: ', plotter.normal_dir_name
        print 'and: ', plotter.modified_dir_name

    def plot_delayed_S2(self):
        experiment_dir = self.root_dir + os.sep + self.experiment_dir
        plotter = ECG_plots.PlotECGs(experiment_dir)
        plotter.normal_dir_name = "run_28"
        plotter.normal_legend = "$I_{Kr} = 0.153$"
        plotter.modified_dir_name = "run_29"
        plotter.modified_legend = "$I_{Kr} = 0.0$"
        plotter.load_data()
        plotter.set_ECG_type(0, -1)
        plotter.overlay_ECG_full("Lead_1_delayed_s2.eps")
        print 'This was for runs: ', plotter.normal_dir_name
        print 'and: ', plotter.modified_dir_name
