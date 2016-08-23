'''
Created on Sep 11, 2012

@author: butler
'''

import numpy as np
import matplotlib.pyplot as plt


class PlotSchema():
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.columnwidth = 431
        self.ptperinch = 72.27
        self.figfractionalwidth = 2.0
        self.fontsize = 10
        self.linewidth = 2
        # setup fonts for figures - the rc interface is really quite clunky

    def set_column_width_pt(self, new_width):
        self.columnwidth = new_width

    def get_column_width_pt(self):
        return self.columnwidth

    def set_column_width_in(self, new_width):
        self.columnwidth = new_width * self.ptperinch

    def get_column_width_in(self):
        return self.columnwidth / self.ptperinch

    def set_frac_text_width(self, new_fraction):
        self.figfractionalwidth = new_fraction

    def get_frac_text_width(self):
        return self.figfractionalwidth

    def apply_figuresize_settings(self, figure):
        (width, height) = figure.get_size_inches()
        new_width = self.columnwidth / self.ptperinch * self.figfractionalwidth
        new_height = height / width * new_width
        figure.set_size_inches((new_width, new_height))

    def set_fontsize(self, fontsize):
        self.fontsize = fontsize
        pass

    def get_fontsize(self):
        return self.fontsize

    def set_linewidth(self, linewidth):
        self.linewidth = linewidth

    def get_linewidth(self):
        return self.linewidth

    def apply_fontsettings(self, pyplot):
        """
        \param
        """
        pyplot.rc('text', usetex=True)
        pyplot.rc('font', family='serif')
        pyplot.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]


    def apply_axes_lables_sizes(self, axes):
        for x_ticl_i in axes.get_xticklabels():
            x_ticl_i.set_fontsize(self.fontsize)
        for y_ticl_i in axes.get_yticklabels():
            y_ticl_i.set_fontsize(self.fontsize)

    def apply_legend_sizes(self, legend):
        for label in legend.get_texts():
            label.set_fontsize(self.fontsize)


def test_figure_class(test_fig_save="test_scaling.png"):
    x = np.linspace(-5, 5, 20)
    y = x ** 2
    figure = plt.figure()
    Alter_scaling = PlotSchema()
    Alter_scaling.apply_figuresize_settings(figure)
    Alter_scaling.apply_fontsettings(plt)
    plt.plot(x, y)
    plt.grid()
    plt.savefig(test_fig_save, bbox_inches='tight')
