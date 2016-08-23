'''
Created on 15/07/2013

@author: butler
'''
import argparse
import priming
import sys
import t_peak

if __name__ == "__main__":
    binary_description = ("Creates plots for the europace parameter space." 
    + "Currently only works for priming parameter space")

    parser = argparse.ArgumentParser(description=binary_description)
    parser.add_argument('path', help="Path to the eurospace directory")
    parser.add_argument('BCL', type=int, help="Basic cycle length (ms)")
    pt_help = "Select the method of plotting: overlay, single or subplot"
    parser.add_argument('plot_type', help=pt_help)
    bt_help = "Optionally specify the beat to be plotted"
    parser.add_argument("--beat", help=bt_help, type=int, default=0)
    parser.add_argument("--tpeak", action='store_true',
                        help="optionally calculate the t-wave peak and print")
    args = parser.parse_args()
    # Setup plotting
    plot_ctl = priming.Priming(args.path)
    assert(args.BCL in [500, 1000, 2000])
    if args.beat > 0:
        if args.plot_type == 'overlay':
            plot_ctl.beat_as_overlay(args.BCL, args.beat)
        elif args.plot_type == 'subplot':
            plot_ctl.beat_as_subplot(args.BCL, args.beat)
        elif args.plot_tupe == 'single':
            plot_ctl.beat_as_single(args.BCL, args.beat)
        else:
            print "Unrecognised plot type"
            sys.exit(1)
    else:
        if args.plot_type == 'overlay':
            plot_ctl.complete_as_overlay(args.BCL)
        elif args.plot_type == 'subplot':
            plot_ctl.complete_as_subplots(args.BCL)
        elif args.plot_tupe == 'single':
            plot_ctl.complete_as_single(args.BCL)
        else:
            print "Unrecognised plot type"
            sys.exit(1)
    if args.tpeak:
        peaker = t_peak.TPeak(args.path)
        peaker.create_list(args.BCL)
        peaker.calc_extrema(7)
