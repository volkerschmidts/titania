#!/usr/bin/env python3
import sys
import os
import subprocess
from select import poll
from time import sleep

import matplotlib.pylab as plt
import matplotlib

monte_carlo_map=["rmsd(S)", "rmsd(sigm[S])", "rmsd(p)", "rmsd(sigm[p])"]

class DataMonitor():
    def __init__(self, sub_p, args):
        self.load_keywords(args)
        self.initialize_plot()
        self.sub_p=sub_p
        self.poll_obj = poll()
        self.poll_obj.register(self.sub_p.stdout)
        self.converged=False

    def load_keywords(self,args):
        self.timeout=args.timeout
        self.timestep=args.timestep
        self.nowait=args.nowait
        self.keywords={"keys": ["Finished iteration"], "position": [4] }
        if args.Soverall:
            self.keywords["keys"].append("S(overall)")
            self.keywords["position"].append(2)
        if args.MonteCarlo:
            self.keywords["keys"].append("Monte Carlo")
            self.keywords["position"].append(4)
        if args.Damping:
            self.keywords["keys"].append("Overall damping")
            self.keywords["position"].append(3)

    def initialize_plot(self):
        self.graph=dict()
        plt.ion()
        for idx, key in enumerate(self.keywords["keys"]):
            if key == "Finished iteration": 
                self.graph[key]={"Figure": None, "Axis": None, "x-value": [], "y-value": [], "index": 0}
            elif key == "Monte Carlo":
                self.graph[key]={"Figure": plt.figure(), "Axis": [], "x-value": [], "y-values": {"rmsd(S)": [], "rmsd(sigm[S])": [], "rmsd(p)": [], "rmsd(sigm[p])": []}, "index": 0}
                curr_g=self.graph[key]
                for i in range(221,225):
                    curr_g["Axis"].append(curr_g["Figure"].add_subplot(i))
                curr_g["Figure"].canvas.set_window_title(key)
            else:
                self.graph[key]={"Figure": plt.figure(), "Axis": None, "x-value": [], "y-value": [], "index": 0}
                curr_g=self.graph[key]
                curr_g["Axis"]=curr_g["Figure"].add_subplot()
                curr_g["Axis"].plot(self.graph["Finished iteration"]["x-value"], curr_g["y-value"])
                curr_g["Figure"].canvas.set_window_title(key)

    def save_data(self, key, value, **kwargs):
        curr_g=self.graph[key]
        if key == "Monte Carlo":
            if kwargs["index"]==0:
                curr_g["index"]=curr_g["index"]+1
            curr_g["y-values"][monte_carlo_map[kwargs["index"]]].append(value)
        else:
            curr_g["y-value"].append(value)
            curr_g["index"]=curr_g["index"]+1
        if ( key == "Finished iteration"):
            curr_g["x-value"].append(value)
            if len(curr_g["x-value"]) >=2 and curr_g["x-value"][-1] < curr_g["x-value"][-2]:
                self.clear_data()
            self.update_plots()

    def update_plots(self):
        if len(plt.get_fignums()) == 0:
            exit()
        for key in self.keywords["keys"]:
            curr_g=self.graph[key]
            if key == "Finished iteration":
                pass
            elif key == "Monte Carlo":
                if len(self.graph["Finished iteration"]["x-value"]) < 2:
                    pass
                elif plt.fignum_exists(curr_g["Figure"].number):
                    if len(self.graph["Finished iteration"]["x-value"]) != len(curr_g["y-values"]["rmsd(S)"]):
                        self.clear_data()
                    for idx in range(4):
                            ax=curr_g["Axis"][idx]
                            ax.cla()
                            ax.semilogy(self.graph["Finished iteration"]["x-value"], curr_g["y-values"][monte_carlo_map[idx]])
                            curr_g["Axis"][idx].title.set_text(monte_carlo_map[idx])
                    curr_g["Figure"].tight_layout()
            else:
                if plt.fignum_exists(curr_g["Figure"].number):
                    if len(self.graph["Finished iteration"]["x-value"]) != len(curr_g["y-value"]):
                        self.clear_data()
                    curr_g["Axis"].cla()
                    curr_g["Axis"].plot(self.graph["Finished iteration"]["x-value"], curr_g["y-value"])
            unpoping_pause(1e-3)

    def clear_data(self):
        self.converged=False
        for key in self.keywords["keys"]:
            curr_g=self.graph[key]
            if key == "Monte Carlo":
                for i in range(4):
                    curr_g["y-values"][monte_carlo_map[i]]=[]
                curr_g["x-value"]=[]
            else:
                curr_g["y-value"]=[]
                curr_g["x-value"]=[]

    def start(self):
        try:
            timer=0
            while True:
                line = self.sub_p.stdout.readline().decode("utf-8")
                if (self.poll_obj.poll(1) or line != "") and not self.converged:
                    if "REACHED" in line or "YES" in line:
                        self.converged = True;
                    for idx, key in enumerate(self.keywords["keys"]):
                        if key == "Monte Carlo":
                            if "rmsd(S)" in line:
                                value=float(line.strip().split()[self.keywords["position"][idx]])
                                self.save_data(key, value, index=0)
                            elif "rmsd(sigm[S])" in line:
                                value=float(line.strip().split()[self.keywords["position"][idx]])
                                self.save_data(key, value, index=1)
                            elif "rmsd(p)" in line:
                                value=float(line.strip().split()[self.keywords["position"][idx]])
                                self.save_data(key, value, index=2)
                            elif "rmsd(sigm[p])" in line:
                                value=float(line.strip().split()[self.keywords["position"][idx]])
                                self.save_data(key, value, index=3)
                        elif key in line:
                            value=float(line.strip().split()[self.keywords["position"][idx]])
                            self.save_data(key, value)
                    timer=0
                elif not self.converged:
                    timer=timer+self.timestep
                    if timer >= self.timeout:
                        print("Time out of monitor_TITANIA.py after {0} seconds...".format(self.timeout))
                        print("Use the argument -t=n (seconds) to change the timeout limit.")
                        exit()
                    sleep(self.timestep)
                elif self.converged:
                    if self.nowait:
                        exit()
                    else:
                        print("Proces finished and converged.")
                        print("Press x to exit or r to start a new run on the same file.")
                        mode = input('Keyboard Input:')
                        if "x" in mode:
                            exit()
                        elif "r" in mode:
                            self.converged=False
                            pass
                        else:
                            self.converged=False
                            pass
        except KeyboardInterrupt:
            exit()

def unpoping_pause(interval):
    backend = plt.rcParams['backend']
    if backend in matplotlib.rcsetup.interactive_bk:
        figManager = matplotlib._pylab_helpers.Gcf.get_active()
        if figManager is not None:
            canvas = figManager.canvas
            if canvas.figure.stale:
                canvas.draw()
            canvas.start_event_loop(interval)
            return

def monitor_argument_parser(argv):
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        "-i",
        nargs=1,
        default=None,
        help="Input file of the TITANIA job to monitor.",
        required=True,
    )
    parser.add_argument(
        "--output",
        "-o", 
        nargs=1,
        default=None, 
        help="Output file name.",
        required=False,
    )
    parser.add_argument(
        "--timeout",
        "-t",
        type=int,
        default=60,
        help="Limit of time waited for file updates.",
        required=False,
    )
    parser.add_argument(
        "--timestep",
        "-s",
        type=int,
        default=1,
        help="Time step to check output file.",
        required=False,
    )
    parser.add_argument(
        "--Soverall",
        action="store_true",
        help="Enables monitoring of S(overall)." 
    )
    parser.add_argument(
        "--MonteCarlo",
        action="store_true",
        help="Enables monitoring of Monte-Carlo stop criteria." 
    )
    parser.add_argument(
        "--Damping",
        action="store_true",
        help="Enables monitoring of redundant internal coordinate damping factors." 
    )
    parser.add_argument(
        "--nowait",
        action="store_true",
        help="Prevents monitor_TITANIA.py from staying open. This is recommenden when calling this from a script." 
    )
    return parser.parse_args(argv)

def main(argv):
    args=monitor_argument_parser(argv)
    f = subprocess.Popen(['tail','-f',args.input[0]], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    dm = DataMonitor(f, args)
    dm.start()


if __name__ == "__main__":
    main(sys.argv[1:])
