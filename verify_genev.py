import argparse
import glob
import itertools
from nuOscillation.tools import cache
import numpy
import ROOT
import math
import nuOscillation.plot.io
from nuOscillation.plot import histogram
import matplotlib.pyplot as plt 
import os
import scipy.stats

###############################################################################

class Record:
    pdg = "pdg"
    reactioncode = "reactioncode"
    enu = "enu"
    x = "x"
    y = "y"
    z = "z"
    
    ALL = (pdg, reactioncode, enu, x, y, z)
    
    def __call__(self, event):
        x, y, z = self._pos(event)
        return [event.StdHepPdg[0],
                event.NEneutmode,
                self._enu(event),
                x,
                y,
                z,
                ]
        
    def _pos(self, event):
        r = event.EvtVtx
        return (r[0],
                r[1],
                r[2],
                )
        
    def _convert2Darray(self, arr, nrows, ncols):
        out = numpy.zeros(shape=(nrows,ncols))
        for i in xrange(nrows):
            for j in xrange(ncols):
                out[i, j] = arr[(i * ncols) + j]
        return out
        
    def _enu(self, event):
        arr = self._convert2Darray(event.StdHepP4, event.StdHepN, 4)
        return arr[0][3]
        
    def indices(self):
        return [Record.pdg,
                Record.reactioncode,
                Record.enu,
                Record.x,
                Record.y,
                Record.z,
                ]

###############################################################################

class Reader:
    def __init__(self, filelist):
        self._filenames = list(itertools.chain.from_iterable(map(glob.glob, filelist)))
        self._func = Record()
        
    def read(self):
        indices = self._func.indices()
        uniquestr = "_".join(indices + self._filenames)
        cachetool = cache.CacheNumpy(uniquestr, prefix="tmp_verify_genev")
        if cachetool.exists():
            data = cachetool.read()
        else:
            data = self._load_data()
            cachetool.write(data)
        indices = dict(zip(indices, xrange(len(indices))))
        return indices, data
    
    def _load_data(self):
        data = []
        for fname in self._filenames:
            tfile = ROOT.TFile(fname, "READ")
            treename = "nRooTracker"
            tree = tfile.Get(treename)
            if not tree:
                raise Exception("Input file does not contain tree with name", treename, fname)
            for event in tree:
                data.append(self._func(event))
        print len(data), data[0]
        data = numpy.array(data)
        return data

###############################################################################

class PlotInteractionRate:
    def __init__(self, indices, data):
        self._data = data
        self._indices = indices
    
    def __iter__(self):
        #for p in self._plot_histograms():
        #    yield p
        for p in self._plot_xy():
            yield p
    
    def _plot_histograms(self):
        for var in Record.ALL:
            fig = plt.figure()
            d = self._data[:, self._indices[var]]
            y, xbinning = numpy.histogram(d)
            ax = fig.add_subplot(1, 1, 1)
            histogram.plot_hist_points(ax, xbinning, y)
            ax.set_xlabel(var)
            ax.set_ylabel("N events")
            name = os.sep.join(("hist_1d", var))
            yield name, fig
    
    def _plot_xy(self):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel(Record.x)
        ax.set_ylabel(Record.y)
        #Get data
        dx = self._data[:, self._indices[Record.x]]
        dy = self._data[:, self._indices[Record.y]]
        #Plot histogram
        heatmap, xedges, yedges = numpy.histogram2d(dx, dy, bins=20)
        extent = [min(xedges), max(xedges), min(yedges), max(yedges)]
        ax.imshow(heatmap, extent=extent, interpolation="nearest")
        dataset = numpy.row_stack((dx, dy))
        x = numpy.linspace(min(dx), max(dx))
        y = numpy.linspace(min(dy), max(dy))
        X, Y = numpy.meshgrid(x, y)
        kernel = scipy.stats.gaussian_kde(dataset)
        Z = kernel(numpy.row_stack((X.ravel(), Y.ravel()))).reshape(X.shape)
        ax.contour(X, Y, Z, color="black", linecolor="black")
        name = os.sep.join(("hist_2d", "xy"))
        yield name, fig

###############################################################################        

def plot_interaction_rate(args):
    reader = Reader(args.patterns)
    indices, data = reader.read()
    out = nuOscillation.plot.io.FigureWriter()
    plotter = PlotInteractionRate(indices, data)
    out(plotter)
    return

def parsecml():
    parser = argparse.ArgumentParser(description="Verify genev.")
    parser.add_argument("patterns", metavar="F", type=str, nargs='+', help="Input file list.")
    args = parser.parse_args()
    return args

def sanitize(args):
    for p in args.patterns:
        if len(glob.glob(p)) < 1:
            raise Exception("Input file does not exist or no matches to the pattern.", p)
    return

def main():
    #ROOT.gROOT.ProcessLine(".L rootrackerhelper.cxx++")
    args = parsecml()
    sanitize(args)
    plot_interaction_rate(args)
    return

if __name__ == "__main__":
    main()

