import argparse
import glob
import itertools
from nuOscillation.tools import cache, progress
import numpy
import ROOT
import math
import nuOscillation.plot.io
import nuOscillation.model
from nuOscillation.plot import histogram
import matplotlib.pyplot as plt 
import os
import scipy.stats
import matplotlib
from nuOscillation.model import constants, runtime

###############################################################################

class Record:
    pdg = "pdg"
    reactioncode = "reactioncode"
    enu = "enu"
    x = "x"
    y = "y"
    z = "z"
    weight = "weight"
    ndid = "ndid"
    
    ALL = (pdg, reactioncode, enu, x, y, z, weight)
    
    _axis_labels = {pdg : "PDG code",
                    reactioncode : "NEUT reaction code",
                    enu : r"$E_{\nu}$ [GeV]",
                    x : "$x$ [m]",
                    y : "$y$ [m]",
                    z : "$z$ [m]",
                    weight : "event weight",
                    ndid : "near detector ID",
                    }

    @classmethod
    def label(cls, field):
        return cls._axis_labels[field]

def _expandfilelist(filelist):
    return list(itertools.chain.from_iterable(map(glob.glob, filelist)))

###############################################################################

class NeutEventReader:
    
    def __call__(self, event):
        x, y, z = self._pos(event)
        return [event.StdHepPdg[0],
                event.NEneutmode,
                self._enu(event),
                x,
                y,
                z,
                1.0,
                0.0,
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
                Record.weight,
                Record.ndid,
                ]


###############################################################################

class FluxEventReader:
    
    def __call__(self, event):
        return [event.pdg,
                0.0,
                event.enu,
                event.xnu,
                event.ynu,
                0.0,
                event.norm,
                event.ndid,
                ]
        
    def indices(self):
        return [Record.pdg,
                Record.reactioncode,
                Record.enu,
                Record.x,
                Record.y,
                Record.z,
                Record.weight,
                Record.ndid,
                ]

###############################################################################

class Reader(object):
    def __init__(self, filelist, treename, func):
        self._filenames = _expandfilelist(filelist)
        self._func = func
        self._treename = treename
        
    def read(self):
        indices = self._func.indices()
        uniquestr = "_".join(indices 
                             + self._filenames
                             + [type(self._func).__name__]
                             )
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
            treename = self._treename
            tree = tfile.Get(treename)
            if not tree:
                raise Exception("Input file does not contain tree with name", treename, fname)
            for event in tree:
                data.append(self._func(event))
        print len(data), data[0]
        data = numpy.array(data)
        return data

###############################################################################

class NeutFileReader(Reader):
    def __init__(self, filelist):
        super(NeutFileReader, self).__init__(filelist, "nRooTracker", NeutEventReader())

###############################################################################

class FluxFileReader(Reader):
    def __init__(self, filelist):
        super(FluxFileReader, self).__init__(filelist, "flux", FluxEventReader())
        
    def _load_data(self):
        data = []
        for fname in self._filenames:
            tfile = ROOT.TFile(fname, "READ")
            treename = self._treename
            tree = tfile.Get(treename)
            if not tree:
                raise Exception("Input file does not contain tree with name", treename, fname)
            prescale = 20
            iterevents = progress.printProgress("FluxFileReader", tree.GetEntries()/prescale, self._prescale(tree, prescale))
            for event in iterevents:
                data.append(self._func(event))
        print len(data), data[0]
        data = numpy.array(data)
        return data
    
    def _prescale(self, tree, prescale):
        n = tree.GetEntries()
        for i in xrange(n):
            if i % prescale == 0:
                tree.GetEntry(i)
                yield tree

###############################################################################

class PlotInteractionRate:
    def __init__(self, outdir, indices, data, scale_to_m=1.0):
        self._scale_to_m = scale_to_m
        self._outdir = outdir
        self._data = data
        self._indices = indices
    
    def __iter__(self):
        for p in self._plot_xy_hist():
            yield p
        for p in self._plot_histograms():
            yield p
    
    def _plot_histograms(self):
        for var in Record.ALL:
            fig = plt.figure()
            d = self._data[:, self._indices[var]]
            weights = self._data[:, self._indices[Record.weight]]
            y, xbinning = numpy.histogram(d, weights=weights)
            ax = fig.add_subplot(1, 1, 1)
            histogram.plot_hist_points(ax, xbinning, y)
            ax.set_xlabel(Record.label(var))
            ax.set_ylabel("N events")
            name = os.sep.join((self._outdir, "hist_1d", var))
            yield name, fig
    
    def _plot_xy_hist(self):
        fig = plt.figure()
        name = os.sep.join((self._outdir, "hist_2d", "xy"))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel(Record.label(Record.x))
        ax.set_ylabel(Record.label(Record.y))
        #Get data
        dx = self._data[:, self._indices[Record.x]]
        dy = self._data[:, self._indices[Record.y]]
        #scale dx and dy
        dx *= self._scale_to_m
        dy *= self._scale_to_m
        weights = self._data[:, self._indices[Record.weight]]
        #Plot data
        ret = ax.hist2d(dx, dy, bins=25, weights=weights, cmap=matplotlib.cm.get_cmap("hot"))
        img = ret[-1]
        fig.colorbar(img)
        yield name, fig
    
#     def _plot_xy_kernel(self):
#         fig = plt.figure()
#         ax = fig.add_subplot(1, 1, 1)
#         ax.set_xlabel(Record.label(Record.x))
#         ax.set_ylabel(Record.label(Record.y))
#         #Get data
#         dx = self._data[:, self._indices[Record.x]]
#         dy = self._data[:, self._indices[Record.y]]
#         #Plot histogram
#         heatmap, xedges, yedges = numpy.histogram2d(dx, dy, bins=20)
#         extent = [min(xedges), max(xedges), min(yedges), max(yedges)]
#         ax.imshow(heatmap, extent=extent, interpolation="nearest")
#         dataset = numpy.row_stack((dx, dy))
#         x = numpy.linspace(min(dx), max(dx))
#         y = numpy.linspace(min(dy), max(dy))
#         X, Y = numpy.meshgrid(x, y)
#         kernel = scipy.stats.gaussian_kde(dataset)
#         Z = kernel(numpy.row_stack((X.ravel(), Y.ravel()))).reshape(X.shape)
#         ax.contour(X, Y, Z, color="black", linecolor="black")
#         name = os.sep.join(("hist_2d", "xy"))
#         yield name, fig

###############################################################################

class PlotFluxIntXsec:
    def __init__(self, neutindices, neutdata, fluxindices, fluxdata, scale_to_m = 1.0/100.0):
        self._neut_indices = neutindices
        self._neut_data = neutdata
        self._fluxindices = fluxindices
        self._fluxdata = fluxdata
        self._scale_to_m = scale_to_m
    
    def __iter__(self):
        return self._plot()
        
    def _plot(self):
        yield figname, fig
            

###############################################################################        

def plot_interaction_rate(args):
    reader = NeutFileReader(args.patterns)
    indices, data = reader.read()
    out = nuOscillation.plot.io.FigureWriter()
    plotter = PlotInteractionRate("evgen", indices, data)
    out(plotter)
    return

###############################################################################

def plot_flux(args):
    ntuple = nuOscillation.model.beam.FluxNtuple(runtime.getcontext().beamcontext)
    reader = FluxFileReader([ntuple.filename()])
    indices, data = reader.read()
    
    DetectorId = nuOscillation.model.beam.constants.DetectorId
    ndid = DetectorId.toint(DetectorId.ND2K)
    #filter out particular near detector
    data = data[data[:, indices[Record.ndid]] == ndid]
    out = nuOscillation.plot.io.FigureWriter()
    cm_to_m = 1.0/100.0
    plotter = PlotInteractionRate("flux", indices, data, scale_to_m = cm_to_m)
    out(plotter)

###############################################################################

def plot_flux_int_xsec(args):
    #load neut data
    reader = NeutFileReader(args.patterns)
    neutindices, neutdata = reader.read()
    #load flux data
    ntuple = nuOscillation.model.beam.FluxNtuple(runtime.getcontext().beamcontext)
    reader = FluxFileReader([ntuple.filename()])
    fluxindices, fluxdata = reader.read()
    DetectorId = nuOscillation.model.beam.constants.DetectorId
    ndid = DetectorId.toint(DetectorId.ND2K)
    #filter out particular near detector
    fluxdata = fluxdata[fluxdata[:, fluxindices[Record.ndid]] == ndid]
    #make plot
    out = nuOscillation.plot.io.FigureWriter()
    cm_to_m = 1.0 / 100.0
    plotter = PlotFluxIntXsec(neutindices, neutdata, fluxindices, fluxdata, scale_to_m = cm_to_m)
    out(plotter) 

###############################################################################

def parsecml():
    parser = argparse.ArgumentParser(description="Verify genev.")
    parser.add_argument("patterns", metavar="M", type=str, nargs='+', help="Input file list.")
    #parser.add_argument("-f", "--flux", metavar="F", dest="fluxpatterns", type=str, nargs='+', help=".")
    parser.add_argument("-p", "--plot", metavar="P", dest="plot", type=str, choices=["flux", "rate"], help="Choose which plot to make.")
    args = parser.parse_args()
    return args

###############################################################################

def sanitize(args):
    for p in args.patterns:
        if len(glob.glob(p)) < 1:
            raise Exception("Input file does not exist or no matches to the pattern.", p)
    return

###############################################################################

def main():
    #ROOT.gROOT.ProcessLine(".L rootrackerhelper.cxx++")
    args = parsecml()
    sanitize(args)
    if args.plot == "flux":
        plot_flux(args)
    elif args.plot == "rate":
        plot_interaction_rate(args)
    return

###############################################################################

if __name__ == "__main__":
    main()

###############################################################################
