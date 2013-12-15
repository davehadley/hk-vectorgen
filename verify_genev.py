import argparse
import glob
import itertools
from nuOscillation.tools import cache, progress
import numpy
import ROOT
import math
import simplot.io
import nuOscillation.model
from simplot import histogram
import matplotlib.pyplot as plt 
import os
import scipy.stats
import matplotlib
from nuOscillation.model import constants, runtime
import collections
import re
import scipy.interpolate

###############################################################################

GenEvFile = collections.namedtuple("GenEvFile", ["polarity", "ndid", "r", "z", "n", "orientation"])
BeamFile = collections.namedtuple("BeamFile", ["polarity"])

###############################################################################

class FileNameParser:
    def __init__(self, patterns):
        self._flist = _expandfilelist(patterns)
        
        
    def __call__(self):
        r = None
        fname = self._flist[0]
        for f in [self._parse_genev, self._parse_beam]:
            r = f(fname)
            if r is not None:
                break
        #check for success
        if r is None:
            raise ValueError("cannot parse ", fname)
        return r
    
    def _parse_genev(self, fname):
        pat = "genev_(.*?)_(.*?)_r(.*?)_z(.*?)_(.)_(.*?).root"
        match = re.search(pat, fname)
        if not match:
            return None
        polarity = match.group(1)
        det = match.group(2)
        ndid = nuOscillation.model.constants.DetectorId.tostring(det)
        r = float(int(match.group(3))) / 100.0
        z = float(int(match.group(4))) / 100.0
        orientation = match.group(5)
        n = match.group(6)
        return GenEvFile(polarity=polarity,
                  ndid=ndid,
                  r=r,
                  z=z,
                  n=n,
                  orientation=orientation,
                  )
        
    def _parse_beam(self, fname):
        pat = "merged_beamfiles_(.*?).root"
        match = re.search(pat, fname)
        if not match:
            return None
        polarity = match.group(1)
        return BeamFile(polarity=polarity)

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
    def __init__(self, filelist, treename, func, prescale=1):
        self._filenames = _expandfilelist(filelist)
        self._func = func
        self._treename = treename
        self._prescale = prescale
        
    def read(self):
        indices = self._func.indices()
        uniquestr = "_".join(indices 
                             + self._filenames
                             + [type(self._func).__name__,
                                str(self._prescale),
                                ]
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
    
    def _iterprescale(self, tree, prescale):
        n = tree.GetEntries()
        for i in xrange(n):
            if i % prescale == 0:
                tree.GetEntry(i)
                yield tree
                
    def _load_data(self):
        data = []
        for fname in self._filenames:
            tfile = ROOT.TFile(fname, "READ")
            treename = self._treename
            tree = tfile.Get(treename)
            if not tree:
                raise Exception("Input file does not contain tree with name", treename, fname)
            prescale = self._prescale
            iterevents = progress.printProgress("Reading " + treename, tree.GetEntries()/prescale, self._iterprescale(tree, prescale))
            for event in iterevents:
                data.append(self._func(event))
        print len(data), data[0]
        data = numpy.array(data)
        return data

###############################################################################

class NeutFileReader(Reader):
    def __init__(self, filelist):
        super(NeutFileReader, self).__init__(filelist, "nRooTracker", NeutEventReader(), prescale=1)

###############################################################################

class FluxFileReader(Reader):
    def __init__(self, filelist):
        super(FluxFileReader, self).__init__(filelist, "flux", FluxEventReader(), prescale=20)
        

###############################################################################

class PlotInteractionRate:
    def __init__(self, outdir, indices, data, scale_to_m=1.0, namemod=None, numbins=25):
        self._scale_to_m = scale_to_m
        self._outdir = outdir
        self._data = data
        self._indices = indices
        self._namemod = namemod
        self._numbins = numbins
    
    def __iter__(self):
        if len(self._data):
            for p in self._plot_xy_hist():
                yield p
#             for p in self._plot_histograms():
#                 yield p
    
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
            if self._namemod is not None:
                name += "_" + self._namemod
            yield name, fig
            plt.close(fig)
    
    def _plot_xy_hist(self):
        fig = plt.figure()
        name = os.sep.join((self._outdir, "hist_2d", "xy"))
        if self._namemod is not None:
            name += "_" + self._namemod
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
        xmin = -4.0
        xmax = 4.0
        prange = ((xmin, xmax), (xmin, xmax))
        ret = ax.hist2d(dx, dy, bins=self._numbins, weights=weights, cmap=matplotlib.cm.get_cmap("hot"), range=prange)
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

class Flux:
    def __init__(self, data, weights):
        #bin the data
        nbins = 100.0
        enurange = (0.0, 5.0)
        hist, edge = numpy.histogram(data, bins=nbins, range=enurange, density=True, weights=weights)
        y = hist
        x = (edge[1:] + edge[:-1] ) / 2.0
        self.f = scipy.interpolate.interp1d(x, y)
    
    def __call__(self, x):
        return self.f(x)

###############################################################################

class Xsec:
    def __init__(self):
        pass
    
    def __call__(self, x):
        return 1.0

###############################################################################

class FluxIntXsec:
    def __init__(self, flux, xsec):
        enu = numpy.linspace(start=0.01, stop=5.0, num=1000)
        ff = numpy.vectorize(flux)
        fa = ff(enu)
        xf = numpy.vectorize(xsec)
        xa = xf(enu)
        avg = numpy.average(xa, weights=fa)
        self._result = avg
    
    def __call__(self):
        return self._result

###############################################################################

class Ratio(object):
    def __init__(self, numerator, denominator):
        self._n = numerator
        self._d = denominator
        try:
            self._r = float(self._n) / float(self._d)
        except ZeroDivisionError:
            self._r = 0.0
    
    def __div__(self, rhs):
        return Ratio(self, rhs)
    
    def __sub__(self, rhs):
        return Difference(self, rhs)
    
    def __float__(self):
        return self._r
        
    def __str__(self):
        #return "{0:.3e} +- {1:.3e}".format(self._r, self.error())
        return "{0:.3f} +- {1:.3f}".format(self._r, self.error())
#         return "\n".join(("Num="+str(self._n),
#                          "Den="+str(self._d),
#                          "R="+str(self._r),
#                          ))
        
    def error(self):
        err = 0.0
        try:
            n = float(self._n)
            d = float(self._d)
            ne = self._n.error()
            de = self._d.error()
            err = self._r * math.sqrt((ne/n)**2 + (de/d)**2) 
        except AttributeError:
            #assume binomial error
            eff = self._r
            err = 0.0
            if self._d > 0:
                var = eff * (1.0 - eff) / float(self._n)
            err = math.sqrt(var)
        return err
    
class Difference(Ratio):
    def __init__(self, numerator, denominator):
        self._n = numerator
        self._d = denominator
        self._r = float(self._n) - float(self._d)
        
    def error(self):
        err = 0.0
        n = float(self._n)
        d = float(self._d)
        ne = self._n.error()
        de = self._d.error()
        err = math.sqrt((ne)**2 + (de)**2)
        return err

###############################################################################


class PlotFluxIntXsec:
    def __init__(self, neutindices, neutdata, fluxindices, fluxdata, scale_to_m = 1.0/100.0):
        self._neut_indices = neutindices
        self._neut_data = neutdata
        self._fluxindices = fluxindices
        self._fluxdata = fluxdata
        self._scale_to_m = scale_to_m
    
    def __iter__(self):
        rad = 4.0
        outer_xmax = rad / math.sqrt(2.0)
        outer_ymax = outer_xmax
        outer_xmin = -outer_xmax
        outer_ymin = -outer_ymax
        #get fiducial flux
        fvflux = _spacecut(self._fluxindices, self._fluxdata, outer_xmin, outer_xmax, outer_ymin, outer_ymax, scale_to_m=self._scale_to_m)
        fvneut = _spacecut(self._neut_indices, self._neut_data, outer_xmin, outer_xmax, outer_ymin, outer_ymax)
        #divide square into quadrants
        rad = 4.0
        outer_xmax = rad / math.sqrt(2.0)
        outer_ymax = outer_xmax
        outer_xmin = -outer_xmax
        outer_ymin = -outer_ymax
        for name, xmin, xmax, ymin, ymax in [("(+x, +y)", 0.0, outer_xmax, 0.0, outer_ymax),
                                             ("(+x, -y)", 0.0, outer_xmax, outer_ymin, 0.0),
                                             ("(-x, +y)", outer_xmin, 0.0, 0.0, outer_ymax),
                                             ("(-x, -y)",outer_xmin, 0.0, outer_ymin, 0.0),
                                       
                                       
                                       ]:
            #print "Plotting ", name
            for p in self._plot(xmin, xmax, ymin, ymax, fvflux, fvneut, self._fluxindices, self._neut_indices):
                yield p
        return
        
    def _plot(self, xmin, xmax, ymin, ymax, fluxdata, neut_data, fluxindices, neut_indices):
        fd = _spacecut(fluxindices, fluxdata, xmin, xmax, ymin, ymax, scale_to_m=self._scale_to_m)
        gd = _spacecut(neut_indices, neut_data, xmin, xmax, ymin, ymax)
#         #calc flux integrated xsec
#         flux = Flux(fd[:, Record.enu], fd[:, Record.weight])
#         xsec = Xsec()
#         fluxintxsec = FluxIntXsec(flux, xsec)
        #yield figname, fig
        expected_frac = Ratio(self._sumflux(fluxindices, fd), self._sumflux(fluxindices, fluxdata))
        observed_frac = Ratio(len(gd), len(neut_data))
        diff = observed_frac - expected_frac
        #print "expected=", str(expected_frac), ", observed=", str(observed_frac),", diff=",(diff)
        l = [(100.0*float(v), 100.0*v.error()) for v in (expected_frac, observed_frac, diff)]
        for t in l: 
            print t,","
        yield None
        
    def _sumflux(self, indices, fd):
        return numpy.sum(fd[:, indices[Record.weight]])
        
def _spacecut(indices, fd, xmin, xmax, ymin, ymax, scale_to_m=1.0):
        xmin, xmax, ymin, ymax = [[v / scale_to_m] for v in [xmin, xmax, ymin, ymax]]
        #x-cut
        ix = indices[Record.x]
        fd = fd[fd[:, ix] > xmin]
        fd = fd[fd[:, ix] < xmax]
        #y cut
        iy = indices[Record.y]
        fd = fd[fd[:, iy] > ymin]
        fd = fd[fd[:, iy] < ymax]
        return fd
            

###############################################################################

def plot_interaction_rate(args):
    parser = FileNameParser(args.patterns)
    f = parser()
    reader = NeutFileReader(args.patterns)
    indices, data = reader.read()
    out = simplot.io.FigureWriter()
    name = "_".join([f.polarity, f.ndid])
    plotter = PlotInteractionRate("evgen", indices, data,
                                  namemod=name,
                                  numbins=10,
                                  )
    out(plotter)
    return

###############################################################################

def plot_flux(args):
    ntuple = nuOscillation.model.beam.FluxNtuple(runtime.getcontext().beamcontext)
    reader = FluxFileReader([ntuple.filename()])
    indices, data = reader.read()
    DetectorId = nuOscillation.model.beam.constants.DetectorId
    out = simplot.io.FigureWriter()
    for polarity in nuOscillation.model.constants.Polarity.ALL:
        for ndid in DetectorId.ALL:
            #filter out particular near detector
            indid = DetectorId.toint(ndid)
            fluxdata = data[data[:, indices[Record.ndid]] == indid]
            if len(fluxdata) > 0:
                if polarity > 0:
                    fluxdata = fluxdata[fluxdata[:, indices[Record.pdg]] > 0]
                else:
                    fluxdata = fluxdata[fluxdata[:, indices[Record.pdg]] < 0]
            detstr = DetectorId.tostring(ndid)
            cm_to_m = 1.0 / 100.0
            name = "_".join((detstr, str(polarity)))
            plotter = PlotInteractionRate("flux", indices, fluxdata, 
                                          scale_to_m=cm_to_m, 
                                          namemod=name,
                                          numbins=10,
                                          )
            out(plotter)

###############################################################################

def plot_flux_int_xsec(args):
    parser = FileNameParser(args.patterns)
    f = parser()
    #load neut data
    reader = NeutFileReader(args.patterns)
    neutindices, neutdata = reader.read()
    #load flux data
    ntuple = nuOscillation.model.beam.FluxNtuple(runtime.getcontext().beamcontext)
    reader = FluxFileReader([ntuple.filename()])
    fluxindices, fluxdata = reader.read()
    DetectorId = nuOscillation.model.beam.constants.DetectorId
    ndid = DetectorId.toint(f.ndid)
    #filter out particular near detector
    fluxdata = fluxdata[fluxdata[:, fluxindices[Record.ndid]] == ndid]
    #make plot
    out = simplot.io.FigureWriter()
    cm_to_m = 1.0 / 100.0
    #cut on flux energy
    plotter = PlotFluxIntXsec(neutindices, neutdata, fluxindices, fluxdata, 
                              scale_to_m=cm_to_m,
                              )
    #out(plotter)
    for p in plotter:
        pass 

def testcut(indices, d, scale_to_m):
        rad = 3.0
        outer_xmax = rad / math.sqrt(2.0)
        outer_ymax = outer_xmax
        outer_xmin = -outer_xmax
        outer_ymin = -outer_ymax
        #get fiducial flux
        return _spacecut(indices, d, outer_xmin, outer_xmax, outer_ymin, outer_ymax, scale_to_m=scale_to_m)

###############################################################################

def parsecml():
    parser = argparse.ArgumentParser(description="Verify genev.")
    parser.add_argument("patterns", metavar="M", type=str, nargs='+', help="Input file list.")
    #parser.add_argument("-f", "--flux", metavar="F", dest="fluxpatterns", type=str, nargs='+', help=".")
    parser.add_argument("-p", "--plot", metavar="P", dest="plot", type=str, choices=["flux", "rate", "ratio"], help="Choose which plot to make.")
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
    elif args.plot == "ratio":
        plot_flux_int_xsec(args)
    return

###############################################################################

if __name__ == "__main__":
    main()

###############################################################################
