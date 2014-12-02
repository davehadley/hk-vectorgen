
import glob
import math
import os
import ROOT
import sys

from argparse import ArgumentParser
from collections import defaultdict

###############################################################################

_INDENTSTR = "--"

_pdgnames = {
    22 : "gamma",
    -13 : "mu+",
     13 : "mu-",
    -14 : "numubar",
     14 : "numu",
    -11 : "e+",
     11 : "e-",
    -12 : "nuebar",
     12 : "nue",
    2212 : "p",
    2112 : "n",
     111 : "pi0",
     211 : "pi+",
    -211 : "pi-",
}

###############################################################################

def pdgstr(pdg):
    try:
        return _pdgnames[pdg]
    except KeyError:
        return str(pdg)

###############################################################################

def dumptracks(tracks, prefix="", iostream=sys.stdout):
    prefix += _INDENTSTR
    for tr in tracks:
        ipnu = pdgstr(tr.GetIpnu())
        flag = tr.GetFlag()
        m = tr.GetM()
        p = tr.GetP()
        e = tr.GetE()
        dx = tr.GetDir(0)
        dy = tr.GetDir(1)
        dz = tr.GetDir(2)
        px = tr.GetPdir(0)
        py = tr.GetPdir(1)
        pz = tr.GetPdir(2)
        bx = tr.GetStart(0)
        by = tr.GetStart(1)
        bz = tr.GetStart(2)
        ex = tr.GetStop(0)
        ey = tr.GetStop(1)
        ez = tr.GetStop(2)
        ptype = tr.GetParenttype()
        row = "id={:>5.0f}, pdg={:>7}, parent={:>5.0f}, flag={:>2.0f}, m={:>5.0f}, p={:>5.0f}, e={:>5.0f}, dir=({:>5.2f}, {:>5.2f}, {:>5.2f}), start=({:>5.0f}, {:>5.0f}, {:>5.0f}), stop=({:>5.0f}, {:>5.0f}, {:>5.0f})".format(tr.GetId(), ipnu, ptype, flag, m, p, e, dx, dy, dz, bx, by, bz, ex, ey, ez)
        print >>iostream, prefix, row
    return

###############################################################################

def dumpevent(event, prefix="", iostream=sys.stdout):
    prefix += _INDENTSTR
    for itrig in xrange(event.GetNumberOfEvents()):
        print >>iostream, prefix, "Trigger", itrig
        trig = event.GetTrigger(itrig)
        tracks = trig.GetTracks()
        dumptracks(tracks, prefix=prefix, iostream=iostream)
    return

###############################################################################

def dumptree(tree, iostream=sys.stdout, selectevent=None, nevents=None):
    prefix = _INDENTSTR
    count = 0
    for ievent, event in enumerate(tree):
        if nevents is not None and count >= nevents:
            break
        if selectevent is None or selectevent(event.wcsimrootevent):
            print prefix, "Event", ievent
            dumpevent(event.wcsimrootevent, prefix=prefix, iostream=iostream)
            count += 1
    return

###############################################################################

def dumpfiles(args, nevents):
    for tree in iterwcsimtrees(args):
        dumptree(tree, nevents=nevents)
        break
    return

###############################################################################

def loadlib():
    wcsimdir = os.environ["WCSIMDIR"]
    lib = os.path.sep.join([wcsimdir, "libWCSimRoot"])
    ROOT.gSystem.Load(lib)
    return

###############################################################################

def iterinputfiles(args):
    for pattern in args.inputfiles:
        for fname in glob.glob(pattern):
            yield fname
    return

###############################################################################

def iterwcsimtrees(args):
    for fname in iterinputfiles(args):
        tfile = ROOT.TFile(fname)
        tree = tfile.Get("wcsimT")
        tree.GetBranch("wcsimrootevent").SetAutoDelete(ROOT.kTRUE) # IMPORTANT: this prevents memory leak somehow
        yield tree
    return

###############################################################################

def iterwcsimevents(args):
    for fname in iterinputfiles(args):
        tfile = ROOT.TFile(fname)
        tree = tfile.Get("wcsimT")
        tree.GetBranch("wcsimrootevent").SetAutoDelete(ROOT.kTRUE) # IMPORTANT: this prevents memory leak somehow
        for entry in tree:
            yield entry
    return

###############################################################################

class PlotVertex:
    def __init__(self):
        self._cm_to_m = 1.0e-2
        self._x = "x"
        self._y = "y"
        self._z = "z"
        self._dim = ["x", "y", "z"]
        self._dimnum = {self._x : 0,
                        self._y : 1,
                        self._z : 2,

        }

        self._proj2d = [(self._x, self._y), (self._x, self._z), (self._y, self._z)]
        self._binning = { self._x : (100, -10., 10.),
                          self._y : (100, -10., 10.),
                          self._z : (100, -15., 15.),
        }
        #create 1D histograms
        self._hist1d = dict()
        for d in self._dim:
            nbins, xlow, xhigh = self._binning[d]
            h = ROOT.TH1F("vertex_1d_" + d, "vertex " + d + " position;" + d + " [m]", nbins, xlow, xhigh)
            self._hist1d[d] = h
        #create 2D histograms
        self._hist2d = dict()
        for dx, dy in self._proj2d:
            xbins, xlow, xhigh = self._binning[dx]
            ybins, ylow, yhigh = self._binning[dy]
            name = "vertex_2d_" + dx + "_" + dy
            title = "vertex " + dx + "-" + dy + " position; " + dx + " [m]; " + dy + " [m]"
            h = ROOT.TH2F(name, title, xbins, xlow, xhigh, ybins, ylow, yhigh)
            self._hist2d[(dx, dy)] = h

    def add(self, event):
        event = event.wcsimrootevent
        vtx = self._getvertex(event)
        self._fillvertex(vtx)
        return

    def _getvertex(self, event):
        metres = self._cm_to_m
        #get first trigger
        for itrig in xrange(event.GetNumberOfEvents()):
            trig = event.GetTrigger(itrig)
            x = trig.GetVtx(0) * metres
            y = trig.GetVtx(1) * metres
            z = trig.GetVtx(2) * metres
            return (x, y, z)
        raise Exception("no vertex in event", event)

    def _fillvertex(self, vtx):
        self._hist1d[self._x].Fill(vtx[0])
        self._hist1d[self._y].Fill(vtx[1])
        self._hist1d[self._z].Fill(vtx[2])
        for xv, yv in self._proj2d:
            xi = self._dimnum[xv]
            yi = self._dimnum[yv]
            self._hist2d[(xv, yv)].Fill(vtx[xi], vtx[yi])
        return

    def draw(self):
        for h in self._hist1d.itervalues():
            canv = ROOT.TCanvas(h.GetName(), h.GetTitle(), 800, 600)
            h.Draw()
            yield canv
        for h in self._hist2d.itervalues():
            canv = ROOT.TCanvas(h.GetName(), h.GetTitle(), 800, 600)
            h.Draw()
            yield canv
        return

###############################################################################

class PlotMuonDirection:
    def __init__(self):
        self._xmin, self._xmax = -1., 1.
        self._nupdg = (12, -12, 14, -14)
        self._leppdg = (11, -11, 13, -13)
        self._hist_mu = ROOT.TH1F("#theta_#mu", "#theta_#mu;#theta_#mu", 100, self._xmin, self._xmax)
        self._hist_nu = ROOT.TH1F("#theta_#nu", "#theta_#nu;#theta_#nu", 100, self._xmin, self._xmax)

    def add(self, event):
        trig = event.wcsimrootevent.GetTrigger(0)
        #get neutrino vector
        nuvec = None
        for tr in trig.GetTracks():
            if tr.GetIpnu() in self._nupdg:
                nuvec = tr
                break
        #get charged lepton vector
        lepvec = None
        for tr in trig.GetTracks():
            if tr.GetIpnu() in self._leppdg:
                lepvec = tr
                break
        if nuvec:
            theta = self._theta(nuvec)
            self._hist_nu.Fill(theta)
        if lepvec:
            theta = self._theta(lepvec)
            self._hist_mu.Fill(theta)
        return

    def _theta(self, track):
        x = track.GetDir(0)
        y = track.GetDir(1)
        z = track.GetDir(2)
        return ROOT.TVector3(x, y, z).Unit().CosTheta()

    def draw(self):
        allhist = [(self._hist_nu, ROOT.kRed), (self._hist_mu, ROOT.kBlue)]
        canv = ROOT.TCanvas("lepton_direction", "lepton_direction", 800, 600)
        xmin = -math.pi
        xmax = math.pi
        ymin = 0.1
        ymax = max([1.1*h.GetMaximum() for h, c in allhist])
        frame = canv.DrawFrame(xmin, ymin, xmax, ymax, ";#theta [radians];events")
        canv.SetLogy()
        leg = ROOT.TLegend(0.8, 0.8, 0.99, 0.99)
        for hist, col in allhist:
            hist.SetLineColor(col)
            hist.SetMarkerColor(col)
            hist.Draw("SAME")
            leg.AddEntry(hist, hist.GetTitle())
        leg.Draw("LP")
        yield canv

###############################################################################

class PlotCapturedNuclei:
    def __init__(self):
        self._counter = defaultdict(int)
        self._neutronpdg = 2112
        self._nucleusname = {1: "H",
                             64 : "Gd",
                             }

    def _getatomname(self, Z):
        try:
            return self._nucleusname[Z]
        except KeyError:
            return "Other"

    def add(self, event):
        track = self._findneutroncapturenucleus(event)
        if track:
            (Z, A) = self._decodenucleuspdg(track.GetIpnu())
            name = self._getatomname(Z)
            self._counter[name] += 1
        return

    def _decodenucleuspdg(self, pdgcode):
        #01 2 345 678 9
        #10 L ZZZ AAA I
        #where A = np + nn
        #L = number of lambdas
        #I = isomer level (>0 is excited states)
        spdg = str(pdgcode)
        L = spdg[3]
        Z = int(spdg[3:6])
        A = int(spdg[6:9])
        I = int(spdg[9])
        return (Z, A)

    def _findneutroncapturenucleus(self, event):
        result = None
        event = event.wcsimrootevent
        for itrig in xrange(event.GetNumberOfEvents()):
            trig = event.GetTrigger(itrig)
            tracks = trig.GetTracks()
            for tr1 in tracks:
                pdg = tr1.GetIpnu()
                if pdg > 1000000:
                    #is nucleus
                    if self._parentpdg(tr1, tracks) == self._neutronpdg:
                        #parent is neutron
                        result = tr1
        return result

    def _parentpdg(self, tr, tracks):
        result = None
        pt = tr.GetParenttype()
        for tr2 in tracks:
            if tr2.GetId() == pt:
                result = tr2.GetIpnu()
        return result

    def draw(self):
        for c in self._drawcounter():
            yield c

    def _drawcounter(self):
        nbins = len(self._counter)
        hist = ROOT.TH1F("ncapture_count", "ncapture_count;capturing nucleus;fraction", nbins, 0, nbins)
        total = sum(self._counter.itervalues())
        for i, (n, c) in enumerate(sorted(self._counter.iteritems())):
            hist.SetBinContent(i + 1, float(c) / float(total))
            hist.GetXaxis().SetBinLabel(i + 1, n)
        canv = ROOT.TCanvas("ncapture_count", "neutron capture nuclei", )
        hist.Draw()
        yield canv

###############################################################################

class PlotPi0:
    def __init__(self):
        self._pi0pdg = 111
        self._gammapdg = 22
        pi0mass = 134.9766
        deltam = 0.001
        self._pi0mass = ROOT.TH1F("pi0mass", "pi0mass;#pi^{0} mass [MeV]", 100, pi0mass - deltam, pi0mass + deltam)
        self._pi0mass.SetBit(ROOT.TH1.kCanRebin)
        ebins = 100
        emin = 0.0
        emax = 1000.0
        self._pi0_g1_e = ROOT.TH1F("pi0_g1_e", "pi0_g1_e;leading photon energy [MeV]", ebins, emin, emax)
        self._pi0_g2_e = ROOT.TH1F("pi0_g2_e", "pi0_g2_e;second photon energy [MeV]", ebins, emin, emax)
        abins = 100
        amin = 0.0 #-math.pi
        amax = math.pi
        self._pi0_g1_theta = ROOT.TH1F("pi0_g1_theta", "pi0_g1_theta;leading photon angle [radians]", abins, amin, amax)
        self._pi0_g2_theta = ROOT.TH1F("pi0_g2_theta", "pi0_g2_theta;second photon angle [radians]", abins, amin, amax)
        self._simplehist = [self._pi0mass, self._pi0_g2_e, self._pi0_g1_e, self._pi0_g2_theta, self._pi0_g1_theta]

    def add(self, event):
        for (pi0, g1, g2) in self._getpi0s(event):
            self._fillpi0(pi0, g1, g2)
        return

    def _fillpi0(self, pi0, g1, g2):
        pg1 = self._p4(g1)
        pg2 = self._p4(g2)
        #fill mass plot
        mass = (pg1 + pg2).M()
        self._pi0mass.Fill(mass)
        #fill photon kinematics
        self._pi0_g1_e.Fill(pg1.E())
        self._pi0_g2_e.Fill(pg2.E())
        self._pi0_g1_theta.Fill(pg1.Theta())
        self._pi0_g2_theta.Fill(pg2.Theta())
        return

    def _p4(self, track):
        p = track.GetP()
        x = track.GetDir(0)
        y = track.GetDir(1)
        z = track.GetDir(2)
        hlv = ROOT.TLorentzVector(p*x, p*y, p*z, track.GetE())
        return hlv

    def _getpi0s(self, event):
        result =  []
        for itrig in xrange(event.wcsimrootevent.GetNumberOfEvents()):
            trig = event.wcsimrootevent.GetTrigger(itrig)
            tracks = trig.GetTracks()
            #find pi0's
            pi0tracks = [tr for tr in tracks if tr.GetIpnu() == self._pi0pdg]
            for pi0 in pi0tracks:
                #get children (only consider pi0 -> gamma gamma channel)
                children = []
                for tr in tracks:
                    if tr.GetParenttype() == pi0.GetId() and tr.GetIpnu() == self._gammapdg:
                        children.append(tr)
                if len(children) == 2:
                    g1, g2 = children
                    if g2.GetE() > g1.GetE():
                        #swap photons
                        g1, g2 = g2, g1
                    result.append( (pi0, g1, g2) )
        return result

    def draw(self):
        for h in self._simplehist:
            canv = ROOT.TCanvas(h.GetName(), h.GetName(), 800, 600)
            h.Draw()
            yield canv


###############################################################################

def makeplots(args):
    ROOT.gROOT.SetBatch(True)
    #define plots
    plots = [PlotVertex(),
             PlotMuonDirection(),
             PlotCapturedNuclei(),
             PlotPi0(),
             ]
    #fill plots
    for event in iterwcsimevents(args):
        for p in plots:
            p.add(event)
    #draw and save plots
    for p in plots:
        for canv in p.draw():
            for ext in ["eps", "png"]:
                fname = os.sep.join(["plots", ext, canv.GetName() + "." + ext])
                if not os.path.exists(os.path.dirname(fname)):
                    os.makedirs(os.path.dirname(fname))
                canv.SaveAs(fname)
    return

###############################################################################

def parsecml():
    parser = ArgumentParser()
    parser.add_argument("inputfiles", help="Input files.", type=str, nargs="+")
    parser.add_argument("-d", "--dump", dest="dump", help="Dump event information.", type=int, default=0)
    parser.add_argument("-p", "--plot", dest="plot", help=".", action="store_true", default=False)
    return parser.parse_args()

###############################################################################

def main():
    loadlib()
    args = parsecml()
    if args.dump:
        dumpfiles(args, args.dump)
    if args.plot:
        makeplots(args)
    return

###############################################################################

if __name__ == "__main__":
    main()
