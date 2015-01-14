#!/bin/env python
import ROOT  # @UnresolvedImport
ROOT.PyConfig.IgnoreCommandLineOptions = True
import pkg_resources
import argparse
import sys
import math

###############################################################################

_pdg_database = ROOT.TDatabasePDG()

###############################################################################

def pdg_to_string(pdg):
    try:
        part = _pdg_database.GetParticle(pdg)
        result = part.GetName()
    except:
        result = str(pdg)
    return result

###############################################################################

def dumpevent(filename, eventnum=None, wait=True):
    #open file and get tree
    tfile = ROOT.TFile(filename)
    tree = tfile.Get("nRooTracker"); # try to get NEUT tree
    if not tree:
        tree = tfile.Get("gRooTracker"); # try to get GENIE tree
    #Loop over events
    nevents = tree.GetEntries();
    if eventnum is not None:
        event_range = xrange(eventnum, eventnum + 1)
    else:
        event_range = xrange(nevents)
    for jj in event_range:
        _dump_single_event(tree, jj)
        if wait:
            raw_input("press enter to continue")
    return

###############################################################################

def _dump_single_event(tree, jj, ostr=sys.stdout):
        tree.GetEntry(jj)
        m_to_cm = 100.0
        GeV_to_MeV = 1000.0
        try:
            mode = tree.NEneutmode
        except AttributeError:
            mode = tree.G2NeutEvtCode
        vtx = tree.EvtVtx
    
        x = vtx[0] * m_to_cm
        y = vtx[1] * m_to_cm
        z = vtx[2] * m_to_cm
        t = vtx[3] * m_to_cm
        r = math.sqrt(x**2 + y**2)

        print >>ostr, "    Event number", jj
        print >>ostr, "Interaction code", mode
        print >>ostr, "          vertex x=%5.0fcm, y=%5.0fcm, r=%5.0fcm, z=%5.0fcm, t=%5.2f" % (x, y, r, z, t)
        #Loop over particles
        npart = tree.StdHepN
        stdhepp4 = ROOT.getRooTrackerHepP4(tree.StdHepN, tree.StdHepP4)
        header = ["index", "pdg", "particle", "status", "E [Mev]", "p [MeV]", "m [MeV]", "status name"]
        headerfmt = " | ".join(["{:<10s}"] * len(header))
        ifmt = "{:>10.0f}"
        ffmt = "{:>10.0f}"
        sfmt = "{:>10s}"
        rowfmt = " | ".join([ifmt, ifmt, sfmt, ifmt, ffmt, ffmt, ffmt, sfmt])
        header = headerfmt.format(*header)
        hline = "-"*len(header)
        print >>ostr,header
        print >>ostr, hline
        for i in xrange(npart):
            #get the information that we need
            pvec = stdhepp4.at(i)
            pid = tree.StdHepPdg[i]
            stdhepstatus = tree.StdHepStatus[i]
            #p = ROOT.TMath.Sqrt(
            #                    pvec.E()*(pvec.E()-part.fMass*part.fMass)
            #                    )
            pvec = ROOT.TLorentzVector(pvec[0] * GeV_to_MeV, 
                                       pvec[1] * GeV_to_MeV, 
                                       pvec[2] * GeV_to_MeV, 
                                       pvec[3] * GeV_to_MeV
                                       )
            if stdhepstatus==0:
                #initial state particle
                status_string = "initial";
            elif stdhepstatus == 1:
                #final state particle
                status_string = "final";
            else:
                #unknown (intermediate?)
                status_string = "inter";
            momentum = pvec.P()
            energy = pvec.E()
            mass = pvec.M()
            #print >>ostr, "particle pdg=%s, status=%s, E=%.0f, p=%.0f, m=%.0f (%s)" % (pid, stdhepstatus, energy, momentum, mass, status_string)
            row = [i, pid, pdg_to_string(pid), stdhepstatus, energy, momentum, mass, status_string]
            print >>ostr, rowfmt.format(*row)
        print >>ostr, hline
        return

###############################################################################

def loadlib():
    srcfile = pkg_resources.resource_filename("vectorgen", "readRooTracker.C+")
    ROOT.gROOT.ProcessLine(".L " + srcfile)
    return

###############################################################################

def parsecml():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=str, help="Input genev file name.")
    parser.add_argument("-n", "--eventnum", type=int, help="Select event number.")
    return parser.parse_args()

###############################################################################

def main():
    loadlib()
    opt = parsecml()
    dumpevent(opt.infile, eventnum=opt.eventnum)
    return

###############################################################################

if __name__ == "__main__":
    main()

