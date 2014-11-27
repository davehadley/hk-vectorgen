import glob
import collections
import itertools
import ROOT

###############################################################################

class Context:
    def __init__(self, beamcontext=None):
        if beamcontext is None:
            beamcontext = BeamContext()
        self.beamcontext = beamcontext

###############################################################################

class BeamContext:
    #TODO : move detector id, plane definitions and flux ntuple into the beam context.
    def __init__(self, jnubeamfiles=None, beamntuplefiles=None, fluxplanes=None):
        if jnubeamfiles is None:
            jnubeamfiles = JnuBeamFiles()
        self._jnubeamfiles = jnubeamfiles
        if beamntuplefiles is None:
            beamntuplefiles = BeamNtupleFiles()
        self._beamntuplefiles = beamntuplefiles
        if fluxplanes is None:
            fluxplanes = FluxPlaneDefinitions()
        self._fluxplanes = fluxplanes
        #get flux for each detector and both nu/antinu running
    
    def jnubeamfiles(self):
        return self._jnubeamfiles
    
    def beam_ntuple_files(self):
        return self._beamntuplefiles
    
    def flux_planes(self):
        return self._fluxplanes

###############################################################################

class BeamNtupleFiles:

    _standard_flux_neardetector_ntuple_file = ""
    _standard_flux_superk_ntuple_file = ""

    def __init__(self, nd_ntuple_name=None, sk_ntuple_name=None):
        if nd_ntuple_name is None:
            nd_ntuple_name = BeamNtupleFiles._standard_flux_neardetector_ntuple_file
        if sk_ntuple_name is None:
            sk_ntuple_name = BeamNtupleFiles._standard_flux_superk_ntuple_file
                
        self.nd_ntuple_name = nd_ntuple_name
        self.far_ntuple_name = sk_ntuple_name 

###############################################################################

class JnuBeamFiles:
    
    _default_nu_flux_files = glob.glob("/home/dave/t2k/data/irods/QMULZone1/home/hyperk/gfluka_flux/flux_320kA/hk_flux_gfluka_*.root")
    _default_antinu_flux_files = glob.glob("/home/dave/t2k/data/irods/QMULZone1/home/hyperk/gfluka_flux/flux_m320kA/hk_flux_gfluka_*.root")
    
    def __init__(self, nu_flux_files=_default_nu_flux_files, antinu_flux_files=_default_antinu_flux_files):
        self.nu_flux_files = list(sorted(nu_flux_files))
        self.antinu_flux_files = list(sorted(antinu_flux_files))
        
    def verify(self):
        #all files should contain the tree "h2000"
        for fname in itertools.chain(self.nu_flux_files, self.antinu_flux_files):
            tfile = ROOT.TFile(fname, "READ")
            if not tfile.IsOpen():
                raise Exception("JnuBeamFiles failed to open file", fname)
            treename = "h2000"
            if not tfile.Get("h2000"):
                raise Exception("JnuBeamFiles missing treename", treename, fname)
        return

###############################################################################

class FluxPlaneDefinitions:
    def __init__(self, fluxplanes=None):
        self._fluxplanes = dict()
        if fluxplanes is not None:
            for p in fluxplanes:
                self.add(p)
    
    def add(self, plane):
        self._fluxplanes[plane.name] = plane
        print "DEBUG adding ", plane.name, plane, self._fluxplanes
        return
    
    def tostring(self, ndid):
        return self._fluxplanes[ndid].name

    def baseline(self, ndid):
        return self._fluxplanes[ndid].baseline
    
    def flukaid(self, ndid):
        return self._fluxplanes[ndid].flukaid

###############################################################################

FluxPlane = collections.namedtuple("FluxPlane", ["name", "baseline", "flukaid"])

###############################################################################

#provide a global instance of the context
_current_context = Context()

def getcontext():
    return _current_context

def setcontext(context):
    global _current_context
    _current_context = context
    
###############################################################################
