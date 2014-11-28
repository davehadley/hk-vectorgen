import ROOT
import shutil
import tempfile
ROOT.PyConfig.IgnoreCommandLineOptions = True

import subprocess
import os
import glob
import argparse
import itertools
import collections

import runtime

#Jobs:
#    (1) Merge flux files.
#    (2) Create geometry.
#    (3) Run event_rate.
#    (4) Run genev.

###############################################################################

def _abspath(path):
    f = os.path.expandvars(os.path.expanduser(path))
    return os.path.abspath(f)

###############################################################################

class IJob(object):
    def __init__(self, rundir=None, test=False):
        self._test = test
        if rundir is None:
            rundir = RunDir()
        self._rundir = rundir
        
    def _tmp_chdir(self, workingdir, func, *args, **kwargs):
        #keep original directory to change back to
        origdir = os.getcwd()
        #change directory
        if workingdir is not None:
            os.chdir(workingdir)
        #actually run function
        try:
            ret = func(*args, **kwargs)
        finally:
            #change back to the original directory
            if workingdir is not None:
                os.chdir(origdir)
        return ret
    
    def _check_call(self, cmd, workingdir=None):
        if workingdir is None:
            workingdir = self._rundir.rundir()
        ret = None
        if self._test:
            print "[TEST]",cmd
        else:
            ret = self._tmp_chdir(workingdir, subprocess.check_call, cmd, shell=True)
        return ret
    
    def verify(self):
        return
    
    def run(self):
        return

###############################################################################

BeamPlane = collections.namedtuple("BeamPlane", ["name", "baseline", "ndcode"])

###############################################################################

def plane_from_ndid(ndid, context):
    name = context.beamcontext.flux_planes().tostring(ndid)
    baseline = context.beamcontext.flux_planes().baseline(ndid)
    code = context.beamcontext.flux_planes().flukaid(ndid)
    plane = BeamPlane(name=name, baseline=baseline, ndcode=code)
    return plane

###############################################################################
        
class RunDir:
    def __init__(self, path=None, card=None):
        if path is None:
            path = tempfile.mkdtemp(prefix="tmp_run_genev_")
        self._rundir = _abspath(path)
        self._card = self._find_card(card)
        self._run_make_links()
        
    def rundir(self):
        return self._rundir
    
    def _find_card(self, card):
        if card is None:
            #use default from NEUTGEOM directory
            card = "".join((os.environ["NEUT_ROOT"], os.sep, "src/neutgeom/neut.card"))
        card = _abspath(card)
        if not os.path.exists(card):
            raise Exception("Cannot find card file", card)
        return card
    
    def _run_make_links(self):
        outdir = self.rundir()
        try:
            os.makedirs(outdir)
        except os.error:
            #ignore as this happens if directory already exists
            pass
        inputdir = "".join((os.environ["NEUT_ROOT"], os.sep, "src", os.sep, "neutsmpl"))
        for fname in os.listdir(inputdir):
            src = "".join((inputdir, os.sep, fname))
            if os.path.islink(src):
                dst = "".join((self.rundir(), os.sep, fname))
                if not os.path.exists(dst):
                    os.symlink(src, dst)
        src = self._card
        dst = "".join((self.rundir(), os.sep, "neut.card"))
        shutil.copyfile(src, dst)
        return

###############################################################################

class MergeFluxJob(IJob):
    def __init__(self, beam_input, rundir=None, test=False):
        super(MergeFluxJob, self).__init__(rundir, test)
        self._beam_input = beam_input
        
    def run(self):
        if not os.path.exists(self._beam_input.filename()):
            self.verify()
            self._run_hadd()
        return
            
    def verify(self):
        self._beam_input.verify()
    
    def _run_hadd(self):
        cmd = " ".join(("hadd",
                       self._beam_input.filename(),
                       " ".join(self._beam_input.filelist()
                        )),
               )
        self._check_call(cmd)

###############################################################################
        
class MakeFluxLinks(IJob):
    def __init__(self, beam_input, rundir=None, test=False, n=None):
        super(MakeFluxLinks, self).__init__(rundir, test)
        self._n = n
        self._beam_input = beam_input
        
    def run(self):
        if not os.path.exists(self._beam_input.filename()):
            self.verify()
            self._run_make_links()
        return
            
    def verify(self):
        self._beam_input.verify()
    
    def _run_make_links(self):
        outdir = "".join((self._rundir.rundir(), os.sep, self._beam_input.linkdir()))
        try:
            os.makedirs(outdir)
        except os.error:
            #ignore as this happens if directory already exists
            pass
        for i, fname in enumerate(self._beam_input.filelist()):
            if self._n is not None and i >= self._n:
                break
            src = fname
            dst = "".join((self._rundir.rundir(), os.sep, self._beam_input.filestem(), ".", str(i), ".root"))
            if not os.path.exists(dst):
                os.symlink(src, dst)
        return

###############################################################################

class Orientation:
    Z = "Z"
    Y = "Y"

###############################################################################

class CylinderGeometry:
    def __init__(self, ndid, radius=4.0, z=8.0, orientation=Orientation.Z, context=None):
        if context is None:
            context = runtime.getcontext()
        self._context = context
        self.ndid = ndid
        self.radius = radius
        self.z = z
        self.orientation = orientation
        self.name = self._uniquestr()
        self._plane = plane_from_ndid(self.ndid, self._context)
        
    def verify(self):
        return
    
    def _uniquestr(self):
        return "_".join((self._context.beamcontext.flux_planes().tostring(self.ndid),
                         "cylinder",
                         self._float_to_string(self.radius, "r"),
                         self._float_to_string(self.z, "z"),
                         self.orientation,
                         ))
        
    def _float_to_string(self, f, prefix):
        return prefix + str(int(round(f * 100.0)))
    
    def filename(self):
        return self.name + ".root"
    
    def volume_name(self):
        return "wc_volume"
    
    def plane(self):
        return self._plane
    
    def build_detector_volume(self):
        #get dimensions
        m_to_mm = 1000.0
        radius_mm = self.radius * m_to_mm
        z_mm = self.z * m_to_mm / 2.0
        #build volume
        vol0 = ROOT.TGeoTube(0, radius_mm, z_mm);
        vol1 = ROOT.TGeoTube(0, radius_mm, z_mm);
        return (vol0, vol1)

###############################################################################

class CuboidGeometry:
    def __init__(self, ndid, radius=4.0, z=8.0, orientation=Orientation.Z, context=None):
        if context is None:
            context = runtime.getcontext()
        self._context = context
        self.ndid = ndid
        self.radius = radius
        self.z = z
        self.orientation = orientation
        self.name = self._uniquestr()
        self._plane = plane_from_ndid(self.ndid, self._context)
        
    def verify(self):
        return
    
    def _uniquestr(self):
        return "_".join((self._context.beamcontext.flux_planes().tostring(self.ndid),
                         "cuboid",
                         self._float_to_string(self.radius, "x"),
                         self._float_to_string(self.z, "z"),
                         self.orientation,
                         ))
        
    def _float_to_string(self, f, prefix):
        return prefix + str(int(round(f * 100.0)))
    
    def filename(self):
        return self.name + ".root"
    
    def volume_name(self):
        return "wc_volume"
    
    def plane(self):
        return self._plane
    
    def build_detector_volume(self):
        #get dimensions
        m_to_mm = 1000.0
        radius_mm = self.radius * m_to_mm
        z_mm = self.z * m_to_mm / 2.0
        x_mm = radius_mm
        y_mm = radius_mm
        #build volume
        vol0 = ROOT.TGeoBBox(x_mm, y_mm, z_mm);
        vol1 = ROOT.TGeoBBox(x_mm, y_mm, z_mm);
        return (vol0, vol1)

###############################################################################

class CreateGeometryJob(IJob):
    def __init__(self, geometry, rundir=None, test=False):
        super(CreateGeometryJob, self).__init__(rundir, test)
        self._geometry = geometry
        
    def run(self):
        if not os.path.exists(self._geometry.filename()):
            self.verify()
            self._run_geometry()
        return self._geometry.filename()
            
    def verify(self):
        self._geometry.verify()
        return
    
    def _run_geometry(self):
        g = self._geometry
        gen = GenWCGeom()
        gen(g)
        return

###############################################################################

class GenWCGeom:
    #according to the manual ROOT geometry distance units are in cm.
    #however looking at the existing nd280geometry.root, it appears to be in mm.    
    def __call__(self, geometry):
        g = geometry
        orientation = g.orientation
        if not orientation == Orientation.Z:
            raise Exception("not implemented") 
        outfilename = g.filename()
        volume_name = g.volume_name()
        #create geometry
        wc_geometry = ROOT.TGeoManager("ND280Geometry","ND280Geometry");
        oxygen = ROOT.TGeoElement("oxygen", "oxygen", 8, 16);
        hydrogen = ROOT.TGeoElement("hydrogen", "hydrogen", 1, 1);
        water = ROOT.TGeoMixture("water", 2, 1);
        water.AddElement(oxygen, 1); 
        water.AddElement(hydrogen, 2); 
        water_med = ROOT.TGeoMedium("water", 1, water);
        vol0, vol1 = g.build_detector_volume()
        t2k = ROOT.TGeoVolume("t2k",vol1);
        wc_volume = ROOT.TGeoVolume(volume_name, vol0, water_med);
        t2k.AddNode(wc_volume, 1);
        wc_geometry.AddVolume(t2k);
        wc_geometry.SetTopVolume(t2k);
        #write geometry to file
        outfile = ROOT.TFile(outfilename, "RECREATE");
        wc_geometry.Write()
#         wc_geometry.Export(outfilename.replace(".root", ".gdml"))
        outfile.Close()
        return
    
###############################################################################

class EventRateJob(IJob):
    def __init__(self, beam_input, geometry, rundir=None, test=False):
        super(EventRateJob, self).__init__(rundir, test)
        self._beam_input = beam_input
        self._geometry = geometry
        
    def filename(self):
        beamname = self._beam_input.name
        geomname = self._geometry.name
        outfilename = "_".join(("eventrate",
                               beamname,
                               geomname,
                               )) + ".root"
        return outfilename
    
    def run(self):
        outfilename = self.filename()
        if not os.path.exists(outfilename):
            self._create_event_rate()
        return
    
    def _create_event_rate(self):
        outfilename = _abspath(self.filename())
        beamfile = _abspath(self._beam_input.filename())
        filestem = self._beam_input.filestem()
        filestem = "".join((self._rundir.rundir(), os.sep, filestem))
        #N = len(self._beam_input.filename()) - 1
        N = len(glob.glob(filestem + "*.root")) - 1
        if N <= 0:
            raise Exception("No flux files matching", filestem)
        geomfile = _abspath(self._geometry.filename())
        volumename = self._geometry.volume_name()
        plane = self._geometry.plane()
        planenum = plane.ndcode
        neutgeompath = os.environ["NEUTGEOM"]
        cmd = " ".join((
                        os.sep.join((neutgeompath, "event_rate")),
                        #"-f", beamfile,
                        "-s", filestem, "0", str(N),
                        "-g",
                        geomfile,
                        "-v",
                        "+" + volumename,
                        "-o",
                        outfilename,
                        "-d",
                        str(planenum),
                        ))
        #setupneutcmd = "source /home/software/neut/setupNeut.sh" # TODO : move this to constants somewhere.
        #cmd = " && ".join((setupneutcmd, cmd))
        self._check_call(cmd)
        return
    
    def verify(self):
        pass

###############################################################################

class GenEvConfig:
    def __init__(self, num_events, nu_pdg):
        self.num_events = num_events
        if nu_pdg is None:
            nu_pdg = 0
        self.nu_pdg = nu_pdg
        self._nu_pdg_names = {0 : "allnuflav",
                              12 : "nue",
                              14 : "numu",
                              -12 : "antinue",
                              -14 : "antinumu",
                              }
        if (self.nu_pdg is not None) and (not self.nu_pdg in self._nu_pdg_names):
            raise Exception("unknown neutrino PDG", self.nu_pdg)
    
    @property
    def name(self):
        nupdgname = self._nu_pdg_names[self.nu_pdg]
        return "_".join((str(self.num_events),
                         nupdgname,
                  ))
    
###############################################################################

class GenEvJob(IJob):
    def __init__(self, gen_config, beam_input, geometry, eventratejob, rundir=None, test=False):
        super(GenEvJob, self).__init__(rundir, test)
        self._gen_config = gen_config
        self._beam_input = beam_input
        self._geometry = geometry
        self._eventrate = eventratejob
        
    def filename(self):
        beamname = self._beam_input.name
        geomname = self._geometry.name
        configname = self._gen_config.name
        outfilename = "_".join(("genev",
                               beamname,
                               geomname,
                               configname,
                               )) + ".root"
        return outfilename
    
    def run(self):
        outfilename = self.filename()
        if not os.path.exists(outfilename):
            self._create_genev()
        return

#genev  -j flux_files.root -g nd280geometry.root -v +Basket -o test.genev.output.2.root -n 10 -f rootracker -i setup_output.root -d 5 2>&1

    def _create_genev(self):
        outfilename = _abspath(self.filename())
        beamfile = _abspath(self._beam_input.filename())
        filestem = self._beam_input.filestem()
        filestem = "".join((self._rundir.rundir(), os.sep, filestem))
        #N = len(self._beam_input.filename()) - 1
        N = len(glob.glob(filestem + "*.root")) - 1
        if N <= 0:
            raise Exception("No flux files matching", filestem) 
        geomfile = _abspath(self._geometry.filename())
        eventratefile = _abspath(self._eventrate.filename())
        volumename = self._geometry.volume_name()
        plane = self._geometry.plane()
        planenum = plane.ndcode
        neutgeompath = os.environ["NEUTGEOM"]
        numevents = self._gen_config.num_events
        nupdg = self._gen_config.nu_pdg
        cmd = " ".join((
                        os.sep.join((neutgeompath, "genev")),
                        #"-j", beamfile,
                        "-s", filestem, "0", str(N),
                        "-g",
                        geomfile,
                        "-v",
                        "+" + volumename,
                        "-o",
                        outfilename,
                        "-d",
                        str(planenum),
                        "-n",
                        str(numevents),
                        "-f rootracker",
                        #"-f neut",
                        "-i",
                        eventratefile,
                        "-w 1 ", #rewind the flux file
                        "-p", 
                        str(nupdg),
                        ))
        #setupneutcmd = "source /home/software/neut/setupNeut.sh" # TODO : move this to constants somewhere.
        #cmd = " && ".join((setupneutcmd, cmd))
        self._check_call(cmd)
        return

###############################################################################

class CompleteJob(IJob):
    def __init__(self, beam_input, geometry, gen_config, rundir=None, test=False):
        super(CompleteJob, self).__init__(rundir, test)
        self._beam_input = beam_input
        self._geometry = geometry
        self._gen_config = gen_config
        
    def run(self):
        beam_input = self._beam_input
        geometry = self._geometry
        gen_config = self._gen_config
        #job_flux = MergeFluxJob(beam_input, test=self._test)
        rundir = self._rundir
        job_flux = MakeFluxLinks(beam_input, test=self._test, rundir=rundir)
        job_creategeometry = CreateGeometryJob(geometry, test=self._test, rundir=rundir)
        job_evrate = EventRateJob(beam_input, geometry, test=self._test, rundir=rundir)
        job_genev = GenEvJob(gen_config, beam_input, geometry, job_evrate, test=self._test, rundir=rundir)
        jobs = [job_flux,
                job_creategeometry,
                job_evrate,
                job_genev,
                ]
        for j in jobs:
            j.verify()
            j.run()
        return

###############################################################################

def str_from_polarity(polarity):
    r = None
    if polarity < 0:
        r = "antinu"
    else: 
        r = "nu"
    return r

###############################################################################

def getjobname(opt):
    return str_from_polarity(opt.polarity)

###############################################################################

def run(opt):
    test = opt.test
    card = opt.card
    ndid = opt.flux
    radius = opt.radius
    polarity = opt.polarity
    z = opt.z
    nevents = opt.n
    nu_pdg = opt.pdg
    jobname = getjobname(opt)
    #beamcontext = runtime.getcontext().beamcontext
    #nu_flux_files = glob.glob(_abspath("~/t2k/data/irods/QMULZone2/home/hyperk/fluxes/fluka_flux/numode/*.root"))
    #antinu_flux_files = glob.glob(_abspath("~/t2k/data/irods/QMULZone2/home/hyperk/fluxes/fluka_flux/anumode/*.root"))
    nu_flux_files = glob.glob(_abspath("~/t2k/data/hk/ryan_flux/numode/*.root"))
    antinu_flux_files = glob.glob(_abspath("~/t2k/data/hk/ryan_flux/antinumode/*.root"))
    fluxplanes = runtime.FluxPlaneDefinitions()
    fluxplanes.add(runtime.FluxPlane(name="nd2k", baseline=2.04, flukaid=1))
    beamcontext = runtime.BeamContext(jnubeamfiles=runtime.JnuBeamFiles(nu_flux_files, antinu_flux_files), fluxplanes=fluxplanes)
    context = runtime.Context(beamcontext=beamcontext)
    jnubeamfiles = beamcontext.jnubeamfiles()
    jnubeamfiles.verify()
    if polarity == 1:
        filelist = jnubeamfiles.nu_flux_files
    elif polarity == -1:
        filelist = jnubeamfiles.antinu_flux_files
    else:
        raise Exception()
    #print "DEBUG speed up process for debugging"
    #filelist = filelist[0:10]
    rundir = RunDir(card=card)
    beam_input = BeamInput(jobname, filelist)
    #geometry = Geometry(ndid=self._context.DetectorId.ND280, radius=2.0, z=4.0, orientation=Orientation.Z)
    if opt.geometry.lower() == "cylinder":
        geometry = CylinderGeometry(ndid=ndid, radius=radius, z=z, orientation=Orientation.Z, context=context)
    else:
        geometry = CuboidGeometry(ndid=ndid, radius=radius, z=z, orientation=Orientation.Z, context=context)
    gen_config = GenEvConfig(num_events=nevents, nu_pdg=nu_pdg)
    job = CompleteJob(beam_input, geometry, gen_config, test=test, rundir=rundir)
    job.run()
    return

###############################################################################

def parsecml():
    parser = argparse.ArgumentParser()
    parser.add_argument("polarity", type=int, choices=[-1, 1], help="+1 to run neutrino, -1 to run anti-neutrino.", default=1)
    parser.add_argument("radius", type=float, help="Set radius of cyclinder in m.")
    parser.add_argument("z", type=float, help="Set z of cyclinder in m.")
    parser.add_argument("flux", type=str, help="choose flux plane.")
    parser.add_argument("--geometry", type=str, choices=["cylinder", "cuboid"], help="choose geoetry type", default="cuboid")
    parser.add_argument("-c", "--card", type=str, default=None)
    parser.add_argument("-n", "--nevents", dest="n", type=int, default=10000)
    parser.add_argument("-p", "--pdg", dest="pdg", type=int, choices=[-14, -12, 12, 14], default=None)
    parser.add_argument("-t", "--test", dest="test", type=bool, default=False)
    return parser.parse_args()

def main():
    opt = parsecml()
    run(opt)
    return

###############################################################################

if __name__ == "__main__":
    main()

###############################################################################
