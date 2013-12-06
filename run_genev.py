import subprocess
import os
import glob
import argparse
import nuOscillation
import itertools
import collections
import ROOT

DetectorId = nuOscillation.model.constants.DetectorId

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
    def __init__(self, test=False):
        self._test = test
        
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

def _default_planes():
    ret = dict()
    for ndid in DetectorId.ALL:
        if DetectorId.isfardetector(ndid):
            plane = plane_from_ndid(ndid)
            ret[key] = plane
    return ret

###############################################################################

def plane_from_ndid(ndid):
    name = DetectorId.tostring(ndid)
    baseline = DetectorId.baseline(ndid)
    code = DetectorId.toint(ndid)
    plane = BeamPlane(name=name, baseline=baseline, ndcode=code)
    return plane

###############################################################################

class BeamInput:
    def __init__(self, name, file_list):
        self.name = name
        self._file_list = file_list

    def filelist(self):
        return list(itertools.chain.from_iterable(glob.glob(f) for f in self._file_list))
    
    def filename(self):
        return "".join(["merged_beamfiles_", self.name, ".root"])
    
    def verify(self):
        if len(self.filelist()) < 1:
            raise Exception("BeamInput has not matching files", self._file_pattern)
        return

###############################################################################

class MergeFluxJob(IJob):
    def __init__(self, beam_input, test=False):
        super(MergeFluxJob, self).__init__(test)
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
                       self._filename(),
                       " ".join(self._beam_input.filelist()
                        )),
               )
        self._check_call(cmd)
        
###############################################################################

class Orientation:
    Z = "Z"
    Y = "Y"

###############################################################################

class Geometry:
    def __init__(self, ndid, radius=4.0, z=8.0, orientation=Orientation.Z):
        self.ndid = ndid
        self.radius = radius
        self.z = z
        self.orientation = orientation
        self.name = self._uniquestr()
        self._plane = plane_from_ndid(self.ndid)
        
    def verify(self):
        return
    
    def _uniquestr(self):
        return "_".join((DetectorId.tostring(self.ndid),
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
    

###############################################################################

class CreateGeometryJob(IJob):
    def __init__(self, geometry, test=False):
        super(CreateGeometryJob, self).__init__(test)
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
        #get dimensions
        m_to_mm = 1000.0
        radius_mm = g.radius * m_to_mm
        z_mm = g.z * m_to_mm
        #create geometry
        wc_geometry = ROOT.TGeoManager("ND280Geometry","ND280Geometry");
        oxygen = ROOT.TGeoElement("oxygen", "oxygen", 8, 16);
        hydrogen = ROOT.TGeoElement("hydrogen", "hydrogen", 1, 1);
        water = ROOT.TGeoMixture("water", 2, 1);
        water.AddElement(oxygen, 1); 
        water.AddElement(hydrogen, 2); 
        water_med = ROOT.TGeoMedium("water", 1, water);
        cylinder = ROOT.TGeoTube(0, radius_mm, z_mm);
        cylinder1 = ROOT.TGeoTube(0, radius_mm, z_mm);
        t2k = ROOT.TGeoVolume("t2k",cylinder1);
        wc_volume = ROOT.TGeoVolume(volume_name, cylinder, water_med);
        t2k.AddNode(wc_volume, 1);
        wc_geometry.AddVolume(t2k);
        wc_geometry.SetTopVolume(t2k);
        #write geometry to file
        outfile = ROOT.TFile(outfilename,"RECREATE");
        wc_geometry.Write()
        outfile.Close()
        return
    
###############################################################################

class EventRateJob(IJob):
    def __init__(self, beam_input, geometry, test=False):
        super(EventRateJob, self).__init__(test)
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
        if os.path.exists(outfilename):
            os.remove(outfilename)
        if not os.path.exists(outfilename):
            self._create_event_rate()
    
    def _create_event_rate(self):
        outfilename = _abspath(self.filename())
        beamfile = _abspath(self._beam_input.filename())
        geomfile = _abspath(self._geometry.filename())
        volumename = self._geometry.volume_name()
        plane = self._geometry.plane()
        planenum = plane.ndcode
        neutgeompath = os.environ["NEUTGEOM"]
        cmd = " ".join((
                        os.sep.join((neutgeompath, "event_rate")),
                        "-f",
                        beamfile,
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
    
    def verify(self):
        pass

###############################################################################

#/home/software/neut/neut_5.1.4.2/src/neutgeom/genev  -j flux_files.root -g nd280geometry.root -v +Basket -o test.genev.output.2.root -n 10 -f rootracker -i setup_output.root -d 5 2>&1

class GenEvJob(IJob):
    def __init__(self, gen_config, beam_input, geometry, test=False):
        super(GenEvJob, self).__init__(test)
        self._gen_config = gen_config
        self._beam_input = beam_input
        self._geometry = geometry

###############################################################################

class CompleteJob(IJob):
    def __init__(self, beam_input, geometry, gen_config, test=False):
        super(CompleteJob, self).__init__(test)
        self._beam_input = beam_input
        self._geometry = geometry
        self._gen_config = gen_config
        
    def run(self):
        beam_input = self._beam_input
        geometry = self._geometry
        gen_config = self._gen_config
        jobs = [MergeFluxJob(beam_input, test=self._test),
                CreateGeometryJob(geometry, test=self._test),
                EventRateJob(beam_input, geometry, test=self._test),
                GenEvJob(gen_config, beam_input, geometry, test=self._test),
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
    jobname = getjobname(opt)
    beamcontext = nuOscillation.model.runtime.getcontext().beamcontext
    jnubeamfiles = beamcontext.jnubeamfiles()
    if opt.polarity == 1:
        filelist = jnubeamfiles.nu_flux_files
    elif opt.polarity == -1:
        filelist = jnubeamfiles.antinu_flux_files
    beam_input = BeamInput(jobname, filelist)
    #geometry = Geometry(ndid=DetectorId.ND280, radius=2.0, z=4.0, orientation=Orientation.Z)
    geometry = Geometry(ndid=DetectorId.ND2K, radius=4.0, z=8.0, orientation=Orientation.Z)
    gen_config = None
    job = CompleteJob(beam_input, geometry, gen_config, test=False)
    job.run()
    return

###############################################################################

def parsecml():
    parser = argparse.ArgumentParser()
    parser.add_argument("polarity", type=int, choices=[-1, 1], help="+1 to run neutrino, -1 to run anti-neutrino.")
    return parser.parse_args()

def main():
    opt = parsecml()
    run(opt)
    return

###############################################################################

if __name__ == "__main__":
    main()

###############################################################################
