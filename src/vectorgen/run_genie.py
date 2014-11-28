import uuid
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

from simplot.batch import warwickcluster

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
    def __init__(self, path=None):
        if path is None:
            path = tempfile.mkdtemp(prefix="tmp_run_genev_")
        self._rundir = _abspath(path)

    def rundir(self):
        return self._rundir

###############################################################################

class BeamInput:
    def __init__(self, name, file_list):
        self.name = name
        self._file_list = file_list

    def filelist(self):
        return sorted(list(itertools.chain.from_iterable(glob.glob(f) for f in self._file_list)))
    
    def filename(self):
        return "".join(["merged_beamfiles_", self.name, ".root"])
    
    def linkdir(self):
        return "flux_links_" + self.filename().replace(".root", "")
    
    def verify(self):
        if len(self.filelist()) < 1:
            raise Exception("BeamInput has not matching files", self._file_list)
        return
    
    def filestem(self):
        return  "".join((self.linkdir(), os.path.sep, "nu.nd280"))

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
    def __init__(self, beam_input, rundir=None, test=False, n=None, docopy=False):
        super(MakeFluxLinks, self).__init__(rundir, test)
        self._n = n
        self._beam_input = beam_input
        self._copy = docopy
        
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
                if self._copy:
                    shutil.copy2(src, dst)
                else:
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
        vacuum = ROOT.TGeoMaterial("vacuum", 0,0,0);
        vacuum_med = ROOT.TGeoMedium("vacuum", 2, vacuum)
        vol0, vol1 = g.build_detector_volume()
        t2k = ROOT.TGeoVolume("t2k",vol1, vacuum_med);
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
    def __init__(self, beam_input, geometry, rundir=None, test=False, maxfiles=None):
        super(EventRateJob, self).__init__(rundir, test)
        self._beam_input = beam_input
        self._geometry = geometry
        self._maxfiles = maxfiles
        
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
        try:
            os.makedirs(os.path.dirname(outfilename))
        except OSError:
            pass
        tmpdir = tempfile.mkdtemp(suffix="run_genie")
        tmpoutfilename = os.sep.join([tmpdir, self.filename()])
        beamfile = _abspath(self._beam_input.filename())
        filestem = self._beam_input.filestem()
        filestem = "".join((self._rundir.rundir(), os.sep, filestem))
        #N = len(self._beam_input.filename()) - 1
        N = len(glob.glob(filestem + "*.root")) - 1
        if N <= 0:
            raise Exception("No flux files matching", filestem)
        if self._maxfiles is not None:
            N = min(N, self._maxfiles)
        geomfile = _abspath(self._geometry.filename())
        volumename = self._geometry.volume_name()
        plane = self._geometry.plane()
        planenum = plane.ndcode
        geniepath = os.environ["GENIE"]
        splinesfile = os.environ["GENIE_SPLINES"]
        cmd = " ".join((
                        os.sep.join((geniepath, "bin", "gevgen_t2k")),
                        #"-n", numevents,
                        #"--run", runnum,
                        "--seed", str(1721827),
                        "--cross-sections", splinesfile,
                        "--message-thresholds ${GENIE}/config/Messenger_laconic.xml",
                        "-g", geomfile,
                        "--event-generator-list DefaultWithMEC",
                        "-f " + "".join([filestem,"@", str(0), "@", str(N), ",nd" + str(planenum)]),
                        "-t", volumename,
                        "-S", tmpoutfilename,
                        ))
        self._check_call(cmd)
        #copy output to destination
        shutil.copy2(tmpoutfilename, outfilename)
        atomicmove(tmpoutfilename, outfilename)
        return
    
    def verify(self):
        pass

###############################################################################

class GenEvConfig:
    def __init__(self, num_events, run_num, nu_pdg=None, startseed=1721827):
        self.num_events = num_events
        self.run_num = run_num
        self._startingseed = startseed
        self._nu_pdg = nu_pdg

    @property
    def seed(self):
        return self._startingseed + self.run_num
    
    @property
    def name(self):
        return "_".join((str(self.nu_pdg_name), str(self._startingseed), str(self.run_num),
        ))

    @property
    def nu_pdg_code(self):
        return self._nu_pdg

    @property
    def nu_pdg_name(self):
        _nu_pdg_names = {None : "allnuflav",
                         0 : "allnuflav",
                         12 : "nue",
                         14 : "numu",
                         -12 : "antinue",
                         -14 : "antinumu",
                              }
        return _nu_pdg_names[self._nu_pdg]
    
###############################################################################

class GenieEvJob(IJob):
    def __init__(self, gen_config, beam_input, geometry, eventratejob, rundir=None, test=False, maxfiles=None):
        super(GenieEvJob, self).__init__(rundir, test)
        self._gen_config = gen_config
        self._beam_input = beam_input
        self._geometry = geometry
        self._eventrate = eventratejob
        self._maxfiles = maxfiles

    def filename(self):
        beamname = self._beam_input.name
        geomname = self._geometry.name
        configname = self._gen_config.name
        outfilename = "_".join(("genie",
                               beamname,
                               geomname,
                               configname,
                               )) + ".0.ghep.root"
        return outfilename
    
    def run(self):
        outfilename = self.filename()
        if not os.path.exists(outfilename):
            self._create_genev()
        return

#genev  -j flux_files.root -g nd280geometry.root -v +Basket -o test.genev.output.2.root -n 10 -f rootracker -i setup_output.root -d 5 2>&1

    def _create_genev(self):
        outfilename = _abspath(self.filename())
        try:
            os.makedirs(os.path.dirname(outfilename))
        except OSError:
            pass
        tmpdir = tempfile.mkdtemp("run_genie")
        tmpoutfilename = os.sep.join([tmpdir, os.path.basename(outfilename)])
        beamfile = _abspath(self._beam_input.filename())
        filestem = self._beam_input.filestem()
        filestem = "".join((self._rundir.rundir(), os.sep, filestem))
        #N = len(self._beam_input.filename()) - 1
        N = len(glob.glob(filestem + "*.root")) - 1
        if N <= 0:
            raise Exception("No flux files matching", filestem)
        if self._maxfiles is not None:
            N = min(N, self._maxfiles)
        geomfile = _abspath(self._geometry.filename())
        eventratefile = _abspath(self._eventrate.filename())
        volumename = self._geometry.volume_name()
        plane = self._geometry.plane()
        planenum = plane.ndcode
        numevents = self._gen_config.num_events
        geniepath = os.environ["GENIE"]
        splinesfile = os.environ["GENIE_SPLINES"]
        outfileprefix = tmpoutfilename.split(".0.ghep.root")[0]
        fluxstr = "-f " + "".join([filestem,"@", str(0), "@", str(N), ",nd" + str(planenum)])
        if self._gen_config.nu_pdg_code is not None:
            fluxstr += "," + str(self._gen_config.nu_pdg_code)
        cmd = " ".join((
                        os.sep.join((geniepath, "bin", "gevgen_t2k")),
                        "-n", str(numevents),
                        "--run", str(self._gen_config.run_num),
                        "--seed", str(self._gen_config.seed),
                        "--cross-sections", splinesfile,
                        "--message-thresholds ${GENIE}/config/Messenger_laconic.xml",
                        "-g", geomfile,
                        "--event-generator-list DefaultWithMEC",
                        fluxstr,
                        "-t", volumename,
                        "-P", eventratefile,
                        "-o", outfileprefix
                        ))
        self._check_call(cmd)
        #copy contents of temporary directory to the destination directory
        for root, dirs, files in os.walk(tmpdir):
            for fname in files:
                f = os.path.join(root, fname)
                atomicmove(f, os.path.dirname(outfilename))
        return

def atomicmove(srcname, dest):
    #determine destination dir and name
    if os.path.exists(dest) and os.path.isdir(dest):
        destdir = dest
        destname = os.path.basename(srcname)
    else:
        destdir = os.path.dirname(dest)
        destname = os.path.basename(dest)
    destname = os.sep.join([destdir, destname])
    #copy file to destination directory
    uniquestr = str(uuid.uuid4()).replace("-", "")
    shutil.move(srcname, destname + uniquestr)
    #rename back to the proper name
    os.rename(destname + uniquestr, destname)
    return

###############################################################################

class ConvertGenieEvJob(IJob):
    def __init__(self, genevjob, rundir=None, test=False):
        super(ConvertGenieEvJob, self).__init__(rundir, test)
        self._infilename = genevjob.filename()
        
    def filename(self):
        outfilename = self._infilename.replace("genie", "genev")
        if outfilename == self._infilename:
            raise Exception(outfilename)
        return outfilename
    
    def run(self):
        outfilename = self.filename()
        if not os.path.exists(outfilename):
            self._convert()
        return

#genev  -j flux_files.root -g nd280geometry.root -v +Basket -o test.genev.output.2.root -n 10 -f rootracker -i setup_output.root -d 5 2>&1

    def _convert(self):
        infilename = _abspath(self._infilename)
        outfilename = _abspath(self.filename())
        cmd = "gntpc -i %s -o %s -f t2k_rootracker" % (infilename, outfilename)
        self._check_call(cmd)
        return

###############################################################################

class CompleteJob(IJob):
    def __init__(self, beam_input, geometry, gen_config, rundir=None, test=False, copyflux=False):
        super(CompleteJob, self).__init__(rundir, test)
        self._beam_input = beam_input
        self._geometry = geometry
        self._gen_config = gen_config
        self._copyflux = copyflux
        
    def run(self):
        beam_input = self._beam_input
        geometry = self._geometry
        gen_config = self._gen_config
        #job_flux = MergeFluxJob(beam_input, test=self._test)
        rundir = self._rundir
        job_flux = MakeFluxLinks(beam_input, test=self._test, rundir=rundir, docopy=self._copyflux)
        job_creategeometry = CreateGeometryJob(geometry, test=self._test, rundir=rundir)
        maxfiles = None
        job_evrate = EventRateJob(beam_input, geometry, test=self._test, rundir=rundir, maxfiles=maxfiles)
        job_genev = GenieEvJob(gen_config, beam_input, geometry, job_evrate, test=self._test, rundir=rundir, maxfiles=maxfiles)
        job_convert = ConvertGenieEvJob(job_genev, test=self._test, rundir=rundir)
        jobs = [job_flux,
                job_creategeometry,
                job_evrate,
                job_genev,
                job_convert,
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
    ndid = opt.flux
    radius = opt.radius
    polarity = opt.polarity
    z = opt.z
    nevents = opt.n
    jobname = getjobname(opt)
    #beamcontext = runtime.getcontext().beamcontext
    nu_flux_files = glob.glob(_abspath("/data/t2k/hk/irods/hk2/home/hyperk/fluxes/flux_2km/plus_minus320kA/t2hk_320a_2km_fluka2011_*.root"))
    antinu_flux_files = glob.glob(_abspath("/data/t2k/hk/irods/hk2/home/hyperk/fluxes/flux_2km/plus_minus320kA/t2hk_m320a_2km_fluka2011_*.root"))
    #nu_flux_files = glob.glob(_abspath("~/t2k/data/irods/hk2/home/hyperk/fluxes/flux_2km/plus_minus320kA/t2hk_m320a_2km_fluka2011_*1.root"))
    #antinu_flux_files = glob.glob(_abspath("~/t2k/data/irods/hk2/home/hyperk/fluxes/flux_2km/plus_minus320kA/t2hk_320a_2km_fluka2011_*1.root"))
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
    rundir = RunDir()
    beam_input = BeamInput(jobname, filelist)
    #geometry = Geometry(ndid=self._context.DetectorId.ND280, radius=2.0, z=4.0, orientation=Orientation.Z)
    if opt.geometry.lower() == "cylinder":
        geometry = CylinderGeometry(ndid=ndid, radius=radius, z=z, orientation=Orientation.Z, context=context)
    else:
        geometry = CuboidGeometry(ndid=ndid, radius=radius, z=z, orientation=Orientation.Z, context=context)
    gen_config = GenEvConfig(num_events=nevents, run_num=opt.runnum, nu_pdg=opt.pdg)
    job = CompleteJob(beam_input, geometry, gen_config, test=test, rundir=rundir, copyflux=opt.copyflux)
    job.run()
    return

###############################################################################

def submitjobs(opt):
    startrunnum = opt.runnum
    jobs = []
    for i in xrange(opt.njobs):
        runnum = startrunnum + i
        cmd = "python -m vectorgen.run_genie " + " ".join([
            str(opt.polarity),
            str(opt.radius),
            str(opt.z),
            str(opt.flux),
            "--geometry=%s" % (opt.geometry),
            "--nevents=%s" % (opt.n),
            "--runnum=%s" % (runnum),
            "--copyflux", #copy flux files to /tmp before running the job.
        ])
        if opt.pdg is not None:
            cmd += " --pdg=" + str(opt.pdg) + " "
        queue = "medium"
        name = "genie" + "_".join([str(opt.polarity), str(runnum)])
        j = warwickcluster.ClusterJob(name, queue, cmd)
        jobs.append(j)
    #run the jobs
    for j in jobs:
        j.submit()
    return

###############################################################################

def parsecml():
    parser = argparse.ArgumentParser()
    parser.add_argument("polarity", type=int, choices=[-1, 1], help="+1 to run neutrino, -1 to run anti-neutrino.", default=1)
    parser.add_argument("radius", type=float, help="Set radius of cyclinder in m.")
    parser.add_argument("z", type=float, help="Set z of cyclinder in m.")
    parser.add_argument("flux", type=str, choices=["nd2k"], help="choose flux plane.")
    parser.add_argument("--geometry", type=str, choices=["cylinder", "cuboid"], help="choose geoetry type", default="cylinder")
    parser.add_argument("-p", "--pdg", dest="pdg", type=int, choices=[-14, -12, 12, 14], default=None)
    parser.add_argument("-n", "--nevents", dest="n", type=int, default=10000)
    parser.add_argument("-t", "--test", dest="test", type=bool, default=False)
    parser.add_argument("-r", "--runnum", dest="runnum", type=int, default=1000)
    parser.add_argument("-b", "--batch", dest="batch", action="store_true", default=False)
    parser.add_argument("-j", "--njobs", dest="njobs", type=int, default=100)
    parser.add_argument("--copyflux", dest="copyflux", action="store_true", default=False)
    return parser.parse_args()

###############################################################################

def main():
    opt = parsecml()
    if opt.batch:
        submitjobs(opt)
    else:
        run(opt)
    return

###############################################################################

if __name__ == "__main__":
    main()

###############################################################################
