import glob
import os
import shutil
import tempfile
import uuid

from vectorgen.jobtools import IJob, abspath

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
        pdgname = self._gen_config.nu_pdg_name
        outfilename = "_".join(("eventrate",
                               beamname,
                               geomname,
                               pdgname,
                               )) + ".root"
        return outfilename

    def run(self):
        outfilename = self.filename()
        if not os.path.exists(outfilename):
            self._create_event_rate()
        return

    def _create_event_rate(self):
        outfilename = abspath(self.filename())
        try:
            os.makedirs(os.path.dirname(outfilename))
        except OSError:
            pass
        tmpdir = tempfile.mkdtemp(suffix="run_genie")
        tmpoutfilename = os.sep.join([tmpdir, self.filename()])
        beamfile = abspath(self._beam_input.filename())
        filestem = self._beam_input.filestem()
        filestem = "".join((self._rundir.rundir(), os.sep, filestem))
        #N = len(self._beam_input.filename()) - 1
        N = len(glob.glob(filestem + "*.root")) - 1
        if N <= 0:
            raise Exception("No flux files matching", filestem)
        if self._maxfiles is not None:
            N = min(N, self._maxfiles)
        geomfile = abspath(self._geometry.filename())
        volumename = self._geometry.volume_name()
        planenum = self._beam_input.planenum()
        geniepath = os.environ["GENIE"]
        splinesfile = os.environ["GENIE_SPLINES"]
        fluxstr = "-f " + "".join([filestem,"@", str(0), "@", str(N), ",nd" + str(planenum)])
        if self._gen_config.nu_pdg_code is not None:
            fluxstr += "," + str(self._gen_config.nu_pdg_code)
        cmd = " ".join((
                        os.sep.join((geniepath, "bin", "gevgen_t2k")),
                        #"-n", numevents,
                        #"--run", runnum,
                        "--seed", str(1721827),
                        "--cross-sections", splinesfile,
                        "--message-thresholds ${GENIE}/config/Messenger_laconic.xml",
                        "-g", geomfile,
                        "--event-generator-list DefaultWithMEC",
                        fluxstr,
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
        outfilename = abspath(self.filename())
        try:
            os.makedirs(os.path.dirname(outfilename))
        except OSError:
            pass
        tmpdir = tempfile.mkdtemp("run_genie")
        tmpoutfilename = os.sep.join([tmpdir, os.path.basename(outfilename)])
        beamfile = abspath(self._beam_input.filename())
        filestem = self._beam_input.filestem()
        filestem = "".join((self._rundir.rundir(), os.sep, filestem))
        #N = len(self._beam_input.filename()) - 1
        N = len(glob.glob(filestem + "*.root")) - 1
        if N <= 0:
            raise Exception("No flux files matching", filestem)
        if self._maxfiles is not None:
            N = min(N, self._maxfiles)
        geomfile = abspath(self._geometry.filename())
        eventratefile = abspath(self._eventrate.filename())
        volumename = self._geometry.volume_name()
        planenum = self._beam_input.planenum()
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
        infilename = abspath(self._infilename)
        outfilename = abspath(self.filename())
        cmd = "gntpc -i %s -o %s -f t2k_rootracker" % (infilename, outfilename)
        self._check_call(cmd)
        return

###############################################################################

