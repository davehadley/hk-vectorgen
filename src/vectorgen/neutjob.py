import glob
import os
import shutil
import tempfile
import uuid

from vectorgen.jobtools import IJob, abspath, mirror_file_with_lock

###############################################################################

class NeutMakeLinks(IJob):
    def __init__(self, rundir, card=None, test=False):
        super(NeutMakeLinks, self).__init__(rundir, test)
        self._card = self._find_card(card)
        self._run_make_links()

    def _find_card(self, card):
        if card is None:
            #use default from NEUTGEOM directory
            card = "".join((os.environ["NEUT_ROOT"], os.sep, "src/neutgeom/neut.card"))
        card = abspath(card)
        if not os.path.exists(card):
            raise Exception("Cannot find card file", card)
        return card

    def _run_make_links(self):
        outdir = self._rundir.rundir()
        try:
            os.makedirs(outdir)
        except os.error:
            #ignore as this happens if directory already exists
            pass
        inputdir = "".join((os.environ["NEUT_ROOT"], os.sep, "src", os.sep, "neutsmpl"))
        for fname in os.listdir(inputdir):
            src = "".join((inputdir, os.sep, fname))
            if os.path.islink(src):
                dst = "".join((outdir, os.sep, fname))
                if not os.path.exists(dst):
                    os.symlink(src, dst)
        src = self._card
        dst = "".join((outdir, os.sep, "neut.card"))
        shutil.copyfile(src, dst)
        return

###############################################################################

class EventRateJob(IJob):
    def __init__(self, beam_input, geometry, gen_config, rundir=None, test=False, maxfiles=None):
        super(EventRateJob, self).__init__(rundir, test)
        self._beam_input = beam_input
        self._geometry = geometry
        self._maxfiles = maxfiles
        self._gen_config = gen_config

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
        tmpdir = tempfile.mkdtemp(suffix="run_neut")
        tmpoutfilename = os.sep.join([tmpdir, self.filename()])
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
        neutgeompath = os.environ["NEUTGEOM"]
        cmd = " ".join((
                        os.sep.join((neutgeompath, "event_rate")),
                        "-s", filestem, "0", str(N),
                        "-g", geomfile,
                        "-v", "+" + volumename,
                        "-o", tmpoutfilename,
                        "-d", str(planenum),
                        ))
        self._check_call(cmd)
        #copy output to destination
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

class NeutEvJob(IJob):
    def __init__(self, gen_config, beam_input, geometry, eventratejob, rundir=None, test=False, maxfiles=None):
        super(NeutEvJob, self).__init__(rundir, test)
        self._gen_config = gen_config
        self._beam_input = beam_input
        self._geometry = geometry
        self._eventrate = eventratejob
        self._maxfiles = maxfiles

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
        outfilename = abspath(self.filename())
        try:
            os.makedirs(os.path.dirname(outfilename))
        except OSError:
            pass
        tmpdir = tempfile.mkdtemp("run_neut")
        tmpoutfilename = os.sep.join([tmpdir, os.path.basename(outfilename)])
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
        numevents = self._gen_config.num_events
        neutgeompath = os.environ["NEUTGEOM"]
        nupdg = self._gen_config.nu_pdg_code
        if nupdg is None:
            nupdg = 0
        #get event rate file and copy to temp
        eventratefile = abspath(self._eventrate.filename())
        eventratefile = mirror_file_with_lock(eventratefile)
        #form genev command
        cmd = " ".join((
                        os.sep.join((neutgeompath, "genev")),
                        "-s", filestem, "0", str(N),
                        "-g", geomfile,
                        "-v", "+" + volumename,
                        "-o", tmpoutfilename,
                        "-d", str(planenum),
                        "-n", str(numevents),
                        "-f rootracker",
                        "-i", eventratefile,
                        "-w 1 ", #rewind the flux file
                        " -p", str(nupdg),
                        " -r", str(self._gen_config.seed),
                        ))
        self._check_call(cmd)
        atomicmove(tmpoutfilename, outfilename)
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



