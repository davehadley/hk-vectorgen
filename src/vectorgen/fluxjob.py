import glob
import itertools
import os
import shutil
import ROOT

from vectorgen.jobtools import IJob, mirror_file
from vectorgen.filelock import FileLock


###############################################################################

class BeamPlane:
    def __init__(self, name, code):
        self.name = name
        self.code = code

###############################################################################

class BeamInput:
    def __init__(self, name, file_list, plane):
        self.name = name
        self._file_list = file_list
        self._plane = plane

    def filelist(self):
        return sorted(list(itertools.chain.from_iterable(glob.glob(f) for f in self._file_list)))

    def filename(self):
        return "".join(["merged_beamfiles_", self.name, ".root"])

    def linkdir(self):
        return "flux_links_" + self.filename().replace(".root", "")

    def verify(self):
        if len(self.filelist()) < 1:
            raise Exception("BeamInput has not matching files", self._file_pattern)
        return

    def filestem(self):
        return  "".join((self.linkdir(), os.path.sep, "nu.nd280"))

    def planenum(self):
        return self._plane.code

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
        filelock = FileLock("/tmp/.lock_file_hk_vectorgen_flux", timeout=60*10)
        filelock.lock()
        try:
            os.makedirs(outdir)
        except os.error:
            #ignore as this happens if directory already exists
            pass
        for i, fname in enumerate(self._beam_input.filelist()):
            if self._n is not None and i >= self._n:
                break
            src = fname
            if self._copy:
                src = self._mirror(src)
            dst = "".join((self._rundir.rundir(), os.sep, self._beam_input.filestem(), ".", str(i), ".root"))
            if not os.path.exists(dst):
                os.symlink(src, dst)
        filelock.release()
        self._checkfiles()
        return

    def _checkfiles(self):
        for fname in self._beam_input.filelist():
            if self._copy:
                fname = self._mirror(fname)
            self._checkbeamfile(fname)
        return

    def _checkbeamfile(self, fname):
        tfile = ROOT.TFile(fname)
        if not tfile.IsOpen():
            raise Exception("failed to open", fname)
        tree = tfile.Get("h2000")
        if not tree:
            raise Exception("failed to get tree", fname)
        if not tree.GetEntries():
            raise Exception("tree has no entries", fname, tree)
        return

    def _mirror(self, fname):
        return mirror_file(fname)
