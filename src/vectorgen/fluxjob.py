import glob
import itertools
import os
import shutil

from vectorgen.jobtools import IJob

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
            raise Exception("BeamInput has not matching files", self._file_pattern)
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

