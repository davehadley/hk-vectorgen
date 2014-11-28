import glob
import subprocess
import tempfile
import os

###############################################################################
import itertools


def abspath(path):
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

class RunDir:
    def __init__(self, path=None):
        if path is None:
            path = tempfile.mkdtemp(prefix="tmp_hk_vectorgen_")
        self._rundir = abspath(path)

    def rundir(self):
        return self._rundir