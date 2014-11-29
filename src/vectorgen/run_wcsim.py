import argparse
import shutil
import tempfile
import os
import glob
import ROOT

from vectorgen import genev_to_nuance_wcsim
from vectorgen.jobtools import IJob, abspath

###############################################################################

_MACRO_TEMPLATE = '''
/run/verbose 0
/tracking/verbose 0
/hits/verbose 0

/WCSim/WCgeom TITUS_2kton_10inch_40perCent
/WCSim/Construct
/WCSim/PMTQEMethod     Stacking_Only
/WCSim/SavePi0 true
/DarkRate/SetDarkRate 4 kHz

/mygen/vecfile {inputfilename}
/WCSimIO/RootFile {outputfilename}
/run/beamOn {numberofevents}
'''

###############################################################################

class WCSimJob(IJob):
    def __init__(self, inputfilename, outputfilename=None, nevents=None):
        self._maxevents = nevents
        #create a temp directory to work in
        self._tmpdir = tempfile.mkdtemp(prefix="tmp_vectorgen_run_wcsim")
        super(WCSimJob, self).__init__(rundir=self._tmpdir)
        #work with absolute path names
        inputfilename = abspath(inputfilename)
        if outputfilename:
            outputfilename = abspath(outputfilename)
        #store in and out file names
        self._originalinputfilename = str(inputfilename)
        self._inputfilename = str(inputfilename)
        if outputfilename is None:
            outputfilename = self._autooutputfilename(inputfilename)
        self._finaloutputfilename = outputfilename
        self._outputfilename = os.path.sep.join([self._tmpdir, os.path.basename(outputfilename)])



    def _autooutputfilename(self, fname):
        dn = os.path.dirname(fname)
        bn = os.path.basename(fname)
        fn = bn.replace("genev", "wcsim")
        if bn == fn:
            raise Exception("unable to automatically determine output file name from input file name", fname)
        result = os.sep.join([dn, fn])
        return abspath(result)

    def run(self):
        self._copyinputtotemp()
        datfile = self._create_datfile()
        macfile = self._create_macro(datfile)
        self._runwcsim(macfile)
        self._copy_original_genev_tree()
        self._copy_result_to_dest()
        return

    def _get_genev_tree(self, tfile):
        tree = None
        for tname in ["nRooTracker", "gRooTracker"]:
            tree = tfile.Get(tname)
            if tree:
                break
        if not tree:
            tfile.ls()
            raise Exception("failed to load tree", self._inputfilename)
        return tree

    def _copyinputtotemp(self):
        src = self._originalinputfilename
        fname = os.path.basename(src)
        dst = os.sep.join([self._tmpdir, fname])
        shutil.copy2(src, dst)
        self._inputfilename = dst
        return

    def _create_datfile(self):
        genev_to_nuance_wcsim.loadlib()
        infile = self._inputfilename
        #get output filename
        fname = os.path.basename(infile).replace(".root", ".dat")
        if fname == infile:
            raise Exception("Error: failed to determine correct filename for dat file.", fname)
        outfile = os.sep.join([self._tmpdir, fname])
        #run conversion
        genev_to_nuance_wcsim.chktrack2(infile, outfile)
        return outfile

    def _create_macro(self, datfile):
        #get macro contents
        inputfilename = datfile
        outputfilename = self._outputfilename
        nevents = self._nevents()
        macstr = _MACRO_TEMPLATE.format(inputfilename=inputfilename, outputfilename=outputfilename, numberofevents=nevents)
        #write macro to a file
        macfname = os.sep.join([self._tmpdir, "run_wcsim.mac"])
        print >>open(macfname, "w"), macstr
        return macfname

    def _nevents(self):
        if self._maxevents is not None:
            return self._maxevents
        else:
            #determine number of events from the input file
            tfile = ROOT.TFile(self._inputfilename)
            if not tfile.IsOpen():
                raise Exception("failed to open ROOT file", self._inputfilename)
            tree = self._get_genev_tree(tfile)
            return tree.GetEntries()

    def _runwcsim(self, macfile):
        wcsimbin = self._findwcsimbinary()
        cmd = " ".join([wcsimbin, macfile])
        self._check_call(cmd, workingdir=self._tmpdir)
        return


    def _findwcsimbinary(self):
        wcsimbin = None
        try:
            wcsimdir = os.environ["WCSIMDIR"]
        except KeyError:
            raise Exception("No WCSIMDIR environment variable set.")
        pattern = os.sep.join([wcsimdir,"exe", "bin", "*", "WCSim"])
        matches = glob.glob(pattern)
        if len(matches) > 0:
            wcsimbin = matches[0]
        else:
            raise Exception("Cannot find WCSim executable. No matches to pattern", pattern)
        return wcsimbin

    def _copy_result_to_dest(self):
        src = self._outputfilename
        dst = self._finaloutputfilename
        shutil.move(src, dst)
        return

    def _copy_original_genev_tree(self):
        oldfile = ROOT.TFile(self._inputfilename)
        if not oldfile.IsOpen():
            raise Exception("failed to open ROOT file", self._inputfilename)
        oldtree = self._get_genev_tree(oldfile)
        newfile = ROOT.TFile(self._outputfilename, "UPDATE")
        newtree = oldtree.CloneTree()
        newtree.Write()
        newfile.Close()
        oldfile.Close()
        return

###############################################################################

def parsecml():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input files.")
    parser.add_argument("-o", help="Output filename", default=None)
    parser.add_argument("--maxevents",  help="Set a maximum number of events to process", type=int, default=None)
    parser.add_argument("--batch" , action="store_true", help="Submit job to batch system.")
    return parser.parse_args()

def run(opt):
    j = WCSimJob(opt.input, opt.o, nevents=opt.maxevents)
    j.run()
    return

def submitjobs(opt):
    pass

def main():
    opt = parsecml()
    if opt.batch:
        submitjobs(opt)
    else:
        run(opt)
    return

if __name__ == "__main__":
    main()
