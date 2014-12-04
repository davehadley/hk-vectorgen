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

from vectorgen.jobtools import IJob, RunDir, abspath
from vectorgen.fluxjob import BeamInput, MakeFluxLinks, BeamPlane
from vectorgen.geometryjob import CylinderGeometry, CuboidGeometry, CreateGeometryJob
from vectorgen.neutjob import EventRateJob, NeutEvJob, GenEvConfig, NeutMakeLinks
from vectorgen.constants import Orientation

import warwickcluster

#Jobs:
#    (1) Merge flux files.
#    (2) Create geometry.
#    (3) Run event_rate.
#    (4) Run genev.

###############################################################################

class CompleteJob(IJob):
    def __init__(self, beam_input, geometry, gen_config, rundir=None, test=False, copyflux=False, card=None):
        super(CompleteJob, self).__init__(rundir, test)
        self._beam_input = beam_input
        self._geometry = geometry
        self._gen_config = gen_config
        self._copyflux = copyflux
        self._card = card

    def run(self):
        beam_input = self._beam_input
        geometry = self._geometry
        gen_config = self._gen_config
        #job_flux = MergeFluxJob(beam_input, test=self._test)
        rundir = self._rundir
        job_neutlinks = NeutMakeLinks(rundir, card=self._card)
        job_flux = MakeFluxLinks(beam_input, test=self._test, rundir=rundir, docopy=self._copyflux)
        job_creategeometry = CreateGeometryJob(geometry, test=self._test, rundir=rundir)
        maxfiles = None
        job_evrate = EventRateJob(beam_input, geometry, gen_config, test=self._test, rundir=rundir, maxfiles=maxfiles)
        job_genev = NeutEvJob(gen_config, beam_input, geometry, job_evrate, test=self._test, rundir=rundir, maxfiles=maxfiles)
        jobs = [job_neutlinks,
                job_flux,
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

def run(opt):
    test = opt.test
    plane = BeamPlane(name="nd2k", code=opt.beamplane)
    radius = opt.radius
    polarity = opt.polarity
    z = opt.z
    nevents = opt.n
    #beamcontext = runtime.getcontext().beamcontext
    nu_flux_files = glob.glob(abspath("/data/t2k/hk/irods/hk2/home/hyperk/fluxes/flux_2km/plus_minus320kA/t2hk_320a_2km_fluka2011_*.root"))
    antinu_flux_files = glob.glob(abspath("/data/t2k/hk/irods/hk2/home/hyperk/fluxes/flux_2km/plus_minus320kA/t2hk_m320a_2km_fluka2011_*.root"))
    #nu_flux_files = glob.glob(abspath("~/t2k/data/irods/hk2/home/hyperk/fluxes/flux_2km/plus_minus320kA/t2hk_m320a_2km_fluka2011_*1.root"))[0:10]
    #antinu_flux_files = glob.glob(abspath("~/t2k/data/irods/hk2/home/hyperk/fluxes/flux_2km/plus_minus320kA/t2hk_320a_2km_fluka2011_*1.root"))[0:10]
    if polarity == 1:
        filelist = nu_flux_files
    elif polarity == -1:
        filelist = antinu_flux_files
    else:
        raise Exception()
    #print "DEBUG speed up process for debugging"
    #filelist = filelist[0:10]
    rundir = RunDir()
    beam_input = BeamInput(str_from_polarity(opt.polarity), filelist, plane)
    #geometry = Geometry(ndid=self._context.DetectorId.ND280, radius=2.0, z=4.0, orientation=Orientation.Z)
    if opt.geometry.lower() == "cylinder":
        geometry = CylinderGeometry(radius=radius, z=z, orientation=Orientation.Z)
    else:
        geometry = CuboidGeometry(radius=radius, z=z, orientation=Orientation.Z)
    gen_config = GenEvConfig(num_events=nevents, run_num=opt.runnum, nu_pdg=opt.pdg)
    job = CompleteJob(beam_input, geometry, gen_config, test=test, rundir=rundir, copyflux=opt.copyflux, card=opt.card)
    job.run()
    return

###############################################################################

def submitjobs(opt):
    startrunnum = opt.runnum
    jobs = []
    for i in xrange(opt.njobs):
        runnum = startrunnum + i
        cmd = "python -m vectorgen.run_neut " + " ".join([
            str(opt.polarity),
            str(opt.radius),
            str(opt.z),
            "--geometry=%s" % (opt.geometry),
            "--nevents=%s" % (opt.n),
            "--runnum=%s" % (runnum),
            "--copyflux", #copy flux files to /tmp before running the job.
            "--beamplane=%s" % opt.beamplane,
        ])
        if opt.pdg is not None:
            cmd += " --pdg=" + str(opt.pdg) + " "
        queue = "long"
        name = "job_run_genie_" + "_".join([str(opt.polarity), str(opt.pdg), str(runnum)])
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
    parser.add_argument("--beamplane", type=int, default=1, help="choose flux plane.")
    parser.add_argument("--geometry", type=str, choices=["cylinder", "cuboid"], help="choose geoetry type", default="cylinder")
    parser.add_argument("-p", "--pdg", dest="pdg", type=int, choices=[-14, -12, 12, 14], default=None)
    parser.add_argument("-n", "--nevents", dest="n", type=int, default=10000)
    parser.add_argument("-t", "--test", dest="test", type=bool, default=False)
    parser.add_argument("-r", "--runnum", dest="runnum", type=int, default=1000)
    parser.add_argument("-b", "--batch", dest="batch", action="store_true", default=False)
    parser.add_argument("-j", "--njobs", dest="njobs", type=int, default=100)
    parser.add_argument("--copyflux", dest="copyflux", action="store_true", default=False)
    parser.add_argument("-c", "--card", type=str, default=None)
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
