import argparse
import math
import StringIO

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

class ParticleType:
    NUMU = "numu"
    MU = "mu-"
    E = "e-"
    PI0 = "pi0"
    PIP = "pi+"
    PIM = "pi-"
    PROTON = "proton"
    NEUTRON = "neutron"

    ALL = [NUMU, MU, E, PI0, PIP, PIM, PROTON, NEUTRON]

    _pdg = {NUMU:14,
            MU:13,
            E:12,
            PI0:111,
            PIP:211,
            PIM:-211,
            PROTON:2212,
            NEUTRON:2112,
    }

    _mass = {NUMU:0.0,
            MU: 105.658e-3,
            E:   0.511e-3,
            PI0: 134.98e-3,
            PIP: 139.57e-3,
            PIM: 139.57e-3,
            PROTON: 938.272e-3,
            NEUTRON: 939.565e-3,
    }

    @classmethod
    def pdg(cls, p):
        return ParticleType._pdg[p]

    @classmethod
    def mass(cls, p):
        return ParticleType._mass[p]

class Vector:
    def __init__(self, xvec, pvec, pdg):
        self.xvec = xvec
        self.pvec = pvec
        self.pdg = pdg

class ParticleGunJob:
    def __init__(self, options):
        self._options = options

    def run(self):
        self._rand = ROOT.TRandom3(self._options.seed)
        with open(self._options.outfilename,"w") as outfile:
            for eventnum in xrange(self._options.n):
                print >>outfile, self._gen_nuance_event(eventnum)

    def _gen_nuance_event(self, eventnum):
        m_to_cm = 100.0
        GeV_to_MeV = 1000.0
        vector = self._gen_vector()
        pvec = vector.pvec
        xvec = vector.xvec
        pdg = vector.pdg
        #format event data into nuance string
        sio = StringIO.StringIO()
        mode = 99
        sio.write("$ begin\n")
        sio.write("$ nuance %d\n" % mode)
        sio.write("$ vertex %5.4f %5.4f %5.4f %5.4f\n" % (m_to_cm*xvec.X(),
                                                          m_to_cm*xvec.Y(),
                                                          m_to_cm*xvec.Z(),
                                                          xvec.T()
                                                        )
        )
        p = pvec.P()
        if p == 0:
            p = 0.000001;
        cosx = pvec.Px() / p;
        cosy = pvec.Py() / p;
        cosz = pvec.Pz() / p;
        #create a fake initial state
        stat = -1
        sio.write("$ track %d %5.4f %5.4f %5.4f %5.4f %d\n" % (ParticleType.pdg(ParticleType.NUMU),
                                                               GeV_to_MeV*pvec.E(),
                                                               cosx,
                                                               cosy,
                                                               cosz,
                                                               stat)
        )
        sio.write("$ track %d %5.4f %5.4f %5.4f %5.4f %d\n" % (ParticleType.pdg(ParticleType.PROTON),
                                                               GeV_to_MeV*ParticleType.mass(ParticleType.PROTON),
                                                               0.0,
                                                               0.0,
                                                               1.0,
                                                               stat)
        )
        #write the final state
        stat = 0
        sio.write("$ info 0 0 %d\n" % eventnum)
        sio.write("$ track %d %5.4f %5.4f %5.4f %5.4f %d\n" % (pdg, pvec.E(), cosx, cosy, cosz, stat))
        #sio.write("$ end %d" % eventnum)
        sio.write("$ end")
        return sio.getvalue()

    def _gen_vector(self):
        xvec = self._gen_position()
        pdg = self._get_pdg()
        pvec = self._gen_momentum()
        return Vector(xvec=xvec, pvec=pvec, pdg=pdg)

    def _get_pdg(self):
        return ParticleType.pdg(self._options.particle)

    def _gen_momentum(self):
        dirvec = self._gen_direction()
        #generate total momentum
        p = self._rand.Uniform(self._options.pmin, self._options.pmax)
        #calculate energy
        m = ParticleType.mass(self._options.particle)
        e = math.sqrt(p**2 + m**2)
        pvec = ROOT.TLorentzVector(p * dirvec.X(),
                                   p * dirvec.Y(),
                                   p * dirvec.Z(),
                                   e,
                                   )
        return pvec

    def _gen_direction(self):
        #generate random direction
        theta = self._rand.Uniform(0.0, 2.0*math.pi)
        z = self._rand.Uniform(-1, 1)
        x = math.sqrt(1-z**2) * math.cos(theta)
        y = math.sqrt(1-z**2) * math.sin(theta)
        dirvec = ROOT.TVector3(x, y, z)
        #guarantee it is a unit vector
        dirvec = dirvec.Unit()
        return dirvec

    def _gen_position(self):
        #generate position
        r = self._rand.Uniform(self._options.rmin, self._options.rmax)
        theta = self._rand.Uniform(0.0, 2.0*math.pi)
        z = self._rand.Uniform(self._options.zmin, self._options.zmax)
        x = r*math.cos(theta)
        y = r*math.sin(theta)
        t = 0.0
        xvec = ROOT.TLorentzVector(x, y, z, t)
        return xvec



def parsecml():
    parser = argparse.ArgumentParser()
    parser.add_argument("outfilename", help="Output file name.", type=str)
    parser.add_argument("-n", default=100, help="number of events", type=int)
    parser.add_argument("--particle", choices=ParticleType.ALL, help="choose what particle type to generate.", default=ParticleType.MU)
    parser.add_argument("--pmin", default=0.0, help="minimum particle momentum in GeV/c.", type=float)
    parser.add_argument("--pmax", default=1000.0, help="maximum particle momentum in GeV/c.", type=float)
    parser.add_argument("--rmin", default=0.0, help="minimum particle radius in metres.", type=float)
    parser.add_argument("--rmax", default=6.0, help="maximum particle radius in metres.", type=float)
    parser.add_argument("--zmin", default=0.0, help="minimum particle z in metres.", type=float)
    parser.add_argument("--zmax", default=5.5, help="maximum particle z in metres.", type=float)
    parser.add_argument("--seed", default=23213, help="Random seed.")
    return parser.parse_args()


def main():
    options = parsecml()
    job = ParticleGunJob(options)
    job.run()
    return

if __name__ == "__main__":
    main()


