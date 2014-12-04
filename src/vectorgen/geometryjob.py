import os
import ROOT

from vectorgen.constants import Orientation
from vectorgen.jobtools import IJob

###############################################################################

class CylinderGeometry:
    def __init__(self, radius=4.0, z=8.0, orientation=Orientation.Z):
        self.radius = radius
        self.z = z
        self.orientation = orientation
        self.name = self._uniquestr()

    def verify(self):
        return

    def _uniquestr(self):
        return "_".join((
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
    def __init__(self, radius=4.0, z=8.0, orientation=Orientation.Z):
        self.radius = radius
        self.z = z
        self.orientation = orientation
        self.name = self._uniquestr()

    def verify(self):
        return

    def _uniquestr(self):
        return "_".join((
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

