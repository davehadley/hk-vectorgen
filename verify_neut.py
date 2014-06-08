import argparse
import ROOT
from simplot.rootplot import drawtools, rootio, drawoptions, style

###############################################################################

_out = rootio.CanvasWriter("plots")

###############################################################################

# def loadgenie():
#     gROOT = ROOT.gROOT
#     gSystem = ROOT.gSystem
#     script_dir = gSystem.Getenv("GENIE")
#     script_dir += "/src/scripts/gcint/"
#     curr_dir = gSystem.pwd()
#     gSystem.cd(script_dir)
#     gROOT.ProcessLine(".x loadincs.C")
#     gROOT.ProcessLine(".x loadlibs.C")
#     gSystem.cd(curr_dir);
#     gROOT.ProcessLine(".x analgenie.C")
#     return

###############################################################################

def parsecml():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_file", type=str, help="select input file")
    return ap.parse_args()

###############################################################################

def gettree(fname):
    tfile = ROOT.TFile(fname, "READ")
    tree = tfile.Get("nRooTracker")
    tree._tfile = tfile # trick to keep file open while tree exists
    return tree

###############################################################################

def run(opt):
    tree = gettree(opt.input_file)
    tree.SetWeight(1.0)
    paint = drawtools.TreePainter(tree=tree)
    x = "EvtVtx[0]"
    y = "EvtVtx[1]"
    z = "EvtVtx[2]"
    e = "StdHepP4[0][3]"
    nu_px = "StdHepP4[0][0]"
    nu_py = "StdHepP4[0][1]"
    nu_pz = "StdHepP4[0][2]"
    pdg = "StdHepPdg[0]"
    pdgsplit = drawoptions.SplitDataset.from_integer_map(pdg, {"#nu_{e}" : 12,
                                                                          "#bar{#nu_{e}}" : -12,
                                                                          "#nu_{#mu}" : 14,
                                                                          "#bar{#nu_{#mu}}" : -14,
                                                                          })
    weight = drawoptions.EventWeight("1.0")
    canv = paint.paint("enu", "enu", "{e}".format(e=e), pdgsplit)
    _out.save(canv)
    canv = paint.paint("nu_px", "nu_px", "{nu_px}".format(nu_px=nu_px), pdgsplit)
    _out.save(canv)
    canv = paint.paint("nu_py", "nu_py", "{nu_py}".format(nu_py=nu_py), pdgsplit)
    _out.save(canv)
    canv = paint.paint("nu_pz", "nu_pz", "{nu_pz}".format(nu_pz=nu_pz), pdgsplit)
    _out.save(canv)
    canv = paint.paint("nupdg", "nupdg", pdg, pdgsplit, drawoptions.UniformBinning(31, -15.5, 15.5))
    _out.save(canv)
    canv = paint.paint("rnu", "rnu", "sqrt(pow({x},2)+pow({y},2))".format(x=x, y=y), weight)
    _out.save(canv)
    xbinning = drawoptions.UniformBinning(100, -6.0, 6.0)
    zbinning = drawoptions.UniformBinning(100, -15.0, 15.0)
    canv = paint.paint("xy", "xy", "{x}:{y}".format(x=x, y=y), xbinning, drawoptions.NDimensions(2), drawoptions.YBinning(xbinning), weight)
    _out.save(canv)
    for name, cmd in [("x", x),
                      ("y", y),
                      ("z", z),
                      ]:
        canv = paint.paint(name, name, cmd, zbinning, weight)
        _out.save(canv)
    return

###############################################################################


def main():
    ROOT.gROOT.SetBatch()
    style.setDefaultStyle()
    opt = parsecml()
    run(opt)

###############################################################################

if __name__ == "__main__":
    main()
