#standard library imports
from argparse import ArgumentParser
#external library imports
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
#imports from this package
import ntuple

################################################################################

def parsecml():
    parser = ArgumentParser()
    parser.add_argument("infile", type=str, help="Input genev file name.")
    parser.add_argument("-o", "--outfile", dest="outfile", type=str, help="Output file name.", default=None)
    return parser.parse_args()

################################################################################

def getfilenames(opt):
    inname = opt.infile
    outname = opt.outfile
    if outname is None:
        #automatically determine output file name
        if not ("genev" in inname and ".root" in inname):
            raise Exception("Cannot automatically determine outputfile name from input. You must specify it manually with the -o option.")
    outname = inname.replace("genev", "valorinput")
    return inname, outname


################################################################################

def _load_input_tree(filename):
    #Load input file and get tree
    infile = ROOT.TFile(inputfile, "READ")
    intree = infile.Get("nRootTracker") # try to get NEUT tree
    if not intree:
        intree = infile.Get("gRooTracker") # try to get GENIE tree
    return infile, intree

################################################################################

class OutputTree:
    def __init__(self, tree):
        self._tree = tree
        #Pdg, Mode, Enu, EnuReco
        Pdg = ntuple.BranchPrimitive("Pdg", tree, 0.0)
        Mode = ntuple.BranchPrimitive("Mode", tree, 0.0)
        Enu = ntuple.BranchPrimitive("Enu", tree, 0.0)
        EnuReco = ntuple.BranchPrimitive("EnuReco", tree, 0.0)
        Weight1RMu = ntuple.BranchPrimitive("Weight1RMu", tree, 0.0)
        Weight1Re = ntuple.BranchPrimitive("Weight1Re", tree, 0.0)

    def Fill(self):
        return self._tree.Fill()

    def Write(self):
        return self._tree.Write()

    def SetValues(self, Pdg, Mode, Enu, EnuReco, Weight1RMu, Weight1Re):
        self.Pdg.setvalue(Pdg)
        self.Mode.setvalue(Mode)
        self.Enu.setvalue(Enu)
        self.EnuReco.setvalue(EnuReco)
        self.Weight1RMu.setvalue(Weight1RMu)
        self.Weight1Re.setvalue(Weight1Re)
        return


################################################################################

def _setup_output_tree(filename):
    tfile = ROOT.TFile(filename, "RECREATE")
    tree = ROOT.TTree("selection_tree", "selection_tree")
    return tfile, tree

################################################################################

def _copy_tree_variables(event, outtree):
    outtree.SetValues(Pdg=,
                      Mode=,
                      Enu=,
                      EnuReco=,
                      Weight1RMu=,
                      Weight1Re=,
                  )
    return

################################################################################

def convert_nd280neut_to_valor(inputfile, outputfile):
    infile, intree = _load_input_tree(inputfilename)
    outfile, outtree = _setup_output_tree(outputfile)
    for event in intree:
        if event.detector = 0:
            _copy_tree_variables(event, outtree)
            outtree.Fill()
    outtree.Write()
    outfile.Close()
    return

################################################################################

def main():
    cml = parsecml()
    inname, outname = getfilenames(opt)
    convert_nd280neut_to_valor(inname, outname)
    return

if __name__ == "__main__":
    main()
