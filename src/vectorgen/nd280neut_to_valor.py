#standard library imports
from argparse import ArgumentParser
#external library imports
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
#imports from this package
import ntuple
from progress import printprogress

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
    infile = ROOT.TFile(filename, "READ")
    #intree = infile.Get("nRooTracker") # try to get NEUT tree
    #if not intree:
    #    intree = infile.Get("gRooTracker") # try to get GENIE tree
    intree = infile.Get("titus_ntuple")
    if not intree:
        raise Exception("could not load TTree from file", filename)
    return infile, intree

################################################################################

class OutputTree:
    def __init__(self, tree):
        self._tree = tree
        #Pdg, Mode, Enu, EnuReco
        self.Pdg = ntuple.BranchPrimitive("Pdg", tree, 0.0)
        self.Mode = ntuple.BranchPrimitive("Mode", tree, 0.0)
        self.Enu = ntuple.BranchPrimitive("Enu", tree, 0.0)
        self.EnuReco = ntuple.BranchPrimitive("EnuReco", tree, 0.0)
        self.Weight1RMu = ntuple.BranchPrimitive("Weight1RMu", tree, 0.0)
        self.Weight1Re = ntuple.BranchPrimitive("Weight1Re", tree, 0.0)

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
    print filename
    tree = ROOT.TTree("selection_tree", "selection_tree")
    tree = OutputTree(tree)
    return tfile, tree

################################################################################

def _copy_tree_variables(event, outtree):
    outtree.SetValues(Pdg=event.true_nu_pdg,
                      Mode=event.true_reaction_code,
                      Enu=event.true_nu_e,
                      EnuReco=event.reco_enuqe_muon,
                      Weight1RMu=event.weight_selection_muon,
                      Weight1Re=event.weight_selection_electron,
                  )
    return

################################################################################

def convert_nd280neut_to_valor(inputfile, outputfile):
    infile, intree = _load_input_tree(inputfile)
    outfile, outtree = _setup_output_tree(outputfile)
    iterevents = printprogress("create_valor", intree.GetEntries(), intree)
    for event in iterevents:
        if event.experiment == 0:
            _copy_tree_variables(event, outtree)
            outtree.Fill()
    outtree.Write()
    outfile.Close()
    return

################################################################################

def main():
    cml = parsecml()
    inname, outname = getfilenames(cml)
    convert_nd280neut_to_valor(inname, outname)
    return

if __name__ == "__main__":
    main()
