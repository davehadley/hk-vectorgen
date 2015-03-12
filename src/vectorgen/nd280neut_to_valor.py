#standard library imports
from argparse import ArgumentParser
#external library imports
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
#imports from this package
import ntuple
from progress import printprogress

################################################################################

class Mode:
    ELECTRON = "electron"
    MUON = "muon"

    ALL = [ELECTRON, MUON]

class Polarity:
    FHC = "FHC"
    RHC = "RHC"
    
    ALL = [FHC, RHC]

################################################################################

def parsecml():
    parser = ArgumentParser()
    parser.add_argument("infile", type=str, help="Input genev file name.")
    parser.add_argument("-o", "--outfile", dest="outfile", type=str, help="Output file name.", default=None)
    parser.add_argument("-m", dest="mode", type=str, choices=Mode.ALL, default=Mode.ELECTRON)
    parser.add_argument("-p", "--polarity", dest="polarity", type=str, choices=Mode.ALL, default=None)
    parser.add_argument("-a", "--all", dest="runall", action="store_true", default=False)
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
    def __init__(self, tree, mode):
        self._tree = tree
        self._mode = mode
        #Pdg, Mode, Enu, EnuReco
        self.Pdg = ntuple.BranchPrimitive("Pdg", tree, 0.0)
        self.Mode = ntuple.BranchPrimitive("Mode", tree, 0.0)
        self.Enu = ntuple.BranchPrimitive("Enu", tree, 0.0)
        self.EnuReco = ntuple.BranchPrimitive("EnuReco", tree, 0.0)
        self.EnuRecoElectron = ntuple.BranchPrimitive("EnuRecoElectron", tree, 0.0)
        self.EnuRecoMuon = ntuple.BranchPrimitive("EnuRecoMuon", tree, 0.0)
        self.Weight = ntuple.BranchPrimitive("Weight", tree, 0.0)
        self.SelectionWeight = ntuple.BranchPrimitive("SelectionWeight", tree, 0.0)
        self.SelectionWeight1RMu = ntuple.BranchPrimitive("SelectionWeight1RMu", tree, 0.0)
        self.SelectionWeight1Re = ntuple.BranchPrimitive("SelectionWeight1Re", tree, 0.0)
        self.BeamPolarity = ntuple.BranchPrimitive("BeamPolarity", tree, 0.0)
        self.PotWeight = ntuple.BranchPrimitive("PotWeight", tree, 0.0)

    def Fill(self):
        return self._tree.Fill()

    def Write(self):
        return self._tree.Write()

    def SetValues(self, Pdg, ReactionCode, Enu, EnuRecoMuon, EnuRecoElectron, Weight1RMu, Weight1Re, BeamPolarity, PotWeight):
        self.Pdg.setvalue(Pdg)
        self.Mode.setvalue(ReactionCode)
        self.Enu.setvalue(Enu)
        if self._mode == Mode.ELECTRON:
            EnuReco = EnuRecoElectron
        else:
            EnuReco = EnuRecoMuon
        self.EnuReco.setvalue(EnuReco)
        self.EnuRecoElectron.setvalue(EnuRecoElectron)
        self.EnuRecoMuon.setvalue(EnuRecoMuon)
        if self._mode == Mode.ELECTRON:
            SelectionWeight = Weight1Re
        else:
            SelectionWeight = Weight1RMu
        self.Weight.setvalue(SelectionWeight*PotWeight)
        self.SelectionWeight.setvalue(SelectionWeight)
        self.SelectionWeight1RMu.setvalue(Weight1RMu)
        self.SelectionWeight1Re.setvalue(Weight1Re)
        self.BeamPolarity.setvalue(BeamPolarity)
        self.PotWeight.setvalue(PotWeight)
        return


################################################################################

def _setup_output_tree(filename, mode):
    tfile = ROOT.TFile(filename, "RECREATE")
    print filename
    tree = ROOT.TTree("selection_tree", "selection_tree")
    tree = OutputTree(tree, mode)
    return tfile, tree

################################################################################

def _copy_tree_variables(event, outtree):
    outtree.SetValues(Pdg=event.true_nu_pdg,
                      ReactionCode=event.true_reaction_code,
                      Enu=event.true_nu_e,
                      #EnuRecoMuon=event.reco_enuqe_muon,
                      #EnuRecoElectron=event.reco_enuqe_electron,
                      EnuRecoMuon=event.smeared_muonring_enuqe,
                      EnuRecoElectron=event.smeared_electronring_enuqe,
                      Weight1RMu=event.weight_selection_muon,
                      Weight1Re=event.weight_selection_electron,
                      BeamPolarity=event.polarity,
                      PotWeight=event.potscale,
                  )
    return

################################################################################

def convert_nd280neut_to_valor(inputfile, outputfile, mode, polarity):
    infile, intree = _load_input_tree(inputfile)
    outfile, outtree = _setup_output_tree(outputfile, mode)
    polarity = {None:None, Polarity.FHC:0, Polarity.RHC:1}[polarity]
    iterevents = printprogress("create_valor", intree.GetEntries(), intree)
    for event in iterevents:
        if polarity is None or polarity == event.polarity:
            if event.experiment == 0:
                _copy_tree_variables(event, outtree)
                outtree.Fill()
    outtree.Write()
    outfile.Close()
    return

################################################################################

def run_all(cml):
    inname, outname = getfilenames(cml)
    jobs = []
    for polarity in Polarity.ALL:
        for mode in Mode.ALL:
            postfix = "_" + "_".join((polarity, mode))
            name = outname.replace(".root", postfix + ".root")
            if name == outname:
                raise Exception("cannot automatically determine output filename")
            def callable(inname=inname, outname=name, mode=mode, polarity=polarity):
                return convert_nd280neut_to_valor(inname, outname, mode, polarity)
            jobs.append((name, callable))
    for n, j in jobs:
        print "processing", n
        j()
    return

################################################################################

def main():
    cml = parsecml()
    if cml.runall:
        run_all(cml)
    else:
        inname, outname = getfilenames(cml)
        convert_nd280neut_to_valor(inname, outname, cml.mode, cml.polarity)
    return

if __name__ == "__main__":
    main()
