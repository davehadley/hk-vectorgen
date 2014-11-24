import ROOT
import argparse

def _parsecml():
    arg = argparse.ArgumentParser()
    arg.add_argument("nevents", type=int)
    arg.add_argument("filename", type=str)
    arg.add_argument("treename", type=str)
    arg.add_argument("-t", dest="test", action="store_true", help="run some tests.")
    return arg.parse_args()

_testfilename = "test.root"

def _testwrite():
    import numpy
    tfile = ROOT.TFile(_testfilename, "RECREATE")
    tree = ROOT.TTree("tree", "tree")
    x = numpy.array([0.0], dtype=float)
    tree.Branch("x", x, "x/D")
    for i in xrange(1000):
        x[0] = float(i)
        tree.Fill()
    tree.Write()
    tfile.Close()
    return

def _testread():
    i = 0
    prediction = 0
    while True:
        fname = "{}_{}.root".format(_testfilename.split(".", 1)[0], i)
        tfile = ROOT.TFile(fname, "READ")
        if not tfile.IsOpen():
            break
        tree = tfile.Get("tree")
        for j, entry in enumerate(tree):
            #print "DEBUG file", i, ", entry",j, ", value", entry.x, ", prediction", prediction
            if not (entry.x == prediction):
                raise Exception("entry does not match prediction", entry.x, prediction, fname)
            prediction += 1
            
        i += 1
    return

def splitfile(nevents, filename, treename):
    tfile = ROOT.TFile(filename, "READ")
    if not tfile.IsOpen():
        raise Exception("failed to open file", filename)
    tree = tfile.Get(treename)
    if not tree:
        raise Exception("failed to retrieve tree", treename)
    ninput = tree.GetEntries()
    start = 0
    i = 0
    while start < ninput:
        end = start + nevents
        outname = "{}_{}.root".format(filename.split(".", 1)[0], i)
        out = ROOT.TFile(outname, "RECREATE")
        outtree = tree.CopyTree("", "", nevents, start)
        outtree.Write()
        start = end
        i += 1
    return

def main():
    opt = _parsecml()
    if opt.test:
        #_testwrite()
        _testread()
    splitfile(opt.nevents, opt.filename, opt.treename)
    return

if __name__ == "__main__":
    main()
