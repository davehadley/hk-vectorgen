import itertools
import ROOT
import numpy
import operator
import os
import glob
import progress
import sys

###############################################################################

class BranchPrimitive:
    def __init__(self, name, tree, start=None):
        """
        :param name: name of the output variable.
        :type name: str
        :param tree: the output TTree.
        """
        self.name = name
        if start is None:
            start = 0.0
        self.value = numpy.array([start], dtype=float)
        self.start = start
        self.tree = tree
        self._setbranch()
        
    def _setbranch(self):
        branch = self.tree.FindBranch(self.name)
        if not branch:
            branch = self.tree.Branch(self.name, self.value, self.name+"/D" )
        else:
            #self.tree.SetBranchAddress(self.name, self.value, self.name+"/D" )
            self.tree.SetBranchAddress(self.name, self.value)
        self.branch = branch
    
    def _clearbranch(self):
        if self.branch:
            self.branch.SetAddress(0)
    
    def setvalue(self, value):
        """Set the value. This should be done prior to calling TTree.Fill().
        
        :param value: value to be set for current entry in TTree.
        :type value: float 
        """
        self.value[0] = value
    
    def getvalue(self):
        """
        :returns: current value.
        """
        return self.value[0]
    
    def reset(self):
        self.setvalue(self.start)
    
    def __str__(self):
        return "(",self.name,self.value[0],")"
    
###############################################################################
    
class BranchObject:
    """Writes C++ objects to a ROOT TTree. 
    
    The object must have a default constructor that takes no arguments.
    The dictionaries for this object must be loaded.
    """
    def __init__(self,ObjectType, name, tree, start=None):
        """
        :param ObjectType: the class of this type eg ROOT.TLorentzVector
        :type ObjectType: class
        :param name: name of the output branch.
        :type name: str
        :param tree: the output TTree.
        :type tree: ROOT.TTree
        """
        self.name = name
        self.ObjectType = ObjectType
        if start is None:
            start = ObjectType() 
        self.value = start
        self.start = start
        self.tree = tree
        self._setbranch()
        
    def _setbranch(self):
        self.branch = self.tree.FindBranch(self.name)
        if not self.branch:
            self.branch = self.tree.Branch(self.name, self.value)
        else:
            #tree.SetBranchAddress(self.name, self.value)
            self.branch.SetAddress(ROOT.AddressOf(self.value))
    
    def _clearbranch(self):
        if self.branch:
            self.branch.SetAddress(0)

    def setvalue(self,value):
        '''Set reference to the output object.
        
        Takes reference of the object to be writen. 
        Preferred method provided value is guaranteed to exist at time of ROOT.TTree.Fill.
        
        :param value: reference to the output object to be filled.
        '''
        self.value = value
        #update branch address
        self.branch.SetAddress(ROOT.AddressOf(self.value))
    
    def copyvalue(self,value):
        '''Tries to copy the value. 
        
        Use in case value is not guaranteed to exist at the time of ROOT.TTree.Fill.
        This is potentially a slow method so set should be used if possible.
        
        :param value: reference to the output object to be copied then filled.
        '''
        self.value = self.ObjectType(value)
        #update branch address
        self.branch.SetAddress(ROOT.AddressOf(self.value))
        
    def reset(self):
        self.setvalue(self.start)
    
    def getvalue(self):
        """
        :returns: reference to the current object.
        """
        return self.value
    
    def __str__(self):
        return "BranchObject("+self.name+")"

###############################################################################

class TreeCopier:
    
    def __init__(self,treeName,inputDataset,outputFileName, outputFile=None):
        self.treeName = treeName
        self.inputDataset = inputDataset
        self.outputFileName = outputFileName
        #create input chain
        self.inputChain = ROOT.TChain(self.treeName,self.treeName)
        for fileName in self.inputDataset:
            self.inputChain.AddFile(fileName)
        #create output file
        if outputFileName:
            self.outputFile = ROOT.TFile(self.outputFileName,"recreate")
        else:
            outputFile.cd()
            self.outputFile = outputFile
        assert self.outputFile.IsOpen()
        self.outputTree = self.inputChain.CloneTree(0)
        self.outputFile.cd()
        
    def __iter__(self):
        for event in self.inputChain:
            yield event
            
    def getOutputTree(self):
        return self.outputTree

    def fill(self):
        self.outputTree.Fill()
    
    def close(self):
        self.outputTree.AutoSave()
        self.outputFile.Close()
        
    def numEntries(self):
        return self.inputChain.GetEntries()

###############################################################################

class Algorithm(object):
    '''Define the interface for algorithms applied to ROOT tree by the ProcessTree runnable.'''
    
    def begin(self):
        '''Called once at the beginning of process.'''
        pass
    
    def file(self, tfile):
        '''Called once each time a new file is loaded.'''
        pass
    
    def event(self, event):
        '''Called once each time a new event is loaded.'''
        pass
    
    def end(self):
        '''Called once at the end of processing.'''
        pass

###############################################################################

class AlgorithmList(Algorithm):
    
    class StopIteration(Exception):
        pass
    
    def __init__(self, algs):
        self._algs = algs
    
    def begin(self):
        for a in self._algs:
            a.begin()
    
    def file(self, tfile):
        for a in self._algs:
            a.file(tfile)
    
    def event(self, event):
        try:
            for a in self._algs:
                a.event(event)
        except AlgorithmList.StopIteration:
            #algorithms can throw this exception to stop iterating
            pass
    
    def end(self):
        for a in self._algs:
            a.end()

###############################################################################

class ClearCaches(Algorithm):
    def __init__(self, caches):
        self._caches = list(caches)
        for c in caches:
            if not callable(getattr(c, "clear", default=None)):
                raise Exception("ClearCaches requires that all input arguments have a method \"clear\".", c)
    
    def event(self, event):
        for c in self._caches:
            c.clear()


###############################################################################

class TreeFillerAlgorithm(Algorithm):
    
    def __init__(self, outfilename, treename, branches, append=False):
        self._outfilename = outfilename
        self._outtreename = treename
        self._branches = branches
        self._outfile = None
        self._outtree = None
        self._append = append

    def begin(self):
        if self._append:
            self._outfile = ROOT.TFile(self._outfilename, "UPDATE")
            self._outtree = self._outfile.Get(self._outtreename)
            if not self._outtree:
                self._outtree = ROOT.TTree(self._outtreename, self._outtreename)
        else:
            self._outfile = ROOT.TFile(self._outfilename, "RECREATE")
            self._outtree = ROOT.TTree(self._outtreename, self._outtreename)
        for b in sorted(self._branches, key=operator.attrgetter("name")):
            b.createbranch(self._outtree)
        return

    def file(self, tfile):
        return

    def event(self, event):
        for br in self._branches:
            br.event(event)
        self._outtree.Fill()
        return

    def end(self):
        self._outtree.AutoSave()
        self._outfile.Close()
        return

###############################################################################

class FileObjectGetter:
    def __init__(self, name):
        self._name = name
    def __call__(self, tfile):
        return tfile.Get(self._name)


###############################################################################

class FriendTreeGetter:
    def __init__(self, treename, friendtrees):
        self._treename = treename
        self._friendtrees = friendtrees
        self._friendlist = []

    def __call__(self, tfile):
        tree = tfile.Get(self._treename)
        self._friendlist = []
        nentries = tree.GetEntries()
        for fn in self._friendtrees:
            friend = tfile.Get(fn)
            if not friend:
                raise Exception("no tree called {}".format(fn))
            if not nentries == friend.GetEntries():
                raise Exception("trees have a different number of entries", self._treename, fn)
            tree.AddFriend(friend, fn)
            self._friendlist.append(friend)
        return tree

###############################################################################

class FileObjectTupleGetter:
    def __init__(self, *args):
        self._names = args
    def __call__(self, tfile):
        return tuple(tfile.Get(name) for name in self._names)

###############################################################################

class ProcessTree(object):
    def __init__(self, infilelist, treename, alg, n_max=None, usefasttree=True, treegetter=None):
        '''Iterates over the input tree and applies the provided algorithm.
        '''
        self._filelist = _expand_file_patterns(infilelist)
        self._alg = alg
        self._treename = treename
        if treegetter is None:
            treegetter = FileObjectGetter(treename)
        self._tree_getter = treegetter
        self._n_max = n_max
        self._usefasttree = usefasttree
        return
    
    def run(self):
        #calculate the total number of events to be processed
        self._total_num_events()
        iterable = progress.printprogress("ProcessTree "+self._treename, 
                                               self._total_num_events(), 
                                               self._iterable(), 
                                               update=True
                                        )
        #only process the required number of events
        n_max = self._n_max
        if n_max is not None:
            if n_max > 0:
                for i, val in enumerate(iterable):
                    #exit if we have done enough events
                    if i >= n_max:
                        break
        else:
            for i in iterable:
                #do every event
                pass
        return
    
    def _iterable(self):
        self._alg.begin()
        for tfile in _iter_root_files(self._filelist):
            self._alg.file(tfile)
            tree = self._tree_getter(tfile)
            if self._usefasttree:
                tree = FastTree(tree)
            for event in tree:
                yield self._alg.event(event)
        self._alg.end()
        return
    
    def _total_num_events(self):
        nevents = 0
        for tfile in _iter_root_files(self._filelist):
            tree = self._tree_getter(tfile)
            nevents += tree.GetEntries()
        return nevents

###############################################################################

class ProcessParallelTrees(ProcessTree):
    def __init__(self, infilelist, treenamelist, alg, n_max=None, usefasttree=True):
        '''Iterates over several input trees within the same file.
        '''
        treename = treenamelist[0]
        super(ProcessParallelTrees, self).__init__(infilelist, treename, alg, n_max=n_max, usefasttree=usefasttree)
        self._treenamelist = treenamelist
        self._tree_getter = FileObjectTupleGetter(*treenamelist)
    
    def _iterable(self):
        self._alg.begin()
        for tfile in _iter_root_files(self._filelist):
            self._alg.file(tfile)
            treetuple = self._tree_getter(tfile)
            if self._usefasttree:
                treetuple = tuple(FastTree(tree) for tree in treetuple)
            for event in itertools.izip(*treetuple):
                yield self._alg.event(event)
        self._alg.end()
        return

    def _total_num_events(self):
        nevents = 0
        for tfile in _iter_root_files(self._filelist):
            tree = self._tree_getter(tfile)[0]
            nevents += tree.GetEntries()
        return nevents

###############################################################################

class ProcessFriendTrees(ProcessTree):
    def __init__(self, inputfiles, alg, n_max=None, usefasttree=True):
        '''Iterates over several input trees in different files.
        '''
        self._filelists = [(tn, flist) for tn, flist in inputfiles]
        treename = self._filelists[0][0]
        infilelist = self._filelists[0][1]
        super(ProcessFriendTrees, self).__init__(infilelist, treename, alg, n_max=n_max, usefasttree=usefasttree)
        self._multi_tree_getter = dict([(tn, FileObjectGetter(tn)) for tn, flist in self._filelists])
    
    def _chain(self, filelist, treename):
        for tfile in _iter_root_files(filelist):
            self._alg.file(tfile)
            tree = self._multi_tree_getter[treename](tfile)
            if self._usefasttree:
                tree = FastTree(tree)
            for event in tree:
                yield event
    
    def _iterable(self):
        self._alg.begin()
        chains = [self._chain(flist, tn) for tn, flist in self._filelists]
        for event in itertools.izip(*chains):
            yield self._alg.event(event)
        self._alg.end()
        return

###############################################################################

class FastTree(object):
    '''Lazilly enable branches as they are used.
    This prevents loading unused data. This will be slower than a regular tree
     in cases where most or all branches are used.
    '''   
    def __init__(self, tree, enablebranches=[]):
        self._tree = tree
        self._tree.SetBranchStatus("*", 0)
        self._encountered = set()
        for name in enablebranches:
            self._encountered.add(name)
            self._enablebranch(name)

    def _enablebranch(self, name):
        self._tree.SetBranchStatus(name, 1)
        #self._tree.ResetBranchAddresses()
        return 

    def __getattr__(self, name):
        if name in self._encountered:
            return getattr(self._tree, name)
        else:
            self._encountered.add(name)
            if self._tree.FindBranch(name):
                self._enablebranch(name)
                self._tree.SetBranchStatus(name, 1)
            return getattr(self._tree, name)

    def __iter__(self):
        i = 0
        tree = self._tree
        n = tree.GetEntries()
        for i in xrange(n):
            tree.GetEntry(i)
            yield self

###############################################################################

class ProcessTreeSubset(ProcessTree):
    def __init__(self, infilelist, treename, alg, cutstr, n_max=None, usefasttree=True, treegetter=None):
        '''Iterates over the input tree and applies the provided algorithm.
        Only events that pass the cut will be processed.
        '''
        super(ProcessTreeSubset, self).__init__(infilelist, treename, alg, n_max=n_max, usefasttree=usefasttree, treegetter=treegetter)
        self._cutstr = cutstr
        return

    def run(self):
        self._alg.begin()
        if len(self._filelist)>1:
            #more than one file, progress display will show an entry for each event.
            iterable = progress.printprogress("ProcessTreeSubset "+self._treename, 
                                               len(self._filelist), 
                                               self._filelist, 
                                               update=True
                                        )
        else:
            #special case for just 1 input file, progress display will show per event information
            iterable = self._filelist
        n_max = self._n_max
        if n_max is None:
            n_max = float("inf")
        count = 0
        self._alg.begin()
        for tfile in _iter_root_files(iterable):
            self._alg.file(tfile)
            tree = self._tree_getter(tfile)
            eventlist = self._create_eventlist(tree)
            if self._usefasttree:
                #use after create event list, as the Draw command may fail otherwise.
                tree = FastTree(tree)
            eventlistiterable = xrange(eventlist.GetN()) 
            if len(self._filelist)<=1:
                eventlistiterable = progress.printprogress("ProcessTreeSubset "+self._treename, 
                                               eventlist.GetN(), 
                                               eventlistiterable, 
                                               update=True
                                        )
            for index_el in eventlistiterable:
                eventnum = eventlist.GetEntry(index_el)
                if eventnum > -1: # -1 is a failure code for TEventList::GetEntry
                    if count >= n_max:
                        break
                    count += 1
                    self._load_event_from_ttree(tree, eventnum)
                    self._alg.event(tree)
            if count >= n_max:
                break
        self._alg.end()
        return
    
    def _load_event_from_ttree(self, tree, eventnum):
        #separating this out into a separate function makes it more 
        #obvious where time is being spent when profiling code.
        return tree.GetEntry(eventnum)
    
    def _create_eventlist(self, tree):
        ret = tree.Draw(">>elist", self._cutstr)
        if ret < 0:
            raise Exception("ProcessTreeSubset failed to get an event list. This is usually because the cut string is invalid.", self._cutstr)
        eventlist = ROOT.gDirectory.Get("elist")
        return eventlist

###############################################################################

class BranchFiller(object):
    def __init__(self, name, function, start_value=0.0, ignore_errors=False):
        '''BranchFiller handles setting branch values on each event and is designed to be provided to a TreeFillerAlgorithm.
        The required arguments are the name of the output branch and a callable object 
        or function that is applied to the event that returns a double.
        Alternatively, function if is a string, a simple attribute getter is constructed with it. 
        Optionally the start_value can be set.
        '''
        self.name = name
        self._start_value = start_value
        self._branch = None
        self._ignore_errors = ignore_errors
        if isinstance(function, str):
            function = operator.attrgetter(function)
        self._function = function
        if not callable(self._function):
            raise Exception("BranchFiller given a non-callable function object.")
    
    def createbranch(self, tree):
        self._branch = BranchPrimitive(name=self.name, tree=tree, start=self._start_value)
        return
    
    def event(self, event):
        try:
            val = self._eval(event)
        except Exception as ex:
            #use default value
            val = self._branch.start
            if not self._ignore_errors:
                self._handle_exception("BranchFiller failed to evaluate function", ex)
        #Sometimes setting fails due to a bug in the above function.
        #for example if the wrong type is returned.
        try:
            self._branch.setvalue(val)
        except Exception as ex:
            if not self._ignore_errors:
                self._handle_exception("BranchFiller failed to set branch value", ex, val)
            
    
    def _handle_exception(self, msg, ex, val=None):
        #get the exception information
        exc_info = sys.exc_info()
        #repackage the 
        if val is None:
            exinst = exc_info[0](msg, 
                             self.name,
                             exc_info[1],
                             )
        else:
            exinst = exc_info[0](msg, 
                             self.name, 
                             val, 
                             type(val), 
                             exc_info[1],
                             )
        #re-raise the re-packaged exception
        raise exinst, None, exc_info[2]
    
    def _eval(self, event):
        return self._function(event)

###############################################################################

class BranchObjectFiller(BranchFiller):
    def __init__(self, name, function, cls, start_value=None, ignore_errors=False):
        super(BranchObjectFiller, self).__init__(name, function, start_value=start_value, ignore_errors=ignore_errors)
        self._cls = cls

    def createbranch(self, tree):
        self._branch = BranchObject(ObjectType=self._cls, name=self.name, tree=tree, start=self._start_value)
        return

###############################################################################

def _find_file_matches(pattern):
    pattern = os.path.expanduser(os.path.expandvars(pattern))
    matches = glob.glob(pattern)
    return matches

###############################################################################

def _expand_file_patterns(pattern_list):
    result = []
    for p in pattern_list:
        result.extend(_find_file_matches(p))
    return result

###############################################################################

def _iter_root_files(filelist):
    for fname in filelist:
        yield ROOT.TFile(fname)

###############################################################################

_testfilename = "unitTestBranchObj.root"
_testtreename = "unitTestBranchObj"

def _unitTestBranchObjWrite():
    print "--- _unitTestBranchObjWrite"
    #generate dictionaries
    ROOT.gInterpreter.GenerateDictionary("std::vector<TLorentzVector>","TLorentzVector.h")
    #create output tree
    tfile = ROOT.TFile(_testfilename,"recreate")
    tree = ROOT.TTree("unitTestBranchObj","unitTestBranchObj")
    #define branches
    x = BranchPrimitive("x", tree)
    y = BranchPrimitive("y", tree)
    z = BranchPrimitive("z", tree)
    testLorentzVec = BranchObject(ROOT.TLorentzVector,"testLorentzVec", tree)
    testStdVec = BranchObject(ROOT.std.vector("double"),"testStdVec", tree)
    #fill tree with some random data
    StdVecHlv = ROOT.std.vector("TLorentzVector")
    testStdVecLorentz = BranchObject(StdVecHlv,"testStdVecLorentz", tree)
    rand = ROOT.TRandom3()
    print "unitTestBranchObj begining event loop"
    for i in xrange(1000):
        x.setvalue(rand.Gaus(1,1))
        y.setvalue(rand.Gaus(10,1))
        z.setvalue(rand.Gaus(20,1))
        e = rand.Gaus(1000,100)
        hlv = ROOT.TLorentzVector(1.0,2.0,3.0,4.0)
        testLorentzVec.setvalue(hlv)
        vec = ROOT.std.vector("double")()
        for i in xrange(10):
            vec.push_back(rand.Gaus(-1,1))
        testStdVec.setvalue(vec)
        vecHlv = StdVecHlv()
        hlv1 = ROOT.TLorentzVector(1.0,2.0,3.0,4.0)
        hlv2 = ROOT.TLorentzVector(5.0,6.0,7.0,8.0)
        vecHlv.push_back(hlv1)
        vecHlv.push_back(hlv2)
        vecHlv.push_back(ROOT.TLorentzVector(9.0,10.0,11.0,12.0))
        testStdVecLorentz.setvalue(vecHlv)
        tree.Fill()
    tree.Write()
    tfile.Close()
    print "--- _unitTestBranchObjWrite... done." 
    return

##################################################

def _unitTestBranchObjRead():
    print "--- _unitTestBranchObjRead"
    tfile = ROOT.TFile(_testfilename,"read")
    tree = tfile.Get(_testtreename)
    print tree
    for i,event in enumerate(tree):
        print "TEST double = ",event.x
        event.testLorentzVec.Print()
        print "TEST testStdVec = ",event.testStdVec.size(),[ x for x in event.testStdVec ]
        print "TEST vec Lorentz Vec:",event.testStdVecLorentz.size(),"entries"
        for hlv in event.testStdVecLorentz:
            hlv.Print()
        if i>=2:
            break
    print "--- _unitTestBranchObjRead ..done"
    return

##################################################

def _unitTestAlgorithm():
    #create and output file
    _unitTestBranchObjWrite()
    branches = [BranchFiller("x", "x", start_value=0.0),
                BranchFiller("x2", lambda x: x.x*2, start_value=0.0),
                ]
    treefiller = TreeFillerAlgorithm(outfilename="test_alg_ntuple.root", treename="test_alg", branches=branches)
    process = ProcessTree(infilelist=[_testfilename], treename=_testtreename, alg=treefiller)
    process.run()
    return

##################################################

def main():
    """Run some unit tests that check reading and writing of a ttree still works."""
    #_unitTestBranchObjWrite()
    #_unitTestBranchObjRead()
    _unitTestAlgorithm()
    return

if __name__ == "__main__":
    main()
