import ROOT, glob, os
from helpers_old import isValidRootFile, getObjFromFile

class Sample:

    def __init__(self, name, path, output, splitJobs, xsec, skipFaultyFiles = False):           #Added last variable due to t2b pnfs problems that took too long to get solved
        self.name               = name
        self.path               = path
        self.isData             = (xsec == 'data')
        self.xsec               = eval(xsec) if not self.isData else None
        self.splitJobs          = splitJobs
        self.hCount             = None
        self.Chain              = None
        self.output             = output
        self.singleFile         = self.path.endswith('.root')
        self.isSkimmed          = not 'pnfs' in self.path
 
    def initTree(self, branches = None, needhCount=True):
        self.Chain              = ROOT.TChain('tOut')
        if self.singleFile:
            listOfFiles         = [self.path]
        else:
            listOfFiles         = sorted(glob.glob(self.path + '/*.root'))
        
        for f in listOfFiles:
            if 'pnfs' in f:
                f = 'root://maite.iihe.ac.be'+f
            self.Chain.Add(f)                                                   
        
        return self.Chain

    def getEventRange(self, subJob):
        limits = [entry*self.Chain.GetEntries()/self.splitJobs for entry in range(self.splitJobs)] + [self.Chain.GetEntries()]
        return xrange(limits[subJob], limits[subJob+1])

def createSampleList(fileName):
    sampleInfos = [line.split('%')[0].strip() for line in open(fileName)]                     # Strip % comments and \n charachters
    sampleInfos = [line.split() for line in sampleInfos if line]                              # Get lines into tuples
    for name, path, output, splitJobs, xsec in sampleInfos:
        yield Sample(name, path, output, int(splitJobs), xsec)

def getSampleFromList(sampleList, name):
  return next((s for s in sampleList if s.name==name), None)
     
