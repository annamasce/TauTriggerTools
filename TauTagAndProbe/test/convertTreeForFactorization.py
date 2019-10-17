from ROOT import *
import numpy as n
from helpers import *

MC = True
DYJets = True   # False means WJet enriched cuts will be used, True means cuts for DYJet enriched samples will be used

Sample2017 = False
Sample2018 = True

saveOnlyOS=False
if MC:
	saveOnlyOS = True # True; save only OS, False; save both and store weight for bkg sub
        if DYJets: tauGenMatching = True
	if not DYJets: tauGenMatching = False
	excludeLumiSections = False
	print "==> OS events are stored and tau gen matching is applied for MC samples! <=="
else:
	saveOnlyOS = False # True; save only OS, False; save both and store weight for bkg sub
	tauGenMatching = False
	excludeLumiSections = False
	print "==> SS events are stored as weights and applied to suppress the bkg for Data samples! <=="

disabledPScolumns = False  # True; to remove the disabled columns, False; to consider all columns

# the hadd of all the output ntuples
path = "/user/lwezenbe/CMSSW_10_2_13/src/TauTriggerTools/TauTagAndProbe/test/"

if MC:
    fname =  path + "crabTuple_v5.root"
else:
	#fname =  path + "Ntuple_SingleMuon_Run2018ABCDReReco17SepPromptRecoD_190121.root"
	#fname =  path + "Ntuple_SingleMuon_Run2018D-PromptReco-v2_190121.root"
	fname = path + "crabTuple_v5.root"
	#fname = "./NTuple_SingleMu_Data_2018_filterTest_Run201.root"



#pt = [20, 26, 30, 34]
pt = [20, 26, 30, 32, 34]
if(Sample2017):
	numberOfHLTTriggers = 23
	numberOfHLTTriggersForFit = 26
if(Sample2018):
	numberOfHLTTriggers = 21
	numberOfHLTTriggersForFit = 21



#######################################################
fIn = TFile.Open(fname)
tIn = fIn.Get('Ntuplizer_noTagAndProbe_multipleTaus/TagAndProbe')
tTriggerNames = fIn.Get("Ntuplizer_noTagAndProbe_multipleTaus/triggerNames")
if MC:
	suppressionType = "OStauGenMatched"
else:
	suppressionType = "SSsubtraction"
if DYJets:
	outname = fname.replace ('.root', '_' + suppressionType + '_VVLooseWP2017v2_forFit_260419_test.root')
else:
	outname = fname.replace ('.root', '_' + suppressionType + '_WjetEnriched_MediumWP2017v2_forFit.root')
fOut = TFile (outname, 'recreate')
tOut = tIn.CloneTree(0)
tOutNames = tTriggerNames.CloneTree(-1) # copy all

briso   = [n.zeros(1, dtype=int) for x in range (0, len(pt))]
brnoiso = [n.zeros(1, dtype=int) for x in range (0, len(pt))]
bkgSubW = n.zeros(1, dtype=float)
bkgSubANDpuW = n.zeros(1, dtype=float)

hltPathTriggered_OS   = [n.zeros(1, dtype=int) for x in range (0, numberOfHLTTriggersForFit+1)]

for i in range (0, len(pt)):
    name = ("hasL1_" + str(pt[i]))
    tOut.Branch(name, brnoiso[i], name+"/I")
    name += "_iso"
    tOut.Branch(name, briso[i], name+"/I")


for i in range (0, numberOfHLTTriggers):
	tTriggerNames.GetEntry(i)
	if(i < numberOfHLTTriggers):
		name = ("hasHLTPath_" + str(i))
	tOut.Branch(name, hltPathTriggered_OS[i], name+"/I")
	
	#tOut.Branch(name, hltPathTriggered_OS[i], name+"/I")

#tOut.Branch("isoHLT", hltPathTriggered_OS[6], name+"/I")

tOut.Branch("bkgSubW", bkgSubW, "bkgSubW/D")
tOut.Branch("bkgSubANDpuW", bkgSubANDpuW, "bkgSubANDpuW/D")

nentries = tIn.GetEntries()
for ev in range (0, nentries):
#    print len(tIn.tauTriggerBitsPerLeg), showBranch(tIn.tauTriggerBitsPerLeg), tIn.tauTriggerBitsPerLeg    
    tIn.GetEntry(ev)
    if (ev%10000 == 0) : print ev, "/", nentries

    ngoodeta = 0
    for n in xrange(tIn.ntau):
        if abs(tIn.tauEta[n]) > 2.1:
            continue
        ngoodeta+=1
    if ngoodeta < 2:    continue

    if saveOnlyOS and not tIn.isOS:
        continue

    for i in range (0, len(pt)):
        briso[i][0] = 0
        brnoiso[i][0] = 0

    for i in range (0, numberOfHLTTriggersForFit):
        hltPathTriggered_OS[i][0] = 0

    L1iso = True if tIn.l1tIso == 1 else False
    L1pt = tIn.l1tPt
    for i in range(0, len(pt)):
        # print L1pt, pt[i]
        #
        if L1pt > pt[i]:
            brnoiso[i][0] = 1
            # print "SUCCESS!! ", brnoiso[i]
            if L1iso:
                briso[i][0] = 1

    HLTpt = tIn.hltPt
    RunNumber = tIn.RunNumber
    lumi = tIn.lumi
    for bitIndex in range(0, numberOfHLTTriggers):
	import itertools as it
#    	if bitIndex in it.chain(range(3, 6), range(20, 23)):   # di-tau paths for non-HPS and HPS ones
        if ((tIn.tauTriggerBitsPerLeg[0] >> bitIndex) & 1) == 1 and ((tIn.tauTriggerBitsPerLeg[1] >> bitIndex) & 1) == 1:
                hltPathTriggered_OS[bitIndex][0] = 1
        else:
                hltPathTriggered_OS[bitIndex][0] = 0

#	if(bitIndex==14):  #mutau
#		if (((triggerBits >> bitIndex) & 1) == 1 and L1pt>=26 and L1iso>0 and HLTpt[bitIndex]>=30):
#			hltPathTriggered_OS[numberOfHLTTriggers][0] = 1	  # this is the path for etau trigger. So (L1iso) should be applied here!
#		else:
#			hltPathTriggered_OS[numberOfHLTTriggers][0] = 0
#
#		if ((triggerBits >> bitIndex) & 1) == 1:
#			hltPathTriggered_OS[numberOfHLTTriggers+1][0] = 1 # this is the path for mutau trigger. So no extra requirement is needed like: L1pt and L1iso and HLTpt
#		else:
#			hltPathTriggered_OS[numberOfHLTTriggers+1][0] = 0
#
#	if(bitIndex==8):  #mutau
#		if (((triggerBits >> bitIndex) & 1) == 1 and L1pt>=26 and L1iso>0 and HLTpt[bitIndex]>=30):
#			hltPathTriggered_OS[numberOfHLTTriggers+3][0] = 1	  # this is the path for etau trigger. So (L1iso) should be applied here!
#		else:
#			hltPathTriggered_OS[numberOfHLTTriggers+3][0] = 0
#
#		if ((triggerBits >> bitIndex) & 1) == 1:
#			hltPathTriggered_OS[numberOfHLTTriggers+4][0] = 1 # this is the path for mutau trigger. So no extra requirement is needed like: L1pt and L1iso and HLTpt
#		else:
#			hltPathTriggered_OS[numberOfHLTTriggers+4][0] = 0

#	# ditau
#    if (((((triggerBits >> 4) & 1) == 1 and HLTpt[4]>=40) or (((triggerBits >> 5) & 1) == 1 and HLTpt[5]>=40) or (((triggerBits >> 6) & 1) == 1))):
#	    hltPathTriggered_OS[numberOfHLTTriggers+5][0] = 1  # this is the path for di-tau trigger. HLTpt cut is required to have the same threshold on tau + L1Pt is needed due to L1 matching differences between MC and Data
#    else:
#	    hltPathTriggered_OS[numberOfHLTTriggers+5][0] = 0
#    if(MC):  #ditau
#            if (((triggerBits >> 15) & 1) == 1):
#                    hltPathTriggered_OS[numberOfHLTTriggers+2][0] = 1     # this is the path for ditau trigger for MC samples
#            else:
#                    hltPathTriggered_OS[numberOfHLTTriggers+2][0] = 0
#
#    elif(not MC):  #ditau
#            if (((triggerBits >> 20) & 1) == 1):
#                    hltPathTriggered_OS[numberOfHLTTriggers+2][0] = 1     # this is the path for ditau trigger for data
#            else:
#                    hltPathTriggered_OS[numberOfHLTTriggers+2][0] = 0
#
#    if(MC):
#            if (((triggerBits >> 14) & 1) == 1 and L1pt>=26 and L1iso>0 and HLTpt[14]>=30):
#                    hltPathTriggered_OS[numberOfHLTTriggers+6][0] = 1	  # this is the path for etau trigger. So (L1iso) should be applied here!
#            else:
#                    hltPathTriggered_OS[numberOfHLTTriggers+6][0] = 0
#
#            if ((triggerBits >> 14) & 1) == 1:
#                    hltPathTriggered_OS[numberOfHLTTriggers+7][0] = 1 # this is the path for mutau trigger. So no extra requirement is needed like: L1pt and L1iso and HLTpt
#            else:
#                    hltPathTriggered_OS[numberOfHLTTriggers+7][0] = 0
#
#            if (((triggerBits >> 15) & 1) == 1):
#                    hltPathTriggered_OS[numberOfHLTTriggers+8][0] = 1
#            else:
#                    hltPathTriggered_OS[numberOfHLTTriggers+8][0] = 0
#    else:
#            if(RunNumber < 317509):
#                    if (((triggerBits >> 8) & 1) == 1 and L1pt>=26 and L1iso>0 and HLTpt[8]>=30):
#                            hltPathTriggered_OS[numberOfHLTTriggers+6][0] = 1	  # this is the path for etau trigger. So (L1iso) should be applied here!
#                    else:
#                            hltPathTriggered_OS[numberOfHLTTriggers+6][0] = 0
#
#                    ie ((triggerBits >> 8) & 1) == 1:
#                            hltPathTriggered_OS[numberOfHLTTriggers+7][0] = 1 # this is the path for mutau trigger. So no extra requirement is needed like: L1pt and L1iso and HLTpt
#                    else:
#                            hltPathTriggered_OS[numberOfHLTTriggers+7][0] = 0
#
#                    if (((((triggerBits >> 4) & 1) == 1 and HLTpt[4]>=40) or (((triggerBits >> 5) & 1) == 1 and HLTpt[5]>=40) or (((triggerBits >> 6) & 1) == 1))):
#                    # this is the path for di-tau trigger. HLTpt cut is required to have the same threshold on tau + L1Pt is needed due to L1 matching differences between MC and Data
#                            hltPathTriggered_OS[numberOfHLTTriggers+8][0] = 1
#                    else:
#                            hltPathTriggered_OS[numberOfHLTTriggers+8][0] = 0
#            else:
#
#                    if (((triggerBits >> 14) & 1) == 1 and L1pt>=26 and L1iso>0 and HLTpt[14]>=30):
#                            hltPathTriggered_OS[numberOfHLTTriggers+6][0] = 1  # this is the path for etau trigger. So (L1iso) should be applied here!
#                    else:
#                            hltPathTriggered_OS[numberOfHLTTriggers+6][0] = 0
#
#                    if ((triggerBits >> 14) & 1) == 1:
#                            hltPathTriggered_OS[numberOfHLTTriggers+7][0] = 1 # this is the path for mutau trigger. So no extra requirement is needed like: L1pt and L1iso and HLTpt
#                    else:
#                            hltPathTriggered_OS[numberOfHLTTriggers+7][0] = 0
#
#                    if (((triggerBits >> 20) & 1) == 1):
#                            hltPathTriggered_OS[numberOfHLTTriggers+8][0] = 1
#                    else:
#                            hltPathTriggered_OS[numberOfHLTTriggers+8][0] = 0

		
    bkgSubW[0] = 1. if tIn.isOS else -1.

#    if not "Run201" in fname:
#        puweight = tIn.puweight
#    else:
#        puweight = 1
    puweight = 1

    bkgSubANDpuW[0] = bkgSubW[0]*puweight

    if(tIn._byVVLooseIsolationMVArun2017v2DBoldDMwLT2017[0] > 0.5 and tIn._byVVLooseIsolationMVArun2017v2DBoldDMwLT2017[1] > 0.5):

#        #Mass cuts, mt and mvis for DY Jets
#	    if DYJets:
#		    if(tIn.mT < 30 and tIn.mVis >40 and tIn.mVis < 80):
#				if(tauGenMatching):  #for tau gen matching
#				    if(tIn.tau_genindex > 0):
#					    tOut.Fill()
#				else:
#				    tOut.Fill()
#        #High mT requirement for WJets
#	    elif not DYJets:
#		    if(tIn.mT > 30):
#			    if(tauGenMatching):    #for tau gen matching
#				    if(tIn.tau_genindex > 0):
#					    tOut.Fill()
#			    else:
#				    tOut.Fill()

        tOut.Fill()

tOutNames.Write()
tOut.Write()
fOut.Close()
