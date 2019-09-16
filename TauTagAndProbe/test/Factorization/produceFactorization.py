from ROOT import *
from array import array
gROOT.SetBatch(True)
from math import sqrt
from histMaps import histMap
import numpy as np
from helpers import progress
from array import array

#Parse arguments
import argparse

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--subJob',              action='store',         default=0)
argParser.add_argument('--isTest',              action='store',         default=False)

args = argParser.parse_args()

#Load in samples
import Sample
sampleList = Sample.createSampleList('/storage_mnt/storage/user/lwezenbe/CMSSW_10_2_13/src/TauTriggerTools/TauTagAndProbe/test/Factorization/data/inputFiles_2018.conf')
sample = Sample.getSampleFromList(sampleList, 'noTagAndProbe')

in_file = TFile.Open(sample.path)
tree = sample.initTree()
#sample.Chain = file.Get('Ntuplizer_noTagAndProbe_multipleTaus/TagAndProbe')
triggerNamesTree = in_file.Get("Ntuplizer_noTagAndProbe_multipleTaus/triggerNames")	
if args.isTest:
    eventRange=xrange(200)
else:
    eventRange = sample.getEventRange(int(args.subJob))


gStyle.SetFrameLineWidth(1)
gStyle.SetPadBottomMargin(0.13)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)

#WPs = ["vvlooseTauMVA", "vlooseTauMVA", "looseTauMVA", "mediumTauMVA", "tightTauMVA", "vtightTauMVA", "vvtightTauMVA"]
WPs = ["vlooseTauMVA", "looseTauMVA", "mediumTauMVA", "tightTauMVA", "vtightTauMVA"]
tauDMs = ["splitPerDM"]

isDMspesific = True

outputname = "/storage_mnt/storage/user/lwezenbe/CMSSW_10_2_13/src/TauTriggerTools/TauTagAndProbe/test/Factorization/data/tauTriggerFactorization2018_"+str(args.subJob)+".root"


# get binning from binning2017.py and binning2018.py files.
from binning2018 import getbinning2018
bins = getbinning2018()
	
bin2D = bins.getBinning()["ditau"]
bin2DDM = bins.getBinning()["ditau"]
bin1D = np.linspace(0, 200, 50)
bin1DDM = np.linspace(0, 200, 50)
	
ditauTriggeredHists_1D = histMap("ditauTriggered_mass", WPs, tauDMs, bin1D, bin1DDM)
weightedDitauHists_1D = histMap("weightedDitau_mass", WPs, tauDMs, bin1D, bin1DDM)
denomHists_1D = histMap('denomHists_mass', WPs, tauDMs, bin1D, bin1DDM)

ditauTriggeredHists_2D = histMap("ditauTriggered_2D", WPs, tauDMs, array('f',bin2D), bin2DDM, isTH2=True)
weightedDitauHists_2D = histMap("weightedDitau_2D", WPs, tauDMs, array('f',bin2D), bin2DDM, isTH2=True)
denomHists_2D = histMap('denomHists_2D', WPs, tauDMs, array('f',bin2D), bin2DDM, isTH2=True)

eff_file = TFile.Open('/storage_mnt/storage/user/lwezenbe/CMSSW_10_2_13/src/TauTriggerTools/TauTagAndProbe/test/Factorization/tauTriggerEfficiencies2018.root')

print "Populating histograms"

Nevts = 0
for iEv in eventRange:
    tree.GetEntry(iEv)
    progress(iEv - eventRange[0], len(eventRange))       
    tauPt = tree.tauPt
    HLTPt = tree.hltPt
    tauEta = tree.tauEta
    tauPhi = tree.tauPhi
    tauE = tree.tauE
    tauDM = tree.tauDecayMode
    Nvtx = tree.Nvtx
    
    # tau Energy Shift (SF) is only applied for 2017 for now!
    if("DYJets" in in_file):
        puweight = tree.puweight
        tauPt_ESshifted = tree.tauPt  
    else:
        puweight = 1.
        tauPt_ESshifted = tree.tauPt 

    #vvlooseWP = tree.byVVLooseIsolationMVArun2v1DBoldDMwLT
    vlooseWP = tree.byVLooseIsolationMVArun2v1DBoldDMwLT
    looseWP = tree.byLooseIsolationMVArun2v1DBoldDMwLT	
    mediumWP = tree.byMediumIsolationMVArun2v1DBoldDMwLT	
    tightWP = tree.byTightIsolationMVArun2v1DBoldDMwLT
    vtightWP = tree.byVTightIsolationMVArun2v1DBoldDMwLT

    p1 = TLorentzVector()
    p1.SetPtEtaPhiE(tauPt[0], tauEta[0], tauPhi[0], tauE[0])
    p2 = TLorentzVector()
    p2.SetPtEtaPhiE(tauPt[1], tauEta[1], tauPhi[1], tauE[1])
    mass = (p1+p2).M()
    pt_values = sorted([p1.Pt(), p2.Pt()], reverse=True)

    passed_ditau_trigger = tree.hasHLTPath_18

    Nevents = tree.EventNumber
    Nevts =Nevts + 1
    
    bkgSubW = 1.
    weight = bkgSubW*puweight

    WPoints = {"vvlooseTauMVA":"vvlooseMVAv2", "vlooseTauMVA":"vlooseMVAv2", "looseTauMVA":"looseMVAv2", "mediumTauMVA":"mediumMVAv2", "tightTauMVA":"tightMVAv2", "vtightTauMVA":"vtightMVAv2", "vvtightTauMVA":"vvtightMVAv2"}
    DMmap = {0:'dm0', 1:'dm1', 10:'dm10'}

    # Filling the histograms
    for WP in WPs:

        eff_curve_leg1 = eff_file.Get('ditau_'+WPoints[WP]+"_"+DMmap[tauDM[0]]+'_MC_fit')
        eff_curve_leg2 = eff_file.Get('ditau_'+WPoints[WP]+"_"+DMmap[tauDM[1]]+'_MC_fit')

        efficiency_leg1_DM = eff_curve_leg1.Eval(tauPt[0])
        efficiency_leg2_DM = eff_curve_leg2.Eval(tauPt[1])
        weightedDitauHists_1D.fillHist(WP, 'splitPerDM', mass, weight*efficiency_leg1_DM*efficiency_leg2_DM)
        weightedDitauHists_2D.fillHist(WP, 'splitPerDM', pt_values, weight*efficiency_leg1_DM*efficiency_leg2_DM)
        denomHists_1D.fillHist(WP, 'splitPerDM', mass, weight)
        denomHists_2D.fillHist(WP, 'splitPerDM', pt_values, weight)
        

        if(passed_ditau_trigger):
            ditauTriggeredHists_1D.fillHist(WP, 'splitPerDM', mass, weight)
            ditauTriggeredHists_2D.fillHist(WP, 'splitPerDM', pt_values, weight)

if not args.isTest: 
    out_file = TFile( outputname, 'recreate')

    for WP in WPs:
        ditauTriggeredHists_1D.returnHist(WP, 'splitPerDM').Write(ditauTriggeredHists_1D.name + "_" +WP + "_inclusive")
        ditauTriggeredHists_2D.returnHist(WP, 'splitPerDM').Write(ditauTriggeredHists_2D.name + "_" +WP + "_inclusive")
        weightedDitauHists_1D.returnHist(WP, 'splitPerDM').Write(weightedDitauHists_1D.name + "_" +WP + "_inclusive")
        weightedDitauHists_2D.returnHist(WP, 'splitPerDM').Write(weightedDitauHists_2D.name + "_" +WP + "_inclusive")
        denomHists_1D.returnHist(WP, 'splitPerDM').Write(denomHists_1D.name + "_" +WP + "_inclusive")
        denomHists_2D.returnHist(WP, 'splitPerDM').Write(denomHists_2D.name + "_" +WP + "_inclusive")

    out_file.Close()
    print "The output ROOT file has been created: " + outputname
