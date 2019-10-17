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
argParser.add_argument('--isTest',              action='store_true')

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
tauDMs = ["inclusive", "dm0", "dm1", "dm10"]

isDMspesific = True

outputname = "/storage_mnt/storage/user/lwezenbe/CMSSW_10_2_13/src/TauTriggerTools/TauTagAndProbe/test/Factorization/data/Output/tauTriggerFactorization2018_"+str(args.subJob)+".root"


# get binning from binning2017.py and binning2018.py files.
from binning2018 import getbinning2018
bins = getbinning2018()
	
#bin2D = bins.getBinning()["ditau"]
#bin2DDM = bins.getBinning()["ditau"]
bin2D = np.linspace(40, 100, 10)
bin2DDM = np.linspace(40, 100, 10)
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
    highestPtIndex = 0
    if p2.Pt() > p1.Pt():
        highestPtIndex = 1
        
    pt_values = sorted([p1.Pt(), p2.Pt()], reverse=True)

    passed_ditau_trigger = tree.hasHLTPath_20

    Nevents = tree.EventNumber
    Nevts =Nevts + 1
    
    bkgSubW = 1.
    weight = bkgSubW*puweight

    WPointNames = {"vvlooseTauMVA":"vvlooseMVAv2", "vlooseTauMVA":"vlooseMVAv2", "looseTauMVA":"looseMVAv2", "mediumTauMVA":"mediumMVAv2", "tightTauMVA":"tightMVAv2", "vtightTauMVA":"vtightMVAv2", "vvtightTauMVA":"vvtightMVAv2"}
    WPoints = {"vvlooseTauMVA":vlooseWP, "vlooseTauMVA":vlooseWP, "looseTauMVA":looseWP, "mediumTauMVA":mediumWP, "tightTauMVA":tightWP, "vtightTauMVA":vtightWP, "vvtightTauMVA":vtightWP}
    DMmap = {0:'dm0', 1:'dm1', 10:'dm10'}

    # Filling the histograms
    for WP in WPs:

        if not WPoints[WP][0] or not WPoints[WP][1]: continue

        eff_curve_leg1 = eff_file.Get('ditau_'+WPointNames[WP]+"_"+DMmap[tauDM[highestPtIndex]]+'_MC_fit')
        eff_curve_leg2 = eff_file.Get('ditau_'+WPointNames[WP]+"_"+DMmap[tauDM[abs(highestPtIndex-1)]]+'_MC_fit')

        efficiency_leg1_DM = eff_curve_leg1.Eval(tauPt[highestPtIndex])
        efficiency_leg2_DM = eff_curve_leg2.Eval(tauPt[abs(highestPtIndex-1)])
        weightedDitauHists_1D.fillHist(WP, ('inclusive', 'inclusive'), mass, weight*efficiency_leg1_DM*efficiency_leg2_DM)
        weightedDitauHists_2D.fillHist(WP, ('inclusive', 'inclusive'), pt_values, weight*efficiency_leg1_DM*efficiency_leg2_DM)
        weightedDitauHists_1D.fillHist(WP, (DMmap[tauDM[highestPtIndex]], DMmap[tauDM[abs(highestPtIndex-1)]]), mass, weight*efficiency_leg1_DM*efficiency_leg2_DM)
        weightedDitauHists_2D.fillHist(WP,  (DMmap[tauDM[highestPtIndex]], DMmap[tauDM[abs(highestPtIndex-1)]]), pt_values, weight*efficiency_leg1_DM*efficiency_leg2_DM)
        denomHists_1D.fillHist(WP,('inclusive', 'inclusive'), mass, weight)
        denomHists_2D.fillHist(WP,('inclusive', 'inclusive'), pt_values, weight)
        denomHists_1D.fillHist(WP,(DMmap[tauDM[highestPtIndex]], DMmap[tauDM[abs(highestPtIndex-1)]]), mass, weight)
        denomHists_2D.fillHist(WP,(DMmap[tauDM[highestPtIndex]], DMmap[tauDM[abs(highestPtIndex-1)]]), pt_values, weight)
        

        if(passed_ditau_trigger):
            ditauTriggeredHists_1D.fillHist(WP, ('inclusive', 'inclusive'), mass, weight)
            ditauTriggeredHists_2D.fillHist(WP, ('inclusive', 'inclusive'), pt_values, weight)
            ditauTriggeredHists_1D.fillHist(WP, (DMmap[tauDM[highestPtIndex]], DMmap[tauDM[abs(highestPtIndex-1)]]), mass, weight)
            ditauTriggeredHists_2D.fillHist(WP, (DMmap[tauDM[highestPtIndex]], DMmap[tauDM[abs(highestPtIndex-1)]]), pt_values, weight)

if not args.isTest: 
    out_file = TFile( outputname, 'recreate')

    for WP in WPs:
        ditauTriggeredHists_1D.returnHist(WP, ('inclusive', 'inclusive')).Write(ditauTriggeredHists_1D.name + "_" +WP + "_inclusive")
        ditauTriggeredHists_2D.returnHist(WP, ('inclusive', 'inclusive')).Write(ditauTriggeredHists_2D.name + "_" +WP + "_inclusive")
        weightedDitauHists_1D.returnHist(WP, ('inclusive', 'inclusive')).Write(weightedDitauHists_1D.name + "_" +WP + "_inclusive")
        weightedDitauHists_2D.returnHist(WP, ('inclusive', 'inclusive')).Write(weightedDitauHists_2D.name + "_" +WP + "_inclusive")
        denomHists_1D.returnHist(WP, ('inclusive', 'inclusive')).Write(denomHists_1D.name + "_" +WP + "_inclusive")
        denomHists_2D.returnHist(WP, ('inclusive', 'inclusive')).Write(denomHists_2D.name + "_" +WP + "_inclusive")
        for tau1DM in tauDMs[1:]:
            for tau2DM in tauDMs[1:]:
                ditauTriggeredHists_1D.returnHist(WP, (tau1DM, tau2DM)).Write(ditauTriggeredHists_1D.name + "_" +WP + "_tau1_" + tau1DM + "_tau2_" + tau2DM)
                ditauTriggeredHists_2D.returnHist(WP, (tau1DM, tau2DM)).Write(ditauTriggeredHists_2D.name + "_" +WP + "_tau1_" + tau1DM + "_tau2_" + tau2DM)
                weightedDitauHists_1D.returnHist(WP, (tau1DM, tau2DM)).Write(weightedDitauHists_1D.name + "_" +WP + "_tau1_" + tau1DM + "_tau2_" + tau2DM)
                weightedDitauHists_2D.returnHist(WP, (tau1DM, tau2DM)).Write(weightedDitauHists_2D.name + "_" +WP + "_tau1_" + tau1DM + "_tau2_" + tau2DM)
                denomHists_1D.returnHist(WP, (tau1DM, tau2DM)).Write(denomHists_1D.name + "_" +WP + "_tau1_" + tau1DM + "_tau2_" + tau2DM)
                denomHists_2D.returnHist(WP, (tau1DM, tau2DM)).Write(denomHists_2D.name + "_" +WP + "_tau1_" + tau1DM + "_tau2_" + tau2DM)

    out_file.Close()
    print "The output ROOT file has been created: " + outputname
