from ROOT import *
from array import array
gROOT.SetBatch(True)
from math import sqrt
from histMaps import histMap
import numpy as np
from helpers_old import progress, makeDirIfNeeded
from array import array

#Parse arguments
import argparse

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--subJob',              action='store',         default=0)
argParser.add_argument('--isTest',              action='store_true')
argParser.add_argument('--ptRegion',              action='store',       default = 'high',       choices=['low', 'high', 'all'])
argParser.add_argument('--ptThreshold',              action='store',       default = '40')
argParser.add_argument('--ptThresholdLeading',              action='store',       default = None)

args = argParser.parse_args()

#Load in samples
import Sample
sampleList = Sample.createSampleList('/storage_mnt/storage/user/lwezenbe/private/CMSSW_10_2_17/src/TauTriggerTools/Factorization/data/inputFiles_2018.conf')
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

if args.ptThresholdLeading:
    thresholdName = args.ptThresholdLeading +'-'+args.ptThreshold
else:
    thresholdName = args.ptThreshold
#outputname = "/storage_mnt/storage/user/lwezenbe/private/CMSSW_10_2_17/src/TauTriggerTools/Factorization/data/Output/test/"+args.ptRegion+"/"+thresholdName+"/tauTriggerFactorization2018_"+str(args.subJob)+".root"
outputname = "/storage_mnt/storage/user/lwezenbe/private/CMSSW_10_2_17/src/TauTriggerTools/Factorization/data/Output/"+args.ptRegion+"/"+thresholdName+"/tauTriggerFactorization2018_"+str(args.subJob)+".root"
makeDirIfNeeded(outputname.rsplit('/', 1)[0])

#
#       Set variable range ranges
#
var = {'high':  {'2D':[40., 45., 50., 60., 75., 90., 115., 175., 200.],
                'leadingpt':[40., 45., 50., 60., 75., 95., 120., 200.],
                'subleadingpt':[40., 45., 50., 60., 75., 95., 120., 200.],
                #'leadingpt':[40.0, 45.0, 50.0, 55.0, 65.0, 85.0, 110.0, 410.0],
                #'subleadingpt':[40.0, 45.0, 50.0, 55.0, 65.0, 85.0, 110.0, 410.0],
                },
    'low':  {'2D':[20., 24., 28., 32., 36., 40.],
            'leadingpt':[20., 24., 28., 32., 36., 40.],
            'subleadingpt':[20., 24., 28., 32., 36., 40.]
            },
    'all': {'2D':[20., 25., 30., 35., 40., 45., 50., 60., 70., 85., 100., 150., 200.],
            'leadingpt':[20., 25., 30., 35., 40., 50., 60., 70., 85., 100., 150., 200.],
            'subleadingpt':[20., 25., 30., 35., 40., 50., 60., 70., 85., 100., 150., 200.]
            }
        }

#make adjustments to these ranges if needed
#Append values below threshold on low and remove from high
for v in var['high'].keys():
    tmp_list = [x for x in var['high'][v]]  #Make temporary copy
    var['high'][v] = [x for x in var['high'][v] if x >= float(args.ptThreshold)]
    var['low'][v].extend([item for item in tmp_list if item not in var['high'][v] and item not in var['low'][v]])
    var['low'][v].append(var['high'][v][0])

ditauTriggeredHists_1D_leading = histMap("ditauTriggered_leadingpt", WPs, tauDMs, np.array(var[args.ptRegion]['leadingpt']), np.array(var[args.ptRegion]['leadingpt']))
weightedDitauHists_1D_leading = histMap("weightedDitau_leadingpt", WPs, tauDMs, np.array(var[args.ptRegion]['leadingpt']), np.array(var[args.ptRegion]['leadingpt']))
denomHists_1D_leading = histMap('denomHists_leadingpt', WPs, tauDMs, np.array(var[args.ptRegion]['leadingpt']), np.array(var[args.ptRegion]['leadingpt']))

ditauTriggeredHists_1D_subleading = histMap("ditauTriggered_subleadingpt", WPs, tauDMs, np.array(var[args.ptRegion]['subleadingpt']), np.array(var[args.ptRegion]['subleadingpt']))
weightedDitauHists_1D_subleading = histMap("weightedDitau_subleadingpt", WPs, tauDMs, np.array(var[args.ptRegion]['subleadingpt']), np.array(var[args.ptRegion]['subleadingpt']))
denomHists_1D_subleading = histMap('denomHists_subleadingpt', WPs, tauDMs, np.array(var[args.ptRegion]['subleadingpt']), np.array(var[args.ptRegion]['subleadingpt']))

ditauTriggeredHists_2D = histMap("ditauTriggered_2D", WPs, tauDMs, array('f',np.array(var[args.ptRegion]['2D'])), np.array(var[args.ptRegion]['2D']), isTH2=True)
weightedDitauHists_2D = histMap("weightedDitau_2D", WPs, tauDMs, array('f', np.array(var[args.ptRegion]['2D'])), np.array(var[args.ptRegion]['2D']), isTH2=True)
denomHists_2D = histMap('denomHists_2D', WPs, tauDMs, np.array(var[args.ptRegion]['2D']), np.array(var[args.ptRegion]['2D']), isTH2=True)

eff_file = TFile.Open('/storage_mnt/storage/user/lwezenbe/private/CMSSW_10_2_17/src/TauTriggerTools/Factorization/tauTriggerEfficiencies2018.root')

print "Populating histograms"

DMmap = {0:'dm0', 1:'dm1', 10:'dm10'}
Nevts = 0

def getYval(graph, xval):

    for i, x in enumerate(np.array(graph.GetX())):
        if xval < x:
            if xval < x-np.array(graph.GetErrorX(i)):     return np.array(graph.GetY())[i-1]
            else:                                       return np.array(graph.GetY())[i]
    else: return 0.


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

    leadingptcut = args.ptThresholdLeading if args.ptThresholdLeading else args.ptThreshold
    if args.ptRegion == 'high' and (tree.tauPt[0] < float(leadingptcut) or tree.tauPt[1] < float(args.ptThreshold)): continue
 
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

    # Filling the histograms
    for WP in WPs:

        if not WPoints[WP][0] or not WPoints[WP][1]: continue

        eff_curve_leg1 = eff_file.Get('ditau_'+WPointNames[WP]+"_"+DMmap[tauDM[highestPtIndex]]+'_MC_fit')
        eff_curve_leg2 = eff_file.Get('ditau_'+WPointNames[WP]+"_"+DMmap[tauDM[abs(highestPtIndex-1)]]+'_MC_fit')
        #eff_curve_leg1 = eff_file.Get('ditau_'+WPointNames[WP]+"_"+DMmap[tauDM[highestPtIndex]]+'_MC_graph')
        #eff_curve_leg2 = eff_file.Get('ditau_'+WPointNames[WP]+"_"+DMmap[tauDM[abs(highestPtIndex-1)]]+'_MC_graph')

        efficiency_leg1_DM = eff_curve_leg1.Eval(tauPt[highestPtIndex])
        efficiency_leg2_DM = eff_curve_leg2.Eval(tauPt[abs(highestPtIndex-1)])
        #efficiency_leg1_DM = getYval(eff_curve_leg1, tauPt[highestPtIndex])
        #efficiency_leg2_DM = getYval(eff_curve_leg2, tauPt[abs(highestPtIndex-1)])
        weightedDitauHists_1D_leading.fillHist(WP, ('inclusive', 'inclusive'), pt_values[0], args.ptRegion, weight*efficiency_leg1_DM*efficiency_leg2_DM)
        weightedDitauHists_1D_subleading.fillHist(WP, ('inclusive', 'inclusive'), pt_values[1],args.ptRegion, weight*efficiency_leg1_DM*efficiency_leg2_DM)
        weightedDitauHists_2D.fillHist(WP, ('inclusive', 'inclusive'), pt_values, args.ptRegion,weight*efficiency_leg1_DM*efficiency_leg2_DM)
        weightedDitauHists_1D_leading.fillHist(WP, (DMmap[tauDM[highestPtIndex]], DMmap[tauDM[abs(highestPtIndex-1)]]), pt_values[0], args.ptRegion,weight*efficiency_leg1_DM*efficiency_leg2_DM)
        weightedDitauHists_1D_subleading.fillHist(WP, (DMmap[tauDM[highestPtIndex]], DMmap[tauDM[abs(highestPtIndex-1)]]), pt_values[1], args.ptRegion,weight*efficiency_leg1_DM*efficiency_leg2_DM)
        weightedDitauHists_2D.fillHist(WP,  (DMmap[tauDM[highestPtIndex]], DMmap[tauDM[abs(highestPtIndex-1)]]), pt_values, args.ptRegion,weight*efficiency_leg1_DM*efficiency_leg2_DM)
        denomHists_1D_leading.fillHist(WP,('inclusive', 'inclusive'), pt_values[0], args.ptRegion,weight)
        denomHists_1D_subleading.fillHist(WP,('inclusive', 'inclusive'), pt_values[1], args.ptRegion,weight)
        denomHists_2D.fillHist(WP,('inclusive', 'inclusive'), pt_values,args.ptRegion, weight)
        denomHists_1D_leading.fillHist(WP,(DMmap[tauDM[highestPtIndex]], DMmap[tauDM[abs(highestPtIndex-1)]]), pt_values[0], args.ptRegion,weight)
        denomHists_1D_subleading.fillHist(WP,(DMmap[tauDM[highestPtIndex]], DMmap[tauDM[abs(highestPtIndex-1)]]), pt_values[1], args.ptRegion,weight)
        denomHists_2D.fillHist(WP,(DMmap[tauDM[highestPtIndex]], DMmap[tauDM[abs(highestPtIndex-1)]]), pt_values, args.ptRegion,weight)
        

        if(passed_ditau_trigger):
            ditauTriggeredHists_1D_leading.fillHist(WP, ('inclusive', 'inclusive'), pt_values[0], args.ptRegion,weight)
            ditauTriggeredHists_1D_subleading.fillHist(WP, ('inclusive', 'inclusive'), pt_values[1], args.ptRegion,weight)
            ditauTriggeredHists_2D.fillHist(WP, ('inclusive', 'inclusive'), pt_values,args.ptRegion, weight)
            ditauTriggeredHists_1D_leading.fillHist(WP, (DMmap[tauDM[highestPtIndex]], DMmap[tauDM[abs(highestPtIndex-1)]]), pt_values[0],args.ptRegion, weight)
            ditauTriggeredHists_1D_subleading.fillHist(WP, (DMmap[tauDM[highestPtIndex]], DMmap[tauDM[abs(highestPtIndex-1)]]), pt_values[1],args.ptRegion, weight)
            ditauTriggeredHists_2D.fillHist(WP, (DMmap[tauDM[highestPtIndex]], DMmap[tauDM[abs(highestPtIndex-1)]]), pt_values, args.ptRegion,weight)

tot = [0., 0.]
for x in DMmap.keys():
    for y in DMmap.keys():
        tot[0] += ditauTriggeredHists_2D.returnHist('mediumTauMVA', (DMmap[x], DMmap[y])).GetBinContent(1, 1)
        tot[1] += weightedDitauHists_2D.returnHist('mediumTauMVA', (DMmap[x], DMmap[y])).GetBinContent(1, 1)

print tot, ditauTriggeredHists_2D.returnHist('mediumTauMVA', ('inclusive', 'inclusive')).GetBinContent(1, 1), weightedDitauHists_2D.returnHist('mediumTauMVA', ('inclusive', 'inclusive')).GetBinContent(1, 1)

if not args.isTest: 
    out_file = TFile( outputname, 'recreate')

    for WP in WPs:
        ditauTriggeredHists_1D_leading.returnHist(WP, ('inclusive', 'inclusive')).Write(ditauTriggeredHists_1D_leading.name + "_" +WP + "_inclusive")
        ditauTriggeredHists_1D_subleading.returnHist(WP, ('inclusive', 'inclusive')).Write(ditauTriggeredHists_1D_subleading.name + "_" +WP + "_inclusive")
        ditauTriggeredHists_2D.returnHist(WP, ('inclusive', 'inclusive')).Write(ditauTriggeredHists_2D.name + "_" +WP + "_inclusive")
        weightedDitauHists_1D_leading.returnHist(WP, ('inclusive', 'inclusive')).Write(weightedDitauHists_1D_leading.name + "_" +WP + "_inclusive")
        weightedDitauHists_1D_subleading.returnHist(WP, ('inclusive', 'inclusive')).Write(weightedDitauHists_1D_subleading.name + "_" +WP + "_inclusive")
        weightedDitauHists_2D.returnHist(WP, ('inclusive', 'inclusive')).Write(weightedDitauHists_2D.name + "_" +WP + "_inclusive")
        denomHists_1D_leading.returnHist(WP, ('inclusive', 'inclusive')).Write(denomHists_1D_leading.name + "_" +WP + "_inclusive")
        denomHists_1D_subleading.returnHist(WP, ('inclusive', 'inclusive')).Write(denomHists_1D_subleading.name + "_" +WP + "_inclusive")
        denomHists_2D.returnHist(WP, ('inclusive', 'inclusive')).Write(denomHists_2D.name + "_" +WP + "_inclusive")
        for tau1DM in tauDMs[1:]:
            for tau2DM in tauDMs[1:]:
                ditauTriggeredHists_1D_leading.returnHist(WP, (tau1DM, tau2DM)).Write(ditauTriggeredHists_1D_leading.name + "_" +WP + "_tau1_" + tau1DM + "_tau2_" + tau2DM)
                ditauTriggeredHists_1D_subleading.returnHist(WP, (tau1DM, tau2DM)).Write(ditauTriggeredHists_1D_subleading.name + "_" +WP + "_tau1_" + tau1DM + "_tau2_" + tau2DM)
                ditauTriggeredHists_2D.returnHist(WP, (tau1DM, tau2DM)).Write(ditauTriggeredHists_2D.name + "_" +WP + "_tau1_" + tau1DM + "_tau2_" + tau2DM)
                weightedDitauHists_1D_leading.returnHist(WP, (tau1DM, tau2DM)).Write(weightedDitauHists_1D_leading.name + "_" +WP + "_tau1_" + tau1DM + "_tau2_" + tau2DM)
                weightedDitauHists_1D_subleading.returnHist(WP, (tau1DM, tau2DM)).Write(weightedDitauHists_1D_subleading.name + "_" +WP + "_tau1_" + tau1DM + "_tau2_" + tau2DM)
                weightedDitauHists_2D.returnHist(WP, (tau1DM, tau2DM)).Write(weightedDitauHists_2D.name + "_" +WP + "_tau1_" + tau1DM + "_tau2_" + tau2DM)
                denomHists_1D_leading.returnHist(WP, (tau1DM, tau2DM)).Write(denomHists_1D_leading.name + "_" +WP + "_tau1_" + tau1DM + "_tau2_" + tau2DM)
                denomHists_1D_subleading.returnHist(WP, (tau1DM, tau2DM)).Write(denomHists_1D_subleading.name + "_" +WP + "_tau1_" + tau1DM + "_tau2_" + tau2DM)
                denomHists_2D.returnHist(WP, (tau1DM, tau2DM)).Write(denomHists_2D.name + "_" +WP + "_tau1_" + tau1DM + "_tau2_" + tau2DM)

    out_file.Close()
    print "The output ROOT file has been created: " + outputname
