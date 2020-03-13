from ROOT import *
from array import array
gROOT.SetBatch(True)
from math import sqrt
from TauTriggerTools.Studies.histMaps import histMap
import numpy as np
from TauTriggerTools.Studies.helpers import progress, makeDirIfNeeded
from TauTriggerTools.Common.AnalysisTypes import *
from TauTriggerTools.Common.AnalysisTools import *

from array import array
import os

#Parse arguments
import argparse

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--subJob',              action='store',         default=0)
argParser.add_argument('--isTest',              action='store_true')
argParser.add_argument('--ptRegion',              action='store',       default = 'high',       choices=['low', 'high', 'all'])
argParser.add_argument('--ptThreshold',              action='store',       default = '40')
argParser.add_argument('--ptThresholdLeading',              action='store',       default = None)
argParser.add_argument('--selection', required=True, type=str, help="tau selection", choices=['DeepTau', 'MVA'])


args = argParser.parse_args()

#Load in samples
import TauTriggerTools.Studies.Sample as Sample
sampleList = Sample.createSampleList(os.path.expandvars('$CMSSW_BASE/src/TauTriggerTools/Studies/data/inputFiles_Factorization_2018.conf'))
sample = Sample.getSampleFromList(sampleList, 'noTagAndProbe_'+args.selection)

in_file = TFile.Open(sample.path)
tree = sample.initTree()
if args.isTest:
    eventRange=xrange(200)
else:
    eventRange = sample.getEventRange(int(args.subJob))

gStyle.SetFrameLineWidth(1)
gStyle.SetPadBottomMargin(0.13)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)

WPs = ['vloose', 'loose', 'medium', 'tight', 'vtight']
tauDMs = ["inclusive", "dm0", "dm1", "dm10"]

isDMspesific = True

if args.ptThresholdLeading:
    thresholdName = args.ptThresholdLeading +'-'+args.ptThreshold
else:
    thresholdName = args.ptThreshold
outputname = os.path.expandvars("$CMSSW_BASE/src/TauTriggerTools/Studies/data/Factorization/"+args.selection+"/"+args.ptRegion+"/"+thresholdName+"/tauTriggerFactorization2018_"+str(args.subJob)+".root")
makeDirIfNeeded(outputname)
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

if args.selection == 'DeepTau':
    eff_file = TFile.Open(os.path.expandvars('$CMSSW_BASE/src/TauTriggerTools/Studies/data/SF/2018_tauTriggerEff_DeepTau2017v2p1.root'))
else:
    eff_file = TFile.Open(os.path.expandvars('$CMSSW_BASE/src/TauTriggerTools/Studies/data/SF/tauTriggerEfficiencies2018.root'))

print "Populating histograms"

DMmap = {0:'dm0', 1:'dm1', 10:'dm10'}
Nevts = 0

WPoints = {'vvvloose' : DiscriminatorWP.VVVLoose,
        'vvloose'  : DiscriminatorWP.VVLoose,
        'vloose'   : DiscriminatorWP.VLoose,
        'loose'    : DiscriminatorWP.Loose,
        'medium'   : DiscriminatorWP.Medium,
        'tight'    : DiscriminatorWP.Tight,
        'vtight'    : DiscriminatorWP.VTight,
        'vvtight'    : DiscriminatorWP.VVTight,
        'vvvtight'    : DiscriminatorWP.VVVTight,
        }
 
def getYval(graph, xval):

    for i, x in enumerate(np.array(graph.GetX())):
        if xval < x:
            if xval < x-np.array(graph.GetErrorX(i)):     return np.array(graph.GetY())[i-1]
            else:                                       return np.array(graph.GetY())[i]
    else: return 0.

test = 0
test1 = 0
test2 = 0
for iEv in eventRange:
    tree.GetEntry(iEv)
    progress(iEv - eventRange[0], len(eventRange))       
    test += 1
    tauPt = tree.tau_pt
    tauEta = tree.tau_eta
    tauPhi = tree.tau_phi
    tauMass = tree.tau_mass
    tauDM = tree.tau_decayMode

    if args.selection == 'DeepTau': tauSel = tree.byDeepTau2017v2p1VSjet
    else: tauSel = tree.byIsolationMVArun2017v2DBoldDMwLT2017
 
    if 11 in tauDM: continue

    # tau Energy Shift (SF) is only applied for 2017 for now!
    if("DYJets" in in_file):
        puweight = tree.puweight
        # tauPt_ESshifted = tree.tauPt  
    else:
        puweight = 1.
        # tauPt_ESshifted = tree.tauPt 

    p1 = TLorentzVector()
    p1.SetPtEtaPhiM(tauPt[0], tauEta[0], tauPhi[0], tauMass[0])
    p2 = TLorentzVector()
    p2.SetPtEtaPhiM(tauPt[1], tauEta[1], tauPhi[1], tauMass[1])
    mass = (p1+p2).M()
    highestPtIndex = 0
    if p2.Pt() > p1.Pt():
        highestPtIndex = 1
        
    pt_values = sorted([p1.Pt(), p2.Pt()], reverse=True)

    leadingptcut = args.ptThresholdLeading if args.ptThresholdLeading else args.ptThreshold
    print pt_values
    if args.ptRegion == 'high' and (pt_values[0] < float(leadingptcut) or pt_values[1] < float(args.ptThreshold)): continue

    passed_ditau_trigger = tree.pass_ditau

    Nevts =Nevts + 1
    
    bkgSubW = 1.
    weight = bkgSubW*puweight

    WPointNames = {'MVA':{"vvloose":"vvlooseMVAv2", "vloose":"vlooseMVAv2", "loose":"looseMVAv2", "medium":"mediumMVAv2", "tight":"tightMVAv2", "vtight":"vtightMVAv2", "vvtight":"vvtightMVAv2"},
                    'DeepTau':{"vvloose":"VVLoose", "vloose":"VLoose", "loose":"Loose", "medium":"Medium", "tight":"Tight", "vtight":"VTight", "vvtight":"VVTight"}}
    # WPoints = {"vvlooseTauMVA":vlooseWP, "vlooseTauMVA":vlooseWP, "looseTauMVA":looseWP, "mediumTauMVA":mediumWP, "tightTauMVA":tightWP, "vtightTauMVA":vtightWP, "vvtightTauMVA":vtightWP}

    test1 += 1
    # Filling the histograms
    for WP in WPs:
        if (int(tauSel[0]) & (1 << WPoints[WP])) == 0 or (int(tauSel[1]) & (1 << WPoints[WP])) == 0: continue
	if WP == 'vloose': test2 += 1
        if args.selection == 'DeepTau':
            eff_curve_leg1 = eff_file.Get('mc_ditau_'+WPointNames[args.selection][WP]+"_"+DMmap[tauDM[highestPtIndex]]+'_fitted')
            eff_curve_leg2 = eff_file.Get('mc_ditau_'+WPointNames[args.selection][WP]+"_"+DMmap[tauDM[abs(highestPtIndex-1)]]+'_fitted')
            efficiency_leg1_DM = eff_curve_leg1.GetBinContent(eff_curve_leg1.FindBin(tauPt[highestPtIndex]))
            efficiency_leg2_DM = eff_curve_leg2.GetBinContent(eff_curve_leg2.FindBin(tauPt[abs(highestPtIndex-1)]))
        else:
            eff_curve_leg1 = eff_file.Get('ditau_'+WPointNames[args.selection][WP]+"_"+DMmap[tauDM[highestPtIndex]]+'_MC_graph')
            eff_curve_leg2 = eff_file.Get('ditau_'+WPointNames[args.selection][WP]+"_"+DMmap[tauDM[abs(highestPtIndex-1)]]+'_MC_graph')
            efficiency_leg1_DM = eff_curve_leg1.Eval(tauPt[highestPtIndex])
            efficiency_leg2_DM = eff_curve_leg2.Eval(tauPt[abs(highestPtIndex-1)])

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

print test, test1, test2

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
