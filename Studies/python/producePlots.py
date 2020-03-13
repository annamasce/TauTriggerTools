from TauTriggerTools.Studies.plottingTools import *
import ROOT
from TauTriggerTools.Studies.helpers import *

import argparse

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--ptRegion',              action='store',       default = 'high',       choices=['low', 'high', 'all'])
argParser.add_argument('--WP',                    action='store',       default = 'tight',       choices=['vvvloose', 'vvloose', 'vloose', 'loose', 'medium', 'tight', 'vtight', 'vvtight'])
argParser.add_argument('--var',                   action='store',       default = 'leadingpt') 
argParser.add_argument('--ptThreshold',           action='store',       default = '40') 
argParser.add_argument('--selection',             action='store',       required=True, choices=['DeepTau', 'MVA']) 

args = argParser.parse_args()

in_file = ROOT.TFile.Open(os.path.expandvars("$CMSSW_BASE/src/TauTriggerTools/Studies/data/Factorization/"+args.selection+"/"+args.ptRegion+"/"+args.ptThreshold+"/tauTriggerFactorization2018_0.root"))
WP = args.WP
tauDMs = ['dm0', 'dm1', 'dm10']

histNames = ['inclusive']
for tau1DM in tauDMs:
    for tau2DM in tauDMs:
        histNames.append('tau1_'+tau1DM+'_tau2_'+tau2DM)

categNames = ['2D', 'leadingpt', 'subleadingpt']

for h in histNames:
    #[text, xpos, ypos, textsize, align] 


    triggeredExtraInfo = []
    WeightedExtraInfo = []
    RatioExtraInfo = []
    triggeredRawExtraInfo = []
    WeightedRawExtraInfo = []
    denomRawExtraInfo = []
    significanceExtraInfo = []

    triggeredExtraInfo.append(extraTextFormat('HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v', 0.2, 0.9, 0.02))
    WeightedExtraInfo.append(extraTextFormat('HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v', 0.2, 0.9, 0.02))
    RatioExtraInfo.append(extraTextFormat('HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v', 0.2, 0.9, 0.02))
    triggeredRawExtraInfo.append(extraTextFormat('HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v', 0.2, 0.9, 0.02))
    WeightedRawExtraInfo.append(extraTextFormat('HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v', 0.2, 0.9, 0.02))
    denomRawExtraInfo.append(extraTextFormat('HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v', 0.2, 0.9, 0.02))
    significanceExtraInfo.append(extraTextFormat('HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v', 0.2, 0.9, 0.02))
    triggeredExtraInfo.append(extraTextFormat('Triggered'))
    WeightedExtraInfo.append(extraTextFormat('Weighted'))
    RatioExtraInfo.append(extraTextFormat('Weighted/Triggered'))
    triggeredRawExtraInfo.append(extraTextFormat('Triggered'))
    significanceExtraInfo.append(extraTextFormat('Significance (|ratio-1|/error)'))
    triggeredExtraInfo.append(extraTextFormat(WP))
    WeightedExtraInfo.append(extraTextFormat(WP))
    RatioExtraInfo.append(extraTextFormat(WP))
    triggeredRawExtraInfo.append(extraTextFormat(WP))
    WeightedRawExtraInfo.append(extraTextFormat(WP))
    denomRawExtraInfo.append(extraTextFormat(WP))
    significanceExtraInfo.append(extraTextFormat(WP))

    triggeredHist = in_file.Get('ditauTriggered_'+args.var+'_'+WP+'_'+h)
    weightedHist = in_file.Get("weightedDitau_"+args.var+"_"+WP+"_"+h)
    denomHist = in_file.Get("denomHists_"+args.var+"_"+WP+"_"+h)
    print 'ditauTriggered_'+args.var+'_'+WP+'_'+h
    triggeredRawHist = triggeredHist.Clone('triggered variable distribution')
    weightedRawHist = weightedHist.Clone('weighted variable distribution')

    triggeredHist.Divide(triggeredHist, denomHist, 1., 1., 'B')
    weightedHist.Divide(weightedHist, denomHist, 1., 1., 'B')

    axis_max = None
    if args.ptRegion == 'high': axis_max = 115.

    destination = os.path.expandvars("$CMSSW_BASE/src/TauTriggerTools/Studies/data/Factorization/Results/Factorization/HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v/"+args.selection+"/"+args.ptRegion+"/"+args.ptThreshold+"/"+args.var+'/'+WP+'_'+h)
    makeDirIfNeeded(destination)
    if args.var == '2D': 
        draw2DHist(triggeredHist, "p_{T}^{#tau_{1}}(leading) [GeV]", "p_{T}^{#tau_{2}}(subleading) [GeV]", destination+"_2DtriggeredEff", option="EtextColz", x_log=False, y_log=False, extraText=(triggeredExtraInfo), axis_max=axis_max)
        draw2DHist(weightedHist, "p_{T}^{#tau_{1}}(leading) [GeV]", "p_{T}^{#tau_{2}}(subleading) [GeV]", destination+"_2DweightedEff", option="EtextColz", x_log=False, y_log=False, extraText=(WeightedExtraInfo), axis_max=axis_max)
        ratio_hist = weightedHist.Clone("ratio plot")
        ratio_hist.Divide(triggeredHist)
        draw2DHist(ratio_hist, "p_{T}^{#tau_{1}}(leading) [GeV]", "p_{T}^{#tau_{2}}(subleading) [GeV]", destination+"_2DratioEff", option="EtextColz", x_log=False, y_log=False, extraText=(RatioExtraInfo))
        significance_hist = ratio_hist.Clone('significance')
        for xbin in xrange(significance_hist.GetXaxis().GetNbins()):
            for ybin in xrange(significance_hist.GetYaxis().GetNbins()):
                if significance_hist.GetBinError(xbin, ybin) != 0: significance_hist.SetBinContent(xbin, ybin, abs((significance_hist.GetBinContent(xbin, ybin) - 1.)/significance_hist.GetBinError(xbin, ybin))) 
        draw2DHist(significance_hist, "p_{T}^{#tau_{1}}(leading) [GeV]", "p_{T}^{#tau_{2}}(subleading) [GeV]", destination+"_2Dsig", option="TextColz", x_log=False, y_log=False, extraText=(significanceExtraInfo))
        
            
    else: 
        #plotClosure(triggeredHist, weightedHist, "p_{T}^{#tau} [GeV]", "Efficiency", "DYJets", destination +'_'+args.var+'_wOverlap', yLog = False, denominator_shape=denomHist)
        plotClosure(triggeredRawHist, weightedRawHist, "p_{T}^{#tau} [GeV]", "Taus", "DYJets", destination +'_'+args.var+'_Events', yLog = False)
        plotClosure(triggeredHist, weightedHist, "p_{T}^{#tau} [GeV]", "Efficiency", "DYJets", destination +'_'+args.var, yLog = False)
        draw1DHist([triggeredRawHist, weightedRawHist], "p_{T}^{#tau} [GeV]", "Events", destination+"_triggeredEvents", legend=['triggered', 'weighted'], x_log=False, y_log=False, extraText=(triggeredRawExtraInfo))
        draw1DHist(denomHist, "p_{T}^{#tau} [GeV]", "Events", destination+"_BaselineEvents", x_log=False, y_log=False, extraText=(denomRawExtraInfo))

    del triggeredHist
    del weightedHist
    del denomHist
    del triggeredRawHist
    del weightedRawHist

