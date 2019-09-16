from plottingTools import *
import ROOT

#in_file = ROOT.TFile.Open("data/tauTriggerFactorization2018-full.root")
in_file = ROOT.TFile.Open("data/tauTriggerFactorization2018-2D-v3.root")
WPs = ['vlooseTauMVA', 'looseTauMVA', 'mediumTauMVA', 'tightTauMVA', 'vtightTauMVA']

#Parse arguments
import argparse

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--isTH2',              action='store',         default=None)

args = argParser.parse_args()


for WP in WPs:
    triggeredHist_2D = in_file.Get('ditauTriggered_2D_'+WP+'_inclusive')
    weightedHist_2D = in_file.Get("weightedDitau_2D_"+WP+"_inclusive")
    denomHist_2D = in_file.Get("denomHists_2D_"+WP+"_inclusive")
    triggeredHist_1D = in_file.Get('ditauTriggered_mass_'+WP+'_inclusive')
    weightedHist_1D = in_file.Get("weightedDitau_mass_"+WP+"_inclusive")
    denomHist_1D = in_file.Get("denomHists_mass_"+WP+"_inclusive")

    triggeredHist_2D.Divide(denomHist_2D)
    weightedHist_2D.Divide(denomHist_2D)
    triggeredHist_1D.Divide(denomHist_1D)
    weightedHist_1D.Divide(denomHist_1D)

    destination = "/user/lwezenbe/CMSSW_10_2_13/src/TauTriggerTools/TauTagAndProbe/test/Factorization/Results/Factorization/"+WP
    
    draw2DHist(triggeredHist_2D, "p_{T}^{#tau_{1}} [GeV]", "p_{T}^{#tau_{2}} [GeV]", destination+"_2Dtriggered", option="Colz", x_log=True, y_log=True)
    draw2DHist(weightedHist_2D, "p_{T}^{#tau_{1}} [GeV]", "p_{T}^{#tau_{2}} [GeV]", destination+"_2Dweighted", option="Colz", x_log=True, y_log=True)
    ratio_hist_2D = weightedHist_2D.Clone("ratio plot 2D")
    ratio_hist_2D.Divide(triggeredHist_2D)
    draw2DHist(ratio_hist_2D, "p_{T}^{#tau_{1}} [GeV]", "p_{T}^{#tau_{2}} [GeV]", destination+"_2Dratio", option="Colz", x_log=True, y_log=True)
        
    plotClosure(triggeredHist_1D, weightedHist_1D, "M(#tau_{h}#tau_{h})", "Events", "DYJets", destination, yLog = False)
