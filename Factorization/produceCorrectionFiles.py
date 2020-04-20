#
#       ProduceCorrectionFiles.py
#
# Code to read files with measured and weighted efficiencies and 
# produce files that can be used for corrections
#

#Parse arguments
import argparse

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--isTest',              action='store_true')
argParser.add_argument('--ptThreshold',              action='store',       default = '50')
argParser.add_argument('--year',              action='store',       default = '2018')

args = argParser.parse_args()

WPs = ["vlooseTauMVA", "looseTauMVA", "mediumTauMVA", "tightTauMVA", "vtightTauMVA"]
tauDMs = ["dm0", "dm1", "dm10"]

weighted_hist = {}
triggered_hist = {}
ratio_hist = {}
correction_hist = {}

from ROOT import TFile
out_file = TFile('data/SF-corrections-2018.root', 'recreate')
from helpers import getObjFromFile
in_path = 'data/all/40/tauTriggerFactorization'+args.year+'-2D-v8.root'
for wp in WPs:
    weighted_hist[wp] = {}
    triggered_hist[wp] = {}
    ratio_hist[wp] = {}
    correction_hist[wp] = {}
    for tau1_dm in tauDMs:
        for tau2_dm in tauDMs:
            weighted_hist[(tau1_dm, tau2_dm)] = getObjFromFile(in_path, 'weightedDitau_2D_'+wp+'_tau1_'+tau1_dm+'_tau2_'+tau2_dm)
            triggered_hist[(tau1_dm, tau2_dm)] = getObjFromFile(in_path, 'ditauTriggered_2D_'+wp+'_tau1_'+tau1_dm+'_tau2_'+tau2_dm)
            out_file.cd()
            ratio_hist[(tau1_dm, tau2_dm)] = triggered_hist[(tau1_dm, tau2_dm)].Clone('ratio_'+wp+'_tau1_'+tau1_dm+'_tau2_'+tau2_dm)
            ratio_hist[(tau1_dm, tau2_dm)].Divide(weighted_hist[(tau1_dm, tau2_dm)])
            correction_hist[(tau1_dm, tau2_dm)] = ratio_hist[(tau1_dm, tau2_dm)].Clone('SFcorrection_'+wp+'_tau1_'+tau1_dm+'_tau2_'+tau2_dm)
            for bx in xrange(correction_hist[(tau1_dm, tau2_dm)].GetXaxis().GetNbins()+1):
                for by in xrange(correction_hist[(tau1_dm, tau2_dm)].GetYaxis().GetNbins()+1):
                    x = correction_hist[(tau1_dm, tau2_dm)].GetXaxis().GetBinUpEdge(bx)
                    y = correction_hist[(tau1_dm, tau2_dm)].GetXaxis().GetBinUpEdge(by)
                    if x > int(args.ptThreshold) and y > int(args.ptThreshold):
                        correction_hist[(tau1_dm, tau2_dm)].SetBinContent(bx, by, 1.)
                        correction_hist[(tau1_dm, tau2_dm)].SetBinError(bx, by, 0.)
            correction_hist[(tau1_dm, tau2_dm)].Write()

out_file.Close()

