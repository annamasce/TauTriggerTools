import ROOT
import numpy as np
class histMap():
    
    def __init__(self, name, WPs, tauDMs, bins, binsDM, isTH2 = False):
        self.name = name
        self.WPs = WPs
        self.tauDMs = tauDMs
        self.hist = {}
        self.isTH2 = isTH2
        
        for WP in WPs:
            for tau1DM in tauDMs:
                for tau2DM in tauDMs:
                    if isTH2:
                        self.hist[(WP, (tau1DM, tau2DM))] = ROOT.TH2D(name + "_" + WP + "_" + tau1DM + '_' + tau2DM, name + "_" + WP + "_" + tau1DM + '_' + tau2DM, len(bins)-1, bins, len(bins)-1, bins)
                        self.hist[(WP, (tau1DM, tau2DM))].Sumw2()
                    else:
                        self.hist[(WP, (tau1DM, tau2DM))] = ROOT.TH1D(name + "_" + WP + "_" + tau1DM + '_' + tau2DM, name + "_" + WP + "_" + tau1DM + '_' + tau2DM, len(bins)-1, bins)
                        self.hist[(WP, (tau1DM, tau2DM))].Sumw2()
    

    def fillHist(self, WP, tauDM, values, weight = 1.):
        if self.isTH2:
            self.hist[(WP, tauDM)].Fill(values[0], values[1], weight)
        else:
            self.hist[(WP, tauDM)].Fill(values, weight)

    def returnHist(self, WP, tauDM):
        return self.hist[(WP, tauDM)]


