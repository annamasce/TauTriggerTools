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
            for tauDM in tauDMs:
                if isTH2:
                    self.hist[(WP, tauDM)] = ROOT.TH2D(name + "_" + WP + "_" + tauDM, name + "_" + WP + "_" + tauDM, len(bins)-1, bins, len(bins)-1, bins)
                    self.hist[(WP, tauDM)].Sumw2()
                else:
                    self.hist[(WP, tauDM)] = ROOT.TH1D(name + "_" + WP + "_" + tauDM, name + "_" + WP + "_" + tauDM, len(bins)-1, bins)
                    self.hist[(WP, tauDM)].Sumw2()
    

    def fillHist(self, WP, tauDM, values, weight = 1.):
        if self.isTH2:
            self.hist[(WP, tauDM)].Fill(values[0], values[1], weight)
        else:
            self.hist[(WP, tauDM)].Fill(values, weight)

    def returnHist(self, WP, tauDM):
        return self.hist[(WP, tauDM)]


