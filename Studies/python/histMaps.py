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

    def overflowvalue(self, WP, tauDM, value, region, isY = False):
        
        if isY:
            last_val = self.hist[(WP, tauDM)].GetYaxis().GetBinCenter(self.hist[(WP, tauDM)].GetYaxis().GetLast())
        else:
            last_val = self.hist[(WP, tauDM)].GetXaxis().GetBinCenter(self.hist[(WP, tauDM)].GetXaxis().GetLast())
       
        val = None 
        if region == 'high' or region == 'all':
            val = min(value, last_val)
        if region == 'low':
            val = value

        return val

    def fillHist(self, WP, tauDM, values, region, weight = 1.):
        if self.isTH2:
            xval = self.overflowvalue(WP, tauDM, values[0], region)
            yval = self.overflowvalue(WP, tauDM, values[1], region, isY=True)
 
            self.hist[(WP, tauDM)].Fill(xval, yval, weight)
        else:
            xval = self.overflowvalue(WP, tauDM, values, region)
            self.hist[(WP, tauDM)].Fill(xval, weight)
        return

    def returnHist(self, WP, tauDM):
        return self.hist[(WP, tauDM)]


