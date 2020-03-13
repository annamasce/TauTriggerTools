import ROOT
from ROOT import TLine, THStack, TColor, TCanvas, TLegend, gROOT, gStyle, TH1, TMultiGraph, TGraphErrors, TGraph, TLatex, TPad, TMath
#import CMS_lumi as cl
import tdrstyle as tdr
import numpy as np
from TauTriggerTools.Studies.helpers import makeDirIfNeeded, sortByOtherList
import os

def GeneralSettings(paintformat = "4.2f"):
    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)
    tdr.setTDRStyle()
    gStyle.SetPaintTextFormat(paintformat)
    gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")

def GetStackColor(index):
    if index == 0:      return "#4C5760"
    if index == 1:      return "#93A8AC"
    if index == 2:      return "#D7CEB2"
    if index == 3:      return "#F4FDD9"
    if index == 4:      return "#AA767C"
    if index == 5:      return "#D6A184"

def GetStackColorTauPOG(index):
    if index == 0:      return "#18252a"
    if index == 1:      return "#de5a6a"
    if index == 2:      return "#9999cc"
    if index == 3:      return "#4496c8"
    if index == 4:      return "#e5b7e5"
    if index == 5:      return "#ffcc66"

def GetStackColorTauPOGbyName(name):
    if 'VVV' in name:           return "#87F1FF"
    elif 'VV' in name:            return "#18252a"
    elif 'TT' in name:          return "#de5a6a"
    elif 'ST' in name:          return "#9999cc"
    elif 'WJets' in name:       return "#4496c8"
    elif 'QCD' in name:         return "#e5b7e5"
    elif 'DY' in name:          return "#ffcc66"
#    elif 'H' in name:           return "#87F1FF"
    else:                       return "#F4F1BB"

def GetHistColor(index):        
    if index == 0:      return "#000000"
    if index == 1:      return "#000075"
    if index == 2:      return "#800000"
    if index == 3:      return "#f58231"
    if index == 4:      return "#3cb44d"
    if index == 5:      return "#ffe119"
    if index == 6:      return "#87F1FF"
    if index == 7:      return "#F4F1BB"

def GetLineColor(index):
    if index == 0:      return "#000000"
    if index == 1:      return "#e6194B"
    if index == 2:      return "#4363d8"
    if index == 3:      return "#B9BAA3"
    if index == 4:      return "#685762"
    if index == 5:      return "#E8C547"

def GetMarker(index):
    if index == 0:      return 20
    if index == 1:      return 21
    if index == 2:      return 22
    if index == 3:      return 23
    if index == 4:      return 24
    if index == 5:      return 25

def GetOverallMaximum(hist):
    currentMax = 0
    for h in hist:
        if h.GetMaximum() > currentMax :    currentMax = h.GetMaximum()
    return currentMax

def GetOverallMinimum(hist, zero_not_allowed = False):
    currentMin = 99999
    for i, h in enumerate(hist):
        if h.GetMinimum() == -1 and zero_not_allowed: continue
        if h.GetMinimum() < currentMin:    currentMin = h.GetMinimum()
    return currentMin

def GetNestedMin(arr):
    localMin = []
    for a in arr:
        localMin.append(np.min(a))
    return np.min(localMin)

def GetNestedMax(arr):
    localMax = []
    for a in arr:
        localMax.append(np.max(a))
    return np.max(localMax)

def GetXMin(graphs):
    xmin = 99999.
    for graph in graphs:
        if TMath.MinElement(graph.GetN(),graph.GetX()) < xmin:
            xmin = TMath.MinElement(graph.GetN(),graph.GetX())
    return xmin

def GetXMax(graphs):
    xmax = -99999.
    for graph in graphs:
        if TMath.MaxElement(graph.GetN(),graph.GetX()) > xmax:
            xmax = TMath.MaxElement(graph.GetN(),graph.GetX())
    return xmax

def GetYMin(graphs):
    ymin = 99999.
    for graph in graphs:
        if TMath.MinElement(graph.GetN(),graph.GetY()) < ymin:
            ymin = TMath.MinElement(graph.GetN(),graph.GetY())
    return ymin

def GetYMax(graphs):
    ymax = -99999.
    for graph in graphs:
        if TMath.MaxElement(graph.GetN(),graph.GetY()) > ymax:
            ymax = TMath.MaxElement(graph.GetN(),graph.GetY())
    return ymax

def orderHist(hist, names, lowestFirst = False):
    weight = -1.
    if lowestFirst: weight = 1.
    sof = [h.GetSumOfWeights()*weight for h in hist]
    return sortByOtherList(hist, sof), sortByOtherList(names, sof)

def extraTextFormat(text, xpos = None, ypos = None, textsize = None, align = 12):
    return [text, xpos, ypos, textsize, align]

def DrawExtraText(pad, additionalInformation):
    pad.cd()
    #Write extra text
    if additionalInformation is not None:
        lastYpos = 0.8
        lastCorrectedYpos = None
        lastXpos = 0.2
        extraText = TLatex()
        for info in additionalInformation:
            try :
                extraTextString = info[0]
                extraTextXpos = info[1]
                extraTextYpos = info[2]
                extraTextSize = info[3]
            except:
                print("Wrong Format for additionalInformation. Stopping")
                pass

            if extraTextSize is None:
                extraTextSize = 0.03
            correction_term = extraTextSize*(5./3.)

            if extraTextXpos is None:
                if extraTextYpos is None:
                    extraTextXpos = lastXpos
                else:
                    extraTextXpos = 0.2
            if extraTextYpos is None:
                if lastYpos is None:
                    extraTextYpos = lastCorrectedYpos - correction_term
                else: extraTextYpos = lastYpos - correction_term
            
            extraText.SetNDC()
            extraText.SetTextAlign(info[4])
            extraText.SetTextSize(extraTextSize)
            extraText.DrawLatex(extraTextXpos, extraTextYpos, extraTextString)
            
            lastXpos = extraTextXpos
            lastYpos = info[2]
            lastCorrectedYpos = extraTextYpos

    pad.Update()

def getUnit(x):
    if x.find('[') == -1:
        return ''
    else:
        return x[x.find('[')+len('['):x.rfind(']')]

def makeList(item):
    if not isinstance(item, (list,)):
        item = [item]
    return item

def savePlots(Canv, destination):

    Canv.SaveAs(destination + ".pdf")
    Canv.SaveAs(destination + ".png")
    Canv.SaveAs(destination + ".root")


def plotClosure(observed, predicted, xtitle, ytitle, DataName, destination, yLog = False, additionalInfo = None, denominator_shape = None):
    GeneralSettings()
    
    Canv = TCanvas("Canv"+destination, "Canv"+destination, 1000, 1000)
 
    observed.SetMarkerColor(ROOT.kBlack)
    observed.SetLineColor(ROOT.kBlack)    
    observed.SetMarkerStyle(20)    
   
    predicted.SetFillColor(TColor.GetColor('#3399ff'))
    predicted.SetLineColor(TColor.GetColor('#3399ff'))
 
    #First pad
    plotpad = TPad("plotpad", "plotpad", 0, .3, 1, 0.98)
    plotpad.SetBottomMargin(0.025)
    plotpad.Draw()
    plotpad.cd()

    predicted.GetXaxis().SetRange(1, predicted.GetXaxis().FindBin(100))
    predicted.Draw("Hist")                                                            #Draw before using GetHistogram, see https://root-forum.cern.ch/t/thstack-gethistogram-null-pointer-error/12892/4
    title = " ; ; "+ytitle+" / " +str(predicted.GetBinWidth(1)) +' '+ getUnit(xtitle)
    predicted.SetTitle(title)
    predicted.GetXaxis().SetLabelOffset(9999999)
    
    predictedError = predicted.Clone("TotalPredictedError")
    predictedError.SetFillStyle(3013)
    predictedError.SetFillColor(ROOT.kGray+2)
    predictedError.SetMarkerStyle(0)
    predictedError.Draw("E2 Same")

    #Set range
    overallMin = GetOverallMinimum([observed, predicted], yLog)
    overallMax = GetOverallMaximum([observed, predicted])

    if yLog:
        plotpad.SetLogy()
        predicted.SetMinimum(0.3*overallMin)    
        predicted.SetMaximum(10*overallMax)    
    else:
        predicted.SetMinimum(0.5*overallMin)    
        predicted.SetMaximum(1.4*overallMax)    
    observed.Draw("EPSame")
     
    if denominator_shape:
        denominator_shape.SetMarkerColor(ROOT.kGray+2)
        denominator_shape.SetLineColor(ROOT.kGray+2)
        denominator_shape.SetMarkerStyle(20)
        denominator_shape.Scale(observed.GetSumOfWeights()/denominator_shape.GetSumOfWeights())
        denominator_shape.Draw('EPSame')

    #Create Legend
    legend = TLegend(0.4, .7, .6, .9)
    legend.AddEntry(observed, DataName + ' (triggered)')  
    legend.AddEntry(predicted, DataName + ' (weighted)')  
    if denominator_shape:       legend.AddEntry(denominator_shape, 'normalized baseline distribution')
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)

    legend.Draw()

    #Draw extra text
    DrawExtraText(plotpad, additionalInfo)

    #Return to canvas
    Canv.cd()

    #Second pad
    ratiopad = TPad("ratiopad", "ratiopad", 0, 0.05, 1, .3)
    ratiopad.SetTopMargin(0.05)
    ratiopad.SetBottomMargin(0.25)
    ratiopad.Draw()
    ratiopad.cd()
   
    ratio = observed.Clone('ratio')
    ratio.Divide(predicted)

    predictedErrorRatio = predictedError.Clone("PredictedErrorsRatio")
    for b in xrange(predictedErrorRatio.GetNbinsX()+1):
        if(predictedErrorRatio.GetBinContent(b) != 0):
            predictedErrorRatio.SetBinError(b, predictedErrorRatio.GetBinError(b)/predictedErrorRatio.GetBinContent(b))
            predictedErrorRatio.SetBinContent(b, 1.)
        else:
            predictedErrorRatio.SetBinContent(b, 0)
    predictedErrorRatio.SetMinimum(0.5)
    predictedErrorRatio.SetMaximum(1.5)
        

    #Set Style for bottom plot
    ratio.SetTitle(";" + xtitle + "; triggered/weighted")
    ratio.GetXaxis().SetTitleSize(.12)
    ratio.GetYaxis().SetTitleSize(.12)
    ratio.GetYaxis().SetTitleOffset(.6)
    ratio.GetXaxis().SetLabelSize(.12)
    ratio.GetYaxis().SetLabelSize(.12)
    ratio.SetMinimum(0.5)
    ratio.SetMaximum(1.5)
    ratio.Draw("EP")
    predictedErrorRatio.Draw("E2 same")
    ratio.Draw("EPsame")

    #Draw a guide for the eye
    line = TLine(ratio.GetXaxis().GetXmin(),1,ratio.GetXaxis().GetXmax(),1)
    line.SetLineColor(ROOT.kRed)
    line.SetLineWidth(1)
    line.SetLineStyle(1)
    line.Draw("Same")
   
    #Throw CMs lumi at it
    # cl.CMS_lumi(Canv, 4, 11, 'Preliminary', True)
 
    #Save everything
    savePlots(Canv, destination)
    ROOT.SetOwnership(Canv,False)               #https://root-forum.cern.ch/t/tlatex-crashing-in-pyroot-after-many-uses/21638/4
    del Canv
    return

def draw2DHist(h, xlabel, ylabel, output, option="ETextColz", x_log = False, y_log = False, extraText=None, axis_max=None):
    
    tdr.setTDRStyle()
    GeneralSettings("4.2f")
    gStyle.SetPalette(ROOT.kIsland)
    gStyle.SetPadRightMargin(0.15)
     
    #Create Canvas
    Canv = TCanvas("Canv"+output, "Canv"+output, 1000, 1000)
    if x_log:
        Canv.SetLogx()
    if y_log:
        Canv.SetLogy()

    h.SetTitle(';'+xlabel+';'+ylabel) 
    h.SetMarkerSize(0.6)
    if axis_max:
        h.GetXaxis().SetRange(1, h.GetXaxis().FindBin(axis_max)) 
        h.GetYaxis().SetRange(1, h.GetXaxis().FindBin(axis_max)) 
    h.Draw(option)

    if extraText:
        DrawExtraText(Canv, extraText)
        
    #Throw CMs lumi at it
    # cl.CMS_lumi(Canv, 4, 0, 'Preliminary', True)

    #Save everything
    savePlots(Canv, output)

    ROOT.SetOwnership(Canv,False)               #https://root-forum.cern.ch/t/tlatex-crashing-in-pyroot-after-many-uses/21638/4
    del Canv
    return

    
def draw1DHist(h, xlabel, ylabel, output, legend = None, x_log = False, y_log = False, extraText=None):
    
    tdr.setTDRStyle()
    GeneralSettings("4.2f")

    h = makeList(h)    
 
    #Create Canvas
    Canv = TCanvas("Canv"+output, "Canv"+output, 1000, 1000)
    if x_log:
        Canv.SetLogx()
    if y_log:
        Canv.SetLogy()

    overallMax = GetOverallMaximum(h)
    h[0].SetTitle(';'+xlabel+';'+ylabel) 
    h[0].SetMaximum(1.2*overallMax)
    for i, hist in enumerate(h):
        hist.SetLineColor(ROOT.TColor.GetColor(GetLineColor(i)))
        hist.SetFillColor(0)
        hist.SetLineWidth(3)
        hist.SetMarkerStyle(8)
        hist.SetMarkerColor(ROOT.TColor.GetColor(GetLineColor(i)))
        if(i == 0) :
            hist.Draw("EP")
        else:
            hist.Draw("EPSAME")

    if extraText:
        DrawExtraText(Canv, extraText)
       
    if legend is not None:
        tlegend = ROOT.TLegend(0.5, .8, .9, .9)
        for hist, n in zip(h, legend):
            tlegend.AddEntry(hist, n)
        tlegend.SetFillStyle(0)
        tlegend.SetBorderSize(0)
        tlegend.Draw()
 
    #Throw CMs lumi at it
    # cl.CMS_lumi(Canv, 4, 0, 'Preliminary', True)

    #Save everything
    savePlots(Canv, output)

    ROOT.SetOwnership(Canv,False)               #https://root-forum.cern.ch/t/tlatex-crashing-in-pyroot-after-many-uses/21638/4
    del Canv
    return

    



