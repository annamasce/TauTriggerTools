import ROOT
from ROOT import TLine, THStack, TColor, TCanvas, TLegend, gROOT, gStyle, TH1, TMultiGraph, TGraphErrors, TGraph, TLatex, TPad, TMath
import CMS_lumi as cl
import tdrstyle as tdr
import numpy as np
from helpers import makeDirIfNeeded, isTimeStampFormat, sortByOtherList
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

def getUnit(x):
    if x.find('[') == -1:
        return ''
    else:
        return x[x.find('[')+len('['):x.rfind(']')]

def savePlots(Canv, destination):
    destination_components = destination.split('/')
    cleaned_components = [x for x in destination_components if not isTimeStampFormat(x)]
    index_for_php = cleaned_components.index('Results')
    php_destination = '/user/lwezenbe/public_html/'
    php_destination += '/'.join(cleaned_components[index_for_php+1:])
    makeDirIfNeeded(php_destination.rsplit('/', 1)[0])    
    os.system('cp -rf /user/lwezenbe/private/PhD/index.php '+ php_destination.rsplit('/', 1)[0]+'/index.php')    

    print destination

    Canv.SaveAs(destination + ".pdf")
    Canv.SaveAs(destination + ".png")
    Canv.SaveAs(destination + ".root")

    #Clean out the php directory you want to write to if it is already filled, otherwise things go wrong with updating the file on the website
    #os.system("rm "+php_destination.rsplit('/')[0]+"/*")

    Canv.SaveAs(php_destination + ".pdf")
    Canv.SaveAs(php_destination + ".png")
    Canv.SaveAs(php_destination + ".root")

def plotClosure(observed, predicted, xtitle, ytitle, DataName, destination, yLog = False):
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

    predicted.Draw("EHist")                                                            #Draw before using GetHistogram, see https://root-forum.cern.ch/t/thstack-gethistogram-null-pointer-error/12892/4
    title = " ; ; "+ytitle+" / " +str(predicted.GetBinWidth(1)) +' '+ getUnit(xtitle)
    predicted.SetTitle(title)
    predicted.GetXaxis().SetLabelOffset(9999999)

    #Set range
    overallMin = GetOverallMinimum([observed, predicted], yLog)
    overallMax = GetOverallMaximum([observed, predicted])

    if yLog:
        plotpad.SetLogy()
        predicted.SetMinimum(0.3*overallMin)    
        predicted.SetMaximum(10*overallMax)    
    else:
        predicted.SetMinimum(0.5*overallMin)    
        predicted.SetMaximum(1.2*overallMax)    
    observed.Draw("EPSame")
     
    #Create Legend
    legend = TLegend(0.7, .7, .9, .9)
    legend.AddEntry(observed, DataName + ' (observed)')  
    legend.AddEntry(predicted, DataName + ' (predicted)')  
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)

    legend.Draw()

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

    #Set Style for bottom plot
    ratio.SetTitle(";" + xtitle + "; Obs./pred.")
    ratio.GetXaxis().SetTitleSize(.12)
    ratio.GetYaxis().SetTitleSize(.12)
    ratio.GetYaxis().SetTitleOffset(.6)
    ratio.GetXaxis().SetLabelSize(.12)
    ratio.GetYaxis().SetLabelSize(.12)
    ratio.SetMinimum(0.5)
    ratio.SetMaximum(1.5)
    ratio.Draw("EP")

    #Draw a guide for the eye
    line = TLine(ratio.GetXaxis().GetXmin(),1,ratio.GetXaxis().GetXmax(),1)
    line.SetLineColor(ROOT.kRed)
    line.SetLineWidth(1)
    line.SetLineStyle(1)
    line.Draw("Same")
   
    #Throw CMs lumi at it
    cl.CMS_lumi(Canv, 4, 11, 'Preliminary', True)
 
    #Save everything
    savePlots(Canv, destination)

def draw2DHist(h, xlabel, ylabel, output, option="ETextColz", x_log = False, y_log = False):
    
    tdr.setTDRStyle()
    GeneralSettings("4.3f")
    gStyle.SetPalette(ROOT.kIsland)
    gStyle.SetPadRightMargin(0.15)
     
    #Create Canvas
    Canv = TCanvas("Canv"+output, "Canv"+output, 1000, 1000)
    if x_log:
        Canv.SetLogx()
    if y_log:
        Canv.SetLogy()

    h.SetTitle(';'+xlabel+';'+ylabel) 
    h.Draw(option)
        
    #Throw CMs lumi at it
    cl.CMS_lumi(Canv, 4, 0, 'Preliminary', True)

    #Save everything
    savePlots(Canv, output)

    ROOT.SetOwnership(Canv,False)               #https://root-forum.cern.ch/t/tlatex-crashing-in-pyroot-after-many-uses/21638/4
    return

    



