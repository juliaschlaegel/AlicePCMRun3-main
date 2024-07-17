# This code was written by Julia Schlägel (July 2024)

import ROOT
import yaml
from ctypes import *
import numpy as np
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullCross, kFullDiamond, kOpenSquare, kOpenTriangleUp, kOpenCircle, kFullCrossX

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
from utility import Utility



class HistogramAcceptance:
    def __init__(self, filename, config, decay): 
        # Open wanted data and load histogramms
        self.decay = decay
        self.filename = TFile.Open(filename, "READ")
        dir_name1 = "pi0eta-to-gammagamma-mc/Generated"
        if self.decay == "pcm":
            dir_name2 = "PCMPCM"
        else:
            dir_name2 = "PCMDalitzEE"
        self.dir = self.filename.Get(dir_name1).FindObject(dir_name2)
        self.config = config    

    def calc_acceptance_vs_pt(self, name_plot, save, utils):
        #rebin histograms
        h1 = utils.rebin_hist("hPt_Pi0_Acc")
        utils.scale_by_Nev(h1)
        print("First Histogram rebinned")
        h2 = utils.rebin_hist("hPt_Pi0")
        utils.scale_by_Nev(h2)
        print("second histogram rebinned")
        if h1 is None:
            print(f"Failed to load histogram: hPt_Pi0")
    #get canvas and pad to plot the acceptance in 
        canv = TCanvas("acc1", "acc1", 0, 0, 900, 900)
        pad = canv.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()

        h_acc = h1.Clone("h_A")   
        
        # calculate acceptance by dividing the events with daughters in acceptance by all events
        h_acc.Divide(h1, h2, 1, 1 )
        # Define marker settings
        h_acc.SetFillColor(kCyan+1)
        h_acc.SetMarkerStyle(kFullCross)
        h_acc.SetMarkerColor(kBlue)
        h_acc.SetLineColor(kBlue)

        #Draw points in the histogramm
        h_acc.Draw("Esame")
       
        # Define y-axis settings
        y = h_acc.GetYaxis() 
        y.SetTitle("Acceptance") 
        y.SetTitleSize(0.048)
        y.SetTitleFont(42)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)

        #Define x-axis settings
        x = h_acc.GetXaxis()
        x.SetTitle("p_{T} [GeV/c]")
        x.SetTitleSize(0.048)
        x.SetTitleFont(42)
        x.SetTitleOffset(0.9)
        x.SetLabelFont(42)
        x.SetLabelSize(0.035)
        x.SetRangeUser(0,14)

        # Define legend settings
        leg = TLegend(0.2, 0.8, 0.4, 0.75)
        leg.SetBorderSize(0)
        leg.SetFillColor(kWhite)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.03)
        if self.decay == "pcm":
            leg.AddEntry(h_acc, "\pi^{0} \\rightarrow \gamma \gamma", "LP")
        else:
            leg.AddEntry(h_acc, "\pi^{0} \\rightarrow e^{+} e^{-} \gamma", "LP")
        leg.Draw("")
        ROOT.SetOwnership(leg,False)
               
        # Add text 
        txt = TPaveText(0.3,0.85,0.3,0.85,"NDC")
        txt.SetFillColor(kWhite)
        txt.SetFillStyle(0)
        txt.SetBorderSize(0)
        txt.SetTextAlign(33);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.02)
        txt.AddText("this thesis")
        txt.Draw()
        ROOT.SetOwnership(txt,False)

        txt3 = TPaveText(0.35,0.82,0.35,0.82,"NDC"); 
        txt3.SetFillColor(kWhite)
        txt3.SetFillStyle(0)
        txt3.SetBorderSize(0)
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.02)
        txt3.AddText("pp, #sqrt{s} = 13.6TeV")
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)


        #save histogram
        canv.Modified()
        canv.Update()
        ROOT.SetOwnership(canv, False)
        if save == "true":
            canv.SaveAs(name_plot.Data())
        if save == "loc":
            canv.SaveAs("Acceptance_with_ut")
        print("acc hist type:", type(h_acc))
        return h_acc
        


