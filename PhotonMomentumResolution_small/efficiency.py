# This code was written by Julia Schlägel (July 2023)

import ROOT
import yaml
import math
import os, sys, shutil
import datetime
from ctypes import *
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullCross, kFullDiamond, kOpenSquare, kOpenTriangleUp, kOpenCircle, kFullCrossX
import numpy as np
from acceptance import HistogramAcceptance
from utility import Utility


gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

class Efficiency:
    def __init__(self, decay): 
        self.decay = decay
       
    def calc_efficiency_vs_pt(self, name_plot, save, raw_mc, utils): 
        h1 = raw_mc
        print("First Histogram rebinned")
        h2 = utils.rebin_hist("hPt_Pi0_Acc")
        utils.scale_by_Nev(h2)
        print("second histogram rebinned")
        if h1 is None:
            print(f"Failed to load histogram: hPt_Pi0")
    #get canvas and pad to plot the acceptance in 
        can = TCanvas("eff1", "eff1", 0, 0, 900, 900)
        pad = can.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()
        pad.SetLogx()
        

        h_eff = h1.Clone("h_E")   
        
        h_eff.Divide(h2)
        # Define marker settings
        h_eff.SetFillColor(kCyan+1)
        h_eff.SetMarkerStyle(kFullCross)
        h_eff.SetMarkerColor(kBlue)
        h_eff.SetMarkerSize(0.85)
        h_eff.SetLineColor(kBlue)

        #Draw points in the histogramm
        h_eff.Draw("Esame")
       
        # Define y-axis settings
        y = h_eff.GetYaxis() 
        y.SetTitle("Efficiency") 
        y.SetTitleSize(0.048)
        y.SetTitleFont(42)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)

        #Define x-axis settings
        x = h_eff.GetXaxis()
        x.SetTitle("p_{T} [GeV/c]")
        x.SetTitleSize(0.048)
        x.SetTitleFont(42)
        x.SetTitleOffset(0.9)
        x.SetLabelFont(42)
        x.SetLabelSize(0.035)
        

        # Define legend settings
        leg = TLegend(0.2, 0.8, 0.4, 0.75)
        leg.SetBorderSize(0)
        leg.SetFillColor(kWhite)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.03)
        if self.decay == "pcm":
            leg.AddEntry(h_eff, "\pi^{0} \\rightarrow \gamma \gamma", "LP")
        else:
            leg.AddEntry(h_eff, "\pi^{0} \\rightarrow e^{+} e^{-} \gamma", "LP")
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

        txt3 = TPaveText(0.35,0.82,0.35,0.82,"NDC"); # txt3 = TPaveText(0.4,0.75,0.4,0.90,"NDC");
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
        can.Modified()
        can.Update()
        ROOT.SetOwnership(can, False)
        if save == "true":
            can.SaveAs(name_plot.Data())
        
        print("eff hist type:", type(h_eff))
        return h_eff
    
        