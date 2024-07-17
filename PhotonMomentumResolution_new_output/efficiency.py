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
from array import array
from acceptance import HistogramAcceptance
from utility import Utility
# from PhotonMomentumResolution import FitInvMassForPt 
#from FitInvMassForPt import PairAnalyzer
#from PlotRawYield import PlotRawYieldInvMass




gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

class Efficiency:
    def __init__(self, decay, meson): 
        self.decay = decay
        self.meson = meson
    
    def eff_and_acc_vs_pt(self, scale_BR, raw_mc, utils, name_plot, save):
        #Get wanted histograms:
        h1 = raw_mc
        print("First Histogram rebinned")
        
        h2 = utils.rebin_hist("hPt")
        utils.scale_by_Nev(h2)
        pt = utils.get_bin_var()
        print("####### pt range used: ", pt)
        
        print("second histogram rebinned")
        if h1 is None:
            print(f"Failed to load histogram: hPt_Pi0")
        

    #get canvas and pad to plot the acceptance in 
        can = TCanvas("eff_acc", "eff_acc", 0, 0, 900, 900)
        pad = can.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()
        #pad.SetLogx()
        

        h_eff_acc = h1.Clone("h_E_A")   
    
        h_eff_acc.Divide(h2)
        if scale_BR == "true":
            utils.scale_by_BR(h_eff_acc)
        
        for i in range(0, len(pt)):
            print("###########Efficiency, values and Errors" , pt[i] ,"########")
            print("Value of raw mc: ", h1.GetBinContent(i), "±", h1.GetBinError(i))
            print("Value of hPt all: ", h2.GetBinContent(i), "±", h2.GetBinError(i))
            print("value of eff x acc x BR: ", h_eff_acc.GetBinContent(i), "±", h_eff_acc.GetBinError(i))
            if h_eff_acc.GetBinError(i) != 0:
                rel = h_eff_acc.GetBinContent(i) / h_eff_acc.GetBinError(i)
                print("relaive error: ", rel)
        
        # Define marker settings
        h_eff_acc.SetFillColor(kCyan+1)
        if self.decay == "pcm":
            h_eff_acc.SetMarkerStyle(34)
        if self.decay == "dalitz":
            h_eff_acc.SetMarkerStyle(20)
        #h_eff_acc.SetMarkerStyle(kFullCross)
        h_eff_acc.SetMarkerColor(kBlue)
        #h_eff_acc.SetMarkerSize(0.85)
        h_eff_acc.SetMarkerSize(1.5)
        h_eff_acc.SetLineColor(kBlue)

        #Draw points in the histogramm
        h_eff_acc.Draw("Esame")
       
        # Define y-axis settings
        y = h_eff_acc.GetYaxis() 
        y.SetTitle("#varepsilon x A x BR") 
        y.SetTitleSize(0.025)
        y.SetTitleFont(42)
        y.SetLabelFont(42)
        y.SetLabelSize(0.025)
        
        #Define x-axis settings
        x = h_eff_acc.GetXaxis()
        x.SetTitle("p_{T} [GeV/c]")
        x.SetTitleSize(0.025)
        x.SetTitleFont(42)
        x.SetTitleOffset(0.9)
        x.SetLabelFont(42)
        x.SetLabelSize(0.02)
        #x.SetRangeUser(0., 15)
        

        # Define legend settings
        leg = TLegend(0.65, 0.3, 0.9, 0.25)
        leg.SetBorderSize(0)
        leg.SetFillColor(kWhite)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.03)
        if self.meson == "pi0":
            if self.decay == "pcm":
                leg.AddEntry(h_eff_acc, "\pi^{0} \\rightarrow \gamma \gamma", "LP")
            else:
                leg.AddEntry(h_eff_acc, "\pi^{0} \\rightarrow e^{+} e^{-} \gamma", "LP")
        if self.meson == "eta":
            if self.decay == "pcm":
                leg.AddEntry(h_eff_acc, "\eta \\rightarrow \gamma \gamma", "LP")
            else:
                leg.AddEntry(h_eff_acc, "\eta \\rightarrow e^{+} e^{-} \gamma", "LP")
        leg.Draw("")
        ROOT.SetOwnership(leg,False)
               
        # Add text 
        txt = TPaveText(0.8,0.25,0.8,0.25,"NDC")
        txt.SetFillColor(kWhite)
        txt.SetFillStyle(0)
        txt.SetBorderSize(0)
        txt.SetTextAlign(33);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.02)
        txt.AddText("this thesis")
        txt.Draw()
        ROOT.SetOwnership(txt,False)

        txt3 = TPaveText(0.85,0.22,0.85,0.22,"NDC"); # txt3 = TPaveText(0.4,0.75,0.4,0.90,"NDC");
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
        print("eff and acc hist type:", type(h_eff_acc))
        return h_eff_acc
    
    def teff_and_acc_vs_pt(self, raw_mc, utils, name_plot, save):
        #h1_b = raw_mc
        h1 = utils.scale_by_BR(raw_mc)
        h2 = utils.rebin_hist("hPt")
        utils.scale_by_Nev(h2)

        shift = 0.1
        def shift_histo(hist, shift):
            nbins = hist.GetNbinsX()
            xbins = [hist.GetBinLowEdge(i)+ shift for i in range(1, nbins+2)]
            new_hist = ROOT.TH1F(hist.GetName() + "_shifted", hist.GetTitle(), nbins, array("d", xbins))
            for i in range(1, nbins + 1):
                new_hist.SetBinContent(i, hist.GetBinContent(i))
                new_hist.SetBinError(i, hist.GetBinError(i))
            return new_hist
        shifted_h1 = shift_histo(h1, shift)
        shifted_h2 = shift_histo(h2, shift)

        
        h_eff_acc = ROOT.TEfficiency(shifted_h1, shifted_h2)
        h_eff_acc.SetStatisticOption(ROOT.TEfficiency.kFCP)
        #h_eff_acc.SetStatisticOption(ROOT.TEfficiency.kBBayesian)

        c = ROOT.TCanvas("teff", "teff", 0, 0, 900, 900)
        pad = c.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()

        h_eff_acc.Draw()
        h_eff_acc.SetTitle("TEfficiency; p_{T} [GeV/c] ; Teff x A x BR")
        c.Update()
        



        #utils.scale_by_BR(h_eff_acc)
        
        # Define marker settings
        h_eff_acc.SetFillColor(kCyan+1)
        h_eff_acc.SetMarkerStyle(kFullCross)
        h_eff_acc.SetMarkerColor(kBlue)
        #h_eff_acc.SetMarkerSize(0.85)
        h_eff_acc.SetMarkerSize(1.5)
        h_eff_acc.SetLineColor(kBlue)

        

        # Define legend settings
        leg = TLegend(0.2, 0.8, 0.4, 0.75)
        leg.SetBorderSize(0)
        leg.SetFillColor(kWhite)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.03)
        if self.meson == "pi0":
            if self.decay == "pcm":
                leg.AddEntry(h_eff_acc, "\pi^{0} \\rightarrow \gamma \gamma", "LP")
            else:
                leg.AddEntry(h_eff_acc, "\pi^{0} \\rightarrow e^{+} e^{-} \gamma", "LP")
        if self.meson == "eta":
            if self.decay == "pcm":
                leg.AddEntry(h_eff_acc, "\eta \\rightarrow \gamma \gamma", "LP")
            else:
                leg.AddEntry(h_eff_acc, "\eta \\rightarrow e^{+} e^{-} \gamma", "LP")
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


        # save histogram
        c.Modified()
        c.Update()
        ROOT.SetOwnership(c, False)
        if save == "true":
            c.SaveAs(name_plot.Data())
        print("eff and acc hist type:", type(h_eff_acc))
        return h_eff_acc