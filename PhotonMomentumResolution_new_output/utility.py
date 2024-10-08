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

class Utility: #in this class all central things are done that are needed in multiple macros like loading a file from a rootfile,...
    def __init__(self, meson, filename, config, decay, typ): 
        self.decay = decay
        self.typ = typ
        self.filename = TFile.Open(filename, "READ")
        self.config = config
        self.meson = meson

    def get_hist (self, hist_name):
        if self.decay == "dalitz":
            dir_name1 = "pi0eta-to-gammagamma-mc-pcmdalitzee/Generated"
        if self.decay == "pcm":
            dir_name1 = "pi0eta-to-gammagamma-mc-pcmpcm/Generated"
        if self.meson == "pi0":
            dir_name2 = "Pi0"
        if self.meson == "eta":
            dir_name2 = "Eta"
        dir = self.filename.Get(dir_name1).Get(dir_name2)
        self.hist_name = dir.Get(hist_name)
        if not self.hist_name:
            raise ValueError("Histogram not found")
        return self.hist_name

    def get_Nev (self): #function to get the number of events from the rootfile
        if self.decay == "dalitz":
            if self.typ =="mc":
                direv = self.filename.Get("pi0eta-to-gammagamma-mc-pcmdalitzee/Event")
            if self.typ == "data":
                direv = self.filename.Get("pi0eta-to-gammagamma-pcmdalitzee/Event")
        if self.decay == "pcm":
            if self.typ == "mc":
                direv = self.filename.Get("pi0eta-to-gammagamma-mc-pcmpcm/Event")
            if self.typ == "data":
                direv = self.filename.Get("pi0eta-to-gammagamma-pcmpcm/Event")
    
        events = direv.Get("after").Get("hCollisionCounter")
        self.Nev = events.GetBinContent(10)
        print("Nev :", self.Nev)
        return self.Nev
    
    def get_bin_var(self): #function to get the pT binning that was set in the config file
            self.bin_var = self.config["common"]["pt_bin"]
            print(self.bin_var)
            return self.bin_var
    
    def scale_by_Nev(self, input_hist): #function to divide a histogram by the number of events
        Nev = self.get_Nev()
        input_hist.Scale(1/Nev)
        return input_hist
    def scale_to_inv_yield(self, input_hist):
        #input_hist.Scale((57.8*1e-3))
        input_hist.Scale(59.4*1e-3)
        input_hist.Scale(1e12)
        return input_hist
    
    def scale_by_pt(self, input_hist): # function to divide a histogram by pT
        p_T = self.get_bin_var()
        for bin in range (1, len(p_T)):
             bin_content = input_hist.GetBinContent(bin)
             bin_center = input_hist.GetBinCenter(bin)
             bin_error = input_hist.GetBinError(bin)

             new_bin_content = bin_content / (bin_center) 
             input_hist.SetBinContent(bin, new_bin_content)
             new_bin_error = bin_error / (bin_center)
             input_hist.SetBinError(bin, new_bin_error)
        return input_hist
     
    def scale_by_BR(self, input_hist): # function to scale a histogram with the branching ratio
        if self.meson == "eta":
            if self.decay == "pcm":
                Br = 0.3936
            if self.decay == "dalitz":
                Br = 0.0069
        if self.meson == "pi0":
            if self.decay == "pcm":
                Br = 0.98823
            if self.decay == "dalitz":
                Br = 0.01174
        input_hist.Scale(1/Br)
        return input_hist



    def rebin_hist(self, hist_name):
        bin_var = self.get_bin_var()
        hist = self.get_hist(hist_name)
        new_hist = TH1D("p_T1", "p_T1", len(bin_var) - 1, np.array(bin_var))
        new_hist.Sumw2()
    
        for bin in range(0, len(bin_var) - 1):
            b1 = bin_var[bin]
            b2 = bin_var[bin + 1]
            bin_b1 = hist.GetXaxis().FindBin(b1 + 1e-6)
            bin_b2 = hist.GetXaxis().FindBin(b2 - 1e-6)
        
            error = c_double(0.0)
            # content = hist.IntegralAndError(bin_b1, bin_b2, error, "")
            content = 0
            print("bin1: ", bin_b1, bin_b2)
            for i in range(bin_b1, bin_b2 +1):
                content = content+ hist.GetBinContent(i)
            error = np.sqrt(content)

            bin_center = (b1 + b2) / 2  # Take the center of the bin to fill
            bin_width = b2-b1
            print("in rebin function error:", error)
            print("content: ", content)
            # 
            new_content = content / bin_width
            # new_hist.Fill(bin_center, new_content)  # Filling with this center of each bin
            new_hist.SetBinContent(new_hist.FindBin(bin_center), new_content)
            #new_hist.SetBinError(bin, error) # / bin_width)  # calculate error to the right bin#
            new_hist.SetBinError(new_hist.FindBin(bin_center), error)


        if new_hist is None:
            print(f"Failed to load histogram: {hist_name}")
        

        # Create canvas and draw histogram
        canvas = TCanvas("acc", "acc", 0, 0, 900, 900)
        pad = canvas.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        new_hist.Draw("ESame")
        ROOT.SetOwnership(canvas, False)
        

        return new_hist