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

class Utility:
    def __init__(self, filename, config, decay, typ): 
        self.decay = decay
        self.typ = typ
        self.filename = TFile.Open(filename, "READ")
        self.config = config

    def get_hist (self, hist_name):
        dir_name1 = "pi0eta-to-gammagamma-mc/Generated"
        
        if self.decay == "dalitz":
            dir_name2 = "PCMDalitzEE"
        if self.decay == "pcm":
            dir_name2 = "PCMPCM"
        self.hist_name = self.filename.Get(dir_name1).FindObject(dir_name2).FindObject(hist_name)
        if self.hist_name:
            print("histogram found")
        if not self.hist_name:
            raise ValueError("Histogram not found")
        return self.hist_name
    
    def get_Nev (self): #dir_name2):
        if self.typ == "mc":
            direv = self.filename.Get("pi0eta-to-gammagamma-mc/Event")
        if self.typ == "data":
            direv = self.filename.Get("pi0eta-to-gammagamma/Event")
        if self.decay == "pcm":
            dir_name2 = "PCMPCM"
        else:
            dir_name2 = "PCMDalitzEE"
        events = direv.FindObject(dir_name2).FindObject("after").FindObject("hCollisionCounter")
        self.Nev = events.GetBinContent(9)
        print("Nev :", self.Nev)
        return self.Nev
    
    def get_bin_var(self):
            self.bin_var = self.config["common"]["pt_bin"]
            print(self.bin_var)
            return self.bin_var
    
    def scale_by_Nev(self, input_hist):
        Nev = self.get_Nev()
        input_hist.Scale(1/Nev)
        return input_hist
    
    def scale_by_pt(self, input_hist):
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



    def rebin_hist(self, hist_name):
        bin_var = self.get_bin_var()
        hist = self.get_hist(hist_name)
        new_hist = TH1D("p_T1", "p_T1", len(bin_var) - 1, np.array(bin_var))
    
        for bin in range(0, len(bin_var) - 1):
            b1 = bin_var[bin]
            b2 = bin_var[bin + 1]
            

            bin_b1 = hist.GetXaxis().FindBin(b1 + 1e-6)
            bin_b2 = hist.GetXaxis().FindBin(b2 - 1e-6)
            

            error = c_double(0.0)
            content = hist.IntegralAndError(bin_b1, bin_b2, error, "")
            bin_center = (b1 + b2) / 2  # Take the center of the bin to fill
            bin_width = b2-b1
            new_content = content/ bin_width
            new_hist.Fill(bin_center, new_content)  # Filling with this center of each bin
            
            new_hist.SetBinError(bin+1, error.value / bin_width)  # calculate error to the right bin#

        
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