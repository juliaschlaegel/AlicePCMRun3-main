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
    def __init__(self, meson, filename, config, decay, typ): #, hist1_name, hist2_name):
        # Open wanted data and load histogramms
        self.decay = decay
        self.typ = typ
        self.filename = TFile.Open(filename, "READ")
        self.config = config
        self.meson = meson

    def get_hist (self, hist_name):
        if self.decay == "dalitz":
            #if self.typ == "mc":
            dir_name1 = "pi0eta-to-gammagamma-mc-pcmdalitzee/Generated"
            #if self.typ == "data":
            #    dir_name1 = ""
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

    def get_Nev (self): #dir_name2):
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
    
        # print(Nev)
    def get_bin_var(self):
         #with open(self.config, "r") as file:
            #config_data = yaml.safe_load(file)
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
             #bin_width = input_hist.GetBinWidth(bin)
             # print("bin width: ", bin_width)
             bin_error = input_hist.GetBinError(bin)

             new_bin_content = bin_content / (bin_center) #* Nev* deltay*2*np.pi)
             input_hist.SetBinContent(bin, new_bin_content)
             new_bin_error = bin_error / (bin_center)
             input_hist.SetBinError(bin, new_bin_error)
        return input_hist
    
    def scale_by_BR(self, input_hist):
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
    
        for bin in range(0, len(bin_var) - 1):
            b1 = bin_var[bin]
            b2 = bin_var[bin + 1]
            # print("b1: ", b1, "b2: ", b2)

            bin_b1 = hist.GetXaxis().FindBin(b1 + 1e-6)
            bin_b2 = hist.GetXaxis().FindBin(b2 - 1e-6)
            # print("bin_b1: ", bin_b1, "bin_b2: ", bin_b2)

            error = c_double(0.0)
            content = hist.IntegralAndError(bin_b1, bin_b2, error, "")
            bin_center = (b1 + b2) / 2  # Take the center of the bin to fill
            bin_width = b2-b1
            new_content = content/ bin_width
            new_hist.Fill(bin_center, new_content)  # Filling with this center of each bin
            
            new_hist.SetBinError(bin+1, error.value / bin_width)  # calculate error to the right bin#

        #NoE = self.get_Nev()
        #new_hist.Scale(1 / NoE)
        if new_hist is None:
            print(f"Failed to load histogram: {hist_name}")
        # print("hist neu type:", type(new_hist))
        # print(f"Histogram {hist_name} loaded with {new_hist.GetEntries()} entries.")
        # return new_hist

        # Create canvas and draw histogram
        canvas = TCanvas("acc", "acc", 0, 0, 900, 900)
        pad = canvas.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        new_hist.Draw("ESame")
        ROOT.SetOwnership(canvas, False)
        # canvas.SaveAs("rebinning_test1")

        return new_hist