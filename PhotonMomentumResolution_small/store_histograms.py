# This code was written by Julia Schlägel (July 2024)

import ROOT
import ctypes
from ctypes import *
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullCross, kFullDiamond, kOpenSquare, kOpenTriangleUp, kOpenCircle, kFullCrossX
import numpy as np
from acceptance import HistogramAcceptance
from efficiency import Efficiency
from corrected_yield import CorrectedYield

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)


class Save_Hist:
    def __init__(self, filename):
        self.filename = filename
        self.histograms = {}

    def add_histogram(self, folder, histogram, new_name=None):
        if new_name is not None:
            histogram.SetName(new_name)
        
        if folder not in self.histograms:
            self.histograms[folder] = []
        self.histograms[folder].append(histogram)

    def save_histograms(self):
        try:
            output_file = ROOT.TFile(self.filename, "RECREATE")
            for folder, hist_list in self.histograms.items():
                dir = output_file.mkdir(folder)
                dir.cd()
                for hist in hist_list:
                    hist.Write()
            output_file.Close()
            print(f"Histograms successfully stored in {self.filename}")
        except Exception as e:
            print(f"Error by storing of the histos:{e}")

    
