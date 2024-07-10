# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import re, os
import numpy as np
import datetime
import math
import ROOT
import ctypes
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaletteAxis, TPaveText, TPython, TMath, TF1, TLine, TPython, TArc
from ROOT import gStyle, gROOT, gSystem, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kViolet
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);
from FomattingMaterialBudget import FrameSettings, ALICEtext, RatioLegendSettings, FrameSettingsRatio, FrameSettings2D, ALICEtext2D

#________________________________________________
def make_common_style(g1,marker,size,color,width=1,fill=0):
    g1.SetMarkerStyle(marker);
    g1.SetMarkerColor(color);
    g1.SetMarkerSize(size);
    g1.SetLineColor(color);
    g1.SetLineWidth(width);
    g1.SetFillColor(color);
    g1.SetFillStyle(fill);

#________________________________________________
def draw_AP_without_cuts(filename_data, suffix, folder, date):
    rootfile_data = TFile.Open(filename_data, "READ");
    rootdir_data = rootfile_data.Get("pcm-qc");
    dir_v0 = rootdir_data.Get("V0")


    h2_data = dir_v0.Get("hAPplot");
    h2_data.Sumw2();
    h2_data.SetXTitle("#it{#alpha} = (p_{L}^{-}-p_{L}^{+})/(p_{L}^{-}+p_{L}^{+})")
    h2_data.SetYTitle("q_{T} (GeV/c)")
    h2_data.GetXaxis().SetTitleOffset(1.5 )
    h2_data.GetZaxis().SetTitleOffset(1.9);
    h2_data.SetDirectory(0);
    h2_data.GetXaxis().SetRangeUser(-100,100);
    h2_data.GetYaxis().SetRangeUser(-100,100);
    ROOT.SetOwnership(h2_data, False);

    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.Divide(1,2,1e-5,1e-5); 
    c1.SetMargin(0.13,0.13,0.13,0.13);
    c1.SetTicks(1,1);

    frame1 = c1.DrawFrame(-1, 0, 1, 0.25);
    FrameSettings2D(frame1)

    gPad.SetLogz()
    h2_data.Draw("COLZ");

    txt = TPaveText(0.0,0.95,1.0,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(22);#centered,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("Amenteros-Podolanski plot")
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    txt = TPaveText(0.52,0.77,0.78,0.85,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(1001);
    txt.SetBorderSize(0);
    # txt.SetTextColor(kRed)
    txt.SetTextAlign(32);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.03);
    # txt.AddText("ALICE");
    txt.AddText("this thesis")
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    outname = os.path.join(folder, "{0}_material_budget_AP_plot_without_cuts_{1}.pdf".format(date, suffix))
    c1.SaveAs(outname);

if __name__ == "__main__":
    filename_data = "/Users/juliaschlagel/Analysis240411/Analysis/HLtrains/data/AnalysisResults_LHC22o_full_statistics.root"
    filename_mc = "/Users/juliaschlagel/Analysis240411/Analysis/HLtrains/MC/AnalysisResults_LHC24b1_full_statistics.root"
    cutname = "qc"
    suffix = "AnyTrack" 
    folder = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/cuts"
    os.makedirs(folder, exist_ok=True);
    date= "this_thesis"
    

    period_data     = "LHC22o";
    period_mc       = "LHC241b";
    
    draw_AP_without_cuts(filename_data, suffix, folder, date);