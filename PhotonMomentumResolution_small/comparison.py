import ROOT
import ctypes
from ctypes import *
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullCross, kFullDiamond, kOpenSquare, kOpenTriangleUp, kOpenCircle, kFullCrossX
import numpy as np
import yaml
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TGaxis, TPad, TLatex, TBox, TString, TLine
# import matplotlib.pyplot as plt
# from hist import Hist
from HistoFormatting import FrameSettings
from PlotRawYield import PlotRawYieldInvMass
from acceptance import HistogramAcceptance
from efficiency import Efficiency
from corrected_yield import CorrectedYield
from utility import Utility

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)


class Compare:
    # Defining a class and load the scaling factors: The kind of meson wanted, the Number of Events,
    # the rapidity delta and the scaling factor
    def __init__(self) -> None:
        pass
    def __init__(self, filename_mc): #, filename, filename1, acc_class, eff_class, corr_class):

        self.filename_mc = TFile.Open(filename_mc, "READ")
    
    def plot_eff_dalitz_vs_pcm(self, filename_mc, config, raw_yield_mc_dalitz, raw_yield_mc_pcm, name_plot, save):
        utils_dalitz = Utility(filename_mc, config, "dalitz", "mc")
        utils_pcm = Utility(filename_mc, config, "pcm", "mc")
        # acc_class_dalitz = HistogramAcceptance(filename_mc, config, "dalitz")
        # acc_class_pcm = HistogramAcceptance(filename_mc, config, "pcm")
        eff_dalitz = Efficiency("dalitz")
        eff_pcm = Efficiency("pcm")
        h1 = eff_dalitz.calc_efficiency_vs_pt("Efficiency_dalitz", "false", raw_yield_mc_dalitz, utils_dalitz)
        h2 = eff_pcm.calc_efficiency_vs_pt("Efficiency_pcm","false", raw_yield_mc_pcm, utils_pcm)
        ratio = h1.Clone("Ratio_eff")
        ratio.Divide(h2)
        #acc_class_dalitz = HistogramAcceptance(filename_mc, config, decay)
        p_T = utils_dalitz.get_bin_var()
        for i in range(1, len(p_T)-1):
            x = h1.GetBinContent(i)
            y= h1.GetBinCenter(i)
            #print("eff_dalitz Bin content: ", y, x )
            x2 = h2.GetBinContent(i)
            y2 = h2.GetBinCenter(i)
            #print("eff_pcm Bin content: ", y2, x2)

        max_h1 = h1.GetMaximum()
        max_h2 = h2.GetMaximum()
        max_y = 1.5 * max(max_h1, max_h2)  # 50% mehr als der höchste Balken

        h1.SetMaximum(max_y)


        # Erstelle einen Canvas
        c1 = ROOT.TCanvas("eff_dal_vs_pcm", "eff_dal_vs_pcm", 800, 600)
        c1.Divide(1, 2)

        p1 = c1.cd(1)
        # p1 = c1.cd(1)
        p1.SetPad(0.0, 0.0, 1, 1)
        p1.SetMargin(0.15,0.02,0.22,0.0)
        p1.SetTicks(1,1)
        p1.SetLogy()
        p1.SetLogx()
        # p2 = c1.cd(2);
        #pad.SetPad(0.0, 0.01, 1, 1)
        #pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        #pad.SetTicks(1,1)
        #pad.SetLogy()

        # Einstellungen für das erste Histogramm
        h1.SetLineColor(ROOT.kOrange)
        h1.SetMarkerColor(ROOT.kOrange)
        h1.SetMarkerStyle(20)  # Kreisförmige Punkte
        h1.Draw("E1P")          # Zeichne das erste Histogramm als Punkte mit Fehlerbalken

        # Einstellungen für das zweite Histogramm
        h2.SetLineColor(ROOT.kViolet)
        h2.SetMarkerColor(ROOT.kViolet)
        h2.SetMarkerStyle(34)  # Dreieckige Punkte
        h2.Draw("E1P same")     # Zeichne das zweite Histogramm über das erste

        # Legende hinzufügen
        legend = ROOT.TLegend(0.2, 0.9, 0.4, 0.8)
        #legend = TLegend(0.6, 0.8, 1.0, 0.9)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)
        legend.AddEntry(h1, "Efficiency for #pi^{0} -> e^{+} e^{-} #gamma", "LP")
        legend.AddEntry(h2, "Efficiency for #pi^{0} -> #gamma #gamma", "LP")
        legend.Draw()

        # Add text 
        txt = TPaveText(0.3, 0.975, 0.3, 0.925,"NDC") #0.875,0.85,0.875,0.85,"NDC")
        txt.SetFillColor(kWhite)
        txt.SetFillStyle(0)
        txt.SetBorderSize(0)
        txt.SetTextAlign(33);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.02)
        txt.AddText("this thesis")
        txt.Draw()
        ROOT.SetOwnership(txt,False)

        txt3 = TPaveText(0.325,0.945,0.325,0.895,"NDC")
        txt3.SetFillColor(kWhite)
        txt3.SetFillStyle(0)
        txt3.SetBorderSize(0)
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.02)
        txt3.AddText("pp, #sqrt{s} = 13.6TeV")
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

        # Achsenbeschriftungen und Titel
        h1.SetTitle("Efficiency comparison dalitz vs. pcm")
        h1.SetXTitle("p_{T} (GeV/c)")
        h1.GetYaxis().SetTitleSize(0.04)
        h1.GetYaxis().SetTitleOffset(1)
        h1.SetYTitle("Efficiency")

        p2 = c1.cd(2)
        p2.SetPad(0,0,1,0.3)
        p2.SetMargin(0.15,0.02,0.22,0.0)
        p2.SetTicks(1,1)
        p2.SetLogx()
        p2.SetLogy()

        ratio.SetLineColor(ROOT.kViolet)
        ratio.SetMarkerColor(ROOT.kViolet)
        ratio.SetMarkerStyle(34)  # Kreisförmige Punkte
        ratio.Draw("E1P")   
         
        legend2 = ROOT.TLegend(0.2,0.9,0.3,0.8)
        legend2.SetBorderSize(0)
        legend2.SetFillColor(kWhite)
        legend2.SetFillStyle(0)
        legend2.SetTextSize(0.05)
        legend2.AddEntry(ratio, "ratio dalitz/pcm", "LP")
        legend2.Draw()  

        ratio.SetTitle("Ratio Dalitz/pcm")
        ratio.SetXTitle("p_{T} (GeV/c)")
        ratio.GetXaxis().SetLabelSize(0.08)
        ratio.GetXaxis().SetTitleSize(0.08)
        ratio.GetXaxis().SetTitleOffset(1.25)
        ratio.GetYaxis().SetLabelSize(0.10)
        
        #ratio.GetYaxis().SetRangeUser(-1, 1.5)
        ratio.GetYaxis().SetTitleSize(0.2)
        ratio.GetYaxis().SetTitleOffset(.4)
        ratio.SetYTitle("Ratio")

        c1.Update()

        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())
        # else:
        #     c1.SaveAs("Comparison_eff_dalitz_pcm")
        #c1.SaveAs("comparison_eff_updated")
    
    # def plot_raw_dalitz_vs_pcm_MC(self, eff_dalitz, eff_pcm):
    #     h1 = eff_dalitz.get_raw_yield(hist_name_raw)
    #     h2 = eff_pcm.get_raw_yield(hist_name_raw)
    #     ratio = h1.Clone("Ratio_raw")
    #     ratio.Divide(h2)
    #     acc_class_dalitz = HistogramAcceptance(filename, dir_name21, dir_name22_dalitz, dir_ev, config)
    #     p_T = acc_class_dalitz.get_bin_var()
    #     for i in range(1, len(p_T)-1):
    #         x = h1.GetBinContent(i)
    #         y= h1.GetBinCenter(i)
    #         # print("raw_dalitz Bin content: ", y, x )
    #         x2 = h2.GetBinContent(i)
    #         y2 = h2.GetBinCenter(i)
    #         # print("raw_pcm Bin content: ", y2, x2)

    #     #h1 = acc_class_dalitz.plot_acceptance_vs_pt_pc(hist_name_acc, hist_name_all)
    #     # h2 = acc_class_pcm.plot_acceptance_vs_pt_pc(hist_name_acc, hist_name_all)

    #     max_h1 = h1.GetMaximum()
    #     max_h2 = h2.GetMaximum()
    #     max_y = 1.5 * max(max_h1, max_h2)  # 5% mehr als der höchste Balken

    #     h1.SetMaximum(max_y)

    #     # Erstelle einen Canvas
    #     c1 = ROOT.TCanvas("raw_dal_vs_pcm", "raw_dal_vs_pcm", 800, 600)
    #     #c1.Divide(1, 2)
    #     pad = c1.cd()
    #     pad.SetPad(0.0, 0.0, 1, 1)
    #     pad.SetMargin(0.15,0.02,0.22,0.0)
    #     pad.SetTicks(1,1)
    #     pad.SetLogy()
    #     pad.SetLogx()

    #     # Einstellungen für das erste Histogramm
    #     h1.SetLineColor(ROOT.kBlue)
    #     h1.SetMarkerColor(ROOT.kBlue)
    #     h1.SetMarkerStyle(21)  # Kreisförmige Punkte
    #     h1.Draw("E1P")          # Zeichne das erste Histogramm als Punkte mit Fehlerbalken

    #     # Einstellungen für das zweite Histogramm
    #     h2.SetLineColor(ROOT.kMagenta)
    #     h2.SetMarkerColor(ROOT.kMagenta)
    #     h2.SetMarkerStyle(22)  # Dreieckige Punkte
    #     h2.Draw("E1P same")     # Zeichne das zweite Histogramm über das erste

    #     # Legende hinzufügen
    #     legend = ROOT.TLegend(0.7, 0.9, 0.9, 0.85)
    #     #legend = TLegend(0.6, 0.8, 1.0, 0.9)
    #     legend.SetBorderSize(0)
    #     legend.SetFillColor(kWhite)
    #     legend.SetFillStyle(0)
    #     legend.SetTextSize(0.03)
    #     legend.AddEntry(h1, "raw yield for #pi^{0} -> e^{+} e^{-} #gamma", "LP")
    #     legend.AddEntry(h2, "raw yield for #pi^{0} -> #gamma #gamma", "LP")
    #     legend.Draw()

    #     # Add text 
    #     txt = TPaveText(0.875,0.96,0.875,0.96,"NDC")
    #     txt.SetFillColor(kWhite)
    #     txt.SetFillStyle(0)
    #     txt.SetBorderSize(0)
    #     txt.SetTextAlign(33);#middle,left
    #     txt.SetTextFont(42);#helvetica
    #     txt.SetTextSize(0.02)
    #     txt.AddText("this thesis")
    #     txt.Draw()
    #     ROOT.SetOwnership(txt,False)

    #     txt3 = TPaveText(0.875,0.93,0.875,0.93,"NDC")
    #     txt3.SetFillColor(kWhite)
    #     txt3.SetFillStyle(0)
    #     txt3.SetBorderSize(0)
    #     txt3.SetTextAlign(33);#middle,left
    #     txt3.SetTextFont(42);#helvetica
    #     txt3.SetTextSize(0.02)
    #     txt3.AddText("pp, #sqrt{s} = 13.6TeV")
    #     txt3.Draw()
    #     ROOT.SetOwnership(txt3,False)

    #     # Achsenbeschriftungen und Titel
    #     h1.SetTitle("Acceptance comparison dalitz vs. pcm")
        
    #     h1.SetXTitle("p_{T} (GeV/c)")
    #     h1.SetYTitle("raw yield")
    #     h1.GetXaxis().SetTitleOffset(1.5)
    
    #     #h1.GetYaxis().SetLabelSize(0.10)
        
    #     # # ratio.GetYaxis().SetRangeUser(0, 15)
    #     # ratio.GetYaxis().SetTitleSize(0.10)
    #     # ratio.GetYaxis().SetTitleOffset(.3)
    #     # p2 = c1.cd(2)
    #     # p2.SetPad(0,0,1,0.3)
    #     # p2.SetMargin(0.15,0.02,0.22,0.0)
    #     # p2.SetTicks(1,1)
    #     # p2.SetLogx()
    #     # # p2.SetLogy()

    #     # ratio.SetLineColor(ROOT.kBlue)
    #     # ratio.SetMarkerColor(ROOT.kBlue)
    #     # ratio.SetMarkerStyle(21)  # Kreisförmige Punkte
    #     # ratio.Draw("E1P")   
         
    #     # legend2 = ROOT.TLegend(0.2,0.9,0.3,0.8)
    #     # legend2.SetBorderSize(0)
    #     # legend2.SetFillColor(kWhite)
    #     # legend2.SetFillStyle(0)
    #     # legend2.SetTextSize(0.03)
    #     # legend2.AddEntry(ratio, "ratio dalitz/pcm", "LP")
    #     # legend2.Draw()  

    #     # ratio.SetTitle("Ratio Dalitz/pcm")
    #     # ratio.SetXTitle("p_{T} (GeV/c)")
    #     # ratio.GetXaxis().SetLabelSize(0.10)
    #     # ratio.GetXaxis().SetTitleSize(0.10)
    #     # ratio.GetXaxis().SetTitleOffset(1.0)
    #     # ratio.GetYaxis().SetLabelSize(0.10)
        
    #     # # ratio.GetYaxis().SetRangeUser(0, 15)
    #     # ratio.GetYaxis().SetTitleSize(0.10)
    #     # ratio.GetYaxis().SetTitleOffset(.3)
    #     # ratio.SetYTitle("Ratio")

    #     c1.Update()

    #     c1.Update()
    #     c1.SaveAs("comparison_raw_yield_dalitz_pcm_MC_updated")

    # def plot_raw_dalitz_vs_pcm_data(self, corr_dalitz, corr_pcm):
    #     h1 = corr_dalitz.get_raw_yield_data(filename1_data_dalitz, dir_name11_dalitz, dir_name12_dalitz, dir_name13_data, dir_name14, hist_name_raw)
    #     h2 = corr_pcm.get_raw_yield_data(filename1_data_pcm, dir_name11_pcm, dir_name12_pcm, dir_name13_data, dir_name14, hist_name_raw)
    #     ratio = h1.Clone("Ratio_raw")
    #     ratio.Divide(h2)
    #     acc_class_dalitz = HistogramAcceptance(filename, dir_name21, dir_name22_dalitz, dir_ev, config)
    #     p_T = acc_class_dalitz.get_bin_var()
    #     for i in range(1, len(p_T)-1):
    #         x = h1.GetBinContent(i)
    #         y= h1.GetBinCenter(i)
    #         # print("raw_dalitz Bin content: ", y, x )
    #         x2 = h2.GetBinContent(i)
    #         y2 = h2.GetBinCenter(i)
    #         # print("raw_pcm Bin content: ", y2, x2)

    #     #h1 = acc_class_dalitz.plot_acceptance_vs_pt_pc(hist_name_acc, hist_name_all)
    #     # h2 = acc_class_pcm.plot_acceptance_vs_pt_pc(hist_name_acc, hist_name_all)

    #     max_h1 = h1.GetMaximum()
    #     max_h2 = h2.GetMaximum()
    #     max_y = 1.5 * max(max_h1, max_h2)  # 50% mehr als der höchste Balken

    #     h1.SetMaximum(max_y)

    #     # Erstelle einen Canvas
    #     c1 = ROOT.TCanvas("raw_dal_vs_pcm_data", "raw_dal_vs_pcm_data", 800, 600)
    #     #c1.Divide(1, 2)
    #     pad = c1.cd()
    #     pad.SetPad(0.0, 0.0, 1, 1)
    #     pad.SetMargin(0.15,0.02,0.22,0.0)
    #     pad.SetTicks(1,1)
    #     pad.SetLogy()
    #     pad.SetLogx()

    #     # Einstellungen für das erste Histogramm
    #     h1.SetLineColor(ROOT.kBlue)
    #     h1.SetMarkerColor(ROOT.kBlue)
    #     h1.SetMarkerStyle(21)  # Kreisförmige Punkte
    #     h1.Draw("E1P")          # Zeichne das erste Histogramm als Punkte mit Fehlerbalken

    #     # Einstellungen für das zweite Histogramm
    #     h2.SetLineColor(ROOT.kMagenta)
    #     h2.SetMarkerColor(ROOT.kMagenta)
    #     h2.SetMarkerStyle(22)  # Dreieckige Punkte
    #     h2.Draw("E1P same")     # Zeichne das zweite Histogramm über das erste

    #     # Legende hinzufügen
    #     legend = ROOT.TLegend(0.7, 0.9, 0.9, 0.85)
    #     #legend = TLegend(0.6, 0.8, 1.0, 0.9)
    #     legend.SetBorderSize(0)
    #     legend.SetFillColor(kWhite)
    #     legend.SetFillStyle(0)
    #     legend.SetTextSize(0.03)
    #     legend.AddEntry(h1, "raw yield for #pi^{0} -> e^{+} e^{-} #gamma", "LP")
    #     legend.AddEntry(h2, "raw yield for #pi^{0} -> #gamma #gamma", "LP")
    #     legend.Draw()

    #     # Add text 
    #     txt = TPaveText(0.875,0.96,0.875,0.96,"NDC")
    #     txt.SetFillColor(kWhite)
    #     txt.SetFillStyle(0)
    #     txt.SetBorderSize(0)
    #     txt.SetTextAlign(33);#middle,left
    #     txt.SetTextFont(42);#helvetica
    #     txt.SetTextSize(0.02)
    #     txt.AddText("this thesis")
    #     txt.Draw()
    #     ROOT.SetOwnership(txt,False)

    #     txt3 = TPaveText(0.875,0.93,0.875,0.93,"NDC")
    #     txt3.SetFillColor(kWhite)
    #     txt3.SetFillStyle(0)
    #     txt3.SetBorderSize(0)
    #     txt3.SetTextAlign(33);#middle,left
    #     txt3.SetTextFont(42);#helvetica
    #     txt3.SetTextSize(0.02)
    #     txt3.AddText("pp, #sqrt{s} = 13.6TeV")
    #     txt3.Draw()
    #     ROOT.SetOwnership(txt3,False)

    #     # Achsenbeschriftungen und Titel
    #     h1.SetTitle("Acceptance comparison dalitz vs. pcm")
    #     h1.SetXTitle("p_{T} (GeV/c)")
    #     h1.SetYTitle("raw yield")
    #     h1.GetXaxis().SetTitleOffset(1.5)
    #     # p2 = c1.cd(2)
    #     # p2.SetPad(0,0,1,0.3)
    #     # p2.SetMargin(0.15,0.02,0.22,0.0)
    #     # p2.SetTicks(1,1)
    #     # p2.SetLogx()
    #     # # p2.SetLogy()

    #     # ratio.SetLineColor(ROOT.kBlue)
    #     # ratio.SetMarkerColor(ROOT.kBlue)
    #     # ratio.SetMarkerStyle(21)  # Kreisförmige Punkte
    #     # ratio.Draw("E1P")   
         
    #     # legend2 = ROOT.TLegend(0.2,0.9,0.3,0.8)
    #     # legend2.SetBorderSize(0)
    #     # legend2.SetFillColor(kWhite)
    #     # legend2.SetFillStyle(0)
    #     # legend2.SetTextSize(0.03)
    #     # legend2.AddEntry(ratio, "ratio dalitz/pcm", "LP")
    #     # legend2.Draw()  

    #     # ratio.SetTitle("Ratio Dalitz/pcm")
    #     # ratio.SetXTitle("p_{T} (GeV/c)")
    #     # ratio.GetXaxis().SetLabelSize(0.10)
    #     # ratio.GetXaxis().SetTitleSize(0.10)
    #     # ratio.GetXaxis().SetTitleOffset(1.0)
    #     # ratio.GetYaxis().SetLabelSize(0.10)
        
    #     # # ratio.GetYaxis().SetRangeUser(0, 15)
    #     # ratio.GetYaxis().SetTitleSize(0.10)
    #     # ratio.GetYaxis().SetTitleOffset(.3)
    #     # ratio.SetYTitle("Ratio")

    #     c1.Update()

    #     c1.Update()
    #     c1.SaveAs("comparison_raw_yield_dalitz_pcm_data_updated")

        
    def plot_acc_dalitz_vs_pcm(self, filename_mc, config_array_dalitz, config_array_pcm, name_plot, save):
        
        config_file_dalitz = config_array_dalitz[0]
        config_file_pcm = config_array_pcm[0]
        with open(config_file_dalitz, "r", encoding="utf-8") as config_dalitz_yml:
            config_dalitz = yaml.safe_load(config_dalitz_yml)
        with open(config_file_pcm, "r", encoding="utf-8") as config_pcm_yml:
            config_pcm = yaml.safe_load(config_pcm_yml)
        utils_pcm = Utility(filename_mc, config_pcm, "pcm", "mc")
        utils_dalitz = Utility(filename_mc, config_dalitz, "dalitz", "mc")
        acc_class_dalitz= HistogramAcceptance(filename_mc, config_dalitz, "dalitz")
        acc_class_pcm = HistogramAcceptance(filename_mc, config_pcm, "pcm")
        h1 = acc_class_dalitz.calc_acceptance_vs_pt("acceptance_dalitz", "false", utils_dalitz)
        h2 = acc_class_pcm.calc_acceptance_vs_pt("acceptance_dalitz", "false", utils_pcm)

        # ratio = h1.Clone("Ratio_acc")
        # ratio.Divide(h2)

        #acc_class_dalitz = HistogramAcceptance(filename, dir_name21, dir_name22_dalitz, dir_ev, config)
        p_T = utils_dalitz.get_bin_var()
        for i in range(1, len(p_T)-1):
            x = h1.GetBinContent(i)
            y= h1.GetBinCenter(i)
            print("acc_dalitz Bin content: ", y, x )
            x2 = h2.GetBinContent(i)
            y2 = h2.GetBinCenter(i)
            print("acc_pcm Bin content: ", y2, x2)


        max_h1 = h1.GetMaximum()
        max_h2 = h2.GetMaximum()
        max_y = 1.5 * max(max_h1, max_h2)  # 5% mehr als der höchste Balken

        h1.SetMaximum(max_y)

        # Erstelle einen Canvas
        c1 = ROOT.TCanvas("acc_dal_vs_pcm", "acc_dal_vs_pcm", 800, 600)
        # c1.Divide(1, 2)
        pad = c1.cd()
        
        pad.SetPad(0.0, 0.0, 1, 1)
        pad.SetMargin(0.15,0.02,0.22,0.0)
        pad.SetTicks(1,1)
        pad.SetLogy()
        #pad.SetLogx()

        # Einstellungen für das erste Histogramm
        h1.SetLineColor(ROOT.kOrange)
        h1.SetMarkerColor(ROOT.kOrange)
        h1.SetMarkerStyle(20)  # Kreisförmige Punkte
        h1.Draw("E1P")          # Zeichne das erste Histogramm als Punkte mit Fehlerbalken

        # Einstellungen für das zweite Histogramm
        h2.SetLineColor(ROOT.kViolet)
        h2.SetMarkerColor(ROOT.kViolet)
        h2.SetMarkerStyle(34)  # Dreieckige Punkte
        h2.Draw("E1P same")     # Zeichne das zweite Histogramm über das erste

        # Legende hinzufügen
        legend = ROOT.TLegend(0.5, 0.8, 0.8, 0.7)
        #legend = TLegend(0.6, 0.8, 1.0, 0.9)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.035)
        legend.AddEntry(h1, "Acceptance x BR for #pi^{0} -> e^{+} e^{-} #gamma", "LP")
        legend.AddEntry(h2, "Acceptance x BR for #pi^{0} -> #gamma #gamma", "LP")
        legend.Draw()

        # Add text 
        txt = TPaveText(0.875,0.85,0.875,0.85,"NDC")
        txt.SetFillColor(kWhite)
        txt.SetFillStyle(0)
        txt.SetBorderSize(0)
        txt.SetTextAlign(33);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.02)
        txt.AddText("this thesis")
        txt.Draw()
        ROOT.SetOwnership(txt,False)

        txt3 = TPaveText(0.875,0.82,0.875,0.82,"NDC")
        txt3.SetFillColor(kWhite)
        txt3.SetFillStyle(0)
        txt3.SetBorderSize(0)
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.02)
        txt3.AddText("pp, #sqrt{s} = 13.6TeV")
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

        # Achsenbeschriftungen und Titel
        h1.SetTitle("Acceptance comparison dalitz vs. pcm")
        h1.SetXTitle("p_{T} (GeV/c)")
        h1.SetYTitle("Acceptance x BR")
        h1.GetXaxis().SetRangeUser(0,12)
        # p2 = c1.cd(2)
        # p2.SetPad(0,0,1,0.3)
        # p2.SetMargin(0.15,0.02,0.22,0.0)
        # p2.SetTicks(1,1)
        # #p2.SetLogx()
        # p2.SetLogy()

        # ratio.SetLineColor(ROOT.kBlue)
        # ratio.SetMarkerColor(ROOT.kBlue)
        # ratio.SetMarkerStyle(21)  # Kreisförmige Punkte
        # ratio.Draw("E1P")   
         
        # legend2 = ROOT.TLegend(0.2,0.25,0.6,0.4)
        # legend2.SetBorderSize(0)
        # legend2.SetFillColor(kWhite)
        # legend2.SetFillStyle(0)
        # legend2.SetTextSize(0.05)
        # legend2.AddEntry(ratio, "ratio dalitz/pcm", "LP")
        # legend2.Draw()  

        # ratio.SetTitle("Ratio Dalitz/pcm")
        # ratio.SetXTitle("p_{T} (GeV/c)")
        # ratio.GetXaxis().SetLabelSize(0.10)
        # ratio.GetXaxis().SetTitleSize(0.10)
        # ratio.GetXaxis().SetTitleOffset(1.0)
        # ratio.GetXaxis().SetRangeUser(0,12)
        # ratio.GetYaxis().SetLabelSize(0.10)
        
        # #ratio.GetYaxis().SetRangeUser(8, 14)
        # ratio.GetYaxis().SetTitleSize(0.10)
        # ratio.GetYaxis().SetTitleOffset(.5)
        # ratio.SetYTitle("Ratio")

        c1.Update()

        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())
        # else:
        #     c1.SaveAs("Comparison_acc_dalitz_pcm")
        #c1.SaveAs("comparison_acc_updated")

    # def plot_calc_yield(self, corr_dalitz, corr_pcm):
    #     h1 = corr_dalitz.calculate_corrected_yield_dalitz(save_name1, acc_class_dalitz, eff_dalitz, hist_name_acc, hist_name_all, hist_name_raw)
    #     h2 = corr_pcm.calculate_corrected_yield_pcm(save_name2, acc_class_pcm, eff_pcm, hist_name_acc, hist_name_all, hist_name_raw)

    #     #acc_class_dalitz = HistogramAcceptance(filename, dir_name21, dir_name22_dalitz, dir_ev, config)
    #     p_T = acc_class_dalitz.get_bin_var()
    #     for i in range(1, len(p_T)-1):
    #         x = h1.GetBinContent(i)
    #         y= h1.GetBinCenter(i)
    #         print("correctedyield_dalitz Bin content: ", y, x )
    #         x2 = h2.GetBinContent(i)
    #         y2 = h2.GetBinCenter(i)
    #         print("correctedyield_pcm Bin content: ", y2, x2)

        

    #     max_h1 = h1.GetMaximum()
    #     max_h2 = h2.GetMaximum()
    #     max_y = 1.05 * max(max_h1, max_h2)  # 5% mehr als der höchste Balken

    #     h1.SetMaximum(max_y)

    #     # Erstelle einen Canvas
    #     c1 = ROOT.TCanvas("acc_dal_vs_pcm", "acc_dal_vs_pcm", 800, 600)
    #     pad = c1.cd()
    #     pad.SetPad(0.0, 0.01, 1, 1)
    #     pad.SetMargin(0.15, 0.1, 0.1, 0.1)
    #     pad.SetTicks(1,1)
    #     pad.SetLogy()

    #     # Einstellungen für das erste Histogramm
    #     h1.SetLineColor(ROOT.kBlue)
    #     h1.SetMarkerColor(ROOT.kBlue)
    #     h1.SetMarkerStyle(21)  # Kreisförmige Punkte
    #     h1.Draw("E1P")          # Zeichne das erste Histogramm als Punkte mit Fehlerbalken

    #     # Einstellungen für das zweite Histogramm
    #     h2.SetLineColor(ROOT.kMagenta)
    #     h2.SetMarkerColor(ROOT.kMagenta)
    #     h2.SetMarkerStyle(22)  # Dreieckige Punkte
    #     h2.Draw("E1P same")     # Zeichne das zweite Histogramm über das erste

    #     # Legende hinzufügen
    #     legend = ROOT.TLegend(0.6, 0.8, 0.8, 0.7)
    #     #legend = TLegend(0.6, 0.8, 1.0, 0.9)
    #     legend.SetBorderSize(0)
    #     legend.SetFillColor(kWhite)
    #     legend.SetFillStyle(0)
    #     legend.SetTextSize(0.03)
    #     legend.AddEntry(h1, "Acceptance for #pi^{0} -> e^{+} e^{-} #gamma", "LP")
    #     legend.AddEntry(h2, "Acceptance for #pi^{0} -> e^{+} e^{-}", "LP")
    #     legend.Draw()

    #     # Add text 
    #     txt = TPaveText(0.875,0.85,0.875,0.85,"NDC")
    #     txt.SetFillColor(kWhite)
    #     txt.SetFillStyle(0)
    #     txt.SetBorderSize(0)
    #     txt.SetTextAlign(33);#middle,left
    #     txt.SetTextFont(42);#helvetica
    #     txt.SetTextSize(0.02)
    #     txt.AddText("this thesis")
    #     txt.Draw()
    #     ROOT.SetOwnership(txt,False)

    #     txt3 = TPaveText(0.875,0.82,0.875,0.82,"NDC")
    #     txt3.SetFillColor(kWhite)
    #     txt3.SetFillStyle(0)
    #     txt3.SetBorderSize(0)
    #     txt3.SetTextAlign(33);#middle,left
    #     txt3.SetTextFont(42);#helvetica
    #     txt3.SetTextSize(0.02)
    #     txt3.AddText("pp, #sqrt{s} = 13.6TeV")
    #     txt3.Draw()
    #     ROOT.SetOwnership(txt3,False)

    #     # Achsenbeschriftungen und Titel
    #     h1.SetTitle("Acceptance comparison dalitz vs. pcm")
    #     h1.SetXTitle("p_{T} (GeV/c)")
    #     h1.SetYTitle("Acceptance [%]")

    #     c1.Update()
    #     c1.SaveAs("comparison_corr_test")

    def compare_corr_dalitz_pcm(self, corr_dalitz, corr_pcm, name_plot, save): #raw_data_dalitz, eff_hist_dalitz, acc_hist_dalitz, raw_data_pcm, eff_hist_pcm, acc_hist_pcm): #,corr_dalitz, corr_pcm):
        # corr_class_dalitz = CorrectedYield(filename_mc, "dalitz")
        # corr_class_pcm = CorrectedYield(filename_mc, "pcm")
        # acc_class_dalitz = HistogramAcceptance(filename_mc, config, "dalitz")
        # acc_class_pcm = HistogramAcceptance(filename_mc, config, "pcm")
        
        #h1 = corr_dalitz.calculate_corrected_yield_dalitz(filename1_data_dalitz, hist_name_acc, hist_name_all, hist_name_raw, acc_class_dalitz, eff_dalitz, dir_name11_dalitz, dir_name12_dalitz, dir_name13_data, dir_name14)
        #h2 = corr_pcm.calculate_corrected_yield_pcm(filename1_data_pcm, hist_name_acc, hist_name_all, hist_name_raw, acc_class_pcm, eff_pcm, dir_name11_pcm, dir_name12_pcm, dir_name13_data, dir_name14)
        #h1 = corr_class_dalitz.calc_corr_yield_vs_pt("corrected_yield_dalitz", "false", raw_data_dalitz, eff_hist_dalitz, acc_hist_dalitz, acc_class_dalitz)
        #h2 = corr_class_pcm.calc_corr_yield_vs_pt("corrected_yield_pcm", "false", raw_data_pcm, eff_hist_pcm, acc_hist_pcm, acc_class_pcm)

        h1 = corr_dalitz 
        h2 = corr_pcm
        print("h1 typ: ", type(h1))
        print("h2 typ: ", type(h2))

        ratio = h1.Clone("Ratio_corr_yield")
        ratio.Divide(h2)

        max_h1 = h1.GetMaximum()
        min_h1 = h1.GetMinimum()
        max_h2 = h2.GetMaximum()
        min_h2 = h2.GetMinimum()
        max_y = 1.5 * max(max_h1, max_h2)  # 5% mehr als der höchste Balken
        min_y =  min(min_h1, min_h2)

        h1.SetMaximum(max_y)
        print("min_y :", min_y)
        #h1.SetMinimum(min_y)
        # y_min = min(min(y1), min(y2))
        # y_max = max(max(y1), max(y2))

        # Setzen der y-Achse mit einem Puffer
        #h1.GetYaxis().SetRangeUser(min_y - 0.1*min_y, max_y + 0.1*max_y)


        # Erstelle einen Canvas
        c1 = ROOT.TCanvas("py_vs_cor", "py_vs_cor", 800, 600)
        c1.Divide(1, 2)

        p1 = c1.cd(1)
        p1.SetPad(0.0, 0.0, 1, 1)
        p1.SetMargin(0.15,0.02,0.22,0.0)
        p1.SetTicks(1,1)
        p1.SetLogy()
        p1.SetLogx()
        
        h1.SetLineColor(ROOT.kBlue)
        h1.SetMarkerColor(ROOT.kBlue)
        h1.SetMarkerStyle(21)
        #h1.GetYaxis().SetRangeUser(1e-5, 0)  # Kreisförmige Punkte
        h1.Draw("E1P")          # Zeichne das erste Histogramm als Punkte mit Fehlerbalken

        # Einstellungen für das zweite Histogramm
        h2.SetLineColor(ROOT.kMagenta)
        h2.SetMarkerColor(ROOT.kMagenta)
        h2.SetMarkerStyle(22)  # Dreieckige Punkte
        #h2.GetYaxis().SetRangeUser(1e-5, 0)
        h2.Draw("E1P same")     # Zeichne das zweite Histogramm über das erste

        # Legende hinzufügen
        legend1 = ROOT.TLegend(0.6, 0.9, 0.95, 0.75)
        legend1.SetBorderSize(0)
        legend1.SetFillColor(kWhite)
        legend1.SetFillStyle(0)
        legend1.SetTextSize(0.03)
        legend1.AddEntry(h1, "corrected yield for #pi^{0} -> e^{+} e^{-} #gamma", "LP")
        legend1.AddEntry(h2, "corrected yield for #pi^{0} -> #gamma #gamma", "LP")
        #legend1.AddEntry(h2, "corrected yield for #pi^{0} -> e^{+} e^{-}", "LP")
        legend1.Draw()

        txt = TPaveText(0.9,0.975,0.95,0.95,"NDC")
        txt.SetFillColor(kWhite)
        txt.SetFillStyle(0)
        txt.SetBorderSize(0)
        txt.SetTextAlign(33);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.02)
        txt.AddText("this thesis")
        txt.Draw()
        ROOT.SetOwnership(txt,False)

        txt3 = TPaveText(0.9,0.945,0.95,0.92,"NDC"); # txt3 = TPaveText(0.4,0.75,0.4,0.90,"NDC");
        txt3.SetFillColor(kWhite)
        txt3.SetFillStyle(0)
        txt3.SetBorderSize(0)
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.02)
        txt3.AddText("pp, #sqrt{s} = 13.6TeV")
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

        # Achsenbeschriftungen und Titel
        h1.SetTitle("Corrected yield")
        h1.SetXTitle("p_{T} (GeV/c)")
        h1.SetYTitle("#frac{1}{2#pi} #frac{d^{2}N}{dp_{T}dy}")

        #Go to second canvas to plot ratio in
        
        p2 = c1.cd(2)
        p2.SetPad(0,0,1,0.3)
        p2.SetMargin(0.15,0.02,0.22,0.0)
        p2.SetTicks(1,1)
        p2.SetLogx()
        # p2.SetLogy()

        ratio.SetLineColor(ROOT.kBlue)
        ratio.SetMarkerColor(ROOT.kBlue)
        ratio.SetMarkerStyle(21)  # Kreisförmige Punkte
        ratio.Draw("E1P")   
         
        legend2 = ROOT.TLegend(0.2,0.9,0.3,0.8)
        legend2.SetBorderSize(0)
        legend2.SetFillColor(kWhite)
        legend2.SetFillStyle(0)
        legend1.SetTextSize(0.03)
        legend2.AddEntry(ratio, "ratio dalitz/pcm", "LP")
        legend2.Draw()  

        ratio.SetTitle("Ratio Dalitz/pcm")
        ratio.SetXTitle("p_{T} (GeV/c)")
        ratio.GetXaxis().SetLabelSize(0.10)
        ratio.GetXaxis().SetTitleSize(0.10)
        ratio.GetXaxis().SetTitleOffset(1.6)
        ratio.GetYaxis().SetLabelSize(0.10)
        
        ratio.GetYaxis().SetRangeUser(-0.5, 2.5)
        ratio.GetYaxis().SetTitleSize(0.10)
        ratio.GetYaxis().SetTitleOffset(.3)
        ratio.SetYTitle("Ratio")

        # Fit Ratio
        fit_result = ratio.Fit("pol0", "S")
        fit = ratio.GetFunction("pol0")
        fit.SetLineColor(ROOT.kCyan)
        fit.SetLineWidth(1)  # Dünne Linie

        # Ausgabe der Fitparameter
        fit_constant = fit_result.Parameter(0)
        fit_error = fit_result.ParError(0)
        fit_info = f"Fit: y = {fit_constant:.3f} #pm {fit_error:.3f}"

        legend2 = ROOT.TLegend(0.15,0.94,0.5,0.75)
        legend2.SetBorderSize(0)
        legend2.SetFillColor(ROOT.kWhite)
        legend2.SetFillStyle(0)
        legend2.SetTextSize(0.05)
        legend2.AddEntry(ratio, "Ratio Dalitz/PCM", "LP")
        # Erstelle einen leeren TGraph
        placeholder_graph = ROOT.TGraph()
        placeholder_graph.SetTitle(fit_info)  # Setze den Titel des Graphen auf die Fit-Informationen
        legend2.AddEntry(fit, fit_info, "l")  # Füge den Graphen mit den Fit-Informationen hinzu
        #legend2.AddEntry(None, fit_info, "")  # Fügt Fit-Informationen hinzu
        legend2.Draw()  

        ratio.SetTitle("Ratio Dalitz/PCM")
        ratio.SetXTitle("p_{T} (GeV/c)")
        ratio.GetXaxis().SetLabelSize(0.09)
        ratio.GetXaxis().SetTitleSize(0.10)
        ratio.GetXaxis().SetTitleOffset(1.0)
        ratio.GetYaxis().SetLabelSize(0.10)
        ratio.GetYaxis().SetRangeUser(0.0, 1.5)
        ratio.GetYaxis().SetTitleSize(0.10)
        ratio.GetYaxis().SetTitleOffset(.3)
        ratio.SetYTitle("Ratio")

        # Gestrichelte Linie bei y=1 hinzufügen
        line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1, ratio.GetXaxis().GetXmax(), 1)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(2)  # Gestrichelte Linie
        line.Draw()

        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())
        #c1.SaveAs("compare_corr_dalitz_vs_pcm_updated")

    def plot_eff_times_acc(self, filename_mc, config, raw_mc_dalitz, raw_mc_pcm, name_plot, save, utils):
        acc_class_dalitz = HistogramAcceptance(filename_mc, config, "dalitz")
        acc_class_pcm = HistogramAcceptance(filename_mc, config, "pcm")
        eff_class_dalitz = Efficiency("dalitz")
        eff_class_pcm = Efficiency("pcm")
        h1_eff = eff_class_dalitz.calc_efficiency_vs_pt("Efficiency_dalitz", "false", raw_mc_dalitz, utils)
        h1_acc = acc_class_dalitz.calc_acceptance_vs_pt("Acceptance_dalitz", "false", utils)
        h2_eff = eff_class_pcm.calc_efficiency_vs_pt("Efficiency_pcm", "false", raw_mc_pcm, utils)
        h2_acc = acc_class_pcm.calc_acceptance_vs_pt("Acceptance_pcm", "false", utils)
        # acc_class_dalitz = HistogramAcceptance(filename, dir_name21, dir_name22_dalitz, dir_ev, config)
        # p_T = acc_class_dalitz.get_bin_var()
        # for i in range(1, len(p_T)-1):
        #     x = h1.GetBinContent(i)
        #     y= h1.GetBinCenter(i)
        #     # print("eff_dalitz Bin content: ", y, x )
        #     x2 = h2.GetBinContent(i)
        #     y2 = h2.GetBinCenter(i)
            # print("eff_pcm Bin content: ", y2, x2)
        # h1_eff = eff_dalitz
        # h1_acc = acc_dalitz
        # h2_eff = eff_pcm
        # h2_acc = acc_pcm
        product_dalitz = h1_eff.Clone("Dalitz_product")
        product_pcm = h2_eff.Clone("PCM_product")
        product_dalitz.Multiply(h1_acc)
        product_pcm.Multiply(h2_acc)

        max_h1 = product_dalitz.GetMaximum()
        max_h2 = product_pcm.GetMaximum()
        max_y = 1.2 * max(max_h1, max_h2)  # 5% mehr als der höchste Balken

        product_dalitz.SetMaximum(max_y)

        # Erstelle einen Canvas
        c1 = ROOT.TCanvas("product_dal_vs_pcm", "product_dal_vs_pcm", 800, 600)
        pad = c1.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()

        # Einstellungen für das erste Histogramm
        product_dalitz.SetLineColor(ROOT.kOrange)
        product_dalitz.SetMarkerColor(ROOT.kOrange)
        product_dalitz.SetMarkerStyle(20)  # Kreisförmige Punkte
        product_dalitz.Draw("E1P")          # Zeichne das erste Histogramm als Punkte mit Fehlerbalken

        # Einstellungen für das zweite Histogramm
        product_pcm.SetLineColor(ROOT.kViolet)
        product_pcm.SetMarkerColor(ROOT.kViolet)
        product_pcm.SetMarkerStyle(34)  # Dreieckige Punkte
        product_pcm.Draw("E1P same")     # Zeichne das zweite Histogramm über das erste

        # Legende hinzufügen
        legend = ROOT.TLegend(0.2, 0.25, 0.4, 0.15)
        #legend = TLegend(0.6, 0.8, 1.0, 0.9)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)
        legend.AddEntry(product_dalitz, "eff x acc x BR for #pi^{0} -> e^{+} e^{-} #gamma", "LP")
        legend.AddEntry(product_pcm, "eff x acc x BRfor #pi^{0} -> #gamma #gamma", "LP")
        legend.Draw()

        # Add text 
        txt = TPaveText(0.875,0.85,0.875,0.85,"NDC")
        txt.SetFillColor(kWhite)
        txt.SetFillStyle(0)
        txt.SetBorderSize(0)
        txt.SetTextAlign(33);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.02)
        txt.AddText("this thesis")
        txt.Draw()
        ROOT.SetOwnership(txt,False)

        txt3 = TPaveText(0.875,0.82,0.875,0.82,"NDC")
        txt3.SetFillColor(kWhite)
        txt3.SetFillStyle(0)
        txt3.SetBorderSize(0)
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.02)
        txt3.AddText("pp, #sqrt{s} = 13.6TeV")
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

        # Achsenbeschriftungen und Titel
        product_dalitz.SetTitle("Efficiency times Acceptance")
        product_dalitz.SetXTitle("p_{T} (GeV/c)")
        product_dalitz.SetYTitle("Efficiency * Acceptance")

        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())
        #c1.SaveAs("eff_times_acc_updated")

    def get_13tev_file(self, filename_comp,hist_comp):
        self.filename_comp = TFile.Open(filename_comp, "READ")
        self.dir_pi0 = self.filename_comp.Get("Pi013TeV") #.Get("AcceptancePi0")
        self.hist_comp = self.dir_pi0.Get(hist_comp) # "AcceptancePi0"
        if not self.hist_comp:
            print ("nicht gefunden")
        if self.hist_comp:
            print("histogram gefunden")
            print("typ von hist: ", type(self.hist_comp))
        return self.hist_comp

        #filename_comp = "/Users/juliaschlagel/Analysis240411/Analysis/Comparison/data_PCMDalitzResultsFullCorrection_PP.root"
    
    def compare_13tev(self, filename_comp, comp_partner, comp_partner_hist, decay, name_plot, save):
        if comp_partner == "acc":
            # h1 = acc_class.calc_acceptance_vs_pt_pc("Acceptance", "false")
            hist_comp = "AcceptancePi0"
            
        if comp_partner == "eff":
            #h1 = eff_class.calc_efficiency_vs_pt("Efficiency_dalitz", "false", raw_mc, acc_class)
            hist_comp = "EfficiencyPi0"
        
        if comp_partner == "corr":
            # h1 = corr_class.calculate_corrected_yield("Corrected_yield_dalitz", filename_mc, config, filename_inv_mass_mc_dalitz, filename_inv_mass_data_dalitz, "false")
            hist_comp = "CorrectedYieldPi0"
        h1 = comp_partner_hist
        if h1: 
            print("Histogram erfolgreich geladen", type(h1))
        
        
        graph = self.get_13tev_file(filename_comp, hist_comp)
        if graph:
            print("graph erfolgreich geladen", type(graph))
            print("histogramm in Graphen: ", hist_comp)
        if not graph:
            print("graph nicht gefunden")

        max_h1 = h1.GetMaximum()
        min_h1 = h1.GetMinimum()
        def get_graph_y_range(graph):
            n_points = graph.GetN()
            y_values = [graph.GetY()[i] for i in range(n_points)]
            y_min = min(y_values)
            y_max = max(y_values)
            return y_min, y_max
        if comp_partner == "corr":
            min_h2 = graph.GetMinimum()
            max_h2 = graph.GetMaximum()
        else: 
            min_h2, max_h2 = get_graph_y_range(graph)

        max_y = 1.5 * max(max_h1, max_h2)  # 50% mehr als der höchste Balken
        #min_y = min(min_h1, min_h2)
        # max_h2 = graph.GetMaximum()
        # min_h2 = graph.GetMinimum()
        # max_y = 1.5 * max(max_h1, max_h2)  # 50% mehr als der höchste Balken
        # min_y =  min(min_h1, min_h2)

        h1.SetMaximum(max_y)
        #print("min_y :", min_y)

        # Erstelle einen Canvas
        if comp_partner == "acc":
            c1 = ROOT.TCanvas("comp_run2_acc", "comp_run2_acc", 800, 600)
        if comp_partner == "eff":
            c1 = ROOT.TCanvas("comp_run2_eff", "comp_run2_eff", 800, 600)
        if comp_partner == "corr":
            c1 = ROOT.TCanvas("comp_run2_corr", "comp_run2_corr", 800, 600)
        #c1.Divide(1, 2)

        p1 = c1.cd()
        p1.SetPad(0.0, 0.0, 1, 1)
        p1.SetMargin(0.15,0.02,0.22,0.0)
        p1.SetTicks(1,1)
        p1.SetLogy()
        #p1.SetLogx()
        
        h1.SetLineColor(ROOT.kBlue)
        h1.SetMarkerColor(ROOT.kBlue)
        h1.SetMarkerStyle(21)
        #h1.GetXaxis().SetRangeUser(0.1, 15)
        #h1.GetYaxis().SetRangeUser(1e-5, 0)  # Kreisförmige Punkte
        h1.Draw("E1P")          # Zeichne das erste Histogramm als Punkte mit Fehlerbalken

        # Einstellungen für das zweite Histogramm
        # graph.SetLineColor(ROOT.kMagenta)
        # graph.SetMarkerColor(ROOT.kMagenta)
        # graph.SetMarkerStyle(22)
        #graph.GetXaxis().SetLimits(0.1, 15)
        graph.SetLineColor(ROOT.kMagenta)
        graph.SetMarkerColor(ROOT.kMagenta)
        graph.SetMarkerStyle(22)  # Dreieckige Punkte
        #h2.GetYaxis().SetRangeUser(1e-5, 0)
        if comp_partner == "corr":
            graph.Draw("E1P Same")
        else:
            graph.Draw("P SAME")
        #h2.Draw("E1P same")     # Zeichne das zweite Histogramm über das erste

        # Legende hinzufügen
        if comp_partner=="acc":
            legend1 = ROOT.TLegend(0.4, 0.4, 0.95, 0.3)
        if comp_partner == "eff":
            legend1 = ROOT.TLegend(0.4, 0.4, 0.95, 0.3)
        if comp_partner == "corr":
            legend1 = ROOT.TLegend(0.4, 0.9, 0.95, 0.75)
        legend1.SetBorderSize(0)
        legend1.SetFillColor(kWhite)
        legend1.SetFillStyle(0)
        legend1.SetTextSize(0.03)
        if comp_partner == "acc":
            if decay == "dalitz":
                legend1.AddEntry(h1, "acceptance for #pi^{0} -> e^{+} e^{-} #gamma at #sqrt{s} = 13.6TeV", "LP")
                legend1.AddEntry(graph, "acceptance for #pi^{0} -> e^{+} e^{-} #gamma at #sqrt{s} = 13.0TeV", "LP")
            if decay == "pcm":
                legend1.AddEntry(h1, "acceptance for #pi^{0} -> #gamma #gamma at #sqrt{s} = 13.6TeV", "LP")
                legend1.AddEntry(graph, "acceptance for #pi^{0} -> #gamma #gamma at #sqrt{s} = 13.0TeV", "LP")
        if comp_partner == "eff":
            if decay == "dalitz":
                legend1.AddEntry(h1, "efficiency for #pi^{0} -> e^{+} e^{-} #gamma at #sqrt{s} = 13.6TeV", "LP")
                legend1.AddEntry(graph, "efficiency for #pi^{0} -> e^{+} e^{-} #gamma at #sqrt{s} = 13.0TeV", "LP")
            if decay == "pcm":
                legend1.AddEntry(h1, "efficiency for #pi^{0} -> #gamma #gamma at #sqrt{s} = 13.6TeV", "LP")
                legend1.AddEntry(graph, "efficiency for #pi^{0} -> #gamma #gamma at #sqrt{s} = 13.0TeV", "LP")
        if comp_partner == "corr":
            if decay == "dalitz":
                legend1.AddEntry(h1, "corrected yield for #pi^{0} -> e^{+} e^{-} #gamma at #sqrt{s} = 13.6TeV", "LP")
                legend1.AddEntry(graph, "corrected yield for #pi^{0} -> e^{+} e^{-} #gamma at #sqrt{s} = 13.0TeV", "LP")
            if decay == "pcm":
                legend1.AddEntry(h1, "corrected yield for #pi^{0} -> #gamma #gamma at #sqrt{s} = 13.6TeV", "LP")
                legend1.AddEntry(graph, "corrected yield for #pi^{0} -> #gamma #gamma at #sqrt{s} = 13.0TeV", "LP")

        #legend1.AddEntry(h2, "corrected yield for #pi^{0} -> e^{+} e^{-}", "LP")
        legend1.Draw()

        txt = TPaveText(0.9,0.975,0.95,0.95,"NDC")
        txt.SetFillColor(kWhite)
        txt.SetFillStyle(0)
        txt.SetBorderSize(0)
        txt.SetTextAlign(33);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.02)
        txt.AddText("this thesis")
        txt.Draw()
        ROOT.SetOwnership(txt,False)

        txt3 = TPaveText(0.9,0.945,0.95,0.92,"NDC"); # txt3 = TPaveText(0.4,0.75,0.4,0.90,"NDC");
        txt3.SetFillColor(kWhite)
        txt3.SetFillStyle(0)
        txt3.SetBorderSize(0)
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.02)
        txt3.AddText("pp, #sqrt{s} = 13.6TeV")
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

        # Achsenbeschriftungen und Titel
        if comp_partner == "acc":
            h1.SetTitle("Acceptance")
            h1.SetXTitle("p_{T} (GeV/c)")
            if decay == "dalitz":
                h1.SetYTitle("Acceptance[%]")
            if decay == "pcm":
                h1.SetYTitle("Acceptance")
        if comp_partner == "eff":
            h1.SetTitle("Efficiency")
            h1.SetXTitle("p_{T} (GeV/c)")
            h1.SetYTitle("Efficiency")
        if comp_partner == "corr":
            h1.SetTitle("Corrected yield")
            h1.SetXTitle("p_{T} (GeV/c)")
            h1.SetYTitle("#frac{1}{2#pip_{T}} #frac{d^{2}N}{dp_{T}dy}")


        c1.Update()
        # if comp_partner == "acc":
        #     c1.SaveAs("compare_13Tev_acceptance")
        # if comp_partner == "eff":
        #     c1.SaveAs("compare_13Tev_efficiency")
        # if comp_partner == "corr":
        #     c1.SaveAs("compare_13Tev_corrected_yield")  
        if save == "true":
            c1.SaveAs(name_plot.Data())

    def compare_cuts(self, fHistolist, cutnames, decay, filename_mc, config, name_plot, comp_partner, utils):
        TGaxis.SetMaxDigits(3);
        acc_class= HistogramAcceptance(filename_mc, config, decay)
        pt = utils.get_bin_var()
        npt = len(pt);

        titlePt = fHistolist[0].GetTitle(); # ACCESS TITLE
        #print("#########fHistoParameter: ", type(fHistoParameter), "#########")
        #print("fHistoParameter",fHistoParameter)
        #print("Länge: ", len(fHistoParameter))
        # y_max = []
        # y_min = []
        # for ic in range(len(fHistolist)):
        #     max_y = max(fHistolist[ic])
        #     y_max.append(max_y)
        #     min_y = min(fHistolist[ic])
        #     y_min.append(min_y)
        # yMax_ = max(y_max)
        # yMin_ = min(y_min)
        if comp_partner == "eff":
            yMin_, yMax_ = 1e-5, 1.4*1e-2
        if comp_partner == "corr":
            yMin_, yMax_ = 0 , 100
        color = [kRed+1, kBlue+1, kGreen+2, kMagenta+2, kCyan+1]

        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);
        p1.SetLogy();


        frame1 = p1.DrawFrame(0., yMin_, 12., yMax_);
        if comp_partner == "eff":
            FrameSettings(frame1, "", "Efficiency",0., 12.)
        if comp_partner == "corr":
            FrameSettings(frame1, "", "#frac{1}{2#pi p_{T}}#frac{d^2 N}{dp_{T} dy}",0., 12.)
        frame1.GetXaxis().SetTitleSize(0.045);
        frame1.GetYaxis().SetTitleSize(0.045);
        frame1.GetXaxis().SetTitleOffset(1.0);
        frame1.GetYaxis().SetTitleOffset(1.52);
        frame1.GetXaxis().SetLabelSize(0.045);
        frame1.GetYaxis().SetLabelSize(0.045);
        frame1.GetYaxis().SetMaxDigits(3);
        frame1.GetXaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetLabelOffset(0.01);

        for icut in range(len(fHistolist)):
            #PlotRawYieldInvMass(meson, filename, dirname)
            #def DrawHisto(self, pad1, histo1, Title, XTitle, YTitle,  markerColor, drawsettings, markerSize =1.):
            fHistolist[icut].GetYaxis().SetLabelSize(0.02);
            fHistolist[icut].GetYaxis().SetTitleSize(0.025);
            fHistolist[icut].GetYaxis().SetDecimals();
            fHistolist[icut].GetXaxis().SetTitleSize(0.025);
            fHistolist[icut].GetXaxis().SetLabelSize(0.02);
            fHistolist[icut].SetMarkerStyle(kFullCircle);
            fHistolist[icut].GetYaxis().SetLabelSize(0.05);
            fHistolist[icut].GetXaxis().SetLabelSize(0.05);
            fHistolist[icut].GetXaxis().SetNdivisions(507, True);

            fHistolist[icut].SetMarkerStyle(kFullCircle);
            fHistolist[icut].SetMarkerColor(color[icut]);
            #fHistolist[icut].SetMarkerSize(markerSize);
            fHistolist[icut].SetLineColor(color[icut]);
            fHistolist[icut].SetLineWidth(1);
            fHistolist[icut].SetFillColor(color[icut]);
            fHistolist[icut].SetFillStyle(0);
            fHistolist[icut].DrawCopy("E1,same");


        p2 = c1.cd(2);
        p2.SetPad(0,0,1,0.3);
        p2.SetMargin(0.15,0.02,0.22,0.0);
        p2.SetTicks(1,1);
        p2.SetLogy();
        frame2 = p2.DrawFrame(0., 5*1e-3, 12., 1e1);
        FrameSettings(frame2, "p_{T} (GeV/c)", "#frac{cut}{analysis}",0., 12.)
        frame2.GetXaxis().SetTitleSize(0.10);
        frame2.GetYaxis().SetTitleSize(0.10);
        frame2.GetXaxis().SetTitleOffset(1.0);
        frame2.GetYaxis().SetTitleOffset(0.7);
        frame2.GetXaxis().SetLabelSize(0.10);
        frame2.GetYaxis().SetLabelSize(0.10);
        frame2.GetYaxis().CenterTitle(True);
        frame2.GetXaxis().SetLabelOffset(0.01);
        frame2.GetYaxis().SetLabelOffset(0.01);
        ROOT.SetOwnership(frame2,False);

        line1 = TLine(0,1,12.,1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(1);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

        for iratio in range(len(fHistolist)-1):
            h1ratio1 = fHistolist[iratio+1].Clone("h1ratio");
            h1ratio1.Reset();
            h1ratio1.Sumw2();
            h1ratio1.Divide(fHistolist[iratio+1], fHistolist[0], 1., 1., "G");
            h1ratio1.DrawCopy("E1,same");
            del h1ratio1

        c1.cd()

        txt = TPaveText(0.,0.9,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);
        if comp_partner == "eff":
            txt.AddText("Efficiency as a function of #it{p}_{T}");
        else:
            txt.AddText("Corrected yield as a function of #it{p}_{T}")
        txt.Draw();
        ROOT.SetOwnership(txt,False);       

        txt = TPaveText(0.65,0.77,0.95,0.87,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(12);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.02);
        txt.AddText("this thesis");
        txt.AddText(decay)
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        leg = TLegend(0.65,0.62,0.95,0.77);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.02);
        leg.SetTextFont(42);#helvetica
        #markersize       = fHistoParameter[0].GetMarkerSize();
        for icut in range(len(fHistolist)):
            #fHistoParameter[icut].SetMarkerSize(markersize);
            leg.AddEntry(fHistolist[icut],"{}".format(cutnames[icut]),"ep");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        c1.SaveAs(name_plot.Data());
        c1.Close();     



# Using the class
# if __name__ == "__main__":
#     filename_comp = "/Users/juliaschlagel/Analysis240411/Comparison/data_PCMDalitzResultsFullCorrection_PP.root"

#     filename_mc = "/Users/juliaschlagel/Analysis240411/Analysis/HLtrains/MC/AnalysisResults_HL197593_LHC24b1.root"
#     dir_name21 = "pi0eta-to-gammagamma-mc/Generated"
#     dir_name22_dalitz = "PCMDalitzEE"              
#     dir_name22_pcm ="PCMPCM"    
#     dir_name_py_1 = "associate-mc-info/Generated"
    
#     #hist = self.filename.Get(dir_name21).FindObject(dir_name22).FindObject(hist_name_all)

#     #filename1_MC_dalitz = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC24b1_invariant_mass_dalitz/this_analysis_LHC24b1_pi0_mc_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
#     #filename1_MC_pcm = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC24b1_invariant_mass_pcm/this_analysis_LHC24b1_pi0_mc_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
#     #filename1_data_dalitz = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC22o_pass6_small_invariant_mass_dalitz/this_analysis_LHC22o_pass6_small_pi0_data_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
#     #filename1_data_pcm = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC22o_pass6_small_invariant_mass_pcm/this_analysis_LHC22o_pass6_small_pi0_data_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"

#     filename_inv_mass_data_dalitz = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC22o_pass6_small_invariant_mass_dalitz_data_MC/this_analysis_LHC22o_pass6_small_pi0_data_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
#     filename_inv_mass_mc_dalitz = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC24b1_invariant_mass_dalitz_data_MC/this_analysis_LHC24b1_pi0_mc_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
#     # filename1_data_dalitz = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC22o_pass6_small_invariant_mass_dalitz_data_MC/this_analysis_LHC22o_pass6_small_pi0_data_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
#     # filename1_MC_dalitz = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC24b1_invariant_mass_dalitz_data_MC/this_analysis_LHC24b1_pi0_mc_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"

#     dir_name11_dalitz = "PCMDalitzEE"
#     dir_name11_pcm = "PCMPCM"
#     dir_name12_dalitz = "qc_mee0_60_minpt100_maxeta09_tpconly_lowB"
#     dir_name12_pcm = "qc_qc"
#     dir_name13_data = "gausexplinear"
#     dir_name13_mc = "gausexp"
#     dir_name14 = "fit_0.04_0.20_GeVc2"
    
#     hist_name_raw = "h1yield_param"
#     hist_name_acc = "hPt_Pi0_Acc"
#     hist_name_all = "hPt_Pi0"
#     dir_ev = "pi0eta-to-gammagamma-mc/Event"
#     config = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC22o_small.yml"
#     acc_class = HistogramAcceptance(filename_mc, config, "dalitz")
#     eff_class = Efficiency("dalitz", acc_class)
#     corr_class = CorrectedYield(filename_mc, filename_inv_mass_data_dalitz, "dalitz")
#     #eff_pcm = Efficiency(acc_class_pcm, "data", "pi0")



    #corr_dalitz = CorrectedYield(filename, filename1_MC_dalitz, filename1_data_dalitz) #, acc_class_dalitz, eff_dalitz)
    #corr_pcm = CorrectedYield(filename, filename1_MC_pcm, filename1_data_dalitz) # acc_class_pcm, eff_pcm)
 #, acc_class_pcm, eff_pcm)
    #corr_pcm = CorrectedYieldDalitz()

    #comp = Compare() #filename, filename1_dalitz, acc_class_dalitz, eff_dalitz, corr_dalitz)
    # comp.plot_raw_dalitz_vs_pcm(eff_dalitz, eff_pcm)

    
    #comp.compare_13tev(filename_comp, "AcceptancePi0", "acc")
    #comp.compare_13tev(filename_comp, "EfficiencyPi0", "eff")
    #comp.compare_13tev(filename_comp, "CorrectedYieldPi0", "corrected_yield")
    
    # comp.plot_calc_yield(corr_dalitz, corr_pcm)
    # comp.compare_corr_dalitz_pcm(corr_dalitz, corr_pcm)
    # comp.plot_eff_times_acc(eff_dalitz, acc_class_dalitz, eff_pcm, acc_class_pcm)
    # comp.plot_acc_dalitz_vs_pcm(acc_class_dalitz, acc_class_pcm)
    # comp.plot_raw_dalitz_vs_pcm_data(corr_dalitz, corr_pcm)
    # comp.plot_raw_dalitz_vs_pcm_MC(eff_dalitz, eff_pcm)
     # comp.plot_raw_dalitz_vs_pcm(eff_dalitz, eff_pcm)
    # comp.plot_acc_dalitz_vs_pcm(acc_class_dalitz, acc_class_pcm)
    # comp.plot_calc_yield(corr_dalitz, corr_pcm)
    #comp.calculate_corrected_yield_dalitz(hist_name_acc, hist_name_all, hist_name_raw)
    # comp.calculate_corrected_yield_pcm(hist_name_acc, hist_name_all, hist_name_raw)
    # comp.calculate_corrected_yield_dalitz(hist_name_acc, hist_name_all, hist_name_raw)
    #comp.calculate_corrected_yield_pcm(save_name2, acc_class_pcm, eff_pcm, hist_name_acc, hist_name_all, hist_name_raw)



