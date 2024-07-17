# This code was written by Julia Schlägel (July 2024)

import ROOT
import ctypes
from ctypes import *
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullCross, kFullDiamond, kOpenSquare, kOpenTriangleUp, kOpenCircle, kFullCrossX
import numpy as np
import yaml
import array
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TGaxis, TPad, TLatex, TBox, TString, TLine
from HistoFormatting import FrameSettings
from PlotRawYield import PlotRawYieldInvMass
from acceptance import HistogramAcceptance
from efficiency import Efficiency
from utility import Utility

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)


class Compare:
    def __init__(self) -> None:
        pass
    def __init__(self, filename_mc, decay, meson): 
        self.meson = meson
        self.decay = decay
        self.filename_mc = TFile.Open(filename_mc, "READ")
    
    def plot_eff_acc_dalitz_vs_pcm(self, filename_mc, config, raw_yield_mc_dalitz, raw_yield_mc_pcm, name_plot, save):
        utils_dalitz = Utility(self.meson, filename_mc, config, "dalitz", "mc")
        utils_pcm = Utility(self.meson, filename_mc, config, "pcm", "mc")
        
        eff_dalitz = Efficiency("dalitz", self.meson)
        eff_pcm = Efficiency("pcm", self.meson)
        h1 = eff_dalitz.eff_and_acc_vs_pt("false", raw_yield_mc_dalitz, utils_dalitz, "eff_acc_dal", "false")
        h2 = eff_pcm.eff_and_acc_vs_pt("false", raw_yield_mc_pcm, utils_dalitz, "eff_acc_pcm", "false")
        ratio = h1.Clone("Ratio_eff")
        ratio.Divide(h2)
        
        p_T = utils_dalitz.get_bin_var()
        for i in range(1, len(p_T)-1):
            x = h1.GetBinContent(i)
            y= h1.GetBinCenter(i)
            
            x2 = h2.GetBinContent(i)
            y2 = h2.GetBinCenter(i)
            

        max_h1 = h1.GetMaximum()
        max_h2 = h2.GetMaximum()
        max_y = 2.0 * max(max_h1, max_h2)  # 2 times higher than the biggest point

        h1.SetMaximum(max_y)


        # Make a Canvas
        c1 = ROOT.TCanvas("eff_acc_dal_vs_pcm", "eff_acc_dal_vs_pcm", 800, 800)
        c1.Divide(1, 2)

        p1 = c1.cd(1)
        p1.SetPad(0.0, 0.0, 1, 1)
        p1.SetMargin(0.15,0.02,0.22,0.02)
        p1.SetTicks(1,1)
        p1.SetLogy()
        p1.SetLogx()


        h1.SetLineColor(ROOT.kOrange)
        h1.GetYaxis().SetRangeUser(5e-8, 5e-4)
        h1.SetMarkerColor(ROOT.kOrange)
        h1.SetMarkerStyle(20)  
        h1.Draw("E1P")          

        
        h2.SetLineColor(ROOT.kViolet)
        h2.SetMarkerColor(ROOT.kViolet)
        h2.SetMarkerStyle(34)  
        h2.Draw("E1P same")     
        
        # Add legend
        legend = ROOT.TLegend(0.65, 0.45, 0.85, 0.35)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)
        if self.meson == "pi0":
            legend.AddEntry(h1, "#pi^{0} -> e^{+} e^{-} #gamma", "LP")
            legend.AddEntry(h2, "#pi^{0} -> #gamma #gamma", "LP")
        if self.meson == "eta":
            legend.AddEntry(h1, "#eta -> e^{+} e^{-} #gamma", "LP")
            legend.AddEntry(h2, "#eta -> #gamma #gamma", "LP")
        legend.Draw()

        # Add text 
        txt = TPaveText(0.3, 0.975, 0.3, 0.925,"NDC") 
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

        
        h1.SetTitle("#varepsilon x A x BR comparison dalitz vs. pcm")
        h1.SetXTitle("p_{T} (GeV/c)")
        h1.GetYaxis().SetTitleSize(0.04)
        h1.GetYaxis().SetTitleOffset(1)
        h1.SetYTitle("#varepsilon x A x BR")

        p2 = c1.cd(2)
        p2.SetPad(0,0,1,0.3)
        p2.SetMargin(0.15,0.02,0.22,0.0)
        p2.SetTicks(1,1)
        p2.SetLogx()
        p2.SetLogy()
        

        ratio.SetLineColor(ROOT.kViolet)
        ratio.SetMarkerColor(ROOT.kViolet)
        ratio.SetMarkerStyle(34)  
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
        
        ratio.GetYaxis().SetTitleSize(0.08)
        ratio.GetYaxis().SetTitleOffset(.6)
        ratio.SetYTitle("Ratio")

        c1.Update()

        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())

        
    def plot_acc_dalitz_vs_pcm(self, filename_mc, config_array_dalitz, config_array_pcm, name_plot, save):
        
        config_file_dalitz = config_array_dalitz[0]
        config_file_pcm = config_array_pcm[0]
        with open(config_file_dalitz, "r", encoding="utf-8") as config_dalitz_yml:
            config_dalitz = yaml.safe_load(config_dalitz_yml)
        with open(config_file_pcm, "r", encoding="utf-8") as config_pcm_yml:
            config_pcm = yaml.safe_load(config_pcm_yml)
        utils_pcm = Utility(filename_mc, config_pcm, "pcm")
        utils_dalitz = Utility(filename_mc, config_dalitz, "dalitz")
        acc_class_dalitz= HistogramAcceptance(filename_mc, config_dalitz, "dalitz")
        acc_class_pcm = HistogramAcceptance(filename_mc, config_pcm, "pcm")
        h1 = acc_class_dalitz.calc_acceptance_vs_pt("acceptance_dalitz", "false", utils_dalitz)
        h2 = acc_class_pcm.calc_acceptance_vs_pt("acceptance_dalitz", "false", utils_pcm)

        ratio = h1.Clone("Ratio_acc")
        ratio.Divide(h2)

        
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
        max_y = 1.5 * max(max_h1, max_h2)  

        h1.SetMaximum(max_y)

        
        c1 = ROOT.TCanvas("acc_dal_vs_pcm", "acc_dal_vs_pcm", 800, 600)
        c1.Divide(1, 2)
        pad = c1.cd(1)
        
        pad.SetPad(0.0, 0.0, 1, 1)
        pad.SetMargin(0.15,0.02,0.22,0.0)
        pad.SetTicks(1,1)
        pad.SetLogy()
        
        h1.SetLineColor(ROOT.kBlue)
        h1.SetMarkerColor(ROOT.kBlue)
        h1.SetMarkerStyle(20)  
        h1.Draw("E1P")         

       
        h2.SetLineColor(ROOT.kMagenta)
        h2.SetMarkerColor(ROOT.kMagenta)
        h2.SetMarkerStyle(34)  
        h2.Draw("E1P same")     

        legend = ROOT.TLegend(0.6, 0.8, 0.8, 0.7)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)
        legend.AddEntry(h1, "Acceptance for #pi^{0} -> e^{+} e^{-} #gamma", "LP")
        legend.AddEntry(h2, "Acceptance for #pi^{0} -> #gamma #gamma", "LP")
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

        
        h1.SetTitle("Acceptance comparison dalitz vs. pcm")
        h1.SetXTitle("p_{T} (GeV/c)")
        h1.SetYTitle("Acceptance")
        h1.GetXaxis().SetRangeUser(0,12)
        p2 = c1.cd(2)
        p2.SetPad(0,0,1,0.3)
        p2.SetMargin(0.15,0.02,0.22,0.0)
        p2.SetTicks(1,1)
        p2.SetLogy()

        ratio.SetLineColor(ROOT.kBlue)
        ratio.SetMarkerColor(ROOT.kBlue)
        ratio.SetMarkerStyle(21)  
        ratio.Draw("E1P")   
         
        legend2 = ROOT.TLegend(0.2,0.25,0.6,0.4)
        legend2.SetBorderSize(0)
        legend2.SetFillColor(kWhite)
        legend2.SetFillStyle(0)
        legend2.SetTextSize(0.05)
        legend2.AddEntry(ratio, "ratio dalitz/pcm", "LP")
        legend2.Draw()  

        ratio.SetTitle("Ratio Dalitz/pcm")
        ratio.SetXTitle("p_{T} (GeV/c)")
        ratio.GetXaxis().SetLabelSize(0.10)
        ratio.GetXaxis().SetTitleSize(0.10)
        ratio.GetXaxis().SetTitleOffset(1.0)
        ratio.GetXaxis().SetRangeUser(0,12)
        ratio.GetYaxis().SetLabelSize(0.10)
        
       
        ratio.GetYaxis().SetTitleSize(0.10)
        ratio.GetYaxis().SetTitleOffset(.5)
        ratio.SetYTitle("Ratio")

        c1.Update()

        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())
    

    def compare_corr_dalitz_pcm(self, corr_dalitz, corr_pcm, name_plot, save): 

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
        max_y = 1.5 * max(max_h1, max_h2)  
        min_y =  min(min_h1, min_h2)

        h1.SetMaximum(max_y)
        print("min_y :", min_y)
        
        c1 = ROOT.TCanvas("py_vs_cor", "py_vs_cor", 800, 600)
        c1.Divide(1, 2)

        p1 = c1.cd(1)
        p1.SetPad(0.0, 0.0, 1, 1)
        p1.SetMargin(0.15,0.02,0.22,0.0)
        p1.SetTicks(1,1)
        p1.SetLogy()
        p1.SetLogx()
        
        h1.SetLineColor(kOrange)
        h1.SetMarkerColor(kOrange)
        h1.SetMarkerStyle(20)
        if self.meson == "eta":
            h1.GetYaxis().SetRangeUser(5e-5, 8e-2)
        h1.Draw("E1P")

        h2.SetLineColor(kViolet)
        h2.SetMarkerColor(kViolet)
        h2.SetMarkerStyle(34)  
        h2.Draw("E1P same")     

        # Add legend
        legend1 = ROOT.TLegend(0.5, 0.75, 0.9, 0.85)
        legend1.SetBorderSize(0)
        legend1.SetFillColor(kWhite)
        legend1.SetFillStyle(0)
        legend1.SetTextSize(0.035)
        if self.meson == "pi0":
            legend1.AddEntry(h1, "corrected yield for #pi^{0} -> e^{+} e^{-} #gamma", "LP")
            legend1.AddEntry(h2, "corrected yield for #pi^{0} -> #gamma #gamma", "LP")
        if self.meson == "eta":
            legend1.AddEntry(h1, "corrected yield for #eta -> e^{+} e^{-} #gamma", "LP")
            legend1.AddEntry(h2, "corrected yield for #eta -> #gamma #gamma", "LP")
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

        txt3 = TPaveText(0.9,0.945,0.95,0.92,"NDC"); 
        txt3.SetFillColor(kWhite)
        txt3.SetFillStyle(0)
        txt3.SetBorderSize(0)
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.02)
        txt3.AddText("pp, #sqrt{s} = 13.6TeV")
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

       
        h1.SetTitle("Corrected yield")
        h1.SetXTitle("p_{T} (GeV/c)")
        h1.SetYTitle("#frac{1}{2#pi} #frac{d^{2}N}{dp_{T}dy} [GeV^{-1} c^{3}]")

        #Go to second canvas to plot ratio in
        p2 = c1.cd(2)
        p2.SetPad(0,0,1,0.3)
        p2.SetMargin(0.15,0.02,0.22,0.0)
        p2.SetTicks(1,1)
        p2.SetLogx()
        p2.SetLogy()
       

        ratio.SetLineColor(ROOT.kViolet)
        ratio.SetMarkerColor(ROOT.kViolet)
        ratio.SetMarkerStyle(34) 
        ratio.GetYaxis().SetNdivisions(507, True)
        ratio.Draw("E1P")   
         
        legend2 = ROOT.TLegend(0.2,0.9,0.3,0.8)
        legend2.SetBorderSize(0)
        legend2.SetFillColor(kWhite)
        legend2.SetFillStyle(0)
        legend2.SetTextSize(0.03)
        legend2.AddEntry(ratio, "ratio dalitz/pcm", "LP")
        legend2.Draw()  

        ratio.SetTitle("Ratio Dalitz/pcm")
        ratio.SetXTitle("p_{T} (GeV/c)")
        ratio.GetXaxis().SetLabelSize(0.10)
        ratio.GetXaxis().SetTitleSize(0.10)
        ratio.GetXaxis().SetTitleOffset(1.6)
        ratio.GetYaxis().SetLabelSize(0.10)
        
        if self.meson == "pi0":
            ratio.GetYaxis().SetRangeUser(-0.5, 2.)
        if self.meson == "eta":
            ratio.GetYaxis().SetRangeUser(1, 1.8)
        ratio.GetYaxis().SetTitleSize(0.10)
        ratio.GetYaxis().SetTitleOffset(.3)
        ratio.SetYTitle("Ratio")

        # Fit Ratio
        fit_result = ratio.Fit("pol0", "S", "", 0.6, 10.)
        fit = ratio.GetFunction("pol0")
        fit.SetLineColor(ROOT.kCyan)
        fit.SetLineWidth(1)  

        
        fit_constant = fit_result.Parameter(0)
        fit_error = fit_result.ParError(0)
        fit_info = f"Fit: y = {fit_constant:.3f} #pm {fit_error:.3f}"

        legend2 = ROOT.TLegend(0.15,0.94,0.5,0.75)
        legend2.SetBorderSize(0)
        legend2.SetFillColor(ROOT.kWhite)
        legend2.SetFillStyle(0)
        legend2.SetTextSize(0.05)
        legend2.AddEntry(ratio, "Ratio Dalitz/PCM", "LP")
        
        placeholder_graph = ROOT.TGraph()
        placeholder_graph.SetTitle(fit_info)  
        legend2.AddEntry(fit, fit_info, "l")  
        legend2.Draw()  

        ratio.SetTitle("Ratio Dalitz/PCM")
        ratio.SetXTitle("p_{T} (GeV/c)")
        ratio.GetXaxis().SetLabelSize(0.09)
        ratio.GetXaxis().SetTitleSize(0.10)
        ratio.GetXaxis().SetTitleOffset(1.0)
        ratio.GetYaxis().SetLabelSize(0.10)
        ratio.GetYaxis().SetRangeUser(0.5, 2.8)
        ratio.GetYaxis().SetTitleSize(0.10)
        ratio.GetYaxis().SetTitleOffset(.3)
        ratio.SetYTitle("Ratio")

        
        line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1, ratio.GetXaxis().GetXmax(), 1)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(2)  
        line.Draw()

        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())
        

    def plot_eff_times_acc(self, filename_mc, config, raw_mc_dalitz, raw_mc_pcm, name_plot, save, utils):
        acc_class_dalitz = HistogramAcceptance(filename_mc, config, "dalitz")
        acc_class_pcm = HistogramAcceptance(filename_mc, config, "pcm")
        eff_class_dalitz = Efficiency("dalitz")
        eff_class_pcm = Efficiency("pcm")
        h1_eff = eff_class_dalitz.calc_efficiency_vs_pt("Efficiency_dalitz", "false", raw_mc_dalitz, utils)
        h1_acc = acc_class_dalitz.calc_acceptance_vs_pt("Acceptance_dalitz", "false", utils)
        h2_eff = eff_class_pcm.calc_efficiency_vs_pt("Efficiency_pcm", "false", raw_mc_pcm, utils)
        h2_acc = acc_class_pcm.calc_acceptance_vs_pt("Acceptance_pcm", "false", utils)
        
        product_dalitz = h1_eff.Clone("Dalitz_product")
        product_pcm = h2_eff.Clone("PCM_product")
        product_dalitz.Multiply(h1_acc)
        product_pcm.Multiply(h2_acc)

        max_h1 = product_dalitz.GetMaximum()
        max_h2 = product_pcm.GetMaximum()
        max_y = 1.2 * max(max_h1, max_h2)  

        product_dalitz.SetMaximum(max_y)

       
        c1 = ROOT.TCanvas("product_dal_vs_pcm", "product_dal_vs_pcm", 800, 600)
        pad = c1.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()

       
        product_dalitz.SetLineColor(ROOT.kBlue)
        product_dalitz.SetMarkerColor(ROOT.kBlue)
        product_dalitz.SetMarkerStyle(21)  
        product_dalitz.Draw("E1P")         

        # second histogram
        product_pcm.SetLineColor(ROOT.kMagenta)
        product_pcm.SetMarkerColor(ROOT.kMagenta)
        product_pcm.SetMarkerStyle(22)  
        product_pcm.Draw("E1P same")    

        
        legend = ROOT.TLegend(0.2, 0.25, 0.4, 0.15)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)
        if self.meson == "pi0":
            legend.AddEntry(product_dalitz, "Efficiency times Acceptance for #pi^{0} -> e^{+} e^{-} #gamma", "LP")
            legend.AddEntry(product_pcm, "Efficiency times Acceptancefor #pi^{0} -> #gamma #gamma", "LP")
        if self.meson == "eta":
            legend.AddEntry(product_dalitz, "Efficiency times Acceptance for #pi^{0} -> e^{+} e^{-} #gamma", "LP")
            legend.AddEntry(product_pcm, "Efficiency times Acceptancefor #pi^{0} -> #gamma #gamma", "LP")
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

        
        product_dalitz.SetTitle("Efficiency times Acceptance")
        product_dalitz.SetXTitle("p_{T} (GeV/c)")
        product_dalitz.SetYTitle("Efficiency * Acceptance")

        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())
        

    def get_13tev_file(self, filename_comp,hist_comp):
        self.filename_comp = TFile.Open(filename_comp, "READ")
        if self.meson == "pi0":
            self.dir_1 = self.filename_comp.Get("Pi013TeV") 
        if self.meson == "eta":
            self.dir_1 = self.filename_comp.Get("Eta13TeV")
        self.hist_comp = self.dir_1.Get(hist_comp) 
        if not self.hist_comp:
            print ("histogram", hist_comp,"not found")
        if self.hist_comp:
            print("histogram found", hist_comp)
            print("typ von hist: ", type(self.hist_comp))
        return self.hist_comp

    
    def compare_13tev(self, filename_comp, comp_partner, comp_partner_hist, decay, name_plot, save):
        if comp_partner == "acc":
            hist_comp = "AcceptancePi0"
            
        if comp_partner == "eff":
            hist_comp = "EfficiencyPi0"
        
        if comp_partner == "corr":
            hist_comp = "CorrectedYieldPi0"

        if comp_partner == "eff_times_acc":
            hist_comp = "EffTimesAccPi0"

        h1 = comp_partner_hist
        if h1: 
            print("histogram loaded successfully", type(h1))
        else:
            print("histogram not found")
        
        graph = self.get_13tev_file(filename_comp, hist_comp)
        if comp_partner == "eff_times_acc":
            if decay == "pcm":
                graph.Scale(0.98823)
                print("typ von graph mit pcm:", type(graph))
            if decay == "dalitz":
                graph.Scale(0.01174)
                print("typ von graph mit dalitz:", type(graph))
                
        if graph:
            print("graph loaded successfully", type(graph))
            print("histogramm in graph: ", hist_comp)
        if not graph:
            print("graph not found")

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

        max_y = 1.5 * max(max_h1, max_h2)  

        h1.SetMaximum(max_y)
        

        
        if comp_partner == "acc":
            c1 = ROOT.TCanvas("comp_run2_acc", "comp_run2_acc", 800, 600)
        if comp_partner == "eff":
            c1 = ROOT.TCanvas("comp_run2_eff", "comp_run2_eff", 800, 600)
        if comp_partner == "corr":
            c1 = ROOT.TCanvas("comp_run2_corr", "comp_run2_corr", 800, 600)
        if comp_partner == "eff_times_acc":
            c1 = ROOT.TCanvas("comp_run2_eff_times_acc", "comp_run2_eff_times_acc", 800, 600)
        

        p1 = c1.cd()
        p1.SetPad(0.0, 0.0, 1, 1)
        p1.SetMargin(0.15,0.02,0.22,0.0)
        p1.SetTicks(1,1)
        p1.SetLogy()
        
        
        h1.SetLineColor(ROOT.kBlue)
        h1.SetMarkerColor(ROOT.kBlue)
        h1.SetMarkerStyle(34)
        h1.Draw("E1P")          

        graph.SetLineColor(ROOT.kMagenta)
        graph.SetMarkerColor(ROOT.kMagenta)
        graph.SetMarkerStyle(34)  
        
        if comp_partner == "corr":
            graph.Draw("E1P Same")
        else:
            graph.Draw("P SAME")
        

        
        if comp_partner=="acc":
            legend1 = ROOT.TLegend(0.4, 0.4, 0.95, 0.3)
        if comp_partner == "eff":
            legend1 = ROOT.TLegend(0.4, 0.4, 0.95, 0.3)
        if comp_partner == "corr":
            legend1 = ROOT.TLegend(0.4, 0.9, 0.95, 0.75)
        if comp_partner=="eff_times_acc":
            legend1 = ROOT.TLegend(0.4, 0.4, 0.95, 0.3)
        legend1.SetBorderSize(0)
        legend1.SetFillColor(kWhite)
        legend1.SetFillStyle(0)
        legend1.SetTextSize(0.03)
        if comp_partner == "acc":
            if decay == "dalitz":
                if self.meson == "pi0":
                    legend1.AddEntry(h1, "acceptance for #pi^{0} -> e^{+} e^{-} #gamma at #sqrt{s} = 13.6TeV", "LP")
                    legend1.AddEntry(graph, "acceptance for #pi^{0} -> e^{+} e^{-} #gamma at #sqrt{s} = 13.0TeV", "LP")
                if self.meson == "eta":
                    legend1.AddEntry(h1, "acceptance for #eta -> e^{+} e^{-} #gamma at #sqrt{s} = 13.6TeV", "LP")
                    legend1.AddEntry(graph, "acceptance for #eta -> e^{+} e^{-} #gamma at #sqrt{s} = 13.0TeV", "LP")
            if decay == "pcm":
                if self.meson == "pi0":
                    legend1.AddEntry(h1, "acceptance for #pi^{0} -> #gamma #gamma at #sqrt{s} = 13.6TeV", "LP")
                    legend1.AddEntry(graph, "acceptance for #pi^{0} -> #gamma #gamma at #sqrt{s} = 13.0TeV", "LP")
                if self.meson == "eta":
                    legend1.AddEntry(h1, "acceptance for #eta -> #gamma #gamma at #sqrt{s} = 13.6TeV", "LP")
                    legend1.AddEntry(graph, "acceptance for #eta -> #gamma #gamma at #sqrt{s} = 13.0TeV", "LP")
        if comp_partner == "eff":
            if decay == "dalitz":
                if self.meson == "pi0":
                    legend1.AddEntry(h1, "efficiency for #pi^{0} -> e^{+} e^{-} #gamma at #sqrt{s} = 13.6TeV", "LP")
                    legend1.AddEntry(graph, "efficiency for #pi^{0} -> e^{+} e^{-} #gamma at #sqrt{s} = 13.0TeV", "LP")
                if self.meson == "eta":
                    legend1.AddEntry(h1, "efficiency for #eta -> e^{+} e^{-} #gamma at #sqrt{s} = 13.6TeV", "LP")
                    legend1.AddEntry(graph, "efficiency for #eta -> e^{+} e^{-} #gamma at #sqrt{s} = 13.0TeV", "LP")
            if decay == "pcm":
                if self.meson == "pi0":
                    legend1.AddEntry(h1, "efficiency for #pi^{0} -> #gamma #gamma at #sqrt{s} = 13.6TeV", "LP")
                    legend1.AddEntry(graph, "efficiency for #pi^{0} -> #gamma #gamma at #sqrt{s} = 13.0TeV", "LP")
                if self.meson == "eta":
                    legend1.AddEntry(h1, "efficiency for #eta -> #gamma #gamma at #sqrt{s} = 13.6TeV", "LP")
                    legend1.AddEntry(graph, "efficiency for #eta -> #gamma #gamma at #sqrt{s} = 13.0TeV", "LP")
        if comp_partner == "corr":
            if decay == "dalitz":
                if self.meson == "pi0":
                    legend1.AddEntry(h1, "corrected yield for #pi^{0} -> e^{+} e^{-} #gamma at #sqrt{s} = 13.6TeV", "LP")
                    legend1.AddEntry(graph, "corrected yield for #pi^{0} -> e^{+} e^{-} #gamma at #sqrt{s} = 13.0TeV", "LP")
                if self.meson == "eta":
                    legend1.AddEntry(h1, "corrected yield for #eta -> e^{+} e^{-} #gamma at #sqrt{s} = 13.6TeV", "LP")
                    legend1.AddEntry(graph, "corrected yield for #eta -> e^{+} e^{-} #gamma at #sqrt{s} = 13.0TeV", "LP")
            if decay == "pcm":
                if self.meson == "pi0":
                    legend1.AddEntry(h1, "corrected yield for #pi^{0} -> #gamma #gamma at #sqrt{s} = 13.6TeV", "LP")
                    legend1.AddEntry(graph, "corrected yield for #pi^{0} -> #gamma #gamma at #sqrt{s} = 13.0TeV", "LP")
                if self.meson == "eta":
                    legend1.AddEntry(h1, "corrected yield for #eta -> #gamma #gamma at #sqrt{s} = 13.6TeV", "LP")
                    legend1.AddEntry(graph, "corrected yield for #eta -> #gamma #gamma at #sqrt{s} = 13.0TeV", "LP")
        if comp_partner == "eff_times_acc":
            if decay == "dalitz":
                if self.meson == "pi0":
                    legend1.AddEntry(h1, "#varepsilon x A x BR for #pi^{0} -> e^{+} e^{-} #gamma at #sqrt{s} = 13.6TeV", "LP")
                    legend1.AddEntry(graph, "#varepsilon x A x BR for #pi^{0} -> e^{+} e^{-} #gamma at #sqrt{s} = 13.0TeV", "LP")
                if self.meson == "eta":
                    legend1.AddEntry(h1, "#varepsilon x A x BR for #eta -> e^{+} e^{-} #gamma at #sqrt{s} = 13.6TeV", "LP")
                    legend1.AddEntry(graph, "#varepsilon x A x BR for #eta -> e^{+} e^{-} #gamma at #sqrt{s} = 13.0TeV", "LP")
            if decay == "pcm":
                if self.meson == "pi0":
                    legend1.AddEntry(h1, "#varepsilon x A x BR for #pi^{0} -> #gamma #gamma at #sqrt{s} = 13.6TeV", "LP")
                    legend1.AddEntry(graph, "#varepsilon x A x BR for #pi^{0} -> #gamma #gamma at #sqrt{s} = 13.0TeV", "LP")
                if self.meson == "eta":
                    legend1.AddEntry(h1, "#varepsilon x A x BR for #eta -> #gamma #gamma at #sqrt{s} = 13.6TeV", "LP")
                    legend1.AddEntry(graph, "#varepsilon x A x BR for #eta -> #gamma #gamma at #sqrt{s} = 13.0TeV", "LP")

        
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

        txt3 = TPaveText(0.9,0.945,0.95,0.92,"NDC"); 
        txt3.SetFillColor(kWhite)
        txt3.SetFillStyle(0)
        txt3.SetBorderSize(0)
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.02)
        txt3.AddText("pp, #sqrt{s} = 13.6TeV")
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

        
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
            h1.SetYTitle("#frac{1}{2#pip_{T}} #frac{d^{2}N}{dp_{T}dy} [GeV^{-2} c^{3}]")
        if comp_partner == "eff_times_acc":
            h1.SetTitle("#varepsilon x A x BR")
            h1.SetXTitle("p_{T} (GeV/c)")
            h1.SetYTitle("#varepsilon x A x BR")
        x = h1.GetXaxis()
        x.SetTitleSize(0.04)
        x.SetTitleFont(42)
        x.SetTitleOffset(0.9)
        x.SetLabelFont(42)
        x.SetLabelSize(0.035)

        y = h1.GetYaxis()
        y.SetTitleSize(0.04)
        y.SetTitleFont(42)
        y.SetTitleOffset(1.3)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)


        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())
    
    

    def compare_13tev_eff_acc_new(self, filename_comp, comp_partner_hist, decay, name_plot, save):
        
        h1 = comp_partner_hist
        if h1: 
            print("histogram loaded successfully", type(h1))
        else:
            print("histogram not found")

        filename_comp = TFile.Open(filename_comp, "READ")
        dir_1 = self.filename_comp.Get("Pi013TeV") 
        
        hist_eff = dir_1.Get("EfficiencyPi0_INT7").Clone("eff_clone")
        hist_acc = dir_1.Get("AcceptancePi0_INT7").Clone("acc_clone")
        
        
        
        
        print( "Type of eff graph: ", type(hist_eff))
        print("Type of acc graph after loading eff: ", type(hist_acc))

        eff_acc_hist = hist_acc.Clone("result_hist")
        eff_acc_hist.Reset()
        for bin in range (1, hist_eff.GetNbinsX()+1):
            content_eff = hist_eff.GetBinContent(bin)
            if self.decay == "pcm":
                content_acc = hist_acc.GetBinContent(bin) 
            if self.decay == "dalitz":
                content_acc = hist_acc.GetBinContent(bin) * (1/100)
            mult_content = content_eff * content_acc
            eff_acc_hist.SetBinContent(bin, mult_content)

    
        max_h1 = h1.GetMaximum()

        min_h2 = eff_acc_hist.GetMinimum()
        max_h2 = eff_acc_hist.GetMaximum()
       
        max_y = 1.5 * max(max_h1, max_h2)  

        h1.SetMaximum(max_y)
        

    
        c1 = ROOT.TCanvas("comp_run2_eff_times_acc", "comp_run2_eff_times_acc", 800, 600)
        

        p1 = c1.cd()
        p1.SetPad(0.0, 0.0, 1, 1)
        p1.SetMargin(0.15,0.02,0.22,0.0)
        p1.SetTicks(1,1)
        p1.SetLogy()
        
        
        h1.SetLineColor(ROOT.kBlue)
        h1.SetMarkerColor(ROOT.kBlue)
        h1.SetMarkerStyle(34)
        h1.Draw("E1P")          

        eff_acc_hist.SetLineColor(ROOT.kMagenta)
        eff_acc_hist.SetMarkerColor(ROOT.kMagenta)
        eff_acc_hist.SetMarkerStyle(34)  
        
        
        eff_acc_hist.Draw("E1P Same")
        
        legend1 = ROOT.TLegend(0.4, 0.4, 0.95, 0.3)
        legend1.SetBorderSize(0)
        legend1.SetFillColor(kWhite)
        legend1.SetFillStyle(0)
        legend1.SetTextSize(0.03)
       
        if decay == "dalitz":
            if self.meson == "pi0":
                legend1.AddEntry(h1, "#varepsilon x A x BR for #pi^{0} -> e^{+} e^{-} #gamma at #sqrt{s} = 13.6TeV", "LP")
                legend1.AddEntry(eff_acc_hist, "#varepsilon x A x BR for #pi^{0} -> e^{+} e^{-} #gamma at #sqrt{s} = 13.0TeV", "LP")
            if self.meson == "eta":
                legend1.AddEntry(h1, "#varepsilon x A x BR for #eta -> e^{+} e^{-} #gamma at #sqrt{s} = 13.6TeV", "LP")
                legend1.AddEntry(eff_acc_hist, "#varepsilon x A x BR for #eta -> e^{+} e^{-} #gamma at #sqrt{s} = 13.0TeV", "LP")
        if decay == "pcm":
            if self.meson == "pi0":
                legend1.AddEntry(h1, "#varepsilon x A x BR for #pi^{0} -> #gamma #gamma at #sqrt{s} = 13.6TeV", "LP")
                legend1.AddEntry(eff_acc_hist, "#varepsilon x A x BR for #pi^{0} -> #gamma #gamma at #sqrt{s} = 13.0TeV", "LP")
            if self.meson == "eta":
                legend1.AddEntry(h1, "#varepsilon x A x BR for #eta -> #gamma #gamma at #sqrt{s} = 13.6TeV", "LP")
                legend1.AddEntry(eff_acc_hist, "#varepsilon x A x BR for #eta -> #gamma #gamma at #sqrt{s} = 13.0TeV", "LP")

        
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

        txt3 = TPaveText(0.9,0.945,0.95,0.92,"NDC"); 
        txt3.SetFillColor(kWhite)
        txt3.SetFillStyle(0)
        txt3.SetBorderSize(0)
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.02)
        txt3.AddText("pp, #sqrt{s} = 13.6TeV")
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

        
        h1.SetTitle("#varepsilon x A x BR")
        h1.SetXTitle("p_{T} (GeV/c)")
        h1.SetYTitle("#varepsilon x A x BR")
        x = h1.GetXaxis()
        x.SetTitleSize(0.04)
        x.SetTitleFont(42)
        x.SetTitleOffset(0.9)
        x.SetLabelFont(42)
        x.SetLabelSize(0.035)

        y = h1.GetYaxis()
        y.SetTitleSize(0.04)
        y.SetTitleFont(42)
        y.SetTitleOffset(1.2)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)


        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())


    def compare_cuts(self, fHistolist, cutnames, decay, filename_mc, config, name_plot, comp_partner, utils):
        TGaxis.SetMaxDigits(3);
        pt = utils.get_bin_var()
        npt = len(pt);

        titlePt = fHistolist[0].GetTitle(); 
        if comp_partner == "eff":
            yMin_, yMax_ = 1e-5, 1.4*1e-2
        if comp_partner == "corr":
            yMin_, yMax_ = 0 , 100
        if comp_partner == "eff_acc":
            yMin_, yMax_ = 1e-7, 1e-4
        color = [kRed+1, kBlue+1, kGreen+2, kMagenta+2, kCyan+1]

        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);
        p1.SetLogy();


        frame1 = p1.DrawFrame(0., yMin_, 12., yMax_);
        if comp_partner == "eff_acc":
            FrameSettings(frame1, "", "#varepsilon x A x BR", 0., 12.)
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
            fHistolist[icut].SetMarkerSize(1.5)
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
        frame2 = p2.DrawFrame(0., 0.5, 12., 1.5);0., 5*1e-3, 12., 1e1
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
        if comp_partner == "eff_acc":
            txt.AddText("#varepsilon x A x BR as a function of #it{p}_{T}")
        if comp_partner == "eff":
            txt.AddText("Efficiency as a function of #it{p}_{T}");
        if comp_partner == "corr":
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

        if comp_partner == "eff_acc":
            leg = TLegend(0.5, 0.35, 0.95, 0.45)
        else:
            leg = TLegend(0.65,0.62,0.95,0.77);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.02);
        leg.SetTextFont(42);#helvetica
        for icut in range(len(fHistolist)):
            leg.AddEntry(fHistolist[icut],"{}".format(cutnames[icut]),"ep");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        c1.SaveAs(name_plot.Data());
        c1.Close();   

    def corr_yield_ratio(self, corr_eta, corr_pi0, name_plot, save):
        print("corrected yield of the eta loaded successfully", type(corr_eta))
        print("corrected yield of the pi0 loaded successfully", type(corr_pi0))

        ratio = corr_eta.Clone("Ratio_corr_eta_pi0")
        ratio.Divide(corr_pi0)

        ratio.SetFillColor(kCyan+1)
        ratio.SetMarkerStyle(kFullCross)
        ratio.SetMarkerColor(kBlue)
        ratio.SetLineColor(kBlue)
        
        

        # new canvas
        canvas = TCanvas("ratio_eta_pi0", "ratio_eta_pi0", 0, 0, 900, 900)
        pad = canvas.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.2, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()
        
        #Draw points in the histogramm
        ratio.Draw("Esame")
       
        # Define y-axis settings
        y = ratio.GetYaxis() 
        y.SetTitle("#eta / #pi^{0}")
        y.SetTitleSize(0.048)
        y.SetTitleFont(42)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)
        

        #Define x-axis settings
        x = ratio.GetXaxis()
        x.SetTitle("p_{T} [GeV/c]")
        x.SetTitleSize(0.048)
        x.SetTitleFont(42)
        x.SetTitleOffset(0.9)
        x.SetLabelFont(42)
        x.SetLabelSize(0.035)


        # Define legend settings
        leg = TLegend(0.4, 0.75, 1.0, 0.7)
        leg.SetBorderSize(0)
        leg.SetFillColor(kWhite)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.03)
        leg.AddEntry(ratio, "corrected yield ratio #eta / #pi^{0}", "LP")
        leg.Draw("")
        ROOT.SetOwnership(leg,False)

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

        txt2 = TPaveText(0.875, 0.79, 0.875, 0.79, "NDC")
        txt2.SetFillColor(kWhite)
        txt2.SetFillStyle(0)
        txt2.SetBorderSize(0)
        txt2.SetTextAlign(33);#middle,left
        txt2.SetTextFont(42);#helvetica
        txt2.SetTextSize(0.02)
        txt2.AddText("pcm")
        txt2.Draw()
        ROOT.SetOwnership(txt2,False)
        

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


        #save histogram
        canvas.Modified()
        canvas.Update()
        ROOT.SetOwnership(canvas, False)
        if save == "true":
            canvas.SaveAs(name_plot.Data())
        return ratio

    

    def comp_eff_Teff(self, utils, raw_mc, h_eff_acc, h_Teff_acc, name_plot, save):
        h1 = h_eff_acc
        h2 = h_Teff_acc

        # for i in range(1,11):
        #     print("## efficiency from TEfficiency: ", i, h2.GetEfficiency(i), h2.GetEfficiencyErrorLow(i), h2.GetEfficiencyErrorUp(i))
        #     print("## efficiency: ", i, h1.GetBinContent(i), h1.GetBinError(i))
        
        p = utils.get_bin_var()
    
        max_h1 = h1.GetMaximum()
        max_y = 1.5 * max_h1

        h1.SetMaximum(max_y)

        
        c1 = ROOT.TCanvas("eff_vs_teff", "eff_vs_teff", 800, 600)
        
        pad = c1.cd()
        
        pad.SetPad(0.0, 0.0, 1, 1)
        pad.SetMargin(0.15,0.02,0.22,0.0)
        pad.SetTicks(1,1)
        pad.SetLogy()
         
        h2.SetLineColor(ROOT.kMagenta)
        h2.SetMarkerColor(ROOT.kMagenta)
        h2.SetMarkerStyle(22)  
        h2.Draw("E1P")  
        

        h1.SetLineColor(ROOT.kBlue)
        h1.SetMarkerColor(ROOT.kBlue)
        h1.SetMarkerStyle(21)  
        h1.Draw("E1P Same")      

        legend = ROOT.TLegend(0.7, 0.8, 0.9, 0.7)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)
        legend.AddEntry(h1, "eff x acc x BR", "LP")
        legend.AddEntry(h2, "Teff x acc x BR", "LP")
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

        
        h1.SetTitle("Comparison efficiency vs Teff")
        h1.SetXTitle("p_{T} (GeV/c)")
        h1.SetYTitle("Efficiency x acc x BR")
        h1.GetXaxis().SetRangeUser(0,12)

        c1.Update()

        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())


    def compare_13tev_etapi0(self, filename_comp, comp_partner_hist, name_plot, save):
        comp_partner= "eta_pi0"
        hist_comp = "EtaToPi0YShiftedStatError"
        h1 = comp_partner_hist
        if h1: 
            print("histogram loaded successfully", type(h1))
        else:
            print("histogram not found")
        
        graph = self.get_13tev_file(filename_comp, hist_comp)
       
                
        if graph:
            print("graph loaded successfully", type(graph))
            print("histogramm in graph: ", hist_comp)
        if not graph:
            print("graph not found")

        max_h1 = h1.GetMaximum()
        min_h1 = h1.GetMinimum()
        
        c1 = ROOT.TCanvas("comp_run2_eta_pi0", "comp_run2_eta_pi0", 800, 600)
        

        p1 = c1.cd()
        p1.SetPad(0.0, 0.0, 1, 1)
        p1.SetMargin(0.15,0.02,0.22,0.0)
        p1.SetTicks(1,1)
        p1.SetLogy()
        
        
        h1.SetLineColor(ROOT.kBlue)
        h1.SetMarkerColor(ROOT.kBlue)
        h1.SetMarkerStyle(34)
        h1.Draw("E1P")          

        graph.SetLineColor(ROOT.kMagenta)
        graph.SetMarkerColor(ROOT.kMagenta)
        graph.SetMarkerStyle(34)  
        graph.SetMarkerSize(1.5)
        
        graph.Draw("E1P Same")
                
        legend1 = ROOT.TLegend(0.4, 0.4, 0.95, 0.3)
        legend1.SetBorderSize(0)
        legend1.SetFillColor(kWhite)
        legend1.SetFillStyle(0)
        legend1.SetTextSize(0.03)
        legend1.AddEntry(h1, "#eta / #pi^{0} ratio at #sqrt{s} = 13.6 TeV", "LP")
        legend1.AddEntry(graph, "#eta /#pi^{0} ratio at #sqrt{s}=13.0 TeV", "LP")
        
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

        txt3 = TPaveText(0.9,0.945,0.95,0.92,"NDC"); 
        txt3.SetFillColor(kWhite)
        txt3.SetFillStyle(0)
        txt3.SetBorderSize(0)
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.02)
        txt3.AddText("pp, #sqrt{s} = 13.6TeV")
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

        
        
        h1.SetTitle("#eta/ #pi^{0}")
        h1.SetXTitle("p_{T} (GeV/c)")
    
        x = h1.GetXaxis()
        x.SetTitleSize(0.04)
        x.SetTitleFont(42)
        x.SetTitleOffset(0.9)
        x.SetLabelFont(42)
        x.SetLabelSize(0.035)

        y = h1.GetYaxis()
        y.SetTitleSize(0.04)
        y.SetTitleFont(42)
        y.SetTitleOffset(1.)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)
        y.SetRangeUser(4e-2, 0.9)


        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())

    def compare_eff_eta_pi0(self, eff_eta, eff_pi0, decay, name_plot, save): 

        h1 = eff_eta 
        h2 = eff_pi0
        print("h1 typ: ", type(h1))
        print("h2 typ: ", type(h2))

        ratio = h1.Clone("Ratio_eff_eta_pi0")
        ratio.Divide(h2)

        max_h1 = h1.GetMaximum()
        min_h1 = h1.GetMinimum()
        max_h2 = h2.GetMaximum()
        min_h2 = h2.GetMinimum()
        max_y = 1.5 * max(max_h1, max_h2)  
        min_y =  min(min_h1, min_h2)

        h1.SetMaximum(max_y)
        print("min_y :", min_y)
        
        c1 = ROOT.TCanvas("eff_eta_pi0", "eff_eta_pi0", 800, 600)
        c1.Divide(1, 2)

        p1 = c1.cd(1)
        p1.SetPad(0.0, 0.0, 1, 1)
        p1.SetMargin(0.15,0.02,0.22,0.0)
        p1.SetTicks(1,1)
        p1.SetLogy()
        
        h1.SetLineColor(kRed)
        h1.SetMarkerColor(kRed)
        h1.SetMarkerStyle(20)
        h1.Draw("E1P")          

        
        
        h2.SetLineColor(kBlack)
        h2.SetMarkerColor(kBlack)
        h2.SetMarkerStyle(34)  
        h2.Draw("E1P same")     

        # Add legend
        legend1 = ROOT.TLegend(0.5, 0.75, 0.9, 0.85)
        legend1.SetBorderSize(0)
        legend1.SetFillColor(kWhite)
        legend1.SetFillStyle(0)
        legend1.SetTextSize(0.035)
        
        if decay == "dalitz":
            legend1.AddEntry(h1, "#varepsilon x A x BR #eta -> e^{+} e^{-} #gamma", "LP")
            legend1.AddEntry(h2, "#varepsilon x A x BR #pi^{0} -> e^{+} e^{-} #gamma", "LP")
        if decay == "pcm":
            legend1.AddEntry(h1, "#varepsilon x A x BR #eta -> #gamma #gamma", "LP")
            legend1.AddEntry(h2, "#varepsilon x A x BR #pi^{0} -> #gamma #gamma", "LP")
        
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

        txt3 = TPaveText(0.9,0.945,0.95,0.92,"NDC"); 
        txt3.SetFillColor(kWhite)
        txt3.SetFillStyle(0)
        txt3.SetBorderSize(0)
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.02)
        txt3.AddText("pp, #sqrt{s} = 13.6TeV")
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

       
        h1.SetTitle("Comparison #varepsilon x A x BR")
        h1.SetXTitle("p_{T} (GeV/c)")
        h1.SetYTitle("#varepsilon x A x BR")

        #Go to second canvas to plot ratio in
        p2 = c1.cd(2)
        p2.SetPad(0,0,1,0.3)
        p2.SetMargin(0.15,0.02,0.22,0.0)
        p2.SetTicks(1,1)
        p2.SetLogy()
       

        ratio.SetLineColor(ROOT.kBlack)
        ratio.SetMarkerColor(ROOT.kBlack)
        ratio.SetMarkerStyle(34) 
        ratio.GetYaxis().SetNdivisions(507, True)
        ratio.Draw("E1P")   
         
        legend2 = ROOT.TLegend(0.2,0.9,0.3,0.8)
        legend2.SetBorderSize(0)
        legend2.SetFillColor(kWhite)
        legend2.SetFillStyle(0)
        legend2.SetTextSize(0.03)
        legend2.AddEntry(ratio, "ratio #eta / #pi^{0}", "LP")
        legend2.Draw()  

        ratio.SetTitle("Ratio #eta/#pi^{0}")
        ratio.SetXTitle("p_{T} (GeV/c)")
        ratio.GetXaxis().SetLabelSize(0.10)
        ratio.GetXaxis().SetTitleSize(0.10)
        ratio.GetXaxis().SetTitleOffset(1.6)
        ratio.GetYaxis().SetLabelSize(0.10)
        
        ratio.GetYaxis().SetTitleSize(0.10)
        ratio.GetYaxis().SetTitleOffset(.3)
        ratio.SetYTitle("Ratio")

        legend2 = ROOT.TLegend(0.15,0.94,0.5,0.75)
        legend2.SetBorderSize(0)
        legend2.SetFillColor(ROOT.kWhite)
        legend2.SetFillStyle(0)
        legend2.SetTextSize(0.05)
        legend2.AddEntry(ratio, "Ratio #eta/#pi^{0}", "LP") 
        legend2.Draw()  

        ratio.SetTitle("Ratio #eta/#pi^{0}")
        ratio.SetXTitle("p_{T} (GeV/c)")
        ratio.GetXaxis().SetLabelSize(0.09)
        ratio.GetXaxis().SetTitleSize(0.10)
        ratio.GetXaxis().SetTitleOffset(1.0)
        ratio.GetYaxis().SetLabelSize(0.10)
        ratio.GetYaxis().SetTitleSize(0.10)
        ratio.GetYaxis().SetTitleOffset(.3)
        ratio.SetYTitle("Ratio")


        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())
    





