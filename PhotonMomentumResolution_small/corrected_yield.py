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
from utility import Utility

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

class CorrectedYield:
    
    def __init__(self, filename_mc, decay): 
        self.filename_mc = TFile.Open(filename_mc, "READ")
        self.decay = decay

    
    
    def calculate_corrected_yield(self, name_plot, filename_mc, config, raw_hist_data, raw_hist_mc, save):
        if self.decay == "pcm":
            acc_class = HistogramAcceptance(filename_mc, config, "pcm")
            eff_class = Efficiency("pcm")
            utils = Utility(filename_mc, config, "pcm")
        else:
            acc_class = HistogramAcceptance(filename_mc, config, "dalitz")
            eff_class = Efficiency("dalitz")
            utils = Utility(filename_mc, config, "dalitz")

        deltaRapid = 1.8

        #Acceptance 
        acc = acc_class.calc_acceptance_vs_pt("Acceptance", "false")
        print("Typ der Akzeptanz: ", type(acc))
        
        # transveres momentum 
        p_T = utils.get_bin_var()

        #Efficiency
        eff = eff_class.calc_efficiency_vs_pt("Efficiency", "false", raw_hist_mc, acc_class)
        
        raw_data = raw_hist_data
        print("Raw data type: ", type(raw_data))
        
        #produce new empty canvas and pad to plot the corrected yield in 
        canvas = TCanvas("corr", "corr", 0, 0, 900, 900)
        pad = canvas.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)

        
        h_corrected_yield = raw_data.Clone("corrected_yield_hist") 
        

        # get sum over squared weights
        h_corrected_yield.Sumw2()


        h_corrected_yield.Scale(1 / (2 * np.pi))
        h_corrected_yield.Scale(1 / deltaRapid)
        h_corrected_yield.Divide(eff)
        h_corrected_yield.Divide(acc)
    
        #return corrected_yield_hist
        canvas = TCanvas("corr", "corr", 0, 0, 900, 900)
        pad = canvas.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()
        pad.SetLogx()

        h_corrected_yield.SetFillColor(kCyan+1)
        h_corrected_yield.SetMarkerStyle(kFullCross)
        h_corrected_yield.SetMarkerColor(kBlue)
        h_corrected_yield.SetMarkerSize(0.85)
        h_corrected_yield.SetLineColor(kBlue)
        
        h_corrected_yield.Draw("Esame")
       
        # Define y-axis settings
        y = h_corrected_yield.GetYaxis() 
        
        y.SetTitle("#frac{1}{2#pip_{T}} #frac{d^{2}N}{dp_{T}dy}")
        y.SetTitleSize(0.048)
        y.SetTitleFont(42)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)

        #Define x-axis settings
        x = h_corrected_yield.GetXaxis()
        x.SetTitle("p_{T} [GeV/c]")
        x.SetTitleSize(0.048)
        x.SetTitleFont(42)
        x.SetTitleOffset(0.9)
        x.SetLabelFont(42)
        x.SetLabelSize(0.035)

        # Define legend settings
        leg = TLegend(0.6, 0.8, 1.0, 0.7)
        leg.SetBorderSize(0)
        leg.SetFillColor(kWhite)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.03)
        
        if self.decay=="pcm":
            leg.AddEntry(h_corrected_yield, "\pi_{0} \\rightarrow \gamma \gamma", "LP")
        else:
            leg.AddEntry(h_corrected_yield, "\pi_{0} \\rightarrow e^{+} e^{-} \gamma", "LP")
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
        if save =="true":
            canvas.SaveAs(name_plot.Data())
        return h_corrected_yield
    
    
    def calc_corr_vs_pt(self, name_plot, save, raw_data, decay, raw_mc, filename_mc, config, eff_hist, acc_hist): # eff_hist, acc_hist, raw_hist_data, filename_mc, config):
         #Get wanted histograms:
        deltay = 1.8
        raw = raw_data
        acc_class = HistogramAcceptance(filename_mc, config, "dalitz")
            
        #get canvas and pad to plot the acceptance in 
        canvas_loc = TCanvas("corr_yield_loc", "corr_yield_loc", 0, 0, 900, 900)
        pad = canvas_loc.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        
        h_corrected_yield = raw.Clone("h_corr")
        
        h_corrected_yield.Sumw2()

        h_corrected_yield.Scale(1/(2*np.pi))
        h_corrected_yield.Scale(1/deltay)
        
    
        # Define marker settings
        h_corrected_yield.SetFillColor(kCyan+1)
        h_corrected_yield.SetMarkerStyle(kFullCross)
        h_corrected_yield.SetMarkerColor(kBlue)
        h_corrected_yield.SetMarkerSize(0.85)
        h_corrected_yield.SetLineColor(kBlue)

        #Draw points in the histogramm
        h_corrected_yield.Draw("ESame")
        
        canvas_loc.Draw()
       
        # Define y-axis settings
        y = h_corrected_yield.GetYaxis() 
        y.SetTitle("Efficiency") 
        y.SetTitleSize(0.048)
        y.SetTitleFont(42)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)

        #Define x-axis settings
        x = h_corrected_yield.GetXaxis()
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
            leg.AddEntry(h_corrected_yield, "\pi^{0} \\rightarrow \gamma \gamma", "LP")
        else:
            leg.AddEntry(h_corrected_yield, "\pi^{0} \\rightarrow e^{+} e^{-} \gamma", "LP")
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
        canvas_loc.Modified()
        canvas_loc.Update()
        ROOT.SetOwnership(canvas_loc, False)
        if save == "true":
            
            canvas_loc.SaveAs(name_plot.Data())
    
        print("corr hist type:", type(h_corrected_yield))
        return h_corrected_yield
    

    def calc_corr_yield_vs_pt(self, name_plot, save, raw_data, eff_hist, acc_hist, acc_class):
        
        #Get wanted histograms:
        deltarap = 1.8 
        h1 = raw_data
        print("First Histogram loaded, raw_data", type(h1))
        eff = eff_hist
        print("second histogram loaded, efficiency", type(eff))
        acc = acc_hist
        print("third histogram loaded, acc", type(acc))
        if h1 is None:
            print(f"Failed to load histogram: hPt_Pi0")

    #get canvas and pad to plot the acceptance in 
        can = TCanvas("correc", "correc", 0, 0, 900, 900)
        pad = can.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()
        pad.SetLogx()
        

        h_corr = h1.Clone("h_c")   
    
        h_corr.Divide(eff)
        h_corr.Divide(acc)
        h_corr.Scale(1/deltarap)
        h_corr.Scale(1/(2* np.pi))
        
        
        # Define marker settings
        h_corr.SetFillColor(kCyan+1)
        h_corr.SetMarkerStyle(kFullCross)
        h_corr.SetMarkerColor(kBlue)
        h_corr.SetMarkerSize(0.85)
        h_corr.SetLineColor(kBlue)

        #Draw points in the histogramm
        h_corr.Draw("Esame")
       
        # Define y-axis settings
        y = h_corr.GetYaxis() 
        y.SetTitle("#frac{1}{2#pip_{T}} #frac{d^{2}N}{dp_{T}dy}") 
        y.SetTitleSize(0.048)
        y.SetTitleFont(42)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)

        #Define x-axis settings
        x = h_corr.GetXaxis()
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
            leg.AddEntry(h_corr, "\pi^{0} \\rightarrow \gamma \gamma", "LP")
        else:
            leg.AddEntry(h_corr, "\pi^{0} \\rightarrow e^{+} e^{-} \gamma", "LP")
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

        txt3 = TPaveText(0.35,0.82,0.35,0.82,"NDC")
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
        print("corrected yield hist type:", type(h_corr))
        return h_corr
    
    
    
    
    def plot_pythia(self, utils, name_plot, save): 
        if self.decay =="pcm":
            dir_name22 = "PCMPCM"    
        else:
            dir_name22 = "PCMDalitzEE"
        hist = self.filename_mc.Get("pi0eta-to-gammagamma-mc/Generated").FindObject(dir_name22).FindObject("hPt_Pi0")
        
        deltay = 1.8
        Nev = utils.get_Nev()
        print("Number of Events: ", Nev)
       
        h1 = hist.Clone("pythia")
        
        for bin in range (1, h1.GetXaxis().GetNbins()):
             bin_content = h1.GetBinContent(bin)
             bin_center = h1.GetBinCenter(bin)
             bin_width = h1.GetBinWidth(bin)
             bin_error = h1.GetBinError(bin)

             new_bin_content = bin_content / (bin_center *bin_width) 
             h1.SetBinContent(bin, new_bin_content)
             new_bin_error = bin_error / (bin_center * bin_width)
             h1.SetBinError(bin, new_bin_error)

        h1.Scale(1 / Nev)
        h1.Scale(1 / deltay)
        h1.Scale(1 / (2 * np.pi))
        

        #plot the pythia corrected yield
        canvas = TCanvas("pyt", "pyt", 0, 0, 900, 900)
        pad = canvas.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()

        #Draw points in the histogramm
        h1.Draw("Esame")
       
        # Define y-axis settings
        y = h1.GetYaxis() 
        y.SetTitle("#frac{1}{2#pip_{T}} #frac{d^{2}N}{dp_{T}dy}")
        y.SetTitleSize(0.048)
        y.SetTitleFont(42)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)

        #Define x-axis settings
        x = h1.GetXaxis()
        x.SetTitle("p_{T} [GeV/c]")
        x.SetTitleSize(0.048)
        x.SetTitleFont(42)
        x.SetTitleOffset(0.9)
        x.SetLabelFont(42)
        x.SetLabelSize(0.035)

        # Define legend settings
        leg = TLegend(0.45, 0.8, 1.0, 0.75)
        leg.SetBorderSize(0)
        leg.SetFillColor(kWhite)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.03)
        leg.AddEntry(hist, "corrected yield Pythia", "LP")
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
        return h1
    
    
    def compare_pythia_corr_yield(self, raw_data, eff_hist, acc_hist, utils, name_plot, save):
        if self.decay == "pcm":
            dir_name22 = "PCMPCM"
        else:
            dir_name22 = "PCMDalitzEE"
        h1 = self.calc_corr_yield_vs_pt("corrected_yield", "false", raw_data, eff_hist, acc_hist, utils)
        h2 = self.plot_pythia(utils, "pythia", "false")

        max_h1 = h1.GetMaximum()
        max_h2 = h2.GetMaximum()
        max_y = 1.05 * max(max_h1, max_h2)  

        h1.SetMaximum(max_y)

        c1 = ROOT.TCanvas("py_vs_cor", "py_vs_cor", 800, 600)
        pad = c1.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()

        
        h1.SetLineColor(ROOT.kBlue)
        h1.SetMarkerColor(ROOT.kBlue)
        h1.SetMarkerStyle(21)  
        h1.Draw("E1P")          

     
        h2.SetLineColor(ROOT.kMagenta)
        h2.SetMarkerColor(ROOT.kMagenta)
        h2.SetMarkerStyle(22)  
        h2.Draw("E1P same")     

        
        legend = ROOT.TLegend(0.5, 0.75, 0.95, 0.65)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)
        legend.AddEntry(h1, "corrected yield from data", "LP")
        legend.AddEntry(h2, "corrected yield from pythia", "LP")
        legend.Draw()

        txt = TPaveText(0.75,0.875,0.85,0.775,"NDC")
        txt.SetFillColor(kWhite)
        txt.SetFillStyle(0)
        txt.SetBorderSize(0)
        txt.SetTextAlign(33);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.02)
        txt.AddText("this thesis")
        txt.Draw()
        ROOT.SetOwnership(txt,False)

        txt3 = TPaveText(0.75,0.845,0.85,0.745,"NDC"); 
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
        h1.SetYTitle("#frac{1}{2#pip_{T}} #frac{d^{2}N}{dp_{T}dy}")

        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())
        