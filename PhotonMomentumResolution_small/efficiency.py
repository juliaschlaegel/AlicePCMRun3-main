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
# from PhotonMomentumResolution import FitInvMassForPt 
#from FitInvMassForPt import PairAnalyzer
#from PlotRawYield import PlotRawYieldInvMass




gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
# gStyle.SetPadTickX(1)
# gStyle.SetPadTickY(1)
# gStyle.SetErrorX(0)
# gStyle.SetEndErrorSize(5)
#from HistoFormatting import FrameSettings, CanvasSettings, PadSettings, DrawHisto, SetTitle, SetStyleTLatex


class Efficiency:
    def __init__(self, decay): 
        self.decay = decay
        #self.utils = Utility(filename_mc, config, decay)
        
        
        
    
    # def get_raw_yield(self, filename_MC_inv):
    
    #     file_MC = TFile.Open(filename_MC_inv, "READ")
    #     if file_MC:
    #         print("file geöffnet")

    #     if self.decay == "pcm":
    #         dir_name11 = "PCMPCM"
    #         dir_name12 = "qc_qc"
    #     else:
    #         dir_name11 = "PCMDalitzEE"
    #         dir_name12 = "qc_mee0_60_minpt100_maxeta09_tpconly_lowB"

    #     raw_yield_MC = file_MC.Get(dir_name11).FindObject(dir_name12).FindObject("gausexp").FindObject("fit_0.04_0.20_GeVc2").FindObject("h1yield_param")
    #     if raw_yield_MC:
    #         print("raw yield gefunden")
    #     if not raw_yield_MC:
    #         raise ValueError("raw yield histogram not found")
    #     return raw_yield_MC
        

    # def plot_combined_hist(self):
    
    #     h1 = self.get_raw_yield(hist_name_raw)
        
    #     # acc = HistogramAcceptance(filename, dir_name21, dir_name22, dir_ev, config)
    #     h2 = self.acc_class.rebin_hist(hist_name_acc)
        

    #     max_h1 = h1.GetMaximum()
    #     max_h2 = h2.GetMaximum()
    #     max_y = 1.05 * max(max_h1, max_h2)  # 5% mehr als der höchste Balken

    #     h1.SetMaximum(max_y)

    #     # Erstelle einen Canvas
    #     c1 = ROOT.TCanvas("c1", "Combined Histogram", 800, 600)
    #     # c1.cd()
    #     pad = c1.cd()
    #     pad.SetPad(0.0, 0.01, 1, 1)
    #     pad.SetMargin(0.15, 0.1, 0.1, 0.1)
    #     pad.SetTicks(1,1)
    #     pad.SetLogy()

    #     # Einstellungen für das erste Histogramm
    #     h1.SetLineColor(ROOT.kBlue)  # Blaue Linie
    #     #h1.SetFillColor(ROOT.kBlue)  # Blaue Füllung
    #     #h1.SetFillStyle(3001)        # Leichte Transparenz
    #     h1.Draw("hist")              # Zeichne das erste Histogramm

    #     # Einstellungen für das zweite Histogramm
    #     h2.SetLineColor(ROOT.kMagenta)   # Rote Linie
    #     #h2.SetFillStyle(3004)        # Andere Art von Transparenz
    #     h2.Draw("histsame")          # Zeichne das zweite Histogramm über das erste

    #     # Legende hinzufügen
    #     legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)  # Position der Legende
    #     legend.AddEntry(h1, "raw yield", "LP")
    #     legend.AddEntry(h2, "daughters in acceptance", "LP")
    #     legend.Draw()

    #     # Achsenbeschriftungen
    #     h1.SetTitle("Combined Histogram")
    #     h1.SetXTitle("p_T")
    #     h1.SetYTitle("rebinned histograms for the efficiency")

    #     c1.Update()
    #     c1.SaveAs("combined_histogram_eff")


    # def calc_efficiency_vs_pt_pc(self, filename_mc, config, filename_MC_inv, name_plot, save): # this definition is only to work on without having all values in percent
    #     #    raise ValueError("Histogram 1 not found")
    #     if self.decay == "pcm":
    #         acc_class= HistogramAcceptance(filename_mc, config, "pcm")
    #     else:
    #         acc_class = HistogramAcceptance(filename_mc, config, "dalitz")
    
    #     h1 = self.get_raw_yield(filename_MC_inv)
       
    #     h2 = acc_class.rebin_hist("hPt_Pi0_Acc")
        

    # #get canvas and pad to plot the efficiency in 
    #     canvas = TCanvas("eff1", "eff1", 0, 0, 900, 900)
    #     pad = canvas.cd()
    #     pad.SetPad(0.0, 0.01, 1, 1)
    #     pad.SetMargin(0.15, 0.1, 0.1, 0.1)
    #     pad.SetTicks(1,1)
    #     pad.SetLogy()
    #     pad.SetLogx()

    #     # Define third histogram to plot acceptance in 
    #     h_eff = h2.Clone("h_eff")

    #     # get sum over squared weights
    #     # h1.Sumw2()
    #     # h2.Sumw2()
    
    #     # Calculate the efficiency by dividing the matched events by the events with daughters in acceptance
    #     # this time in percent 
    #     h_eff.Divide(h1, h2, 100, 1 )
    #     # Define marker settings
    #     h_eff.SetFillColor(kCyan+1)
    #     h_eff.SetMarkerStyle(kFullCross)
    #     h_eff.SetMarkerColor(kBlue)
    #     h_eff.SetLineColor(kBlue)

    #     h_eff.Draw("Esame")

        
       
    #     # Define y-axis settings
    #     y = h_eff.GetYaxis() #.SetRangeUser(min(bin_var), max(bin_var))
    #     y.SetTitle("Efficiency")
    #     y.SetTitleSize(0.048)
    #     y.SetTitleFont(42)
    #     y.SetLabelFont(42)
    #     y.SetLabelSize(0.035)
        
    #     #y.SetRangeUser(0, 0.4)
    #     #ratio.GetYaxis().SetRangeUser(-0.5, 2.5)

    #     #Define x-axis settings
    #     x = h_eff.GetXaxis()
    #     x.SetTitle("p_{T} [GeV/c]")
    #     x.SetTitleSize(0.048)
    #     x.SetTitleFont(42)
    #     x.SetTitleOffset(0.9)
    #     x.SetLabelFont(42)
    #     x.SetLabelSize(0.035)

    #     # Define legend settings
    #     leg = TLegend(0.6, 0.7, 1.0, 0.7)
    #     leg.SetBorderSize(0)
    #     leg.SetFillColor(kWhite)
    #     leg.SetFillStyle(0)
    #     leg.SetTextSize(0.03)
    #     if self.decay == "pcm":
    #         leg.AddEntry(h_eff, "\pi^{0} \\rightarrow \gamma \gamma", "LP")
    #     else:
    #         leg.AddEntry(h_eff, "\pi^{0} \\rightarrow e^{+} e^{-} \gamma", "LP")
    #     leg.Draw("")
    #     ROOT.SetOwnership(leg,False)
        
        
    #     # Add text 
    #     txt = TPaveText(0.85,0.85,0.85,0.85,"NDC")
    #     txt.SetFillColor(kWhite)
    #     txt.SetFillStyle(0)
    #     txt.SetBorderSize(0)
    #     txt.SetTextAlign(33);#middle,left
    #     txt.SetTextFont(42);#helvetica
    #     txt.SetTextSize(0.02)
    #     txt.AddText("this thesis")
    #     txt.Draw()
    #     ROOT.SetOwnership(txt,False)

    #     txt3 = TPaveText(0.85,0.82,0.85,0.82,"NDC"); # txt3 = TPaveText(0.4,0.75,0.4,0.90,"NDC");
    #     txt3.SetFillColor(kWhite)
    #     txt3.SetFillStyle(0)
    #     txt3.SetBorderSize(0)
    #     txt3.SetTextAlign(33);#middle,left
    #     txt3.SetTextFont(42);#helvetica
    #     txt3.SetTextSize(0.02)
    #     txt3.AddText("pp, #sqrt{s} = 13.6TeV")
    #     txt3.Draw()
    #     ROOT.SetOwnership(txt3,False)

    #     #save histogram
    #     canvas.Modified()
    #     canvas.Update()
    #     ROOT.SetOwnership(canvas, False)
    #     if save =="true":
    #         canvas.SaveAs(name_plot.Data())
    #     # else:
    #     #     if self.decay == "pcm":
    #     #         canvas.SaveAs("Efficiency_pcm")
    #     #     else:
    #     #         canvas.SaveAs("Efficiency_dalitz")
    #     print("eff hist type:", type(h_eff))
    #     return h_eff
    
    
    # def efficiency_vs_pt_default(self, filename_mc, config, raw_yield_mc, name_plot, save): 
    #     if self.decay == "pcm":
    #         acc_class= HistogramAcceptance(filename_mc, config, "pcm")
    #     else:
    #         acc_class = HistogramAcceptance(filename_mc, config, "dalitz")
    #     h1 = raw_yield_mc
    #     print("####Raw yield: ", type(h1), "####")
    #     h2 = acc_class.rebin_hist("hPt_Pi0_Acc")

    #     #get canvas and pad to plot the efficiency in 
    #     canvas = TCanvas("eff1", "eff1", 0, 0, 900, 900)
    #     pad = canvas.cd()
    #     pad.SetPad(0.0, 0.01, 1, 1)
    #     pad.SetMargin(0.15, 0.1, 0.1, 0.1)
    #     pad.SetTicks(1,1)
    #     pad.SetLogy()
    #     pad.SetLogx()

    #     # Define third histogram to plot acceptance in 
    #     h_eff = h2.Clone("h_eff")
            
    #     # Calculate the efficiency by dividing the matched events by the events with daughters in acceptance
    #     h_eff.Divide(h1, h2, 1, 1 )
    #     # Define marker settings
    #     h_eff.SetFillColor(kCyan+1)
    #     h_eff.SetMarkerStyle(kFullCross)
    #     h_eff.SetMarkerColor(kBlue)
    #     h_eff.SetLineColor(kBlue)

    #     h_eff.Draw("Esame")

    #     # Define y-axis settings
    #     y = h_eff.GetYaxis() #.SetRangeUser(min(bin_var), max(bin_var))
    #     y.SetTitle("Efficiency")
    #     y.SetTitleSize(0.048)
    #     y.SetTitleFont(42)
    #     y.SetLabelFont(42)
    #     y.SetLabelSize(0.035)
            

    #     #Define x-axis settings
    #     x = h_eff.GetXaxis()
    #     x.SetTitle("p_{T} [GeV/c]")
    #     x.SetTitleSize(0.048)
    #     x.SetTitleFont(42)
    #     x.SetTitleOffset(0.9)
    #     x.SetLabelFont(42)
    #     x.SetLabelSize(0.035)

    #     # Define legend settings
    #     leg = TLegend(0.2, 0.8, 0.4, 0.75)
    #     leg.SetBorderSize(0)
    #     leg.SetFillColor(kWhite)
    #     leg.SetFillStyle(0)
    #     leg.SetTextSize(0.03)
    #     if self.decay == "pcm":
    #         leg.AddEntry(h_eff, "\pi^{0} \\rightarrow \gamma \gamma", "LP")
    #     else:
    #         leg.AddEntry(h_eff, "\pi^{0} \\rightarrow e^{+} e^{-} \gamma", "LP")
    #     leg.Draw("")
    #     ROOT.SetOwnership(leg,False)
            
            
    #     # Add text 
    #     txt = TPaveText(0.3,0.85,0.3,0.85,"NDC")
    #     txt.SetFillColor(kWhite)
    #     txt.SetFillStyle(0)
    #     txt.SetBorderSize(0)
    #     txt.SetTextAlign(33);#middle,left
    #     txt.SetTextFont(42);#helvetica
    #     txt.SetTextSize(0.02)
    #     txt.AddText("this thesis")
    #     txt.Draw()
    #     ROOT.SetOwnership(txt,False)

    #     txt3 = TPaveText(0.35,0.82,0.35,0.82,"NDC"); # txt3 = TPaveText(0.4,0.75,0.4,0.90,"NDC");
    #     txt3.SetFillColor(kWhite)
    #     txt3.SetFillStyle(0)
    #     txt3.SetBorderSize(0)
    #     txt3.SetTextAlign(33);#middle,left
    #     txt3.SetTextFont(42);#helvetica
    #     txt3.SetTextSize(0.02)
    #     txt3.AddText("pp, #sqrt{s} = 13.6TeV")
    #     txt3.Draw()
    #     ROOT.SetOwnership(txt3,False)

    #     #save histogram
    #     canvas.Modified()
    #     canvas.Update()
    #     ROOT.SetOwnership(canvas, False)
    #     if save =="true":
    #         canvas.SaveAs(name_plot.Data())
    #         print("eff hist type:", type(h_eff))
    #         print("Return h_eff")
    #         return h_eff
    #     # else:
    #     #     if self.decay == "pcm":
    #     #         canvas.SaveAs("Efficiency_pcm")
    #     #         #return h_eff
    #     #     else:
    #     #         canvas.SaveAs("Efficiency_dalitz")
    #         print("eff hist type:", type(h_eff))
    #         print("Return h_eff")
    #         return(h_eff)
                    
    #     # print("eff hist type:", type(h_eff))
    #     # print("Return h_eff")
        # return h_eff
    def calc_efficiency_vs_pt(self, name_plot, save, raw_mc, utils): #, acc_class):
        
        #Get wanted histograms:
        h1 = raw_mc
        print("First Histogram rebinned")
        #utils = Utility(filename_mc, config, decay)
        h2 = utils.rebin_hist("hPt_Pi0_Acc")
        utils.scale_by_Nev(h2)
        #h2 = acc_class.rebin_hist("hPt_Pi0")
        print("second histogram rebinned")
        if h1 is None:
            print(f"Failed to load histogram: hPt_Pi0")
        # print("h1 type:", type(h1))

    #get canvas and pad to plot the acceptance in 
        can = TCanvas("eff1", "eff1", 0, 0, 900, 900)
        pad = can.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()
        pad.SetLogx()
        

        h_eff = h1.Clone("h_E")   
        
        
        # calculate acceptance by dividing the events with daughters in acceptance by all events
        # this time not in percent
        #h_eff.Divide(h1, h2, 1, 1 )#, option="B") ###Fragen: Welcher Fehler? und in prozent?
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
        #x.SetRangeUser(0,14)

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
        # else:
        #     if self.decay =="pcm":
        #         canvas.SaveAs("Acceptance_pcm")
        #     else:
        #         canvas.SaveAs("Acceptance_dalitz")
        print("eff hist type:", type(h_eff))
        return h_eff
    
        

# using the class
if __name__ == "__main__":
    filename_mc = "/Users/juliaschlagel/Analysis240411/Analysis/HLtrains/MC/AnalysisResults_HL197593_LHC24b1.root"
    #dir_name = "pi0eta-to-gammagamma"
    filename_MC = "/Users/juliaschlagel/Analysis240411/Analysis/HLtrains/MC/AnalysisResults_HL197593_LHC24b1.root"
    

    
    
    filename1_MC_dalitz = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC24b1_invariant_mass_dalitz/this_analysis_LHC24b1_pi0_mc_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"

    config = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC22o_small.yml"

    eff = Efficiency("pi0", "dalitz")
    filename_MC_inv = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC24b1_invariant_mass_dalitz_data_MC/this_analysis_LHC24b1_pi0_mc_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
    #eff.get_raw_yield(filename_MC)
    #eff.calc_efficiency_vs_pt(filename_MC_inv, "Efficiency")

    # eff.get_raw_yield(hist_name_raw)
    # eff.plot_efficiency_vs_pt(filename, dir_name21, dir_name22, dir_ev, config, hist_name_acc)
    # eff.hist_op()
    # eff.test()
    # eff.use_histogram()

#PCMDalitzEE
# hist1 = filename2.Get("PCMPCM").FindObject("qc_qc").FindObject("gausexp").FindObject("fit_0.04_0.20_GeVc2").FindObject("h1yield_param")
# if not hist1:
