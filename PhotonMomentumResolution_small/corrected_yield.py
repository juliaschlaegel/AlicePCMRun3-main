import ROOT
import ctypes
from ctypes import *
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullCross, kFullDiamond, kOpenSquare, kOpenTriangleUp, kOpenCircle, kFullCrossX
import numpy as np
# import matplotlib.pyplot as plt
# from hist import Hist
from acceptance import HistogramAcceptance
from efficiency import Efficiency
from utility import Utility

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
# gStyle.SetErrorX(0)
# gStyle.SetEndErrorSize(5)

class CorrectedYield:
    # Defining a class and load the scaling factors: The kind of meson wanted, the Number of Events,
    # the rapidity delta and the scaling factor
    def __init__(self, filename_mc, decay): #, acc_class, eff_class):
        # self.acc_class = acc_class
        # self.eff_class = eff_class
        # self.BR = BR
        # self.deltarapid = deltarapid
        self.filename_mc = TFile.Open(filename_mc, "READ")
        #self.filename1_MC = TFile.Open(filename1_MC, "READ")
        #self.filename_inv_mass_data_dalitz = TFile.Open(filename_inv_mass_data_dalitz)
        #self.filename1_data = TFile.Open(filename1_data, "READ")
        self.decay = decay

    
    # def get_raw_yield_data(self, filename_inv_mass_data):
    #     #if self.filename1_data:
    #     filename = TFile.Open(filename_inv_mass_data, "READ")
    #     if filename:
    #         print("file wurde geladen")
    #     if self.decay == "pcm":
    #         dir_name11 = "PCMPCM"
    #         dir_name12 = "qc_qc"
    #     else:
    #         dir_name11 = "PCMDalitzEE"
    #         dir_name12 = "qc_mee0_60_minpt100_maxeta09_tpconly_lowB"
    #     raw_1 = filename.Get(dir_name11).FindObject(dir_name12)
    #     if raw_1:
    #       print("raw yield data gefunden")
    #     raw_yield_data = raw_1.FindObject("gausexplinear").FindObject("fit_0.04_0.20_GeVc2").FindObject("h1yield_param")
    #     if raw_yield_data:
    #         print("alles geladen")
    #     return raw_yield_data
    

    
    #def calculate_corrected_yield(self, name_plot, filename_mc, config, raw_hist_data, raw_hist_mc, save): # filename_inv_mass_mc, filename_inv_mass_data, save): #, hist_name_acc, hist_name_all, hist_name_raw):
    def calculate_corrected_yield(self, name_plot, filename_mc, config, raw_hist_data, raw_hist_mc, save):
        if self.decay == "pcm":
            acc_class = HistogramAcceptance(filename_mc, config, "pcm")
            eff_class = Efficiency("pcm")
            utils = Utility(filename_mc, config, "pcm")
        else:
            acc_class = HistogramAcceptance(filename_mc, config, "dalitz")
            eff_class = Efficiency("dalitz")
            utils = Utility(filename_mc, config, "dalitz")

        #Zuerst alle nötigen Histogramme und Werte laden (Efficiency, Acceptance, alle caling faktoren)
        deltaRapid = 1.8

        #Acceptance 
        acc = acc_class.calc_acceptance_vs_pt("Acceptance", "false")
        print("Typ der Akzeptanz: ", type(acc))
        
        # transveres momentum 
        p_T = utils.get_bin_var()
        #print("p_T: ", p_T, "typ: ", type(p_T), "length: ", len(p_T))

        #Efficiency
        eff = eff_class.calc_efficiency_vs_pt("Efficiency", "false", raw_hist_mc, acc_class)
        #eff = eff_hist
        print("Typ der Effizienz: ", type(eff))
        
        #raw yield
        #raw_data = self.get_raw_yield_data(filename_inv_mass_data)
        raw_data = raw_hist_data
        print("Raw data type: ", type(raw_data))
        
        #produce new empty canvas and pad to plot the corrected yield in 
        canvas = TCanvas("corr", "corr", 0, 0, 900, 900)
        pad = canvas.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)

        
        h_corrected_yield = raw_data.Clone("corrected_yield_hist") 
        print("clone von raw data type: ", type(h_corrected_yield))

        #h_corrected_yield = ROOT.TH1D("h_corrected_yield", "h_corrected_yield", (len(p_T)-1), 0, (max(p_T)))
        #h_corrected_yield.Add(raw_data)

    
        # get sum over squared weights
        h_corrected_yield.Sumw2()

        # x = h_corrected_yield.GetBinContent(5)
        #print("bin content before scaling:", x)
        
        # for i in range (1, len(p_T) - 1):
        #      bin_content = h_corrected_yield.GetBinContent(i)
        #      print("bin_content_raw: ", bin_content)
        #      bin_center = h_corrected_yield.GetBinCenter(i)
        #      bin_width = h_corrected_yield.GetBinWidth(i)
        #      bin_error = h_corrected_yield.GetBinError(i)
             
        #      new_bin_content = bin_content / (bin_width * bin_center) 
             
        #      h_corrected_yield.SetBinContent(i, new_bin_content)

        #      new_bin_error = bin_error / (bin_width * bin_center) 
        #      h_corrected_yield.SetBinError(i, new_bin_error)


        h_corrected_yield.Scale(1 / (2 * np.pi))
        h_corrected_yield.Scale(1 / deltaRapid)
        # y= h_corrected_yield.GetBinContent(5)
        # print("bin content after scaling: ", y)


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
        #Draw points in the histogramm
        h_corrected_yield.Draw("Esame")
       
        # Define y-axis settings
        y = h_corrected_yield.GetYaxis() 
        # y.SetLog()
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
        #if save_name == "Corrected_yield_dalitz":
        #    leg.AddEntry(corrected_yield_hist, "\pi_{0} \\rightarrow e^{+} e^{-} \gamma", "LP")
        #else: 
            #leg.AddEntry(corrected_yield_hist, "\pi_{0} \\rightarrow e^{+} e^{-}", "LP")
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
        # else:
        #     if self.decay == "pcm":
        #         canvas.SaveAs("Corrected_yield_pcm")
        #     else:
        #         canvas.SaveAs("Corrected_yield_dalitz")
        print("type corrected yield: ", type(h_corrected_yield))
        # h3 = TH1D("pT3", "pT3", (len(p_T)-1),min(p_T), max(p_T)) 
        # h3.Add(h_corrected_yield)

        return h_corrected_yield
    
    
    def calc_corr_vs_pt(self, name_plot, save, raw_data, decay, raw_mc, filename_mc, config, eff_hist, acc_hist): # eff_hist, acc_hist, raw_hist_data, filename_mc, config):
         #Get wanted histograms:
        deltay = 1.8
        raw = raw_data
        #     raw_h1d = ROOT.TH1D(raw.GetName() + "_TH1D", raw.GetTitle(), raw.GetNbinsX(), raw.GetXaxis().GetXmin(), raw.GetXaxis().GetXmax())

        # # Übertrage die Daten von h1f zu h1d
        #     for i in range(1, raw.GetNbinsX() + 1):
        #         raw_h1d.SetBinContent(i, raw.GetBinContent(i))
        #         raw_h1d.SetBinError(i, raw.GetBinError(i))
        acc_class = HistogramAcceptance(filename_mc, config, "dalitz")
        
        # p_T = acc_class.get_bin_var()
        # eff = eff_hist
        # for i in range (1, len(p_T)):
        #      bin_content = eff.GetBinContent(i)
        #      #print("bin_content_raw: ", bin_content)
        #      bin_center = eff.GetBinCenter(i)
        #      bin_width = eff.GetBinWidth(i)
        #      bin_error = eff.GetBinError(i)
        #      print("Bin center eff", bin_center, "Bin content eff", bin_content)
        # acc = acc_hist
        # for i in range (1, len(p_T)):
        #      bin_content = acc.GetBinContent(i)
        #      #print("bin_content_: ", bin_content)
        #      bin_center = acc.GetBinCenter(i)
        #      bin_width = acc.GetBinWidth(i)
        #      bin_error = acc.GetBinError(i)
        #      print("Bin center acc", bin_center, "Bin content acc", bin_content)
        #h_corrected_yield = raw.Clone("corr_yield_1")
        #print("h_corr_yield: ", type(h_corrected_yield))
            
            
        #get canvas and pad to plot the acceptance in 
        canvas_loc = TCanvas("corr_yield_loc", "corr_yield_loc", 0, 0, 900, 900)
        pad = canvas_loc.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        #pad.SetLogy()
        #pad.SetLogx()
        print("typ des canvas: ", type(canvas_loc))
        
        print("Typ von raw: ", type(raw))
        #h_raw = raw_h1d.Clone("h_corr")  
        #h_eff = eff.Clone("h_eff") 
        h_corrected_yield = raw.Clone("h_corr")
        print("typ von corr yield vor eff und acc", type(h_corrected_yield))
        
        h_corrected_yield.Sumw2()

        h_corrected_yield.Scale(1/(2*np.pi))
        h_corrected_yield.Scale(1/deltay)
        #eff_class = Efficiency(decay)
        #acc_class = HistogramAcceptance(filename_mc, config, "dalitz")
        
        #p_T = acc_class.get_bin_var()
        # print(p_T)
        
        #eff = eff_class.calc_efficiency_vs_pt("Efficiency", "false", raw_mc, acc_class)
        # for i in range (1, len(p_T)):
        #      bin_content = eff.GetBinContent(i)
        #      print("bin_content_raw: ", bin_content)
        #      bin_center = eff.GetBinCenter(i)
        #      bin_width = eff.GetBinWidth(i)
        #      bin_error = eff.GetBinError(i)
        #      print("Bin center eff", bin_center, "Bin content eff", bin_content)
        #h_corrected_yield.Divide(eff)
        # for i in range (1, len(p_T)):
        #      bin_content = h_corrected_yield.GetBinContent(i)
        #      print("bin_content_raw: ", bin_content)
        #      bin_center = h_corrected_yield.GetBinCenter(i)
        #      #bin_width = eff.GetBinWidth(i)
        #      #bin_error = eff.GetBinError(i)
        #      print("Bin center eff", bin_center, "Bin content eff", bin_content)
        
        #acc = acc_class.calc_acceptance_vs_pt("Acceptance", "false")
        #acc_loc = acc.Clone("Acc_loc")
        # for i in range (1, len(p_T) - 1):
        #      bin_content = acc.GetBinContent(i)
        #      print("bin_content_raw: ", bin_content)
        #      bin_center = acc.GetBinCenter(i)
        #      bin_width = acc.GetBinWidth(i)
        #      bin_error = acc.GetBinError(i)
        #      print("Bin center acc", bin_center, "Bin content acc", bin_content)
        #h_corrected_yield.Divide(acc)
        # for i in range (1, len(p_T) - 1):
        #      bin_content = h_corrected_yield.GetBinContent(i)
        #      print("bin_content_raw: ", bin_content)
        #      bin_center = h_corrected_yield.GetBinCenter(i)
        #      bin_width = h_corrected_yield.GetBinWidth(i)
        #      bin_error = h_corrected_yield.GetBinError(i)
        #      print("Bin center corrected yield", bin_center, "Bin content corrected yield", bin_content)
        # for i in range (1, len(p_T) - 1):
        #      bin_content = h_corrected_yield.GetBinContent(i)
        #      print("bin_content_raw: ", bin_content)
        #      bin_center = h_corrected_yield.GetBinCenter(i)
        #      bin_width = h_corrected_yield.GetBinWidth(i)
        #      bin_error = h_corrected_yield.GetBinError(i)
             
        #      new_bin_content = bin_content / (bin_width * bin_center) 
             
        #      h_corrected_yield.SetBinContent(i, new_bin_content)

        #      new_bin_error = bin_error / (bin_width * bin_center) 
        #      h_corrected_yield.SetBinError(i, new_bin_error)
        
        print("Typ von korrigiertem Histogram: ", type(h_corrected_yield))
    
        # Define marker settings
        h_corrected_yield.SetFillColor(kCyan+1)
        h_corrected_yield.SetMarkerStyle(kFullCross)
        h_corrected_yield.SetMarkerColor(kBlue)
        h_corrected_yield.SetMarkerSize(0.85)
        h_corrected_yield.SetLineColor(kBlue)

        #Draw points in the histogramm
        h_corrected_yield.Draw("ESame")
        #eff.Draw("ESame")
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
        #x.SetRangeUser(0,14)

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
        # else:
        #     if self.decay =="pcm":
        #         canvas.SaveAs("Acceptance_pcm")
        #     else:
        #         canvas.SaveAs("Acceptance_dalitz")
        print("corr hist type:", type(h_corrected_yield))
        return h_corrected_yield
    

    def calc_corr_yield_vs_pt(self, name_plot, save, raw_data, eff_hist, acc_hist, acc_class):
        
        #Get wanted histograms:
        deltarap = 1.8 
        #acc_class = HistogramAcceptance(filename_mc, config, decay)
        #p_T = utils.get_bin_var()
        #print("Binning edges: ", p_T)
        #NeV = utils.get_Nev()
        #print("NUMBER OF EVENTS: ", NeV)
        h1 = raw_data
        print("First Histogram loaded, raw_data", type(h1))
        eff = eff_hist
        print("second histogram loaded, efficiency", type(eff))
        acc = acc_hist
        print("third histogram loaded, acc", type(acc))
        if h1 is None:
            print(f"Failed to load histogram: hPt_Pi0")
        # print("h1 type:", type(h1))

    #get canvas and pad to plot the acceptance in 
        can = TCanvas("correc", "correc", 0, 0, 900, 900)
        pad = can.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()
        pad.SetLogx()
        

        h_corr = h1.Clone("h_c")   
        
        
        # calculate acceptance by dividing the events with daughters in acceptance by all events
        # this time not in percent
        #h_eff.Divide(h1, h2, 1, 1 )#, option="B") ###Fragen: Welcher Fehler? und in prozent?
        h_corr.Divide(eff)
        h_corr.Divide(acc)
        # for i in range (1, len(p_T) - 1):
        #     bin_content = h_corr.GetBinContent(i)
        #     print("bin_content_raw: ", bin_content)
        #     bin_center = h_corr.GetBinCenter(i)
        #     bin_width = h_corr.GetBinWidth(i)
        #     bin_error = h_corr.GetBinError(i)
             
        #     new_bin_content = bin_content / (bin_width * bin_center) 
             
        #     h_corr.SetBinContent(i, new_bin_content)

        #     new_bin_error = bin_error / (bin_width * bin_center) 
        #     h_corr.SetBinError(i, new_bin_error)

        h_corr.Scale(1/deltarap)
        h_corr.Scale(1/(2* np.pi))
        
        #h_corr.Scale(1/NeV)
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
        #x.SetRangeUser(0,14)

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
        print("corrected yield hist type:", type(h_corr))
        return h_corr
    
    
    
    
    def plot_pythia(self, utils, name_plot, save): #This is the corrected yield we would expect from MC
        if self.decay =="pcm":
            dir_name22 = "PCMPCM"    
            #hist = self.filename.Get(dir_name21).FindObject(dir_name22).FindObject(hist_name_all)
        else:
            dir_name22 = "PCMDalitzEE"
        hist = self.filename_mc.Get("pi0eta-to-gammagamma-mc/Generated").FindObject(dir_name22).FindObject("hPt_Pi0")
        #hist = raw
        deltay = 1.8
        if hist:
            print("histogram importiert")
        # hist.Scale(1 / (2*np.pi))
        Nev = utils.get_Nev()
        print("Number of Events: ", Nev)
       
        h1 = hist.Clone("pythia")
        
        for bin in range (1, h1.GetXaxis().GetNbins()):
             bin_content = h1.GetBinContent(bin)
             bin_center = h1.GetBinCenter(bin)
             bin_width = h1.GetBinWidth(bin)
             # print("bin width: ", bin_width)
             bin_error = h1.GetBinError(bin)

             new_bin_content = bin_content / (bin_center *bin_width) #* Nev* deltay*2*np.pi)
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
        # y.SetLog()
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
        # else:
        #     if self.decay =="pcm":
        #         canvas.SaveAs("pythia_pcm")
        #     else:
        #         canvas.SaveAs("pythia_dalitz")
        #canvas.SaveAs("pythia_test")
        #canvas.SaveAs("pythia.root")
        return h1
    
    
    def compare_pythia_corr_yield(self, raw_data, eff_hist, acc_hist, utils, name_plot, save):
        # h1 = self.calculate_corrected_yield_dalitz(filename1_data_dalitz, hist_name_acc, hist_name_all, hist_name_raw, acc_class_dalitz, eff_dalitz, dir_name11_dalitz, dir_name12_dalitz, dir_name13_data, dir_name14)
        #h1 = self.calculate_corrected_yield_dalitz( filename1_data_dalitz, hist_name_acc, hist_name_all, hist_name_raw, acc_class_dalitz, eff_dalitz, dir_name11_dalitz, dir_name12_dalitz, dir_name13_data, dir_name14)
        if self.decay == "pcm":
            dir_name22 = "PCMPCM"
        else:
            dir_name22 = "PCMDalitzEE"
        h1 = self.calc_corr_yield_vs_pt("corrected_yield", "false", raw_data, eff_hist, acc_hist, utils)
        h2 = self.plot_pythia(utils, "pythia", "false")

        max_h1 = h1.GetMaximum()
        max_h2 = h2.GetMaximum()
        max_y = 1.05 * max(max_h1, max_h2)  # 5% mehr als der höchste Balken

        h1.SetMaximum(max_y)

        # Erstelle einen Canvas
        c1 = ROOT.TCanvas("py_vs_cor", "py_vs_cor", 800, 600)
        pad = c1.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()

        # Einstellungen für das erste Histogramm
        h1.SetLineColor(ROOT.kBlue)
        h1.SetMarkerColor(ROOT.kBlue)
        h1.SetMarkerStyle(21)  # Kreisförmige Punkte
        h1.Draw("E1P")          # Zeichne das erste Histogramm als Punkte mit Fehlerbalken

        # Einstellungen für das zweite Histogramm
        h2.SetLineColor(ROOT.kMagenta)
        h2.SetMarkerColor(ROOT.kMagenta)
        h2.SetMarkerStyle(22)  # Dreieckige Punkte
        h2.Draw("E1P same")     # Zeichne das zweite Histogramm über das erste

        # Legende hinzufügen
        # legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
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

        txt3 = TPaveText(0.75,0.845,0.85,0.745,"NDC"); # txt3 = TPaveText(0.4,0.75,0.4,0.90,"NDC");
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
        h1.SetYTitle("#frac{1}{2#pip_{T}} #frac{d^{2}N}{dp_{T}dy}")

        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())
        # else:
        #     if self.decay =="pcm":
        #         c1.SaveAs("comp_pythia_corr_pcm")
        #     else:
        #         c1.SaveAs("comp_pythia_corr_dalitz")

    # def compare_corr_dalitz_pcm(self,corr_dalitz, corr_pcm):
    #     h1 = corr_dalitz.calculate_corrected_yield_dalitz()

    #     h2 = corr_pcm.calculate_corrected_yield("Corrected_yield_pcm")
    #     print("h1 Bin content:")
    #     h1.Print()
    #     p_T = self.acc_class.get_bin_var()
    #     for i in range(1, len(p_T)-1):
    #         x = h1.GetBinContent(i)
    #         y= h1.GetBinCenter(i)
    #         print("h1 Bin content: ", y, x )
    #         x2 = h2.GetBinContent(i)
    #         y2 = h2.GetBinCenter(i)
    #         print("h2 Bin content: ", y2, x2)

        
    #     print("h2 Bin content")
    #     h2.Print()
    #     ratio = h1.Clone("Ratio_corr_yield")
    #     ratio.Divide(h2)


    #     max_h1 = h1.GetMaximum()
    #     min_h1 = h1.GetMinimum()
    #     max_h2 = h2.GetMaximum()
    #     min_h2 = h2.GetMinimum()
    #     max_y = 1.05 * max(max_h1, max_h2)  # 5% mehr als der höchste Balken
    #     min_y = 1.05 * min(min_h1, min_h2)

    #     h1.SetMaximum(max_y)
    #     # h2.SetMinimum(min_y)

    #     # Erstelle einen Canvas
    #     c1 = ROOT.TCanvas("py_vs_cor", "py_vs_cor", 800, 600)
    #     c1.Divide(1, 2)

    #     p1 = c1.cd(1)
    #     p1.SetPad(0.0, 0.0, 1, 1)
    #     p1.SetMargin(0.15,0.02,0.22,0.0)
    #     p1.SetTicks(1,1)
    #     p1.SetLogy()
    #     # p2 = c1.cd(2);
    #     # p2.SetPad(0,0,1,0.3);
    #     # p2.SetMargin(0.15,0.02,0.22,0.0);
    #     # p2.SetTicks(1,1);
    #     # p2.SetLogy();
    #     # Margin: 0.15, 0.1, 0.1, 0.1
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
    #     legend1 = ROOT.TLegend(0.6, 0.9, 0.95, 0.75)
    #     legend1.SetBorderSize(0)
    #     legend1.SetFillColor(kWhite)
    #     legend1.SetFillStyle(0)
    #     legend1.SetTextSize(0.03)
    #     legend1.AddEntry(h1, "corrected yield for #pi^{0} -> e^{+} e^{-} #gamma", "LP")
    #     legend1.AddEntry(h2, "corrected yield for #pi^{0} -> e^{+} e^{-}", "LP")
    #     legend1.Draw()

    #     txt = TPaveText(0.9,0.975,0.95,0.95,"NDC")
    #     txt.SetFillColor(kWhite)
    #     txt.SetFillStyle(0)
    #     txt.SetBorderSize(0)
    #     txt.SetTextAlign(33);#middle,left
    #     txt.SetTextFont(42);#helvetica
    #     txt.SetTextSize(0.02)
    #     txt.AddText("this thesis")
    #     txt.Draw()
    #     ROOT.SetOwnership(txt,False)

    #     txt3 = TPaveText(0.9,0.945,0.95,0.92,"NDC"); # txt3 = TPaveText(0.4,0.75,0.4,0.90,"NDC");
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
    #     h1.SetTitle("Corrected yield")
    #     h1.SetXTitle("p_{T} (GeV/c)")
    #     h1.SetYTitle("#frac{1}{2#pip_{T}} #frac{d^{2}N}{dp_{T}dy}")

    #     #Go to second canvas to plot ratio in
        
    #     p2 = c1.cd(2)
    #     p2.SetPad(0,0,1,0.3)
    #     p2.SetMargin(0.15,0.02,0.22,0.0)
    #     p2.SetTicks(1,1)
    #     # p2.SetLogy()

    #     ratio.SetLineColor(ROOT.kBlue)
    #     ratio.SetMarkerColor(ROOT.kBlue)
    #     ratio.SetMarkerStyle(21)  # Kreisförmige Punkte
    #     ratio.Draw("E1P")   
         
    #     legend2 = ROOT.TLegend(0.2,0.9,0.3,0.8)
    #     legend2.SetBorderSize(0)
    #     legend2.SetFillColor(kWhite)
    #     legend2.SetFillStyle(0)
    #     legend1.SetTextSize(0.03)
    #     legend2.AddEntry(ratio, "ratio dalitz/pcm", "LP")
    #     legend2.Draw()  

    #     ratio.SetTitle("Ratio Dalitz/pcm")
    #     ratio.SetXTitle("p_{T} (GeV/c)")
    #     ratio.GetXaxis().SetLabelSize(0.10)
    #     ratio.GetXaxis().SetTitleSize(0.10)
    #     ratio.GetXaxis().SetTitleOffset(1.0)
    #     ratio.GetYaxis().SetLabelSize(0.10)
        
    #     ratio.GetYaxis().SetRangeUser(0.8, 1.2)
    #     ratio.GetYaxis().SetTitleSize(0.10)
    #     ratio.GetYaxis().SetTitleOffset(.3)
    #     ratio.SetYTitle("Ratio")

    #     c1.Update()
    #     c1.SaveAs("compare_corr_dalitz_vs_pcm_test")





            

        
        
        

# Using the class
if __name__ == "__main__":
    filename_mc = "/Users/juliaschlagel/Analysis240411/Analysis/HLtrains/MC/AnalysisResults_HL197593_LHC24b1.root"
    # dir_name21 = "pi0eta-to-gammagamma-mc/Generated"
    # dir_name22_dalitz = "PCMDalitzEE"              
    # dir_name22_pcm ="PCMPCM"    
    # dir_name_py_1 = "associate-mc-info/Generated"
    
    
   
    filename_inv_mass_data_dalitz = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC22o_pass6_small_invariant_mass_dalitz_data_MC/this_analysis_LHC22o_pass6_small_pi0_data_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
    filename_inv_mass_data = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC22o_pass6_small_invariant_mass_pcm_data_MC/this_analysis_LHC22o_pass6_small_pi0_data_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
    filename_inv_mass_mc_dalitz = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC24b1_invariant_mass_dalitz_data_MC/this_analysis_LHC24b1_pi0_mc_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
    filename_inv_mass_mc = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC24b1_invariant_mass_pcm_data_MC/this_analysis_LHC24b1_pi0_mc_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"

    # dir_name11_dalitz = "PCMDalitzEE"
    # dir_name11_pcm = "PCMPCM"
    # dir_name12_dalitz = "qc_mee0_60_minpt100_maxeta09_tpconly_lowB"
    # dir_name12_pcm = "qc_qc"
    # dir_name13_data = "gausexplinear"
    # dir_name13_mc = "gausexp"
    # dir_name14 = "fit_0.04_0.20_GeVc2"
    
    # hist_name_raw = "h1yield_param"
    # hist_name_acc = "hPt_Pi0_Acc"
    # hist_name_all = "hPt_Pi0"
    # dir_ev = "pi0eta-to-gammagamma-mc/Event"
    config = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC22o_small.yml"
    # acc_class_dalitz = HistogramAcceptance(filename, dir_name21, dir_name22_dalitz, dir_ev, config)
    # # acc_class_pcm = HistogramAcceptance(filename, dir_name21, dir_name22_pcm, dir_ev, config)
    # eff_dalitz = Efficiency(filename1_MC_dalitz, filename, dir_name11_dalitz, dir_name12_dalitz, dir_name13_mc, dir_name14, dir_name21, dir_name22_dalitz, acc_class_dalitz)
    # eff_pcm = Efficiency(filename1_MC_pcm, filename, dir_name11_pcm, dir_name12_pcm, dir_name13_mc, dir_name14, dir_name21, dir_name22_pcm, acc_class_pcm)



    #corr_dalitz = CorrectedYield(filename, filename1_MC_dalitz, filename1_data_dalitz) #, acc_class_dalitz, eff_dalitz)
    # corr_pcm = CorrectedYield(filename, filename1_MC_pcm, filename1_data_dalitz) # acc_class_pcm, eff_pcm)
    # corr_dalitz.calculate_corrected_yield_dalitz_crosscheck(hist_name_acc, hist_name_all, hist_name_raw, acc_class_dalitz, eff_dalitz)
    # corr_dalitz.compare_pythia_corr_yield(dir_name22_dalitz)
    #corr_dalitz.calculate_corrected_yield("dalitz", hist_name_acc, hist_name_all, hist_name_raw, dir_name13_data, dir_name14)
    #corr_dalitz.compare_pythia_corr_yield("dalitz")
    #corr_dalitz.get_raw_yield_data(dir_name11_dalitz, dir_name12_dalitz, dir_name13_data, dir_name14, hist_name_raw)
    # corr_dalitz.get_raw_yield_data(filename1_data_dalitz, dir_name11_dalitz, dir_name12_dalitz, dir_name13_data, dir_name14, hist_name_raw)
    # corr_pcm.get_raw_yield_data(filename1_data_pcm, dir_name11_pcm, dir_name12_pcm, dir_name13_data, dir_name14, hist_name_raw)
    #corr_dalitz.calculate_corrected_yield_dalitz(filename1_data_dalitz, hist_name_acc, hist_name_all, hist_name_raw, acc_class_dalitz, eff_dalitz, dir_name11_dalitz, dir_name12_dalitz, dir_name13_data, dir_name14)
    #corr_dalitz.calculate_corrected_yield_pcm(filename1_data_pcm, hist_name_acc, hist_name_all, hist_name_raw, acc_class_pcm, eff_pcm, dir_name11_pcm, dir_name12_pcm, dir_name13_data, dir_name14)
    



    # corr_dalitz.plot_eff_dalitz_vs_pcm(eff_dalitz, eff_pcm)
    # corr_dalitz.plot_acc_dalitz_vs_pcm(acc_class_dalitz, acc_class_pcm)
   
    