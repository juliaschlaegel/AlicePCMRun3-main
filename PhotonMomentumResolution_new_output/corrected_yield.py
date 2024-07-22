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
from comparison import Compare

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

class CorrectedYield:
    def __init__(self, filename_mc, decay, meson): 
        self.filename_mc = TFile.Open(filename_mc, "READ")
        self.decay = decay
        self.meson = meson

    
    def calc_corr_yield_new(self, utils, raw_data, eff_acc_hist, decay, scale_pT ,name_plot, save):
        print("Scaling by pT:", scale_pT)
        # the option scale_pt is just "true" or "false" which means if you want your corrected yield to be scaled with pT or not.
        #get wanted histograms:
        raw = raw_data
        if raw:
            print("raw yield loaded")
        eff_acc = eff_acc_hist
        if eff_acc:
            print("eff_acc_loaded")
        deltarap = 1.6
    

    #get canvas and pad to plot the corrected yield in 
        can = TCanvas("corr", "corr", 0, 0, 900, 900)
        pad = can.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.2, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()
        
        h_corrected = raw.Clone("h_corrected") 

        h_corrected.Divide(eff_acc)
        h_corrected.Scale(1/deltarap)
        h_corrected.Scale(1/(2* np.pi))
        
        if scale_pT == "true":
            utils.scale_by_pt(h_corrected)

        # Define marker settings
        h_corrected.SetFillColor(kCyan+1)
        if self.decay == "pcm":
            h_corrected.SetMarkerStyle(34)
        if self.decay== "dalitz":
            h_corrected.SetMarkerStyle(20)
        h_corrected.SetMarkerColor(kBlue)
        h_corrected.SetMarkerSize(1.5)
        h_corrected.SetLineColor(kBlue)

        #Draw points in the histogramm
        h_corrected.Draw("Esame")
       
        # Define y-axis settings
        y = h_corrected.GetYaxis() 
        if scale_pT == "true":
            y.SetTitle("#frac{1}{2#pi} #frac{1}{p_{T}} #frac{d^{2}N}{dp_{T}dy} [GeV^{-2} c^{3}]")
        else:
            y.SetTitle("#frac{1}{2#pi} #frac{d^{2}N}{dp_{T}dy} [GeV^{-1} c^{3}]") 
        y.SetTitleSize(0.048)
        y.SetTitleFont(42)
        y.SetTitleOffset(1.7)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)

        #Define x-axis settings
        x = h_corrected.GetXaxis()
        x.SetTitle("p_{T} [GeV/c]")
        x.SetTitleSize(0.048)
        x.SetTitleFont(42)
        x.SetTitleOffset(0.9)
        x.SetLabelFont(42)
        x.SetLabelSize(0.035)
        x.SetRangeUser(0.0, 15.0)

        # Define legend settings
        leg = TLegend(0.7, 0.75, 0.9, 0.75)
        leg.SetBorderSize(0)
        leg.SetFillColor(kWhite)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.03)
        if self.decay == "pcm":
            if self.meson == "pi0":
                leg.AddEntry(h_corrected, "\pi^{0} \\rightarrow \gamma \gamma", "LP")
            if self.meson == "eta":
                leg.AddEntry(h_corrected, "\eta \\rightarrow \gamma \gamma", "LP")
        else:
            if self.meson == "pi0":
                leg.AddEntry(h_corrected, "\pi^{0} \\rightarrow e^{+} e^{-} \gamma", "LP")
            if self.meson == "eta":
                leg.AddEntry(h_corrected, "\eta \\rightarrow e^{+} e^{-} \gamma", "LP")
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

        txt3 = TPaveText(0.875,0.82,0.875,0.82,"NDC"); 
        txt3.SetFillColor(kWhite)
        txt3.SetFillStyle(0)
        txt3.SetBorderSize(0)
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.02)
        txt3.AddText("pp, #sqrt{s} = 13.6TeV")
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

        txt2 = TPaveText(0.875,0.79,0.875,0.79,"NDC")
        txt2.SetFillColor(kWhite)
        txt2.SetFillStyle(0)
        txt2.SetBorderSize(0)
        txt2.SetTextAlign(33);#middle,left
        txt2.SetTextFont(42);#helvetica
        txt2.SetTextSize(0.02)
        txt2.AddText(decay)
        txt2.Draw()
        ROOT.SetOwnership(txt2,False)


        #save histogram
        can.Modified()
        can.Update()
        ROOT.SetOwnership(can, False)
        if save == "true":
            can.SaveAs(name_plot.Data())
        
        print("corrected yield hist type:", type(h_corrected))
        return h_corrected
    
    
    
    
    def plot_pythia(self, utils, name_plot, save): #This is the PYTHIA 8 model prediction
        hist =  utils.get_hist("hPt")
        deltay = 1.8
        if hist:
            print("histogram imported")
        
        Nev = utils.get_Nev()
        print("Number of Events: ", Nev)
       
        h1 = hist.Clone("pythia")
        
        for bin in range (1, h1.GetXaxis().GetNbins()):
            bin_content = h1.GetBinContent(bin)
            bin_width = h1.GetBinWidth(bin)
            bin_error = h1.GetBinError(bin)
            if bin_width != 0:
                new_bin_content = bin_content / (bin_width) 
                h1.SetBinContent(bin, new_bin_content)
                new_bin_error = bin_error / (bin_width)
                h1.SetBinError(bin, new_bin_error)

        

        h1.Scale(1 / Nev)
        h1.Scale(1 / deltay)
        h1.Scale(1 / (2 * np.pi))
        
        #This was just made, because there seems to be a bugg: the value when the binning changes is set to 0 so at 
        # pT = 5 the value of the histogram is 0
        # this part of the code just cuts out this one strange bin
        x_values = []
        y_values = []
        bin_change = h1.FindBin(5)
        print("++++++ bin change: ", bin_change, h1.GetBinCenter(bin_change))
        for bin in range(1, h1.GetNbinsX() + 1):
        #for bin in range(1, bin_change +1):
            bin_content = h1.GetBinContent(bin)
            x_values.append(h1.GetBinCenter(bin))
            y_values.append(bin_content)
        print("x_values: ", x_values, "y_values: ", y_values)
        del x_values[bin_change-1]
        del y_values[bin_change-1]
        print("x_values after removing: ", x_values, "y_values: ", y_values)
        graph = ROOT.TGraph(len(x_values), np.array(x_values), np.array(y_values))
        
        graph.SetLineWidth(3)
        graph.SetLineColor(8)
        
        

        #plot the pythia corrected yield
        canvas = TCanvas("pyt", "pyt", 0, 0, 900, 900)
        pad = canvas.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.2, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()
        
        graph.Draw("AC") #C is with a smooth line
        graph.GetYaxis().SetRangeUser(5e-7, 1.1)
       
        # Define y-axis settings
        y = graph.GetYaxis() 
        y.SetTitle("#frac{1}{2#pi} #frac{d^{2}N}{dp_{T}dy} [GeV^{-1} c^{3}]")
        y.SetTitleSize(0.048)
        y.SetTitleFont(42)
        y.SetTitleOffset(1.4)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)

        #Define x-axis settings
        x = graph.GetXaxis()
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
        return graph
    
    
    def compare_pythia_corr_yield(self, corrected_hist, pythia_hist, decay, name_plot, save):
        corr = corrected_hist
        pythia = pythia_hist

        max_h1 = corr.GetMaximum()
        max_y = 2 * max_h1  

        corr.SetMaximum(max_y)

        
        c1 = ROOT.TCanvas("py_vs_cor", "py_vs_cor", 800, 600)
        pad = c1.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()

        pythia.SetLineColor(8)
        pythia.SetMarkerColor(8)
        pythia.SetMarkerStyle(47)  
        pythia.Draw("AL")
        
        corr.SetLineColor(ROOT.kBlue)
        corr.SetMarkerColor(ROOT.kBlue)
        if self.decay == "pcm":
            corr.SetMarkerStyle(34)
        if self.decay == "dalitz":
            corr.SetMarkerStyle(20)  
        corr.Draw("E1P Same")          

        
        
        
        legend = ROOT.TLegend(0.575, 0.75, 0.8, 0.65)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)
        legend.AddEntry(corr, "data", "LP")
        legend.AddEntry(pythia, "PYTHIA", "LP")
        legend.Draw()

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

        txt3 = TPaveText(0.875,0.82,0.875,0.82,"NDC"); 
        txt3.SetFillColor(kWhite)
        txt3.SetFillStyle(0)
        txt3.SetBorderSize(0)
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.02)
        txt3.AddText("pp, #sqrt{s} = 13.6TeV")
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

        txt2 = TPaveText(0.875,0.79,0.875,0.79,"NDC")
        txt2.SetFillColor(kWhite)
        txt2.SetFillStyle(0)
        txt2.SetBorderSize(0)
        txt2.SetTextAlign(33);#middle,left
        txt2.SetTextFont(42);#helvetica
        txt2.SetTextSize(0.02)
        txt2.AddText(decay)
        txt2.Draw()
        ROOT.SetOwnership(txt2,False)

        
        corr.SetTitle("Corrected yield")
        corr.SetXTitle("p_{T} (GeV/c)")
        corr.SetYTitle("#frac{1}{2#pi} #frac{d^{2}N}{dp_{T}dy} [GeV^{-1} c^{3}]")
        corr.GetYaxis().SetTitleOffset(1.4)

        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())

    
    def plot_fit_func(self, corr_hist_scaled, decay, name_plot, save):
        #This is the fit with the Two-Component Model to the corrected yield

        c1 = ROOT.TCanvas("c_fit", "c_fit", 900,900)
        pad = c1.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.2, 0.1, 0.1, 0.1)
        pad.SetLogy()
        
        h_corr = corr_hist_scaled.Clone("corr_fit")
        h_corr.Scale(59.4*1e-3) #getting the differential invariant yield
        h_corr.Scale(1e12) #scaling to pico barn
        
        #Giving the start parameters
        A_e = 244 * 1e9
        M_pi = 0.135
        T_e = 0.157
        A_init = 27 * 1e9 
        T_init = 0.6  
        n = 2.96
    
        #defining the fit function
        f1TCM = TF1("TCM_fit", "[0]* (TMath::Exp(-(TMath::Power(TMath::Power(x,2)+TMath::Power([1],2), 0.5)-[1])/([2])))+ [3] * (TMath::Power((1+(TMath::Power(x,2))/(TMath::Power([4],2)*[5])), -[5]))", 0.2, 15.,6)
        f1TCM.SetParameters(A_e, M_pi, T_e, A_init, T_init, n)
        f1TCM.SetLineColor(ROOT.kBlack)
            
        f1TCM.SetParLimits(0, 0.5*A_e, 1.5*A_e)
        f1TCM.SetParLimits(1, M_pi, M_pi)
        f1TCM.SetParLimits(2, 0.5 * T_e, 1.5*T_e);
        f1TCM.SetParLimits(3,0.5*A_init, 1.5*A_init)
        f1TCM.SetParLimits(4, 0.5 * T_init, 1.5 * T_init)
        f1TCM.SetParLimits(5, 0.5 * n,  1.5*n)

        h_corr.Fit(f1TCM,"TCM","",0.2, 15);
        y = h_corr.GetYaxis() 
        y.SetTitle("E #frac{d^{3}#sigma}{dp^{3}} (pb GeV^{-2} c^{3})")
        y.SetTitleSize(0.048)
        y.SetTitleFont(42)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)
        
        legend = ROOT.TLegend(0.575, 0.75, 0.8, 0.65)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.02)
        legend.AddEntry(h_corr, "corrected yield from data", "LP")
        legend.AddEntry(f1TCM, "Fit with two component model", "LP")
        legend.Draw()

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

        txt3 = TPaveText(0.875,0.82,0.875,0.82,"NDC"); 
        txt3.SetFillColor(kWhite)
        txt3.SetFillStyle(0)
        txt3.SetBorderSize(0)
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.02)
        txt3.AddText("pp, #sqrt{s} = 13.6TeV")
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

        txt2 = TPaveText(0.875,0.79,0.875,0.79,"NDC")
        txt2.SetFillColor(kWhite)
        txt2.SetFillStyle(0)
        txt2.SetBorderSize(0)
        txt2.SetTextAlign(33);#middle,left
        txt2.SetTextFont(42);#helvetica
        txt2.SetTextSize(0.02)
        txt2.AddText(decay)
        txt2.Draw()
        ROOT.SetOwnership(txt2,False)


        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())
        
        params = f1TCM.GetParameters()
        params_err = f1TCM.GetParErrors()
        for i in range(6):
            print(f"Parameter {i}: {params[i]}", "±", params_err[i])


    def plot_fit_func_13Tev(self, filename_mc, filename_comp, decay, name_plot, save):
        # this function applies the TCM fit to the yield from run 2
        c1 = ROOT.TCanvas("c_fit", "c_fit", 900,900)
        pad = c1.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.2, 0.1, 0.1, 0.1)
        pad.SetLogy()
        pad.SetLogx()
        comp = Compare(filename_mc, self.decay, self.meson)
        
        h_corr_13 = comp.get_13tev_file(filename_comp, "CorrectedYieldPi0")

        h_corr_13.Scale(57.8*1e-3) #getting the differential invariant yield
        h_corr_13.Scale(1e12) #scaling to pico barn
        
        A_e = 536 * 1e9
        M_pi = 0.135
        T_e = 0.142
        A_init = 30 *1e9
        T_init = 0.63 
        n = 2.96
        
    
        f1TCM_13 = TF1("TCM_fit", "[0]* (TMath::Exp(-(TMath::Power(TMath::Power(x,2)+TMath::Power([1],2), 0.5)-[1])/([2])))+ [3] * (TMath::Power((1+(TMath::Power(x,2))/(TMath::Power([4],2)*[5])), -[5]))", 0.2, 15.,6)
        f1TCM_13.SetParameters(A_e, M_pi, T_e, A_init, T_init, n)
        parameters = ["A_e", "M_pi", "T_e", "A", "T", "n"]
        f1TCM_13.SetLineColor(ROOT.kOrange)
        
        

        
        f1TCM_13.SetParLimits(0, 0.5*A_e, 1.5*A_e)
        f1TCM_13.SetParLimits(1, M_pi, M_pi)
        f1TCM_13.SetParLimits(2, 0.5 * T_e, 1.5*T_e);
        f1TCM_13.SetParLimits(3,0.5*A_init, 1.5*A_init)
        f1TCM_13.SetParLimits(4, 0.5 * T_init, 1.5 * T_init)
        f1TCM_13.SetParLimits(5, 0.5 * n,  1.5*n)

        h_corr_13.Fit(f1TCM_13,"TCM","",0.2, 15);
        y = h_corr_13.GetYaxis() 
        y.SetTitle("E #frac{d^{3}#sigma}{dp^{3}} (pb GeV^{-2} c^{3})")
        y.SetTitleSize(0.048)
        y.SetTitleFont(42)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)    
        
        legend = ROOT.TLegend(0.575, 0.75, 0.8, 0.65)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.02)
        legend.AddEntry(h_corr_13, "corrected yield for 13TeV", "LP")
        legend.AddEntry(f1TCM_13, "Fit with two component model", "LP")
        legend.Draw()

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

        txt3 = TPaveText(0.875,0.82,0.875,0.82,"NDC"); 
        txt3.SetFillColor(kWhite)
        txt3.SetFillStyle(0)
        txt3.SetBorderSize(0)
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.02)
        txt3.AddText("pp, #sqrt{s} = 13.6TeV")
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

        txt2 = TPaveText(0.875,0.79,0.875,0.79,"NDC")
        txt2.SetFillColor(kWhite)
        txt2.SetFillStyle(0)
        txt2.SetBorderSize(0)
        txt2.SetTextAlign(33);#middle,left
        txt2.SetTextFont(42);#helvetica
        txt2.SetTextSize(0.02)
        txt2.AddText(decay)
        txt2.Draw()
        ROOT.SetOwnership(txt2,False)


        c1.Update()
        if save == "true":
            c1.SaveAs(name_plot.Data())
        
        params = f1TCM_13.GetParameters()
        params_err = f1TCM_13.GetParErrors()
        for i in range(6):
            print(f"Parameter {i} {parameters[i]}: {params[i]}", "±", params_err[i])






                


    