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
# #gStyle.SetPadTickX(1)
# gStyle.SetPadTickY(1)
# gStyle.SetErrorX(0)
# gStyle.SetEndErrorSize(5)
#from HistoFormatting import FrameSettings, CanvasSettings, PadSettings, DrawHisto, SetTitle, SetStyleTLatex
from utility import Utility



class HistogramAcceptance:
    def __init__(self, filename, config, decay): #, hist1_name, hist2_name):
        # Open wanted data and load histogramms
        self.decay = decay
        self.filename = TFile.Open(filename, "READ")
        dir_name1 = "pi0eta-to-gammagamma-mc/Generated"
        if self.decay == "pcm":
            dir_name2 = "PCMPCM"
        else:
            dir_name2 = "PCMDalitzEE"
        #self.dir_name1 = dir_name1
        #self.dir_name2 = dir_name2
        self.dir = self.filename.Get(dir_name1).FindObject(dir_name2)
        #self.dir_ev = dir_ev
        self.config = config
        #self.utils = Utility(filename, config, decay)
    # def __init(self, utils):
    #     self.utils = Utility(filename, config, decay)

        
    # def get_hist (self, hist_name):
    #     self.hist_name = self.dir.FindObject(hist_name)
    #     if not self.hist_name:
    #         raise ValueError("Histogram not found")
    #     return self.hist_name

    # def get_Nev (self): #dir_name2):
    #     # dir_ev = "pi0eta-to-gammagamma-mc/Event"
    #     direv = self.filename.Get("pi0eta-to-gammagamma-mc/Event")
    #     if self.decay == "pcm":
    #         dir_name2 = "PCMPCM"
    #     else:
    #         dir_name2 = "PCMDalitzEE"
    #     events = direv.FindObject(dir_name2).FindObject("after").FindObject("hCollisionCounter")
    #     self.Nev = events.GetBinContent(9)
    #     print("Nev :", self.Nev)
    #     return self.Nev
    
    #     # print(Nev)
    # def get_bin_var(self):
    #      #with open(self.config, "r") as file:
    #         #config_data = yaml.safe_load(file)
    #         self.bin_var = self.config["common"]["pt_bin"]
    #         print(self.bin_var)
    #         return self.bin_var


    # def rebin_hist(self, hist_name):
    #     bin_var = self.get_bin_var()
    #     hist = self.get_hist(hist_name)
    #     new_hist = TH1D("p_T1", "p_T1", len(bin_var) - 1, np.array(bin_var))
    
    #     for bin in range(0, len(bin_var) - 1):
    #         b1 = bin_var[bin]
    #         b2 = bin_var[bin + 1]
    #         # print("b1: ", b1, "b2: ", b2)

    #         bin_b1 = hist.GetXaxis().FindBin(b1 + 1e-6)
    #         bin_b2 = hist.GetXaxis().FindBin(b2 - 1e-6)
    #         # print("bin_b1: ", bin_b1, "bin_b2: ", bin_b2)

    #         error = c_double(0.0)
    #         content = hist.IntegralAndError(bin_b1, bin_b2, error, "")
    #         bin_center = (b1 + b2) / 2  # Take the center of the bin to fill
    #         bin_width = b2-b1
    #         new_content = content/ bin_width
    #         new_hist.Fill(bin_center, new_content)  # Filling with this center of each bin
            
    #         new_hist.SetBinError(bin+1, error.value / bin_width)  # calculate error to the right bin#

    #     NoE = self.get_Nev()
    #     new_hist.Scale(1 / NoE)
    #     if new_hist is None:
    #         print(f"Failed to load histogram: {hist_name}")
    #     # print("hist neu type:", type(new_hist))
    #     # print(f"Histogram {hist_name} loaded with {new_hist.GetEntries()} entries.")
    #     # return new_hist

    #     # Create canvas and draw histogram
    #     canvas = TCanvas("acc", "acc", 0, 0, 900, 900)
    #     pad = canvas.cd()
    #     pad.SetPad(0.0, 0.01, 1, 1)
    #     pad.SetMargin(0.15, 0.1, 0.1, 0.1)
    #     pad.SetTicks(1,1)
    #     new_hist.Draw("ESame")
    #     ROOT.SetOwnership(canvas, False)
    #     # canvas.SaveAs("rebinning_test1")

    #     return new_hist



    # def plot_combined_hist(self):
    #     h1 = self.rebin_hist(hist_name_acc)
    #     h2 = self.rebin_hist(hist_name_all)

    #     max_h1 = h1.GetMaximum()
    #     max_h2 = h2.GetMaximum()
    #     max_y = 1.05 * max(max_h1, max_h2)  # 20% mehr als der höchste Balken

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
    #     h1.Draw("hist")              # Zeichne das erste Histogramm

    #     # Einstellungen für das zweite Histogramm
    #     h2.SetLineColor(ROOT.kMagenta)   # Rote Linie
    #     h2.Draw("histsame")          # Zeichne das zweite Histogramm über das erste

    #     # Legende hinzufügen
    #     legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)  # Position der Legende
    #     legend.AddEntry(h1, "daughters in acceptance", "LP")
    #     legend.AddEntry(h2, "all", "LP")
    #     legend.Draw()

    #     # Achsenbeschriftungen
    #     h1.SetTitle("Combined Histogram")
    #     h1.SetXTitle("p_T")
    #     h1.SetYTitle("rebinned histograms for the acceptance")

    #     c1.Update()
    #     c1.SaveAs("combined_histogram_acc")

        

    # #def calc_acceptance_vs_pt(self, name_plot, save):
    #     #Get wanted histograms:
    #     # Re-bin histograms
    #     h1 = self.rebin_hist("hPt_Pi0_Acc")
    #     # Error_1 = h1.GetBinError()
    #     # print("Bin error hPt_acc: ", Error_1)
    #     print("First Histogram rebinned")
    #     h2 = self.rebin_hist("hPt_Pi0")
    #     #Error_2 = h2.GetBinError()
    #     print("Bin error hPt all: ", Error_2)
    #     print("second histogram rebinned")
    #     if h1 is None:
    #         print(f"Failed to load histogram: hPt_Pi0")
    #     # print("h1 type:", type(h1))

    # #get canvas and pad to plot the acceptance in 
    #     canv = TCanvas("acc1", "acc1", 0, 0, 900, 900)
    #     pad = canv.cd()
    #     pad.SetPad(0.0, 0.01, 1, 1)
    #     pad.SetMargin(0.15, 0.1, 0.1, 0.1)
    #     pad.SetTicks(1,1)

    #     h_acc = h1.Clone("h_A")   
        
    #     # calculate acceptance by dividing the events with daughters in acceptance by all events
    #     # this time in percent
    #     h_acc.Divide(h2)
    #     #h_acc.Divide(h1, h2, 1, 1 )#, option="B") ###Fragen: Welcher Fehler? und in prozent?
    #     # Define marker settings
    #     h_acc.SetFillColor(kCyan+1)
    #     h_acc.SetMarkerStyle(kFullCross)
    #     h_acc.SetMarkerColor(kBlue)
    #     h_acc.SetLineColor(kBlue)

    #     #Draw points in the histogramm
    #     h_acc.Draw("Esame")
       
    #     # Define y-axis settings
    #     y = h_acc.GetYaxis() #.SetRangeUser(min(bin_var), max(bin_var))
    #     y.SetTitle("Acceptance") #_{e^{+}e^{-}\gamma}^{\pi^0}
    #     y.SetTitleSize(0.048)
    #     y.SetTitleFont(42)
    #     y.SetLabelFont(42)
    #     y.SetLabelSize(0.035)

    #     #Define x-axis settings
    #     x = h_acc.GetXaxis()
    #     x.SetTitle("p_{T} [GeV/c]")
    #     x.SetTitleSize(0.048)
    #     x.SetTitleFont(42)
    #     x.SetTitleOffset(0.9)
    #     x.SetLabelFont(42)
    #     x.SetLabelSize(0.035)

    #     # Define legend settings
    #     leg = TLegend(0.6, 0.5, 1.0, 0.75)
    #     leg.SetBorderSize(0)
    #     leg.SetFillColor(kWhite)
    #     leg.SetFillStyle(0)
    #     leg.SetTextSize(0.03)
    #     if self.decay == "pcm":
    #         leg.AddEntry(h_acc, "\pi^{0} \\rightarrow \gamma \gamma", "LP")
    #     else:
    #         leg.AddEntry(h_acc, "\pi^{0} \\rightarrow e^{+} e^{-} \gamma", "LP")
    #     leg.Draw("")
    #     ROOT.SetOwnership(leg,False)
               
    #     # Add text 
    #     txt = TPaveText(0.875,0.2,0.875,0.2,"NDC")
    #     txt.SetFillColor(kWhite)
    #     txt.SetFillStyle(0)
    #     txt.SetBorderSize(0)
    #     txt.SetTextAlign(33);#middle,left
    #     txt.SetTextFont(42);#helvetica
    #     txt.SetTextSize(0.02)
    #     txt.AddText("this thesis")
    #     txt.Draw()
    #     ROOT.SetOwnership(txt,False)

    #     txt3 = TPaveText(0.875,0.17,0.875,0.17,"NDC"); # txt3 = TPaveText(0.4,0.75,0.4,0.90,"NDC");
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
    #     canv.Modified()
    #     canv.Update()
    #     ROOT.SetOwnership(canv, False)
    #     if save == "true":
    #         canv.SaveAs(name_plot.Data())
    #     # else:
    #     #     if self.decay =="pcm":
    #     #         canvas.SaveAs("Acceptance_pcm")
    #     #     else:
    #     #         canvas.SaveAs("Acceptance_dalitz")
    #     # print("acc hist type:", type(h_acc))
    #     # return h_acc
    

    def calc_acceptance_vs_pt(self, name_plot, save, utils):
        #Get wanted histograms:
        # Re-bin histograms
        h1 = utils.rebin_hist("hPt_Pi0_Acc")
        utils.scale_by_Nev(h1)
        #h1 = self.rebin_hist("hPt_Pi0_Acc")
        print("First Histogram rebinned")
        h2 = utils.rebin_hist("hPt_Pi0")
        utils.scale_by_Nev(h2)
        #h2 = self.rebin_hist("hPt_Pi0")
        print("second histogram rebinned")
        if h1 is None:
            print(f"Failed to load histogram: hPt_Pi0")
        # print("h1 type:", type(h1))

    #get canvas and pad to plot the acceptance in 
        canv = TCanvas("acc1", "acc1", 0, 0, 900, 900)
        pad = canv.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()
        #pad.SetLogx()
        

        h_acc = h1.Clone("h_A")   
        
        # calculate acceptance by dividing the events with daughters in acceptance by all events
        # this time not in percent
        h_acc.Divide(h1, h2, 1, 1 )#, option="B") ###Fragen: Welcher Fehler? und in prozent?
        # Define marker settings
        h_acc.SetFillColor(kCyan+1)
        h_acc.SetMarkerStyle(kFullCross)
        h_acc.SetMarkerColor(kBlue)
        h_acc.SetLineColor(kBlue)

        #Draw points in the histogramm
        h_acc.Draw("Esame")
       
        # Define y-axis settings
        y = h_acc.GetYaxis() #.SetRangeUser(min(bin_var), max(bin_var))
        y.SetTitle("Acceptance") #_{e^{+}e^{-}\gamma}^{\pi^0}
        y.SetTitleSize(0.048)
        y.SetTitleFont(42)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)

        #Define x-axis settings
        x = h_acc.GetXaxis()
        x.SetTitle("p_{T} [GeV/c]")
        x.SetTitleSize(0.048)
        x.SetTitleFont(42)
        x.SetTitleOffset(0.9)
        x.SetLabelFont(42)
        x.SetLabelSize(0.035)
        x.SetRangeUser(0,14)

        # Define legend settings
        leg = TLegend(0.2, 0.8, 0.4, 0.75)
        leg.SetBorderSize(0)
        leg.SetFillColor(kWhite)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.03)
        if self.decay == "pcm":
            leg.AddEntry(h_acc, "\pi^{0} \\rightarrow \gamma \gamma", "LP")
        else:
            leg.AddEntry(h_acc, "\pi^{0} \\rightarrow e^{+} e^{-} \gamma", "LP")
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
        canv.Modified()
        canv.Update()
        ROOT.SetOwnership(canv, False)
        if save == "true":
            canv.SaveAs(name_plot.Data())
        if save == "loc":
            canv.SaveAs("Acceptance_with_ut")
        # else:
        #     if self.decay =="pcm":
        #         canvas.SaveAs("Acceptance_pcm")
        #     else:
        #         canvas.SaveAs("Acceptance_dalitz")
        print("acc hist type:", type(h_acc))
        return h_acc
        

# using the class
if __name__ == "__main__":
    filename = "/Users/juliaschlagel/Analysis240411/Analysis/HLtrains/MC/AnalysisResults_HL197593_LHC24b1.root"
    #dir_name1 = "pi0eta-to-gammagamma-mc/Generated"
    # dir_name2 = "PCMPCM"
    #dir_name2 = "PCMDalitzEE"
    #hist_name_all = "hPt_Pi0"
    #hist_name_acc = "hPt_Pi0_Acc"
    #dir_ev = "pi0eta-to-gammagamma-mc/Event"
    config_file = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC22o_small_dalitz.yml"
    with open(config_file, "r", encoding="utf-8") as config_yml:
                config = yaml.safe_load(config_yml)
    acc = HistogramAcceptance(filename, config, "dalitz")
    #acc.calc_acceptance_vs_pt("Accceptance", "loc")
