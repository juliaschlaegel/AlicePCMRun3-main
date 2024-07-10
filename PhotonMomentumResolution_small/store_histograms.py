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

    
# if __name__ == "__main__":
#     filename = "/Users/juliaschlagel/Analysis240411/Analysis/HLtrains/MC/AnalysisResults_HL197593_LHC24b1.root"
#     dir_name21 = "pi0eta-to-gammagamma-mc/Generated"
#     dir_name22_dalitz = "PCMDalitzEE"              
#     dir_name22_pcm ="PCMPCM"    
#     dir_name_py_1 = "associate-mc-info/Generated"
    
#     #hist = self.filename.Get(dir_name21).FindObject(dir_name22).FindObject(hist_name_all)

#     filename1_MC_dalitz = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC24b1_invariant_mass_dalitz/this_analysis_LHC24b1_pi0_mc_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
#     filename1_MC_pcm = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC24b1_invariant_mass_pcm/this_analysis_LHC24b1_pi0_mc_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
#     filename1_data_dalitz = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC22o_pass6_small_invariant_mass_dalitz/this_analysis_LHC22o_pass6_small_pi0_data_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
#     filename1_data_pcm = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/this_analysis_LHC22o_pass6_small_invariant_mass_pcm/this_analysis_LHC22o_pass6_small_pi0_data_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"

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
#     acc_class_dalitz = HistogramAcceptance(filename, dir_name21, dir_name22_dalitz, dir_ev, config)
#     acc_class_pcm = HistogramAcceptance(filename, dir_name21, dir_name22_pcm, dir_ev, config)
#     eff_dalitz = Efficiency(filename1_MC_dalitz, filename, dir_name11_dalitz, dir_name12_dalitz, dir_name13_mc, dir_name14, dir_name21, dir_name22_dalitz, acc_class_dalitz)
#     eff_pcm = Efficiency(filename1_MC_pcm, filename, dir_name11_pcm, dir_name12_pcm, dir_name13_mc, dir_name14, dir_name21, dir_name22_pcm, acc_class_pcm)



#     #corr_dalitz = CorrectedYieldDalitz(filename, filename1_MC_dalitz, filename1_data_dalitz) #, acc_class_dalitz, eff_dalitz)
#     #corr_pcm = CorrectedYieldDalitz(filename, filename1_MC_pcm, filename1_data_dalitz)
#     # Erstellen einer Instanz der zentralen Manager-Klasse
#     save_hist = Save_Hist("results.root")

    # Erstellen von Histogrammen mit verschiedenen Klassen
    # eff = Efficiency(filename1_MC_dalitz, filename, dir_name11_dalitz, dir_name12_dalitz, dir_name13_mc, dir_name14, dir_name21, dir_name22, acc_class)
    # acc = HistogramAcceptance(filename, dir_name21, dir_name22_dalitz, dir_ev, config)
    # corr = CorrectedYieldDalitz(filename, filename1_MC_dalitz, filename1_data_dalitz)
    # creator_a = HistogramCreatorA()
    # creator_b = HistogramCreatorB()
    #save_hist = Save_Hist("output.root")
# save_hist.add_histogram("Folder1", hist1, "new_hist1_name")
# save_hist.add_histogram("Folder1", hist2)  # Ohne Namensänderung
# save_hist.save_histograms()

    # efficiency_dalitz = eff_dalitz.calc_efficiency_vs_pt(hist_name_acc, hist_name_raw)
    # efficiency_pcm = eff_pcm.calc_efficiency_vs_pt(hist_name_acc, hist_name_raw)
    # acceptance_dalitz = acc_class_dalitz.calc_acceptance_vs_pt(hist_name_acc, hist_name_all)
    # acceptance_pcm = acc_class_pcm.calc_acceptance_vs_pt(hist_name_acc, hist_name_all)
    # raw_yield_MC_dalitz = eff_dalitz.get_raw_yield(hist_name_raw)
    # raw_yield_MC_pcm = eff_pcm.get_raw_yield(hist_name_raw)
    # raw_yield_data_dalitz = corr_dalitz.get_raw_yield_data(filename1_data_dalitz, dir_name11_dalitz, dir_name12_dalitz, dir_name13_data, dir_name14, hist_name_raw)
    # raw_yield_data_pcm = corr_pcm.get_raw_yield_data(filename1_data_pcm, dir_name11_pcm, dir_name12_pcm, dir_name13_data, dir_name14, hist_name_raw)
    # corrected_yield_dalitz = corr_dalitz.calculate_corrected_yield_dalitz(filename1_data_dalitz, hist_name_acc, hist_name_all, hist_name_raw, acc_class_dalitz, eff_dalitz, dir_name11_dalitz, dir_name12_dalitz, dir_name13_data, dir_name14)
    # corrected_yield_pcm = corr_pcm.calculate_corrected_yield_pcm(filename1_data_pcm, hist_name_acc, hist_name_all, hist_name_raw, acc_class_pcm, eff_pcm, dir_name11_pcm, dir_name12_pcm, dir_name13_data, dir_name14)

    #histB = creator_b.create_histogram()

    # Hinzufügen der Histogramme zum Manager
    # save_hist.add_histogram("PCMDalitzEE", efficiency_dalitz)
    # save_hist.add_histogram("PCMDalitzEE", acceptance_dalitz)
    # save_hist.add_histogram("PCMDalitzEE", raw_yield_data_dalitz, "h_raw_yield_data")
    # save_hist.add_histogram("PCMDalitzEE", raw_yield_MC_dalitz, "h_raw_yield_MC")
    # save_hist.add_histogram("PCMDalitzEE", corrected_yield_dalitz)

    # save_hist.add_histogram("PCMPCM", efficiency_pcm)
    # save_hist.add_histogram("PCMPCM", acceptance_pcm)
    # save_hist.add_histogram("PCMPCM", raw_yield_data_pcm, "h_raw_yield_data")
    # save_hist.add_histogram("PCMPCM", raw_yield_MC_pcm, "h_raw_yield_MC")
    # save_hist.add_histogram("PCMPCM", corrected_yield_pcm)
    

    # Speichern der Histogramme in einer ROOT-Datei
    #save_hist.save_histograms()

    # Laden der Histogramme zur Überprüfung
    # loaded_histA = manager.load_histogram("histA")
    # loaded_histB = manager.load_histogram("histB")
