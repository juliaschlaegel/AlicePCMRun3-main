import ROOT
from utility import Utility
import yaml
import numpy as np
from math import sqrt
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency

ROOT.TH1.SetDefaultSumw2()
ROOT.gROOT.SetBatch(True)

config_file ="/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC24b1_full_dalitz.yml" 
with open(config_file, "r", encoding="utf-8") as config_yml:
    config = yaml.safe_load(config_yml)

filename_mc = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/analysis_results_final_pi0_begin_0/inv_yield_calculation/Dalitz/this_analysis_pi0eta-to-gammagamma-pcmdalitzee_results.root"
rootfile_mc = TFile.Open(filename_mc, "READ")
rootdir = rootfile_mc.Get("PCMDalitzEE")
raw = rootdir.Get("h1_yield_mc")
if raw:
    print("raw yield loaded")

filename = "/Users/juliaschlagel/Analysis240411/Analysis/HLtrains/MC/AnalysisResults_LHC24b1_full_statistics.root"
rootfile = TFile.Open(filename, "READ")

dir_name1 = rootfile.Get("pi0eta-to-gammagamma-mc-pcmdalitzee/Generated")
dir_name2 = dir_name1.Get("Pi0")

histo = dir_name2.Get("hPt")
if histo:
    print("hPt loaded")
utils = Utility("pi0", filename, config, "dalitz", "mc")
h_2 = utils.rebin_hist("hPt")
# if h2:
#     print("histogram rebinned")

pt = utils.get_bin_var()
print(np.array(pt))

h1 = ROOT.TH1F("h1", "h1", len(pt)-1, np.array(pt))
h2 = ROOT.TH1F("h2", "h2", len(pt)-1, np.array(pt))

raw_hist = raw.Clone("raw_hist")

# h1.SetBinContent(1, 24)
# h2.SetBinContent(1, 155651658.0)
# h1.SetBinContent(2, 112)
# h2.SetBinContent(2, 59532536.0)
h1.SetBinContent(1, 0)
h2.SetBinContent(1,0)
h1.SetBinContent(2, 12)
h2.SetBinContent(2, 77825829.0)
h1.SetBinContent(3, 28)
h2.SetBinContent(3, 14883134.0)
h1.SetBinContent(4, 116)
h2.SetBinContent(4, 38272084.0)
h1.SetBinContent(5, 68)
h2.SetBinContent(5, 15074734.0)
h1.SetBinContent(6, 68)
h2.SetBinContent(6, 6334180.0)
h1.SetBinContent(7, 50)
h2.SetBinContent(7, 2893210.0)
h1.SetBinContent(8, 15)
h2.SetBinContent(8, 1081112.0)
h1.SetBinContent(9, 7)
h2.SetBinContent(9, 325950.0)
h1.SetBinContent(10, 2.5)
h2.SetBinContent(10, 87452.5)
h1.SetBinContent(11, 0)
h2.SetBinContent(11, 4942.538461538462)


# h1.SetBinError(1, sqrt(24))
# h2.SetBinError(1, sqrt(155651658.0))
# h1.SetBinError(2, sqrt(112))
# h2.SetBinError(2, sqrt(59532536.0))
h1.SetBinError(1, 0)
h2.SetBinError(1,0)
h1.SetBinError(2, sqrt(12))
h2.SetBinError(2, sqrt(77825829.0))
h1.SetBinError(3, sqrt(28))
h2.SetBinError(3, sqrt(14883134.0))
h1.SetBinError(4, sqrt(116))
h2.SetBinError(4, sqrt(38272084.0))
h1.SetBinError(5, sqrt(68))
h2.SetBinError(5, sqrt(15074734.))
h1.SetBinError(6, sqrt(68))
h2.SetBinError(6, sqrt(6334180.0))
h1.SetBinError(7, sqrt(50))
h2.SetBinError(7, sqrt(2893210.0))
h1.SetBinError(8, sqrt(15))
h2.SetBinError(8, sqrt(1081112.0))
h1.SetBinError(9, sqrt(7))
h2.SetBinError(9, sqrt(325950.0))
h1.SetBinError(10, sqrt(2.5))
h2.SetBinError(10, sqrt(87452.5))
h1.SetBinError(11, sqrt(0))
h2.SetBinError(11, sqrt(4942.538461538462))


teff = ROOT.TEfficiency(h1, h2)
teff.SetStatisticOption(ROOT.TEfficiency.kFCP)
#teff.SetStatisticOption(ROOT.TEfficiency.kBBayesian)

teff_raw = ROOT.TEfficiency(raw_hist, h_2)
teff_raw.SetStatisticOption(ROOT.TEfficiency.kFCP)

# eff = ROOT.TGraphAsymmErrors(raw, h2)
# eff.Print()

for i in range (1,len(pt)-1):
    print("Content raw yield: ", raw.GetBinContent(i), "±", raw.GetBinError(i), raw.GetBinCenter(i))
    #print("sqrt error: ", sqrt(raw.GetBinContent(i)))
#     if raw.GetBinContent(i) != 0:
#         print("relative: ", raw.GetBinError(i)/raw.GetBinContent(i))
    print("Content von rebinned hPt: ", h_2.GetBinContent(i), "±", h_2.GetBinError(i), h_2.GetBinCenter(i))
    print("sqrt error rebinned hPt: ", sqrt(h_2.GetBinContent(i)))
    # print("Content von hPt: ", histo.GetBinContent(i), "±", histo.GetBinError(i))
    # print("sqrt error hPt: ", sqrt(histo.GetBinContent(i)))

#     if h2.GetBinContent(i) != 0:
#         print("relative: ", h2.GetBinError(i)/h2.GetBinContent(i))
    print("Teff: ", teff.GetEfficiency(i),teff.GetEfficiencyErrorLow(i),teff.GetEfficiencyErrorUp(i))
    print("Teff with raw: ", teff_raw.GetEfficiency(i), teff_raw.GetEfficiencyErrorLow(i), teff_raw.GetEfficiencyErrorUp(i))
    # if teff.GetEfficiency(i) != 0:
    #     print("relative: ", teff.GetEfficiencyErrorLow(i)/teff.GetEfficiency(i), teff.GetEfficiencyErrorUp(i)/teff.GetEfficiency(i))
        
#     # print("eff: ", eff.GetBinContent(i), "±", eff.GetBinError(i))
#     # if eff.GetBinContent(i) != 0:
#     #     print("relative: ", eff.GetBinError(i)/eff.GetBinContent(i))
    

