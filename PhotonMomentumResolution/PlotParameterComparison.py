# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3
# This code was written by Alica Enderich (Febuary 2024)
# This code was modified by Julia Schlägel (July 2024)

import numpy as np
import yaml
import datetime
import ROOT
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TGaxis, TPad, TLatex, TBox, TString, TLine
from ROOT import gStyle, gROOT, gSystem
from ROOT import TFile, TDirectory, THashList, TF1, TH1F, TCanvas, TH1
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kFullCircle, kOrange, kTRUE
import re
import numpy as np
import datetime
import math
import ctypes
import os, sys, shutil
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);
gStyle.SetPadTickX(1);
gStyle.SetPadTickY(1);
gStyle.SetErrorX(0)
gStyle.SetEndErrorSize(5)
ROOT.TGaxis.SetMaxDigits(2)

from HistoFormatting import FrameSettings, CanvasSettings, PadSettings, DrawHisto, SetTitle, SetStyleTLatex

def make_common_style(g1,marker,size,color,width=1,fill=0):
    g1.SetMarkerStyle(marker);
    g1.SetMarkerColor(color);
    g1.SetMarkerSize(size);
    g1.SetLineColor(color);
    g1.SetLineWidth(width);
    g1.SetFillColorAlpha(color, 0.65);
    g1.SetFillStyle(fill);

def draw_comparison_parameters(meson, mean, list_data, list_mc, index, cutname, date, period_str, decay):
        ymin = [6*1e-3, 130*1e-3,  10*1e-4, 2e-3, 1e-8]
        ymax = [22*1e-3, 140*1e-3 , 110*1e-4, 14e-3, 1e-3]
        for i in range(len(list_data)):
            print("MAX", list_data[i].GetMaximum(),list_mc[i].GetMaximum())
    #canvas plotting
        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);
        if list_data[0].GetTitle() == "raw yield":
            p1.SetLogy()

        frame1 = p1.DrawFrame(0, ymin[index], 10., ymax[index]);
        frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/c)");
        if cutname == "mean": 
            unit = "[MeV/c^{2}]"
        if cutname == "FWHM/2.36":
            unit = "[GeV/c^{2}]"
        if cutname == "lambda":
            unit = "[GeV/c^{2}]"
        else:
            unit = "[GeV/c^{2}]"

        frame1.GetYaxis().SetTitle("{0} {1}".format(cutname, unit));
        frame1.GetXaxis().SetTitleSize(0.045);
        frame1.GetYaxis().SetTitleSize(0.045);
        frame1.GetXaxis().SetTitleOffset(1.0);
        frame1.GetYaxis().SetTitleOffset(1.52);
        frame1.GetXaxis().SetLabelSize(0.045);
        frame1.GetYaxis().SetLabelSize(0.045);
        frame1.GetYaxis().SetMaxDigits(3);
        frame1.GetXaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetNdivisions(5, 5, 0, kTRUE)
        ROOT.SetOwnership(frame1,False);

        color = [kRed+1, kBlue+1, kGreen+2, kMagenta+2, kCyan+1]
        for i in range(len(list_data)):
                print(i, list_data[i].GetTitle(), len(list_data))
                make_common_style(list_data[i], 20, 0.9, color[i], 1, 0);
                make_common_style(list_mc[i], 25, 0.9, color[i], 1, 0);
            
                list_data[i].Draw("E0same");
                list_mc[i].Draw("E0same");
        if mean == "true":
            line2 = TLine(0.,0.13497,10.,0.13497);
            line2.SetLineColor(kRed);
            line2.SetLineStyle(2);
            line2.SetLineWidth(1);
            line2.Draw("");
            ROOT.SetOwnership(line2,False)

        txt = TPaveText(0.,0.9,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);
        txt.AddText("{0} for different cuts in data and MC".format(cutname));
        txt.Draw();
        ROOT.SetOwnership(txt,False);


        txt = TPaveText(0.50,0.7,0.88,0.82,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(32);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.03);
        txt.AddText("this thesis");
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
        txt.AddText(decay)
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        leg = TLegend(0.5,0.65,0.9,0.7);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.03);
        leg.SetTextAlign(32);
        leg.SetTextFont(42);#helvetica
        
        if decay == "pcm":
            cut = ["pi0eta-to-gamma-gamma-pcmpcm"]
        if decay == "dalitz":
            cut = ["pi0eta-to-gamma-gamma-pcmdalitzee", "pi0eta-to-gamma-gamma-pcmdalitzee-itsibany"]
        for i in range(len(list_data)):
            leg.AddEntry(list_data[i], "{0}".format(cut[i]),"LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        leg = TLegend(0.2,0.72,0.50,0.82);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.03);
        leg.SetTextAlign(12);
        leg.SetTextFont(42);#helvetica
        leg.AddEntry(list_data[0], "Data #gamma candidates (LHC220)","LP");
        leg.AddEntry(list_mc[0] , "MC rec. #gamma (LHC24b1)","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        p2 = c1.cd(2);
        p2.SetPad(0,0,1,0.3);
        p2.SetMargin(0.15,0.02,0.22,0.0);
        p2.SetTicks(1,1);
        cut_in_ratio = 0.5
        ymin, ymax = 0.5, 2.0
        if index ==1:
            ymin, ymax = 0.95, 1.05

        frame2 = p2.DrawFrame(0,ymin,10.,ymax);
        frame2.GetXaxis().SetTitle("#it{p}_{T} (GeV/c)");
        frame2.GetYaxis().SetTitle("#frac{Data}{MC}");
        frame2.GetXaxis().SetTitleSize(0.10);
        frame2.GetYaxis().SetTitleSize(0.10);
        frame2.GetXaxis().SetTitleOffset(1.0);
        frame2.GetYaxis().SetTitleOffset(0.7);
        frame2.GetXaxis().SetLabelSize(0.10);
        frame2.GetYaxis().SetLabelSize(0.10);
        frame2.GetYaxis().CenterTitle(True);
        frame2.GetXaxis().SetLabelOffset(0.01);
        frame2.GetYaxis().SetLabelOffset(0.01);
        frame2.GetYaxis().SetNdivisions(5, 5, 0, kTRUE)
        ROOT.SetOwnership(frame2,False);

        color = [kRed+1, kBlue+1, kGreen+2, kMagenta+2, kCyan+1]
        for i in range(len(list_data)):
                h1ratio = list_data[i].Clone("h1ratio");
                h1ratio.Reset();
                h1ratio.Sumw2();
                h1ratio.Divide(list_data[i], list_mc[i], 1., 1., "G");
                h1ratio.Draw("E0same");
                h1ratio.SetDirectory(0);
                ROOT.SetOwnership(h1ratio,False);

        line1 = TLine(0.,1,10.,1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        filepath = os.path.join(folder, "{0}_Parameter_Comparison_{1}{2}.pdf".format(date, period_str, cutname));    
        c1.SaveAs(filepath);

def run(filename_data, filename_mc,config,type,folder, period_str, decay):
    print(sys._getframe().f_code.co_name);
    arr_pt = np.array(config["common"]["pt_bin"],dtype=float);
    print("pT binning = ",arr_pt);
    print("type = ",type);
    print("input = ",filename_data);
    rootfile_data = TFile.Open(filename_data,"READ");
    rootfile_mc = TFile.Open(filename_mc,"READ");
    meson = config["common"]["meson"];

    list_fit_func = config[type]["fit_func"];

    list_fit_min = config["common"]["fit_min"];
    list_fit_max = config["common"]["fit_max"];
    if len(list_fit_min) != len(list_fit_max):
        return;

    list_integral_min = config["common"]["integral_min"];
    list_integral_max = config["common"]["integral_max"];
    if len(list_integral_min) != len(list_integral_max):
        return;

    list_yield_min = config["common"]["yield_min"];
    list_yield_max = config["common"]["yield_max"];
    if len(list_yield_min) != len(list_yield_max):
        return;

    nsys = len(config[type]['subsystems']);
    print("nsys",nsys); 

    def loop_data(rootfile): 
            for isys in range(0,nsys):
                ssname = config[type]['subsystems'][isys]['name']; #subsystem name
                print("plot subsystem", ssname);
                
                nfit = len(list_fit_func)
                list_parameters_comparison = [];
                yield_list = [];
                list_ss   = rootfile.Get(ssname);     
                for ifit in range(0, nfit): 
                    fitname = "gausexplinear";
                    list_fitname  = list_ss.FindObject(fitname);
                    print("fitname ", fitname)

                    list_plot = list_fitname.FindObject("fit_0.04_0.20_GeVc2")
                    yield_list.append(list_plot.FindObject("h1yield_param"))

                    parameter_list = [];
                    it = ROOT.TIter(list_plot)
                    obj = it.Next()
                    while obj:
                        objName = obj.GetName()
                        if "param" in objName:
                            parameter_list.append(obj)
                        obj = it.Next()
                    list_parameters_comparison.append(parameter_list);
            return list_parameters_comparison, yield_list  
    def loop_mc(rootfile): 
            for isys in range(0,nsys):
                ssname = config[type]['subsystems'][isys]['mc_name']; #subsystem name
                print("plot subsystem", ssname);
                nfit = len(list_fit_func)
                list_parameters_comparison = [];
                yield_list = [];
                list_ss   = rootfile.Get(ssname);
                print("ss list: ", list_ss)
             
                for ifit in range(0, nfit): 
                    fitname = "gausexp";
                    list_fitname  = list_ss.FindObject("gausexp");
                    print("fitname ", fitname)

                    list_plot = list_fitname.FindObject("fit_0.04_0.20_GeVc2")
                    yield_list.append(list_plot.FindObject("h1yield_param"))

                    parameter_list = [];
                    same_list = [];
                    mixed_list = [];
                    parameter_list = [];
                    it = ROOT.TIter(list_plot)
                    obj = it.Next()
                    while obj:
                        objName = obj.GetName()
                        if "param" in objName:
                            parameter_list.append(obj)
                        obj = it.Next()
                    list_parameters_comparison.append(parameter_list);
            return list_parameters_comparison, yield_list     
    list_parameters_data, yield_list_data = loop_data(rootfile_data)
    
    
    

    list_parameters_mc, yield_list_mc = loop_mc(rootfile_mc)
    
    for i in range(len(list_parameters_mc)):
         print(len(list_parameters_mc[0]))
         print(list_parameters_mc[i][0].GetTitle())
         print(list_parameters_mc[i][1].GetTitle())
         print(list_parameters_mc[i][2].GetTitle())
         print(list_parameters_mc[i][3].GetTitle())
         print(list_parameters_mc[i][4].GetTitle())
         print(list_parameters_mc[i][5].GetTitle())

    for i in range(len(list_parameters_data)):
        for j in range(len(list_parameters_data[0])):
            print(list_parameters_data[i][j].GetTitle())

    list_data_mean = []
    list_data_lambda = []
    list_data_sigma = []
    list_data_fwhm = []
    for i in range(len(list_parameters_data)):
        for j in range(len(list_parameters_data[0])):
            if list_parameters_data[i][j].GetTitle() == "mean":
                list_data_mean.append(list_parameters_data[i][j])
            if list_parameters_data[i][j].GetTitle() == "lambda":
                list_data_lambda.append(list_parameters_data[i][j])
            if list_parameters_data[i][j].GetTitle() == "sigma":
                list_data_sigma.append(list_parameters_data[i][j])
            if list_parameters_data[i][j].GetTitle() == "fwhm/2.36":
                 list_data_fwhm.append(list_parameters_data[i][j])

    list_mc_mean = []
    list_mc_lambda = []
    list_mc_sigma = []
    list_mc_fwhm = []
    for i in range(len(list_parameters_mc)):
        for j in range(len(list_parameters_mc[0])):
            if list_parameters_mc[i][j].GetTitle() == "mean":
                list_mc_mean.append(list_parameters_mc[i][j])
            if list_parameters_mc[i][j].GetTitle() == "lambda":
                list_mc_lambda.append(list_parameters_mc[i][j])
            if list_parameters_mc[i][j].GetTitle() == "sigma":
                list_mc_sigma.append(list_parameters_mc[i][j])
            if list_parameters_mc[i][j].GetTitle() == "fwhm/2.36":
                 list_mc_fwhm.append(list_parameters_mc[i][j])
    
    print(list_data_lambda)
    print(list_data_mean)
    print(list_data_sigma)
    print(list_data_fwhm)
    print(list_mc_lambda)
    print(list_mc_mean)
    print(list_mc_sigma)
    print(list_mc_fwhm)

    list_data = [list_data_lambda, list_data_mean, list_data_sigma, list_data_fwhm, yield_list_data]
    list_mc = [list_mc_lambda, list_mc_mean, list_mc_sigma, list_mc_fwhm, yield_list_mc]

    print("length", len(list_data), (len(list_data_lambda,)), (len(list_data_mean)), (len(list_data_sigma)), len(list_data_fwhm), len(list_parameters_data))
    print("length", len(list_mc), (len(list_mc_lambda,)), (len(list_mc_mean)), (len(list_mc_sigma)), len(list_mc_fwhm), len(list_parameters_mc))
    for i in range(5):
        list_mc[i][0].SetBinContent(1, 0)
        list_mc[i][0].SetBinError(1,0)
        list_data[i][0].SetBinContent(1, 0)
        list_data[i][0].SetBinError(1,0)

        print(list_data[i], list_mc[i], i, list_data[i][0].GetTitle())
        if i == 1:
            draw_comparison_parameters("pi0", "true",list_data[i], list_mc[i], i, list_data[i][0].GetTitle(), date, period_str, decay)
        else:
            draw_comparison_parameters("pi0", "false", list_data[i], list_mc[i], i, list_data[i][0].GetTitle(), date, period_str, decay)
   
    for i in range(len(list_data[3][0])):
        bc = list_data[3][0].GetBinContent(i)
        center = list_data[3][0].GetBinCenter(i)
        print("##bin content: ", bc, center)
    draw_comparison_parameters("pi0", "false", list_data[3], list_mc[3], 3, "FWHM div 2.36", date, period_str, decay)

if __name__ == "__main__":
   
    period_array = ["LHC22o", "LHC24b1"]
    filename_array = ["/Users/juliaschlagel/Analysis240411/Analysis/HLtrains/data/AnalysisResults_LHC22o_full_statistics.root", "/Users/juliaschlagel/Analysis240411/Analysis/HLtrains/MC/AnalysisResults_LHC24b1_full_statistics.root"]
    type_array = ["data", "mc"]
    decay = "dalitz"
    config_array = ["/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC22o_full_dalitz.yml", "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC24b1_full_dalitz.yml"] #dalitz
    # config_array = ["/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC22o_full_pcm.yml", "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC24b1_full_pcm.yml"] #pcm

    period_str = ""
    for i in range(len(period_array)):
        period_str += period_array[i]
        period_str += "_"
        print("period_str", period_str)

    for i in range(len(period_array)):
        period = period_array[i]
        filename = filename_array[i]
        cutname = "qc"
        suffix = "AnyTrack";
        type = type_array[i]
        config_file = config_array[i]
        with open(config_file, "r", encoding="utf-8") as config_yml:
            config = yaml.safe_load(config_yml)
        # Date or prefix "this_thesis"
        date = "this_thesis" #datetime.date.today().strftime("%Y%m%d");
        folder = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/param_comp_{0}_small_binning".format(decay)
        os.makedirs(folder, exist_ok=True);
        # filename_data = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/analysis_results_pi0_begin_1/inv_mass_analysis/PCM/this_analysis_pcm_LHC22o_pi0_data_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root" #pcm
        filename_data = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/analysis_results_pi0_begin_1/inv_mass_analysis/Dalitz/this_analysis_dalitz_LHC22o_pi0_data_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
        # filename_mc = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/analysis_results_pi0_begin_1/inv_mass_analysis/PCM/this_analysis_pcm_LHC24b1_pi0_mc_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root" #pcm
        filename_mc = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/analysis_results_pi0_begin_1/inv_mass_analysis/Dalitz/this_analysis_dalitz_LHC24b1_pi0_mc_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
        run(filename_data, filename_mc,config,type,folder, period_str, decay)
        print("input period_str", period_str)