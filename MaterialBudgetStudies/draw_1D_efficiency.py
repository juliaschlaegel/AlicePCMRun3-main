# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import numpy as np
# import datetime
import math
import ROOT
import os
import yaml
from ctypes import *
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython
from ROOT import gStyle, gROOT, gSystem, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kTRUE
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);
from FomattingMaterialBudget import ALICEtext, FrameSettingsRatio, FrameSettings, RatioLegendSettings

#_____________________________________________________________________
def make_common_style(g1,marker,size,color,width=1,fill=0):
    g1.SetMarkerStyle(marker);
    g1.SetMarkerColor(color);
    g1.SetMarkerSize(size);
    g1.SetLineColor(color);
    g1.SetLineWidth(width);
    g1.SetFillColor(color);
    g1.SetFillStyle(fill);

class draw_1D_efficiency:
    def __init__(self):
        print("default constructor is called");
    def __init__(self, filename_data, filename_mc, folder, period, suffix):
        #self.cutname = cutname;
        self.folder = folder;
        self.period = period;
        self.suffix = suffix;

        self.rootfile_data = TFile.Open(filename_data, "READ");
        self.rootfile_mc   = TFile.Open(filename_mc  , "READ");

        # self.rootdir_mc_gen  = self.rootfile_mc.Get("material-budget-mc")
        # self.list_gen        = self.rootdir_mc_gen.Get("Generated");
        # self.list_ev_gen     = self.rootdir_mc_gen.Get("Event");
        # self.list_ev_mc_gen  = self.list_ev_gen.FindObject("PCMDalitzEE");
        # self.rootdir_mc_rec  = self.rootfile_mc.Get("material-budget-mc");
        # self.list_v0_mc_rec  = self.rootdir_mc_rec.Get("V0");
        # self.list_ev_rec     = self.rootdir_mc_rec.Get("Event")
        # self.list_ev_mc_rec  = self.list_ev_rec.FindObject("PCMDalitzEE");
        # self.list_cut_mc_rec = self.list_v0_mc_rec.FindObject(cutname);


        # self.rootdir_mc_pcmqc    = self.rootfile_mc.Get("pcm-qc-mc");
        # self.list_gen_pcmqc        = self.rootdir_mc_pcmqc.Get("Generated");
        # self.list_ev_mc_pcm  = self.rootdir_mc_pcmqc.Get("Event");
        # self.h1nch_mc_gen    = self.list_ev_mc_pcm.FindObject("hMultNTracksPV").Clone("h1mult");
        # self.nev_gen         = self.h1nch_mc_gen.GetEntries();
        # self.nch_gen         = self.h1nch_mc_gen.GetMean();

        # self.h1nch_mc_rec    = self.list_ev_mc_pcm.FindObject("hMultNTracksPV");
        # self.nch_rec         = self.h1nch_mc_rec.GetMean();
        # self.nev_rec         = self.h1nch_mc_rec.GetEntries();
    
        # self.rootdir_data    = self.rootfile_data.Get("material-budget");
        # self.list_v0_data    = self.rootdir_data.Get("V0");
        # self.list_ev_data_1  = self.rootdir_data.Get("Event");
        # self.list_ev_data    = self.list_ev_data_1.FindObject("PCMDalitzEE");
        # self.list_cut_data   = self.list_v0_data.FindObject(cutname);
        self.rootdir_data_pcmqc    = self.rootfile_data.Get("pcm-qc");
        self.list_v0 = self.rootdir_data_pcmqc.Get("V0")
        

        self.list_ev_data_pcm  = self.rootdir_data_pcmqc.Get("Event");
        self.list_ev_data_pcmqc_after = self.list_ev_data_pcm.Get("after")

        self.h1nch_data = self.list_ev_data_pcmqc_after.Get("hMultNTracksPV");
        self.nch_data = self.h1nch_data.GetMean();
        print("nch data: ", self.nch_data)

        self.h1nch_data      = self.list_ev_data_pcmqc_after.FindObject("hMultNTracksPV");
        if self.h1nch_data:
            print("h1nch found")
        #self.nev_data        = self.h1nch_data.GetEntries();
        #self.nch_data        = self.h1nch_data.GetMean();
        self.coll = self.list_ev_data_pcmqc_after.Get("hCollisionCounter")# .Clone("h1ev");
        print("type h1ev_data: ", type(self.coll))
        self.nev_data = self.coll.GetBinContent(10);
        print("NEV data: ", self.nev_data)

        self.suffix = suffix;
        self.folder = folder;
        #self.cutname = cutname
        self.arr_rxy = np.array([0,1,2,3,4,5], dtype=float);

    # def __del__(self):
    #     if self.rootfile_data.IsOpen():
    #         print("close input data root file.");
    #         self.rootfile_data.Close();
    #     if self.rootfile_mc.IsOpen():
    #         print("close input mc root file.");
    #         self.rootfile_mc.Close();
   
    #_____________________________________________________________________
    def draw_eff(self, cuts, generated, date):

        # hs_mc_rec = self.list_cut_mc_rec.FindObject("hs_conv_point").Clone("hs_rec");
        # h1_mc_rec = hs_mc_rec.Projection(1 ,"");
        # h1_mc_rec.SetDirectory(0);
        # ROOT.SetOwnership(h1_mc_rec, False);
        # h1_mc_rec.Sumw2();
        # h1_mc_rec.Scale(1,"width");
        # h1_mc_rec.Scale(1/self.nev_rec);
        # h1_mc_rec.Scale(1/self.nch_rec);#nch
        # make_common_style(h1_mc_rec, 20, 1.0, kRed+1, 1, 0)

        # h2_mc_gen = self.list_gen_pcmqc.FindObject("hPhotonRZ");
        # h2_mc_gen.Sumw2();

        # h1_mc_gen = h2_mc_gen.ProjectionY("h1_mc_gen");
        # h2_mc_gen.SetDirectory(0);
        # h1_mc_gen.SetDirectory(0);
        # h1_mc_gen.RebinX(10);
        # ROOT.SetOwnership(h2_mc_gen,False);
        # ROOT.SetOwnership(h1_mc_gen,False);
        # # print("BINS", h1_mc_gen.GetXaxis().GetNbins())

        # mc_rec_list = self.rootdir_mc_pcmqc.Get("V0")
        # qc_mc_rec = mc_rec_list.FindObject("qc")
        # h2_mc_rec_to_MC = qc_mc_rec.FindObject("hRZ_Photon_Primary_MC")
        # h1_mc_rec_to_MC = h2_mc_rec_to_MC.ProjectionY("h1_mc_recMC");
        # print("BINS", h1_mc_rec_to_MC.GetXaxis().GetNbins())

        # h1_mc_rec_to_MC.Scale(1,"width");
        # h1_mc_rec_to_MC.Scale(1/self.nev_rec);
        # h1_mc_rec_to_MC.Scale(1/self.nch_rec);
        # make_common_style(h1_mc_rec_to_MC, 20, 1.0, kGreen+2, 1, 0);
        # h1_mc_gen.Scale(1,"width");
        # h1_mc_gen.Scale(1/self.nev_gen);
        # h1_mc_gen.Scale(1/self.nch_gen);
        # make_common_style(h1_mc_gen, 20, 1.0, kMagenta+2, 1, 0);
        
        outname = os.path.join(self.folder, "{0}_pp_13.6TeV_{1}_material_budget_efficiency_gen_histo.root".format(date, self.period, ))

        outfile = TFile(outname,"RECREATE");
        outlist = THashList();
        outlist.SetName("outlist");
        arr_gen = np.arange(91).astype(np.float64)
        npt = len(arr_gen);
        h1amplitude     = TH1F("h1gen","gen",npt-1, arr_gen);
        # for i in range(0, npt-1):
        #     value = h1_mc_gen.GetBinContent(i)
        #     error = h1_mc_gen.GetBinError(i)
        #     h1amplitude.SetBinContent(i, value)
        #     h1amplitude.SetBinError(i, error)
        # outlist.Add(h1amplitude);

        h1recToMC     = TH1F("h1rectToMC","rec",npt-1, arr_gen);
        # for i in range(0, npt-1):
        #     value = h1_mc_rec_to_MC.GetBinContent(i)
        #     error = h1_mc_rec_to_MC.GetBinError(i)
        #     h1recToMC.SetBinContent(i, value)
        #     h1recToMC.SetBinError(i, error)
        # outlist.Add(h1recToMC);

        outfile.WriteTObject(outlist);
        outlist.Clear();

        hs_data = self.list_v0.Get("hRadius").Clone("hs_data") #FindObject("hs_conv_point").Clone("hs_data");
        h1_data = hs_data.ProjectionY("")
        #h1_data = hs_data.Projection(1 ,"");
        h1_data.SetDirectory(0);
        ROOT.SetOwnership(h1_data, False);
        h1_data.Scale(1,"width");
        h1_data.Scale(1/self.nev_data);
        h1_data.Scale(1/self.nch_data);
        make_common_style(h1_data, 20, 1.0, kBlue+1, 1, 0);    

    #_________________________________________________________________________________________    

        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);
        p1.SetLogy(1);

        frame1 = p1.DrawFrame(0,2e-7, 90,4*1e-1);  #previously 1e-2 as ymax
        frame1.GetXaxis().SetTitle("conversion radius #it{R}_{xy} (cm)");
        frame1.GetYaxis().SetTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}_{#gamma}}{d#it{r}_{xy} d#it{#eta}} (cm)^{#minus1}");
        FrameSettings(frame1)

        # if generated == True:
        #     h1_mc_gen.Draw("E0h,same");
        # h1_mc_rec.Draw("E0h,same");
        h1_data.Draw("E0h,same");
        # h1_mc_rec_to_MC.Draw("E0h,same");

        if cuts == True:
            line2 = TLine(42,1e-9,42,1e-1);
            line2.SetLineColor(kMagenta);
            line2.SetLineStyle(2);
            line2.SetLineWidth(2);
            line2.Draw("");
            ROOT.SetOwnership(line2,False);

            line3 = TLine(58,1e-9,58,1e-1);
            line3.SetLineColor(kMagenta);
            line3.SetLineStyle(2);
            line3.SetLineWidth(2);
            line3.Draw("");
            ROOT.SetOwnership(line3,False);

            line7 = TLine(14,1e-9,14,1e-1);
            line7.SetLineColor(kMagenta);
            line7.SetLineStyle(2);
            line7.SetLineWidth(2);
            line7.Draw("");
            ROOT.SetOwnership(line7,False);

            line8 = TLine(30,1e-9,30,1e-1);
            line8.SetLineColor(kMagenta);
            line8.SetLineStyle(2);
            line8.SetLineWidth(2);
            line8.Draw("");
            ROOT.SetOwnership(line8,False);

            line10 = TLine(69,1e-9,69,1e-1);
            line10.SetLineColor(kMagenta);
            line10.SetLineStyle(2);
            line10.SetLineWidth(2);
            line10.Draw("");
            ROOT.SetOwnership(line10,False);

        txt = TPaveText(0.,0.9,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);
        txt.AddText("Converted / Reconstructed photons as a function of #it{R_{xy}}, |#it{#eta}_{#gamma}| < 0.9 ");
        #txt.AddText("");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        ALICEtext("thesis")
        leg = TLegend(0.17,0.67,0.35,0.82);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.045);
        # if generated == True:
        #     leg.AddEntry(h1_mc_gen ,"M.C. gen. (LHC23d1k)","LP");
        # leg.AddEntry(h1_mc_rec ,"M.C. rec. primary #gamma (LHC23d1k)","LP");
        # leg.AddEntry(h1_mc_rec_to_MC, "M.C. rec. (MCTruth conversion point)", "LP")
        leg.AddEntry(h1_data   ,"Data #gamma candidates (LHC22o)","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        p2 = c1.cd(2);
        p2.SetPad(0,0,1,0.3);
        p2.SetMargin(0.15,0.02,0.22,0.0);
        p2.SetTicks(1,1);
        p2.SetLogy();

        frame2 = p2.DrawFrame(0,1e-3,90,20);
        frame2.GetXaxis().SetTitle("#it{R}_{xy} (cm)");
        frame2.GetYaxis().SetTitle("ratio");
        FrameSettingsRatio(frame2)

        line1 = TLine(0,1,90,1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

    # 5% lines:
        line2 = TLine(0,1.05,90,1.05);
        line2.SetLineColor(kBlack);
        line2.SetLineStyle(2);
        line2.SetLineWidth(2);
        line2.Draw("");
        ROOT.SetOwnership(line2,False); 

        line2 = TLine(0,0.95,90,0.95);
        line2.SetLineColor(kBlack);
        line2.SetLineStyle(2);
        line2.SetLineWidth(2);
        line2.Draw("");
        ROOT.SetOwnership(line2,False);  

        if cuts == True:
            line4 = TLine(42,0,42,4);
            line4.SetLineColor(kMagenta);
            line4.SetLineStyle(2);
            line4.SetLineWidth(2);
            line4.Draw("");
            ROOT.SetOwnership(line4,False);

            line5 = TLine(58,0,58,4);
            line5.SetLineColor(kMagenta);
            line5.SetLineStyle(2);
            line5.SetLineWidth(2);
            line5.Draw("");
            ROOT.SetOwnership(line5,False);

            line6 = TLine(14,0,14,4);
            line6.SetLineColor(kMagenta);
            line6.SetLineStyle(2);
            line6.SetLineWidth(2);
            line6.Draw("");
            ROOT.SetOwnership(line6,False);

            line9 = TLine(30,0,30,4);
            line9.SetLineColor(kMagenta);
            line9.SetLineStyle(2);
            line9.SetLineWidth(2);
            line9.Draw("");
            ROOT.SetOwnership(line9,False);

            line11 = TLine(69,0,69,4);
            line11.SetLineColor(kMagenta);
            line11.SetLineStyle(2);
            line11.SetLineWidth(2);
            line11.Draw("");
            ROOT.SetOwnership(line11,False);

        # h1ratio = h1_data.Clone("h1ratio");
        # make_common_style(h1ratio, 20, 1.0, kRed+1, 1, 0);
        # h1ratio.Reset();
        # h1ratio.Sumw2();
        # h1ratio.Divide(h1_data, h1_mc_rec, 1., 1., "G");
        # h1ratio.Draw("E0h,same");
        # h1ratio.SetDirectory(0);
        # ROOT.SetOwnership(h1ratio,False);

        # h1ratio1 = h1_mc_gen.Clone("h1ratio1");
        # make_common_style(h1ratio1, 20, 1.0, kGreen+2, 1, 0);
        # h1ratio1.Reset();
        # h1ratio1.Sumw2();
        # h1ratio1.Divide(h1_mc_rec_to_MC,h1_mc_gen, 1., 1., "G");
        # h1ratio1.Draw("E0h,same");
        # h1ratio1.SetDirectory(0);
        # ROOT.SetOwnership(h1ratio1,False);

        # leg = RatioLegendSettings()
        # leg.AddEntry(h1ratio   ,"Data / M.C. rec.","LP");
        # leg.AddEntry(h1ratio1 ,"M.C. rec / M.C. gen.","LP");
        # leg.Draw("");
        # ROOT.SetOwnership(leg,False);

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        if cuts == True:
            self.suffix += "_with_cuts"
        if generated == True:
            self.suffix += "_with_generated"
        filepath = os.path.join(self.folder, "{0}_pp_13.6TeV_{1}_material_budget_efficiency.pdf".format(date, self.period))
        c1.SaveAs(filepath);

        self.rootfile_data.Close();
        self.rootfile_mc  .Close();
        c1.Close();

        if cuts == True:
            self.suffix = self.suffix.replace("_with_cuts", "");
        if generated == True:
            self.suffix = self.suffix.replace("_with_generated", "");

    #_____________________________________________________________________
    
# if __name__ == "__main__":
#     cutname = "qc"
#     period_mc = "LHC23d1k";
#     period_data = "LHC22f"
#     suffix = "AnyTrack";
#     filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_155278_LHC23d1k.root"
#     filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_155756_LHC22f_pass4.root"
#     config_file = "/Users/alicamarieenderich/202312_material_budget_code/config_pp_13.6TeV_LHC22f_material.yml"
#     with open(config_file, "r", encoding="utf-8") as config_yml:
#         config = yaml.safe_load(config_yml)
#     date = "this_thesis" # datetime.date.today().strftime("%Y%m%d");
#     folder = "/Users/alicamarieenderich/{0}_efficiency_plots/".format(date);  
#     os.makedirs(folder, exist_ok=True);
#     cuts = False
#     generated = True
#     draw = draw_1D_efficiency(filename_data, filename_mc, cutname, folder, period_data, suffix)
#     draw.draw_eff(cuts, generated, date);

if __name__ == "__main__":
    filename_data = "/Users/juliaschlagel/Analysis240411/Analysis/HLtrains/data/AnalysisResults_LHC22o_full_statistics.root"
    filename_mc = "/Users/juliaschlagel/Analysis240411/Analysis/HLtrains/MC/AnalysisResults_LHC24b1_full_statistics.root"
    cutname = "qc"
    period_data= "LHC22o"
    period_mc = "LHC24b1"
    suffix = "AnyTrack" 
    folder = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/cuts"
    os.makedirs(folder, exist_ok=True);
    date= "this_thesis"
    generated = False
    cuts = False
    draw = draw_1D_efficiency(filename_data, filename_mc, folder, period_data, suffix)
    draw.draw_eff(cuts, generated, date)