# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3

import os, sys, shutil
import math
import argparse
import pickle
import copy
import numpy as np
import ctypes
import yaml
from copy import deepcopy
import ROOT
import datetime
ROOT.gROOT.SetBatch(True);
from ROOT import gStyle
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullCross, kFullDiamond, kOpenSquare, kOpenTriangleUp, kOpenCircle, kFullCrossX
from ROOT import TFile, THashList, TF1, TString, TCanvas, TLegend, TPaveText
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
from FitInvMassForPt import PairAnalyzer
from PlotInvMass import PlotInvMass
from PlotParameterHistoInvMass import PlotHistoInvMass
from PlotParametersCombined import PlotHistoParametersCombined
from PlotRawYield import PlotRawYieldInvMass
from PlotSameAndMixed import PlotSameMixed
from acceptance import HistogramAcceptance
from efficiency import Efficiency
from corrected_yield import CorrectedYield
from comparison import Compare
from utility import Utility
from store_histograms import Save_Hist
#_________________________________________________________________________________________
def run(start_pt, filename, config, typ, suffix, folder):
    
    print("##############", decay, "###########")
    print(sys._getframe().f_code.co_name);
    arr_pt = np.array(config["common"]["pt_bin"],dtype=float);
    print("pT binning = ",arr_pt);
    print("type = ",typ);
    print("input = ",filename);
    rootfile = TFile.Open(filename,"READ");
    meson = config["common"]["meson"];

    list_fit_func = config[typ]["fit_func"];

    list_fit_min = config["common"]["fit_min"];
    list_fit_max = config["common"]["fit_max"];
    if len(list_fit_min) != len(list_fit_max):
        return;

    list_integral_min = config["common"]["integral_min"];
    print("Integral min: ", list_integral_min)
    list_integral_max = config["common"]["integral_max"];
    print("Integral max: ", list_integral_max)
    if len(list_integral_min) != len(list_integral_max):
        return;

    list_yield_min = config["common"]["yield_min"];
    list_yield_max = config["common"]["yield_max"];
    if len(list_yield_min) != len(list_yield_max):
        return;
    subsystems = config[typ]["subsystems"]
    
    nsys = len(config[typ]['subsystems']);
    print("Number of subsystems: ", nsys); 

    # start_pt = config["common"]["start_pT_bin"]
    # print("start_pT: ", start_pt)
    

    meson = config["common"]["meson"];

    if config["common"]["do_ptspectrum"] == True:
        if run_type == "one":
            if decay =="dalitz":
                outname = os.path.join(folder,"inv_mass_analysis", "Dalitz","{0}_{1}_{2}_{3}_{4}_ptspectrum_{5}_{6}TeV_{7}_{8}.root".format(date, decay, period, meson, typ, config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix))
            else:
                outname = os.path.join(folder,"inv_mass_analysis", "PCM","{0}_{1}_{2}_{3}_{4}_ptspectrum_{5}_{6}TeV_{7}_{8}.root".format(date, decay, period, meson, typ, config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix));
            print("output file name = ",outname);
            outfile = TFile(outname,"RECREATE");
            
            for isys in range(0,nsys):
                if typ == "data":
                    ssname = config[typ]['subsystems'][isys]['name']; #subsystem name
                if typ == "mc":
                    ssname = config[typ]['subsystems'][isys]['mc_name']
                if typ == "mc":
                    if decay == "dalitz":
                        ana_pi0 = PairAnalyzer(meson, filename, ssname, typ, filename_mc, config, decay);
                    if decay == "pcm":
                        ana_pi0 = PairAnalyzer(meson, filename, ssname, typ, filename_mc, config, decay)
                else:
                    if decay == "dalitz":
                        ana_pi0 = PairAnalyzer(meson, filename, ssname, typ, filename_mc, config, decay);
                    if decay == "pcm":
                        ana_pi0 = PairAnalyzer(meson, filename, ssname, typ, filename_mc, config, decay);
                ana_pi0.set_arr_pt(arr_pt)
                #dirname_ss_mc = "pi0eta-to-gammagamma-mc-pcmdalitzee"
                ana_pi0.set_subsystem(ssname);
                ana_pi0.set_xtitle("#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})");
                ana_pi0.set_ytitle("#it{p}_{T,#gamma#gamma} (GeV/#it{c}^{2})");
                print("analyze subsystem", ssname);
                
                fit_parameters = subsystems[isys].get('fit_parameters', [])
                fit_limit_min = subsystems[isys].get('fit_limit_min', [])
                fit_limit_max = subsystems[isys].get('fit_limit_max', [])
                ana_pi0.set_fit_params(fit_parameters, fit_limit_min, fit_limit_max)
                
                outlist_ss = THashList()
                outlist_ss.SetName(ssname);
                outlist_ss.SetOwner(True);
                
                
                for ifunc in list_fit_func:
                    ana_pi0.set_fit_function(typ, ifunc);
                    outlist_func = THashList();
                    outlist_func.SetName(ifunc);
                    outlist_func.SetOwner(True);
                    outlist_ss.Add(outlist_func);
                    for ir in range(0, len(list_fit_min)):
                        fit_min = list_fit_min[ir];
                        fit_max = list_fit_max[ir];
                        ana_pi0.set_fit_range(fit_min, fit_max);
                        integral_min = list_integral_min[0];
                        integral_max = list_integral_max[0];
                        ana_pi0.set_integral_range(integral_min, integral_max);
                        yield_min = list_yield_min[0];
                        yield_max = list_yield_max[0];
                        ana_pi0.set_yield_range(yield_min, yield_max);
                        if start_pt == 0:
                            outlist_fit_range_zer = ana_pi0.analyze_ptspectrum(0)
                            outlist_fit_range_zer.SetName("fit_{0:3.2f}_{1:3.2f}_GeVc2".format(fit_min, fit_max))
                            outlist_func.Add(outlist_fit_range_zer)
                        else:
                            outlist_fit_range = ana_pi0.analyze_ptspectrum(start_pt);
                            outlist_fit_range.SetName("fit_{0:3.2f}_{1:3.2f}_GeVc2".format(fit_min, fit_max));
                            outlist_func.Add(outlist_fit_range);
                outfile.WriteTObject(outlist_ss);
                outlist_ss.Clear();
                print("Fitinvmass script beendet")
            del ana_pi0;

            #plot_pi0 = PlotInvMass(meson, outname, "pi0eta-to-gammagamma");
            #plot_pi0.set_arr_pt(arr_pt);
            yield_list = []
            ssname_list = []
            list_parameters_comparison = []
            for isys in range(0,nsys):
                if typ == "data":
                    ssname = config[typ]['subsystems'][isys]['name']; #subsystem name
                if typ == "mc":
                    ssname = config[typ]['subsystems'][isys]['mc_name']
                plot_pi0 = PlotInvMass(meson, outname, ssname)
                plot_pi0.set_arr_pt(arr_pt);
                #ssname = config[typ]['subsystems'][isys]['name']; #subsystem name
                ssname_list.append(ssname)
                plot_pi0.set_subsystem(ssname);
                print("plot subsystem", ssname)
                nfit = len(list_fit_func)
                # list_parameters_comparison = [];
                # yield_list = [];           
                for ifit in range(0, nfit): 
                    fitname = list_fit_func[ifit];
                    print("fitname ", fitname)
                    plot_pi0.set_fitname(fitname);
                    for ir in range(len(list_fit_min)):
                        fit_min = list_fit_min[ir];
                        fit_max = list_fit_max[ir];
                        plot_pi0.set_fit_range(fit_min, fit_max);
                        yield_min = list_yield_min[0];
                        yield_max = list_yield_max[0];
                        plot_pi0.set_yield_range(yield_min, yield_max);
                        if meson == "pi0":
                            plottingRange = [0., 0.3];
                        if meson == "eta":
                            plottingRange = [0.3, 0.8]
                        list_plot = plot_pi0.set_fit_list();
                        yield_list.append(list_plot.FindObject("h1yield_param"))
                        function_list = [];
                        histogram_list = [];
                        parameter_list = [];
                        same_list = [];
                        mixed_list = [];
                        it = ROOT.TIter(list_plot)
                        obj = it.Next()
                        while obj:
                            objName = obj.GetName()
                            if "param" in objName:
                                parameter_list.append(obj)
                            obj = it.Next()

            # pdf output of all parameters for each combination of cut and fit                 
                        # if typ == "mc":
                        #     plot_histo = PlotHistoInvMass(meson, outname, "pi0eta-to-gammagamma-mc");
                        # else:
                        #     plot_histo = PlotHistoInvMass(meson, outname, "pi0eta-to-gammagamma");  
                        
                        plot_histo = PlotHistoInvMass(meson, outname, ssname)
                        if decay=="dalitz":
                            outname_histo = os.path.join(folder, "inv_mass_analysis","Dalitz","{0}_{1}_{2}_InvMass_Fitparameters_{3}.pdf".format(date, decay, period, fitname))
                        else:
                            outname_histo = os.path.join(folder, "inv_mass_analysis","PCM","{0}_{1}_{2}_InvMass_Fitparameters_{3}.pdf".format(date, decay, period, fitname))
                        plot_histo.PlotHistoParameters(parameter_list, TString(outname_histo), period, 2,3,"#gamma#gamma", "{0}".format(typ))

            # pdf output of all fitted histograms for each combination of cut and fit
                        for ipt in range(start_pt, len(arr_pt)):
                            histo_ipt =list_plot.FindObject("h1mgg_pt{0}".format(ipt));
                            histogram_list.append(histo_ipt)
                            #print("##Histogram list: ", histogram_list)
                            function_ipt = list_plot.FindObject("f1total_pt{0}".format(ipt));
                            function_list.append(function_ipt)
                            #print("###Function list: ", function_list)
                        if decay =="dalitz":
                            output_name_with_yield = os.path.join(folder,"inv_mass_analysis", "Dalitz","{0}_{1}_{2}_InvMass_Overview_{3}_with_yield.pdf".format(date, decay, period, fitname))
                            #print("outname", output_name_with_yield)
                        else:
                            output_name_with_yield = os.path.join(folder,"inv_mass_analysis","PCM", "{0}_{1}_{2}_InvMass_Overview_{3}_with_yield.pdf".format(date, decay, period, fitname))
                        if typ == "mc":
                            utils = Utility(meson, filename_mc, config, decay, "mc")
                        if typ == "data":
                            utils = Utility(meson, filename_data, config, decay, "data")
                        #nev = utils.get_Nev()
                        #print("###Nev : ", nev)
                        plot_pi0.PlotInvMassInPtBins(utils, start_pt, histogram_list, function_list, TString(output_name_with_yield), "", "", 
                                                    plottingRange, period, 5, 7, 0, len(arr_pt),len(arr_pt), "#gamma#gamma", 
                                                    typ, "yield", decay) #for pi0 binning: 5, 6 (for eta binning: 3, 4)
                                
            # pdf output of all same and mixed scaled histograms for each cut
                        if typ == "data":
                            plot_same_mixed = PlotSameMixed(meson, outname, ssname)
                            #plot_same_mixed = PlotSameMixed(meson, outname,"pi0eta-to-gammagamma" )
                            plot_same_mixed.set_arr_pt(arr_pt)
                            plot_same_mixed.set_fit_range(fit_min, fit_max)
                            plot_same_mixed.set_integral_range(integral_min, integral_max)
                            for ipt in range(start_pt,len(arr_pt)):
                                histo_ipt_same =list_plot.FindObject("h1mgg_same_pt{0}".format(ipt));
                                same_list.append(histo_ipt_same)
                                histo_ipt_mixed = list_plot.FindObject("h1mgg_mix_scaled_pt{0}".format(ipt));
                                mixed_list.append(histo_ipt_mixed)
                            if decay =="dalitz":
                                output_name = os.path.join(folder, "inv_mass_analysis","Dalitz", "{0}_{1}_{2}_InvMass_Scaled_Mixed_{3}.pdf".format(date, decay, period, fitname))
                            else:
                                output_name = os.path.join(folder, "inv_mass_analysis","PCM", "{0}_{1}_{2}_InvMass_Scaled_Mixed_{3}.pdf".format(date, decay, period, fitname))
                            plot_same_mixed.PlotSameMixedInPtBins(start_pt, same_list, mixed_list, TString(output_name), "", "", plottingRange, period, 5, 7, 0, 
                                                                len(arr_pt),len(arr_pt), "#gamma#gamma", "{0}".format(typ), decay) #for pi0 binning: 5, 6 (for eta: 3,4)

            # pdf output of mass, amplitude and width for each fit and comparison of all cuts
                        parameter_comparison = [];
                        for i_parameter in range(5):
                            parameter_comparison.append(parameter_list[i_parameter]);
                        list_parameters_comparison.append(parameter_comparison);

            print("ssname_list: ", ssname_list)
            plot_raw_yield = PlotRawYieldInvMass(decay, meson, outname, "pi0eta-to-gammagamma") 
            if decay == "dalitz":    
                outname_raw_yield = os.path.join(folder, "inv_mass_analysis", "Dalitz","{0}_{1}_{2}_InvMass_RawYield_{3}.pdf".format(date, decay, period, fitname))
            else:
                outname_raw_yield = os.path.join(folder, "inv_mass_analysis", "PCM","{0}_{1}_{2}_InvMass_RawYield_{3}.pdf".format(date, decay, period, fitname));  
            plot_raw_yield.PlotHistoYield(yield_list, TString(outname_raw_yield), "", "", period, 1,1, "#gamma#gamma", "data", ssname_list)                   
            plot_parameters = PlotHistoParametersCombined(meson, outname,"pi0eta-to-gammagamma")
            if decay =="dalitz":
                outname_histo_linear = os.path.join(folder, "inv_mass_analysis", "Dalitz", "{0}_{1}_{2}_InvMass_Parameters_Combined_{3}_{4}.pdf".format(date, decay, period, fitname, ssname));
            else:
                outname_histo_linear = os.path.join(folder, "inv_mass_analysis", "PCM","{0}_{1}_{2}_InvMass_Parameters_Combined_{3}_{4}.pdf".format(date, decay, period, fitname, ssname));
            print("##################### länge histo param: ", len(parameter_comparison))
            plot_parameters.PlotHistoParametersCombined(list_parameters_comparison, TString(outname_histo_linear), "", "", period, 2, 3, 
                                                        "#gamma#gamma", "{0}".format(typ), ssname_list)
                    
            for iss in range (len(yield_list)):
                if typ == "mc":
                    yield_used = yield_list[iss]
                    raw_mc = yield_used
                    print("#### Raw yield mc:")
                    for i in range(1, raw_mc.GetNbinsX() + 1):
                        bin_content_m = raw_mc.GetBinContent(i)
                        bin_center_m = raw_mc.GetBinCenter(i)
                        print(f"Bin {i}: Center = {bin_center_m}, Content = {bin_content_m}")
                    
                    filename_mc_pcm = f'th1d_raw_mc_pcm_{iss}.pkl'
                    filename_mc_dalitz = f"th1d_raw_mc_dalitz{iss}.pkl"

                    if decay == "pcm":
                        with open(filename_mc_pcm, 'wb') as raw_mc_file_pcm:
                            pickle.dump(raw_mc, raw_mc_file_pcm)
                    else:
                        with open(filename_mc_dalitz, 'wb') as raw_mc_file_dalitz:
                            pickle.dump(raw_mc, raw_mc_file_dalitz)

                elif typ == "data":
                    yield_used = yield_list[iss]
                    raw_data = yield_used
                    print("#### Raw yield data:")
                    for i in range(1, raw_data.GetNbinsX() + 1):
                        bin_content_d = raw_data.GetBinContent(i)
                        bin_center_d = raw_data.GetBinCenter(i)
                        print(f"Bin {i}: Center = {bin_center_d}, Content = {bin_content_d}")
                    
                    filename_data_pcm = f'th1d_raw_data_pcm_{iss}.pkl'
                    filename_data_dalitz = f"th1d_raw_data_dalitz{iss}.pkl"
                    if decay == "pcm":
                        with open(filename_data_pcm, 'wb') as raw_data_file_pcm:
                            pickle.dump(raw_data, raw_data_file_pcm)
                    else:
                        with open(filename_data_dalitz, 'wb') as raw_data_file_dalitz:
                            pickle.dump(raw_data, raw_data_file_dalitz)
            
            del plot_pi0
            del plot_histo
            del parameter_list
            del histogram_list
            del function_list
            outfile.Close()
                
        if run_type == "both_types":
            print("run type: ", run_type)
            eff_acc_list_dalitz = []
            eff_acc_list_pcm = []
            corr_list_dalitz = []
            corr_list_pcm_pi0 = []
            corr_list_pcm_eta = []

            utils_mc = Utility(meson, filename_mc, config, decay, "mc")
            utils_data = Utility(meson, filename_data, config, decay, "data")

            for isys in range(0, nsys):
                ssname = config[typ]["subsystems"][isys]["name"]
                print("Analyse Subsystem: ", ssname)
                print("Analyse SS for decay: ", decay,)
                if decay == "dalitz":
                    filename_rootfile = os.path.join(folder,"inv_yield_calculation", "Dalitz","{0}_{1}_results.root".format(date, ssname))
                if decay == "pcm":
                    filename_rootfile = os.path.join(folder, "inv_yield_calculation", "PCM","{0}_{1}results.root".format(date, ssname))
            
                # save_hist_dal = Save_Hist(filename_rootfile_dal)
                save_hist = Save_Hist(filename_rootfile)
            
                filename_mc_pcm = f'th1d_raw_mc_pcm_{isys}.pkl'
                filename_data_pcm = f'th1d_raw_data_pcm_{isys}.pkl'
                if os.path.exists(filename_data_pcm):
                    with open(filename_data_pcm, 'rb') as raw_data_file_pcm:
                        raw_data_pcm = pickle.load(raw_data_file_pcm)
                        #if decay == "pcm":
                        #outlist_ss_both.Add(raw_data_pcm)
                        if decay == "pcm":
                            save_hist.add_histogram("PCMPCM", raw_data_pcm, "h1_yield_data")
                            save_hist.save_histograms()
                        print("####raw data file pcm geladen: ", raw_data_pcm)
                        # for i in range(1, raw_data_pcm.GetNbinsX() + 1):
                        #     bin_content_d = raw_data_pcm.GetBinContent(i)
                        #     bin_center_d = raw_data_pcm.GetBinCenter(i)
                        #     print(f"Bin {i}: Center = {bin_center_d}, Content = {bin_content_d}")
                if os.path.exists(filename_mc_pcm):
                    with open(filename_mc_pcm, 'rb') as raw_mc_file_pcm:
                        raw_mc_pcm = pickle.load(raw_mc_file_pcm)
                        #outlist_ss_both.Add(raw_data_pcm)
                        if decay == "pcm":
                            save_hist.add_histogram("PCMPCM", raw_mc_pcm, "h1_yield_mc")
                            save_hist.save_histograms()
                        print("####raw mc file pcm geladen: ", raw_mc_pcm)
                        print("typ von raw mc: ", type(raw_mc_pcm))
                        print("####raw mc file pcm geladen: ", raw_data_pcm)
                        # for i in range(1, raw_mc_pcm.GetNbinsX() + 1):
                        #     bin_content_d = raw_mc_pcm.GetBinContent(i)
                        #     bin_center_d = raw_mc_pcm.GetBinCenter(i)
                        #     print(f"Bin {i}: Center = {bin_center_d}, Content = {bin_content_d}")
                    
                filename_mc_dalitz = f"th1d_raw_mc_dalitz{isys}.pkl"
                filename_data_dalitz = f"th1d_raw_data_dalitz{isys}.pkl"
                if os.path.exists(filename_data_dalitz):
                    with open(filename_data_dalitz, 'rb') as raw_data_file_dalitz:
                        raw_data_dalitz = pickle.load(raw_data_file_dalitz)
                        if decay == "dalitz":
                            save_hist.add_histogram("PCMDalitzEE", raw_data_dalitz,"h1_yield_data")
                            save_hist.save_histograms()
                        print("####raw data file dalitz geladen: ", raw_data_dalitz)
                        # for i in range(1, raw_data_dalitz.GetNbinsX() + 1):
                        #     bin_content_d = raw_data_dalitz.GetBinContent(i)
                        #     bin_center_d = raw_data_dalitz.GetBinCenter(i)
                            # print(f"Bin {i}: Center = {bin_center_d}, Content = {bin_content_d}")
                if os.path.exists(filename_mc_dalitz):
                    with open(filename_mc_dalitz, 'rb') as raw_mc_file_dalitz:
                        raw_mc_dalitz = pickle.load(raw_mc_file_dalitz)
                        # outlist_ss_both.Add(raw_data_pcm)
                        if decay == "dalitz":
                            save_hist.add_histogram("PCMDalitzEE", raw_mc_dalitz, "h1_yield_mc")
                            save_hist.save_histograms()

                        print("####raw mc file dalitz geladen: ", raw_mc_dalitz)
                        print("typ von raw mc: ", type(raw_mc_dalitz))
                        # for i in range(1, raw_mc_dalitz.GetNbinsX() + 1):
                        #     bin_content_d = raw_mc_dalitz.GetBinContent(i)
                        #     bin_center_d = raw_mc_dalitz.GetBinCenter(i)
                        #     print(f"Bin {i}: Center = {bin_center_d}, Content = {bin_content_d}")
                # outfile.WriteTObject(outlist_ss_both);
                # outlist_ss_both.Clear();
                #outfile.Close()

                if decay == "pcm":
                    raw_mc = raw_mc_pcm
                    raw_data = raw_data_pcm
                else:
                    raw_mc = raw_mc_dalitz
                    raw_data = raw_data_dalitz
                
                        
                eff_class = Efficiency(decay, meson)
                if decay == "dalitz":
                    eff_acc_name = os.path.join(folder, "inv_yield_calculation", "Dalitz", "Efficiency_Acceptance_{0}".format(ssname))
                    eff_acc_name_scaled = os.path.join(folder, "inv_yield_calculation", "Dalitz", "Efficiency_Acceptance_scaled_by_BR_{0}".format(ssname))
                    Teff_acc_name = os.path.join(folder, "inv_yield_calculation", "Dalitz", "T_Efficiency_Acceptance_shifted_{0}".format(ssname))
                else:
                    eff_acc_name = os.path.join(folder, "inv_yield_calculation", "PCM", "Efficiency_Acceptance_{0}".format(ssname))
                    eff_acc_name_scaled = os.path.join(folder, "inv_yield_calculation", "PCM", "Efficiency_Acceptance_scaled_by_BR_{0}".format(ssname))
                    Teff_acc_name = os.path.join(folder, "inv_yield_calculation", "PCM", "T_Efficiency_Acceptance_shifted_{0}".format(ssname))
                eff_acc = eff_class.eff_and_acc_vs_pt("false", raw_mc, utils_mc, TString(eff_acc_name), "true")
                eff_acc_scaled = eff_class.eff_and_acc_vs_pt("true", raw_mc, utils_mc, TString(eff_acc_name_scaled), "true")
                teff_acc = eff_class.teff_and_acc_vs_pt(raw_mc, utils_mc, TString(Teff_acc_name), "true")
                
                for i in range(1, eff_acc.GetNbinsX() + 1):
                    bin_content_d = eff_acc.GetBinContent(i)
                    bin_center_d = eff_acc.GetBinCenter(i)
                    print(f"Bin {i}: Center = {bin_center_d}, Content = {bin_content_d}")
                if decay == "pcm":
                    eff_acc_list_pcm.append(copy.deepcopy(eff_acc_scaled))
                    save_hist.add_histogram("PCMPCM", eff_acc)
                    save_hist.save_histograms()
                if decay == "dalitz":
                    eff_acc_list_dalitz.append(copy.deepcopy(eff_acc_scaled))
                    save_hist.add_histogram("PCMDalitzEE", eff_acc, "h_eff_acc")
                    save_hist.save_histograms()
                print("List length eff and acc: ", len(eff_acc_list_dalitz), eff_acc_list_dalitz, "++++++++++++")


                corr_class = CorrectedYield(filename_mc, decay, meson)
                if decay == "dalitz":
                    pythia_name = os.path.join(folder, "inv_yield_calculation", "Dalitz", "{0}_pythia".format(date))
                    corr_name= os.path.join(folder, "inv_yield_calculation", "Dalitz", "{0}_corrected_yield".format(date))
                    corr_name_scaled = os.path.join(folder, "inv_yield_calculation", "Dalitz", "{0}_corrected_yield_scaled_by_pT".format(date))
                    comp_pythia_corr_name = os.path.join(folder, "inv_yield_calculation",  "Dalitz", "{0}_comp_pythia_corr_yield".format(date))
                else:                    
                    pythia_name = os.path.join(folder, "inv_yield_calculation", "PCM", "{0}_pythia".format(date))
                    corr_name = os.path.join(folder, "inv_yield_calculation", "PCM", "{0}_corrected_yield".format(date))
                    corr_name_scaled = os.path.join(folder, "inv_yield_calculation", "PCM", "{0}_corrected_yield_scaled_by_pT".format(date))
                    comp_pythia_corr_name = os.path.join(folder, "inv_yield_calculation", "PCM", "{0}_comp_pythia_corr_yield".format(date))
                    
                pyth = corr_class.plot_pythia(utils_mc, TString(pythia_name), "true")
                if decay == "pcm":
                    corr = corr_class.calc_corr_yield_new(utils_data, raw_data_pcm, eff_acc, decay, "false",TString(corr_name), "true")
                    corr_scaled = corr_class.calc_corr_yield_new(utils_data, raw_data_pcm, eff_acc, decay, "true", TString(corr_name_scaled), "true")
                else:
                    corr =corr_class.calc_corr_yield_new(utils_data, raw_data_dalitz, eff_acc, decay, "false", TString(corr_name), "true")
                    corr_scaled = corr_class.calc_corr_yield_new(utils_data, raw_data_dalitz, eff_acc, decay, "true", TString(corr_name_scaled), "true")

                if decay == "dalitz":
                    print("######corrected yield dalitz:")
                else:
                    print("######corrected yield pcm:")
                for i in range(1, corr.GetNbinsX() + 1):
                        bin_content_d = raw_mc_dalitz.GetBinContent(i)
                        bin_center_d = raw_mc_dalitz.GetBinCenter(i)
                        print(f"Bin {i}: Center = {bin_center_d}, Content = {bin_content_d}")

                corr_class.compare_pythia_corr_yield(corr, pyth, decay,TString(comp_pythia_corr_name), "true")
                if decay == "pcm":
                    fit_name = os.path.join(folder, "inv_yield_calculation", "PCM", "TCM_fit")
                if decay == "dalitz":
                    fit_name = os.path.join(folder, "inv_yield_calculation", "Dalitz", "TCM_fit")
                
                corr_class.plot_fit_func(corr_scaled, decay, TString(fit_name), "true")
                if decay == "pcm":
                    fit_13Tev_pcm = os.path.join(folder, "inv_yield_calculation", "PCM", "TCM_13_Tev")
                    if meson == "pi0":
                        corr_class.plot_fit_func_13Tev(filename_mc, filename_comp_pcm, "pcm", TString(fit_13Tev_pcm), "true")
                if decay == "dalitz":
                    fit_13Tev_dalitz = os.path.join(folder, "inv_yield_calculation", "Dalitz", "TCM_13_Tev")
                    if meson == "pi0":
                        corr_class.plot_fit_func_13Tev(filename_mc, filename_comp_dalitz, "dalitz", TString(fit_13Tev_dalitz), "true")

                if decay == "pcm":
                    if meson == "pi0":
                        corr_list_pcm_pi0.append(copy.deepcopy(corr))
                    if meson == "eta":
                        corr_list_pcm_eta.append(copy.deepcopy(corr))
                    save_hist.add_histogram("PCMPCM", pyth)
                    save_hist.save_histograms()
                    save_hist.add_histogram("PCMPCM", corr)
                    save_hist.save_histograms()
                    print("meson: ", meson)
                    print("List corr_yield_pcm_pi0: ", corr_list_pcm_pi0)
                    print("list corr_yield_pcm_eta: ", corr_list_pcm_eta)

                if decay == "dalitz":
                    corr_list_dalitz.append(copy.deepcopy(corr))
                    save_hist.add_histogram("PCMDalitzEE", pyth, "pythia")
                    save_hist.save_histograms()
                    save_hist.add_histogram("PCMDalitzEE", corr, "h_corr")
                    save_hist.save_histograms()
                    print("###################Typ des in der Liste gespeicherten Histogramms: ", type(corr))
                    print("Liste corr_yield_dalitz: ", corr_list_dalitz)
                    
            
                if decay == "dalitz":
                    print("#######Länge Liste corr_dalitz", len(corr_list_dalitz),corr_list_dalitz)
                    
                    corr_yield_dal = corr_list_dalitz[isys]
                    filename_corr_y_dal = f"corr_dalitz_{isys}"
                    with open(filename_corr_y_dal,"wb") as corr_yield_dalitz_file:
                        pickle.dump(corr_yield_dal, corr_yield_dalitz_file)
                if decay =="pcm":
                    if meson == "pi0":
                        print("#######Länge Liste corr_pcm pi0: ", len(corr_list_pcm_pi0),corr_list_pcm_pi0)
                        corr_yield_pcm_pi0 = corr_list_pcm_pi0[isys]
                        filename_corr_y_pcm_pi0 = f"corr_pcm_pi0_{isys}"
                        with open(filename_corr_y_pcm_pi0, "wb") as corr_yield_pcm_file_pi0:
                            pickle.dump(corr_yield_pcm_pi0, corr_yield_pcm_file_pi0)
                    #print("##############list corr_pcm_eta: ", len(corr_list_pcm_eta))
                    if meson == "eta":
                        print("#######Länge Liste corr_pcm eta: ", len(corr_list_pcm_eta),corr_list_pcm_eta)
                        corr_yield_pcm_eta = corr_list_pcm_eta[isys]
                        filename_corr_y_pcm_eta = f"corr_pcm_eta_{isys}"
                        with open(filename_corr_y_pcm_eta, "wb") as corr_yield_pcm_file_eta:
                            pickle.dump(corr_yield_pcm_eta, corr_yield_pcm_file_eta)
                
                        

                comp_corr_13Tev_dalitz_name = os.path.join(folder, "Comparison", "{0}_Comp_13TeV_corr_yield_dalitz".format(date))
                comp_corr_13Tev_pcm_name = os.path.join(folder, "Comparison", "{0}_Comp_13TeV_corr_yield_pcm".format(date))
                

                comp_class = Compare(filename_mc, decay, meson)
                ssname_list_dalitz = []
                ssname_list_pcm = []
                for id in range(0, len(eff_acc_list_dalitz)):
                    ssname_cut = config[typ]["subsystems"][id]["name"]
                    ssname_list_dalitz.append(ssname)
                    print(ssname_cut)
                    comp_eff_times_acc_13Tev_dalitz_name = os.path.join(folder, "Comparison", "Comp_13TeV_eff_times_acc_dalitz_{0}".format(ssname_cut))
                    if meson == "pi0":
                        comp_class.compare_13tev(filename_comp_dalitz, "eff_times_acc", eff_acc_list_dalitz[id], "dalitz", TString(comp_eff_times_acc_13Tev_dalitz_name), "true")

                for ip in range (0, len(eff_acc_list_pcm)):
                    ssname_cut = config[typ]["subsystems"][ip]["name"]
                    print(ssname_cut)
                    ssname_list_pcm.append(ssname_cut)
                    comp_eff_times_acc_13Tev_pcm_name = os.path.join(folder, "Comparison", "Comp_13TeV_eff_times_acc_pcm_{0}".format(ssname_cut))
                    if meson == "pi0":
                        comp_class.compare_13tev(filename_comp_pcm, "eff_times_acc", eff_acc_list_pcm[ip], "pcm", TString(comp_eff_times_acc_13Tev_pcm_name), "true")

                if decay == "dalitz":
                    comp_teff_eff = os.path.join(folder, "inv_yield_calculation", "Dalitz", "Comp_eff_vs_Teff")
                    
                    comp_name_cuts_eff_acc = os.path.join(folder, "inv_yield_calculation", "Dalitz", "Comp_cuts_eff_times_acc_{0}".format(date))
                    comp_class.compare_cuts(eff_acc_list_dalitz, ssname_list_dalitz, "dalitz", filename_mc, config, TString(comp_name_cuts_eff_acc), "eff_acc", utils_mc)
                if decay == "pcm":
                    comp_teff_eff = os.path.join(folder, "inv_yield_calculation", "PCM", "Comp_eff_vs_Teff")
                    comp_name_cuts_eff_acc = os.path.join(folder, "inv_yield_calculation", "PCM", "Comp_cuts_eff_times_acc_{0}".format(date))
                    comp_class.compare_cuts(eff_acc_list_pcm, ssname_list_pcm, "pcm", filename_mc, config, TString(comp_name_cuts_eff_acc), "eff_acc", utils_mc)
                utils = Utility(meson, filename, config, decay, typ)
                comp_class.comp_eff_Teff(utils, raw_mc, eff_acc_scaled, teff_acc, TString(comp_teff_eff), "true")

                

                
                if decay == "dalitz":
                    if meson == "pi0":
                        comp_class.compare_13tev(filename_comp_dalitz, "corr", corr, "dalitz", TString(comp_corr_13Tev_dalitz_name), "true")
                        comp_class.compare_13tev(filename_comp_dalitz,"eff_times_acc", eff_acc, "dalitz",TString(comp_eff_times_acc_13Tev_dalitz_name), "true")
                if decay == "pcm":
                    if meson == "pi0":
                        comp_class.compare_13tev(filename_comp_pcm, "corr", corr, "pcm", TString(comp_corr_13Tev_pcm_name), "true")
                        comp_class.compare_13tev(filename_comp_pcm, "eff_times_acc", eff_acc, "pcm", TString(comp_eff_times_acc_13Tev_pcm_name), "true")
                    
                
        if run_type == "both_decays":
            print("run type: ", run_type)
            final_dec = decay_array[(len(decay_array)-1)]
            print("final decay: ", final_dec)
            config_file_dalitz = config_array_dalitz[0]
            print("config file dal: ", config_file_dalitz)
            config_file_pcm = config_array_pcm[0]
            with open(config_file_dalitz, "r", encoding="utf-8") as config_dalitz_yml:
                config_dalitz = yaml.safe_load(config_dalitz_yml)
            with open(config_file_pcm, "r", encoding="utf-8") as config_pcm_yml:
                config_pcm = yaml.safe_load(config_pcm_yml)

            nsys_dalitz = len(config_dalitz[typ]['subsystems'])
            nsys_pcm = len(config_pcm[typ]['subsystems'])

            for isys_d in range(0,nsys_dalitz):
                filename_corr_y_dal = f"corr_dalitz_{isys_d}"
                with open(filename_corr_y_dal, "rb") as corr_dalitz_file:
                    corr_yield_dalitz = pickle.load(corr_dalitz_file)
                print("corrected_yield_dalitz geladen: ", corr_yield_dalitz, isys_d)
            for isys_p in range(0, nsys_pcm):
                print("######MESON: ", meson)
                filename_corr_y_pcm = f"corr_pcm_{meson}_{isys_p}"
                with open(filename_corr_y_pcm, "rb") as corr_pcm_file:
                    corr_yield_pcm = pickle.load(corr_pcm_file)
                print("corrected_yield_pcm loaded: ", corr_yield_pcm, isys_p)
            comp_class = Compare(filename_mc, decay, meson)
            comp_corr_name = os.path.join(folder, "Comparison", "{0}_Comp_corr_yield_dalitz_pcm".format(date))
            comp_class.compare_corr_dalitz_pcm(corr_yield_dalitz, corr_yield_pcm, TString(comp_corr_name), "true")

        if run_type == "both_mesons":
            print("run type: ", run_type)
            final_meson = meson_array[1]
            print("final meson: ", final_meson)
            if meson == final_meson:
                filename_corr_pcm_pi0 = f"corr_pcm_pi0_0"
                with open(filename_corr_pcm_pi0, "rb") as corr_pcm_pi0_file:
                    corr_yield_pcm_pi0 = pickle.load(corr_pcm_pi0_file)
                    print("corrected yield pcm pi0 loaded: ", corr_yield_pcm_pi0, type(corr_yield_pcm_pi0))
                filename_corr_pcm_eta = f"corr_pcm_eta_0"
                with open(filename_corr_pcm_eta, "rb") as corr_pcm_eta_file:
                    corr_yield_pcm_eta = pickle.load(corr_pcm_eta_file)
                    print("corrected yield pcm eta loaded: ", corr_yield_pcm_eta, type(corr_yield_pcm_eta))

                comp_class = Compare(filename_mc, decay, meson)
                comp_class_eta = Compare(filename_mc, decay, "eta")
                name_ratio = os.path.join(folder, "Comp_corr_eta_pi0")
                name_comp_13Tev = os.path.join(folder, "Comp_corr_eta_pi0_13_TeV")
                h_eta_pi0 = comp_class.corr_yield_ratio(corr_yield_pcm_eta, corr_yield_pcm_pi0, TString(name_ratio), "true")
                print(" typ von eta to pi0 ratio hist: ", type(h_eta_pi0))
                comp_class_eta.compare_13tev_etapi0(filename_comp_pcm, h_eta_pi0, TString(name_comp_13Tev), "true")


            
           

            
            
    else:
        print("please check what to do in",config)

        rootfile.Close()

    
#_________________________________________________________________________________________

########################################
#       Run over a single Dataset      #
########################################

# period = "LHC22f"
# filename = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_124838_LHC22f_pass4.root"
# type = "data"
# config_file = "configs/config_pp_13.6TeV_pi0_LHC22f.yml" #"/Users/alicamarieenderich/202312_invariant_mass/invariant_mass_code/configs/config_pp_13.6TeV_pi0_LHC22f.yml"

# period = "LHC22o_small"
# filename = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_129486_LHC22o_pass4_small.root"
# type = "data"
# config_file = "/Users/alicamarieenderich/202312_invariant_mass/invariant_mass_code/configs/config_pp_13.6TeV_pi0_LHC22o_small.yml"

# period = "LHC22o_minBias"
# filename = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_135860_LHC22o_pass4_minBias_medium.root"
# type = "data"
# config_file = "/Users/alicamarieenderich/202312_invariant_mass_code/configs/config_pp_13.6TeV_pi0_LHC22o_medium.yml"
# cutname = "qc"
# suffix = "AnyTrack";


# period = "LHC23d1k"
# filename = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_124837_LHC23d1k.root"
# type = "mc"
# config_file ="/Users/alicamarieenderich/202312_invariant_mass_code/configs/config_pp_13.6TeV_pi0_LHC23d1k.yml"
# cutname = "qc"
# suffix = "AnyTrack";

# config_file = "/Users/alicamarieenderich/202312_invariant_mass/invariant_mass_code/configs/config_pp_13.6TeV_pi0.yml"
# with open(config_file, "r", encoding="utf-8") as config_yml:
#     config = yaml.safe_load(config_yml)
# date = datetime.date.today().strftime("%Y%m%d"); #"this_thesis" #
# folder = "/Users/alicamarieenderich/{0}_{1}_invariant_mass_plots_new_config/".format(date, period);  
# os.makedirs(folder, exist_ok=True);

# run(filename,config,type,suffix, folder);

########################################
#        Run over all Datasets         #
########################################

myDir = "/Users/juliaschlagel/Analysis240411/Analysis/";

os.makedirs(myDir, exist_ok=True);


period_array = [ "LHC22o", "LHC24b1"]

filename_array = [os.path.join(myDir,"HLtrains/data/AnalysisResults_LHC22o_full_statistics.root"),os.path.join(myDir,"HLtrains/MC/AnalysisResults_LHC24b1_full_statistics.root")]

filename_data = os.path.join(myDir,"HLtrains/data/AnalysisResults_LHC22o_full_statistics.root")
filename_mc = os.path.join(myDir,"HLtrains/MC/AnalysisResults_LHC24b1_full_statistics.root")

filename_comp_dalitz = "/Users/juliaschlagel/Analysis240411/Comparison/data_PCMDalitzResultsFullCorrection_PP.root"
filename_comp_pcm = "/Users/juliaschlagel/Analysis240411/Comparison/data_PCMResultsFullCorrection_PP.root"


type_array = ["data", "mc"]
decay_array = ["dalitz", "pcm"]
run_type_array = ["one", "both_types", "both_decays", "both_mesons"]
meson_array = ["eta", "pi0"]

start_pt_array = ["false"] #, "true"]


for mes in range(len(meson_array)):
    meson = meson_array[mes]
    if meson == "pi0":
        config_array_pcm = [os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC22o_full_pcm.yml"), os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC24b1_full_pcm.yml")]
        config_array_dalitz = [os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC22o_full_dalitz.yml"), os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC24b1_full_dalitz.yml")]
    if meson == "eta":  
        config_array_pcm = [os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_eta_LHC22o_full_pcm.yml"), os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_eta_LHC24b1_full_pcm.yml")]
        config_array_dalitz = [os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_eta_LHC22o_full_dalitz.yml"), os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_eta_LHC24b1_full_dalitz.yml")]

    for start in range(len(start_pt_array)):
        start_pt_ = start_pt_array[start]

        for dec in range(len(decay_array)):
            decay = decay_array[dec]
            if decay == "dalitz":
                config_array = config_array_dalitz
            else:
                config_array = config_array_pcm

            for r in range(len(run_type_array)):
                run_type = run_type_array[r]
                for i in range(len(period_array)):
                    period = period_array[i]
                    filename = filename_array[i]
                    cutname = "qc"
                    suffix = "AnyTrack";
                    typ = type_array[i]
                    config_file = config_array[i]
                    
                    if config_file:
                        print("config file: ", config_file)
                    with open(config_file, "r", encoding="utf-8") as config_yml:
                        config = yaml.safe_load(config_yml)
                    if start_pt_ == "false":
                        start_pt = 0
                    else:
                        start_pt = config["common"]["start_pT_bin"]
                    # start_pt = config["common"]["start_pT_bin"]
                    date = "this_analysis" #datetime.date.today().strftime("%Y%m%d");
                    #if meson == "eta":
                    folder = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/analysis_results_{0}_begin_{1}".format(meson, start_pt)
                    subfolders = ["inv_mass_analysis", "inv_yield_calculation", "Comparison"]
                    subsubfolders1 = ["Dalitz", "PCM"]
                    subsubfolders2 = ["Dalitz", "PCM"]
                    os.makedirs(folder, exist_ok=True)
                    for subfolder in subfolders:
                        os.makedirs(os.path.join(folder, subfolder), exist_ok=True)
            
                    for subsubfolder in subsubfolders1:
                        subsubfolder_path = os.path.join(folder, "inv_mass_analysis", subsubfolder)
                        os.makedirs(subsubfolder_path, exist_ok=True)
                    for subsub in subsubfolders2:
                        subsubfolder_path_2 = os.path.join(folder, "inv_yield_calculation", subsub)
                        os.makedirs(subsubfolder_path_2, exist_ok=True)

                    run(start_pt, filename,config,typ,suffix, folder)
            