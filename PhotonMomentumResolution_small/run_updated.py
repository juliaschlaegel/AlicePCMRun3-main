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
def run(filename, config, typ, suffix, folder):
    
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
    list_integral_max = config["common"]["integral_max"];
    if len(list_integral_min) != len(list_integral_max):
        return;

    list_yield_min = config["common"]["yield_min"];
    list_yield_max = config["common"]["yield_max"];
    if len(list_yield_min) != len(list_yield_max):
        return;

    nsys = len(config[typ]['subsystems']);
    print(nsys); 

    meson = config["common"]["meson"];

    if config["common"]["do_ptspectrum"] == True:
        if run_type == "one":
            if decay =="dalitz":
                outname = os.path.join(folder,"inv_mass_analysis", "Dalitz","{0}_{1}_{2}_{3}_{4}_ptspectrum_{5}_{6}TeV_{7}_{8}.root".format(date, decay, period, meson, typ, config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix))
            else:
                outname = os.path.join(folder,"inv_mass_analysis", "PCM","{0}_{1}_{2}_{3}_{4}_ptspectrum_{5}_{6}TeV_{7}_{8}.root".format(date, decay, period, meson, typ, config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix));
            print("output file name = ",outname);
            outfile = TFile(outname,"RECREATE");

            if typ == "mc":
                ana_pi0 = PairAnalyzer(meson, filename, "pi0eta-to-gammagamma-mc", "mc", filename_mc, config, decay);
            else:
                ana_pi0 = PairAnalyzer(meson, filename, "pi0eta-to-gammagamma", "data", filename_mc, config, decay);

            ana_pi0.set_arr_pt(arr_pt);
            for isys in range(0,nsys):
                ssname = config[typ]['subsystems'][isys]['name']; #subsystem name
                ana_pi0.set_subsystem(ssname);
                ana_pi0.set_xtitle("#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})");
                ana_pi0.set_ytitle("#it{p}_{T,#gamma#gamma} (GeV/#it{c}^{2})");
                print("analyze subsystem", ssname);
                
                cuts = config[typ]["subsystems"][isys]['cuts']
                cutnames = [cut['name'] for cut in cuts]
                print("cutnames", cutnames); 
                nc = len(cutnames);
                outlist_ss = THashList();
                outlist_ss.SetName(ssname);
                outlist_ss.SetOwner(True);
                for ic in range(0,nc):
                    cutname = cutnames[ic];
                    ana_pi0.set_cutname(cutname);
                    fit_parameters = cuts[ic].get('fit_parameters', [])
                    fit_limit_min = cuts[ic].get('fit_limit_min', [])
                    fit_limit_max = cuts[ic].get('fit_limit_max', [])
                    ana_pi0.set_fit_params(fit_parameters, fit_limit_min, fit_limit_max)
                    outlist_cut = THashList();
                    outlist_cut.SetName(cutname);
                    outlist_cut.SetOwner(True);
                    outlist_ss.Add(outlist_cut);
                    for ifunc in list_fit_func:
                        ana_pi0.set_fit_function(typ, ifunc);
                        outlist_func = THashList();
                        outlist_func.SetName(ifunc);
                        outlist_func.SetOwner(True);
                        outlist_cut.Add(outlist_func);
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
                            outlist_fit_range = ana_pi0.analyze_ptspectrum();
                            outlist_fit_range.SetName("fit_{0:3.2f}_{1:3.2f}_GeVc2".format(fit_min, fit_max));
                            outlist_func.Add(outlist_fit_range);
                outfile.WriteTObject(outlist_ss);
                outlist_ss.Clear();
            del ana_pi0;

            plot_pi0 = PlotInvMass(meson, outname, "pi0eta-to-gammagamma");
            plot_pi0.set_arr_pt(arr_pt);
            for isys in range(0,nsys):
                ssname = config[typ]['subsystems'][isys]['name']; #subsystem name
                plot_pi0.set_subsystem(ssname);
                print("plot subsystem", ssname);
                cuts = config[typ]["subsystems"][isys]['cuts']
                cutnames = [cut['name'] for cut in cuts]
                print("cutnames", cutnames); 
                nc = len(cutnames);
                nfit = len(list_fit_func)
                list_parameters_comparison = [];
                yield_list = [];
                for ic in range(0,nc):
                    cutname = cutnames[ic];
                    print("cutname ", cutname)
                    plot_pi0.set_cutname(cutname);               
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
                            plottingRange = [0., 0.3];
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
                            if typ == "mc":
                                plot_histo = PlotHistoInvMass(meson, outname, "pi0eta-to-gammagamma-mc");
                            else:
                                plot_histo = PlotHistoInvMass(meson, outname, "pi0eta-to-gammagamma");  
                            if decay=="dalitz":
                                outname_histo = os.path.join(folder, "inv_mass_analysis","Dalitz","{0}_{1}_{2}_InvMass_Fitparameters_{3}_{4}.pdf".format(date, decay, period, cutname, fitname))
                            else:
                                outname_histo = os.path.join(folder, "inv_mass_analysis","PCM","{0}_{1}_{2}_InvMass_Fitparameters_{3}_{4}.pdf".format(date, decay, period, cutname, fitname))
                            plot_histo.PlotHistoParameters(parameter_list, TString(outname_histo), period, 2,3,"#gamma#gamma", "{0}, cut: {1}".format(typ, cutname))

            # pdf output of all fitted histograms for each combination of cut and fit
                            for ipt in range(len(arr_pt)):
                                histo_ipt =list_plot.FindObject("h1mgg_pt{0}".format(ipt));
                                histogram_list.append(histo_ipt)
                                function_ipt = list_plot.FindObject("f1total_pt{0}".format(ipt));
                                function_list.append(function_ipt)
                            if decay =="dalitz":
                                output_name_with_yield = os.path.join(folder,"inv_mass_analysis", "Dalitz","{0}_{1}_{2}_InvMass_Overview_{3}_{4}_with_yield.pdf".format(date, decay, period, cutname, fitname))
                            else:
                                output_name_with_yield = os.path.join(folder,"inv_mass_analysis","PCM", "{0}_{1}_{2}_InvMass_Overview_{3}_{4}_with_yield.pdf".format(date, decay, period, cutname, fitname))
                            plot_pi0.PlotInvMassInPtBins(histogram_list, function_list, TString(output_name_with_yield), "", "", 
                                                        plottingRange, period, 4,3 , 0, len(arr_pt),len(arr_pt), "#gamma#gamma", 
                                                        "{0}, cut: {1}".format(typ, cutname), "yield")
                                
            # pdf output of all same and mixed scaled histograms for each cut
                            if typ == "data":
                                plot_same_mixed = PlotSameMixed(meson, outname,"pi0eta-to-gammagamma" )
                                plot_same_mixed.set_fit_range(fit_min, fit_max)
                                plot_same_mixed.set_integral_range(integral_min, integral_max)
                                for ipt in range(len(arr_pt)):
                                    histo_ipt_same =list_plot.FindObject("h1mgg_same_pt{0}".format(ipt));
                                    same_list.append(histo_ipt_same)
                                    histo_ipt_mixed = list_plot.FindObject("h1mgg_mix_scaled_pt{0}".format(ipt));
                                    mixed_list.append(histo_ipt_mixed)
                                if decay =="dalitz":
                                    output_name = os.path.join(folder, "inv_mass_analysis","Dalitz", "{0}_{1}_{2}_InvMass_Scaled_Mixed_{3}_{4}.pdf".format(date, decay, period, cutname, fitname))
                                else:
                                    output_name = os.path.join(folder, "inv_mass_analysis","PCM", "{0}_{1}_{2}_InvMass_Scaled_Mixed_{3}_{4}.pdf".format(date, decay, period, cutname, fitname))
                                plot_same_mixed.PlotSameMixedInPtBins(same_list, mixed_list, TString(output_name), "", "", plottingRange, period, 4,3 , 0, 
                                                                    len(arr_pt),len(arr_pt), "#gamma#gamma", "{0}, cut: {1}".format(typ, cutname))

            # pdf output of mass, amplitude and width for each fit and comparison of all cuts
                            parameter_comparison = [];
                            for i_parameter in range(5):
                                parameter_comparison.append(parameter_list[i_parameter]);
                            list_parameters_comparison.append(parameter_comparison);

                plot_raw_yield = PlotRawYieldInvMass(meson, outname, "pi0eta-to-gammagamma") 
                if decay == "dalitz":    
                    outname_raw_yield = os.path.join(folder, "inv_mass_analysis", "Dalitz","{0}_{1}_{2}_InvMass_RawYield_{3}.pdf".format(date, decay, period, fitname))
                else:
                    outname_raw_yield = os.path.join(folder, "inv_mass_analysis", "PCM","{0}_{1}_{2}_InvMass_RawYield_{3}.pdf".format(date, decay, period, fitname));  
                plot_raw_yield.PlotHistoYield(yield_list, TString(outname_raw_yield), "", "", period, 1,1, "#gamma#gamma", "data", cutnames)                   
                plot_parameters = PlotHistoParametersCombined(meson, outname,"pi0eta-to-gammagamma")
                if decay =="dalitz":
                    outname_histo_linear = os.path.join(folder, "inv_mass_analysis", "Dalitz", "{0}_{1}_{2}_InvMass_Parameters_Combined_{3}.pdf".format(date, decay, period, fitname));
                else:
                    outname_histo_linear = os.path.join(folder, "inv_mass_analysis", "PCM","{0}_{1}_{2}_InvMass_Parameters_Combined_{3}.pdf".format(date, decay, period, fitname));
                plot_parameters.PlotHistoParametersCombined(list_parameters_comparison, TString(outname_histo_linear), "", "", period, 2, 3, 
                                                            "#gamma#gamma", "{0}".format(typ), cutnames)
                    
                for icut in range (len(yield_list)):
                    if typ == "mc":
                        yield_used = yield_list[icut]
                        #raw_mc = plot_raw_yield.get_raw_yield(yield_used)
                        raw_mc = yield_used
                        print("#### Raw yield mc:")
                        for i in range(1, raw_mc.GetNbinsX() + 1):
                            bin_content_m = raw_mc.GetBinContent(i)
                            bin_center_m = raw_mc.GetBinCenter(i)
                            print(f"Bin {i}: Center = {bin_center_m}, Content = {bin_content_m}")
                    # for icut in range (len(yield_list)):
                        filename_mc_pcm = f'th1d_raw_mc_pcm_{icut}.pkl'
                        filename_mc_dalitz = f"th1d_raw_mc_dalitz{icut}.pkl"

                        if decay == "pcm":
                            with open(filename_mc_pcm, 'wb') as raw_mc_file_pcm:
                                pickle.dump(raw_mc, raw_mc_file_pcm)
                        else:
                            with open(filename_mc_dalitz, 'wb') as raw_mc_file_dalitz:
                                pickle.dump(raw_mc, raw_mc_file_dalitz)

                    elif typ == "data":
                        yield_used = yield_list[icut]
                        #raw_data = plot_raw_yield.get_raw_yield(yield_used)
                        raw_data = yield_used
                        #raw_data_list.append(raw_data)
                        print("#### Raw yield data:")
                        for i in range(1, raw_data.GetNbinsX() + 1):
                            bin_content_d = raw_data.GetBinContent(i)
                            bin_center_d = raw_data.GetBinCenter(i)
                            print(f"Bin {i}: Center = {bin_center_d}, Content = {bin_content_d}")
                    # for icut in range (len(yield_list)):
                        filename_data_pcm = f'th1d_raw_data_pcm_{icut}.pkl'
                        filename_data_dalitz = f"th1d_raw_data_dalitz{icut}.pkl"
                        if decay == "pcm":
                            with open(filename_data_pcm, 'wb') as raw_data_file_pcm:
                                pickle.dump(raw_data, raw_data_file_pcm)
                        else:
                            with open(filename_data_dalitz, 'wb') as raw_data_file_dalitz:
                                pickle.dump(raw_data, raw_data_file_dalitz)
                #print("############ Länge yield list: ", len(yield_list), "###############################")
                                
                # acc_class = HistogramAcceptance(filename_mc, config, decay)
                # if decay =="dalitz":
                #     acc_name = os.path.join(folder, "Dalitz", "Acceptance")
                # else:
                #     acc_name = os.path.join(folder, "PCM", "Accetpance")
                # acc = acc_class.calc_acceptance_vs_pt(TString(acc_name), "true")
                # print("Akzeptanz: ", acc)

                # if type =="mc":
                #     eff_class = Efficiency(decay)
                #     if decay =="dalitz":
                #         eff_name = os.path.join(folder, "Dalitz", "Efficiency")
                #     else:
                #         eff_name = os.path.join(folder, "PCM", "Efficiency")
                #     eff = eff_class.calc_efficiency_vs_pt(filename_mc, config, raw_mc, TString(eff_name), "true")
                    # print("################Effizienz: ", type(eff))
            del plot_pi0
            del plot_histo
            del parameter_list
            del histogram_list
            del function_list
            outfile.Close()
                
        if run_type == "both_types":
            print("run type: ", run_type)
            eff_list_dalitz = []
            eff_list_pcm = []
            acc_list_dalitz = []
            acc_list_pcm = []
            corr_list_dalitz = []
            corr_list_pcm = []
            filename_rootfile_dal = os.path.join(folder,"Dalitz","{0}_results.root".format(date))
            filename_rootfile_pcm = os.path.join(folder, "PCM","{0}_results.root".format(date))
            
            save_hist_dal = Save_Hist(filename_rootfile_dal)
            save_hist_pcm = Save_Hist(filename_rootfile_pcm)
            


            acc_class = HistogramAcceptance(filename_mc, config, decay)
            if decay =="dalitz":
                acc_name = os.path.join(folder, "Dalitz", "{0}_acceptance".format(date))
            else:
                acc_name = os.path.join(folder, "PCM", "{0}_acceptance".format(date))
            utils_mc = Utility(filename_mc, config, decay, "mc")
            utils_data = Utility(filename_data, config, decay, "data")
            nev_data = utils_data.get_Nev()
            nev_mc = utils_mc.get_Nev()
            print("##################################Number of events: ", nev_data, nev_mc)
            acc = acc_class.calc_acceptance_vs_pt(TString(acc_name), "true", utils_mc)
            if decay =="dalitz":
                acc_list_dalitz.append(acc)
                save_hist_dal.add_histogram("PCMDalitzEE", acc)
                save_hist_dal.save_histograms()
            else:
                acc_list_pcm.append(acc)
                save_hist_pcm.add_histogram("PCMPCM", acc)
                save_hist_pcm.save_histograms()

            print("Akzeptanz: ", acc, "Listenlänge: ")

            for isys in range(0, nsys):
                ssname = config[typ]["subsystems"][isys]["name"]
                print("Analyse Subsystem: ", ssname)

                cuts = config[typ]["subsystems"][isys]["cuts"]
                cutnames = [cut["name"] for cut in cuts]
                print("cutnames: ", cutnames)
                nc = len(cutnames)
                print("#########", decay, "#######", nc)
            
                for ic in range(0,nc):
                    cutname = cutnames[ic]
                    print("cutname: ", cutname)
                    #if decay == "pcm":
                    filename_mc_pcm = f'th1d_raw_mc_pcm_{ic}.pkl'
                    filename_data_pcm = f'th1d_raw_data_pcm_{ic}.pkl'
                    if os.path.exists(filename_data_pcm):
                        with open(filename_data_pcm, 'rb') as raw_data_file_pcm:
                            raw_data_pcm = pickle.load(raw_data_file_pcm)
                            save_hist_pcm.add_histogram("PCMPCM", raw_data_pcm, "h1_yield_data")
                            save_hist_pcm.save_histograms()
                        # print("####raw data file pcm geladen: ", raw_data_pcm, cutname[ic])
                        # for i in range(1, raw_data_pcm.GetNbinsX() + 1):
                        #     bin_content_d = raw_data_pcm.GetBinContent(i)
                        #     bin_center_d = raw_data_pcm.GetBinCenter(i)
                        #     print(f"Bin {i}: Center = {bin_center_d}, Content = {bin_content_d}")
                    if os.path.exists(filename_mc_pcm):
                        with open(filename_mc_pcm, 'rb') as raw_mc_file_pcm:
                            raw_mc_pcm = pickle.load(raw_mc_file_pcm)
                            save_hist_pcm.add_histogram("PCMPCM", raw_mc_pcm, "h1_yield_mc")
                            save_hist_pcm.save_histograms()
                        # print("####raw mc file pcm geladen: ", raw_mc_pcm)
                        # print("typ von raw mc: ", type(raw_mc_pcm))
                        # print("####raw mc file pcm geladen: ", raw_data_pcm, cutname[ic])
                        # for i in range(1, raw_mc_pcm.GetNbinsX() + 1):
                        #     bin_content_d = raw_mc_pcm.GetBinContent(i)
                        #     bin_center_d = raw_mc_pcm.GetBinCenter(i)
                        #     print(f"Bin {i}: Center = {bin_center_d}, Content = {bin_content_d}")
                    
                    filename_mc_dalitz = f"th1d_raw_mc_dalitz{ic}.pkl"
                    filename_data_dalitz = f"th1d_raw_data_dalitz{ic}.pkl"
                    if os.path.exists(filename_data_dalitz):
                        with open(filename_data_dalitz, 'rb') as raw_data_file_dalitz:
                            raw_data_dalitz = pickle.load(raw_data_file_dalitz)
                            save_hist_dal.add_histogram("PCMDalitzEE", raw_data_dalitz,"h1_yield_data")
                            save_hist_dal.save_histograms()
                        # print("####raw data file dalitz geladen: ", raw_data_dalitz, cutname[ic])
                        # for i in range(1, raw_data_dalitz.GetNbinsX() + 1):
                        #     bin_content_d = raw_data_dalitz.GetBinContent(i)
                        #     bin_center_d = raw_data_dalitz.GetBinCenter(i)
                        #     print(f"Bin {i}: Center = {bin_center_d}, Content = {bin_content_d}")
                    if os.path.exists(filename_mc_dalitz):
                        with open(filename_mc_dalitz, 'rb') as raw_mc_file_dalitz:
                            raw_mc_dalitz = pickle.load(raw_mc_file_dalitz)
                            save_hist_dal.add_histogram("PCMDalitzEE", raw_mc_dalitz, "h1_yield_mc")
                            save_hist_dal.save_histograms()

                        # print("####raw mc file dalitz geladen: ", raw_mc_dalitz,cutname[ic])
                        # print("typ von raw mc: ", type(raw_mc_dalitz))
                        # for i in range(1, raw_mc_dalitz.GetNbinsX() + 1):
                        #     bin_content_d = raw_mc_dalitz.GetBinContent(i)
                        #     bin_center_d = raw_mc_dalitz.GetBinCenter(i)
                        #     print(f"Bin {i}: Center = {bin_center_d}, Content = {bin_content_d}")

                    if decay == "pcm":
                        raw_mc = raw_mc_pcm
                        raw_data = raw_data_pcm
                    else:
                        raw_mc = raw_mc_dalitz
                        raw_data = raw_data_dalitz
                
                        
                    eff_class = Efficiency(decay)
                    if decay == "dalitz":
                        eff_acc_name = os.path.join(folder, "Dalitz", "Efficiency_Acceptance")
                    else:
                        eff_acc_name = os.path.join(folder, "PCM", "Efficiency_Acceptance")
                    #eff_class.eff_and_acc_vs_pt(raw_mc, utils, TString(eff_acc_name), "true")
                    if decay =="dalitz":
                        eff_name = os.path.join(folder, "Dalitz", "{0}_Efficiency_{1}".format(date, cutname)) 
                               
                    else:
                        eff_name = os.path.join(folder, "PCM", "{0}_Efficiency_{1}".format(date, cutname))
                    #utils = Utility(filename_mc, config, decay)
                    eff = eff_class.calc_efficiency_vs_pt(TString(eff_name), "true", raw_mc, utils_mc)
                    if decay =="dalitz":
                        eff_list_dalitz.append(copy.deepcopy(eff))
                        save_hist_dal.add_histogram("PCMDalitzEE", eff)
                        save_hist_dal.save_histograms()
                    else:
                        eff_list_pcm.append(copy.deepcopy(eff))
                        save_hist_pcm.add_histogram("PCMPCM", eff)
                        save_hist_pcm.save_histograms()
                    #eff = eff_class.eff_vs_pt(filename_mc, config, raw_mc, TString(eff_name), "true")
                    print("Efficiency geladen: ", type(eff), cutname)
                    #save_hist.save_histograms()
                        

                    corr_class = CorrectedYield(filename_mc, decay)
                    if decay == "dalitz":
                        pythia_name = os.path.join(folder, "Dalitz", "{0}_pythia".format(date))
                        corr_name= os.path.join(folder, "Dalitz", "{0}_corrected_yield_{1}".format(date, cutname))
                        comp_pythia_corr_name = os.path.join(folder, "Dalitz", "{0}_comp_pythia_corr_yield_{1}".format(date, cutname))
                    else:                    
                        pythia_name = os.path.join(folder, "PCM", "{0}_pythia".format(date))
                        corr_name = os.path.join(folder, "PCM", "{0}_corrected_yield_{1}".format(date, cutname))
                        comp_pythia_corr_name = os.path.join(folder, "PCM", "{0}_comp_pythia_corr_yield_{1}".format(date, cutname))
                    
                    pyth = corr_class.plot_pythia(utils_mc, TString(pythia_name), "true")
                    if decay == "pcm":
                        save_hist_pcm.add_histogram("PCMPCM", pyth)
                        save_hist_pcm.save_histograms()
                        #save_hist.save_histograms()
                    if decay == "dalitz":
                        print("dalitz")
                        if ic == (nc-1):
                            save_hist_dal.add_histogram("PCMDalitzEE", pyth)
                            save_hist_dal.save_histograms()
                        # save_hist_dal.save_histograms()
                    #corr = corr_class.calculate_corrected_yield(TString(corr_name), filename_mc, config, raw_data, raw_mc,  "true")
                    # corr = corr_class.calc_corr_vs_pt(TString(corr_name), "true", raw_data, decay, raw_mc, filename_mc, config, eff, acc)
                    corr = corr_class.calc_corr_yield_vs_pt(TString(corr_name), "true", raw_data, eff, acc, utils_mc)
                    if decay == "pcm":
                        corr_list_pcm.append(copy.deepcopy(corr))
                        save_hist_pcm.add_histogram("PCMPCM", corr)
                        save_hist_pcm.save_histograms()

                    else:
                        corr_list_dalitz.append(copy.deepcopy(corr))
                        if ic == (nc-1):
                            save_hist_dal.add_histogram("PCMDalitzEE", corr)
                            save_hist_dal.save_histograms()
                        print("###################Typ des in der Liste gespeicherten Histogramms: ", type(corr))
                        print("Liste corr_yield_dalitz: ", corr_list_dalitz)
                    print("Corrected_yield importiert: ", type(corr))
                    corr_class.compare_pythia_corr_yield(raw_data, eff, acc, utils_mc, TString(comp_pythia_corr_name), "true")

            
            
                    if ic == (len(cutnames)-1):
                        if decay == "dalitz":
                            print("#######Länge Liste corr_dalitz", len(corr_list_dalitz),corr_list_dalitz)
                            for ncut in range(len(cutnames)):
                                corr_yield_dal = corr_list_dalitz[ncut]
                                filename_corr_y_dal = f"corr_dalitz_{ncut}"
                                with open(filename_corr_y_dal,"wb") as corr_yield_dalitz_file:
                                    pickle.dump(corr_yield_dal, corr_yield_dalitz_file)
                        if decay =="pcm":
                            print("#######Länge Liste corr_pcm", len(corr_list_pcm),corr_list_pcm)
                            for ncut in range(len(cutnames)):
                                corr_yield_pcm_ = corr_list_pcm[ncut]
                                filename_corr_y_pcm = f"corr_pcm_{ncut}"
                                with open(filename_corr_y_pcm, "wb") as corr_yield_pcm_file:
                                    pickle.dump(corr_yield_pcm_, corr_yield_pcm_file)
                        
                        comp_class = Compare(filename_mc)
                        if decay == "dalitz":
                            comp_name_cuts_eff = os.path.join(folder, "Dalitz", "comp_cuts_eff")
                            comp_name_cuts_corr = os.path.join(folder, "Dalitz", "comp_cuts_corr")
                            comp_eff_cuts = comp_class.compare_cuts(eff_list_dalitz, cutnames, "dalitz", filename_mc, config, TString(comp_name_cuts_eff), "eff", utils_mc)
                            # save_hist.add_histogram("PCMDalitzEE", comp_eff_cuts)
                            # save_hist.save_histograms()
                            #comp_class.compare_cuts(corr_list_dalitz, cutnames, "dalitz", filename_mc, config, TString(comp_name_cuts_corr),"corr")
                        else:
                            comp_name_cuts_eff = os.path.join(folder, "PCM", "comp_cuts_eff")
                            comp_name_cuts_corr = os.path.join(folder, "PCM", "comp_cuts_corr")
                            comp_eff_cuts = comp_class.compare_cuts(eff_list_pcm, cutnames, "pcm", filename_mc, config, TString(comp_name_cuts_eff), "eff", utils_mc)
                            # save_hist.add_histogram("PCMPCM", comp_eff_cuts)
                            # save_hist.save_histograms()
                            #comp_class.compare_cuts(corr_list_pcm, cutnames, "pcm", filename_mc, config, TString(comp_name_cuts_corr),"corr")
            
                    eff_pc = eff * 100
                    acc_pc = acc * 100   
                    acc_name_comp = os.path.join(folder, "Comparison", "{0}_Comp_acc_dalitz_pcm".format(date))
                    eff_name_comp = os.path.join(folder, "Comparison", "{0}_Comp_eff_dalitz_pcm_{1}".format(date, cutname))
                    comp_acc_13Tev_dalitz_name = os.path.join(folder, "Comparison", "{0}_Comp_13TeV_acc_dalitz".format(date))
                    comp_eff_13Tev_dalitz_name = os.path.join(folder, "Comparison", "{0}_Comp_13TeV_eff_dalitz_{1}".format(date, cutname))
                    comp_corr_13Tev_dalitz_name = os.path.join(folder, "Comparison", "{0}_Comp_13TeV_corr_yield_dalitz".format(date))
                    comp_acc_13Tev_pcm_name = os.path.join(folder, "Comparison", "{0}_Comp_13TeV_acc_pcm".format(date))
                    comp_eff_13Tev_pcm_name = os.path.join(folder, "Comparison", "{0}_Comp_13TeV_eff_pcm".format(date))
                    comp_corr_13Tev_pcm_name = os.path.join(folder, "Comparison", "{0}_Comp_13TeV_corr_yield_pcm".format(date))
                    comp_class = Compare(filename_mc)
                    comp_acc_dal_pcm =comp_class.plot_acc_dalitz_vs_pcm(filename_mc, config_array_dalitz, config_array_pcm, TString(acc_name_comp), "true")
                    comp_eff_dal_pcm = comp_class.plot_eff_dalitz_vs_pcm(filename_mc, config, raw_mc_dalitz, raw_mc_pcm, TString(eff_name_comp),"true")
                    # save_hist.add_histogram("Comparison", comp_acc_dal_pcm)
                    # save_hist.add_histogram("Comparison", comp_eff_dal_pcm)
                    # save_hist.save_histograms()
                    #comp_class.compare_corr_dalitz_pcm(raw_data_dalitz, eff_list_dalitz[ic], acc_list_dalitz[0], raw_data_pcm, eff_list_pcm[0], acc_list_pcm[0])
                    #comp_class.compare_corr_dalitz_pcm(raw_data_dalitz, eff_list_dalitz[ic], acc_list_dalitz[0], raw_data_pcm, eff_list_pcm[icut], acc_list_pcm[0])
                    if decay == "dalitz":
                        comp_class.compare_13tev(filename_comp_dalitz, "acc", acc_pc, "dalitz", TString(comp_acc_13Tev_dalitz_name) , "true")
                        comp_class.compare_13tev(filename_comp_dalitz, "eff", eff, "dalitz",TString(comp_eff_13Tev_dalitz_name), "true")
                        comp_class.compare_13tev(filename_comp_dalitz, "corr", corr, "dalitz", TString(comp_corr_13Tev_dalitz_name), "true")
                    if decay == "pcm":
                        comp_class.compare_13tev(filename_comp_pcm, "acc", acc, "pcm", TString(comp_acc_13Tev_pcm_name) , "true")
                        comp_class.compare_13tev(filename_comp_pcm, "eff", eff, "pcm", TString(comp_eff_13Tev_pcm_name), "true")
                        comp_class.compare_13tev(filename_comp_pcm, "corr", corr, "pcm", TString(comp_corr_13Tev_pcm_name), "true")
                    
                    comp_eff_times_acc_name = os.path.join(folder, "Comparison", "{0}Comp_eff_times_acc".format(date))
                    comp_class.plot_eff_times_acc(filename_mc, config, raw_mc_dalitz, raw_mc_pcm, TString(comp_eff_times_acc_name), "true", utils_mc)
            #save_hist.save_histograms()
                
        if run_type == "both_decays":
            print("run type: ", run_type)
            final_dec = decay_array[(len(decay_array)-1)]
            print("final decay: ", final_dec)
            config_file_dalitz = config_array_dalitz[0]
            config_file_pcm = config_array_pcm[0]
            with open(config_file_dalitz, "r", encoding="utf-8") as config_dalitz_yml:
                config_dalitz = yaml.safe_load(config_dalitz_yml)
            with open(config_file_pcm, "r", encoding="utf-8") as config_pcm_yml:
                config_pcm = yaml.safe_load(config_pcm_yml)
            
            for isys in range(0,nsys):
                ssname_dalitz = config_dalitz[typ]['subsystems'][isys]['name']; #subsystem name
                print("analyze subsystem", ssname_dalitz);
                ssname_pcm = config_pcm[typ]['subsystems'][isys]['name']; #subsystem name
                print("analyze subsystem", ssname_pcm);
                cuts_dalitz = config_dalitz["data"]["subsystems"][isys]['cuts']
                cutnames_dalitz = [cut['name'] for cut in cuts_dalitz]
                print("cutnames dalitz", cutnames_dalitz); 
                cuts_pcm = config_pcm["data"]["subsystems"][isys]['cuts']
                cutnames_pcm = [cut['name'] for cut in cuts_pcm]
                print("cutnames pcm", cutnames_pcm)
                nc_dalitz = len(cutnames_dalitz);
                #print("nc_dalitz", nc_dalitz, range(0,nc_dalitz-1), "cutname 0", cutnames_dalitz[0], "cutname 1", cutnames_dalitz[1])
                nc_pcm = len(cutnames_pcm)
                print("nc_pcm: ", nc_pcm)
                
                for ic in range(0,nc_dalitz):
                    cutname = cutnames_dalitz[ic]
                    print("cutname dalitz: ", cutname)
                    filename_corr_y_dal = f"corr_dalitz_{ic}"
                    with open(filename_corr_y_dal, "rb") as corr_dalitz_file:
                        corr_yield_dalitz = pickle.load(corr_dalitz_file)
                    print("corrected_yield_dalitz geladen: ", corr_yield_dalitz, ic)
                for nc in range(0, nc_pcm):
                    filename_corr_y_pcm = f"corr_pcm_{nc}"
                    with open(filename_corr_y_pcm, "rb") as corr_pcm_file:
                        corr_yield_pcm = pickle.load(corr_pcm_file)
                    print("corrected_yield_pcm geladen: ", corr_yield_pcm)
                comp_class = Compare(filename_mc)
                comp_corr_name = os.path.join(folder, "Comparison", "{0}_Comp_corr_yield_dalitz_pcm".format(date))
                comp_class.compare_corr_dalitz_pcm(corr_yield_dalitz, corr_yield_pcm, TString(comp_corr_name), "true")
            
            # if decay == final_dec:
            #     for isys in range(0, nsys):
            #         ssname = config[typ]["subsystems"][isys]["name"]
            #         print("Analyse Subsystem: ", ssname)

            #         cuts = config[typ]["subsystems"][isys]["cuts"]
            #         cutnames = [cut["name"] for cut in cuts]
            #         print("cutnames: ", cutnames)
            #         nc = len(cutnames)
                
            #         for ic in range(0,nc):
            #             cutname = cutnames[ic]
            #             print("cutname: ", cutname)
            #             filename_corr_y_dal = f"corr_dalitz_{ic}"
            #             with open(filename_corr_y_dal, "rb") as corr_dalitz_file:
            #                 corr_yield_dalitz = pickle.load(corr_dalitz_file)
            #             filename_corr_y_pcm = f"corr_pcm_{ic}"
            #             with open(filename_corr_y_pcm, "rb") as corr_pcm_file:
            #                 corr_yield_pcm = pickle.load(corr_pcm_file)
                    
                    
                    # filename_mc_pcm = f'th1d_raw_mc_pcm_{ic}.pkl'
                    # filename_data_pcm = f'th1d_raw_data_pcm_{ic}.pkl'
                    # if os.path.exists(filename_data_pcm):
                    #     with open(filename_data_pcm, 'rb') as raw_data_file_pcm:
                    #         raw_data_pcm = pickle.load(raw_data_file_pcm)



            
            
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
#LHC22o_pass6_small, HLtrains/data/AnalysisResults_HL197579_LHC22o_pass6_small.root

filename_array = [os.path.join(myDir,"HLtrains/data/AnalysisResults_HL197579_LHC22o_pass6_small.root"),os.path.join(myDir,"HLtrains/MC/AnalysisResults_HL197593_LHC24b1.root")]

filename_data = os.path.join(myDir,"HLtrains/data/AnalysisResults_HL197579_LHC22o_pass6_small.root")
filename_mc = os.path.join(myDir,"HLtrains/MC/AnalysisResults_HL197593_LHC24b1.root")

filename_comp_dalitz = "/Users/juliaschlagel/Analysis240411/Comparison/data_PCMDalitzResultsFullCorrection_PP.root"
filename_comp_pcm = "/Users/juliaschlagel/Analysis240411/Comparison/data_PCMResultsFullCorrection_PP.root"
#config_1 = os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC22o_small_pcm.yml")


type_array = ["data", "mc"]
decay_array = ["dalitz", "pcm"]
run_type_array = ["one", "both_types", "both_decays"]
#decay = "pcm"



config_array_pcm = [os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution_small/configs/config_pp_13.6TeV_pi0_LHC22o_small_pcm.yml"), os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution_small/configs/config_pp_13.6TeV_pi0_LHC24b1_pcm.yml")]
config_array_dalitz = [os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution_small/configs/config_pp_13.6TeV_pi0_LHC22o_small_dalitz.yml"), os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution_small/configs/config_pp_13.6TeV_pi0_LHC24b1_dalitz.yml")]
# config_array_pcm = [os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC22o_full_pcm.yml"), os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC24b1_full_pcm.yml")]
# config_array_dalitz = [os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC22o_full_dalitz.yml"), os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC24b1_full_dalitz.yml")]

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
            with open(config_file, "r", encoding="utf-8") as config_yml:
                config = yaml.safe_load(config_yml)
            date = "this_analysis" #datetime.date.today().strftime("%Y%m%d");
            folder = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/analysis_results_small"
            subfolders = ["inv_mass_analysis", "Dalitz", "PCM", "Comparison"]
            subsubfolders = ["Dalitz", "PCM"]
            #folder = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/{0}_{1}_invariant_mass_plots/".format(date, period);  
            os.makedirs(folder, exist_ok=True)
            for subfolder in subfolders:
                os.makedirs(os.path.join(folder, subfolder), exist_ok=True)
    
            for subsubfolder in subsubfolders:
                subsubfolder_path = os.path.join(folder, "inv_mass_analysis", subsubfolder)
                os.makedirs(subsubfolder_path, exist_ok=True)

            run(filename,config,typ,suffix, folder)
    