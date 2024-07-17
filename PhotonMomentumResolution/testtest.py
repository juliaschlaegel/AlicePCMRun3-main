import os 
from ROOT import TFile
import yaml
from utility import Utility

def run(filename_mc, config):
    rootfile_mc = TFile.Open(filename_mc,"READ")
    nsys = len(config[type]['subsystems']);
    print("nsys",nsys); 
    for isys in range(0, nsys):
        ssname = config[type]['subsystems'][isys]['mc_name']; #subsystem name
        print("plot subsystem", ssname)
        list_ss   = rootfile_mc.Get(ssname);
        print("ss list: ", list_ss)
                #for ic in range(0,nc):
                #    cutname = cutnames[ic];
                #    list_ss_cut = list_ss.FindObject(cutname);
                #    print("cutname ", cutname)             
                #for ifit in range(0, nfit): 
        fitname = "gausexp";
        list_fitname  = list_ss.FindObject("gausexp");
        print("fitname ", fitname)


if __name__ == "__main__":
    
    period_array = ["LHC24b1"]
    filename_array = [ "/Users/juliaschlagel/Analysis240411/Analysis/HLtrains/MC/AnalysisResults_LHC24b1_full_statistics.root"]
    type_array = ["mc"]
    
    config_array = ["/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC24b1_full_pcm.yml"]
    #config_array_dalitz = [os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC22o_full_dalitz.yml"), os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC24b1_full_dalitz.yml")]

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
        
        folder = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/param_comp"
        os.makedirs(folder, exist_ok=True);
        filename_data = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/analysis_results_subtr_pi0/inv_mass_analysis/PCM/this_analysis_pcm_LHC22o_pi0_data_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
        #filename_data = "/Users/alicamarieenderich/this_thesis_LHC22f_invariant_mass_plots/this_thesis_LHC22f_pi0_data_ptspectrum_pp_13.6TeV_LHC22f_AnyTrack.root"
        filename_mc = "/Users/juliaschlagel/Analysis240411/Analysis/AlicePCMRun3-main/analysis_results_subtr_pi0/inv_mass_analysis/PCM/this_analysis_pcm_LHC24b1_pi0_mc_ptspectrum_pp_13.6TeV_LHC22o_AnyTrack.root"
        #filename_mc = "/Users/alicamarieenderich/this_thesis_LHC23d1k_invariant_mass_plots/this_thesis_LHC23d1k_pi0_mc_ptspectrum_pp_13.6TeV_LHC22f_AnyTrack.root"
        run(filename_mc, config)
        print("input period_str", period_str)





#class test:
    
    # def __init__(self, filename):
    #     self.rootfile= TFile.Open(filename, "READ")

    # def aufruf_NEV(self, decay):
    #     utils_mc = Utility(filename_mc, decay, "mc")
    #     utils_data = Utility(filename_data, decay, "data")
    #     Nev_data = utils_data.get_Nev()

    #     Nev_mc = utils_mc.get_Nev()
    #     print("Number of events, mc: ", Nev_mc)
    #     print("Number of events, data: ", Nev_data)
#     def __init__(self, filename, config, decay, typ): 
#         self.decay = decay
#         self.typ = typ
#         self.filename = TFile.Open(filename, "READ")
#         # dir_name1 = "pi0eta-to-gammagamma-mc/Generated"
#         # if self.decay == "dalitz":
#         #     #if self.typ == "mc":
#         #     dir_name1 = "pi0eta-to-gammagamma-mc-pcmdalitzee/Generated"
#         # if self.decay == "pcm":
#         #     dir_name1 = "pi0eta-to-gammagamma-mc-pcmpcm/Generated"
#         # dir_name2 = "Pi0"
#         # self.dir = self.filename.Get(dir_name1).Get(dir_name2)
#         self.config = config

#     def get_hist (self, hist_name):
#         dir_name1 = "pi0eta-to-gammagamma-mc/Generated"
        
#         if self.decay == "dalitz":
#             dir_name2 = "PCMDalitzEE"
#         if self.decay == "pcm":
#             dir_name2 = "PCMPCM"
#         #dir = self.filename.Get(dir_name1).Get(dir_name2)
#         #self.hist_name = dir.Get(hist_name)
#         #print(hist_name)
#         self.hist_name = self.filename.Get(dir_name1).FindObject(dir_name2).FindObject(hist_name)
#         if self.hist_name:
#             print("histogram found")
#         if not self.hist_name:
#             raise ValueError("Histogram not found")
#         return self.hist_name
    

#     def get_Nev (self): #dir_name2):
#         if self.typ == "mc":
#             direv = self.filename.Get("pi0eta-to-gammagamma-mc/Event")
#         if self.typ == "data":
#             direv = self.filename.Get("pi0eta-to-gammagamma/Event")
#         if self.decay == "pcm":
#             dir_name2 = "PCMPCM"
#         else:
#             dir_name2 = "PCMDalitzEE"
#         events = direv.FindObject(dir_name2).FindObject("after").FindObject("hCollisionCounter")
#         self.Nev = events.GetBinContent(9)
#         print("Nev :", self.Nev)
#         return self.Nev


#     # def __init__(self, filename, dirname):
#     #     self.rootfile = TFile.Open(filename, "READ");
#     #     self.rootdir = self.rootfile.Get(dirname);
#     #     #self.list_ev = self.rootdir.Get("Event")
#     #     self.list_pair = self.rootdir.Get("Pair")
#     # def get_events(self):
#     #     # if self.list_ev:
#     #     #     print("list_ev: ", self.list_ev)
#     #     # self.list_ev_af = self.list_ev.Get("after")
#     #     # if self.list_ev_af:
#     #     #     print("after gefunden")
#     #     if self.list_pair:
#     #         print("list_pair: ", self.list_pair)
#     #     #if self.list_pair:
#     #     h2mc_help = self.list_pair.Get("Pi0")
#     #     if h2mc_help:
#     #         print("Pi0 gefunden")
#     #         h2mc = h2mc_help.Get("hMggPt_Primary") #.Clone("h2mc");
#     #         if h2mc:
#     #             print("hMggPt gefunden")

        
        

#             #h2same_help = self.list_pair.Get("same")
#             #h2same = self.list_pair.FindObject("hMggPt_Same").Clone("h2same")
#             # if h2same_help:
#             #     print("same geladen")
#             #     h2same = h2same_help.Get("hMggPt").Clone("h2same")
#             #     if h2same:
#             #         print("h2same loaded", type(h2same))



#     # def get_Nev (self): #dir_name2):
#     #     #if decay == "dalitz":
#     #     direv = self.rootfile.Get("pi0eta-to-gammagamma-mc-pcmdalitzee/Event")
#     #     #if decay == "pcm":
#     #     #direv = self.filename.Get("pi0eta-to-gammagamma-mc-pcmpcm/Event")
#     #     # direv = self.filename.Get("pi0eta-to-gammagamma-mc/Event")
#     #     # if self.decay == "pcm":
#     #     #     dir_name2 = "PCMPCM"
#     #     # else:
#     #     #     dir_name2 = "PCMDalitzEE"
#     #     # events = direv.FindObject(dir_name2).FindObject("after").FindObject("hCollisionCounter")
#     #     events = direv.Get("after").Get("hCollisionCounter")
#     #     self.Nev = events.GetBinContent(10)
#     #     print("Nev :", self.Nev)
#     #     return self.Nev
        
        
        
#         # self.list_ev.ls()
#         # list_ev_af = self.list_ev.Get("after")
#         # if list_ev_af:
#         #     print("after gefunden")
#         # else:
#         #     print("nicht gefunden")
#         # h1ev   = list_ev_af.FindObject("hCollisionCounter");
#         # print(type(h1ev))
#         # nentries = h1ev.GetEntries()
#         # print(nentries)




# if __name__ == "__main__":

#     myDir = "/Users/juliaschlagel/Analysis240411/Analysis/";

#     os.makedirs(myDir, exist_ok=True);


#     # filename_mc = os.path.join(myDir,"HLtrains/MC/AnalysisResults_LHC24b1_full_statistics.root")
#     # filename_data = os.path.join(myDir,"HLtrains/data/AnalysisResults_LHC22o_full_statistics.root")

#     filename_data = os.path.join(myDir, "HLtrains/data/AnalysisResults_HL197579_LHC22o_pass6_small.root")
#     filename_mc = os.path.join(myDir,"HLtrains/MC/AnalysisResults_HL197593_LHC24b1.root")

#     #config_array_pcm = [os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC22o_full_pcm.yml"), os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC24b1_full_pcm.yml")]
#     #config_array_dalitz = [os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC22o_full_dalitz.yml"), os.path.join(myDir,"AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC24b1_full_dalitz.yml")]
#     config_data_dalitz = os.path.join(myDir, "AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC22o_small_dalitz.yml")
#     config_mc_dalitz = os.path.join(myDir, "AlicePCMRun3-main/PhotonMomentumResolution/configs/config_pp_13.6TeV_pi0_LHC24b1_dalitz.yml")
#     with open(config_mc_dalitz, "r", encoding="utf-8") as config_yml:
#         config_mc = yaml.safe_load(config_yml)
#     with open(config_data_dalitz, "r", encoding="utf-8") as config_yml:
#         config_data = yaml.safe_load(config_yml)
#     #test_cl = test(filename_mc)
#     #test_cl.aufruf_NEV("pcm")
#     #test_cl = test(filename_mc, "pi0eta-to-gammagamma-mc-pcmdalitzee")
#     #test_cl = test(filename_data, "pi0eta-to-gammagamma-pcmdalitzee")
#     #test_cl.get_events()
#     #test_cl.get_Nev()
#     test_cl = test(filename_mc, config_mc, "dalitz", "mc")
#     test_cl.get_Nev()

