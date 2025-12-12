from analysis_tools import ObjectCollection, Category, Process, Dataset, Feature, Systematic
from analysis_tools.utils import DotDict
from analysis_tools.utils import join_root_selection as jrs
from plotting_tools import Label
from collections import OrderedDict

from config.hh_2024 import Config as base_config


class Config(base_config):
    def add_processes(self):
        processes = [
            Process("htt", Label("H->tautau"), color=(0, 0, 0), isSignal=True),
            Process("hh4b", Label("HH->4b"), color=(0, 0, 0), isSignal=True),
        ]

        process_group_names = {
            "default": [
                
            ],
        }

        process_training_names = {}

        # adding reweighed processes
        processes = ObjectCollection(processes)

        return ObjectCollection(processes), process_group_names, process_training_names

    def add_datasets(self):
        datasets = [
            Dataset("htt",
                dataset="/GluGluHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8_PU200/"
                    "jleonhol-Phase2Nano15_wjets_reg_new-00000000000000000000000000000000/USER",
                process=self.processes.get("htt"),
                xs=0.02964,
                tags=["ul"]),
            Dataset("hh4b",
                dataset="/GluGluToHHTo4B_node_SM_TuneCP5_14TeV-amcatnlo-pythia8/"
                    "jleonhol-Phase2Nano15_wjets_reg_new-00000000000000000000000000000000/USER",
                process=self.processes.get("hh4b"),
                xs=0.02964,
                tags=["ul"]),
            Dataset("hh4b_puppi",
                dataset="/GluGluToHHTo4B_node_SM_TuneCP5_14TeV-amcatnlo-pythia8/"
                    "jleonhol-Phase2Nano15_wpfpuppi-00000000000000000000000000000000/USER",
                process=self.processes.get("hh4b"),
                xs=0.02964,
                tags=["ul"]),
        ]
        return ObjectCollection(datasets)


config = Config("hh_phase2", year="Phase 2", ecm=14, lumi_pb=3000)