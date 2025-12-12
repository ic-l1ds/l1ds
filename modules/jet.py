from analysis_tools.utils import import_root

ROOT = import_root()

class AntiKtFastJetProducer():
    def __init__(self, *args, **kwargs):
        ROOT.gInterpreter.Declare("""
            #include "fastjet/ClusterSequence.hh"
            #include <iostream>
            #include <cstring>
            #include <array>

            using Vd = const ROOT::RVec<double>&;
            using Vint = const ROOT::RVec<int>&;
            
            const int max_jets = 8;
            const int max_const = 16;

            struct output_fastjet {
                std::vector<std::vector<double>> jets;                          // [max_jets]
                std::vector<std::vector<std::vector<double>>> constituents;     // [max_jets][max_const]

                output_fastjet() :
                    jets(max_jets, {0.0, 0.0, 0.0}),
                    constituents(max_jets, std::vector<std::vector<double>>(max_const, {0.0,0.0,0.0,0.0,0.0,0.0,0.0}))
                {}
            };

                                  

            output_fastjet run_antikt(
                    int nL1BarrelExtPuppi,
                    Vd L1BarrelExtPuppi_pt, Vd L1BarrelExtPuppi_eta, Vd L1BarrelExtPuppi_phi,
                    Vd L1BarrelExtPuppi_dxy, Vd L1BarrelExtPuppi_z0, Vd L1BarrelExtPuppi_puppiWeight) {
                output_fastjet out;
                std::vector<fastjet::PseudoJet> particles;
                for (size_t i = 0; i < nL1BarrelExtPuppi; i++) {
                    particles.push_back(fastjet::PseudoJet(
                        L1BarrelExtPuppi_pt[i] * cos(L1BarrelExtPuppi_phi[i]),
                        L1BarrelExtPuppi_pt[i] * sin(L1BarrelExtPuppi_phi[i]),
                        L1BarrelExtPuppi_pt[i] * sinh(L1BarrelExtPuppi_eta[i]),
                        L1BarrelExtPuppi_pt[i]
                    ));
                    particles.back().set_user_index(i);
                }
                double R = 0.4;
                fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);

                fastjet::ClusterSequence cs(particles, jet_def);
                std::vector<fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

                for (size_t i = 0; i < jets.size(); i++) {
                    if (i >= max_jets)
                        break;
                    if (jets[i].pt() < 5)
                        continue;
                    out.jets[i] = {jets[i].pt(), jets[i].eta(), jets[i].phi()};

                    for (unsigned j = 0; j < jets[i].constituents().size(); j++) {
                        if (j >= max_const)
                            break;
                        out.constituents[i][j] =
                            {
                                jets[i].constituents()[j].pt(),
                                jets[i].constituents()[j].eta(),
                                jets[i].constituents()[j].phi(),
                                0,
                                L1BarrelExtPuppi_dxy[jets[i].constituents()[j].user_index()],
                                L1BarrelExtPuppi_z0[jets[i].constituents()[j].user_index()],
                                L1BarrelExtPuppi_puppiWeight[jets[i].constituents()[j].user_index()],
                            };
                    }
                }
                return out;
            }
        """)

    def run(self, df):
        df = df.Define("tmp",
            "run_antikt(nL1BarrelExtPuppi, L1BarrelExtPuppi_pt, L1BarrelExtPuppi_eta, L1BarrelExtPuppi_phi, L1BarrelExtPuppi_dxy, L1BarrelExtPuppi_z0, L1BarrelExtPuppi_puppiWeight)"
        ).Define("jets", "tmp.jets").Define("constituents", "tmp.constituents")

        return df, ["jets", "constituents"]


def AntiKtFastJet(*args, **kwargs):
    return lambda: AntiKtFastJetProducer()


class GenJetMatchingProducer():
    def __init__(self, *args, **kwargs):
        ROOT.gInterpreter.Declare("""
            using Vd = const ROOT::RVec<double>&;
            using Vfloat = const ROOT::RVec<float>&;
            using Vint = const ROOT::RVec<int>&;
            #include "DataFormats/Math/interface/deltaR.h"
            
            std::vector<int> jet_matches_bjet(
                std::vector<std::vector<double>> jets,
                int nGenPart,
                const Vfloat& GenPart_eta,
                const Vfloat& GenPart_phi,
                const Vint& GenPart_pdgId,
                const Vint& GenPart_status
            ) {
                std::vector<int> jet_isb(8, 0);
                for (size_t i = 0; i < 8; i++) {
                    if (jets[i][0] == 0)
                        break;
                    for (size_t igen = 0; igen < nGenPart; igen++) {
                        if (abs(GenPart_pdgId[igen]) != 5 || GenPart_status[igen] != 23)
                            continue;
                        if (reco::deltaR(GenPart_eta[igen], GenPart_phi[igen], jets[i][1], jets[i][2]) < 0.4) {
                            jet_isb[i] = 1;
                            break;
                        }
                    }
                }
                return jet_isb;
            }
        """)

    def run(self, df):
        df = df.Define(
            "jet_isb", "jet_matches_bjet(jets, nGenPart, GenPart_eta, GenPart_phi, GenPart_pdgId, GenPart_status)"
        )
        return df, ["jet_isb"]


def GenJetMatching(*args, **kwargs):
    return lambda: GenJetMatchingProducer()