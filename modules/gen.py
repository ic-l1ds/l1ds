from analysis_tools.utils import import_root

ROOT = import_root()

class GenTauHadProducer():
    def __init__(self, *args, **kwargs):
        ROOT.gInterpreter.Declare("""
            using Vd = const ROOT::RVec<double>&;
            using Vint = const ROOT::RVec<int>&;
            ROOT::RVec<int> get_had_tau(int nGenPart, Vint GenPart_pdgId, Vint GenPart_genPartIdxMother) {
                ROOT::RVec<int> is_had_tau(nGenPart, 0);
                for (size_t i = 0; i < nGenPart; i++) {
                    if (abs(GenPart_pdgId[i]) != 15)
                        continue;
                    is_had_tau[i] = 1;
                    for (size_t j = 0; j < nGenPart; j++) {
                        if (GenPart_genPartIdxMother[j] != i)
                            continue;
                        if (abs(GenPart_pdgId[j]) == 11 || abs(GenPart_pdgId[j]) == 13 || abs(GenPart_pdgId[j]) == 15) { // including 15 to avoid multiple copies
                            is_had_tau[i] = 0;
                            break;
                        }
                    }
                }
                return is_had_tau;
            }

        """)

    def run(self, df):
        df = df.Define("GenPart_isHadronicTau",
            "get_had_tau(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother)")

        return df, ["GenPart_isHadronicTau"]


def GenTauHad(*args, **kwargs):
    return lambda: GenTauHadProducer()
