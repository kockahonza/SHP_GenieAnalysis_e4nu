#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>
#include <tuple>

#include "GenieAnalysis/GenieAnalysisDeltaStudies.h"
#include "GenieAnalysis/misc.h"

class GenieAnalysis1Proton : public GenieAnalysisDeltaStudies {
  public:
    GenieAnalysis1Proton(const char *filename, const char *output_filename,
                         // Specify the analysis - which stages, properties and types to do histograms for
                         const vector<string> &stages = {}, const vector<string> &properties = {},
                         const vector<string> &types = {},
                         const vector<GenieAnalysisAutoHistograms::AutoVsPlot> &vs_property_plots = {},
                         // Select run
                         const Target &target = Target::C12, const BeamEnergy &beam_energy = BeamEnergy::MeV_2261,
                         // Pass this to GenieAnalysis
                         const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name),
          GenieAnalysisDeltaStudies(filename, output_filename, stages, properties, types, vs_property_plots, target,
                                    beam_energy) {
        m_known_properties.insert(m_new_known_properties.begin(), m_new_known_properties.end());
    }

    Double_t passesCuts() override {
        Double_t weight{GenieAnalysisDeltaStudies::passesCuts()};

        if (m_fs_number_of_protons != 1) {
            return 0;
        }

        return weight;
    }
};

class GenieAnalysisOptimalWCut : public GenieAnalysisDeltaStudiesCuts {
  public:
    const Double_t m_Wmin;
    const Double_t m_Wmax;

    int m_number_of_delta_events;
    double m_result;

  public:
    GenieAnalysisOptimalWCut(const char *filename, const Double_t Wmin, const Double_t Wmax,
                             const Target &target = Target::C12, const BeamEnergy &beam_energy = BeamEnergy::MeV_2261,
                             const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name),
          GenieAnalysisDeltaStudiesCuts(filename, target, beam_energy), m_Wmin{Wmin}, m_Wmax{Wmax} {}

    Double_t passesCuts() override {
        Double_t weight{GenieAnalysisDeltaStudiesCuts::passesCuts()};

        if ((m_reconstructed_W < m_Wmin) || (m_reconstructed_W > m_Wmax)) {
            return 0;
        }

        if ((m_fs_number_of_neutrons == 1) && (m_passed_pi_plus.size() == 1) ||
            (m_fs_number_of_protons == 1) && (m_passed_pi_minus.size() == 1)) {
            return weight;
        }

        return 0;
    }

    void runPreAnalysis() override { m_number_of_delta_events = 0; }

    void useEntry(const Double_t &weight) override { m_number_of_delta_events += 1; }

    void runPostAnalysis(const Long64_t &number_of_entries) override {
        m_result = static_cast<double>(m_number_of_delta_events) / number_of_entries;
    }
};

void runWScan() {
    using std::cout, std::endl;

    vector<Double_t> Wmins{0, 0.4, 0.8, 1.2};
    vector<Double_t> Wmaxs{1.3, 1.4, 1.5, 1.6};

    vector<tuple<Double_t, Double_t, double>> results;

    GenieAnalysisOptimalWCut *ga;

    for (auto Wmin : Wmins) {
        for (auto Wmax : Wmaxs) {
            ga = new GenieAnalysisOptimalWCut(
                "/home/honza/Sync/University/CurrentCourses/SHP/data/Genie_gst_2000000.root", Wmin, Wmax);
            ga->runAnalysis();
            results.push_back({Wmin, Wmax, ga->m_result});
            delete ga;
        }
    }

    cout << "--------------------" << endl;
    for (auto &[Wmin, Wmax, result] : results) {
        cout << Wmin << "," << Wmax << "," << result << endl;
    }
    cout << "--------------------" << endl;
}

void scratch(int argc, char *argv[]) {
    string input_file, output_file;
    if (argc == 2) {
        std::string arg{argv[1]};
        if (arg == "local") {
            input_file = "/home/honza/Sync/University/CurrentCourses/SHP/data/Genie_gst_2000000.root";
            output_file = "output_local";
        } else if (arg == "full") {
            input_file = "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_2261GeV/"
                         "apapadop_SuSav2_C12_2261GeV_master.root";
            output_file = "output_full";
        }
    } else {
        std::cout << "Needs an argument to specify which way to run (should be \"local\" or \"full\"), check the code"
                  << std::endl;
    }

    GenieAnalysis1Proton ga1{input_file.c_str(),
                             (output_file + "_1p.root").c_str(),
                             {},
                             {
                                 "W",
                                 "Ws",
                                 "resc",
                                 "reco_W",
                                 "bjorken_x",
                                 "el_p",
                                 "passed_num_pip",
                                 "passed_num_pim",
                                 "ps_num_pip",
                                 "ps_num_pim",
                                 "ps_num_protons",
                                 "ps_num_neutrons",
                                 "fs_num_pip",
                                 "fs_num_pim",
                                 "fs_num_protons",
                                 "fs_num_neutrons",
                             },
                             {},
                             {
                                 {"W", "Ws", {}},
                                 {"W", "el_p", {}},
                                 {"W", "bjorken_x", {}},
                                 {"Ws", "bjorken_x", {}},
                                 {"passed_num_pip", "passed_num_pim", {}},
                                 {"fs_num_pim", "passed_num_pim", {}},
                                 {"fs_num_pip", "passed_num_pip", {}},
                             }};
    ga1.runAnalysis();
}

int main(int argc, char *argv[]) {

    runWScan();

    /* scratch(argc, argv); */

    return 0;
}
