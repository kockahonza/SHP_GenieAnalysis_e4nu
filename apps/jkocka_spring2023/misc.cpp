#include <iostream>

#include "TLorentzVector.h"

#include "GenieAnalysis/GAAutoHistograms.h"

class Faff : public GenieAnalysis {
  public:
    std::map<int, int> all_resc_counts;
    std::map<int, int> pip_resc_counts;

    using GenieAnalysis::GenieAnalysis;

    Double_t passesCuts() override { return (m_ge.res && (m_ge.resid == 0)); }

    void useEntry(const Double_t &weight = 1) override {
        for (int i = 0; i < m_ge.ni; i++) {
            all_resc_counts[m_ge.resc[i]] += 1;

            if (m_ge.pdgi[i] == 211) {
                pip_resc_counts[m_ge.resc[i]] += 1;
            }
        }
    }
};

class Faff2 : public GAAutoHistograms {
  public:
    const double m_beam_energy_val{2.261};
    const TLorentzVector m_beam_V4{0, 0, sqrt(m_beam_energy_val *m_beam_energy_val - mass_electron * mass_electron),
                                   m_beam_energy_val};

    double m_reco_W;
    TLorentzVector m_el_V4;
    double m_reco_Q2;
    double m_reco_x;

    map<string, AutoProperty> m_new_known_properties{
        // Reconstructed properties of the event and then pi
        {"reco_W", {"Reconstructed W [GeV]", {500, 0, 4}, [this]() { return m_reco_W; }}},
        {"reco_Q2", {"Reconsotructed Q^2 [GeV]", {500, 0, 4}, [this]() { return m_reco_Q2; }}},
        {"reco_x", {"Reconstructed Bjorken x", {500, 0, 1.01}, [this]() { return m_reco_x; }}},
    };

    Faff2(const char *filename, const char *output_filename, const vector<string> &stages = {},
          const vector<string> &properties = {}, const vector<string> &types = {},
          const vector<AutoVsPlot> vs_property_plots = {}, const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name), GAAutoHistograms{filename,   output_filename, stages,
                                                                    properties, types,           vs_property_plots} {
        m_known_properties.insert(m_new_known_properties.begin(), m_new_known_properties.end());
    }

    Double_t passesCuts() override {
        m_el_V4.SetPxPyPzE(m_ge.pxl, m_ge.pyl, m_ge.pzl, m_ge.El);
        m_el_V4.SetPhi(m_el_V4.Phi() + TMath::Pi());

        const TLorentzVector el_change{m_el_V4 - m_beam_V4};
        m_reco_Q2 = -el_change.Mag2();
        const double nu = -el_change.E();
        m_reco_x = m_reco_Q2 / (2 * mass_proton * nu);
        m_reco_W = TMath::Sqrt((mass_proton + nu) * (mass_proton + nu) - el_change.Vect().Mag2());

        return m_ge.wght;
    }
};

int main(int argc, char *argv[]) {
    /* Faff f{"../data/Genie_gst_2000000.root"}; */
    /* f.runAnalysis(); */

    /* std::cout << "All rescs:" << std::endl; */
    /* for (auto const& [key, val] : f.all_resc_counts) { */
    /*     std::cout << "\t" << key << "\t" << val << std::endl; */
    /* } */

    /* std::cout << "Pip rescs:" << std::endl; */
    /* for (auto const& [key, val] : f.pip_resc_counts) { */
    /*     std::cout << "\t" << key << "\t" << val << std::endl; */
    /* } */

    /* return 0; */

    // Old stuff
    string input_file, output_file;
    if (argc == 2) {
        std::string arg{argv[1]};
        if (arg == "local") {
            input_file = "/home/honza/Sync/University/CurrentCourses/SHP/data/Genie_gst_2000000.root";
            output_file = "local";
        } else if (arg == "full") {
            input_file = "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_2261GeV/"
                         "apapadop_SuSav2_C12_2261GeV_master.root";
            output_file = "full";
        }
    } else {
        std::cout << "Needs an argument to specify which way to run (should be \"local\" or \"full\"), check the code"
                  << std::endl;
        return -1;
    }

    Faff2 ga{
        input_file.c_str(), ("outs/misc_" + output_file + ".root").c_str(), {}, {}, {"ALL", "QE", "RES", "DELTA1232", "RES_OTHER", "DIS"}};
    ga.runAnalysis();

    return 0;
}
