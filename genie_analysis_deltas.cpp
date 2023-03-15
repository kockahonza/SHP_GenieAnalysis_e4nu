#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>
#include <tuple>

#include "GenieAnalysis/GenieAnalysisDeltaStudies.h"
#include "GenieAnalysis/misc.h"

/**
 * This is for studies that have exactly one charged Pion, does all the stuff from DeltaStudies and add's the pion
 * filter, can check for nucleon numbers and uses pion acceptance for the weight, also adds the pions properties
 */
class GenieAnalysis1Pion : public GenieAnalysisDeltaStudies {
  public:
    enum class PionType { Minus, Plus, Either };

    const PionType m_pion_type;
    const optional<int> m_proton_count;
    const optional<int> m_neutron_count;

  protected:
    TLorentzVector m_pion_V4;
    TVector3 m_pion_V3;
    Double_t m_pion_acceptance;

  protected:
    map<string, AutoProperty> m_new_known_properties{
        {"pi_phi",
         {"Pion phi [°]",
          {720, -30, 330},
          [this]() {
              double phi_deg{m_pion_V3.Phi() * TMath::RadToDeg()};
              return (phi_deg < -30) ? phi_deg + 360 : phi_deg;
          }}},
        {"pi_cos_theta", {"Pion cos theta", {720, -1, 1}, [this]() { return m_pion_V3.CosTheta(); }}},
        {"pi_p", {"Pion momentum [GeV/c]", {720, 0, 3}, [this]() { return m_pion_V4.P(); }}},
        {"pi_E", {"Pion energy [GeV]", {720, 0, 3}, [this]() { return m_pion_V4.E(); }}},
        {"pi_acceptance", {"Pion acceptance weight", {100, 0, 1}, [this]() { return m_pion_acceptance; }}}};

  public:
    GenieAnalysis1Pion(const char *filename, const char *output_filename, PionType pion_type = PionType::Either,
                       optional<int> proton_count = {}, optional<int> neutron_count = {},
                       const vector<string> &stages = {}, const vector<string> &properties = {},
                       const vector<string> &types = {},
                       const vector<GenieAnalysisAutoHistograms::AutoVsPlot> &vs_property_plots = {},
                       const Target &target = GenieAnalysisDeltaStudies::Target::C12,
                       const BeamEnergy &beam_energy = GenieAnalysisDeltaStudies::BeamEnergy::MeV_2261)

        : GenieAnalysis(filename), GenieAnalysisDeltaStudies(filename, output_filename, stages, properties, types,
                                                             vs_property_plots, target, beam_energy),
          m_pion_type{pion_type}, m_proton_count{proton_count}, m_neutron_count{neutron_count} {
        m_known_properties.insert(m_new_known_properties.begin(), m_new_known_properties.end());
    }

    Double_t passesCuts() override {
        Double_t weight{GenieAnalysisDeltaStudies::passesCuts()};
        if (weight == 0) {
            return 0;
        }

        if ((m_pion_type != PionType::Plus) && (m_passed_pi_plus.size() == 0) && (m_passed_pi_minus.size() == 1)) {
            std::tie(m_pion_V4, m_pion_V3, m_pion_acceptance) = m_passed_pi_minus[0];
        } else if ((m_pion_type != PionType::Minus) && (m_passed_pi_minus.size() == 0) &&
                   (m_passed_pi_plus.size() == 1)) {
            std::tie(m_pion_V4, m_pion_V3, m_pion_acceptance) = m_passed_pi_plus[0];
        } else {
            return 0;
        }

        if (m_proton_count && (m_proton_count != m_fs_number_of_protons)) {
            return 0;
        }

        if (m_neutron_count && (m_neutron_count != m_fs_number_of_neutrons)) {
            return 0;
        }

        if (m_do_pion_acceptance) {
            return weight * m_pion_acceptance;
        } else {
            return weight;
        }
    }
};

int main(int argc, char *argv[]) {
    vector<string> properties{"W",
                              "wght",
                              "el_phi",
                              "el_cos_theta",
                              "el_p",
                              "el_E",
                              "el_acceptance",
                              "pi_phi",
                              "pi_cos_theta",
                              "pi_p",
                              "pi_E",
                              "pi_acceptance",
                              "reco_W",
                              "bjorken_x",
                              "fs_num_protons",
                              "fs_num_neutrons",
                              "fs_num_pip",
                              "fs_num_pim"};
    vector<string> types{"ALL", "QE", "RES_ALL", "DELTA1232", "DIS"};

    string input_file, output_file;

    if (argc == 2) {
        std::string arg{argv[1]};
        if (arg == "local") {
            input_file = "/home/honza/Sync/University/CurrentCourses/SHP/data/Genie_gst_2000000.root";

            GenieAnalysis1Pion ga{input_file.c_str(),
                                  "output_local_pip.root",
                                  GenieAnalysis1Pion::PionType::Plus,
                                  {},
                                  {},
                                  {"nocut"},
                                  properties,
                                  types};
            ga.runAnalysis();

            GenieAnalysis1Pion ga2{input_file.c_str(),
                                   "output_local_pim.root",
                                   GenieAnalysis1Pion::PionType::Minus,
                                   {},
                                   {},
                                   {"nocut"},
                                   properties,
                                   types};
            ga2.runAnalysis();

            /* GenieAnalysis1Pion ga3{input_file.c_str(), */
            /*                        "output_local_pie.root", */
            /*                        {"nocut"}, */
            /*                        properties, */
            /*                        types, */
            /*                        GenieAnalysis1Pion::PionType::Either}; */
            /* ga3.runAnalysis(); */

        } else if (arg == "full") {
            input_file = "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_2261GeV/"
                         "apapadop_SuSav2_C12_2261GeV_master.root";

            /* GenieAnalysis1Pion gap{ */
            /*     input_file.c_str(), "output_full_pip.root", {}, properties, types,
             * GenieAnalysis1Pion::PionType::Plus}; */
            /* gap.runAnalysis(); */

            /* GenieAnalysis1Pion gam{ */
            /*     input_file.c_str(), "output_full_pim.root", {}, properties, types,
             * GenieAnalysis1Pion::PionType::Minus}; */
            /* gam.runAnalysis(); */

            /* GenieAnalysis1Pion gae{input_file.c_str(), */
            /*                        "output_full_pie.root", */
            /*                        {}, */
            /*                        properties, */
            /*                        types, */
            /*                        GenieAnalysis1Pion::PionType::Either}; */
            /* gae.runAnalysis(); */
        }
    } else {
        std::cout << "Needs an argument to specify which way to run (should be \"local\" or \"full\"), check the code"
                  << std::endl;
        return -1;
    }

    /* ga.m_do_precuts = false; */
    /* ga.m_do_electron_fiducials = false; */
    /* ga.m_do_pion_fiducials = false; */
    /* ga.m_do_photon_fiducials = false; */

    /* ga.m_p_pion_momentum_threshold = 0; */

    return 0;
}
