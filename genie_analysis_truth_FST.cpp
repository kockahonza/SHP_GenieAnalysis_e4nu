#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>
#include <limits>
#include <tuple>

#include "GenieAnalysis/GAAutoHistograms.h"
#include "GenieAnalysis/GACLAS6MC.h"
#include "GenieAnalysis/misc.h"

class GA1Pion1Nucleon : public GAAutoHistograms {
  public:
    enum class RunType { Truth, DetectorLike };

    enum class PionType { Plus, Minus };

    enum class NucleonType { Proton, Neutron };

  private:
    const RunType m_run_type;
    const PionType m_pi_type;
    const NucleonType m_nuc_type;

  protected:
    TLorentzVector m_pion_V4;
    TVector3 m_pion_V3;
    Double_t m_pion_acceptance;

    TLorentzVector m_nucleon_V4;
    TVector3 m_nucleon_V3;
    Double_t m_nucleon_acceptance;

    map<string, AutoProperty> m_new_known_properties{};

  public:
    GA1Pion1Nucleon(const char *filename, const char *output_filename, const RunType &run_type, const PionType &pi_type,
                    const NucleonType &nuc_type, const vector<string> &stages = {},
                    const vector<string> &properties = {}, const vector<string> &types = {},
                    const vector<GAAutoHistograms::AutoVsPlot> &vs_property_plots = {},
                    /* const Target &target = Target::C12, const BeamEnergy &beam_energy = BeamEnergy::MeV_2261, */
                    const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name), GAAutoHistograms{filename, output_filename},
          m_run_type{run_type}, m_pi_type{pi_type}, m_nuc_type{nuc_type} {
        m_known_properties.insert(m_new_known_properties.begin(), m_new_known_properties.end());
    }

    Double_t passesCuts() override {
        /* Double_t weight{GACLAS6MC::passesCuts()}; */

        /* if (weight == 0) { */
        /*     return 0; */
        /* } */

        /* if ((m_reconstructed_W < m_p_Wcut_min) || (m_reconstructed_W > m_p_Wcut_max)) { */
        /*     return 0; */
        /* } */

        /* size_t pip_num, pim_num, p_num, n_num; */
        /* if (m_run_type == RunType::Detector) { */
        /*     pip_num = m_passed_pi_plus.size(); */
        /*     pim_num = m_passed_pi_minus.size(); */
        /*     p_num = m_passed_protons.size(); */
        /*     n_num = m_fs_neutrons.size(); */
        /* } else if (m_run_type == RunType::FinalState) { */
        /*     pip_num = m_fs_pi_plus.size(); */
        /*     pim_num = m_fs_pi_minus.size(); */
        /*     p_num = m_fs_protons.size(); */
        /*     n_num = m_fs_neutrons.size(); */
        /* } */

        /* /1* std::cout << pip_num << "," << pim_num << "," << p_num << "," << n_num << std::endl; *1/ */

        /* if ((m_pi_type == PionType::Plus) && (pip_num == 1) && (pim_num == 0)) { */
        /*     if (m_run_type == RunType::Detector) { */
        /*         std::tie(m_pion_V4, m_pion_V3, m_pion_acceptance) = m_passed_pi_plus[0]; */
        /*     } else if (m_run_type == RunType::FinalState) { */
        /*         std::tie(m_pion_V4, m_pion_V3) = m_fs_pi_plus[0]; */
        /*         m_pion_acceptance = 1; */
        /*     } */
        /* } else if ((m_pi_type == PionType::Minus) && (pip_num == 0) && (pim_num == 1)) { */
        /*     if (m_run_type == RunType::Detector) { */
        /*         std::tie(m_pion_V4, m_pion_V3, m_pion_acceptance) = m_passed_pi_minus[0]; */
        /*     } else if (m_run_type == RunType::FinalState) { */
        /*         std::tie(m_pion_V4, m_pion_V3) = m_fs_pi_minus[0]; */
        /*         m_pion_acceptance = 1; */
        /*     } */
        /* } else { */
        /*     return 0; */
        /* } */

        /* if ((m_nuc_type == NucleonType::Proton) && (p_num == 1) && (n_num == 0)) { */
        /*     if (m_run_type == RunType::Detector) { */
        /*         std::tie(m_nucleon_V4, m_nucleon_V3, m_nucleon_acceptance) = m_passed_protons[0]; */
        /*     } else if (m_run_type == RunType::FinalState) { */
        /*         std::tie(m_nucleon_V4, m_nucleon_V3) = m_fs_protons[0]; */
        /*         m_nucleon_acceptance = 1; */
        /*     } */
        /* } else if ((m_nuc_type == NucleonType::Neutron) && (p_num == 0) && (n_num == 1)) { */
        /*     // Neutrons aren't detected so always use fs neutrons */
        /*     std::tie(m_nucleon_V4, m_nucleon_V3) = m_fs_neutrons[0]; */
        /*     m_nucleon_acceptance = 1; */
        /* } else { */
        /*     return 0; */
        /* } */

        /* if (m_use_acceptances) { */
        /*     return weight * m_pion_acceptance * m_nucleon_acceptance; */
        /* } else { */
        /*     return weight; */
        /* } */
        return 1;
    }
};

int main(int argc, char *argv[]) {
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
        return -1;
    }

    vector<GAAutoHistograms::AutoVsPlot> vs_property_plots{
        {"el_p", "pi_p", {}},
        {"el_p", "nuc_p", {}},
        {"pi_p", "nuc_p", {}},
        {"el_phi", "pi_phi", {}},
        {"el_phi", "nuc_phi", {}},
        {"pi_phi", "nuc_phi", {}},
        {"el_cos_theta", "pi_cos_theta", {}},
        {"el_cos_theta", "nuc_cos_theta", {}},
        {"pi_cos_theta", "nuc_cos_theta", {}},

        {"el_p", "ptest", {}},
        {"el_phi", "ptest", {}},
        {"el_cos_theta", "ptest", {}},

        {"pi_p", "ptest", {}},
        {"pi_phi", "ptest", {}},
        {"pi_cos_theta", "ptest", {}},

        {"ni", "resc_same", {}},
        {"nfp", "resc_same", {}},
        {"nfpip", "resc_same", {}},
        {"nfpim", "resc_same", {}},
        {"nf", "resc_same", {}},
    };

    GA1Pion1Nucleon gapnlc{input_file.c_str(),
                           (output_file + "_n_1pip0pim0p1n_Wlucas.root").c_str(),
                           GA1Pion1Nucleon::RunType::DetectorLike,
                           GA1Pion1Nucleon::PionType::Plus,
                           GA1Pion1Nucleon::NucleonType::Neutron,
                           {},
                           {},
                           {},
                           vs_property_plots};
    gapnlc.runAnalysis();
}
