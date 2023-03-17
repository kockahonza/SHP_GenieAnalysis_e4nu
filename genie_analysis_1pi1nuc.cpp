#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>
#include <tuple>

#include "GenieAnalysis/GenieAnalysisDeltaStudies.h"
#include "GenieAnalysis/misc.h"

class GenieAnalysis1Pion1Nucleon : public GenieAnalysisDeltaStudies {
  public:
    enum class RunType { Detector, FinalState };

    enum class PionType { Plus, Minus };

    enum class NucleonType { Proton, Neutron };

  private:
    const RunType m_run_type;
    const PionType m_pi_type;
    const NucleonType m_nuc_type;

  protected:
    TLorentzVector m_pion_V4;
    TVector3 m_pion_V3;
    /* Double_t m_pion_acceptance; */

    TLorentzVector m_nucleon_V4;
    TVector3 m_nucleon_V3;
    /* Double_t m_nucleon_acceptance; */

  public:
    GenieAnalysis1Pion1Nucleon(const char *filename, const char *output_filename, const RunType &run_type,
                               const PionType &pi_type, const NucleonType &nuc_type, const vector<string> &stages = {},
                               const vector<string> &properties = {}, const vector<string> &types = {},
                               const vector<GenieAnalysisAutoHistograms::AutoVsPlot> &vs_property_plots = {},
                               const Target &target = Target::C12, const BeamEnergy &beam_energy = BeamEnergy::MeV_2261,
                               const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name), GenieAnalysisDeltaStudies{filename, output_filename,
                                                                             stages,   properties,
                                                                             types,    vs_property_plots,
                                                                             target,   beam_energy},
          m_run_type{run_type}, m_pi_type{pi_type}, m_nuc_type{nuc_type} {
        m_gather_fs_particles = true;
    }

    Double_t passesCuts() override {
        Double_t weight{GenieAnalysisDeltaStudies::passesCuts()};

        if (weight == 0) {
            return 0;
        }

        size_t pip_num, pim_num, p_num, n_num;
        if (m_run_type == RunType::Detector) {
            pip_num = m_passed_pi_plus.size();
            pim_num = m_passed_pi_minus.size();
            p_num = m_passed_protons.size();
            n_num = m_fs_neutrons.size();
        } else if (m_run_type == RunType::FinalState) {
            pip_num = m_fs_pi_plus.size();
            pim_num = m_fs_pi_minus.size();
            p_num = m_fs_protons.size();
            n_num = m_fs_neutrons.size();
        }

        if ((m_pi_type == PionType::Plus) && (pip_num == 1) && (pim_num == 0)) {
            if (m_run_type == RunType::Detector) {
                std::tie(m_pion_V4, m_pion_V3, std::ignore) = m_passed_pi_plus[0];
            } else if (m_run_type == RunType::FinalState) {
                std::tie(m_pion_V4, m_pion_V3) = m_fs_pi_plus[0];
            }
        } else if ((m_pi_type == PionType::Minus) && (pip_num == 0) && (pim_num == 1)) {
            if (m_run_type == RunType::Detector) {
                std::tie(m_pion_V4, m_pion_V3, std::ignore) = m_passed_pi_minus[0];
            } else if (m_run_type == RunType::FinalState) {
                std::tie(m_pion_V4, m_pion_V3) = m_fs_pi_minus[0];
            }
        } else {
            return 0;
        }

        if ((m_nuc_type == NucleonType::Proton) && (p_num == 1) && (n_num == 0)) {
            if (m_run_type == RunType::Detector) {
                std::tie(m_nucleon_V4, m_nucleon_V3, std::ignore) = m_passed_protons[0];
            } else if (m_run_type == RunType::FinalState) {
                std::tie(m_nucleon_V4, m_nucleon_V3) = m_fs_protons[0];
            }
        } else if ((m_nuc_type == NucleonType::Neutron) && (p_num == 0) && (n_num == 1)) {
            // Neutrons aren't detected so always use fs neutrons
            std::tie(m_nucleon_V4, m_nucleon_V3) = m_fs_neutrons[0];
        } else {
            return 0;
        }

        return weight;
    }
};

int simple_pi_nucleon_counts(int argc, char *argv[]) {
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

    GenieAnalysis1Pion1Nucleon gapp{input_file.c_str(), (output_file + "_n_1pip0pim1p0n.root").c_str(),
                                    GenieAnalysis1Pion1Nucleon::RunType::Detector,
                                    GenieAnalysis1Pion1Nucleon::PionType::Plus,
                                    GenieAnalysis1Pion1Nucleon::NucleonType::Proton};
    gapp.runAnalysis();

    GenieAnalysis1Pion1Nucleon gapn{input_file.c_str(), (output_file + "_n_1pip0pim0p1n.root").c_str(),
                                    GenieAnalysis1Pion1Nucleon::RunType::Detector,
                                    GenieAnalysis1Pion1Nucleon::PionType::Plus,
                                    GenieAnalysis1Pion1Nucleon::NucleonType::Neutron};
    gapn.runAnalysis();

    GenieAnalysis1Pion1Nucleon gamp{input_file.c_str(), (output_file + "_n_0pip1pim1p0n.root").c_str(),
                                    GenieAnalysis1Pion1Nucleon::RunType::Detector,
                                    GenieAnalysis1Pion1Nucleon::PionType::Minus,
                                    GenieAnalysis1Pion1Nucleon::NucleonType::Proton};
    gamp.runAnalysis();

    GenieAnalysis1Pion1Nucleon gamn{input_file.c_str(), (output_file + "_n_0pip1pim0p1n.root").c_str(),
                                    GenieAnalysis1Pion1Nucleon::RunType::Detector,
                                    GenieAnalysis1Pion1Nucleon::PionType::Minus,
                                    GenieAnalysis1Pion1Nucleon::NucleonType::Neutron};
    gamn.runAnalysis();

    return 0;
}

int main(int argc, char *argv[]) { return simple_pi_nucleon_counts(argc, argv); }
