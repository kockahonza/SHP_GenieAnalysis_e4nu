#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <utility>

#include <TH3D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector3.h>

#include "GenieAnalysis/GACLAS6Common.h"
#include "GenieAnalysis/GACLAS6FinalFST.h"

int do_long(
    string local_or_full, Target target, BeamEnergy beam_energy, vector<string> types,
    vector<tuple<string, NuclearTransparencyStudies::PionType, NuclearTransparencyStudies::NucleonType>> variants,
    bool doWcut) {
    string input_file{get_data_filename(target, beam_energy, local_or_full)};

    string desc = doWcut ? "Wcut_" + local_or_full : local_or_full;
    if (target == Target::C12) {
        desc += "_C12";
    } else if (target == Target::Fe56) {
        desc += "_Fe56";
    } else {
        return -1;
    }

    // ps
    for (auto const &[ext, pi_t, nuc_t] : variants) {
        NuclearTransparencyStudies ga{input_file.c_str(),
                                      ("outs/final_ps_" + desc + ext + ".root").c_str(),
                                      pi_t,
                                      nuc_t,
                                      NuclearTransparencyStudies::RunType::PrimaryState,
                                      {},
                                      {},
                                      types,
                                      {
                                          {"pi_phi", "reco_pi_phi", {}},
                                          {"pi_ct", "reco_pi_ct", {}},
                                          {"pi_p", "reco_pi_p", {}},
                                          {"pi_E", "reco_pi_E", {}},
                                          {"pi_phi", "reco_pi_p_test", {}},
                                          {"pi_ct", "reco_pi_p_test", {}},
                                          {"pi_p", "reco_pi_p_test", {}},
                                          {"pi_E", "reco_pi_p_test", {}},
                                          {"pi_resc", "nuc_resc", {}},
                                      }};
        if (doWcut) {
            ga.m_Wcut = {0, 1.4};
        }
        ga.runAnalysis();
    }

    // fsr
    for (auto const &[ext, pi_t, nuc_t] : variants) {
        NuclearTransparencyStudies ga{input_file.c_str(),
                                      ("outs/final_fsr_" + desc + ext + ".root").c_str(),
                                      pi_t,
                                      nuc_t,
                                      NuclearTransparencyStudies::RunType::FinalStateResc1,
                                      {},
                                      {},
                                      types,
                                      {
                                          {"pi_phi", "reco_pi_phi", {}},
                                          {"pi_ct", "reco_pi_ct", {}},
                                          {"pi_p", "reco_pi_p", {}},
                                          {"pi_E", "reco_pi_E", {}},
                                          {"pi_phi", "reco_pi_p_test", {}},
                                          {"pi_ct", "reco_pi_p_test", {}},
                                          {"pi_p", "reco_pi_p_test", {}},
                                          {"pi_E", "reco_pi_p_test", {}},
                                      }};
        if (doWcut) {
            ga.m_Wcut = {0, 1.4};
        }
        ga.runAnalysis();
    }

    // fs
    for (auto const &[ext, pi_t, nuc_t] : variants) {
        NuclearTransparencyStudies ga{input_file.c_str(),
                                      ("outs/final_fs_" + desc + ext + ".root").c_str(),
                                      pi_t,
                                      nuc_t,
                                      NuclearTransparencyStudies::RunType::FinalState,
                                      {},
                                      {},
                                      types,
                                      {
                                          {"pi_phi", "reco_pi_phi", {}},
                                          {"pi_ct", "reco_pi_ct", {}},
                                          {"pi_p", "reco_pi_p", {}},
                                          {"pi_E", "reco_pi_E", {}},
                                          {"pi_phi", "reco_pi_p_test", {}},
                                          {"pi_ct", "reco_pi_p_test", {}},
                                          {"pi_p", "reco_pi_p_test", {}},
                                          {"pi_E", "reco_pi_p_test", {}},
                                          {"resc_same", "reco_pi_p_test", {}},
                                          {"resc_same", "p_change_mag", {}},
                                          {"resc_same", "E_change", {}},
                                      }};
        if (doWcut) {
            ga.m_Wcut = {0, 1.4};
        }
        ga.runAnalysis();
    }

    // dl
    for (auto const &[ext, pi_t, nuc_t] : variants) {
        NuclearTransparencyStudies ga{input_file.c_str(),
                                      ("outs/final_dl_" + desc + ext + ".root").c_str(),
                                      pi_t,
                                      nuc_t,
                                      NuclearTransparencyStudies::RunType::DetectorLike,
                                      {},
                                      {},
                                      types,
                                      {
                                          {"pi_phi", "reco_pi_phi", {}},
                                          {"pi_ct", "reco_pi_ct", {}},
                                          {"pi_p", "reco_pi_p", {}},
                                          {"pi_E", "reco_pi_E", {}},
                                          {"pi_phi", "reco_pi_p_test", {}},
                                          {"pi_ct", "reco_pi_p_test", {}},
                                          {"pi_p", "reco_pi_p_test", {}},
                                          {"pi_E", "reco_pi_p_test", {}},
                                          {"resc_same", "reco_pi_p_test", {}},
                                          {"resc_same", "p_change_mag", {}},
                                          {"resc_same", "E_change", {}},
                                      }};
        if (doWcut) {
            ga.m_Wcut = {0, 1.4};
        }
        ga.runAnalysis();
    }

    return 0;
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        throw std::runtime_error("Needs arguments, check the code");
    }
    string local_or_full{argv[1]};
    string target_str{argv[2]};
    string do_Wcut_str{argv[3]};

    Target target;
    if (target_str == "C") {
        target == Target::C12;
    } else if (target_str == "Fe") {
        target == Target::Fe56;
    } else {
        throw std::runtime_error("Target specifier not valid");
    }

    bool doWcut;
    if (do_Wcut_str == "0") {
        doWcut = false;
    } else if (do_Wcut_str == "1") {
        doWcut = true;
    } else {
        throw std::runtime_error("Wcut flag not valid");
    }

    BeamEnergy beam_energy{BeamEnergy::MeV_2261};
    vector<string> types{"ALL", "QE", "DELTA1232", "RES_OTHER", "DIS"};

    vector<tuple<string, NuclearTransparencyStudies::PionType, NuclearTransparencyStudies::NucleonType>> variants{
        {"_0pip1pim1p0n", NuclearTransparencyStudies::PionType::Minus, NuclearTransparencyStudies::NucleonType::Proton},
        {"_1pip0pim0p1n", NuclearTransparencyStudies::PionType::Plus, NuclearTransparencyStudies::NucleonType::Neutron},

        /* {"_1pip0pim1p0n", NuclearTransparencyStudies::PionType::Plus,
           NuclearTransparencyStudies::NucleonType::Proton}, */
        /* {"_0pip1pim0p1n", NuclearTransparencyStudies::PionType::Minus,
         * NuclearTransparencyStudies::NucleonType::Neutron},
         */
    };

    return do_long(local_or_full, target, beam_energy, types, variants, doWcut);
}
