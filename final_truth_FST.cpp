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

int main(int argc, char *argv[]) {
    if (argc != 2) {
        throw std::runtime_error(
            "Needs an argument to specify which way to run (should be \"local\" or \"full\"), check the code");
    }
    Target target{Target::Fe56};
    BeamEnergy beam_energy{BeamEnergy::MeV_2261};

    string output_file{argv[1]};
    string input_file{get_data_filename(target, beam_energy, output_file)};

    /* Final1Pion1NucleonTruth ga{input_file.c_str(), */
    /*                            ("outs/final_test_" + output_file + ".root").c_str(), */
    /*                            Final1Pion1NucleonTruth::PionType::Minus, */
    /*                            Final1Pion1NucleonTruth::NucleonType::Proton, */
    /*                            Final1Pion1NucleonTruth::RunType::PrimaryState, */
    /*                            {}, */
    /*                            {}, */
    /*                            {"ALL", "QE", "DELTA1232", "RES_OTHER", "DIS"}, */
    /*                            { */
    /*                                {"pi_phi", "reco_pi_phi", {}}, */
    /*                                {"pi_ct", "reco_pi_ct", {}}, */
    /*                                {"pi_p", "reco_pi_p", {}}, */
    /*                                {"pi_E", "reco_pi_E", {}}, */
    /*                                {"pi_resc", "nuc_resc", {}}, */
    /*                            }, */
    /*                            ElectronFiducials::UseFiducials::Option1, */
    /*                            target, */
    /*                            beam_energy}; */
    /* ga.runAnalysis(); */


    bool doWcut{false};

    vector<tuple<string, Final1Pion1NucleonTruth::PionType, Final1Pion1NucleonTruth::NucleonType>> variants{
        {"_0pip1pim1p0n", Final1Pion1NucleonTruth::PionType::Minus, Final1Pion1NucleonTruth::NucleonType::Proton},
        {"_1pip0pim0p1n", Final1Pion1NucleonTruth::PionType::Plus, Final1Pion1NucleonTruth::NucleonType::Neutron},

        /* {"_1pip0pim1p0n", Final1Pion1NucleonTruth::PionType::Plus, Final1Pion1NucleonTruth::NucleonType::Proton}, */
        /* {"_0pip1pim0p1n", Final1Pion1NucleonTruth::PionType::Minus, Final1Pion1NucleonTruth::NucleonType::Neutron},
         */
    };

    string desc = doWcut ? "Wcut_" + output_file : output_file;
    if (target == Target::C12) {
        desc += "C12_";
    } else if (target == Target::Fe56) {
        desc += "Fe56_";
    } else {
        return -1;
    }

    // ps
    for (auto const &[ext, pi_t, nuc_t] : variants) {
        Final1Pion1NucleonTruth ga{input_file.c_str(),
                                   ("outs/final_ps_" + desc + ext + ".root").c_str(),
                                   pi_t,
                                   nuc_t,
                                   Final1Pion1NucleonTruth::RunType::PrimaryState,
                                   {},
                                   {},
                                   {"ALL", "QE", "DELTA1232", "RES_OTHER", "DIS"},
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
        Final1Pion1NucleonTruth ga{input_file.c_str(),
                                   ("outs/final_fsr_" + desc + ext + ".root").c_str(),
                                   pi_t,
                                   nuc_t,
                                   Final1Pion1NucleonTruth::RunType::FinalStateResc1,
                                   {},
                                   {},
                                   {"ALL", "QE", "DELTA1232", "RES_OTHER", "DIS"},
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
        Final1Pion1NucleonTruth ga{input_file.c_str(),
                                   ("outs/final_fs_" + desc + ext + ".root").c_str(),
                                   pi_t,
                                   nuc_t,
                                   Final1Pion1NucleonTruth::RunType::FinalState,
                                   {},
                                   {},
                                   {"ALL", "QE", "DELTA1232", "RES_OTHER", "DIS"},
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
