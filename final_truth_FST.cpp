#include <iostream>
#include <limits>
#include <optional>
#include <tuple>

#include <TH3D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector3.h>

#include "GenieAnalysis/GAAutoHistograms.h"
#include "GenieAnalysis/GACLAS6Common.h"
#include "GenieAnalysis/misc.h"

using namespace GACLAS6Common;

/*
 * This is largely copied from GACLAS6MC, needed to get a blank as I want to alter some of the processing in passesCuts
 * without compromising the previous stuff
 */
class Final1Pion1Nucleon : public GAAutoHistograms {
  public:
    enum class RunType { Truth, DetectorLike };
    enum class PionType { Plus, Minus };
    enum class NucleonType { Proton, Neutron };

    const Target m_target;
    const BeamEnergy m_beam_energy;

    const RunType m_run_type;
    const PionType m_pi_type;
    const NucleonType m_nuc_type;

  public:
    double m_p_pion_momentum_threshold{0.15};
    double m_p_photon_momentum_threshold{0.3};
    double m_p_proton_momentum_threshold{0.3};

    // Electron, pion and nucleon smearing, photons aren't smeared
    double m_smearing_reso_electron{0.005};
    double m_smearing_reso_pion{0.007};
    double m_smearing_reso_proton{0.01};

  protected:
    // Automatically determined parameters -- should be essentially const, but leaving mutable for easy initialization
    // and maybe someone would like to change them
    double m_precut_parameter1;
    double m_precut_parameter2;

    // Initialized in constructor and then const
    Double_t m_beam_energy_val;

  private:
    // Paths and other system stuff
    const unique_ptr<FiducialWrapper> m_fiducials;

    // e2a acceptance maps
    const unique_ptr<TFile> m_el_acceptance_file;  // electron
    const unique_ptr<TFile> m_p_acceptance_file;   // proton
    const unique_ptr<TFile> m_pip_acceptance_file; // pi plus
    const unique_ptr<TFile> m_pim_acceptance_file; // pi minus
    // It seems to be necessary for I suppose performance reason? that these are found beforehand and not at each call
    // to acceptanceJoined
    const unique_ptr<TH3D> m_acc_el_gen, m_acc_el_acc, m_acc_p_gen, m_acc_p_acc, m_acc_pip_gen, m_acc_pip_acc,
        m_acc_pim_gen, m_acc_pim_acc;

  protected:
    // Physical properties of the event and system to be used in passesCuts and for various properties
    TLorentzVector m_el_V4;

    bool m_found_pi;
    TLorentzVector m_pi_V4;
    Int_t m_pi_resc;

    bool m_found_nuc;
    TLorentzVector m_nuc_V4;
    Int_t m_nuc_resc;

    Double_t m_reco_W;
    Double_t m_reco_Q2;
    Double_t m_reco_x;

  public:
    Final1Pion1Nucleon(const char *filename, const char *output_filename,
                       // Params for this analysis
                       const PionType &pi_type, const NucleonType &nuc_type, const RunType &run_type = RunType::Truth,
                       // Setup for AutoHistograms
                       const vector<string> &stages = {}, const vector<string> &properties = {},
                       const vector<string> &types = {},
                       const vector<GAAutoHistograms::AutoVsPlot> &vs_property_plots = {},
                       // CLAS6 run params
                       const Target &target = Target::C12, const BeamEnergy &beam_energy = BeamEnergy::MeV_2261,
                       // Pass this to GenieAnalysis
                       const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name), GAAutoHistograms{filename,   output_filename, stages,
                                                                    properties, types,           vs_property_plots},
          m_target{target}, m_beam_energy{beam_energy}, m_run_type{run_type}, m_pi_type{pi_type}, m_nuc_type{nuc_type},
          m_fiducials{std::make_unique<FiducialWrapper>(m_target, m_beam_energy)},

          // Initializing acceptance map files and TH3Ds
          m_el_acceptance_file{getAcceptanceMapFile(m_target, m_beam_energy, "")},
          m_p_acceptance_file{getAcceptanceMapFile(m_target, m_beam_energy, "_p")},
          m_pip_acceptance_file{getAcceptanceMapFile(m_target, m_beam_energy, "_pip")},
          m_pim_acceptance_file{getAcceptanceMapFile(m_target, m_beam_energy, "_pim")},
          m_acc_el_gen{(TH3D *)m_el_acceptance_file->Get("Generated Particles")},
          m_acc_el_acc{(TH3D *)m_el_acceptance_file->Get("Accepted Particles")},
          m_acc_p_gen{(TH3D *)m_p_acceptance_file->Get("Generated Particles")},
          m_acc_p_acc{(TH3D *)m_p_acceptance_file->Get("Accepted Particles")},
          m_acc_pip_gen{(TH3D *)m_pip_acceptance_file->Get("Generated Particles")},
          m_acc_pip_acc{(TH3D *)m_pip_acceptance_file->Get("Accepted Particles")},
          m_acc_pim_gen{(TH3D *)m_pim_acceptance_file->Get("Generated Particles")},
          m_acc_pim_acc{(TH3D *)m_pim_acceptance_file->Get("Accepted Particles")} {
        // Initialize precut parameters
        if (beam_energy == BeamEnergy::MeV_1161) {
            m_precut_parameter1 = 17;
            m_precut_parameter2 = 7;
            m_beam_energy_val = 1.161;
        } else if (beam_energy == BeamEnergy::MeV_2261) {
            m_precut_parameter1 = 16;
            m_precut_parameter2 = 10.5;
            m_beam_energy_val = 2.261;
        } else if (beam_energy == BeamEnergy::MeV_4461) {
            m_precut_parameter1 = 13.5;
            m_precut_parameter2 = 15;
            m_beam_energy_val = 4.461;
        }
    }

    Double_t passesCuts() {
        // Electron 4 momentum, rotate to be the same as CLAS6
        m_el_V4.SetPxPyPzE(m_ge.pxl, m_ge.pyl, m_ge.pzl, m_ge.El);
        m_el_V4.SetPhi(m_el_V4.Phi() + TMath::Pi());

        // Electron theta and momentum fiducial (essentially I think) cut, the values are specifically for C12 2.261GeV
        // set by inspecting
        // https://docs.google.com/presentation/d/1ghG08JfCYXRXh6O8hcXKrhJOFxkAs_9i5ZfoIkiiEHU/edit?usp=sharing and
        // previous values
        const double theta_min{TMath::DegToRad() * (m_precut_parameter1 + m_precut_parameter2 / m_el_V4.P())};
        if ((m_el_V4.Theta() < theta_min) || (m_el_V4.Theta() > 80 * TMath::DegToRad())) {
            return 0;
        }
        // Electron fiducials
        if (!m_fiducials->electronCut(m_el_V4.Vect())) {
            return 0;
        }

        // Calculation of kinematic quantities (nu, Q2, x bjorken, q and W) -- literally taken from original though
        const TLorentzVector el_change{m_el_V4 - TLorentzVector{0, 0, m_beam_energy_val, m_beam_energy_val}};
        m_reco_Q2 = -el_change.Mag2();
        const double nu = -el_change.E();
        m_reco_x = m_reco_Q2 / (2 * mass_proton * nu);
        m_reco_W = TMath::Sqrt((mass_proton + nu) * (mass_proton + nu) - el_change.Vect().Mag2());

        m_found_pi = false;
        m_found_nuc = false;
        // Temp variables for the final state hadrons(and photon) loop, declared here for performance, these are the
        // most interesting hadrons for comparison to experiment
        double smearing;
        TVector3 V3;
        for (int i{0}; i < m_ge.nf; i++) {
            if (m_ge.pdgf[i] == 2212) { // proton
                V3.SetXYZ(m_ge.pxf[i], m_ge.pyf[i], m_ge.pzf[i]);
                V3.SetPhi(V3.Phi() + TMath::Pi());

                if (V3.Mag() < m_p_proton_momentum_threshold) {
                    continue;
                }

                if (!m_fiducials->protonCut(V3)) {
                    continue;
                }

                // If we got here, it means it qualifies as signal, if we are still looking for a proton, make this
                // found, if we are looking for neutrons or have already found a proton, the event is invalid so return
                // 0
                if ((m_nuc_type == NucleonType::Proton) && (!m_found_nuc)) {
                    m_found_nuc = true;
                    m_nuc_V4 = TLorentzVector(V3, m_ge.Ef[i]);
                    m_nuc_resc = m_ge.resc[i];
                } else {
                    return 0;
                }

            } else if (m_ge.pdgf[i] == 2112) { // neutron
                V3.SetXYZ(m_ge.pxf[i], m_ge.pyf[i], m_ge.pzf[i]);
                V3.SetPhi(V3.Phi() + TMath::Pi());

                // If we got here, it means it qualifies as signal, if we are still looking for a neutron, make this
                // found, if we are looking for protons or have already found a nucleon, the event is invalid so return
                // 0
                if ((m_nuc_type == NucleonType::Neutron) && (!m_found_nuc)) {
                    m_found_nuc = true;
                    m_nuc_V4 = TLorentzVector(V3, m_ge.Ef[i]);
                    m_nuc_resc = m_ge.resc[i];
                } else {
                    return 0;
                }

            } else if (m_ge.pdgf[i] == 211) { // pi+
                V3.SetXYZ(m_ge.pxf[i], m_ge.pyf[i], m_ge.pzf[i]);
                V3.SetPhi(V3.Phi() + TMath::Pi());

                if (V3.Mag() < m_p_pion_momentum_threshold) {
                    continue;
                }

                if (!m_fiducials->piAndPhotonCuts(V3, FiducialWrapper::PiPhotonId::Plus)) {
                    continue;
                }

                // If we got here, it means it qualifies as signal, if we are still looking for a pi+, make this found,
                // if we are looking for pi- or have already found a suitable pion, the event is invalid so return 0
                if ((m_pi_type == PionType::Plus) && (!m_found_pi)) {
                    m_found_pi = true;
                    m_pi_V4 = TLorentzVector(V3, m_ge.Ef[i]);
                    m_pi_resc = m_ge.resc[i];
                } else {
                    return 0;
                }

            } else if (m_ge.pdgf[i] == -211) { // pi-
                V3.SetXYZ(m_ge.pxf[i], m_ge.pyf[i], m_ge.pzf[i]);
                V3.SetPhi(V3.Phi() + TMath::Pi());

                if (V3.Mag() < m_p_pion_momentum_threshold) {
                    continue;
                }

                if (!m_fiducials->piAndPhotonCuts(V3, FiducialWrapper::PiPhotonId::Minus)) {
                    continue;
                }

                // If we got here, it means it qualifies as signal, if we are still looking for a pi-, make this found,
                // if we are looking for pi+ or have already found a suitable pion, the event is invalid so return 0
                if ((m_pi_type == PionType::Minus) && (!m_found_pi)) {
                    m_found_pi = true;
                    m_pi_V4 = TLorentzVector(V3, m_ge.Ef[i]);
                    m_pi_resc = m_ge.resc[i];
                } else {
                    return 0;
                }

            } else if (m_ge.pdgf[i] == 22) { // photon
                // Nothing now
            } else {
                /* std::cout << "extra FS particle -- " << m_ge.pdgf[i] << std::endl; */
            }
        }

        if ((!m_found_pi) || (!m_found_nuc)) {
            return 0;
        }

        return m_ge.wght;
    };

    double acceptanceJoined(const double &p, const double &cos_theta, double phi, const unique_ptr<TH3D> &generated,
                            const unique_ptr<TH3D> &accepted);

    double electronAcceptance(const double &p, const double &cos_theta, const double &phi) {
        return acceptanceJoined(p, cos_theta, phi, m_acc_el_gen, m_acc_el_acc);
    }

    double protonAcceptance(const double &p, const double &cos_theta, const double &phi) {
        return acceptanceJoined(p, cos_theta, phi, m_acc_p_gen, m_acc_p_acc);
    }

    double piPlusAcceptance(const double &p, const double &cos_theta, const double &phi) {
        return acceptanceJoined(p, cos_theta, phi, m_acc_pip_gen, m_acc_pip_acc);
    }

    double piMinusAcceptance(const double &p, const double &cos_theta, const double &phi) {
        return acceptanceJoined(p, cos_theta, phi, m_acc_pim_gen, m_acc_pim_acc);
    }

    static unique_ptr<TFile> getAcceptanceMapFile(const Target &target, const BeamEnergy &beam_energy,
                                                  const string &ending) {
        string target_str, beam_energy_str;
        if (target == Target::C12) {
            target_str = "12C";
        } else if (target == Target::Fe56) {
            target_str = "12C"; // There's no dedicated Fe file and original used 12C for anything except He
        }

        if (beam_energy == BeamEnergy::MeV_1161) {
            beam_energy_str = "1_161";
        } else if (beam_energy == BeamEnergy::MeV_2261) {
            beam_energy_str = "2_261";
        } else if (beam_energy == BeamEnergy::MeV_4461) {
            beam_energy_str = "4_461";
        }

        return unique_ptr<TFile>(TFile::Open(
            ("original/e2a_maps/e2a_maps_" + target_str + "_E_" + beam_energy_str + ending + ".root").c_str(), "READ"));
    };
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

    Final1Pion1Nucleon gamp{input_file.c_str(), ("final_" + output_file + "0pip1pim1p0n.root").c_str(),
                            Final1Pion1Nucleon::PionType::Minus, Final1Pion1Nucleon::NucleonType::Proton};
    gamp.runAnalysis();
}
