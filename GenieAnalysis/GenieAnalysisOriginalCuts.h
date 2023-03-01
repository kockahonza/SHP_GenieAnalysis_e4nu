#ifndef GENIE_ANALYSIS_ORIGINAL_CUTS_H
#define GENIE_ANALYSIS_ORIGINAL_CUTS_H

#include <TH3D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>
#include <memory>

#include "GenieAnalysisAuto.h"
#include "misc.h"

class FiducialWrapper;

using std::unique_ptr;

// So far I'm only focusing on exactly the stuff used with the command line arguments in the README,
// specifically this means 2.261GeV beam on C12 target and more -- not quite anymore
class GenieAnalysisOriginalCuts : public GenieAnalysisAutoTH1Fs {
  public:
    // Configuration options for the major target/energy runs
    enum class Target { C12, Fe56 };
    enum class BeamEnergy { MeV_1161, MeV_2261, MeV_4461 };

    const Target m_target;
    const BeamEnergy m_beam_energy;

    const bool m_do_precuts;
    const bool m_do_electron_fiducials;
    const bool m_do_sectors;

  protected:
    // Parameters -- should be essentially const, but leaving mutable for easy initialization and maybe someone would
    // like to change them
    double m_smearing_reso_el; // smearing for the electrons
    double m_smearing_reso_pi; // smearing for pions, executive decision by Larry (28.08.19)
    double m_smearing_reso_p;  // smearing for the proton -- comments from original

    double m_precut_parameter1;
    double m_precut_parameter2;

  private:
    // Paths and other system stuff

    const unique_ptr<FiducialWrapper> m_fiducials;

    // e2a acceptance maps
    const unique_ptr<TFile> m_el_acceptance_file;
    const unique_ptr<TFile> m_p_acceptance_file;
    const unique_ptr<TFile> m_pip_acceptance_file;
    const unique_ptr<TFile> m_pim_acceptance_file;
    // It seems to be necessary for I suppose performance reason? that these are found beforehand and not at each call
    // to acceptanceJoined
    const unique_ptr<TH3D> m_acc_el_gen, m_acc_el_acc, m_acc_p_gen, m_acc_p_acc, m_acc_pip_gen, m_acc_pip_acc,
        m_acc_pim_gen, m_acc_pim_acc;

  protected:
    // Properties of the loaded event to be set and used in passesCuts and new properties
    // electron
    TVector3 m_smeared_el_V3; // Smeared and rotated by pi
    TLorentzVector m_smeared_el_V4;
    double m_electron_acceptance_weight;

    // hadrons -- all of these only contain information on particles passing relevant cuts (pions need momentum above
    // 0.15 for example)
    vector<tuple<TLorentzVector, TVector3, double>>
        m_passed_pi_plus; // tuple has the smeared 4 momentum, smeared 3 momentum and the pions calculated acceptance
    vector<tuple<TLorentzVector, TVector3, double>>
        m_passed_pi_minus; // tuple has the smeared 4 momentum, smeared 3 momentum and the pions calculated acceptance
    vector<tuple<TLorentzVector, TVector3, double>>
        m_passed_photons; // tuple has the smeared 4 momentum, smeared 3 momentum and the pions calculated acceptance

    // Extensions to automatic TH1Fs
    map<string, AutoProperty> m_new_known_properties{
        {"el_phi",
         {"Out electron phi",
          {720, -30, 330},
          [this]() {
              double phi_deg{m_smeared_el_V3.Phi() * TMath::RadToDeg()};
              return (phi_deg < -30) ? phi_deg + 360 : phi_deg;
          }}},
        {"el_cos_theta", {"Out electron cos theta", {720, -1, 1}, [this]() { return m_smeared_el_V3.CosTheta(); }}},
        {"el_p", {"Out electron momentum", {720, 0, 3}, [this]() { return m_smeared_el_V4.P(); }}},
        {"el_E", {"Out electron energy", {720, 0, 3}, [this]() { return m_smeared_el_V4.Energy(); }}},
        {"el_acceptance",
         {"Out electron acceptance weight", {100, 0, 1}, [this]() { return m_electron_acceptance_weight; }}}};

  public:
    GenieAnalysisOriginalCuts(const char *filename, const char *output_filename,
                              // Select run
                              const Target &target, const BeamEnergy &beam_energy,
                              // Specify the analysis - which stages, properties and types to do histograms for
                              const vector<string> &stages, const vector<string> &properties,
                              const vector<string> &types,
                              // Some flags about which cuts to use
                              const bool &do_precuts = true, const bool &do_electron_fiducials = true,
                              const bool &do_sectors = false,
                              // Smearing parameters
                              const double &smearing_reso_el = 0.005, const double &smearing_reso_pi = 0.007,
                              const double &smearing_reso_p = 0.01,
                              // Pass this to GenieAnalysis
                              const char *gst_ttree_name = "gst")
        : GenieAnalysisAutoTH1Fs(filename, output_filename, stages, properties, types, gst_ttree_name),
          m_target{target}, m_beam_energy{beam_energy},

          m_do_precuts{do_precuts}, m_do_electron_fiducials{do_electron_fiducials}, m_do_sectors{do_sectors},
          m_smearing_reso_el{smearing_reso_el}, m_smearing_reso_pi{smearing_reso_pi},
          m_smearing_reso_p{smearing_reso_p},

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
          m_acc_pim_acc{(TH3D *)m_pim_acceptance_file->Get("Accepted Particles")}

    {
        m_known_properties.insert(m_new_known_properties.begin(), m_new_known_properties.end());

        // Initialize precut parameters and also multiply smearing resos by 3 for the ~1GeV as it was in the original
        if (beam_energy == BeamEnergy::MeV_1161) {
            m_precut_parameter1 = 17;
            m_precut_parameter2 = 7;

            m_smearing_reso_el *= 3;
            m_smearing_reso_pi *= 3;
            m_smearing_reso_p *= 3;
        } else if (beam_energy == BeamEnergy::MeV_2261) {
            m_precut_parameter1 = 16;
            m_precut_parameter2 = 10.5;
        } else if (beam_energy == BeamEnergy::MeV_4461) {
            m_precut_parameter1 = 13.5;
            m_precut_parameter2 = 15;
        }
    }

    Double_t passesCuts();

    double acceptanceJoined(const double &p, const double &cos_theta, double phi, const unique_ptr<TH3D> &generated,
                            const unique_ptr<TH3D> &accepted);

    double electronAcceptance(const double &p, const double &cos_theta, const double &phi) {
        return acceptanceJoined(p, cos_theta, phi, m_acc_el_gen, m_acc_el_acc);
    }

    double photonAcceptance(const double &p, const double &cos_theta, const double &phi) {
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

// A very simple wrapper around the Fiducial class, taken exactly as in original, so that it can be referred to in a
// more modern way
class FiducialWrapper {
  public:
    enum class PiPhotonId : int { Minus = -1, Photon = 0, Plus = 1 };

    const GenieAnalysisOriginalCuts::Target m_target;
    const GenieAnalysisOriginalCuts::BeamEnergy m_beam_energy;
    const string m_target_str;
    const string m_beam_energy_str;

  private:
    Fiducial m_fiducial;

  public:
    FiducialWrapper(const GenieAnalysisOriginalCuts::Target &target,
                    const GenieAnalysisOriginalCuts::BeamEnergy &beam_energy)
        : m_target(target), m_beam_energy(beam_energy), m_target_str{targetStr(m_target)},
          m_beam_energy_str{beamEnergyStr(m_beam_energy)}, m_fiducial{} {

        m_fiducial.InitPiMinusFit(m_beam_energy_str);
        m_fiducial.InitEClimits();

        // The first value is torusCurrent, values taken from original
        m_fiducial.SetConstants(m_beam_energy == GenieAnalysisOriginalCuts::BeamEnergy::MeV_1161 ? 750 : 2250,
                                m_target_str, {{"1161", 1.161}, {"2261", 2.261}, {"4461", 4.461}});
        m_fiducial.SetFiducialCutParameters(m_beam_energy_str);
    }

    bool electronCut(const TVector3 &momentum_V3) { return m_fiducial.EFiducialCut(m_beam_energy_str, momentum_V3); }

    bool piAndPhotonCuts(const TVector3 &momentum_V3, const PiPhotonId &which_particle) {
        return m_fiducial.Pi_phot_fid_united(m_beam_energy_str, momentum_V3, static_cast<int>(which_particle));
    }

    static string targetStr(const GenieAnalysisOriginalCuts::Target &target) {
        if (target == GenieAnalysisOriginalCuts::Target::C12) {
            return "12C";
        } else if (target == GenieAnalysisOriginalCuts::Target::Fe56) {
            return "12C"; // There's no dedicated Fe file and original used 12C for anything except He
        }
        throw "This should not happen";
    }

    static string beamEnergyStr(const GenieAnalysisOriginalCuts::BeamEnergy &beam_energy) {
        if (beam_energy == GenieAnalysisOriginalCuts::BeamEnergy::MeV_1161) {
            return "1161";
        } else if (beam_energy == GenieAnalysisOriginalCuts::BeamEnergy::MeV_2261) {
            return "2261";
        } else if (beam_energy == GenieAnalysisOriginalCuts::BeamEnergy::MeV_4461) {
            return "4461";
        }
        throw "This should not happen";
    }
};

#endif
