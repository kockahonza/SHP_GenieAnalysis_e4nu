#ifndef GACLAS6FINALFST_H
#define GACLAS6FINALFST_H

#include <iostream>
#include <limits>
#include <optional>
#include <tuple>
#include <utility>

#include <TH3D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector3.h>

#include "GAAutoHistograms.h"
#include "GACLAS6Common.h"
#include "misc.h"

using namespace GACLAS6Common;

/*
 * This is largely copied from GACLAS6MC, needed to get a blank as I want to alter some of the processing in passesCuts
 * without compromising the previous stuff, it may be that this is redundant and can be replaced by GACLAS6MC with the
 * right flags.
 */
class NuclearTransparencyStudies : public GAAutoHistograms {
  public:
    enum class PionType { Plus, Minus };
    enum class NucleonType { Proton, Neutron };

    enum class RunType { PrimaryState, FinalStateResc1, FinalState, DetectorLike };

    const PionType m_pi_type;
    const NucleonType m_nuc_type;
    const RunType m_run_type;

    const Target m_target;
    const BeamEnergy m_beam_energy;

    // This is a part of my "extra" fiducials for the report
    double m_p_pion_momentum_threshold{0.15};
    double m_p_photon_momentum_threshold{0.3};
    double m_p_proton_momentum_threshold{0.3};

    // Electron, pion and nucleon smearing, photons aren't smeared
    double m_smearing_reso_electron{0.005};
    double m_smearing_reso_pion{0.007};
    double m_smearing_reso_proton{0.01};

    optional<pair<Double_t, Double_t>> m_Wcut{};

  protected:
    // Initialized in constructor and then const
    Double_t m_beam_energy_val;
    TLorentzVector m_beam_V4;
    // The 4 momentum of the nucleon of the given type before the collision assuming it at rest
    TLorentzVector m_rest_nuc_V4;

    double m_el_precut_parameter1;
    double m_el_precut_parameter2;

    // Fiducials
    const unique_ptr<FiducialWrapper> m_fiducials;

  private:
    // e2a acceptance maps
    const unique_ptr<TFile> m_el_acceptance_file;  // electron
    const unique_ptr<TFile> m_p_acceptance_file;   // proton
    const unique_ptr<TFile> m_pip_acceptance_file; // pi plus
    const unique_ptr<TFile> m_pim_acceptance_file; // pi minus
    // It seems to be necessary for I suppose performance reason? that these are found beforehand and not at each call
    // to acceptanceJoined
    const unique_ptr<TH3D> m_acc_el_gen, m_acc_el_acc, m_acc_p_gen, m_acc_p_acc, m_acc_pip_gen, m_acc_pip_acc,
        m_acc_pim_gen, m_acc_pim_acc;

  private:
    // Simple flags
    bool m_found_pi;
    bool m_found_nuc;

    // Internl storage for DetectorLike particle counting
    vector<TVector3> m_det_pips;
    vector<TVector3> m_det_pims;
    vector<TVector3> m_det_protons;
    vector<TVector3> m_det_neutrons;

  protected:
    // Physical properties of the event and system to be used in passesCuts and for various properties
    TLorentzVector m_el_V4;
    TLorentzVector m_pi_V4;
    TLorentzVector m_nuc_V4;

    // Reconstructed purely from electron data
    Double_t m_reco_W;
    Double_t m_reco_Q2;
    Double_t m_reco_x;

    // Kinematics reconstruction
    // The reconstructed pion 4momentum
    TLorentzVector m_reco_pi_V4;
    // The difference of the actual and the reconstructed pion 4 momentum
    TLorentzVector m_reco_pi_diff_V4;
    // The measured 4 momentum of the original nucleon calculated by assuming 0 4 momentum change, see how it's
    // calculated dat the end of passesCuts
    TLorentzVector m_original_nuc_V4;
    // The total 4 momentum change, it is the sum of the incident electron and the hit nucleon assumed at rest, minus
    // the total 4 momentum of the electron, pion and nucleon
    TLorentzVector m_total_p_change_V4;

    // Mainly intende for final state run types
    optional<Int_t> m_resc_same;

    // These are only available with primary state runs (RunType::PrimaryState)
    Int_t m_ps_pi_resc;
    Int_t m_ps_nuc_resc;

    // Extensions to automatic TH1Fs
    map<string, AutoProperty> m_new_known_properties{
        // Electrons properties
        {"el_phi",
         {"Out electron phi [°]",
          {720, -30, 330},
          [this]() {
              double phi_deg{m_el_V4.Phi() * TMath::RadToDeg()};
              return (phi_deg < -30) ? phi_deg + 360 : phi_deg;
          }}},
        {"el_ct", {"Out electron cos theta", {500, -1, 1}, [this]() { return m_el_V4.CosTheta(); }}},
        {"el_p", {"Out electron momentum [GeV/c]", {500, 0, 3}, [this]() { return m_el_V4.P(); }}},
        {"el_E", {"Out electron energy [GeV]", {500, 0, 3}, [this]() { return m_el_V4.E(); }}},
        // Pion
        {"pi_phi",
         {"Pion phi [°]",
          {720, -30, 330},
          [this]() {
              double phi_deg{m_pi_V4.Phi() * TMath::RadToDeg()};
              return (phi_deg < -30) ? phi_deg + 360 : phi_deg;
          }}},
        {"pi_ct", {"Pion cos theta", {500, -1, 1}, [this]() { return m_pi_V4.CosTheta(); }}},
        {"pi_p", {"Pion momentum [GeV/c]", {500, 0, 3}, [this]() { return m_pi_V4.P(); }}},
        {"pi_E", {"Pion energy [GeV]", {500, 0, 3}, [this]() { return m_pi_V4.E(); }}},
        // Nucleon
        {"nuc_phi",
         {"Nucleon phi [°]",
          {720, -30, 330},
          [this]() {
              double phi_deg{m_nuc_V4.Phi() * TMath::RadToDeg()};
              return (phi_deg < -30) ? phi_deg + 360 : phi_deg;
          }}},
        {"nuc_ct", {"Nucleon cos theta", {500, -1, 1}, [this]() { return m_nuc_V4.CosTheta(); }}},
        {"nuc_p", {"Nucleon momentum [GeV/c]", {500, 0, 3}, [this]() { return m_nuc_V4.P(); }}},
        {"nuc_E", {"Nucleon energy [GeV]", {500, 0, 3}, [this]() { return m_nuc_V4.E(); }}},

        // Reconstructed properties of the event and then pi
        {"reco_W", {"Reconstructed W [GeV]", {500, 0, 4}, [this]() { return m_reco_W; }}},
        {"reco_Q2", {"Reconsotructed Q^2 [GeV]", {500, 0, 4}, [this]() { return m_reco_Q2; }}},
        {"reco_x", {"Reconstructed Bjorken x", {500, 0, 1.01}, [this]() { return m_reco_x; }}},

        {"reco_pi_phi",
         {"Kinematically reconstructed pion phi [°]",
          {720, -30, 330},
          [this]() {
              double phi_deg{m_reco_pi_V4.Phi() * TMath::RadToDeg()};
              return (phi_deg < -30) ? phi_deg + 360 : phi_deg;
          }}},
        {"reco_pi_ct",
         {"Kinematically reconstructed pion cos theta", {500, -1, 1}, [this]() { return m_reco_pi_V4.CosTheta(); }}},
        {"reco_pi_p",
         {"Kinematically reconstructed pion momentum [GeV/c]", {500, 0, 3}, [this]() { return m_reco_pi_V4.P(); }}},
        {"reco_pi_E",
         {"Kinematically reconstructed pion energy [GeV]", {500, 0, 3}, [this]() { return m_reco_pi_V4.E(); }}},

        /* {"reco_pi_diff_phi", */
        /*  {"Error of the reconstructed pion - phi [°]", */
        /*   {720, -30, 330}, */
        /*   [this]() { */
        /*       double phi_deg{m_reco_pi_diff_V4.Phi() * TMath::RadToDeg()}; */
        /*       return (phi_deg < -30) ? phi_deg + 360 : phi_deg; */
        /*   }}}, */
        /* {"reco_pi_diff_ct", */
        /*  {"Error of the reconstructed pion - cos theta", {500, -1, 1}, [this]() { return
           m_reco_pi_diff_V4.CosTheta(); }}}, */
        /* {"reco_pi_diff_p", */
        /*  {"Error of the reconstructed pion - momentum [GeV/c]", {500, 0, 3}, [this]() { return m_reco_pi_diff_V4.P();
           }}}, */
        /* {"reco_pi_diff_E", */
        /*  {"Error of the reconstructed pion - energy [GeV]", {500, -2, 2}, [this]() { return m_reco_pi_diff_V4.E();
           }}}, */

        {"reco_pi_p_test",
         {"Difference of the 3 momentum magnitudes of the reconstructed and actual pions",
          {500, -2, 2},
          [this]() { return m_reco_pi_V4.P() - m_pi_V4.P(); }}},

        {"org_nuc_phi",
         {"Measured pre-collision nucleon phi [°]",
          {720, -30, 330},
          [this]() {
              double phi_deg{m_original_nuc_V4.Phi() * TMath::RadToDeg()};
              return (phi_deg < -30) ? phi_deg + 360 : phi_deg;
          }}},
        {"org_nuc_ct",
         {"Measured pre-collision nucleon cos theta", {500, -1, 1}, [this]() { return m_original_nuc_V4.CosTheta(); }}},
        {"org_nuc_mag",
         {"Measured pre-collision nucleon 3 momentum magnitude [GeV/c]",
          {500, -1, 2},
          [this]() { return m_original_nuc_V4.P(); }}},
        {"org_nuc_E",
         {"Measured pre-collision nucleon energy [GeV/c^2]", {500, -1, 2}, [this]() { return m_original_nuc_V4.E(); }}},

        {"p_change_phi",
         {"Total momentum change phi [°]",
          {720, -30, 330},
          [this]() {
              double phi_deg{m_total_p_change_V4.Phi() * TMath::RadToDeg()};
              return (phi_deg < -30) ? phi_deg + 360 : phi_deg;
          }}},
        {"p_change_ct",
         {"Total momentum change cos theta", {500, -1, 1}, [this]() { return m_total_p_change_V4.CosTheta(); }}},
        {"p_change_mag",
         {"Total momentum change magnitude [GeV/c]", {500, -1, 2}, [this]() { return m_total_p_change_V4.P(); }}},
        {"E_change", {"Total energy change [GeV/c^2]", {500, -1, 2}, [this]() { return m_total_p_change_V4.E(); }}},

        // Rescattering
        {"resc_same",
         {"If all resc values are the same, it is that value, 0 if not",
          {11, -2.5, 8.5},
          [this]() { return m_resc_same.value_or(0); }}},
    };

  public:
    NuclearTransparencyStudies(const char *filename, const char *output_filename,
                               // Params for this analysis
                               const PionType &pi_type, const NucleonType &nuc_type,
                               const RunType &run_type = RunType::PrimaryState,
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
          m_pi_type{pi_type}, m_nuc_type{nuc_type}, m_run_type{run_type}, m_target{target}, m_beam_energy{beam_energy},
          m_fiducials{std::make_unique<FiducialWrapper>(m_beam_energy)},

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
        m_known_properties.insert(m_new_known_properties.begin(), m_new_known_properties.end());

        // Initialize precut parameters
        if (beam_energy == BeamEnergy::MeV_1161) {
            m_el_precut_parameter1 = 17;
            m_el_precut_parameter2 = 7;
            m_beam_energy_val = 1.161;
        } else if (beam_energy == BeamEnergy::MeV_2261) {
            m_el_precut_parameter1 = 16;
            m_el_precut_parameter2 = 10.5;
            m_beam_energy_val = 2.261;
        } else if (beam_energy == BeamEnergy::MeV_4461) {
            m_el_precut_parameter1 = 13.5;
            m_el_precut_parameter2 = 15;
            m_beam_energy_val = 4.461;
        }
        m_beam_V4 = TLorentzVector(0, 0, sqrt(m_beam_energy_val * m_beam_energy_val - mass_electron * mass_electron),
                                   m_beam_energy_val);

        if (m_run_type == RunType::PrimaryState) {
            m_known_properties["pi_resc"] = {
                "Pion resc code from gst", {11, -2.5, 8.5}, [this]() { return m_ps_pi_resc; }};
            m_known_properties["nuc_resc"] = {
                "Nucleon resc code from gst", {11, -2.5, 8.5}, [this]() { return m_ps_nuc_resc; }};
        }

        if (m_nuc_type == NucleonType::Proton) {
            m_rest_nuc_V4 = TLorentzVector(0, 0, 0, mass_proton);
        } else if (m_nuc_type == NucleonType::Neutron) {
            m_rest_nuc_V4 = TLorentzVector(0, 0, 0, mass_neutron);
        }
    }

    Double_t passesCuts() override;

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

#endif
