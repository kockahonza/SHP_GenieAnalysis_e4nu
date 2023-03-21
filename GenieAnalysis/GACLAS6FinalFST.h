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
class ElectronFiducials : public GAAutoHistograms {
  public:
    enum class PionType { Plus, Minus };
    enum class NucleonType { Proton, Neutron };

    const Target m_target;
    const BeamEnergy m_beam_energy;

  public:
    optional<pair<Double_t, Double_t>> m_Wcut{};

  protected:
    // Automatically determined parameters -- should be essentially const, but leaving mutable for easy initialization
    // and maybe someone would like to change them
    double m_el_precut_parameter1;
    double m_el_precut_parameter2;

    // Initialized in constructor and then const
    Double_t m_beam_energy_val;

    // Fiducials, accesible to inheriting classes
    const unique_ptr<FiducialWrapper> m_fiducials;

  protected:
    // Physical properties of the event and system to be used in passesCuts and for various properties
    TLorentzVector m_el_V4;

    Double_t m_reco_W;
    Double_t m_reco_Q2;
    Double_t m_reco_x;

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
        {"el_ct", {"Out electron cos theta", {720, -1, 1}, [this]() { return m_el_V4.CosTheta(); }}},
        {"el_p", {"Out electron momentum [GeV/c]", {720, 0, 3}, [this]() { return m_el_V4.P(); }}},
        {"el_E", {"Out electron energy [GeV]", {720, 0, 3}, [this]() { return m_el_V4.Energy(); }}},
    };

  public:
    ElectronFiducials(const char *filename, const char *output_filename,
                      // Setup for AutoHistograms
                      const vector<string> &stages = {}, const vector<string> &properties = {},
                      const vector<string> &types = {},
                      const vector<GAAutoHistograms::AutoVsPlot> &vs_property_plots = {},
                      // CLAS6 run params
                      const Target &target = Target::C12, const BeamEnergy &beam_energy = BeamEnergy::MeV_2261,
                      // Pass this to GenieAnalysis
                      const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name),
          GAAutoHistograms{filename, output_filename, stages, properties, types, vs_property_plots}, m_target{target},
          m_beam_energy{beam_energy}, m_fiducials{std::make_unique<FiducialWrapper>(m_target, m_beam_energy)} {
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
    }

    Double_t passesCuts() {
        // Electron 4 momentum, rotate to be the same as CLAS6
        m_el_V4.SetPxPyPzE(m_ge.pxl, m_ge.pyl, m_ge.pzl, m_ge.El);
        m_el_V4.SetPhi(m_el_V4.Phi() + TMath::Pi());

        // Electron theta and momentum fiducial (essentially I think) cut, the values are specifically for C12 2.261GeV
        // set by inspecting
        // https://docs.google.com/presentation/d/1ghG08JfCYXRXh6O8hcXKrhJOFxkAs_9i5ZfoIkiiEHU/edit?usp=sharing and
        // previous values
        const double theta_min{TMath::DegToRad() * (m_el_precut_parameter1 + m_el_precut_parameter2 / m_el_V4.P())};
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

        if (m_Wcut) {
            if ((m_reco_W < m_Wcut->first) || (m_reco_W > m_Wcut->second)) {
                std::cout << "Wcut trigerred" << std::endl;
                return 0;
            }
        }

        return m_ge.wght;
    };
};

class Final1Pion1NucleonTruth : public ElectronFiducials {
  public:
    enum class PionType { Plus, Minus };
    enum class NucleonType { Proton, Neutron };

    const PionType m_pi_type;
    const NucleonType m_nuc_type;

  public:
    double m_p_pion_momentum_threshold{0.15};
    double m_p_photon_momentum_threshold{0.3};
    double m_p_proton_momentum_threshold{0.3};

  protected:
    // Physical properties of the event and system to be used in passesCuts and for various properties
    bool m_found_pi;
    TLorentzVector m_pi_V4;
    Int_t m_pi_resc;

    bool m_found_nuc;
    TLorentzVector m_nuc_V4;
    Int_t m_nuc_resc;

    // Extensions to automatic TH1Fs
    map<string, AutoProperty> m_new_known_properties{
        {"pi_phi",
         {"Pion phi [°]",
          {720, -30, 330},
          [this]() {
              double phi_deg{m_pi_V4.Phi() * TMath::RadToDeg()};
              return (phi_deg < -30) ? phi_deg + 360 : phi_deg;
          }}},
        {"pi_ct", {"Pion cos theta", {720, -1, 1}, [this]() { return m_pi_V4.CosTheta(); }}},
        {"pi_p", {"Pion momentum [GeV/c]", {720, 0, 3}, [this]() { return m_pi_V4.P(); }}},
        {"pi_E", {"Pion energy [GeV]", {720, 0, 3}, [this]() { return m_pi_V4.Energy(); }}},
        {"pi_resc", {"Pion resc code from gst", {11, -2.5, 8.5}, [this]() { return m_pi_resc; }}},

        {"nuc_phi",
         {"Nucleon phi [°]",
          {720, -30, 330},
          [this]() {
              double phi_deg{m_nuc_V4.Phi() * TMath::RadToDeg()};
              return (phi_deg < -30) ? phi_deg + 360 : phi_deg;
          }}},
        {"nuc_ct", {"Nucleon cos theta", {720, -1, 1}, [this]() { return m_nuc_V4.CosTheta(); }}},
        {"nuc_p", {"Nucleon momentum [GeV/c]", {720, 0, 3}, [this]() { return m_nuc_V4.P(); }}},
        {"nuc_E", {"Nucleon energy [GeV]", {720, 0, 3}, [this]() { return m_nuc_V4.Energy(); }}},
        {"nuc_resc", {"Nucleon resc code from gst", {11, -2.5, 8.5}, [this]() { return m_nuc_resc; }}},

        // Physical properties of the event as a whole
        {"reco_W", {"Reconstructed W [GeV]", {1000, 0, 4}, [this]() { return m_reco_W; }}},
        {"reco_Q2", {"Reconsotructed Q^2 [GeV]", {1000, 0, 4}, [this]() { return m_reco_Q2; }}},
        {"reco_x", {"Reconstructed Bjorken x", {1000, 0, 1.01}, [this]() { return m_reco_x; }}},
    };

  public:
    Final1Pion1NucleonTruth(const char *filename, const char *output_filename,
                            // Params for this analysis
                            const PionType &pi_type, const NucleonType &nuc_type,
                            // Setup for AutoHistograms
                            const vector<string> &stages = {}, const vector<string> &properties = {},
                            const vector<string> &types = {},
                            const vector<GAAutoHistograms::AutoVsPlot> &vs_property_plots = {},
                            // CLAS6 run params
                            const Target &target = Target::C12, const BeamEnergy &beam_energy = BeamEnergy::MeV_2261,
                            // Pass this to GenieAnalysis
                            const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name), ElectronFiducials{filename, output_filename,   stages, properties,
                                                                     types,    vs_property_plots, target, beam_energy},
          m_pi_type{pi_type}, m_nuc_type{nuc_type} {
        m_known_properties.insert(m_new_known_properties.begin(), m_new_known_properties.end());
    }

    Double_t passesCuts() {
        Double_t weight{ElectronFiducials::passesCuts()};
        if (weight == 0) {
            return 0;
        }

        m_found_pi = false;
        m_found_nuc = false;
        // Temp variables for the primary state hadrons(and photon) loop, declared here for performance, these are the
        // most interesting hadrons for comparison to experiment
        double smearing;
        TVector3 V3;
        for (int i{0}; i < m_ge.ni; i++) {
            if (m_ge.pdgi[i] == 2212) { // proton
                V3.SetXYZ(m_ge.pxi[i], m_ge.pyi[i], m_ge.pzi[i]);
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
                    m_nuc_V4 = TLorentzVector(V3, m_ge.Ei[i]);
                    m_nuc_resc = m_ge.resc[i];
                } else {
                    return 0;
                }

            } else if (m_ge.pdgi[i] == 2112) { // neutron
                V3.SetXYZ(m_ge.pxi[i], m_ge.pyi[i], m_ge.pzi[i]);
                V3.SetPhi(V3.Phi() + TMath::Pi());

                // If we got here, it means it qualifies as signal, if we are still looking for a neutron, make this
                // found, if we are looking for protons or have already found a nucleon, the event is invalid so return
                // 0
                if ((m_nuc_type == NucleonType::Neutron) && (!m_found_nuc)) {
                    m_found_nuc = true;
                    m_nuc_V4 = TLorentzVector(V3, m_ge.Ei[i]);
                    m_nuc_resc = m_ge.resc[i];
                } else {
                    return 0;
                }

            } else if (m_ge.pdgi[i] == 211) { // pi+
                V3.SetXYZ(m_ge.pxi[i], m_ge.pyi[i], m_ge.pzi[i]);
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
                    m_pi_V4 = TLorentzVector(V3, m_ge.Ei[i]);
                    m_pi_resc = m_ge.resc[i];
                } else {
                    return 0;
                }

            } else if (m_ge.pdgi[i] == -211) { // pi-
                V3.SetXYZ(m_ge.pxi[i], m_ge.pyi[i], m_ge.pzi[i]);
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
                    m_pi_V4 = TLorentzVector(V3, m_ge.Ei[i]);
                    m_pi_resc = m_ge.resc[i];
                } else {
                    return 0;
                }

            } else if (m_ge.pdgi[i] == 22) { // photon
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
};

#endif
