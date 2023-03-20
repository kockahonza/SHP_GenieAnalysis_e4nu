#ifndef GENIE_ANALYSIS_GACLAS6DATA_H
#define GENIE_ANALYSIS_GACLAS6DATA_H

#include <optional>

#include <TH3D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>

#include "GAAutoHistograms.h"
#include "misc.h"

class FiducialWrapper;

using std::unique_ptr, std::optional;

class GACLAS6DataPrepare : public virtual GenieAnalysis {
  public:
    // Configuration options for the major target/energy runs
    enum class Target { C12, Fe56 };
    // If the ~1GeV energy is to be added, smearing resolutions were tripled!, do't forget to add that
    enum class BeamEnergy { MeV_1161, MeV_2261, MeV_4461 };

    // These should not be changed once initialized
    const Target m_target;
    const BeamEnergy m_beam_energy;

  protected:
    // Mostly physical properties of the event and system to be used in passesCuts and for various properties
    // will be set in constructor and then const
    Double_t m_beam_energy_val;

    TVector3 m_el_V3; // Smeared and rotated by pi
    TLorentzVector m_el_V4;

    Double_t m_reconstructed_W;
    Double_t m_reconstructed_Q2;
    Double_t m_bjorken_x;

    vector<tuple<TLorentzVector, TVector3>> m_pi_plus;
    vector<tuple<TLorentzVector, TVector3>> m_pi_minus;
    vector<tuple<TLorentzVector, TVector3>> m_protons;

  public:
    GACLAS6DataPrepare(const char *filename,
                       // Select run
                       const Target &target = Target::C12, const BeamEnergy &beam_energy = BeamEnergy::MeV_2261,
                       // Pass this to GenieAnalysis
                       const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name), m_target{target}, m_beam_energy{beam_energy} {}

    Double_t passesCuts();
};

/*
 */
class GACLAS6Data : public GACLAS6DataPrepare, public GAAutoHistograms {
  protected:
    // Extensions to automatic TH1Fs
    map<string, AutoProperty> m_new_known_properties{
        // Electron properties
        {"el_phi",
         {"Out electron phi [°]",
          {720, -30, 330},
          [this]() {
              double phi_deg{m_el_V3.Phi() * TMath::RadToDeg()};
              return (phi_deg < -30) ? phi_deg + 360 : phi_deg;
          }}},
        {"el_cos_theta", {"Out electron cos theta", {720, -1, 1}, [this]() { return m_el_V3.CosTheta(); }}},
        {"el_p", {"Out electron momentum [GeV/c]", {720, 0, 3}, [this]() { return m_el_V4.P(); }}},
        {"el_E", {"Out electron energy [GeV]", {720, 0, 3}, [this]() { return m_el_V4.Energy(); }}},

        // Physical properties of the event as a whole
        {"reco_W", {"Reconstructed W [GeV]", {1000, 0, 4}, [this]() { return m_reconstructed_W; }}},
        {"reco_Q2", {"Reconsotructed Q^2 [GeV]", {1000, 0, 4}, [this]() { return m_reconstructed_Q2; }}},
        {"bjorken_x", {"Bjorken x", {1000, 0, 1.01}, [this]() { return m_bjorken_x; }}},

        {"num_pip", {"Number of detected pi plus", {6, -0.5, 5.5}, [this]() { return m_pi_plus.size(); }}},
        {"num_pim", {"Number of detected pi minus", {6, -0.5, 5.5}, [this]() { return m_pi_minus.size(); }}},
        {"num_protons", {"Number of detected protons", {6, -0.5, 5.5}, [this]() { return m_protons.size(); }}},
    };

  public:
    GACLAS6Data(const char *filename, const char *output_filename,
                // Specify the analysis - which stages, properties and types to do histograms for
                const vector<string> &stages = {}, const vector<string> &properties = {},
                const vector<string> &types = {}, const vector<GAAutoHistograms::AutoVsPlot> &vs_property_plots = {},
                // Select run
                const Target &target = Target::C12, const BeamEnergy &beam_energy = BeamEnergy::MeV_2261,
                // Pass this to GenieAnalysis
                const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name), GACLAS6DataPrepare(filename, target, beam_energy),
          GAAutoHistograms(filename, output_filename, stages, properties, types, vs_property_plots) {
        m_known_properties.insert(m_new_known_properties.begin(), m_new_known_properties.end());
    }
};

#endif
