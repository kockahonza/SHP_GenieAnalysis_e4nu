#ifndef GENIE_ANALYSIS_ORIGINAL_CUTS_H
#define GENIE_ANALYSIS_ORIGINAL_CUTS_H

#include <optional>

#include <TH3D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>

#include "GenieAnalysisAuto.h"
#include "misc.h"

class FiducialWrapper;

using std::unique_ptr, std::optional;

/**
 * Main class for Delta studies analysis, mainly tailored for Genie data tuned for the CLAS6 runs. It takes care of the
 * electron smearing, cuts and so on. It also gathers all pions, gets the numbers of selected hadrons (pi+-, p, n) in
 * the primary and final states. The electron properties are made available here any single pion properties should be
 * done in inheriting classes.
 *
 * The class collects no data, `GenieAnalysisDeltaStudies` also adds the auto histogram functionality and many available
 * properties, most other classes then inherit from that.
 */
class GenieAnalysisDeltaStudiesCuts : public virtual GenieAnalysis {
  public:
    // Configuration options for the major target/energy runs
    enum class Target { C12, Fe56 };
    // If the ~1GeV energy is to be added, smearing resolutions were tripled!, do't forget to add that
    enum class BeamEnergy { MeV_1161, MeV_2261, MeV_4461 };

    // These should not be changed once initialized
    const Target m_target;
    const BeamEnergy m_beam_energy;

  public:
    // Some parameters to be set by directly accessing them, otherwise the constructor would be massive,
    // the values here are used for the CLAS6 C12 2GeV GENIE data.
    bool m_do_precuts{true};
    bool m_do_sectors{false};
    bool m_do_radiation_check{false};

    bool m_do_electron_fiducials{true};
    bool m_do_pion_fiducials{true};
    bool m_do_proton_fiducials{true};
    bool m_do_photon_fiducials{true};

    bool m_do_electron_acceptance{true};
    bool m_do_pion_acceptance{true};
    bool m_do_proton_acceptance{true};

    bool m_gather_fs_particles{false};

    double m_p_pion_momentum_threshold{0.15};
    double m_p_photon_momentum_threshold{0.3};
    double m_p_proton_momentum_threshold{0.3};
    double m_p_radiation_photon_angle{40};
    double m_p_radiation_photon_phi_diff{30};

    // Electron, pion and nucleon smearing, photons aren't smeared
    double m_smearing_reso_electron{0.005};
    double m_smearing_reso_pion{0.007};
    double m_smearing_reso_proton{0.01};

  protected:
    // Automatically determined parameters -- should be essentially const, but leaving mutable for easy initialization
    // and maybe someone would like to change them
    double m_precut_parameter1;
    double m_precut_parameter2;

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
    // Mostly physical properties of the event and system to be used in passesCuts and fot
    // various properties
    // will be set in constructor and then const
    Double_t m_beam_energy_val;

    TVector3 m_smeared_el_V3; // Smeared and rotated by pi
    TLorentzVector m_smeared_el_V4;
    Double_t m_electron_acceptance_weight;

    Double_t m_reconstructed_W;
    Double_t m_bjorken_x;

    // Collection of particles detectable by the CLAS6 detector -- these come from final state particles, passing
    // fiducials and momentum thresholds and have acceptances (except photons) Only charged pions, protons and photons
    // are detected as that's all the particles we have fiducials for. The tuples contain the smeared 4 momentum and 3
    // momentum and calculated acceptance.
    vector<tuple<TLorentzVector, TVector3, double>> m_passed_pi_plus;
    vector<tuple<TLorentzVector, TVector3, double>> m_passed_pi_minus;
    vector<tuple<TLorentzVector, TVector3, double>> m_passed_protons;
    vector<tuple<TLorentzVector, TVector3>> m_passed_photons; // photons don't have acceptances
                                                              //

    // Also gather charged pions and nucleons in the final state, to perform ?truth level? studies and also we don't
    // have detected neutrons so this way we can still look at all 1pion1nucleon events not just 1pim1p. (these
    // naturally don't have acceptances)
    vector<tuple<TLorentzVector, TVector3>> m_fs_pi_plus;
    vector<tuple<TLorentzVector, TVector3>> m_fs_pi_minus;
    vector<tuple<TLorentzVector, TVector3>> m_fs_protons;
    vector<tuple<TLorentzVector, TVector3>> m_fs_neutrons;

    // "Primary" state properties, used as in the gst documentation, these are particles after interaction but before
    // FSI
    Int_t m_ps_number_of_pi_plus;
    Int_t m_ps_number_of_pi_minus;
    Int_t m_ps_number_of_protons;
    Int_t m_ps_number_of_neutrons;
    Int_t m_ps_number_of_photons;

    // The genie final state particle counters, I think this is often referred to as truth
    Int_t m_fs_number_of_pi_plus;
    Int_t m_fs_number_of_pi_minus;
    Int_t m_fs_number_of_protons;
    Int_t m_fs_number_of_neutrons;
    Int_t m_fs_number_of_photons;

  public:
    GenieAnalysisDeltaStudiesCuts(const char *filename,
                                  // Select run
                                  const Target &target = Target::C12,
                                  const BeamEnergy &beam_energy = BeamEnergy::MeV_2261,
                                  // Pass this to GenieAnalysis
                                  const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name), m_target{target}, m_beam_energy{beam_energy},

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

    Double_t passesCuts();

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

/**
 * A very simple wrapper around the Fiducial class, which is taken exactly as in original, so that it can be referred
 * to in a more modern way. No need to pas strings to specify what is to be done in each call.
 */
class FiducialWrapper {
  public:
    enum class PiPhotonId : int { Minus = -1, Photon = 0, Plus = 1 };

    const GenieAnalysisDeltaStudiesCuts::Target m_target;
    const GenieAnalysisDeltaStudiesCuts::BeamEnergy m_beam_energy;
    const string m_target_str;
    const string m_beam_energy_str;

  private:
    Fiducial m_fiducial;

  public:
    FiducialWrapper(const GenieAnalysisDeltaStudiesCuts::Target &target,
                    const GenieAnalysisDeltaStudiesCuts::BeamEnergy &beam_energy)
        : m_target(target), m_beam_energy(beam_energy), m_target_str{targetStr(m_target)},
          m_beam_energy_str{beamEnergyStr(m_beam_energy)}, m_fiducial{} {

        m_fiducial.InitPiMinusFit(m_beam_energy_str);
        m_fiducial.InitEClimits();

        // The first value is torusCurrent, values taken from original
        m_fiducial.SetConstants(m_beam_energy == GenieAnalysisDeltaStudiesCuts::BeamEnergy::MeV_1161 ? 750 : 2250,
                                m_target_str, {{"1161", 1.161}, {"2261", 2.261}, {"4461", 4.461}});
        m_fiducial.SetFiducialCutParameters(m_beam_energy_str);
    }

    bool electronCut(const TVector3 &momentum_V3) { return m_fiducial.EFiducialCut(m_beam_energy_str, momentum_V3); }

    bool protonCut(const TVector3 &momentum_V3) { return m_fiducial.PFiducialCut(m_beam_energy_str, momentum_V3); }

    bool piAndPhotonCuts(const TVector3 &momentum_V3, const PiPhotonId &which_particle) {
        return m_fiducial.Pi_phot_fid_united(m_beam_energy_str, momentum_V3, static_cast<int>(which_particle));
    }

    static string targetStr(const GenieAnalysisDeltaStudiesCuts::Target &target) {
        if (target == GenieAnalysisDeltaStudiesCuts::Target::C12) {
            return "12C";
        } else if (target == GenieAnalysisDeltaStudiesCuts::Target::Fe56) {
            return "12C"; // There's no dedicated Fe file and original used 12C for anything except He
        }
        throw "This should not happen";
    }

    static string beamEnergyStr(const GenieAnalysisDeltaStudiesCuts::BeamEnergy &beam_energy) {
        if (beam_energy == GenieAnalysisDeltaStudiesCuts::BeamEnergy::MeV_1161) {
            return "1161";
        } else if (beam_energy == GenieAnalysisDeltaStudiesCuts::BeamEnergy::MeV_2261) {
            return "2261";
        } else if (beam_energy == GenieAnalysisDeltaStudiesCuts::BeamEnergy::MeV_4461) {
            return "4461";
        }
        throw "This should not happen";
    }
};

/*
 * Adds automatic histogram generation from `GenieAnalysisAutoHistograms` and a bunch of suitable properties to
 * `GenieAnalysisDeltaStudiesCuts`.
 */
class GenieAnalysisDeltaStudies : public GenieAnalysisDeltaStudiesCuts, public GenieAnalysisAutoHistograms {
  protected:
    // Extensions to automatic TH1Fs
    map<string, AutoProperty> m_new_known_properties{
        // Electrons properties
        {"el_phi",
         {"Out electron phi [°]",
          {720, -30, 330},
          [this]() {
              double phi_deg{m_smeared_el_V3.Phi() * TMath::RadToDeg()};
              return (phi_deg < -30) ? phi_deg + 360 : phi_deg;
          }}},
        {"el_cos_theta", {"Out electron cos theta", {720, -1, 1}, [this]() { return m_smeared_el_V3.CosTheta(); }}},
        {"el_p", {"Out electron momentum [GeV/c]", {720, 0, 3}, [this]() { return m_smeared_el_V4.P(); }}},
        {"el_E", {"Out electron energy [GeV]", {720, 0, 3}, [this]() { return m_smeared_el_V4.Energy(); }}},
        {"el_acceptance",
         {"Out electron acceptance weight", {100, 0, 1}, [this]() { return m_electron_acceptance_weight; }}},

        // Physical properties of the event as a whole
        {"reco_W", {"Reconstructed W [GeV]", {1000, 0, 4}, [this]() { return m_reconstructed_W; }}},
        {"bjorken_x", {"Bjorken x", {1000, 0, 1.01}, [this]() { return m_bjorken_x; }}},

        // Passed particle properties, this mainly applies to pi+- and protons which there may be multiple of and we
        // have acceptances fo5
        {"passed_num_pip",
         {"Number of detected pi plus", {6, -0.5, 5.5}, [this]() { return m_passed_pi_plus.size(); }}},
        {"passed_num_pim",
         {"Number of detected pi minus", {6, -0.5, 5.5}, [this]() { return m_passed_pi_minus.size(); }}},
        {"passed_num_protons",
         {"Number of detected protons", {6, -0.5, 5.5}, [this]() { return m_passed_protons.size(); }}},
        {"passed_num_photons",
         {"Number of detected photons", {6, -0.5, 5.5}, [this]() { return m_passed_photons.size(); }}},

        // "Truth" data follows, there's a lot of it and it's quite repetitive
        // Primary state (pre FSI) data
        {"ps_num_pip",
         {"Number of pi plus in the primary state", {6, -0.5, 5.5}, [this]() { return m_ps_number_of_pi_plus; }}},
        {"ps_num_pim",
         {"Number of pi minus in the primary state", {6, -0.5, 5.5}, [this]() { return m_ps_number_of_pi_minus; }}},
        {"ps_num_protons",
         {"Number of protons in the primary state", {6, -0.5, 5.5}, [this]() { return m_ps_number_of_protons; }}},
        {"ps_num_neutrons",
         {"Number of neutrons in the primary state", {6, -0.5, 5.5}, [this]() { return m_ps_number_of_neutrons; }}},
        {"ps_num_photons",
         {"Number of photons in the primary state", {6, -0.5, 5.5}, [this]() { return m_ps_number_of_photons; }}},
        // Final state (post FSI) data, truth - before accounting for detector with fiducials and acceptances
        {"fs_num_pip",
         {"Number of pi plus in the final state", {6, -0.5, 5.5}, [this]() { return m_fs_number_of_pi_plus; }}},
        {"fs_num_pim",
         {"Number of pi minus in the final state", {6, -0.5, 5.5}, [this]() { return m_fs_number_of_pi_minus; }}},
        {"fs_num_protons",
         {"Number of protons in the final state", {6, -0.5, 5.5}, [this]() { return m_fs_number_of_protons; }}},
        {"fs_num_neutrons",
         {"Number of neutrons in the final state", {6, -0.5, 5.5}, [this]() { return m_fs_number_of_neutrons; }}},
        {"fs_num_photons",
         {"Number of photons in the final state", {6, -0.5, 5.5}, [this]() { return m_fs_number_of_photons; }}},
    };

  public:
    GenieAnalysisDeltaStudies(const char *filename, const char *output_filename,
                              // Specify the analysis - which stages, properties and types to do histograms for
                              const vector<string> &stages = {}, const vector<string> &properties = {},
                              const vector<string> &types = {},
                              const vector<GenieAnalysisAutoHistograms::AutoVsPlot> &vs_property_plots = {},
                              // Select run
                              const Target &target = Target::C12, const BeamEnergy &beam_energy = BeamEnergy::MeV_2261,
                              // Pass this to GenieAnalysis
                              const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name), GenieAnalysisDeltaStudiesCuts(filename, target, beam_energy),
          GenieAnalysisAutoHistograms(filename, output_filename, stages, properties, types, vs_property_plots) {
        m_known_properties.insert(m_new_known_properties.begin(), m_new_known_properties.end());
    }
};

/**
 * Continues on from DeltaStudies, adds specific cuts for the number of protons, neutrons and charged pions.
 * It is a very simple wrapper that might be good for quick tests, that's why I keep it here.
 */
class GenieAnalysisPiNucleonCounts : public GenieAnalysisDeltaStudies {
  public:
    const optional<int> m_pi_plus_count;
    const optional<int> m_pi_minus_count;
    const optional<int> m_proton_count;
    const optional<int> m_neutron_count;

  public:
    GenieAnalysisPiNucleonCounts(const char *filename, const char *output_filename, optional<int> pi_plus_count = {},
                                 optional<int> pi_minus_count = {}, optional<int> proton_count = {},
                                 optional<int> neutron_count = {}, const vector<string> &stages = {},
                                 const vector<string> &properties = {}, const vector<string> &types = {},
                                 const vector<GenieAnalysisAutoHistograms::AutoVsPlot> &vs_property_plots = {},
                                 const Target &target = GenieAnalysisDeltaStudies::Target::C12,
                                 const BeamEnergy &beam_energy = GenieAnalysisDeltaStudies::BeamEnergy::MeV_2261)

        : GenieAnalysis(filename), GenieAnalysisDeltaStudies(filename, output_filename, stages, properties, types,
                                                             vs_property_plots, target, beam_energy),
          m_pi_plus_count{pi_plus_count}, m_pi_minus_count{pi_minus_count}, m_proton_count{proton_count},
          m_neutron_count{neutron_count} {
        /* m_known_properties.insert(m_new_known_properties.begin(), m_new_known_properties.end()); */
    }

    Double_t passesCuts() override {
        Double_t weight{GenieAnalysisDeltaStudies::passesCuts()};
        if (weight == 0) {
            return 0;
        }

        if (m_pi_plus_count && (m_pi_plus_count != m_passed_pi_plus.size())) {
            return 0;
        }
        if (m_pi_minus_count && (m_pi_minus_count != m_passed_pi_minus.size())) {
            return 0;
        }
        if (m_proton_count && (m_proton_count != m_fs_number_of_protons)) {
            return 0;
        }
        if (m_neutron_count && (m_neutron_count != m_fs_number_of_neutrons)) {
            return 0;
        }

        return weight;
    }
};

#endif
