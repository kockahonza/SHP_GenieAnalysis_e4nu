#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>
#include <tuple>

#include "GenieAnalysis/GenieAnalysisOriginalCuts.h"
#include "GenieAnalysis/misc.h"

// And specifically only pi- and pi+ for now, same as original code
class GenieAnalysis1Pion : public GenieAnalysisOriginalCuts {
  private:
    int m_pion_charge;
    /* bool m_pion_plus; */
    /* bool m_pion_minus; */
    TLorentzVector m_pion_V4;
    TVector3 m_pion_V3;
    Double_t m_pion_acceptance;

  protected:
    map<string, AutoProperty> m_new_known_properties{
        {"pi_phi",
         {"Pion phi [°]",
          {720, -30, 330},
          [this]() {
              double phi_deg{m_pion_V3.Phi() * TMath::RadToDeg()};
              return (phi_deg < -30) ? phi_deg + 360 : phi_deg;
          }}},
        {"pi_cos_theta", {"Pion cos theta", {720, -1, 1}, [this]() { return m_pion_V3.CosTheta(); }}},
        {"pi_p", {"Pion momentum [GeV/c]", {720, 0, 3}, [this]() { return m_pion_V4.P(); }}},
        {"pi_E", {"Pion energy [GeV]", {720, 0, 3}, [this]() { return m_pion_V4.E(); }}},
        {"pi_acceptance", {"Pion acceptance weight", {100, 0, 1}, [this]() { return m_pion_acceptance; }}}};

    /* map<string, AutoType> m_new_known_types{ */
    /*     {"PIP", {"All events", [this]() { return m_pion_plus; }}}, */
    /*     {"QE_PIP", {"Quasi-Elastic events", [this]() { return m_ge.qel && m_pion_plus; }}}, */
    /*     {"RES_ALL_PIP", {"Resonant events", [this]() { return m_ge.res && m_pion_plus; }}}, */
    /*     {"DELTA1232_PIP", {"Resonant events with a Delta1232", [this]() { return (m_ge.res && (m_ge.resid == 0) &&
     * m_pion_plus); }}}, */
    /*     {"DIS_PIP", {"Deep-inelastic events", [this]() { return m_ge.dis && m_pion_plus; }}}, */
    /*     {"PIM", {"All events", [this]() { return m_pion_minus; }}}, */
    /*     {"QE_PIM", {"Quasi-Elastic events", [this]() { return m_ge.qel && m_pion_minus; }}}, */
    /*     {"RES_ALL_PIM", {"Resonant events", [this]() { return m_ge.res && m_pion_minus; }}}, */
    /*     {"DELTA1232_PIM", {"Resonant events with a Delta1232", [this]() { return (m_ge.res && (m_ge.resid == 0) &&
     * m_pion_minus); }}}, */
    /*     {"DIS_PIM", {"Deep-inelastic events", [this]() { return m_ge.dis && m_pion_minus; }}}, */
    /* }; */

  public:
    GenieAnalysis1Pion(const char *filename, const char *output_filename, const Target &target,
                       const BeamEnergy &beam_energy, const vector<string> &stages, const vector<string> &properties,
                       const vector<string> &types, const char *gst_ttree_name = "gst")
        : GenieAnalysisOriginalCuts(filename, output_filename, target, beam_energy, stages, properties, types) {
        m_known_properties.insert(m_new_known_properties.begin(), m_new_known_properties.end());
        /* m_known_types.insert(m_new_known_types.begin(), m_new_known_types.end()); */
    }

    Double_t passesCuts() override {
        Double_t weight{GenieAnalysisOriginalCuts::passesCuts()};
        if (weight == 0) {
            return 0;
        }

        const size_t num_pi_minus{m_passed_pi_minus.size()};
        const size_t num_pi_plus{m_passed_pi_plus.size()};

        if ((num_pi_minus + num_pi_plus) == 1) {
            if (num_pi_minus == 1) {
                m_pion_charge = -1;
                /* m_pion_plus = false; */
                /* m_pion_minus = true; */
                std::tie(m_pion_V4, m_pion_V3, m_pion_acceptance) = m_passed_pi_minus[0];
                useEntryAtStage("PIP", weight * m_pion_acceptance);
            } else if (num_pi_plus == 1) {
                m_pion_charge = +1;
                /* m_pion_plus = true; */
                /* m_pion_minus = false; */
                std::tie(m_pion_V4, m_pion_V3, m_pion_acceptance) = m_passed_pi_plus[0];
                useEntryAtStage("PIM", weight * m_pion_acceptance);
            }

            /* std::cout << weight << ", " << m_pion_acceptance << ", " << m_pion_charge << std::endl; */
            return weight * m_pion_acceptance;
        } else {
            return 0;
        }
    }
};

int main(int argc, char *argv[]) {
    GenieAnalysis1Pion ga{gst_path_jan,
                          "output_jan.root",
                          GenieAnalysisOriginalCuts::Target::C12,
                          GenieAnalysisOriginalCuts::BeamEnergy::MeV_2261,
                          {"nocut", "π+", "π-"},
                          {"W", "wght", "el_phi", "el_cos_theta", "el_p", "el_E", "el_acceptance", "pi_phi",
                           "pi_cos_theta", "pi_p", "pi_E", "pi_acceptance"},
                          {"ALL", "QE", "RES_ALL", "DELTA1232", "DIS"}};
    /* "PIP", "QE_PIP", "RES_ALL_PIP", "DELTA1232_PIP", "DIS_PIP", */
    /* "PIM", "QE_PIM", "RES_ALL_PIM", "DELTA1232_PIM", "DIS_PIM"}, true, true, false}; */

    /* GenieAnalysisOriginalCuts ga{gst_path_jan, */
    /*                              "output_jan.root", */
    /*                              GenieAnalysisOriginalCuts::Target::C12, */
    /*                              GenieAnalysisOriginalCuts::BeamEnergy::MeV_2261, */
    /*                              {"nocut", "p_gdoc", "p_efid"}, */
    /*                              {"W", "el_phi", "el_cos_theta", "el_p", "el_E", "el_acceptance"}, */
    /*                              {"ALL", "QE", "RES_ALL", "DELTA1232", "DIS"}, false, true, false}; */

    /* GenieAnalysis1Pion ga{genie_machine_path, */
    /*                       "output_jan_full.root", */
    /*                       GenieAnalysisOriginalCuts::Target::C12, */
    /*                       GenieAnalysisOriginalCuts::BeamEnergy::MeV_2261, */
    /*                       {"PIP", "PIM"}, */
    /*                       {"W", "wght", "el_phi", "el_cos_theta", "el_p", "el_E", "el_acceptance", "pi_phi", */
    /*                        "pi_cos_theta", "pi_p", "pi_E", "pi_acceptance"}, */
    /*                       {"ALL", "QE", "RES_ALL", "DELTA1232", "DIS"}}; */

    ga.runAnalysis();

    return 0;
}
