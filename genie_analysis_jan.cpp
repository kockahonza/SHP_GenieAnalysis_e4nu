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
    TLorentzVector m_pion_V4;
    TVector3 m_pion_V3;
    Double_t m_pion_acceptance;

  protected:
    map<string, AutoProperty> m_new_known_properties{
        {"pi_phi",
         {"Pion phi",
          {720, -30, 330},
          [this]() {
              double phi_deg{m_pion_V3.Phi() * TMath::RadToDeg()};
              return (phi_deg < -30) ? phi_deg + 360 : phi_deg;
          }}},
        {"pi_cos_theta", {"Pion cos theta", {720, -1, 1}, [this]() { return m_pion_V3.CosTheta(); }}},
        {"pi_mag", {"Pion momentum", {720, 0, 3}, [this]() { return m_pion_V3.Mag(); }}},
        {"pi_acceptance", {"Pion acceptance weight", {100, 0, 1}, [this]() { return m_pion_acceptance; }}}};

  public:
    GenieAnalysis1Pion(const char *filename, const char *output_filename, const vector<string> &stages,
                       const vector<string> &properties, const vector<string> &types,
                       const char *gst_ttree_name = "gst")
        : GenieAnalysisOriginalCuts(filename, output_filename, stages, properties, types) {
        m_known_properties.insert(m_new_known_properties.begin(), m_new_known_properties.end());
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
                std::tie(m_pion_V4, m_pion_V3, m_pion_acceptance) = m_passed_pi_minus[0];
            } else if (num_pi_plus == 1) {
                m_pion_charge = +1;
                std::tie(m_pion_V4, m_pion_V3, m_pion_acceptance) = m_passed_pi_plus[0];
            }

            if (m_pion_acceptance != TMath::Abs(m_pion_acceptance)) {
                // Seems the pion maps have a lot of invalid values so don't use the events that fail, but don't throw
                // and error, it is expected
                return 0;
            }

            /* std::cout << weight << ", " << m_pion_acceptance << ", " << m_pion_charge << std::endl; */
            return weight * m_pion_acceptance;
        } else {
            return weight;
        }
    }
};

int main(int argc, char *argv[]) {
    GenieAnalysis1Pion ga{gst_path_jan,
                          "output_jan.root",
                          {"nocut"},
                          {"W", "el_smeared_mag", "el_smeared_phi", "el_smeared_cos_theta", "el_acceptance", "pi_mag",
                           "pi_phi", "pi_cos_theta", "pi_acceptance"},
                          {"ALL", "QE", "RES_ALL", "DELTA1232", "DIS"}};

    ga.runAnalysis();

    return 0;
}
