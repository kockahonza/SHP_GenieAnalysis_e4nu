#include "GACLAS6Data.h"
#include <TMath.h>

Double_t GACLAS6DataPrepare::passesCuts() {
    m_el_V3.SetXYZ(m_ge.pxl, m_ge.pyl, m_ge.pzl);
    m_el_V4.SetPxPyPzE(m_el_V3.X(), m_el_V3.Y(), m_el_V3.Z(), sqrt(m_el_V3.Mag2() + mass_electron * mass_electron));

    // Calculation of kinematic quantities (nu, Q2, x bjorken, q and W) -- literally taken from original though
    const TLorentzVector el_change{m_el_V4 - TLorentzVector{0, 0, m_beam_energy_val, m_beam_energy_val}};
    m_reconstructed_Q2 = -el_change.Mag2();

    const double nu = -el_change.E();
    m_bjorken_x = m_reconstructed_Q2 / (2 * mass_proton * nu);

    m_reconstructed_W = TMath::Sqrt((mass_proton + nu) * (mass_proton + nu) - el_change.Vect().Mag2());

    // Hadron loop, first clear the vectors though
    m_pi_minus.clear();
    m_pi_plus.clear();
    m_protons.clear();

    TVector3 V3;
    for (int i{0}; i < m_ge.nf; i++) {
        if (m_ge.pdgf[i] == 2212) { // proton
            V3.SetXYZ(m_ge.pxf[i], m_ge.pyf[i], m_ge.pzf[i]);

            m_protons.push_back({{V3, sqrt(V3.Mag2() + mass_proton * mass_proton)}, V3});
        } else if (m_ge.pdgf[i] == -211) {
            V3.SetXYZ(m_ge.pxf[i], m_ge.pyf[i], m_ge.pzf[i]);

            m_pi_minus.push_back({{V3, sqrt(V3.Mag2() + mass_pion * mass_pion)}, V3});

        } else if (m_ge.pdgf[i] == 211) {
            V3.SetXYZ(m_ge.pxf[i], m_ge.pyf[i], m_ge.pzf[i]);

            m_pi_plus.push_back({{V3, sqrt(V3.Mag2() + mass_pion * mass_pion)}, V3});

        } else if (m_ge.pdgf[i] == 22) {

        } else {
            std::cout << "extra Data particle -- " << m_ge.pdgf[i] << std::endl;
        }
    }

    return m_ge.wght;
}
