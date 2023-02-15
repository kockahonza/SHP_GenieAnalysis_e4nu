#ifndef MISC_H
#define MISC_H

#include "GenieAnalysis.h"

// Not sure if this is a useful thing or not, but keeping it for now, it's the same values as in the original code
enum class InteractionType {
    AllInteractions = 0,
    QuasiElastic = 1,
    MesonExchangeCurrent = 2,
    Resonance = 3,
    DeepInelastic = 4,
    CoherentMesonProduction = 5
};

constexpr double m_pimi = 0.139570, m_pipl = 0.139570, m_pion = 0.139570;
constexpr double m_prot = 0.9382720813, m_neut = 0.939565;
constexpr double H3_bind_en = 0.008481, He4_bind_en = 0.0283, C12_bind_en = 0.09215, B_bind_en = 0.0762;
constexpr double He3_bind_en = 0.0077, D2_bind_en = 0.00222, Fe_bind_en = 0.49226, Mn_bind_en = 0.4820764;
constexpr double e_mass = 0.000510998;

constexpr double fine_struc_constexpr = 0.007297;
constexpr double ns_to_s = 1.0E-9;
constexpr Double_t c = 2.99792E+10;

// smithja: Note these variables are only half of the Delta-phi cuts made in the
// genie_analysis code, e.g. PhiOpeningAngleEl = 6 --> Delta-phi for the electron is 12.
constexpr double PhiOpeningAngleEl = 6;
constexpr double PhiOpeningAngleProt = 22.5;
constexpr double PhiOpeningAngle = 6.;

static constexpr Float_t par_EcUVW[6][3] = {{60, 360, 400}, {55, 360, 400}, {50, 363, 400},
                                            {52, 365, 396}, {60, 360, 398}, {50, 362, 398}};

constexpr double MinThetaProton = 12.;
constexpr double MinThetaPiPlus = 12.;
constexpr double MinThetaPiMinus = 0.;
constexpr double MinThetaGamma = 8.;

constexpr double CenterFirstSector = 30;
constexpr double CenterSecondSector = 90;
constexpr double CenterThirdSector = 150;
constexpr double CenterFourthSector = 210;
constexpr double CenterFifthSector = 270;
constexpr double CenterSixthSector = 330;

constexpr int RotCounterLimit = 100;

#endif
