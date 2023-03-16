#ifndef MISC_H
#define MISC_H

#include <TMath.h>

// Useful environment constants

// Constants used within GenieAnalysis
constexpr double mass_pion{0.139570};
constexpr double mass_electron{0.000510998};
constexpr double mass_proton{0.9382720813};

// Parameters for Fiducial I expect, they are all used in Fiducial.cpp
constexpr double PhiOpeningAngleEl = 6;
constexpr double PhiOpeningAngleProt = 22.5;
constexpr double PhiOpeningAngle = 6.;

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

#endif
