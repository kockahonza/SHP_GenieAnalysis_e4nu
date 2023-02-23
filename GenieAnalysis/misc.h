#ifndef MISC_H
#define MISC_H

#include <TMath.h>

// Useful environment constants
constexpr const char *gst_path_jan{"/home/honza/Sync/University/CurrentCourses/SHP/data/Genie_gst_2000000.root"};
constexpr const char *genie_machine_path{"/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/"
                                         "C12_2261GeV/apapadop_SuSav2_C12_2261GeV_master.root"};

// Constants used within GenieAnalysis
constexpr double m_pion{0.139570};
constexpr double e_mass{0.000510998};

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
