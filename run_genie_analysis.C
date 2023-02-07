#ifndef rungenieanalysis
#define rungenieanalysis
#include "genie_analysis.h"

#include <iostream>

using namespace std;

int main(int argc, char **argv) {
  std::cout << std::endl; // separate this initial output from the last command typed in the terminal
  std::cout << "/--------------------------------------------------------\\" << std::endl;
  std::cout << "| Precision of kinematic cuts AS DOCUMENTED IN THE FILE  |" << std::endl;
  std::cout << "| NAMES truncates at two decimal places, but the code    |" << std::endl;
  std::cout << "| still makes the complete cut, e.g. 1.657 cut is used   |" << std::endl;
  std::cout << "| in its entirety, but the file name labels it as 1pt65. |" << std::endl;
  std::cout << "| You can change that precision in genie_analysis.C (see |" << std::endl;
  std::cout << "| 'CAUTION' statements) or you can manually rename the   |" << std::endl;
  std::cout << "| file after processing. If editing the precision in     |" << std::endl;
  std::cout << "| genie_analysis.C, be sure to update this statement in  |" << std::endl;
  std::cout << "| run_genie_analysis.C. - J. Smith (smithja) on 10/27/21 |" << std::endl;
  std::cout << "\\--------------------------------------------------------/" << std::endl;
  std::cout << std::endl; // separate this output from the rest of the output

  // TODO: This message doesn't mention the proton number and acceptance parameters and it should
  if (argc < 15) {
    std::cout << "/-------------------------------------------------------\\" << std::endl;
    std::cout << "| Please specify the target (3He, 56Fe, C12, 4He), the  |" << std::endl;
    std::cout << "| beam energy (1161, 2261 or 4461), the data type       |" << std::endl;
    std::cout << "| (CLAS=0 or SuSav2 simulation=1 or G18_10a_02_11a      |" << std::endl;
    std::cout << "| simulation=2) and the number of rotations.            |" << std::endl;
    std::cout << "|                                                       |" << std::endl;
    std::cout << "| Note not all target and beam energy combinations are  |" << std::endl;
    std::cout << "| valid.                                                |" << std::endl;
    std::cout << "|    	       	       	       	       	       	          |" << std::endl;
    std::cout << "| Please also specify the kinematic cuts you want to    |" << std::endl;
    std::cout << "| make. The order of those cuts is as follows:          |" << std::endl;
    std::cout << "|  1.  electron sectors (sectors 1 through 6 = 0,       |" << std::endl;
    std::cout << "|      sectors 1, 2, and 6 = 1, sectors 1 and 6 = 2)    |" << std::endl;
    std::cout << "|  2.  flag to make cut on Delta-phi for the electron   |" << std::endl;
    std::cout << "|      (false = 0 or true = 1)                          |" << std::endl;
    std::cout << "|  3.  electron theta lower bound (in degrees, 0 to 180)|" << std::endl;
    std::cout << "|  4.  electron theta upper bound (in degrees, 0 to 180)|" << std::endl;
    std::cout << "|  5.  electron momentum lower bound (in GeV, greater   |" << std::endl;
    std::cout << "|      than or equal to 0)                              |" << std::endl;
    std::cout << "|  6.  flag to require proton in sector m to have a     |" << std::endl;
    std::cout << "|      corresponding electron in sector m-3             |" << std::endl;
    std::cout << "|  7.  flag to make cut on Delta-phi for the proton     |" << std::endl;
    std::cout << "|      (false = 0, true = 1)                            |" << std::endl;
    std::cout << "|  8.  proton theta lower bound (in degrees, 0 to 180)  |" << std::endl;
    std::cout << "|  9.  proton theta upper bound (in degrees, 0 to 180)  |" << std::endl;
    std::cout << "|  10. proton momentum lower bound (in GeV, greater than|" << std::endl;
    std::cout << "|      or equal to 0)                                   |" << std::endl;
    std::cout << "|  11. proton momentum upper bound (in GeV, greater than|" << std::endl;
    std::cout << "|      or equal to 0)                                   |" << std::endl;
    std::cout << "\\-------------------------------------------------------/" << std::endl;
    std::cout << std::endl; // separate output
    std::cout << "====== Usage ======" << std::endl;
    std::cout << "./genie_analysis target beam_energy 0/1 #rot 0/1/2 0/1 thetaEl_lb thetaEl_ub elMom_lb 0/1/2 0/1 thetaProt_lb "
                 "thetaProt_ub protMom_lb protMom_ub"
              << std::endl;
    exit(1);
  }

  std::string target = argv[1];
  std::string beam_en = argv[2];
  int choice = atoi(argv[3]);
  int rotations = atoi(argv[4]);
  int elSectors_flag = atoi(argv[5]);
  int deltaPhiEl = atoi(argv[6]);
  int thetaEl_lb = stod(argv[7]);
  int thetaEl_ub = stod(argv[8]);
  double elMom_lb = stod(argv[9]);
  int protSectors_flag = atoi(argv[10]);
  int deltaPhiProt = atoi(argv[11]);
  int thetaProt_lb = stod(argv[12]);
  int thetaProt_ub = stod(argv[13]);
  double protMom_lb = stod(argv[14]);
  double protMom_ub = stod(argv[15]);

  // vchand - added extra parameters which are missing, corresponds to detector acceptance and number of protons
  int detector_acceptance = atoi(argv[16]);
  int NumOfProton = atoi(argv[17]);

  if (choice != 2 && choice != 1 && choice != 0) {
    std::cout << "Unknown option for data type (parameter 3). It should be either 0 or 1 (or 2). The given value is " << choice
              << std::endl;
    return 0;
  }

  if (rotations <= 0) {
    std::cout << "Not a valid number for the number of rotations (parameter 4). It should be greater than 0. The given value is "
              << rotations << std::endl;
    return 0;
  }

  if (elSectors_flag != 2 && elSectors_flag != 1 && elSectors_flag != 0) {
    std::cout << "Unknown flag for electron sectors (parameter 5). It should be either 0, 1, or 2. The given value is " << elSectors_flag
              << std::endl;
    std::cout << "If you are wanting to add a new option for this parameter, make sure you propogate all changes in genie_analysis.h, "
                 "genie_analysis.C, and run_genie_analysis.C."
              << std::endl;
    return 0;
  }

  if (deltaPhiEl != 1 && deltaPhiEl != 0) {
    std::cout << "Not a valid number for Delta-phi for the electron (parameter 6). It should be either 0 or 1. The given value is "
              << deltaPhiEl << std::endl;
    return 0;
  }

  if (thetaEl_lb < 0 || thetaEl_lb > 180) {
    std::cout << "Not a valid number for the electron theta lower bound (parameter 7). It should be between 0 and 180 (inclusive). The "
                 "given value is "
              << thetaEl_lb << std::endl;
    return 0;
  }

  if (thetaEl_ub < 0 || thetaEl_ub > 180) {
    std::cout << "Not a valid number for the electron theta upper bound (parameter 8). It should be between 0 and 180 (inclusive). The "
                 "given value is "
              << thetaEl_ub << std::endl;
    return 0;
  }

  if (thetaEl_lb > thetaEl_ub) {
    std::cout << "The electron theta lower bound (parameter 7) must be less than or equal to the electron theta upper bound (parameter 8)."
              << std::endl;
    return 0;
  }

  if (elMom_lb < 0) {
    std::cout << "Not a valid number for the electron momentum lower bound (parameter 9). It should be greater than or equal to 0. The "
                 "given value is "
              << elMom_lb << std::endl;
    return 0;
  }

  if (protSectors_flag != 1 && protSectors_flag != 0) {
    std::cout << "Unknown flag for proton sectors (parameter 10). It should be either 0 or 1. The given value is " << protSectors_flag
              << std::endl;
    std::cout << "If you are wanting to add a new option for this parameter, make sure you propogate all changes in genie_analysis.h, "
                 "genie_analysis.C, and run_genie_analysis.C."
              << std::endl;
    return 0;
  }

  if (deltaPhiProt != 1 && deltaPhiProt != 0) {
    std::cout << "Not a valid number for Delta-phi for the proton (parameter 11). It should be either 0 or 1. The given value is "
              << deltaPhiProt << std::endl;
    return 0;
  }

  if (thetaProt_lb < 0 || thetaProt_lb > 180) {
    std::cout << "Not a valid number for the proton theta lower bound (parameter 12). It should be between 0 and 180 (inclusive). The "
                 "given value is "
              << thetaProt_lb << std::endl;
    return 0;
  }

  if (thetaProt_ub < 0 || thetaProt_ub > 180) {
    std::cout << "Not a valid number for the proton theta upper bound (parameter 13). It should be between 0 and 180 (inclusive). The "
                 "given value is "
              << thetaProt_ub << std::endl;
    return 0;
  }

  if (thetaProt_lb > thetaProt_ub) {
    std::cout << "The proton theta lower bound (parameter 12) must be less than or equal to the proton theta upper bound (parameter 13)."
              << std::endl;
    return 0;
  }

  if (protMom_lb < 0) {
    std::cout << "Not a valid number for the proton momentum lower bound (parameter 14). It should be greater than or equal to 0. The "
                 "given value is "
              << protMom_lb << std::endl;
    return 0;
  }
  if (protMom_ub < 0) {
    std::cout << "Not a valid number for the proton momentum upper bound (parameter 15). It should be greater than or equal to 0. The "
                 "given value is "
              << protMom_ub << std::endl;
    return 0;
  }
  if (protMom_lb > protMom_ub) {
    std::cout
        << "The proton momentum lower bound (parameter 14) must be less than or equal to the proton momentum upper bound (parameter 15)."
        << std::endl;
    return 0;
  }

  genie_analysis t(target, beam_en, rotations, choice, elSectors_flag, deltaPhiEl, thetaEl_lb, thetaEl_ub, elMom_lb, protSectors_flag,
                   deltaPhiProt, thetaProt_lb, thetaProt_ub, protMom_lb, protMom_ub, NumOfProton, detector_acceptance);
  t.Loop();

  return 0;
}
#endif
