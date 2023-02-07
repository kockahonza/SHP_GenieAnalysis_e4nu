#include "genie_analysis.h"

#include <iostream>

int main(int argc, char *argv[]) {
    genie_analysis t(target, beam_en, rotations, choice, elSectors_flag, deltaPhiEl, thetaEl_lb, thetaEl_ub, elMom_lb,
                     protSectors_flag, deltaPhiProt, thetaProt_lb, thetaProt_ub, protMom_lb, protMom_ub, NumOfProton,
                     detector_acceptance);

    return 0;
}
