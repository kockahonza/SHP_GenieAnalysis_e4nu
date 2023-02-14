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

bool isQuasiElastic(GenieEvent ge) { return ge.qel; }

bool isResonant(GenieEvent ge) { return ge.res; }

bool isDeepInelastic(GenieEvent ge) { return ge.dis; }

#endif
