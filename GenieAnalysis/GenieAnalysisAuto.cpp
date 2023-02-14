#include "GenieAnalysisAuto.h"

bool GenieAnalysisAutoTH1Fs::isObservable(const string &observable, const GenieEvent &ge) {
    if (observable == "QE") {
        return ge.qel;
    } else if (observable == "RES") {
        return ge.res;
    } else if (observable == "DIS") {
        return ge.dis;
    } else {
        // TODO: improve errors
        throw -1;
    }
}

Double_t GenieAnalysisAutoTH1Fs::getProperty(const string &property, const GenieEvent &ge) {
    if (property == "W") {
        return ge.W;
    } else {
        // TODO: improve errors
        throw -1;
    }
}
