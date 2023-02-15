#include "GenieAnalysisAuto.h"

bool GenieAnalysisAutoTH1Fs::isType(const string &type, const GenieEvent &ge) {
    if (type == "QE") {
        return ge.qel;
    } else if (type == "RES") {
        return ge.res;
    } else if (type == "DIS") {
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
