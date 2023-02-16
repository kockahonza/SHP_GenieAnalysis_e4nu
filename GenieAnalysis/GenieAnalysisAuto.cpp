#include "GenieAnalysisAuto.h"

bool GenieAnalysisAutoTH1Fs::isType(const string &type, const GenieEvent &ge) {
    if (type == "QE") {
        return ge.qel;
    } else if (type == "RES") {
        return ge.res;
    } else if (type == "DIS") {
        return ge.dis;
    } else if (type == "ALL") {
        return true;
    } else {
        // TODO: improve errors
        throw -1;
    }
}

Double_t GenieAnalysisAutoTH1Fs::getProperty(const string &property, const GenieEvent &ge) {
    if (property == "W") {
        return ge.W;
    } else if (property == "wght") {
        return ge.wght;
    } else if (property == "el_phi") {
        double phi_deg{TVector3(ge.pxl, ge.pyl, ge.pzl).Phi() * TMath::RadToDeg()};
        return (phi_deg < -30) ? phi_deg + 360 : phi_deg;
    } else {
        // TODO: improve errors
        throw -1;
    }
}
