#ifndef GENIE_ANALYSIS_H
#define GENIE_ANALYSIS_H

#include <memory>
#include <string_view>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>


// Very simple struct to have all the fields stored in the genie summary table (gst) format
struct GenieEvent {
    Int_t iev;
    Int_t neu;
    Int_t fspl;
    Int_t tgt;
    Int_t Z;
    Int_t A;
    Int_t hitnuc;
    Int_t hitqrk;
    Int_t resid;
    Bool_t sea;
    Bool_t qel;
    Bool_t mec;
    Bool_t res;
    Bool_t dis;
    Bool_t coh;
    Bool_t dfr;
    Bool_t imd;
    Bool_t imdanh;
    Bool_t singlek;
    Bool_t nuel;
    Bool_t em;
    Bool_t cc;
    Bool_t nc;
    Bool_t charm;
    Int_t neut_code;
    Int_t nuance_code;
    Double_t wght;
    Double_t xs;
    Double_t ys;
    Double_t ts;
    Double_t Q2s;
    Double_t Ws;
    Double_t x;
    Double_t y;
    Double_t t;
    Double_t Q2;
    Double_t W;
    Double_t EvRF;
    Double_t Ev;
    Double_t pxv;
    Double_t pyv;
    Double_t pzv;
    Double_t En;
    Double_t pxn;
    Double_t pyn;
    Double_t pzn;
    Double_t El;
    Double_t pxl;
    Double_t pyl;
    Double_t pzl;
    Double_t pl;
    Double_t cthl;
    Int_t nfp;
    Int_t nfn;
    Int_t nfpip;
    Int_t nfpim;
    Int_t nfpi0;
    Int_t nfkp;
    Int_t nfkm;
    Int_t nfk0;
    Int_t nfem;
    Int_t nfother;
    Int_t nip;
    Int_t nin;
    Int_t nipip;
    Int_t nipim;
    Int_t nipi0;
    Int_t nikp;
    Int_t nikm;
    Int_t nik0;
    Int_t niem;
    Int_t niother;
    Int_t ni;
    Int_t pdgi[2];   //[ni]
    Int_t resc[1];   //[ni]
    Double_t Ei[2];  //[ni]
    Double_t pxi[2]; //[ni]
    Double_t pyi[2]; //[ni]
    Double_t pzi[2]; //[ni]
    Int_t nf;
    Int_t pdgf[120];    //[nf]
    Double_t Ef[120];   //[nf]
    Double_t pxf[120];  //[nf]
    Double_t pyf[120];  //[nf]
    Double_t pzf[120];  //[nf]
    Double_t pf[120];   //[nf]
    Double_t cthf[120]; //[nf]
    Double_t vtxx;
    Double_t vtxy;
    Double_t vtxz;
    Double_t vtxt;
    Double_t sumKEf;
    Double_t calresp0;
};


class GenieAnalysis {
public:
    GenieEvent m_loaded_event;
    const std::unique_ptr<TFile> m_genie_data_file;
    const std::unique_ptr<TTree> m_genie_data;
    

    // Set all the gst field branches in m_genie_data to point at the relevant fields in m_loaded_event
    void point_branches_at_event();

public:
    GenieAnalysis(const char* filename, const char* gst_ttree_name="gst") :
        m_genie_data_file(TFile::Open(filename, "READ")),
        m_genie_data((TTree*)m_genie_data_file->Get(gst_ttree_name))
    {
    };
};

#endif
