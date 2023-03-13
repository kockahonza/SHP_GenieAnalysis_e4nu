#include <iostream>
#include <string>

#include <TCanvas.h>
#include <TFile.h>
#include <THStack.h>
#include <TLegend.h>

// General setup
using std::string, std::cout, std::endl;

const char *default_hstack_opt{"nostack,hist"};
const Int_t canvas_ww{1000};
const Int_t canvas_wh{618};

TCanvas cc{"canvas_i", "Interactive canvas"};

// Some basic plotting used throughout
void draw_stack(TFile *tf, string name) {
    tf->Get<THStack>(name.c_str())->Draw(default_hstack_opt);
    tf->Get<TLegend>((name + "_legend").c_str())->Draw();
}

void draw_rebinned(TFile *tf, string name, Int_t rebin = 2) {
    THStack *hs{tf->Get<THStack>(name.c_str())};
    THStack *rebinhs{new THStack()};
    TH1F *thist;
    for (auto h : *hs->GetHists()) {
        thist = (TH1F *)h;
        thist->Rebin(rebin);
        rebinhs->Add(thist);
    }
    rebinhs->Draw(default_hstack_opt);
    tf->Get<TLegend>((name + "_legend").c_str())->Draw();
}

// Specific data parts
//
// 1Pion data from 0313 either with or without using pion acceptance
/* auto
 * napief{TFile::Open("/home/honza/Sync/University/CurrentCourses/SHP/data/230310/full_piacceptance_test/nopiacceptance/output_full_pie.root",
 * "READ")}; */
/* auto
 * napipf{TFile::Open("/home/honza/Sync/University/CurrentCourses/SHP/data/230310/full_piacceptance_test/nopiacceptance/output_full_pip.root",
 * "READ")}; */
/* auto
 * napimf{TFile::Open("/home/honza/Sync/University/CurrentCourses/SHP/data/230310/full_piacceptance_test/nopiacceptance/output_full_pim.root",
 * "READ")}; */
auto pipf{TFile::Open("output_local_pip.root", "READ")};
auto pimf{TFile::Open("output_local_pim.root", "READ")};

void plots() { /* gStyle->SetOptStat(""); */ }
