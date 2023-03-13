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

/* TCanvas cc{"canvas_i", "Interactive canvas"}; */

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
auto fstpip{TFile::Open("output_local_FSt_pip.root", "READ")};
auto fstpim{TFile::Open("output_local_FSt_pim.root", "READ")};
auto fstex{TFile::Open("output_local_FSt_e.root", "READ")};
auto fstex_noa{TFile::Open("output_local_FSt_e_noa.root", "READ")};

void fs_ps_p_n(TFile *dataf) {
    TCanvas *tc{new TCanvas("fs_ps_p_n", "ps,fs vs p,n", canvas_ww, canvas_wh)};
    tc->Divide(2, 2);

    THStack *hs;

    tc->cd(1);
    hs = dataf->Get<THStack>("ps_num_protons");
    hs->Draw("nostack");
    dataf->Get<TLegend>("ps_num_protons_legend")->Draw();

    tc->cd(2);
    hs = dataf->Get<THStack>("ps_num_neutrons");
    hs->Draw("nostack");
    dataf->Get<TLegend>("ps_num_neutrons_legend")->Draw();

    tc->cd(3);
    hs = dataf->Get<THStack>("fs_num_protons");
    hs->Draw("nostack");
    dataf->Get<TLegend>("fs_num_protons_legend")->Draw();

    tc->cd(4);
    hs = dataf->Get<THStack>("fs_num_neutrons");
    hs->Draw("nostack");
    dataf->Get<TLegend>("fs_num_neutrons_legend")->Draw();

    /* tc->Print("draft_plots/pions_ps_to_fs.pdf"); */
}

void pions_ps_to_fs(TFile *dataf) {
    TCanvas *tc{new TCanvas("pion_ps_to_fs", "pion_ps_to_fs", canvas_ww, canvas_wh)};
    tc->Divide(2, 2);

    THStack *hs;

    tc->cd(1);
    hs = dataf->Get<THStack>("ps_num_pip");
    hs->Draw("nostack");
    dataf->Get<TLegend>("ps_num_pip_legend")->Draw();

    tc->cd(2);
    hs = dataf->Get<THStack>("ps_num_pim");
    hs->Draw("nostack");
    dataf->Get<TLegend>("ps_num_pim_legend")->Draw();

    tc->cd(3);
    hs = dataf->Get<THStack>("fs_num_pip");
    hs->Draw("nostack");
    dataf->Get<TLegend>("fs_num_pip_legend")->Draw();

    tc->cd(4);
    hs = dataf->Get<THStack>("fs_num_pim");
    hs->Draw("nostack");
    dataf->Get<TLegend>("fs_num_pim_legend")->Draw();

    /* tc->Print("../plotting/draft_plots/pions_ps_to_fs.pdf"); */
}


void plots() { /* gStyle->SetOptStat(""); */
}
