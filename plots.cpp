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
// Final state transparency studies
auto fstpip{TFile::Open("output_local_FSt_pip.root", "READ")};
auto fstpim{TFile::Open("output_local_FSt_pim.root", "READ")};
auto fstex{TFile::Open("output_local_FSt_e.root", "READ")};
auto fstex_noa{TFile::Open("output_local_FSt_e_noa.root", "READ")};

TCanvas *fs_ps_p_n(TFile *dataf) {
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

    return tc;
}

TCanvas *pions_ps_to_fs(TFile *dataf) {
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

    return tc;
}

// 1Pion 1 nucleon studies, or the beginnings of
auto opop_1pip0pim1p0n{TFile::Open("output_local_1pip0pim1p0n.root", "READ")};
auto opop_0pip1pim1p0n{TFile::Open("output_local_0pip1pim1p0n.root", "READ")};

void plots() {
    pions_ps_to_fs(fstex_noa)->Print("../plotting/draft_plots/fst_pions_ps_to_fs.pdf");
    fs_ps_p_n(fstpip)->Print("../plotting/draft_plots/fst_pip_psfs_pn.pdf");
    fs_ps_p_n(fstpim)->Print("../plotting/draft_plots/fst_pim_psfs_pn.pdf");
}
