// SEARCH THROUGH THIS FILE TO FIND ALL THE 'CODE WARNING' AND 'CODE CAUTION' STATEMENTS BEFORE RUNNING THIS CODE.
// LATALING: NOTE THAT PROTONS SHOULD HAVE pf[i] > 0.3 GeV AND PIONS pf[i] > 0.15GeV!!!!
// You can remove them but do so at your own peril (the "peril" being that Fe56 gives a segmentation violation)
// ctrl+f and look for comments that say "here" which will tell you where they are (there are only 3 places)

// Authors: Afroditi Papadopoulou (apapadop), Graham Chambers-Wall (gchamber), Jacob Smith (smithja)
// Date	of Creation: N/A (you could ask Afroditi about this)

#define GENIE_ANALYSIS_C

#include "genie_analysis.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH3.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TVectorT.h>

#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

using namespace std;

// -----------------------------------------------------------------------------------------------------------------------------------------------

// Not used anywhere
vector<double> CalculateCalKineVars(double ECal, TLorentzVector FSElectron) {

    vector<double> CalKineVars;
    CalKineVars.clear();

    TLorentzVector V4_beam_Cal(0, 0, ECal, ECal);
    double nu_Cal = -(FSElectron - V4_beam_Cal).E();
    double Q2_Cal = -(FSElectron - V4_beam_Cal).Mag2();
    double x_bjk_Cal = Q2_Cal / (2 * m_prot * nu_Cal);
    TVector3 V3_q_Cal = (V4_beam_Cal - FSElectron).Vect();
    double W_var_Cal = TMath::Sqrt((m_prot + nu_Cal) * (m_prot + nu_Cal) - V3_q_Cal * V3_q_Cal);

    CalKineVars.push_back(nu_Cal);    // 0-th element: energy transfer using Ecal
    CalKineVars.push_back(Q2_Cal);    // 1st element: Q2 using Ecal
    CalKineVars.push_back(x_bjk_Cal); // 2nd element: xB using Ecal
    CalKineVars.push_back(W_var_Cal); // 3rd element: invariant mass using Ecal

    return CalKineVars;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Loading all the constants from Constant.h ( e_mass, m_prot, m_pimi, m_pipl, m_pion, m_neut = 0.939565,
// H3_bind_en, He4_bind_en, C12_bind_en, B_bind_en, He3_bind_en, D2_bind_en, Fe_bind_en, Mn_bind_en)

void genie_analysis::Loop() {

    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int NSectors = 6;

    const int NInt = 6; // All Interactions = 0, QE = 1, MEC = 2, RES = 3, DIS = 4, COH = 5
    const int NumInteractions =
        7; // Lataling: All Interactions = 0, QE = 1, MEC = 2, OtherRES = 3, DIS = 4, 1232RES = 5, COH = 6

    // ---------------------------------------------------------------------------------------------------------------

    // Apply Cuts

    // Do we want to apply fiducials & the acceptance map weights
    // Do we want a truth level study ? if so, stop ditching sectors

    bool UseAllSectors = false;
    bool ApplyFiducials = true;
    bool ApplyAccWeights = true;
    bool ApplyReso = true;

    // Lataling: turn these off to remove acceptances
    if (detector_acceptance == 0) {
        UseAllSectors = true;
        ApplyFiducials = false;
        ApplyAccWeights = false;
        ApplyReso = false; // Make this true if you want to make cuts on signal -- since this changes the kinematics of
                           // the event
    }

    // Lataling: use this to turn the threshold momentum restrictions on or off
    bool Applymomthresh = true;

    bool TruthLevel1p0piSignalStudy = false;
    bool TruthLevel0piSignalStudy = false;

    bool UsePhiThetaBand = false;

    double PtMax = 0.2; // gchamber: max pt cut (1e1p spectrum)

    // ---------------------------------------------------------------------------------------------------------------

    // fchoice = 0 is for analysis of CLAS data while fchoice = 1 is for the analysis of GENIE simulation
    if (fchoice != 4 && fchoice != 3 && fchoice != 2 && fchoice != 1 &&
        fchoice != 0) { // if the user has specified an invalid fchoice, let them know
        std::cout << "This parameter value is not implemented in genie_analysis::Loop(). It should be either 0 or 1. "
                     "The given value is "
                  << fchoice << std::endl;
        std::exit(0);
    }

    std::map<std::string, double> bind_en;
    std::map<std::string, double> target_mass;
    std::map<std::string, double> residual_target_mass;
    std::map<std::string, double> Ecal_offset; // that might not be necessary for simulation data

    target_name = ftarget; // std string for target name
    en_beam["1161"] = 1.161;
    en_beam["2261"] = 2.261;
    en_beam["4461"] = 4.461;

    en_beam_Ecal["1161"] = 1.161;
    en_beam_Ecal["2261"] = 2.261;
    en_beam_Ecal["4461"] = 4.461;

    en_beam_Eqe["1161"] = 1.161;
    en_beam_Eqe["2261"] = 2.261;
    en_beam_Eqe["4461"] = 4.461;

    if (fChain == 0)
        return;

    Long64_t nentries = fChain->GetEntriesFast();
    //	nentries = 1000;
    nentries = 10000000; // smithja: max for C-12 simulation to not crash
                         //	nentries = 200000000;

    // smithja: Tells the user if nentries is maually set and if the
    //          manually set value is valid given the amount of
    //          actual entries in a given data file.
    if (nentries < fChain->GetEntriesFast()) {
        std::cout << "nentries IS HARD CODED AND LESS THAN THE ACTUAL NUMBER OF EVENTS IN THIS FILE!!!" << std::endl
                  << std::endl;
    }
    if (nentries > fChain->GetEntriesFast()) {
        std::cout << "nentries is greater than number of actual entries in file." << std::endl;
        std::cout << "To not cause any weird errors in the loops, the program" << std::endl;
        std::cout << "reset nentries to be equal to the number of actual" << std::endl;
        std::cout << "entries in this file." << std::endl;
        nentries = fChain->GetEntriesFast();
    }

    // Resolutions for Smearing for GENIE simulation data
    double reso_p = 0.01;   // smearing for the proton
    double reso_e = 0.005;  // smearing for the electrons
    double reso_pi = 0.007; // smearing for pions, executive decision by Larry (28.08.19)

    if (!ApplyReso) {
        reso_p = 0.;
        reso_e = 0.;
        reso_pi = 0;
    }

    // Resolution defined above seems to be insufficient at 1.1 GeV -> tripled it for all particles
    if (fbeam_en == "1161") {
        reso_p = 3 * reso_p;
        reso_e = 3 * reso_e;
        reso_pi = 3 * reso_pi;
    }

    // Lataling: I do not know why these cuts are here, so I have commented them out for now
    // double Wcut = 2; // cut for all beam energies < 2
    // double Q2cut = 0; // cut for 1.1 GeV > 0.1, for 2.2 GeV > 0.4 and 4.4 GeV > 0.8

    const int n_slice = 3; // Stick to the 3 slices
    const double pperp_min[n_slice] = {0., 0.2, 0.4};
    const double pperp_max[n_slice] = {0.2, 0.4, 10.};

    TVector3 V3_rotprot1, V3_rotprot2, V3_rotprot3, V3_rot_pi, V3_rotprot;

    TString E_acc_file;

    // Lataling: The Q2 cuts have been removed here as well
    if (en_beam[fbeam_en] > 1. && en_beam[fbeam_en] < 2.) { // 1.1 GeV Configuration parameters and cuts
        E_acc_file = "1_161";
        // Q2cut = 0.1;
    }

    if (en_beam[fbeam_en] > 2. && en_beam[fbeam_en] < 3.) { // 2.2 GeV Configuration parameters and cuts
        E_acc_file = "2_261";
        // Q2cut = 0.4;
    }

    if (en_beam[fbeam_en] > 4. && en_beam[fbeam_en] < 5.) { // 4.4 GeV Configuration parameters and cuts
        E_acc_file = "4_461";
        // Q2cut = 0.8;
    }

    // Further constants for binding energies and target masses

    Ecal_offset["3He"] = 0.004;
    Ecal_offset["4He"] = 0.005;
    Ecal_offset["C12"] = 0.005;
    Ecal_offset["56Fe"] = 0.011;

    bind_en["3He"] = He3_bind_en - D2_bind_en + Ecal_offset["3He"]; // the offset is used to shift the peak to be at 0
    bind_en["4He"] = He4_bind_en - H3_bind_en + Ecal_offset["4He"];
    bind_en["C12"] = C12_bind_en - B_bind_en + Ecal_offset["C12"];
    bind_en["56Fe"] = Fe_bind_en - Mn_bind_en + Ecal_offset["56Fe"];
    bind_en["CH2"] = C12_bind_en - B_bind_en;

    target_mass["3He"] = 2 * m_prot + m_neut - He3_bind_en;
    target_mass["4He"] = 2 * m_prot + 2 * m_neut - He4_bind_en;
    target_mass["C12"] = 6 * m_prot + 6 * m_neut - C12_bind_en;
    target_mass["56Fe"] = 26 * m_prot + 30 * m_neut - Fe_bind_en;
    target_mass["CH2"] = 6 * m_prot + 6 * m_neut - C12_bind_en;

    residual_target_mass["3He"] = m_prot + m_neut - D2_bind_en;
    residual_target_mass["4He"] = m_prot + 2 * m_neut - H3_bind_en;
    residual_target_mass["C12"] = 5 * m_prot + 6 * m_neut - B_bind_en;
    residual_target_mass["56Fe"] = 25 * m_prot + 30 * m_neut - Mn_bind_en;
    residual_target_mass["CH2"] = 25 * m_prot + 30 * m_neut - Mn_bind_en;

    // ----------------------------------------------------------------------------------------
    // NOAH :: Remove everything with E_tot in it (calorimetric reconstruction, do not need.
    // Care about data driven background subtraction (funfction of p_p or p_e))

    TH2F *h2_N_pi_phot[20];

    gRandom = new TRandom3();
    gRandom->SetSeed(10);

    TLorentzVector V4_beam(0, 0, en_beam[fbeam_en], en_beam[fbeam_en]);
    TLorentzVector V4_target(0, 0, 0, target_mass[ftarget]);

    // Acceptance Maps

    TString WhichMap = "e2a_maps";
    TFile *file_acceptance;
    TFile *file_acceptance_p;
    TFile *file_acceptance_pip;
    TFile *file_acceptance_pim;

    TString Target = "12C";
    if (ftarget.c_str() == "3He") {
        Target = "3He";
    }
    if (ftarget.c_str() == "4He") {
        Target = "4He";
    }

    if (fchoice > 0) { // Only need acceptance maps for GENIE simulation
        file_acceptance = TFile::Open(WhichMap + "/" + WhichMap + "_" + Target + "_E_" + E_acc_file + ".root");
        file_acceptance_p = TFile::Open(WhichMap + "/" + WhichMap + "_" + Target + "_E_" + E_acc_file + "_p.root");
        file_acceptance_pip = TFile::Open(WhichMap + "/" + WhichMap + "_" + Target + "_E_" + E_acc_file + "_pip.root");
        file_acceptance_pim = TFile::Open(WhichMap + "/" + WhichMap + "_" + Target + "_E_" + E_acc_file + "_pim.root");
    }

    // ---------------------------------------------------------------------------------------------------------------

    double XSecScale = 1.;
    //	TFile* XSecFile = TFile::Open("/uboone/app/users/apapadop/R-3_0_6/mySplines/xsec_gxspl-FNALbig.root");

    //	TGraph* gr = NULL;

    //	if (XSecFile) {
    //		TDirectory* dir = (TDirectory*)(XSecFile->Get("nu_mu_C12"));
    //		gr = (TGraph*)dir->Get("tot_cc");
    //	}

    // ---------------------------------------------------------------------------------------------------------------

    // Output file definition

    // CODE WARNING: I, Jacob Smith (smithja), added the tag 'MottXSecEq1' to all of the GENIE file names
    //          to go along with the analysis I have been doing as part of the 2021 Fall SULI program.
    //          The appropriate variables are set throughout this code to document that change. Future
    //          users may want to remove this tag. To change the MottXSecEq1 tag, look at the FileName
    //          statements below.
    // CAUTION: All kinematic cuts in the file names are given to a set precision. However, all the
    //          cuts made in this analysis use the full cut you specified when doing ./genie_analysis.
    //          To change the precison of how the cuts are labeled in the output file name, change the
    //          second argument of the .substr() function in lines relating to 'deci_num' variables.

    TFile *file_out;
    TString FileName = "";

    t_Run->SetTitle((std::string(fchoice > 0 ? "Genie" : "Data") + "_" + (fchoice > 0 ? std::to_string(fchoice) : "") +
                     "_" + std::string(detector_acceptance == 0 ? "noAcc" : "Acc") + "_" + std::to_string(NumOfProton) +
                     "P")
                        .c_str());
    t_target->SetTitle(ftarget.c_str());
    t_beam_en->SetVal(en_beam[fbeam_en]);
    // Save these options
    my_options->Add(t_Run);
    my_options->Add(t_target);
    my_options->Add(t_beam_en);

    // Number of Entries
    t_nentries->SetVal(nentries);
    std::cout << ""
              << "\n";
    FileName = ("./tmp/" + std::string(t_target->GetTitle()) + "/" + std::string(t_Run->GetTitle()) + "_" +
                std::string(t_target->GetTitle()) + "_" + std::to_string(t_beam_en->GetVal()) + ".root")
                   .c_str();
    file_out = new TFile(FileName, "Recreate");

    // Write out TList of run options see
    // (https://root.cern/doc/master/classTCollection.html#a3992401270fb2383d6d6003b60d81146)
    my_options->Write("Run_Info", TObject::kSingleKey);

    std::cout << "el theta lb: " << t_thetaEl_lb->GetVal() << ", el theta ub: " << t_thetaEl_ub->GetVal() << "\n";
    std::cout << "prot theta lb: " << t_thetaProt_lb->GetVal() << ", prot theta ub: " << t_thetaProt_ub->GetVal()
              << "\n";
    std::cout << "Apply theta slice prot = " << fApplyThetaSliceProt << "\n";
    std::cout << "prot mom lb: " << t_ProtMom_lb->GetVal() << ",prot mom ub: " << t_ProtMom_ub->GetVal() << "\n";
    std::cout << "Apply mom slice prot = " << fApplyProtMomCut << "\n";
    std::cout << "Apply phi opening angle cut on electrons: " << fApplyPhiOpeningAngleEl << "\n";
    std::cout << "Apply phi opening angle cut on protons: " << fApplyPhiOpeningAngleProt << "\n";
    // ---------------------------------------------------------------------------------------------------------------
    // NOAH ::
    fiducialcut->InitPiMinusFit(fbeam_en);

    // initialize Fiducial functions for EC limits
    fiducialcut->InitEClimits();
    std::cout << " Test InitEClimits Loop " << fiducialcut->up_lim1_ec->Eval(60) << std::endl;

    // Lataling: Initialising the histograms & PiType variable
    int PiType = 3; // Values 0,1,2 where 0=pi-,1=pi0,2=pi+
                    // TH2F *h2_pion_mom_theta[NumInteractions][PiType];
    // TH2F *h2_pion_mom_phi[NumInteractions][PiType];
    TH1F *h1_pion_momentum[NumInteractions][PiType];
    TH1F *h1_pion_theta[NumInteractions][PiType];
    TH1F *h1_pion_phi[NumInteractions][PiType];
    TH1F *h1_pion_omega[NumInteractions][PiType];
    TH1F *h1_pion_electrontheta[NumInteractions][PiType];
    TH1F *h1_pion_electronmomentum[NumInteractions][PiType];
    TH1F *h1_pion_w[NumInteractions][PiType];
    TH1F *h1_pion_q2[NumInteractions][PiType];
    TH1F *h1_pion_x[NumInteractions][PiType];
    TH1F *h1_pion_p_npi[NumInteractions][PiType];
    TH1F *h1_pion_p_1pi[NumInteractions][PiType];
    TH2I *h2_threshold_passing[NumInteractions];

    // Looping over both NumInteractions and PiType
    for (int WhichInt = 0; WhichInt < NumInteractions; WhichInt++) {
        h2_threshold_passing[WhichInt] =
            new TH2I("h2_threshold_passing_" + TString(std::to_string(WhichInt)), "", 3, 1, 4, 3, 1, 4);
        for (int WhichType = 0; WhichType < PiType; WhichType++) {
            // h2_pion_mom_theta[WhichInt][WhichType] = new
            // TH2F("h2_pion_mom_theta_interaction_"+TString(std::to_string(WhichInt))+"_type_"+TString(std::to_string(WhichType)),"",500,0,5,360,0,360);
            // h2_pion_mom_phi[WhichInt][WhichType] = new
            // TH2F("h2_pion_mom_phi_interaction_"+TString(std::to_string(WhichInt))+"_type_"+TString(std::to_string(WhichType)),"",500,0,5,360,0,360);
            h1_pion_momentum[WhichInt][WhichType] =
                new TH1F("h1_pion_momentum_interaction_" + TString(std::to_string(WhichInt)) + "_type_" +
                             TString(std::to_string(WhichType)),
                         "", 500, 0, 5);
            h1_pion_theta[WhichInt][WhichType] =
                new TH1F("h1_pion_theta_interaction_" + TString(std::to_string(WhichInt)) + "_type_" +
                             TString(std::to_string(WhichType)),
                         "", 360, 0, 360);
            h1_pion_phi[WhichInt][WhichType] = new TH1F("h1_pion_phi_interaction_" + TString(std::to_string(WhichInt)) +
                                                            "_type_" + TString(std::to_string(WhichType)),
                                                        "", 360, 0, 360);
            h1_pion_omega[WhichInt][WhichType] =
                new TH1F("h1_pion_omega_interaction_" + TString(std::to_string(WhichInt)) + "_type_" +
                             TString(std::to_string(WhichType)),
                         "", 500, 0, 5);
            h1_pion_electrontheta[WhichInt][WhichType] =
                new TH1F("h1_pion_electrontheta_interaction_" + TString(std::to_string(WhichInt)) + "_type_" +
                             TString(std::to_string(WhichType)),
                         "", 360, 0, 360);
            h1_pion_electronmomentum[WhichInt][WhichType] =
                new TH1F("h1_pion_electronmomentum_interaction_" + TString(std::to_string(WhichInt)) + "_type_" +
                             TString(std::to_string(WhichType)),
                         "", 500, 0, 5);
            h1_pion_w[WhichInt][WhichType] = new TH1F("h1_pion_w_interaction_" + TString(std::to_string(WhichInt)) +
                                                          "_type_" + TString(std::to_string(WhichType)),
                                                      "", 500, 0, 5);
            h1_pion_q2[WhichInt][WhichType] = new TH1F("h1_pion_q2_interaction_" + TString(std::to_string(WhichInt)) +
                                                           "_type_" + TString(std::to_string(WhichType)),
                                                       "", 500, 0, 5);
            h1_pion_x[WhichInt][WhichType] = new TH1F("h1_pion_x_interaction_" + TString(std::to_string(WhichInt)) +
                                                          "_type_" + TString(std::to_string(WhichType)),
                                                      "", 100, 0, 1);
            h1_pion_p_npi[WhichInt][WhichType] =
                new TH1F("h1_Npion_p_interaction_" + TString(std::to_string(WhichInt)) + "_type_" +
                             TString(std::to_string(WhichType)),
                         "", 400, 0, 4);
            h1_pion_p_1pi[WhichInt][WhichType] =
                new TH1F("h1_1pion_p_interaction_" + TString(std::to_string(WhichInt)) + "_type_" +
                             TString(std::to_string(WhichType)),
                         "", 400, 0, 4);
        }
    }

    // Vector containing kinematic variables using Ecal
    vector<double> CalKineVars{};
    // Weight to fill the plots mentioned above
    // double LocalWeight;

    // Signal Event Counter -> 1e1p0pi events (everything lese is bkg)
    int SignalEventsplus = 0;
    int SignalEventszero = 0;
    int SignalEventsminus = 0;
    int eventremoved = 0; // Lataling: Delete this
    int totalpi0 = 0;     // also delete this
    // int QESignalEvents = 0;
    // int MECSignalEvents = 0;
    int RESSignalEventsplus = 0;
    int DISSignalEventsplus = 0;
    int OtherSignalEventsplus = 0;
    int RESSignalEventszero = 0;
    int DISSignalEventszero = 0;
    int OtherSignalEventszero = 0;
    int RESSignalEventsminus = 0;
    int DISSignalEventsminus = 0;
    int OtherSignalEventsminus = 0;

    // Lataling: define more Counters for truth/reconstructed study (for piminus)
    int tyrnDISm = 0; // This means, truth = yes (ty), reconstructed = no (rn), DIS mechanism
    int tyryDISm = 0; // This means, truth = yes (ty), reconstructed = yes (ry), DIS mechanism
    int tnrnDISm = 0; // etc..
    int tnryDISm = 0;
    int tyrnRESm = 0;
    int tyryRESm = 0;
    int tnrnRESm = 0;
    int tnryRESm = 0;
    int tyrnOTHm = 0;
    int tyryOTHm = 0;
    int tnrnOTHm = 0;
    int tnryOTHm = 0;

    // Lataling: define more Counters for truth/reconstructed study (for piplus)
    int tyrnDISp = 0; // This means, truth = yes (ty), reconstructed = no (rn), DIS mechanism
    int tyryDISp = 0; // This means, truth = yes (ty), reconstructed = yes (ry), DIS mechanism
    int tnrnDISp = 0; // etc..
    int tnryDISp = 0;
    int tyrnRESp = 0;
    int tyryRESp = 0;
    int tnrnRESp = 0;
    int tnryRESp = 0;
    int tyrnOTHp = 0;
    int tyryOTHp = 0;
    int tnrnOTHp = 0;
    int tnryOTHp = 0;

    // ---------------------------------------------------------------------------------------------------------------

    // Get the number of events to run overall

    //	int Nentries = TMath::Min(Ntfileentries,NtweightsEntries);

    // ---------------------------------------------------------------------------------------------------------------

    // Justification for the parameter fchoice
    // https://docs.google.com/presentation/d/1ghG08JfCYXRXh6O8hcXKrhJOFxkAs_9i5ZfoIkiiEHU/edit?usp=sharing
    // NOAH :: Need this (signal definition)
    TF1 *myElectronFit = new TF1("myElectronFit", "[0]+[1]/x", 0., 5.);

    if (en_beam[fbeam_en] == 1.161) {
        myElectronFit->SetParameters(17, 7);
    }
    if (en_beam[fbeam_en] == 2.261) {
        myElectronFit->SetParameters(16, 10.5);
    }
    if (en_beam[fbeam_en] == 4.461) {
        myElectronFit->SetParameters(13.5, 15);
    }

    // ---------------------------------------------------------------------------------------------------------------

    // Keeping track of the energy

    // ---------------------------------------------------------------------------------------------------------------

    /** Beginning of Event Loop **/

    int TotalCounter = 0;

    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        // Read Entry
        int nb = GetEntry(jentry);
        if (nb == 0) {
            std::cout << "Event loop: 0 byte read for entry " << jentry << ". Indicate failure in reading the file"
                      << std::endl;
        }

        if (jentry % 10000 == 0) {
            std::cout << jentry / 1000 << " k " << std::setprecision(3) << double(jentry) / double(nentries) * 100.
                      << " %" << std::endl;
        }

        /*
        if( jentry%200000 == 0 )
        {
                gDirectory->Write("hist_Files", TObject::kOverwrite);
                //cout<<jentry<<endl;
        }
        */

        TotalCounter++;

        // ---------------------------------------------------------------------------------------------------------------

        std::string StoreEnergy = fbeam_en;
        if (UsePhiThetaBand) {
            StoreEnergy = "";
        }

        // -------------------------------------------------------------------------------------------------------------------------

        // For GENIE samples, identify the interaction type
        // Lataling: copied this but with 1232 resonances separated out

        int Interaction = -1;
        int Separate_Interaction = -1;

        if (fchoice > 0) {
            if (qel) {
                Interaction = 1;
            }
            if (mec) {
                Interaction = 2;
            }
            if (res) {
                Interaction = 3;
            }
            if (dis) {
                Interaction = 4;
            }
        } else if (fchoice == 0) {
            Interaction = 0;
        } // smithja: added when I started running CLAS data on 10/29/2021

        if (fchoice > 0) {
            if (qel) {
                Separate_Interaction = 1;
            }
            if (mec) {
                Separate_Interaction = 2;
            }
            if (res && resid != 0) {
                Separate_Interaction = 3;
            }
            if (dis) {
                Separate_Interaction = 4;
            }
            if (res && resid == 0) {
                Separate_Interaction = 5;
            } // If type is delta 1232 resonance
        } else if (fchoice == 0) {
            Separate_Interaction = 0;
        } // smithja: added when I started running CLAS data on 10/29/2021

        // if (qel == false) continue;
        //  ---------------------------------------------------------------------------------------------------------------

        // For GENIE samples, identify the value for resc

        int resc_val = -99;

        if (fchoice > 0) {
            resc_val = resc[0]; // smithja: index always zero; resc looks to be single-element array in genie_analysis.h
        }

        // if (resc_val > 1) continue;

        // ---------------------------------------------------------------------------------------------------------------

        if (jentry == 0) { // first entry to initialize TorusCurrent, Fiducials and Subtraction classes

            // The TorusField has to be set before the Fiducialcut parameters are initialized
            if (en_beam[fbeam_en] > 1. &&
                en_beam[fbeam_en] < 2.) // 1.1 GeV, we are not using the 1.1 GeV data with 1500 current field
            {
                fTorusCurrent = 750;
            } else if ((en_beam[fbeam_en] > 2. && en_beam[fbeam_en] < 3.) ||
                       (en_beam[fbeam_en] > 4. && en_beam[fbeam_en] < 5.)) // 2.2 GeV	or 4.4 GeV
            {
                fTorusCurrent = 2250;
            } else {
                std::cout << "genie_analysis::Loop(): fTorusCurrent could not be assigned" << std::endl;
            }

            fiducialcut->SetConstants(fTorusCurrent, target_name, en_beam);
            fiducialcut->SetFiducialCutParameters(fbeam_en);
            std::cout << " EventLoop: Finished setting up fiducial cut class " << std::endl;
            //			rotation->InitSubtraction(fbeam_en, target_name, bind_en, N_tot, fiducialcut);
            rotation->InitSubtraction(StoreEnergy, target_name, bind_en, N_tot, fiducialcut);
            std::cout << " EventLoop: Finished setting up rotation initialize " << std::endl;
        }

        // Resets q vector to (0,0,0)
        rotation->ResetQVector();

        // -----------------------------------------------------------------------------------------------------------------------------------------------------------

        // Counters for truth level studies

        int TrueElectronsAboveThreshold = 0;
        int TrueProtonsAboveThreshold = 0;
        int TrueChargedPionsAboveThreshold = 0;
        int TruePiPlusAboveThreshold = 0;
        int TruePiMinusAboveThreshold = 0;
        int TrueGammasAboveThreshold = 0;

        // Define the counters
        int pionabove = 0;
        int pionbelow = 0;

        // Define temporary counters for how many pions there are in an event
        int true_pimin = 0;
        int reco_pimin = 0;
        int true_piplus = 0;
        int reco_piplus = 0;

        // -----------------------------------------------------------------------------------------------------------------------------------------------------------

        double SmearedPe;
        double SmearedEe;
        double e_acc_ratio = 1.; // will be 1 for CLAS data
        // Outgoing e',	Uncorr and corrected are the same read from root file.
        // V4_el and V3_el will be changed by smearing for GENIE simulation data
        TLorentzVector V4_el(pxl, pyl, pzl, El);
        TVector3 V3_el(pxl, pyl, pzl);

        double el_momentum = V3_el.Mag();
        double el_theta = V3_el.Theta();

        // ----------------------------------------------------------------------------------------------------------------------

        if (fchoice > 0) { // smearing, fiducials and acceptance ratio for GENIE simulation data

            // Smearing of Electron Vector from Simulation
            SmearedPe = gRandom->Gaus(pl, reso_e * pl);
            SmearedEe = sqrt(SmearedPe * SmearedPe + e_mass * e_mass);
            V3_el.SetXYZ(SmearedPe / pl * pxl, SmearedPe / pl * pyl, SmearedPe / pl * pzl);
            V4_el.SetPxPyPzE(V3_el.X(), V3_el.Y(), V3_el.Z(), SmearedEe);
            double phi_ElectronOut = V3_el.Phi(); // in Radians

            V3_el.SetPhi(
                phi_ElectronOut +
                TMath::Pi()); // Vec.Phi() is between (-180,180), GENIE coordinate system flipped with respect to CLAS

            //			//Fiducial Cuts with the smeared values // moved it further down after W & Q2 cuts
            //			if (ApplyFiducials)  { if (!EFiducialCut(fbeam_en,V3_el) ) continue;} // Electron theta & phi
            //fiducial cuts

            phi_ElectronOut += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS
            el_momentum = V3_el.Mag();      // Momentum after smearing
            el_theta = V3_el.Theta();       // Angle after smearing

            // acceptance_c takes phi in radians and here unmodified by 30 degree.

            e_acc_ratio = acceptance_c(el_momentum, cos(el_theta), phi_ElectronOut, 11, file_acceptance, ApplyAccWeights);
            if (fabs(e_acc_ratio) != e_acc_ratio) {
                cout << "HAHA" << endl;
                continue;
            }

            // --------------------------------------------------------------------------------------------------

            // GENIE Systematic Uncertainties

            //			tweights->GetEntry(jentry);
            //			float* ArrayWeights = weights->GetArray();

            ////			double TuningWeight = ArrayWeights[0]; // - 1 sigma variation
            ////			double TuningWeight = ArrayWeights[1]; // 0 sigma variation
            //			double TuningWeight = ArrayWeights[2]; // + 1 sigma variation

            //			e_acc_ratio = e_acc_ratio * TuningWeight;

            // --------------------------------------------------------------------------------------------------
        }

        // ----------------------------------------------------------------------------------------------------------------------

        double theta_min = myElectronFit->Eval(el_momentum); // in deg lataling: removed too
        if (el_theta * 180. / TMath::Pi() < theta_min) {
            continue;
        }

        if (fApplyThetaSliceEl) { // hard coded range for now
            if (el_theta * 180. / TMath::Pi() < t_thetaEl_lb->GetVal()) {
                continue;
            }
            if (el_theta * 180. / TMath::Pi() > t_thetaEl_ub->GetVal()) {
                continue;
            }
        }
        // ----------------------------------------------------------------------------------------------------------------------

        // Explicit cuts on electron momentum
        // Lataling: There are also some electron cuts here, which I'm commenting out for now

        // if (fbeam_en=="1161" && (el_momentum < 0.4 || (el_momentum < t_elMom_lb->GetVal() && fApplyElMomCut ==
        // true))) { continue; } if (fbeam_en=="2261" && (el_momentum < 0.55 || (el_momentum < t_elMom_lb->GetVal() &&
        // fApplyElMomCut == true))) { continue; } if (fbeam_en=="4461" && (el_momentum < 1.1 || (el_momentum <
        // t_elMom_lb->GetVal() && fApplyElMomCut == true))) { continue; }

        // Definition as for data. It is also correct for GENIE simulation data since V3_el is rotated above by 180
        // degree in phi
        double el_phi_mod = V3_el.Phi() * TMath::RadToDeg() + 30; // Add 30 degree for plotting and photon phi cut
        if (el_phi_mod < 0)
            el_phi_mod = el_phi_mod + 360; // Add 360 so that electron phi is between 0 and 360 degree

        if (fApplyPhiOpeningAngleEl) {
            if (!(TMath::Abs(el_phi_mod - 30) <= PhiOpeningAngleEl ||
                  TMath::Abs(el_phi_mod - 90) <= PhiOpeningAngleEl ||
                  TMath::Abs(el_phi_mod - 150) <= PhiOpeningAngleEl ||
                  TMath::Abs(el_phi_mod - 210) <= PhiOpeningAngleEl ||
                  TMath::Abs(el_phi_mod - 270) <= PhiOpeningAngleEl ||
                  TMath::Abs(el_phi_mod - 330) <= PhiOpeningAngleEl)) {
                continue;
            }
        }

        int ElectronSector = el_phi_mod / 60.;

        if (fApplyPhiSliceEl_Sectors16) {
            if ((ElectronSector != 5 && ElectronSector != 0))
                continue;
        }
        if (fApplyPhiSliceEl_Sectors126) {
            if (ElectronSector != 5 && ElectronSector != 1 && ElectronSector != 0)
                continue;
        }

        // Calculated Mott Cross Section and Weights for Inclusive Histograms
        // Wght and e_acc_ratio is 1 for CLAS data
        // double Mott_cross_sec = (
        // pow(fine_struc_const,2.)*(cos(el_theta)+1))/(2*pow(El,2.)*pow((1-cos(el_theta)),2.));

        double reco_Q2 = -(V4_el - V4_beam).Mag2();
        //		double Q4 = reco_Q2 * reco_Q2;
        double Q4 = 1.; // smithja: to eliminate Q4 dependence as part of re-doing Graham's Procedure 2 analysis, I
                        //          have set Q4 to simply be 1. The only instances where Q4 pops up are when weighting
                        //          some histograms, i.e. in the second argument of some Fill() statements. To put back
                        //          in the Q4 dependence in these weights, simply comment out this instance of Q4 and
                        //          uncomment the original one directly above.
        //		double Mott_cross_sec = (1./Q4) * XSecScale;
        double Mott_cross_sec =
            1; // smithja: this is the weight Dr. Betancourt said to do away with from Graham's analysis.
               //          the object of setting this to 1 is to remove the Q4 dependence from the
               //          Mott cross section. You can see that XSecScale is also part of the instance of
               //          Mott_cross_sec directly above. However, when I received this script, XSecScale
               //          was set to equal 1. While the above instance where Q4 = 1 technically does the
               //          same thing as this line, I am changing Mott_cross_sec here in case some analysis
               //          down the line wants to make use of XSecScale != 1.
        // ---------------------------------------------------------------------------------------------------------------------

        // ---------------------------------------------------------------------------------------------------------------------

        // Sanity check, especially for radiation
        if (wght < 0 || wght > 10) {
            std::cout << "Something is really wrong with your weights !!!" << std::endl;
        }

        double WeightIncl = wght * e_acc_ratio / Mott_cross_sec;

        // Securing ourselves against infinities
        if (fabs(WeightIncl) != WeightIncl) {
            continue;
        }

        // Calculation of Reconstructed Energy from ELectron only
        // using the same value of single nucleon separation E Ecal and Eqe
        /* double E_rec = */
            /* (m_prot * bind_en[ftarget] + m_prot * V4_el.E()) / (m_prot - V4_el.E() + V4_el.Rho() * cos(el_theta)); */
        /* double EQE_Reso = (E_rec - en_beam_Ecal[fbeam_en]) / en_beam_Ecal[fbeam_en]; */

        // Calculation of kinematic quantities (nu, Q2, x bjorken, q and W)
        double nu = -(V4_el - V4_beam).E();
        double x_bjk = reco_Q2 / (2 * m_prot * nu);

        // QE selection
        // if ( fabs(x_bjk - 1.) > 0.2) { continue; }

        // ---------------------------------------------------------------------------------------------------------------------
        TVector3 V3_q = (V4_beam - V4_el).Vect();
        double W_var = TMath::Sqrt((m_prot + nu) * (m_prot + nu) - V3_q * V3_q);

        // converting theta to degrees
        el_theta = el_theta * TMath::RadToDeg();

        // Lataling: commented this out to successfully remove the cut
        // Cuts on Q2 and W, only keep events with Q2 > Q2cut and W < Wcut
        // if ( reco_Q2 < Q2cut || W_var > Wcut) continue;

        // ---------------------------------------------------------------------------------------------------------------------

        // apapadop Nov 4 2020: true electron counter for truth level studies
        TrueElectronsAboveThreshold++;

        if (fchoice > 0) {

            // Fiducial Cuts with the smeared values
            if (ApplyFiducials) {
                if (!EFiducialCut(fbeam_en, V3_el))
                    continue;
            } // Electron theta & phi fiducial cuts
        }

        // ---------------------------------------------------------------------------------------------------------------------

        // Set q vector for the following rotations for the subtraction procedure
        rotation->SetQVector(V3_q);
        //		rotation->PrintQVector();
        // lataling: these are removed too
        // h2_el_theta_phi->Fill(el_phi_mod,el_theta,WeightIncl);

        // if (el_phi_mod > 0 && el_phi_mod < 60) {
        // h2_Electron_Theta_Momentum_FirstSector->Fill(V4_el.Rho(),V3_el.Theta()*180./TMath::Pi(),e_acc_ratio); } if
        // (el_phi_mod > 60 && el_phi_mod < 120) {
        // h2_Electron_Theta_Momentum_SecondSector->Fill(V4_el.Rho(),V3_el.Theta()*180./TMath::Pi(),e_acc_ratio); } if
        // (el_phi_mod > 120 && el_phi_mod < 180) {
        // h2_Electron_Theta_Momentum_ThirdSector->Fill(V4_el.Rho(),V3_el.Theta()*180./TMath::Pi(),e_acc_ratio); } if
        // (el_phi_mod > 180 && el_phi_mod < 240) {
        // h2_Electron_Theta_Momentum_FourthSector->Fill(V4_el.Rho(),V3_el.Theta()*180./TMath::Pi(),e_acc_ratio); } if
        // (el_phi_mod > 240 && el_phi_mod < 300) {
        // h2_Electron_Theta_Momentum_FifthSector->Fill(V4_el.Rho(),V3_el.Theta()*180./TMath::Pi(),e_acc_ratio); } if
        // (el_phi_mod > 300 && el_phi_mod < 360) {
        // h2_Electron_Theta_Momentum_SixthSector->Fill(V4_el.Rho(),V3_el.Theta()*180./TMath::Pi(),e_acc_ratio); }

        // ---------------------------------------------------------------------------------------------------------------------

        // apapadop: Oct 8 2020: ditching bad sectors
        // Counting sectors from 0 to 5

        if (!UseAllSectors) {
            if ((ElectronSector == 2 || ElectronSector == 4) && fbeam_en == "1161") {
                continue;
            }
            if ((ElectronSector == 2 || ElectronSector == 3 || ElectronSector == 4) && fbeam_en == "2261") {
                continue;
            }
        }

        // ---------------------------------------------------------------------------------------------------------------------

        // Fully inclusive plots & counters shown first below


        // Index variables for hadrons (p and pions)
        int index_pi[20];    // index for each pion
        int charge_pi[20];   // Charge for the pions and photons
        // Smeared Momentum and Energy values for GENIE (simulation) data
        double Smeared_Ppi[20]; // smeared momentum values for pions

        // Number of hadrons
        int num_p = 0;
        int num_pi = 0;
        int num_pi_phot = 0; // couting all pions and photons
        // Index and number variables for neutral particles
        bool ec_radstat_n[20];

        // Array initialize to -1 or false
        for (int i = 0; i < 20; i++) {
            index_pi[i] = -1;
            ec_radstat_n[i] = false;
            charge_pi[i] = -2; // default number should be not a possible real charge
            Smeared_Ppi[i] = 0;
        }

        const double phot_rad_cut = 40;
        const double phot_e_phidiffcut = 30; // electron - photon phi difference cut

        // Creating vectors to store id of particles in the array
        vector<int> ProtonID;
        vector<int> PiPlusID;
        vector<int> PiMinusID;
        vector<int> PhotonID;
        ProtonID.clear();
        PiPlusID.clear();
        PiMinusID.clear();
        PhotonID.clear();

        // Loop for Hadrons
        for (int i = 0; i < nf; i++) {
            // lataling: add basic pion counter here:
            if (pdgf[i] == 211) {
                true_piplus = true_piplus + 1;
            }
            if (pdgf[i] == -211) {
                true_pimin = true_pimin + 1;
            }

            if (Applymomthresh ? pdgf[i] == -211 && pf[i] > 0.15
                               : pdgf[i] == -211) { // Pi minus Lataling: here (0.15 GeV)
                double PiMinusWeight = 1.;
                double PiMinusPhi_Deg = -999.;
                double PiMinusTheta_Deg = -999.;
                double PiMinusCosTheta = -999.;
                double PiMinusMag = -999.;

                if (fchoice > 0) { // GENIE data
                    // Smearing of pi minus
                    double temp_smear_P = gRandom->Gaus(pf[i], reso_pi * pf[i]);
                    double temp_smear_E = sqrt(temp_smear_P * temp_smear_P + m_pion * m_pion);

                    TVector3 V3_pi_corr(temp_smear_P / pf[i] * pxf[i], temp_smear_P / pf[i] * pyf[i],
                                        temp_smear_P / pf[i] * pzf[i]);
                    double phi_pion = V3_pi_corr.Phi();
                    V3_pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)

                    // apapadop Nov 4 2020: true charged pion counter for truth level studies above a min theta
                    // threshold given by a functional form A + B / P

                    if (PimiFiducialCutExtra(StoreEnergy, V3_pi_corr)) {

                        TrueChargedPionsAboveThreshold++;
                        TruePiMinusAboveThreshold++;
                    }

                    // Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus
                    //					if (ApplyFiducials) { if ( !Pi_phot_fid_united(fbeam_en, V3_pi_corr, -1) ) {  continue;
                    //} }
                    if (ApplyFiducials) {
                        if (!Pi_phot_fid_united(StoreEnergy, V3_pi_corr, -1)) {
                            pionbelow = pionbelow + 1;
                            continue;
                        }
                    } // LAtaling: if fail fiducial then pion_below adds 1

                    num_pi = num_pi + 1;
                    num_pi_phot = num_pi_phot + 1;
                    // index_pi[num_pi_phot - 1] = i; Lataling: This is to stop pi0 taking the place of a charged pion
                    // at index 0
                    index_pi[num_pi - 1] = i;
                    PiMinusID.push_back(i);
                    // charge_pi[num_pi_phot - 1] = 1; Same reason as above
                    charge_pi[num_pi - 1] = -1;
                    Smeared_Ppi[num_pi - 1] = temp_smear_P;

                    phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS
                    PiMinusCosTheta = V3_pi_corr.CosTheta();
                    PiMinusMag = V3_pi_corr.Mag();

                    // acceptance_c takes phi in radians and here unmodified by 30 degree.
                    //					PiMinusWeight = wght * acceptance_c(PiMinusMag,PiMinusCosTheta,
                    //phi_pion, -211,file_acceptance,ApplyAccWeights);

                    // UNTIL AXEL CREATES THE CORRECT PIMINUS MAP, WE SET THE PIMINUS ACCEPTANCE TO BE 1
                    //					PiMinusWeight = wght * acceptance_c(PiMinusMag,PiMinusCosTheta,
                    //phi_pion, -211,file_acceptance,false);
                    PiMinusWeight =
                        wght *
                        acceptance_c(
                            PiMinusMag, PiMinusCosTheta, phi_pion, -211, file_acceptance_pim,
                            ApplyAccWeights); // lataling: this was set to true, for some reason. I have changed it

                    if (fabs(PiMinusWeight) != PiMinusWeight) {
                        continue;
                    }

                    PiMinusPhi_Deg = V3_pi_corr.Phi() * 180. / TMath::Pi() + 180. + 30.;
                    if (PiMinusPhi_Deg > 360.) {
                        PiMinusPhi_Deg -= 360.;
                    }
                    PiMinusTheta_Deg = V3_pi_corr.Theta() * 180. / TMath::Pi();
                } else { // CLAS data does not need Fiducial Cut again
                    num_pi = num_pi + 1;
                    num_pi_phot = num_pi_phot + 1;
                    index_pi[num_pi_phot - 1] = i;
                    PiMinusID.push_back(i);
                    charge_pi[num_pi_phot - 1] = -1;

                    TVector3 V3_pi_corr(pxf[i], pyf[i], pzf[i]);

                    PiMinusPhi_Deg = V3_pi_corr.Phi() * 180. / TMath::Pi() + 180. + 30.;
                    if (PiMinusPhi_Deg > 360.) {
                        PiMinusPhi_Deg -= 360.;
                    }
                    PiMinusTheta_Deg = V3_pi_corr.Theta() * 180. / TMath::Pi();
                    PiMinusMag = V3_pi_corr.Mag();
                }
                // lataling: these things removed -- useful to make acceptance maps
                // h1_PiMinus_AccMapWeights->Fill(PiMinusWeight);
                // h1_PiMinus_Momentum->Fill(PiMinusMag,PiMinusWeight);

                // h3_PiMinus_Mom_Theta_Phi->Fill(PiMinusMag,PiMinusTheta_Deg,PiMinusPhi_Deg,PiMinusWeight);
            }

            // -------------------------------------------------------------------------------------------------------------------

            if (Applymomthresh ? pdgf[i] == 211 && pf[i] > 0.15
                               : pdgf[i] == 211) { // Lataling: Here as well
                double PiPlusWeight = 1.;
                double PiPlusPhi_Deg = -999.;
                double PiPlusTheta_Deg = -999.;
                double PiPlusCosTheta = -999.;
                double PiPlusMag = -999.;

                if (fchoice > 0) { // GENIE data
                    // Smearing of pi plus
                    double temp_smear_P = gRandom->Gaus(pf[i], reso_pi * pf[i]); // Smearing is momentum dependent
                    double temp_smear_E = sqrt(temp_smear_P * temp_smear_P + m_pion * m_pion);

                    TVector3 V3_pi_corr(temp_smear_P / pf[i] * pxf[i], temp_smear_P / pf[i] * pyf[i],
                                        temp_smear_P / pf[i] * pzf[i]);
                    double phi_pion = V3_pi_corr.Phi();
                    V3_pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)

                    // apapadop Nov 4 2020: true charged pion counter for truth level studies with min theta threshold
                    // (12 deg)

                    if (PiplFiducialCutExtra(StoreEnergy, V3_pi_corr)) {

                        TrueChargedPionsAboveThreshold++;
                        TruePiPlusAboveThreshold++;
                    }

                    // Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus
                    //					if (ApplyFiducials) { if ( !Pi_phot_fid_united(fbeam_en, V3_pi_corr, 1) )     {
                    //continue; }
                    //}
                    if (ApplyFiducials) {
                        if (!Pi_phot_fid_united(StoreEnergy, V3_pi_corr, 1)) {
                            pionbelow = pionbelow + 1;
                            continue;
                        }
                    } // LAtaling: if fail fiducial then pion_below adds 1

                    num_pi = num_pi + 1;
                    num_pi_phot = num_pi_phot + 1;
                    // index_pi[num_pi_phot - 1] = i; Lataling: This is to stop pi0 taking the place of a charged pion
                    // at index 0
                    index_pi[num_pi - 1] = i;
                    PiPlusID.push_back(i);
                    // charge_pi[num_pi_phot - 1] = 1; Same reason as above
                    charge_pi[num_pi - 1] = 1;
                    Smeared_Ppi[num_pi - 1] = temp_smear_P;

                    phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS
                    PiPlusCosTheta = V3_pi_corr.CosTheta();
                    PiPlusMag = V3_pi_corr.Mag();

                    // acceptance_c takes phi in radians and here unmodified by 30 degree.
                    PiPlusWeight = wght * acceptance_c(PiPlusMag, PiPlusCosTheta, phi_pion, 211, file_acceptance_pip,
                                                       ApplyAccWeights);
                    if (fabs(PiPlusWeight) != PiPlusWeight) {
                        continue;
                    }

                    PiPlusPhi_Deg = V3_pi_corr.Phi() * 180. / TMath::Pi() + 180. + 30.;
                    if (PiPlusPhi_Deg > 360.) {
                        PiPlusPhi_Deg -= 360.;
                    }
                    PiPlusTheta_Deg = V3_pi_corr.Theta() * 180. / TMath::Pi();

                } else { // CLAS data does not need Fiducial Cut again
                    num_pi = num_pi + 1;
                    num_pi_phot = num_pi_phot + 1;
                    PiPlusID.push_back(i);
                    charge_pi[num_pi_phot - 1] = 1;

                    TVector3 V3_pi_corr(pxf[i], pyf[i], pzf[i]);

                    PiPlusPhi_Deg = V3_pi_corr.Phi() * 180. / TMath::Pi() + 180. + 30.;
                    if (PiPlusPhi_Deg > 360.) {
                        PiPlusPhi_Deg -= 360.;
                    }
                    PiPlusTheta_Deg = V3_pi_corr.Theta() * 180. / TMath::Pi();
                    PiPlusMag = V3_pi_corr.Mag();
                }

                // lataling: likewise, these too are removed
                // h1_PiPlus_AccMapWeights->Fill(PiPlusWeight);
                // h1_PiPlus_Momentum->Fill(PiPlusMag,PiPlusWeight);

                // h3_PiPlus_Mom_Theta_Phi->Fill(PiPlusMag,PiPlusTheta_Deg,PiPlusPhi_Deg,PiPlusWeight);
            }

            // ---------------------------------------------------------------------------------------------------------------------------

            // lataling: here too (0.3) -- photons have to be > 0.3 GeV
            if (Applymomthresh ? pdgf[i] == 22 && pf[i] > 0.3
                               : pdgf[i] == 22) {
                // Determine photon vector for the cut on radiation photon via angle with respect to the electron
                TVector3 V3_phot_angles(pxf[i], pyf[i], pzf[i]);
                if (fchoice > 0) { // GENIE data
                    // no smearing of GENIE photons
                    double phi_photon = V3_phot_angles.Phi();
                    V3_phot_angles.SetPhi(phi_photon + TMath::Pi()); // Vec.Phi() is between (-180,180)

                    // apapadop Nov 4 2020: true photon counter for truth level studies

                    if (Phot_fidExtra(V3_phot_angles)) {

                        TrueGammasAboveThreshold++;
                    }

                    if (ApplyFiducials) {
                        if (!Pi_phot_fid_united(fbeam_en, V3_phot_angles, 0)) {
                            continue;
                        }
                    }
                }

                double neut_phi_mod = V3_phot_angles.Phi() * TMath::RadToDeg() + 30; // Add 30 degree
                if (neut_phi_mod < 0)
                    neut_phi_mod = neut_phi_mod + 360; // Neutral particle is between 0 and 360 degree

                num_pi_phot = num_pi_phot + 1;
                PhotonID.push_back(i);

                Smeared_Ppi[num_pi_phot - 1] = V3_phot_angles.Mag();

                // lataling: removed too
                // CosDeltaThetaElectronPhotonAboveThreshold->Fill( cos( V3_phot_angles.Angle(V3_el) ) );
                // CosDeltaPhiElectronPhotonAboveThreshold->Fill( cos( neut_phi_mod-el_phi_mod*TMath::Pi()/180. ) );

                // within 40 degrees in theta and 30 degrees in phi. Electron phi has already added 30 degree and
                // between 0 to 360

                if (V3_phot_angles.Angle(V3_el) * TMath::RadToDeg() < phot_rad_cut &&
                    fabs(neut_phi_mod - el_phi_mod) < phot_e_phidiffcut) {

                    ec_radstat_n[num_pi_phot - 1] = true; // select radiation photons

                    // Skip event if there is at least one radiation photon
                    continue;
                }

                // not used anymore
                if (!ec_radstat_n[num_pi_phot - 1]) {
                    double GammaPhi_Deg = neut_phi_mod;
                    double GammaTheta_Deg = V3_phot_angles.Theta() * TMath::RadToDeg();
                    double GammaMag = V3_phot_angles.Mag();
                    double GammaWeight = wght;
                }
            }
        } // end of hadron loop

        // ----------------------------------------------------------------------------------------------------------------------------

        // TODO: Keep this!
        // apapadop: executive decision Dec 3
        // given that genie has much higher proton multiplicities than we observe in data
        // we ignore the num_p > 4 cases

        // if (num_p > 4) { continue; }



        // ---------------------------------------------------------------------------------------------------------------------------------------

        // lataling: this has also been removed

        // Lataling: for counting events only, this statement looks for events in which any number of pions are produced
        // if (num_pi != 0 || pionbelow != 0){
        if (num_pi != 0) {
            // Define the counters
            // int pionabove = 0;
            // int pionbelow = 0;

            // Count the pions and make temporary vector
            for (int i = 0; i < num_pi; i++) {
                // if (pf[index_pi[i]] > 0.15 && charge_pi[i] == -1) {pionabove = pionabove + 1;}
                // else if (pf[index_pi[i]] < 0.15 && charge_pi[i] == -1) {pionbelow = pionbelow + 1;}
                // else if (pf[index_pi[i]] > 0.15 && charge_pi[i] == 1) {pionabove = pionabove + 1;}
                // else if (pf[index_pi[i]] < 0.15 && charge_pi[i] == 1) {pionbelow = pionbelow + 1;}
                if (charge_pi[i] == -1) {
                    reco_pimin = reco_pimin + 1;
                } // Adds to the pi minus counter
                else if (charge_pi[i] == 1) {
                    reco_piplus = reco_piplus + 1;
                } // Adds to the pi plus counter
                  // h1_pion_p_npi[Separate_Interaction][charge_pi[i]+1]->Fill(pf[index_pi[i]]); //fill the histogram
            }

            // Change counters to 3 if greater than 3
            if (pionabove > 3) {
                pionabove = 3;
            }
            if (pionbelow > 3) {
                pionbelow = 3;
            }

            // Fill the histogram
            h2_threshold_passing[Separate_Interaction]->Fill(pionabove + 1, pionbelow + 1);
        }

        //----------------------------- e- ,1pi  -----------------------------------------
        // With NumOfProton == 0, this is 1e and 1pi ONLY (i.e., fully exclusive) [at the moment: maybe other exotic
        // particles?] With NumOfProton != 0, this is 1e, 1pi and and the selected number of protons [at the moment:
        // maybe other exotic particles?] If detector_acceptance == 0, then the unchanged values are used, and if det...
        // == 1, then the detector idiosynchracies are included This does re-invent the wheel somewhat as events such as
        // 2p1pi events are already counted above. But the code is so long and most of it irrelevant, that it's easiest
        // to just re-do it here.

        // Starting off by choosing how many protons (if number of protons is set to -1, then all protons are included)
        // Note: this used say if (num_pi_phot == 1), but now I'm only looking for charged pions as pi0 aren't decayed
        // by GENIE Thus, I have changed it to num_pi instead, which only looks at the number of pions
        if (NumOfProton != -1 ? num_pi == 1 && num_p == NumOfProton : num_pi == 1) {
            TVector3 V3_pi_corr;
            double P_undet = 0;
            double pion_acc_ratio = 1;

            int pluscount = 0;
            int mincount = 0;
            int zerocount = 0;

            // SignalEvents++;
            //  std::cout << "SignalEvents: " << SignalEvents << std::endl; // smithja: used to track the number of
            //  events in the 1e1p0pi case

            // if (Separate_Interaction == 5) { RESSignalEvents++; }
            // else if (Separate_Interaction == 4) { DISSignalEvents++; }
            // else { OtherSignalEvents++; }

            // for (int i = 0; i < num_pi_phot; i++){
            if (charge_pi[0] == -1) {
                mincount = mincount + 1;
            }
            // else if (charge_pi[0] == 0){zerocount = zerocount + 1;}
            else if (charge_pi[0] == 1) {
                pluscount = pluscount + 1;
            }
            //}

            if (pluscount > 0) {
                SignalEventsplus = SignalEventsplus + 1;
            }
            if (mincount > 0) {
                SignalEventsminus = SignalEventsminus + 1;
            }
            // if (zerocount > 0){SignalEventszero = SignalEventszero + 1;}

            if (pluscount > 0 && Separate_Interaction == 5) {
                RESSignalEventsplus = RESSignalEventsplus + 1;
            }
            if (mincount > 0 && Separate_Interaction == 5) {
                RESSignalEventsminus = RESSignalEventsminus + 1;
            }
            // if (zerocount > 0 && Separate_Interaction == 5){RESSignalEventszero = RESSignalEventszero + 1;}

            if (pluscount > 0 && Separate_Interaction == 4) {
                DISSignalEventsplus = DISSignalEventsplus + 1;
            }
            if (mincount > 0 && Separate_Interaction == 4) {
                DISSignalEventsminus = DISSignalEventsminus + 1;
            }
            // if (zerocount > 0 && Separate_Interaction == 4){DISSignalEventszero = DISSignalEventszero + 1;}

            if (pluscount > 0 && Separate_Interaction != 4 && Separate_Interaction != 5) {
                OtherSignalEventsplus = OtherSignalEventsplus + 1;
            }
            if (mincount > 0 && Separate_Interaction != 4 && Separate_Interaction != 5) {
                OtherSignalEventsminus = OtherSignalEventsminus + 1;
            }
            // if (zerocount > 0 && Separate_Interaction != 4 && Separate_Interaction != 5){OtherSignalEventszero =
            // OtherSignalEventszero + 1;}

            if (fchoice == 0) { // CLAS data
                V3_pi_corr.SetXYZ(pxf[index_pi[0]], pyf[index_pi[0]], pzf[index_pi[0]]);
                pion_acc_ratio = 1; // acceptance is 1 for CLAS data
            }

            if (fchoice > 0) {      // GENIE data
                pion_acc_ratio = 0; // reset just to be sure
                // Lataling: If & else statements toggle detector induced smearing on or off
                if (detector_acceptance == 1) {
                    V3_pi_corr.SetXYZ(Smeared_Ppi[0] / pf[index_pi[0]] * pxf[index_pi[0]],
                                      Smeared_Ppi[0] / pf[index_pi[0]] * pyf[index_pi[0]],
                                      Smeared_Ppi[0] / pf[index_pi[0]] * pzf[index_pi[0]]);
                } else {
                    V3_pi_corr.SetXYZ(pxf[index_pi[0]], pyf[index_pi[0]], pzf[index_pi[0]]);
                }
                double phi_pion = V3_pi_corr.Phi();        // in Radians
                V3_pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)
                phi_pion += TMath::Pi();                   // GENIE coordinate system flipped with respect to CLAS

                double pion_theta = V3_pi_corr.Theta();
                double pion_mom_corr = V3_pi_corr.Mag();

                if (charge_pi[0] == 1) { // acceptance for pi plus and introducing the cuts
                    pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip,
                                                  ApplyAccWeights);
                    if (fabs(pion_acc_ratio) != pion_acc_ratio) {
                        continue;
                    }
                } else if (charge_pi[0] == -1) { // acceptance for pi minus. using electron acceptance map

                    // pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211,
                    // file_acceptance,ApplyAccWeights);

                    // UNTIL AXEL CREATES THE CORRECT PIMINUS MAP, WE SET THE PIMINUS ACCEPTANCE TO BE 1
                    //					pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion,
                    //-211, file_acceptance,false);
                    pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance_pim,
                                                  ApplyAccWeights);

                    if (fabs(pion_acc_ratio) != pion_acc_ratio) {
                        continue;
                    }
                } else if (charge_pi[0] == 0) { // acceptance for photon/pi0 is 1 for now F.H. 09/24/19
                    pion_acc_ratio = 1;
                } else {
                    std::cout << "WARNING: 1 Pion Events. pion_acc_ratio is still 0. Continue with next event "
                              << std::endl;
                    continue;
                }

                // histoweights is 1/Mott_cross_sec for CLAS data -- redefined here by lataling as a "fudge factor"
                // pion_acc_ratio = 1; //delete when needed
                double histoweights = pion_acc_ratio * WeightIncl;
                float noweight = 1.;

                /*
                //lataling: These cuts have been determined with Optimal_Cut.C
                if (charge_pi[0] == 1 && W_var < 1.405){
                histoweights = pion_acc_ratio * WeightIncl;
                noweight = 1.;}

                else if (charge_pi[0] == -1 && W_var < 1.445){
                histoweights = pion_acc_ratio * WeightIncl;
                noweight = 1.;}

                //else if (charge_pi[0] == 0 && W_var < 1.315){
                //histoweights = pion_acc_ratio * WeightIncl;
                //noweight = 1.;}

                else {
                        histoweights = 0;
                        noweight = 0;
                }
                */

                if (detector_acceptance == 1) {
                    // Lataling: doing the same as above but better (unless it doesn't work... in which case, it is
                    // worse)
                    // h2_pion_mom_theta[Separate_Interaction][charge_pi[0]+1]->Fill(pion_mom_corr,pion_theta*(180/TMath::Pi()),histoweights);
                    // h2_pion_mom_phi[Separate_Interaction][charge_pi[0]+1]->Fill(pion_mom_corr,phi_pion*(180/TMath::Pi()),histoweights);
                    h1_pion_momentum[Separate_Interaction][charge_pi[0] + 1]->Fill(pion_mom_corr, histoweights);
                    h1_pion_theta[Separate_Interaction][charge_pi[0] + 1]->Fill(pion_theta * (180 / TMath::Pi()),
                                                                                histoweights);
                    h1_pion_phi[Separate_Interaction][charge_pi[0] + 1]->Fill(phi_pion * (180 / TMath::Pi()),
                                                                              histoweights);
                    h1_pion_omega[Separate_Interaction][charge_pi[0] + 1]->Fill(nu, histoweights);
                    h1_pion_electrontheta[Separate_Interaction][charge_pi[0] + 1]->Fill(el_theta, histoweights);
                    h1_pion_electronmomentum[Separate_Interaction][charge_pi[0] + 1]->Fill(el_momentum, histoweights);
                    h1_pion_q2[Separate_Interaction][charge_pi[0] + 1]->Fill(reco_Q2, histoweights);
                    h1_pion_w[Separate_Interaction][charge_pi[0] + 1]->Fill(W_var, histoweights);
                    h1_pion_x[Separate_Interaction][charge_pi[0] + 1]->Fill(x_bjk, histoweights);
                } else {

                    // h2_pion_mom_theta[Separate_Interaction][charge_pi[0]+1]->Fill(pion_mom_corr,pion_theta*(180/TMath::Pi()),1);
                    // h2_pion_mom_phi[Separate_Interaction][charge_pi[0]+1]->Fill(pion_mom_corr,phi_pion*(180/TMath::Pi()),1);
                    h1_pion_momentum[Separate_Interaction][charge_pi[0] + 1]->Fill(pion_mom_corr, noweight);
                    h1_pion_theta[Separate_Interaction][charge_pi[0] + 1]->Fill(pion_theta * (180 / TMath::Pi()),
                                                                                noweight);
                    h1_pion_phi[Separate_Interaction][charge_pi[0] + 1]->Fill(phi_pion * (180 / TMath::Pi()), noweight);
                    h1_pion_omega[Separate_Interaction][charge_pi[0] + 1]->Fill(nu, noweight);
                    h1_pion_electrontheta[Separate_Interaction][charge_pi[0] + 1]->Fill(el_theta, noweight);
                    h1_pion_electronmomentum[Separate_Interaction][charge_pi[0] + 1]->Fill(el_momentum, noweight);
                    h1_pion_q2[Separate_Interaction][charge_pi[0] + 1]->Fill(reco_Q2, noweight);
                    h1_pion_w[Separate_Interaction][charge_pi[0] + 1]->Fill(W_var, noweight);
                    h1_pion_x[Separate_Interaction][charge_pi[0] + 1]->Fill(x_bjk, noweight);
                    h1_pion_p_1pi[Separate_Interaction][charge_pi[0] + 1]->Fill(pf[index_pi[0]]); // fill the histogram
                }
            }

            rotation->pi1_rot_func(V3_pi_corr, charge_pi[0], &P_undet);

            // histoweight is 1/Mott_cross_sec for CLAS data
            double histoweight = pion_acc_ratio * WeightIncl;
        }

        //----------------------------- e- ,2pi  -----------------------------------------

        // Add the heatmap values
        if (Separate_Interaction == 4) { // DIS
            if (true_pimin == 0 &&
                true_piplus != 0) { // Everything in this if statement is added exclusively to the pi+ group
                if (true_piplus != 1 && reco_piplus != 1) {
                    tnrnDISp = tnrnDISp + 1;
                }
                if (true_piplus != 1 && reco_piplus == 1) {
                    tnryDISp = tnryDISp + 1;
                }
                if (true_piplus == 1 && reco_piplus != 1) {
                    tyrnDISp = tyrnDISp + 1;
                }
                if (true_piplus == 1 && reco_piplus == 1) {
                    tyryDISp = tyryDISp + 1;
                }
            } else if (true_piplus == 0 &&
                       true_pimin != 0) { // Everything in this if statement is added exclusively to the pi- group
                if (true_pimin != 1 && reco_pimin != 1) {
                    tnrnDISm = tnrnDISm + 1;
                }
                if (true_pimin != 1 && reco_pimin == 1) {
                    tnryDISm = tnryDISm + 1;
                }
                if (true_pimin == 1 && reco_pimin != 1) {
                    tyrnDISm = tyrnDISm + 1;
                }
                if (true_pimin == 1 && reco_pimin == 1) {
                    tyryDISm = tyryDISm + 1;
                }
            } else if (true_piplus != 0 &&
                       true_pimin != 0) { // Everything in this if statement is added to both the pi+ and pi- groups
                if ((reco_pimin + reco_piplus == 1) && reco_pimin == 1) {
                    tnryDISm = tnryDISm + 1;
                } else if ((reco_pimin + reco_piplus == 1) && reco_piplus == 1) {
                    tnryDISp = tnryDISp + 1;
                } else {
                    tnrnDISm = tnrnDISm + 1;
                    tnrnDISp = tnrnDISp + 1;
                }
            }
        }

        else if (Separate_Interaction == 5) { // RES
            if (true_pimin == 0 &&
                true_piplus != 0) { // Everything in this if statement is added exclusively to the pi+ group
                if (true_piplus != 1 && reco_piplus != 1) {
                    tnrnRESp = tnrnRESp + 1;
                }
                if (true_piplus != 1 && reco_piplus == 1) {
                    tnryRESp = tnryRESp + 1;
                }
                if (true_piplus == 1 && reco_piplus != 1) {
                    tyrnRESp = tyrnRESp + 1;
                }
                if (true_piplus == 1 && reco_piplus == 1) {
                    tyryRESp = tyryRESp + 1;
                }
            } else if (true_piplus == 0 &&
                       true_pimin != 0) { // Everything in this if statement is added exclusively to the pi- group
                if (true_pimin != 1 && reco_pimin != 1) {
                    tnrnRESm = tnrnRESm + 1;
                }
                if (true_pimin != 1 && reco_pimin == 1) {
                    tnryRESm = tnryRESm + 1;
                }
                if (true_pimin == 1 && reco_pimin != 1) {
                    tyrnRESm = tyrnRESm + 1;
                }
                if (true_pimin == 1 && reco_pimin == 1) {
                    tyryRESm = tyryRESm + 1;
                }
            } else if (true_piplus != 0 &&
                       true_pimin != 0) { // Everything in this if statement is added to both the pi+ and pi- groups
                if (reco_pimin + reco_piplus == 1) {
                    tnryRESm = tnryRESm + 1;
                    tnryRESp = tnryRESp + 1;
                } else {
                    tnrnRESm = tnrnRESm + 1;
                    tnrnRESp = tnrnRESp + 1;
                }
            }
        }

        else { // other
            if (true_pimin == 0 &&
                true_piplus != 0) { // Everything in this if statement is added exclusively to the pi+ group
                if (true_piplus != 1 && reco_piplus != 1) {
                    tnrnOTHp = tnrnOTHp + 1;
                }
                if (true_piplus != 1 && reco_piplus == 1) {
                    tnryOTHp = tnryOTHp + 1;
                }
                if (true_piplus == 1 && reco_piplus != 1) {
                    tyrnOTHp = tyrnOTHp + 1;
                }
                if (true_piplus == 1 && reco_piplus == 1) {
                    tyryOTHp = tyryOTHp + 1;
                }
            } else if (true_piplus == 0 &&
                       true_pimin != 0) { // Everything in this if statement is added exclusively to the pi- group
                if (true_pimin != 1 && reco_pimin != 1) {
                    tnrnOTHm = tnrnOTHm + 1;
                }
                if (true_pimin != 1 && reco_pimin == 1) {
                    tnryOTHm = tnryOTHm + 1;
                }
                if (true_pimin == 1 && reco_pimin != 1) {
                    tyrnOTHm = tyrnOTHm + 1;
                }
                if (true_pimin == 1 && reco_pimin == 1) {
                    tyryOTHm = tyryOTHm + 1;
                }
            } else if (true_piplus != 0 &&
                       true_pimin != 0) { // Everything in this if statement is added to both the pi+ and pi- groups
                if (reco_pimin + reco_piplus == 1) {
                    tnryOTHm = tnryOTHm + 1;
                    tnryOTHp = tnryOTHp + 1;
                } else {
                    tnrnOTHm = tnrnOTHm + 1;
                    tnrnOTHp = tnrnOTHp + 1;
                }
            }
        }

    } // end of event loop (jentry)

    gStyle->SetOptFit(1);

    //------------------------------------fractional energy reconstruction plots --------------------------------------

    gDirectory->Write("hist_Files", TObject::kOverwrite);
    // skim_tree->AutoSave();

    // --------------------------------------------------------------------------------------------------------

    // Lataling: This has also been changed
    std::cout << std::endl
              << "-----------------------------------------------------------------------------------------------------"
              << std::endl;
    std::cout << std::endl << "# Processed Events = " << TotalCounter << std::endl;
    std::cout << std::endl << "1Pi plus Signal # Events = " << SignalEventsplus << std::endl;
    std::cout << std::endl << "1Pi zero Signal # Events = " << SignalEventszero << std::endl;
    std::cout << std::endl << "1Pi minus Signal # Events = " << SignalEventsminus << std::endl;
    // std::cout << std::endl << "Passing Rate = " << int(double(SignalEvents) / double(TotalCounter)*100.) << " \%"<<
    // std::endl << std::endl;

    if (fchoice > 0) {

        // std::cout << std::endl << "QE Fractional Contribution = " << int(double(QESignalEvents) /
        // double(SignalEvents)*100.) << " \%" << std::endl; std::cout << std::endl << "MEC Fractional Contribution = "
        // << int(double(MECSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl; std::cout << std::endl <<
        // "RES Fractional Contribution = " << int(double(RESSignalEvents) / double(SignalEvents)*100.) << " \%" <<
        // std::endl; std::cout << std::endl << "DIS Fractional Contribution = " << int(double(DISSignalEvents) /
        // double(SignalEvents)*100.) << " \%" << std::endl;
        std::cout << std::endl << "1232 RES plus Contribution = " << RESSignalEventsplus << std::endl;
        std::cout << std::endl << "DIS plus Contribution = " << DISSignalEventsplus << std::endl;
        std::cout << std::endl << "Other plus Contribution = " << OtherSignalEventsplus << std::endl;
        std::cout << std::endl << "1232 RES zero Contribution = " << RESSignalEventszero << std::endl;
        std::cout << std::endl << "DIS zero Contribution = " << DISSignalEventszero << std::endl;
        std::cout << std::endl << "Other zero Contribution = " << OtherSignalEventszero << std::endl;
        std::cout << std::endl << "1232 RES minus Contribution = " << RESSignalEventsminus << std::endl;
        std::cout << std::endl << "DIS minus Contribution = " << DISSignalEventsminus << std::endl;
        std::cout << std::endl << "Other minus Contribution = " << OtherSignalEventsminus << std::endl;
        std::cout
            << std::endl
            << "-----------------------------------------------------------------------------------------------------"
            << std::endl;
    }

    ofstream output_file("./" + (fchoice > 0 ? std::to_string(fchoice) : "") + "_" +
                         std::string(detector_acceptance == 0 ? "noAcc" : "Acc") + "1piSelection.txt");

    output_file
        << std::endl
        << "-----------------------------------------------------------------------------------------------------"
        << std::endl;
    output_file << std::endl << "# Processed Events = " << TotalCounter << std::endl;
    output_file << std::endl << "1Pi plus Signal # Events = " << SignalEventsplus << std::endl;
    output_file << std::endl << "1Pi zero Signal # Events = " << SignalEventszero << std::endl;
    output_file << std::endl << "1Pi minus Signal # Events = " << SignalEventsminus << std::endl;
    // output_file << std::endl << "Passing Rate = " << int(double(SignalEvents) / double(TotalCounter)*100.) << " \%"<<
    // std::endl << std::endl;

    if (fchoice > 0) {

        // output_file << std::endl << "QE Fractional Contribution = " << int(double(QESignalEvents) /
        // double(SignalEvents)*100.) << " \%" << std::endl; output_file << std::endl << "MEC Fractional Contribution =
        // " << int(double(MECSignalEvents) / double(SignalEvents)*100.)
        // << " \%" << std::endl; output_file << std::endl << "RES Fractional Contribution = " <<
        // int(double(RESSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl; output_file << std::endl <<
        // "DIS Fractional Contribution = " << int(double(DISSignalEvents) / double(SignalEvents)*100.) << " \%" <<
        // std::endl;
        output_file << std::endl << "1232 RES plus Contribution = " << RESSignalEventsplus << std::endl;
        output_file << std::endl << "DIS plus Contribution = " << DISSignalEventsplus << std::endl;
        output_file << std::endl << "Other plus Contribution = " << OtherSignalEventsplus << std::endl;
        output_file << std::endl << "1232 RES zero Contribution = " << RESSignalEventszero << std::endl;
        output_file << std::endl << "DIS zero Contribution = " << DISSignalEventszero << std::endl;
        output_file << std::endl << "Other zero Contribution = " << OtherSignalEventszero << std::endl;
        output_file << std::endl << "1232 RES minus Contribution = " << RESSignalEventsminus << std::endl;
        output_file << std::endl << "DIS minus Contribution = " << DISSignalEventsminus << std::endl;
        output_file << std::endl << "Other minus Contribution = " << OtherSignalEventsminus << std::endl;
        // output_file << std::endl << "Total Number of Pi0 = " << totalpi0 << std::endl;
        // output_file << std::endl << "Signal Events Removed Because of Pi0 = " << eventremoved << std::endl;
        output_file
            << std::endl
            << "-----------------------------------------------------------------------------------------------------"
            << std::endl;
    }

    ofstream output_file2("./TrueFalseMap.txt");

    output_file2 << std::endl << "min DIS; Truth: yes; Reco: no  = " << tyrnDISm << std::endl;
    output_file2 << std::endl << "min DIS; Truth: yes; Reco: yes = " << tyryDISm << std::endl;
    output_file2 << std::endl << "min DIS; Truth: no; Reco: no   = " << tnrnDISm << std::endl;
    output_file2 << std::endl << "min DIS; Truth: no; Reco: yes  = " << tnryDISm << std::endl;
    output_file2 << std::endl << "min RES; Truth: yes; Reco: no  = " << tyrnRESm << std::endl;
    output_file2 << std::endl << "min RES; Truth: yes; Reco: yes = " << tyryRESm << std::endl;
    output_file2 << std::endl << "min RES; Truth: no; Reco: no   = " << tnrnRESm << std::endl;
    output_file2 << std::endl << "min RES; Truth: no; Reco: yes  = " << tnryRESm << std::endl;
    output_file2 << std::endl << "min OTH; Truth: yes; Reco: no  = " << tyrnOTHm << std::endl;
    output_file2 << std::endl << "min OTH; Truth: yes; Reco: yes = " << tyryOTHm << std::endl;
    output_file2 << std::endl << "min OTH; Truth: no; Reco: no   = " << tnrnOTHm << std::endl;
    output_file2 << std::endl << "min OTH; Truth: no; Reco: yes  = " << tnryOTHm << std::endl;
    output_file2 << std::endl << "---------------------------------------------" << std::endl;
    output_file2 << std::endl << "plus DIS; Truth: yes; Reco: no  = " << tyrnDISp << std::endl;
    output_file2 << std::endl << "plus DIS; Truth: yes; Reco: yes = " << tyryDISp << std::endl;
    output_file2 << std::endl << "plus DIS; Truth: no; Reco: no   = " << tnrnDISp << std::endl;
    output_file2 << std::endl << "plus DIS; Truth: no; Reco: yes  = " << tnryDISp << std::endl;
    output_file2 << std::endl << "plus RES; Truth: yes; Reco: no  = " << tyrnRESp << std::endl;
    output_file2 << std::endl << "plus RES; Truth: yes; Reco: yes = " << tyryRESp << std::endl;
    output_file2 << std::endl << "plus RES; Truth: no; Reco: no   = " << tnrnRESp << std::endl;
    output_file2 << std::endl << "plus RES; Truth: no; Reco: yes  = " << tnryRESp << std::endl;
    output_file2 << std::endl << "plus OTH; Truth: yes; Reco: no  = " << tyrnOTHp << std::endl;
    output_file2 << std::endl << "plus OTH; Truth: yes; Reco: yes = " << tyryOTHp << std::endl;
    output_file2 << std::endl << "plus OTH; Truth: no; Reco: no   = " << tnrnOTHp << std::endl;
    output_file2 << std::endl << "plus OTH; Truth: no; Reco: yes  = " << tnryOTHp << std::endl;

    std::cout << "File " << FileName << " created" << std::endl << std::endl;

} // End of program

// End Loop function

// -------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------

double genie_analysis::acceptance_c(double p, double cost, double phi, int particle_id, TFile *file_acceptance,
                                    bool ApplyAccWeights) {

    if (ApplyAccWeights) {

        // Redefinition of the phi angle
        // because the acceptance maps are defined between (-30,330)

        // Check that phi is between (0,360)

        // int redef = -30;
        int redef = 0;

        TH3D *acc;
        TH3D *gen;

        acc = (TH3D *)file_acceptance->Get("Accepted Particles");
        gen = (TH3D *)file_acceptance->Get("Generated Particles");

        // map 330 till 360 to [-30:0] for the acceptance map histogram
        if (phi > (2 * TMath::Pi() - TMath::Pi() / 6.)) {
            phi -= 2 * TMath::Pi();
        }
        // Find number of generated events

        double pbin_gen = gen->GetXaxis()->FindBin(p);
        double tbin_gen = gen->GetYaxis()->FindBin(cost);
        double phibin_gen = gen->GetZaxis()->FindBin(phi * 180 / TMath::Pi() + redef);
        double num_gen = gen->GetBinContent(pbin_gen, tbin_gen, phibin_gen);

        // Find number of accepted events

        double pbin_acc = acc->GetXaxis()->FindBin(p);
        double tbin_acc = acc->GetYaxis()->FindBin(cost);
        double phibin_acc = acc->GetZaxis()->FindBin(phi * 180 / TMath::Pi() + redef);
        double num_acc = acc->GetBinContent(pbin_acc, tbin_acc, phibin_acc);

        double acc_ratio = num_acc / num_gen;
        double acc_err = sqrt(acc_ratio * (1 - acc_ratio)) / num_gen;

        return acc_ratio;

    }

    else {
        return 1.;
    }
}
