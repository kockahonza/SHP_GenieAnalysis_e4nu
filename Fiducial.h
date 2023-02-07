#ifndef FIDUCIAL_H
#define FIDUCIAL_H

#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TVectorT.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include "Constants.h"

struct Fiducial {

    void InitPiMinusFit(std::string beam_en);
    void InitEClimits();
    void SetConstants(int in_TorusCurrent, std::string in_target_name, std::map<std::string, double> in_en_beam);
    void SetFiducialCutParameters(std::string beam_en);
    Bool_t GetEPhiLimits(std::string beam_en, Float_t momentum, Float_t theta, Int_t sector, Float_t *EPhiMin,
                         Float_t *EPhiMax);
    Bool_t EFiducialCut(std::string beam_en, TVector3 momentum);
    Bool_t PFiducialCut(std::string beam_en, TVector3 momentum);
    Bool_t PiplFiducialCut(std::string beam_en, TVector3 momentum, Float_t *philow, Float_t *phiup);
    Bool_t PimiFiducialCut(std::string beam_en, TVector3 momentum, Float_t *pimi_philow, Float_t *pimi_phiup);
    bool Phot_fid(TVector3 V3_phot);
    bool Pi_phot_fid_united(std::string beam_en, TVector3 V3_pi_phot, int q_pi_phot);

    // --------------------------------------------------------------------------

    // apapadop // Nov 11 2020 // Narrow band 30 deg in phi and either accepting ALL theta or theta_pos > 12 deg (piplus
    // & protons) and theta_pi- > 30

    double GetPhi(TVector3 momentum);
    double GetTheta(TVector3 momentum);

    // --------------------------------------------------------------------------

    // apapadop // Nov 23 2020 // Narrow band 30 deg in phi and either accepting ALL theta or theta_pos > 12 deg (piplus
    // & protons) and theta_pi- > 30

    Bool_t PFiducialCutExtra(std::string beam_en, TVector3 momentum);
    Bool_t PiplFiducialCutExtra(std::string beam_en, TVector3 momentum);
    Bool_t PimiFiducialCutExtra(std::string beam_en, TVector3 momentum);
    bool Phot_fidExtra(TVector3 V3_phot);
    bool Pi_phot_fid_unitedExtra(std::string beam_en, TVector3 V3_pi_phot, int q_pi_phot);

    // --------------------------------------------------------------------------

    int fTorusCurrent;
    std::string target_name;
    std::map<std::string, double> en_beam;

    TF1 *myPiMinusFit;

    TF1 *up_lim1_ec, *up_lim2_ec, *up_lim3_ec, *up_lim4_ec, *up_lim5_ec, *up_lim6_ec, *low_lim1_ec, *low_lim2_ec,
        *low_lim3_ec, *low_lim4_ec, *low_lim5_ec, *low_lim6_ec;
    TF1 *leftside_lim1_ec, *leftside_lim2_ec, *leftside_lim3_ec, *leftside_lim4_ec, *leftside_lim5_ec,
        *leftside_lim6_ec, *rightside_lim1_ec, *rightside_lim2_ec, *rightside_lim3_ec, *rightside_lim4_ec,
        *rightside_lim5_ec, *rightside_lim6_ec;

    // 4.4 GeV parameters
    // e- parameters
    Float_t fgPar_4Gev_2250_Efid_t0_p[6][2]; // 4GeV e- fiducial cut parameters
    Float_t fgPar_4Gev_2250_Efid_t1_p[6][6];
    Float_t fgPar_4Gev_2250_Efid_b_p[6][2][6];
    Float_t fgPar_4Gev_2250_Efid_a_p[6][2][6];
    Float_t fgPar_4Gev_2250_Efid_Theta_S5_extra[8][4];
    Float_t fgPar_4Gev_2250_Efid_Theta_S4_extra[2][4];
    Float_t fgPar_4Gev_2250_Efid_Theta_S3_extra[4][4];
    Float_t fgPar_4Gev_2250_Efid_Theta_S5[8][8];
    Float_t fgPar_4Gev_2250_Efid_Theta_S4[2][8];
    Float_t fgPar_4Gev_2250_Efid_Theta_S3[4][8];
    // proton parameter
    Float_t fgPar_4Gev_2250_Pfid_ScpdS2[2][6];
    Float_t fgPar_4Gev_2250_Pfid_ScpdS3[8][6];
    Float_t fgPar_4Gev_2250_Pfid_ScpdS4[4][6];
    Float_t fgPar_4Gev_2250_Pfid_ScpdS5[8][6];
    Float_t fgPar_4Gev_2250_Pfid_ScpdS2_extra[2][4];
    Float_t fgPar_4Gev_2250_Pfid_ScpdS3_extra[8][4];
    Float_t fgPar_4Gev_2250_Pfid_ScpdS4_extra[4][4];
    Float_t fgPar_4Gev_2250_Pfid_ScpdS5_extra[8][4];
    Float_t fgPar_4Gev_2250_Pfidft1l[6][6];
    Float_t fgPar_4Gev_2250_Pfidft1r[6][6];
    Float_t fgPar_4Gev_2250_Pfidft2l[6][6];
    Float_t fgPar_4Gev_2250_Pfidft2r[6][6];
    Float_t fgPar_4Gev_2250_Pfidbt1l[6][6];
    Float_t fgPar_4Gev_2250_Pfidbt1r[6][6];
    Float_t fgPar_4Gev_2250_Pfidbt2l[6][6];
    Float_t fgPar_4Gev_2250_Pfidbt2r[6][6];
    Float_t fgPar_4Gev_2250_Pfidbl[6][6];
    Float_t fgPar_4Gev_2250_Pfidbr[6][6];

    // 1.1 GeV parameters
    // electron parameters
    //  Parameters for 1.1 GeV electron momentum corrections
    //  Steven McLauchlan 22/6/00
    //  750A Torus Current
    const Double_t fgPar_1Gev_750[6][6][4] = {
        {{0.898182469399665, 0.009415220076383, -0.000321553327295, 0.000003360299566},
         {0.001433232973487, -0.000338746882960, 0.000014468017093, -0.000000163991932},
         {0.000244813389886, -0.000022627173189, 0.000000607066532, -0.000000005038009},
         {0.000037945227448, -0.000003357138839, 0.000000094832472, -0.000000000888288},
         {0.000003459429765, -0.000000272550830, 0.000000007449915, -0.000000000068594},
         {0.000000141526247, -0.000000011239286, 0.000000000308219, -0.000000000002808}},
        {{0.853575021628482, 0.013322226005995, -0.000475465884494, 0.000005233361275},
         {-0.000898871223165, -0.000025733753966, 0.000002903086982, -0.000000036882883},
         {0.000255677579516, -0.000029557989606, 0.000001033241381, -0.000000011429397},
         {-0.000010492173422, 0.000000857171206, -0.000000024089042, 0.000000000201169},
         {0.000002408872495, -0.000000203171222, 0.000000005773277, -0.000000000053976},
         {0.000000202338739, -0.000000018563944, 0.000000000564077, -0.000000000005575}},
        {{0.911770651472992, 0.007889821170976, -0.000243754875987, 0.000002347995903},
         {-0.003133249707879, 0.000372790714701, -0.000011783918257, 0.000000110900393},
         {0.000221894300060, -0.000021371378287, 0.000000679044607, -0.000000007081231},
         {-0.000010543974814, 0.000001049148338, -0.000000031334725, 0.000000000302148},
         {0.000000070833462, -0.000000007103670, 0.000000000202768, -0.000000000001513},
         {0.000000043515251, -0.000000004389163, 0.000000000136925, -0.000000000001374}},
        {{0.900135854789733, 0.008977098989975, -0.000305003333843, 0.000003336049834},
         {-0.001117976109835, 0.000159882041378, -0.000005825403903, 0.000000060340853},
         {0.000039567581002, -0.000002981619404, 0.000000051083716, -0.000000000166212},
         {-0.000014748831605, 0.000001520252189, -0.000000048731285, 0.000000000501377},
         {0.000001631829839, -0.000000145173506, 0.000000004318539, -0.000000000042152},
         {-0.000000072869112, 0.000000006405672, -0.000000000189076, 0.000000000001833}},
        {{0.907900355079541, 0.008068953843190, -0.000262565305558, 0.000002728907009},
         {-0.001364036882089, 0.000126141191203, -0.000004431528252, 0.000000050451361},
         {-0.000229998089014, 0.000024863166816, -0.000000876826425, 0.000000009694051},
         {0.000028442568480, -0.000002742455419, 0.000000085045282, -0.000000000861494},
         {0.000003017482254, -0.000000268590765, 0.000000007948062, -0.000000000076899},
         {-0.000000226776950, 0.000000020641958, -0.000000000610890, 0.000000000005897}},
        {{0.923053551059967, 0.006404299573232, -0.000206431010254, 0.000002238151922},
         {0.000309738865566, -0.000035186041814, 0.000001580561451, -0.000000032799694},
         {0.000194357506809, -0.000019027756245, 0.000000553587088, -0.000000005329462},
         {0.000038419512443, -0.000003258534473, 0.000000089783112, -0.000000000807361},
         {0.000002417465412, -0.000000178513924, 0.000000004525706, -0.000000000038157},
         {-0.000000012031345, 0.000000002470352, -0.000000000096234, 0.000000000001073}}};

    // 1500A Torus Current
    const Double_t fgPar_1Gev_1500[6][6][4] = {
        {{0.927522250807456, 0.006838093705586, -0.000216634914610, 0.000002126491398},
         {-0.001065424181991, 0.000012118819419, 0.000001252854372, -0.000000020379573},
         {0.000287399937244, -0.000028481681938, 0.000000844998746, -0.000000007970496},
         {-0.000005701504251, 0.000000474052216, -0.000000014666966, 0.000000000148530},
         {0.000003021883145, -0.000000224346212, 0.000000005678879, -0.000000000047999},
         {0.000000229918096, -0.000000018427429, 0.000000000494174, -0.000000000004374}},
        {{0.914601170548483, 0.007251908801882, -0.000230573292148, 0.000002287085867},
         {-0.002089055438004, 0.000154186060907, -0.000003918714172, 0.000000035014463},
         {0.000790794703066, -0.000076159473051, 0.000002316285702, -0.000000022628480},
         {-0.000000708021955, -0.000000347377892, 0.000000018081857, -0.000000000233126},
         {-0.000001213272731, 0.000000135645188, -0.000000004427876, 0.000000000045009},
         {0.000000063276196, -0.000000004082613, 0.000000000088234, -0.000000000000617}},
        {{0.935631222686180, 0.005585404612784, -0.000156376049159, 0.000001379236192},
         {-0.004877851226350, 0.000504997794802, -0.000015818237893, 0.000000153744479},
         {0.000490900090576, -0.000044390335634, 0.000001288625471, -0.000000012208528},
         {0.000004566247909, -0.000000541218941, 0.000000021937466, -0.000000000265330},
         {-0.000001194006097, 0.000000111576754, -0.000000003288688, 0.000000000031335},
         {0.000000041007057, -0.000000003568147, 0.000000000093856, -0.000000000000784}},
        {{0.931040799286932, 0.006371487609753, -0.000201715932102, 0.000002050468907},
         {-0.000489244108847, 0.000045653853665, -0.000001203796876, 0.000000007136121},
         {0.000461822839512, -0.000044509197550, 0.000001362037039, -0.000000013428053},
         {0.000012999567648, -0.000000870321095, 0.000000018181961, -0.000000000115335},
         {0.000000642661575, -0.000000038084641, 0.000000000689551, -0.000000000003303},
         {-0.000000143596871, 0.000000011525046, -0.000000000303449, 0.000000000002618}},
        {{0.928191885851796, 0.006645080881866, -0.000204271556022, 0.000002008077759},
         {-0.002053902117661, 0.000179241092848, -0.000005343840064, 0.000000051303231},
         {0.000371986998586, -0.000033350741883, 0.000000943897251, -0.000000008689088},
         {0.000035581588774, -0.000003077572346, 0.000000084399090, -0.000000000751597},
         {0.000000945636957, -0.000000071310579, 0.000000001864456, -0.000000000016189},
         {-0.000000220399980, 0.000000018723362, -0.000000000512754, 0.000000000004567}},
        {{0.933503294604856, 0.006054244382619, -0.000184762370931, 0.000001846020607},
         {-0.000448143402542, 0.000063933289754, -0.000002321140845, 0.000000018965335},
         {0.000519400090511, -0.000050090643825, 0.000001500055721, -0.000000014413091},
         {-0.000002207853107, 0.000000192320612, -0.000000005418561, 0.000000000044185},
         {0.000001217729185, -0.000000066145585, 0.000000001150279, -0.000000000005763},
         {0.000000138636968, -0.000000010625772, 0.000000000273367, -0.000000000002318}}};

    double fgPar_1gev_750_Efid[6][5][6];
    double fgPar_1gev_750_Efid_Theta_S3[4][8];
    double fgPar_1gev_750_Efid_Theta_S4[2][8];
    double fgPar_1gev_750_Efid_Theta_S5[8][8];
    double fgPar_1gev_1500_Efid[6][5][6];
    double fgPar_1gev_1500_Efid_Theta_S3[4][8];
    double fgPar_1gev_1500_Efid_Theta_S4[2][8];
    double fgPar_1gev_1500_Efid_Theta_S5[8][8];

    // proton parameters
    double fgPar_1gev_750_Pfid[6][5][6];
    double fgPar_1gev_750_Pfid_ScpdS2[2][6];
    double fgPar_1gev_750_Pfid_ScpdS3[8][6];
    double fgPar_1gev_750_Pfid_ScpdS4[4][6];
    double fgPar_1gev_750_Pfid_ScpdS5[8][6];
    double fgPar_1gev_1500_Pfid[6][5][6];
    double fgPar_1gev_1500_Pfid_ScpdS2[2][6];
    double fgPar_1gev_1500_Pfid_ScpdS3[8][6];
    double fgPar_1gev_1500_Pfid_ScpdS4[4][6];
    double fgPar_1gev_1500_Pfid_ScpdS5[8][6];

    // Pion Plus parameter?
    double fgPar_1gev_750_Piplfid[6][5][6];
    double fgPar_1gev_750_Piplfid_ScpdS2[2][6];
    double fgPar_1gev_750_Piplfid_ScpdS3[8][6];
    double fgPar_1gev_750_Piplfid_ScpdS4[4][6];
    double fgPar_1gev_750_Piplfid_ScpdS5[8][6];
    double fgPar_1gev_1500_Piplfid[6][5][6];
    double fgPar_1gev_1500_Piplfid_ScpdS2[2][6];
    double fgPar_1gev_1500_Piplfid_ScpdS3[8][6];
    double fgPar_1gev_1500_Piplfid_ScpdS4[4][6];
    double fgPar_1gev_1500_Piplfid_ScpdS5[8][6];

    // Pion Minus parameter?
    double fgPar_1gev_750_Pimfid[6][5][6];
    double fgPar_1gev_750_Pimfid_Theta_S3[4][8];
    double fgPar_1gev_750_Pimfid_Theta_S4[2][8];
    double fgPar_1gev_750_Pimfid_Theta_S5[8][8];
    double fgPar_1gev_750_Pimfid_Theta_S3_extra[4][4];
    double fgPar_1gev_750_Pimfid_Theta_S4_extra[2][4];
    double fgPar_1gev_750_Pimfid_Theta_S5_extra[8][4];
    double fgPar_1gev_1500_Pimfid[6][5][6];
    double fgPar_1gev_1500_Pimfid_Theta_S3[4][8];
    double fgPar_1gev_1500_Pimfid_Theta_S4[2][8];
    double fgPar_1gev_1500_Pimfid_Theta_S5[8][8];
    double fgPar_1gev_1500_Pimfid_Theta_S3_extra[4][4];
    double fgPar_1gev_1500_Pimfid_Theta_S4_extra[2][4];
    double fgPar_1gev_1500_Pimfid_Theta_S5_extra[8][4];

    // 1GeV e- and pi- fiducial parameters (theta gaps)
    const Double_t fid_1gev_750_efid_S5[3][2][4] = {
        {{111.088, 7.4495, -0.810668, 0.0314985}, {119.088, 7.4495, -0.810668, 0.0314985}},
        {{19.2082, 1.65144, 0.661954, -0.165823}, {18.2403, 8.87948, -5.21333, 1.28453}},
        {{34.099, -16.9152, 13.2937, -2.83509}, {23.2426, 10.3194, -5.30272, 1.35695}}};
    // sector 4 has 2 gaps
    const Double_t fid_1gev_750_efid_S4[2][2][4] = {
        {{30.5039, 5.98362, -12.825, 7.04185}, {25.4698, 13.3482, -7.32756, 1.61262}},
        {{99.758, 1.35143, -0.115901, 0.0033244}, {111.575, 0.00204503, 0.122899, -0.00687251}}};
    // sector 3 has 4 gaps
    const Double_t fid_1gev_750_efid_S3[4][2][4] = {
        {{17.1438, 5.82672, -1.33914, 0.370319}, {8.71859, 23.4486, -10.8172, 1.97916}},
        {{34.8286, 4.57815, -0.259835, 0.00424679}, {39.4437, 4.52068, -0.78184, 0.0893456}},
        {{109.431, 5.24502, -0.439109, 0.0104981}, {127.462, 1.91687, -0.174431, 0.0118813}},
        {{97.15, 1.94998, 0, 0}, {109.858, 1.19822, 0, 0}}};
    // Original sector 3 parameters, not defined to low momentum
    /*   const Double_t fid_1gev_750_efid_S3[4][2][4] = {{{17.1438 , 5.82672 , -1.33914 , 0.370319}, */
    /*               {8.71859 , 23.4486 , -10.8172 , 1.97916}}, */
    /*               {{37.5065 , 1.33304 , 0.319901 , 0.181076}, */
    /*               {38.9433 , 4.26507 , 0.260988 , -0.296786}}, */
    /*               {{109.431 , 5.24502 , -0.439109 , 0.0104981}, */
    /*               {127.462,  1.91687,  -0.174431,  0.0118813}}, */
    /*               {{97.15,  1.94998 , 0 , 0},{109.858, 1.19822 ,0, 0}}}; */
    // sector 2 has 1 gap1
    const Double_t fid_1gev_750_efid_S2[2][4] = {{19.6927, 8.30968, -2.7451, 0.457127}, {24.816, 3.00121, 0, 0}};

    // pi minus theta gap parameters
    const Double_t fid_1gev_750_pimifid_S5[2][2][4] = {
        {{111.088, 7.4495, -0.810668, 0.0314985}, {119.088, 7.4495, -0.810668, 0.0314985}},
        //{{19.2082, 1.65144, 0.661954, -0.165823}, {18.2403, 8.87948 ,-5.21333, 1.28453}},
        {{82.588, 7.4495, -0.810668, 0.0314985}, {90.588, 7.4495, -0.810668, 0.0314985}}};
    // sector 4 has 2 gaps
    const Double_t fid_1gev_750_pimifid_S4[2][2][4] = {
        {{35.8238, 6.86212, -0.686454, 0.0249186}, {38.5397, 8.72133, -1.18367, 0.0582613}},
        {{99.758, 1.35143, -0.115901, 0.0033244}, {111.575, 0.00204503, 0.122899, -0.00687251}}};
    // sector 3 has 4 gaps
    const Double_t fid_1gev_750_pimifid_S3[3][2][4] =
        { //{{17.1438 , 5.82672 , -1.33914 , 0.370319},  {8.71859 , 23.4486 , -10.8172 , 1.97916}},
            {{34.8286, 4.57815, -0.259835, 0.00424679}, {39.4437, 4.52068, -0.78184, 0.0893456}},
            {{109.431, 5.24502, -0.439109, 0.0104981}, {127.462, 1.91687, -0.174431, 0.0118813}},
            {{97.15, 1.94998, 0, 0}, {109.858, 1.19822, 0, 0}}};
    // sector 2 has 1 gap1
    const Double_t fid_1gev_750_pimifid_S2[2][4] = {{19.6927, 8.30968, -2.7451, 0.457127}, {24.816, 3.00121, 0, 0}};

    // sector 1 has 1 gap
    const Double_t fid_1gev_750_pimifid_S1[2][4] = {{113.259, 10.763, -2.45037, 0.201482},
                                                    {135.89, 1.31017, -0.27533, 0.0355386}};

    // 1GeV p and pi+ fiducial parameters (theta gaps)
    // sector 1 has 2 gaps
    const Double_t fid_1gev_750_pfid_S1[2][2][4] = {
        {{120.448, 1.10613, -0.794723, 0.0656712}, {147.559, -11.4625, 2.2709, -0.158945}},
        {{67.1972, -25.2405, 7.16785, -0.748601}, {50.2896, -0.809075, -1.79123, 0.334088}}};
    // sector 2 has 1 gap1
    const Double_t fid_1gev_750_pfid_S2[2][4] = {{47.7655, -4.02984, 0, 0}, {49.2938, -2.62524, 0, 0}};
    // sector 3 has 4 gaps, first argument is the # of gaps, second corresponds to upper and lower limits of gap, the
    // last one is the number of parameters
    const Double_t fid_1gev_750_pfid_S3[4][2][4] = {
        {{118.212, -2.17485, 0.0211904, 0.00847772}, {137.918, -11.0922, 2.63701, -0.227367}},
        {{99.4546, -1.72803, 0, 0}, {108.326, -1.95039, 0.0688691, -0.00303378}},
        {{48.1261, -4.15847, 0, 0}, {49.2546, -2.54427, 0, 0}},
        {{35.373, -0.974947, -1.47085, 0.164571}, {43.4268, -6.57691, 0.940029, -0.134398}}};
    // sector 4 has 2 gaps
    const Double_t fid_1gev_750_pfid_S4[2][2][4] = {
        {{104.013, -0.902977, 0, 0}, {110.284, -0.663729, 0, 0}},
        {{47.9585, -4.17542, 0, 0}, {56.2454, -12.0738, 4.15965, -0.586082}}};
    // sector 5 has 4 gaps, first argument is the # of gaps, second corresponds to upper and lower limits of gap, the
    // last one is the number of parameters
    const Double_t fid_1gev_750_pfid_S5[4][2][4] = {
        {{122.754, -1.38182, 0, 0}, {150.113, -15.231, 3.59444, -0.294359}},
        {{46.5518, -3.76049, 0, 0}, {54.3328, -8.0828, 1.8003, -0.212943}},
        {{24.943, -4.91541, 0, 0}, {4.61651, 43.3407, -32.1895, 6.87298}},
        {{97.1144, -1.62607, 0, 0}, {121.22, -12.8209, 1.93743, -0.108376}}};
    // sector 1 has 2 gaps
    const Double_t fid_1gev_750_pfid_S6[2][2][4] = {
        {{64.3241, -25.288, 8.33577, -1.01602}, {37.2644, 10.6393, -4.63834, 0.507171}},
        {{45.6131, 7.2492, -1.9949, 0.200048}, {57.2786, 1.79723, -0.506705, 0.117659}}};

    // 2.2 GeV parameters
    const Float_t fgPar_2GeV_2250_Efid[6][6][9] = {
        {{62.2935, -92.5133, 87.0360, -38.4696, 6.3177, 0, 0, 0, 0},
         {78.5134, -58.5975, 3.30928, 77.4749, -64.3984, 14.4860, 0, 0, 0},
         {-140.845, 1381.30, -4499.99, 7557.27, -7140.27, 3828.75, -1086.21, 126.468, 0},
         {497.951, -1846.42, 2759.58, -1634.71, 345.006, 0, 0, 0, 0},
         {9.40986, 180.752, -646.771, 1055.14, -909.094, 424.435, -99.8368, 9.02086, 0},
         {288.485, -1016.03, 1463.72, -859.231, 185.976, 0, 0, 0, 0}},
        {{61.1474, -88.768, 82.6446, -36.2780, 5.92310, 0, 0, 0, 0},
         {78.5134, -58.5975, 3.30928, 77.4749, -64.3984, 14.4860, 0, 0, 0},
         {21.3087, 138.975, -672.710, 1324.20, -1326.12, 714.866, -197.531, 21.9144, 0},
         {375.091, -1411.50, 2082.58, -1192.17, 239.685, 0, 0, 0, 0},
         {-121.816, 1182.59, -3800.98, 6319.82, -5937.33, 3179.37, -903.954, 105.764, 0},
         {-4781.96, 43165.9, -159567, 318502, -376469, 271207, -116893, 27698.9, -2775.61}},
        {{61.1474, -88.7680, 82.6446, -36.2780, 5.92310, 0, 0, 0, 0},
         {73.7620, -34.6321, -41.8796, 117.543, -81.2043, 17.1718, 0, 0, 0},
         {157.046, -765.472, 1735.21, -2053.86, 1371.34, -515.214, 101.081, -8.07402, 0},
         {-608.740, 4827.18, -13239.6, 17742.4, -12420.0, 4369.11, -607.877, 0, 0},
         {-274.278, 2380.63, -7560.19, 12582.3, -11924.5, 6464.66, -1863.44, 221.134, 0},
         {-1240.72, 8096.04, -19407.0, 23942.9, -16052.3, 5559.32, -776.123, 0, 0}},
        {{61.1474, -88.7680, 82.6446, -36.2780, 5.92310, 0, 0, 0, 0},
         {78.5134, -58.5975, 3.30928, 77.4749, -64.3984, 14.4860, 0, 0, 0},
         {-71.2528, 879.668, -3027.37, 5226.61, -4999.19, 2689.35, -761.206, 88.1242, 0},
         {-1269.89, 9486.25, -26103.8, 35581.2, -25373.0, 9062.87, -1277.60, 0, 0},
         {-186.640, 1811.85, -6032.01, 10283.3, -9808.11, 5285.35, -1501.87, 174.799, 0},
         {-530.826, 4643.56, -13864.2, 20580.2, -15898.0, 6106.69, -916.365, 0, 0}},
        {{61.6665, -90.4268, 84.5606, -37.2240, 6.09207, 0, 0, 0, 0},
         {78.5134, -58.5975, 3.30928, 77.4749, -64.3984, 14.4860, 0, 0, 0},
         {-1.53910, 216.936, -701.057, 1167.26, -1111.92, 615.364, -183.854, 22.8595, 0},
         {-19.7415, 454.317, -1250.51, 1512.52, -762.408, 137.695, 0, 0, 0},
         {-55.9612, 657.449, -2049.73, 3295.30, -2995.85, 1553.68, -427.764, 48.4324, 0},
         {-522.682, 3356.77, -7535.50, 8756.49, -5518.61, 1795.60, -235.144, 0, 0}},
        {{61.1474, -88.7680, 82.6446, -36.2780, 5.92310, 0, 0, 0, 0},
         {73.7620, -34.6321, -41.8796, 117.543, -81.2043, 17.1718, 0, 0, 0},
         {-82.0368, 883.261, -2828.84, 4621.53, -4223.56, 2185.52, -598.218, 67.2908, 0},
         {608.323, -2743.56, 4942.01, -4045.58, 1558.07, -226.240, 0, 0, 0},
         {4.07203, 138.882, -321.983, 282.702, -12.9566, -129.159, 74.5884, -12.9994, 0},
         {-866.737, 5984.13, -15129.6, 19134.6, -12757.7, 4276.79, -566.056, 0, 0}}};

    const Float_t fgPar_2GeV_2250_EfidTheta_S3[4][8] = {
        {74.4893, -158.720, 251.241, -200.000, 52.2984, 25.4188, -18.8692, 3.27217},
        {90.8413, -226.800, 358.487, -259.260, 30.9359, 68.7248, -38.9760, 6.47933},
        {117.102, -429.455, 1208.29, -1922.72, 1791.40, -965.135, 277.459, -32.8536},
        {55.0676, 91.0959, -444.252, 791.284, -717.492, 350.325, -87.3235, 8.68087}};

    const Float_t fgPar_2GeV_2250_EfidTheta_S4[2][8] = {
        {77.7940, -192.492, 361.852, -394.127, 246.499, -84.6133, 13.9182, -0.713846},
        {41.2902, 110.603, -586.690, 1130.70, -1137.27, 633.345, -185.038, 22.1482}};

    const Float_t fgPar_2GeV_2250_EfidTheta_S5[8][8] = {
        {-12998.3, 57694.0, -109085, 114102, -71303.6, 26616.9, -5494.45, 483.756},
        {17842.9, -74659.9, 133869, -133170, 79380.0, -28352.3, 5617.88, -476.314},
        {65.5364, -99.8689, 88.5645, 24.8299, -121.327, 102.818, -37.8275, 5.28492},
        {66.4049, -76.4096, -20.8674, 230.072, -318.905, 206.721, -66.3286, 8.48753},
        {100.262, -358.882, 957.267, -1495.42, 1396.73, -765.881, 226.791, -27.9341},
        {50.4447, 48.3032, -315.976, 580.141, -525.583, 252.075, -59.9294, 5.34805},
        {78.5845, -155.728, 320.528, -420.296, 341.899, -164.626, 42.5274, -4.50224},
        {95.9430, -221.787, 391.495, -350.033, 131.391, 13.2965, -24.0460, 4.92253}};

    const Float_t fgPar_2GeV_2250_Pfid_For[6][4][7] = {
        {{60.2165, -189.720, 446.990, -523.122, 320.721, -97.8518, 11.5258},
         {-1457.16, 13814.2, -43182.7, 66646.0, -54355.1, 22423.5, -3683.76},
         {17.1086, 54.2974, -103.464, 111.325, -70.7673, 27.2551, -5.02858},
         {-2547.86, 22143.1, -66326.6, 101105.0, -82187.8, 33959.7, -5607.59}},
        {{65.7242, -246.922, 759.745, -1198.32, 1007.05, -428.060, 72.2644},
         {3384.16, -19353.1, 54083.5, -79843.4, 63870.2, -26079.2, 4250.29},
         {85.2489, -441.821, 1327.52, -1978.53, 1567.84, -633.530, 102.928},
         {411.998, -533.572, 599.925, 2099.52, -5061.48, 3701.58, -891.843}},
        {{110.022, -558.044, 1512.96, -2098.53, 1579.55, -613.478, 96.3279},
         {3937.29, -23745.1, 59651.0, -76988.6, 54276.0, -19900.2, 2974.95},
         {35.8488, -46.9595, 107.492, -93.9141, 10.5845, 26.1910, -9.89460},
         {-326.838, 4634.99, -11155.2, 11811.4, -5405.80, 554.030, 175.526}},
        {{38.9338, -62.8663, 118.218, -56.6953, -40.5083, 46.1782, -11.5822},
         {1864.83, -11735.6, 34175.4, -48928.5, 37315.8, -14496.1, 2254.05},
         {23.6892, 9.69854, 94.4521, -270.119, 288.132, -140.031, 25.9272},
         {-261.086, 4863.13, -11760.4, 13791.1, -8983.19, 3136.52, -457.183}},
        {{-11.0252, 348.901, -1172.63, 1980.73, -1759.08, 786.043, -139.299},
         {-2231.41, 23477.1, -78229.3, 129238.0, -111761.0, 48561.4, -8370.65},
         {104.415, -548.464, 1506.70, -2064.10, 1507.55, -561.677, 83.9247},
         {1402.87, -9008.78, 25660.0, -37543.3, 29860.8, -12238.4, 2019.03}},
        {{20.4577, 66.1373, -205.218, 372.864, -366.625, 177.596, -33.1168},
         {2059.77, -14468.3, 46492.9, -72168.2, 58275.9, -23615.8, 3800.60},
         {-18.9897, 392.519, -1234.31, 1950.24, -1623.01, 681.260, -113.806},
         {-3478.50, 32840.9, -104381.0, 167656.0, -143070.0, 61909.3, -10690.1}}};

    const Float_t fgPar_2GeV_2250_Pfid_Bak[6][4][7] = {
        {{110.007, 121.302, 97.8380, -1679.71, 4022.73, -3973.09, 1422.42},
         {69.7305, 359.843, -876.383, 649.612, 600.059, -1155.43, 472.866},
         {13.9334, -236.587, 810.783, -1614.65, 1851.97, -1125.48, 280.069},
         {10.1644, 51.7943, -527.843, 2071.12, -3480.34, 2663.52, -768.498}},
        {{161.555, -263.801, 770.924, -902.814, 503.641, -319.619, 171.147},
         {154.660, -619.711, 3444.65, -8994.29, 12253.9, -8439.82, 2321.14},
         {117.461, -1429.96, 6117.79, -13492.3, 16142.2, -9965.40, 2490.47},
         {7.77411, -17.3501, 279.462, -876.326, 1398.82, -1137.49, 365.383}},
        {{-31.1460, 1942.49, -9193.97, 21731.0, -26961.3, 16701.7, -4067.85},
         {154.660, -654.420, 3774.08, -9920.36, 13333.7, -8953.68, 2386.32},
         {63.2709, -867.859, 4000.97, -9557.57, 12215.1, -7926.91, 2052.90},
         {-28.1127, 484.636, -2665.71, 7484.94, -10740.7, 7561.79, -2076.70}},
        {{172.853, -656.312, 3768.76, -10243.0, 14600.3, -10616.3, 3095.27},
         {270.076, -1938.46, 9276.01, -21861.1, 27363.7, -17479.9, 4490.05},
         {32.2327, -432.593, 1666.57, -3491.43, 4031.58, -2406.30, 579.944},
         {-44.9153, 638.112, -2971.77, 7223.13, -9328.99, 6080.46, -1576.13}},
        {{45.7403, 875.133, -3646.85, 7848.52, -8905.36, 4914.78, -1010.91},
         {138.000, -449.485, 2806.13, -7725.44, 10777.3, -7482.95, 2056.80},
         {72.7551, -944.002, 4200.92, -9776.76, 12316.6, -7955.78, 2066.50},
         {-9.59531, 180.519, -795.797, 2124.85, -2978.29, 2040.14, -541.811}},
        {{77.5100, 494.571, -1625.99, 2397.48, -1177.99, -574.604, 530.446},
         {117.869, -56.8761, 330.252, -715.276, 807.257, -497.124, 133.989},
         {7.66164, -208.001, 996.883, -2772.33, 4100.81, -3008.90, 864.126},
         {-25.3497, 346.501, -1458.46, 3513.62, -4625.70, 3088.01, -818.696}}};

    const Float_t fgPar_2GeV_2250_Pfid_ScpdS2[2][6] = {{-28.1486, 425.124, -935.693, 1065.39, -608.526, 137.658},
                                                       {-15.2084, 345.466, -697.657, 751.738, -419.288, 95.2206}};

    const Float_t fgPar_2GeV_2250_Pfid_ScpdS3[8][6] = {{17.1490, 294.605, -640.590, 707.758, -386.730, 83.2529},
                                                       {35.9318, 204.580, -404.489, 413.240, -209.580, 41.7819},
                                                       {47.6825, 274.777, -754.725, 1117.80, -846.816, 255.607},
                                                       {44.7484, 344.543, -872.200, 1113.89, -694.736, 168.061},
                                                       {-205.978, 828.617, -1199.65, 875.482, -317.846, 45.6938},
                                                       {-240.595, 961.068, -1370.34, 977.625, -345.743, 48.3834},
                                                       {-136.104, 479.276, -593.135, 374.730, -118.350, 14.7923},
                                                       {-196.773, 700.974, -894.540, 577.460, -185.690, 23.6201}};

    const Float_t fgPar_2GeV_2250_Pfid_ScpdS4[4][6] = {{81.8115, 139.810, -445.130, 804.212, -821.194, 364.924},
                                                       {79.5053, 317.287, -1582.80, 3987.05, -4880.55, 2305.63},
                                                       {-137.480, 633.288, -954.383, 721.057, -269.140, 39.4822},
                                                       {-145.605, 697.662, -1088.74, 853.855, -330.883, 50.3421}};

    const Float_t fgPar_2GeV_2250_Pfid_ScpdS5[8][6] = {{-29.9426, 370.963, -714.697, 707.343, -348.995, 67.7647},
                                                       {-27.4173, 372.536, -693.341, 652.792, -302.559, 54.7761},
                                                       {-47.1617, 132.967, -104.776, 41.7673, -7.68238, 0.404311},
                                                       {-54.5895, 149.685, -111.590, 41.2556, -6.93943, 0.301087},
                                                       {-79.1386, 275.678, -341.972, 218.907, -69.5520, 8.66381},
                                                       {-97.5794, 352.616, -468.487, 322.829, -111.159, 15.0975},
                                                       {22.5823, -182.064, 365.317, -294.653, 108.779, -15.2712},
                                                       {-7.59521, 2.91795, 31.6773, -28.3085, 10.5943, -1.57966}};

    // 2 and 4GeV analysis use the same fiducial functions for pimi, the parameters for theta vs phi outline cuts are
    // shown below
    //
    const Double_t fid_2gev_2250_pimifid_outline[6][6][5] = {{{437.992, -6758.33, 44725, -130667, 140000},
                                                              {-108.89, 2441.34, -16273.9, 46629.1, -48613.3},
                                                              {23.4972, -22.8102, 58.051, 455.812, -1355.76},
                                                              {-1282.61, 26039.3, -176254, 497349, -504262},
                                                              {79.2767, 10.9846, -1097.42, -2505.95, 11295.1},
                                                              {-116.656, 5464.17, -39666.7, 121333, -133333}},
                                                             {{294.906, -3993.33, 25700, -74666.7, 80000},
                                                              {-1.73619, 230.33, -230.96, -2041.14, 4247.97},
                                                              {-8.35081, 294.696, -351.216, -2783.33, 6284.44},
                                                              {-654.396, 13447.8, -95677.4, 289413, -316333},
                                                              {-1087.77, 21543.4, -150270, 446839, -481346},
                                                              {-12.4219, 3796.67, -30450, 101333, -120000}},
                                                             {{169.867, -1621.67, 9791.67, -29333.3, 33333.3},
                                                              {1.31563, 228.273, -277.831, -2151.54, 4681.77},
                                                              {4.99614, 143.415, -99.2912, -1182.49, 2131.09},
                                                              {-442.257, 10171.2, -74092.5, 222973, -240250},
                                                              {286.575, -4807.23, 31987.4, -94543, 102888},
                                                              {-80.793, 4620.91, -32275.5, 94668.1, -100002}},
                                                             {{171.781, -1604.17, 9133.31, -25333.3, 26666.7},
                                                              {15.4017, 65.7045, -29.3075, -398.534, 391.892},
                                                              {-105.713, 2167.31, -13426.9, 36592.6, -36900.5},
                                                              {-117.86, 4749.21, -44679.4, 161213, -199605},
                                                              {-1228.95, 24386.5, -171732, 516890, -564141},
                                                              {314.227, -2770.83, 17208.3, -46666.6, 46666.7}},
                                                             {{227.727, -2592.5, 15608.3, -44000, 46666.7},
                                                              {-108.819, 2413.21, -16243.6, 47773.3, -51709},
                                                              {-1.34202, 206.441, -194.54, -1788.54, 3644.16},
                                                              {-2105.78, 42101.1, -299023, 908108, -1.00031e+06},
                                                              {-542.533, 11447.1, -83967.3, 260870, -292292},
                                                              {198.746, -468.327, 479.183, 5155.48, -11578.4}},
                                                             {{205.403, -1966.22, 9777.77, -22675.4, 20010.2},
                                                              {5.0732, 177.775, -195.135, -1565.26, 3219.17},
                                                              {114.9, -2022.82, 15044.2, -46493.2, 51247.4},
                                                              {214.243, -4966.33, 42750.9, -149352, 180657},
                                                              {1230.36, -24868.8, 180507, -555189, 613732},
                                                              {166.471, -85.4706, -94.7685, 1150.18, -2073.65}}};

    // 2 and 4 GeV pi- fiducial parameters (theta gaps)
    // sector 5 has 2 gaps
    const Double_t fid_2gev_2250_efid_S5[2][2][4] = {
        {{88.912, 39.7076, -11.7886, 1.219}, {102.112, 37.1564, -12.6383, 1.49529}},
        {{100.903, -0.978802, 1.99633, -0.285675}, {103.959, 2.03364, -0.0566607, 0.11236}}};

    // sector 3 has 2 gaps
    const Double_t fid_2gev_2250_efid_S3[2][2][4] = {
        {{112.584, -8.50953, 3.65541, -0.321358}, {99.7915, 11.1871, -1.96114, 0.167723}},
        {{15.4065, 35.6826, -11.4842, 1.35174}, {-8.2081, 44.1171, -9.23191, 0.696626}}};

    // sector 1 has 2 gaps
    const Double_t fid_2gev_2250_efid_S1[2][2][4] = {
        {{113.956, 0.363286, 0.198558, -0.0144791}, {124.963, -1.25147, 0.390258, -0.0194753}},
        {{172.863, -31.0937, 6.70311, -0.4403}, {135.614, -2.55523, 0.881528, -0.0632905}}};

    // sector 4 has 1 gap
    const Double_t fid_2gev_2250_efid_S4[2][4] = {{99.1183, 3.87431, -0.651923, 0.0541116},
                                                  {114.798, -1.74102, 0.670956, -0.0506822}};

    // sector 6 has 1 gap
    const Double_t fid_2gev_2250_efid_S6[2][4] = {{80.5638, -15.6315, 3.64252, -0.267854},
                                                  {51.8427, 10.8915, -2.38221, 0.177976}};

    // 2 and 4 GeV pi- fiducial parameters (the pol3 theta gaps, with a sensible variable name, redefined some electron
    // gaps here) sector 1 has 2 gaps
    const Double_t fid_2gev_2250_pimifid_S1[2][2][4] = {
        {{113.956, 0.363286, 0.198558, -0.0144791}, {124.963, -1.25147, 0.390258, -0.0194753}},
        {{124.435, 2.93323, -0.655971, 0.0667888}, {135.614, -2.55523, 0.881528, -0.0632905}}};

    // sector 2 currently has no gaps
    // const Double_t fid_2gev_2250_pimifid_S2

    // sector 3 has 3 gaps
    const Double_t fid_2gev_2250_pimifid_S3[3][2][4] = {
        {{112.584, -8.50953, 3.65541, -0.321358}, {99.7915, 11.1871, -1.96114, 0.167723}},
        {{37.6064, 7.00954, -0.489086, 0.0444082}, {41.3276, 4.96432, 0.660286, -0.0701486}},
        {{15.6704, 15.8788, -2.70209, 0.203812}, {17.941, 15.2645, -2.67317, 0.220213}}};

    // sector 4 has 2 gaps
    const Double_t fid_2gev_2250_pimifid_S4[2][2][4] = {
        {{99.1183, 3.87431, -0.651923, 0.0541116}, {114.798, -1.74102, 0.670956, -0.0506822}},
        {{17.4548, 12.0517, -1.60906, 0.108019}, {20.2729, 8.8091, -0.0491563, -0.0567015}}};

    // sector 5 has 4 gaps
    const Double_t fid_2gev_2250_pimifid_S5[4][2][4] = {
        {{99.9404, 24.7386, -5.60018, 0.431221}, {119.312, 13.4537, -2.63984, 0.190205}},
        {{97.5935, 4.44785, -0.611382, 0.0758808}, {103.94, 1.76833, 0.32453, -0.00709944}},
        {{28.009, 12.8943, -2.13127, 0.166442}, {30.6787, 13.2269, -2.2456, 0.210455}},
        {{18.7992, 8.15369, -2.39897, 0.235203}, {21.175, 6.73867, -1.91665, 0.180353}}};
    // sector 6 has 1 gap
    const Double_t fid_2gev_2250_pimifid_S6[2][4] = {{80.5638, -15.6315, 3.64252, -0.267854},
                                                     {51.8427, 10.8915, -2.38221, 0.177976}};

    // The momentum dependent theta upper limits for pi- at 4GeV and all particles at 1GeV
    Float_t pimi_thetamax2and4[5] = {148.8, 64.7164, -234.766, 140.342, -25.4952};
    Float_t pimi_thetamax1[5] = {180.316, -510.125, 2592.66, -5389.62, 3393.97};
    Float_t pipl_thetamax1[5] = {143.52, -114.506, 409.901, -461.16, 97.7215};
    Float_t el_thetamax1[5] = {105.51, -262.424, 469.016, -365.019, 102.453};
};
#endif
