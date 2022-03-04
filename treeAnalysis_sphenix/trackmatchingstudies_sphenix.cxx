#include <algorithm>
// ANCHOR debug output verbosity
int verbosityTM = 1;
bool initializeHistsTMS = false;
//************************************************************************************************************
//************************************************************************************************************
// declaration of histogramas and global variables
//************************************************************************************************************
//************************************************************************************************************
TH2F*  h_TMstudies_2D_clus[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_2D_true[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_2D_true_electron[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_2D_true_hadron[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_2D_track[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_2D_projection[_active_calo]; // [calorimeter_enum][algorithm_enum]

TH1F*  h_TMstudies_clusSpec[_active_calo]; // [calorimeter_enum][algorithm_enum]-
TH1F*  h_TMstudies_clusSpec_Eta[_active_calo][nEta+1]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_Matched[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_Matched_Eta[_active_calo][nEta+1]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_MatchedCalo[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_MatchedCalo_EDiff1[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_MatchedCalo_EDiff2[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_special[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_special_Matched[_active_calo]; // [calorimeter_enum][algorithm_enum]

const int diffCases_TMStudies = 1;
TString strCases_TMStudies[diffCases_TMStudies] = {"all"};//, "leading"};

TH2F*  h_TMstudies_2D_delta_track[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_2D_delta_projection[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_2D_delta_projection_special[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_2D_delta_TTLhits[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_2D_delta_calomatch[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_E1_E2_calomatch[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dx_vsClsEta_track[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dx_vsClsEta_projection[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dx_vsTrkEta_track[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dx_vsTrkEta_projection[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dx_vsClsE_projection[diffCases_TMStudies][_active_calo][_active_algo][nEta+1]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dy_vsClsE_projection[diffCases_TMStudies][_active_calo][_active_algo][nEta+1]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_isMatched[diffCases_TMStudies][_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_isMatched_E[diffCases_TMStudies][_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_isMatched_PDG[diffCases_TMStudies][_active_calo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_isMatched_findable[diffCases_TMStudies][_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_isMatched_findable_E[diffCases_TMStudies][_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_isMatched_findable_PDG[diffCases_TMStudies][_active_calo]; // [calorimeter_enum][algorithm_enum]

//************************************************************************************************************
//************************************************************************************************************
// ANCHOR track/projection matching function
// TODO at some point we might want E/p
//************************************************************************************************************
//************************************************************************************************************
bool trackmatchingstudies_sphenix(bool runSpecialCuts = false){


  // boolean for single event histo?
  bool b_TMstudies_2D_track_filled[_active_calo] = {false};

  //************************************************************************************************************
  //*********************************  create histograms  ******************************************************
  //************************************************************************************************************
  if (!initializeHistsTMS){
    for(int icalo=0;icalo<_active_calo;icalo++){
      if (!caloEnabled[icalo]) continue;

      Int_t caloDir           = GetCaloDirection(icalo);
      Float_t minEtaCurCalo   = partEtaCalo[minEtaBinCalo[caloDir]];
      Float_t maxEtaCurCalo   = partEtaCalo[maxEtaBinCalo[caloDir]];

    //  bool isFwd = IsForwardCalorimeter(icalo);
      float nbins2D =300;
      float min2Dhist = -2.0;
      float max2Dhist = 2.0;
      float min2Dhist2 =  -TMath::Pi();
      float max2Dhist2 = TMath::Pi();
      TString str2DhistoName =  "eta_phi";
      if(!h_TMstudies_2D_projection[icalo])
        h_TMstudies_2D_projection[icalo]        = new TH2F(Form("h_TMstudies_%s_projection_%s",str2DhistoName.Data(),str_calorimeter[icalo].Data()), "", nbins2D, min2Dhist, max2Dhist, nbins2D, min2Dhist2, max2Dhist2);
      if(!h_TMstudies_2D_track[icalo])
        h_TMstudies_2D_track[icalo]             = new TH2F(Form("h_TMstudies_%s_track_%s",str2DhistoName.Data(),str_calorimeter[icalo].Data()), "", nbins2D, min2Dhist, max2Dhist, nbins2D, min2Dhist2, max2Dhist2);
      if(!h_TMstudies_2D_true[icalo])
        h_TMstudies_2D_true[icalo]              = new TH2F(Form("h_TMstudies_%s_true_%s",str2DhistoName.Data(),str_calorimeter[icalo].Data()), "", nbins2D, min2Dhist, max2Dhist, nbins2D, min2Dhist2, max2Dhist2);
      if(!h_TMstudies_2D_true_electron[icalo])
        h_TMstudies_2D_true_electron[icalo]     = new TH2F(Form("h_TMstudies_%s_true_electron_%s",str2DhistoName.Data(),str_calorimeter[icalo].Data()), "", nbins2D, min2Dhist, max2Dhist, nbins2D, min2Dhist2, max2Dhist2);
      if(!h_TMstudies_2D_true_hadron[icalo])
        h_TMstudies_2D_true_hadron[icalo]       = new TH2F(Form("h_TMstudies_%s_true_hadron_%s",str2DhistoName.Data(),str_calorimeter[icalo].Data()), "", nbins2D, min2Dhist, max2Dhist, nbins2D, min2Dhist2, max2Dhist2);
      if(!h_TMstudies_clusSpec[icalo])
        h_TMstudies_clusSpec[icalo]             = new TH1F(Form("h_TMstudies_clusSpec_%s",str_calorimeter[icalo].Data()), "", 120,0,60);
      if(!h_TMstudies_clusSpec_Matched[icalo])
        h_TMstudies_clusSpec_Matched[icalo]     = new TH1F(Form("h_TMstudies_clusSpec_Matched_%s",str_calorimeter[icalo].Data()), "", 120,0,60);

      // add histos for HCal cluster to ECal cluster matching
      if (IsHCALCalorimeter(icalo)){
        if(!h_TMstudies_clusSpec_MatchedCalo[icalo])
          h_TMstudies_clusSpec_MatchedCalo[icalo] = new TH1F(Form("h_TMstudies_clusSpec_MatchedCalo_%s",str_calorimeter[icalo].Data()), "", 120,0,60);
        if(!h_TMstudies_clusSpec_MatchedCalo_EDiff1[icalo])
          h_TMstudies_clusSpec_MatchedCalo_EDiff1[icalo]  = new TH1F(Form("h_TMstudies_clusSpec_MatchedCalo_EDiff1_%s",str_calorimeter[icalo].Data()), "", 120,0,60);
        if(!h_TMstudies_clusSpec_MatchedCalo_EDiff2[icalo])
          h_TMstudies_clusSpec_MatchedCalo_EDiff2[icalo]  = new TH1F(Form("h_TMstudies_clusSpec_MatchedCalo_EDiff2_%s",str_calorimeter[icalo].Data()), "", 120,0,60);
      }

      // add histos with special cuts to clusters
      if (runSpecialCuts){
        if(!h_TMstudies_clusSpec_special[icalo])
          h_TMstudies_clusSpec_special[icalo]     = new TH1F(Form("h_TMstudies_clusSpec_special_%s",str_calorimeter[icalo].Data()), "", 120,0,60);
        if(!h_TMstudies_clusSpec_special_Matched[icalo])
          h_TMstudies_clusSpec_special_Matched[icalo]     = new TH1F(Form("h_TMstudies_clusSpec_special_Matched_%s",str_calorimeter[icalo].Data()), "", 120,0,60);
      }

      // add eta dependent hists
      for (Int_t et = minEtaBinCalo[caloDir]; et<maxEtaBinCalo[caloDir]; et++){
        Double_t etaMin = partEtaCalo[0];
        Double_t etaMax = partEtaCalo[nEta];
        if (et < nEta){
          etaMin = partEtaCalo[et];
          etaMax = partEtaCalo[et+1];
        }
        // reconstructed  particles
        if(!h_TMstudies_clusSpec_Eta[icalo][et])
          h_TMstudies_clusSpec_Eta[icalo][et]         = new TH1F(Form("h_TMstudies_clusSpec_Eta_%1.1f_%1.1f_%s",etaMin,etaMax,str_calorimeter[icalo].Data()), "", 120,0,60);
        // reconstructed & matched particles
        if(!h_TMstudies_clusSpec_Matched_Eta[icalo][et])
          h_TMstudies_clusSpec_Matched_Eta[icalo][et] = new TH1F(Form("h_TMstudies_clusSpec_Matched_Eta_%1.1f_%1.1f_%s",etaMin,etaMax,str_calorimeter[icalo].Data()), "", 120,0,60);
      }

      // bookkeeping hists for simple effi estimates
      for(int icase=0;icase<diffCases_TMStudies;icase++){
        if(!h_TMstudies_isMatched[icase][icalo]){
          h_TMstudies_isMatched[icase][icalo]               = new TH1F(Form("h_TMstudies_%s_isMatched_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data()), "", _active_algo*5,0,_active_algo*5);
          for(int ialgo=0;ialgo<_active_algo;ialgo++){
            h_TMstudies_isMatched[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+1, Form("%s clusters", str_clusterizer[ialgo].Data()));
            h_TMstudies_isMatched[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+2, "track matched");
            h_TMstudies_isMatched[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+3, "proj matched");
            h_TMstudies_isMatched[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+4, "hit matched");
            h_TMstudies_isMatched[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+5, "calo matched");
          }
        }
        if(!h_TMstudies_isMatched_PDG[icase][icalo]){
          h_TMstudies_isMatched_PDG[icase][icalo]           = new TH2F(Form("h_TMstudies_%s_isMatched_PDG_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data()), "", _active_algo*5,0,_active_algo*5,2500,0,2500);
          for(int ialgo=0;ialgo<_active_algo;ialgo++){
            h_TMstudies_isMatched_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+1, Form("%s clusters", str_clusterizer[ialgo].Data()));
            h_TMstudies_isMatched_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+2, "track matched");
            h_TMstudies_isMatched_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+3, "proj matched");
            h_TMstudies_isMatched_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+4, "hit matched");
            h_TMstudies_isMatched_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+5, "calo matched");
          }
        }
        if(!h_TMstudies_isMatched_E[icase][icalo]){
          h_TMstudies_isMatched_E[icase][icalo]             = new TH2F(Form("h_TMstudies_%s_isMatched_E_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data()), "", _active_algo*5,0,_active_algo*5,100,0,100);
          for(int ialgo=0;ialgo<_active_algo;ialgo++){
            h_TMstudies_isMatched_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+1, Form("%s clusters", str_clusterizer[ialgo].Data()));
            h_TMstudies_isMatched_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+2, "track matched");
            h_TMstudies_isMatched_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+3, "proj matched");
            h_TMstudies_isMatched_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+4, "hit matched");
            h_TMstudies_isMatched_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+5, "calo matched");
          }
        }

        if(!h_TMstudies_isMatched_findable[icase][icalo]){
          h_TMstudies_isMatched_findable[icase][icalo]      = new TH1F(Form("h_TMstudies_%s_isMatched_findable_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data()), "", _active_algo*5,0,_active_algo*5);
          for(int ialgo=0;ialgo<_active_algo;ialgo++){
            h_TMstudies_isMatched_findable[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+1, Form("%s clusters", str_clusterizer[ialgo].Data()));
            h_TMstudies_isMatched_findable[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+2, "track matched");
            h_TMstudies_isMatched_findable[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+3, "proj matched");
            h_TMstudies_isMatched_findable[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+4, "hit matched");
            h_TMstudies_isMatched_findable[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+5, "calo matched");
          }
        }
        if(!h_TMstudies_isMatched_findable_E[icase][icalo]){
          h_TMstudies_isMatched_findable_E[icase][icalo]    = new TH2F(Form("h_TMstudies_%s_isMatched_findable_E_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data()), "", _active_algo*5,0,_active_algo*5,100,0,100);
          for(int ialgo=0;ialgo<_active_algo;ialgo++){
            h_TMstudies_isMatched_findable_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+1, Form("%s clusters", str_clusterizer[ialgo].Data()));
            h_TMstudies_isMatched_findable_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+2, "track matched");
            h_TMstudies_isMatched_findable_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+3, "proj matched");
            h_TMstudies_isMatched_findable_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+4, "hit matched");
            h_TMstudies_isMatched_findable_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*5+4, "calo matched");
          }
        }
        if(!h_TMstudies_isMatched_findable_PDG[icase][icalo]){
          h_TMstudies_isMatched_findable_PDG[icase][icalo]  = new TH2F(Form("h_TMstudies_%s_isMatched_findable_PDG_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data()), "", _active_algo*5,0,_active_algo*5,2500,0,2500);
          for(int ialgo=0;ialgo<_active_algo;ialgo++){
            h_TMstudies_isMatched_findable_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+1, Form("%s clusters", str_clusterizer[ialgo].Data()));
            h_TMstudies_isMatched_findable_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+2, "track matched");
            h_TMstudies_isMatched_findable_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+3, "proj matched");
            h_TMstudies_isMatched_findable_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+4, "hit matched");
            h_TMstudies_isMatched_findable_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*4+5, "hit matched");
          }
        }
      }
      for(int ialgo=0;ialgo<_active_algo;ialgo++){
        for(int icase=0;icase<diffCases_TMStudies;icase++){
          float nbins2Ddelta =  200;
          float min2Ddeltahist =  -0.25;
          float max2Ddeltahist =  0.25;

          h_TMstudies_2D_delta_track[icase][icalo][ialgo]               = new TH2F(Form("h_TMstudies_%s_2D_delta_track_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "",
                                                                                  nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, nbins2Ddelta, min2Ddeltahist, max2Ddeltahist);
          h_TMstudies_2D_delta_projection[icase][icalo][ialgo]          = new TH2F(Form("h_TMstudies_%s_2D_delta_projection_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "",
                                                                                    nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, nbins2Ddelta, min2Ddeltahist, max2Ddeltahist);
          if (runSpecialCuts){
              h_TMstudies_2D_delta_projection_special[icase][icalo][ialgo]  = new TH2F(Form("h_TMstudies_%s_2D_delta_projection_special_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "",
                                                                                    nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, nbins2Ddelta, min2Ddeltahist, max2Ddeltahist);
          }
          h_TMstudies_2D_delta_TTLhits[icase][icalo][ialgo]               = new TH2F(Form("h_TMstudies_%s_2D_delta_TTLhits_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "",
                                                                                    nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, nbins2Ddelta, min2Ddeltahist, max2Ddeltahist);
          if (IsHCALCalorimeter(icalo)){
            h_TMstudies_2D_delta_calomatch[icase][icalo][ialgo]           = new TH2F(Form("h_TMstudies_%s_2D_delta_calomatch_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "",
                                                                                    nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, nbins2Ddelta, min2Ddeltahist, max2Ddeltahist);
            h_TMstudies_E1_E2_calomatch[icase][icalo][ialgo]              = new TH2F(Form("h_TMstudies_%s_E1_E2_calomatch_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 100,0, 50.,100,0, 50.);
          }
          h_TMstudies_dx_vsClsEta_track[icase][icalo][ialgo]            = new TH2F(Form("h_TMstudies_%s_dx_vsClsEta_track_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "",
                                                                                  nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, (Int_t)(TMath::Abs(minEtaCurCalo-maxEtaCurCalo)*10), minEtaCurCalo, maxEtaCurCalo);
          h_TMstudies_dx_vsClsEta_projection[icase][icalo][ialgo]       = new TH2F(Form("h_TMstudies_%s_dx_vsClsEta_projection_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "",
                                                                                  nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, (Int_t)(TMath::Abs(minEtaCurCalo-maxEtaCurCalo)*10), minEtaCurCalo, maxEtaCurCalo);
          h_TMstudies_dx_vsTrkEta_track[icase][icalo][ialgo]            = new TH2F(Form("h_TMstudies_%s_dx_vsTrkEta_track_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "",
                                                                                  nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, (Int_t)(TMath::Abs(minEtaCurCalo-maxEtaCurCalo)*10), minEtaCurCalo, maxEtaCurCalo);
          h_TMstudies_dx_vsTrkEta_projection[icase][icalo][ialgo]       = new TH2F(Form("h_TMstudies_%s_dx_vsTrkEta_projection_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "",
                                                                                    nbins2Ddelta, min2Ddeltahist, max2Ddeltahist, (Int_t)(TMath::Abs(minEtaCurCalo-maxEtaCurCalo)*10), minEtaCurCalo, maxEtaCurCalo);
          for (Int_t et = minEtaBinCalo[caloDir]; et<maxEtaBinCalo[caloDir]; et++){
            Double_t etaMin = partEtaCalo[0];
            Double_t etaMax = partEtaCalo[nEta];
            if (et < nEta){
              etaMin = partEtaCalo[et];
              etaMax = partEtaCalo[et+1];
            }
            h_TMstudies_dx_vsClsE_projection[icase][icalo][ialgo][et] = new TH2F(Form("h_TMstudies_%s_dx_vsClsEta_projection_E_%1.1f_%1.1f_%s_%s",strCases_TMStudies[icase].Data(),etaMin,etaMax, str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP,
                                                                                  nbins2Ddelta, min2Ddeltahist, max2Ddeltahist);
            h_TMstudies_dy_vsClsE_projection[icase][icalo][ialgo][et] = new TH2F(Form("h_TMstudies_%s_dy_vsClsEta_projection_E_%1.1f_%1.1f_%s_%s",strCases_TMStudies[icase].Data(),etaMin,etaMax, str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP,
                                                                                  nbins2Ddelta, min2Ddeltahist, max2Ddeltahist);
          }
        }
        if(!h_TMstudies_2D_clus[icalo][ialgo])
          h_TMstudies_2D_clus[icalo][ialgo]       = new TH2F(Form("h_TMstudies_%s_clus_%s_%s",str2DhistoName.Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nbins2D, min2Dhist, max2Dhist, nbins2D, min2Dhist2, max2Dhist2);
      }
    }
    initializeHistsTMS = true;
  }

  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    Int_t caloDir           = GetCaloDirection(icalo);
    Float_t minEtaCurCalo   = partEtaCalo[minEtaBinCalo[caloDir]];
    Float_t maxEtaCurCalo   = partEtaCalo[maxEtaBinCalo[caloDir]];
  //  bool isFwd = IsForwardCalorimeter(icalo);
    bool b_TMstudies_2D_projection_filled = false;

    // loop over all different algorithms
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(!loadClusterizerInput( ialgo, icalo)) continue;

      int projectionlayer = ReturnProjectionIndexForCalorimeter(icalo, true);

      float matchingwdw = ReturnTrackMatchingWindowForCalo(icalo);

      std::vector<int> usedIDs;
      std::vector<int> currentIDs;
      std::vector<int> usedMClabels;
      for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[ialgo][icalo].size(); iclus++){
        currentIDs.clear();
        bool ismatched_track      = false;
        bool ismatched_projection = false;
        bool ismatched_hit        = false;
        bool ismatched_largest_track      = false;
        bool ismatched_largest_projection = false;
        bool ismatched_calo       = false;
        int mcID                = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_trueID;

        if(((_clusters_calo[ialgo][icalo].at(iclus)).cluster_M02<0.1))continue;
        h_TMstudies_clusSpec[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E);

        // determine eta bin
        Int_t et = 0;
        while (partEtaCalo[et+1] < (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta && et < nEta) et++;

        if(et >= minEtaBinCalo[caloDir] && et<maxEtaBinCalo[caloDir]){
          h_TMstudies_clusSpec_Eta[icalo][et]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E);
        }

        // if(!((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta>3.0))continue;
        // if(((_clusters_calo[ialgo][icalo].at(iclus)).cluster_NTowers<2))continue;
        // if(((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta>3.2))continue;

        // find largest cluster for current MC label and make matching plots for only the largest cluster
        usedIDs.push_back(iclus);
        currentIDs.push_back(iclus);
        for(Int_t icluslbl=0; icluslbl<(Int_t)_clusters_calo[ialgo][icalo].size(); icluslbl++){
          if((int)(mcID)==(int)((_clusters_calo[ialgo][icalo].at(icluslbl)).cluster_trueID)){

            if ( !(std::find(usedIDs.begin(), usedIDs.end(), icluslbl) != usedIDs.end()) ){
              usedIDs.push_back(icluslbl);
              currentIDs.push_back(icluslbl);
            }
          }
        }

        int largestClusterID = -1;
        float largestClusterE = -1;
        for (int id = 0; id < (int)currentIDs.size() ; id++){
          if ( (_clusters_calo[ialgo][icalo].at(currentIDs.at(id))).cluster_E > largestClusterE )
            largestClusterE = (_clusters_calo[ialgo][icalo].at(currentIDs.at(id))).cluster_E;
            largestClusterID = currentIDs.at(id);
        }
        bool fillcurrentlargest = false;
        if ( diffCases_TMStudies>1 && !(std::find(usedMClabels.begin(), usedMClabels.end(), (int)(mcID)) != usedMClabels.end()) ){
          usedMClabels.push_back((int)(mcID));
          fillcurrentlargest = true;
        }

        // bool matchfound = false;
        h_TMstudies_isMatched[0][icalo]->Fill(4*ialgo);
        h_TMstudies_isMatched_E[0][icalo]->Fill(4*ialgo,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E );
        h_TMstudies_isMatched_PDG[0][icalo]->Fill(4*ialgo,abs(_mcpart_PDG[mcID]) );
        if(fillcurrentlargest){
          h_TMstudies_isMatched[1][icalo]->Fill(4*ialgo);
          h_TMstudies_isMatched_E[1][icalo]->Fill(4*ialgo,(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_E );
          h_TMstudies_isMatched_PDG[1][icalo]->Fill(4*ialgo,abs(_mcpart_PDG[(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_trueID]) );
        }
        bool isfindable = true;
        if(abs(_mcpart_PDG[mcID])==22 || abs(_mcpart_PDG[mcID])==2112 || abs(_mcpart_PDG[mcID])==130 || abs(_mcpart_PDG[mcID])==13) isfindable = false;
        if(isfindable){
          h_TMstudies_isMatched_findable[0][icalo]->Fill(4*ialgo);
          h_TMstudies_isMatched_findable_E[0][icalo]->Fill(4*ialgo,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E );
          h_TMstudies_isMatched_findable_PDG[0][icalo]->Fill(4*ialgo,abs(_mcpart_PDG[mcID]) );
          if(fillcurrentlargest){
            h_TMstudies_isMatched_findable[1][icalo]->Fill(4*ialgo);
            h_TMstudies_isMatched_findable_E[1][icalo]->Fill(4*ialgo,(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_E );
            h_TMstudies_isMatched_findable_PDG[1][icalo]->Fill(4*ialgo,abs(_mcpart_PDG[(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_trueID]) );
          }
        }
        h_TMstudies_2D_clus[icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta, (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi);
        for(Int_t iproj=0; iproj<_nProjections; iproj++){
          // perform matching only if projection layer is the correct one for the current calorimeter
          if(_track_ProjLayer[iproj]!=projectionlayer) continue;
          // if timing is -9000 then projection failed
          if(_track_Proj_t[iproj]<-9000) continue;
          // if x and y in forward direction is 0, then projection can not be used
          if( _track_Proj_x[iproj]==0 && _track_Proj_y[iproj]==0) continue;
          // create projection vector from x,y,z positions
          TVector3 projvec(_track_Proj_x[iproj],_track_Proj_y[iproj],_track_Proj_z[iproj]);
          float projeta = projvec.Eta();
          float projphi = projvec.Phi();
          if( projphi==0 && projeta==0) continue;
          // this condition is called only once to fill 2D histograms with all projections
          if(!b_TMstudies_2D_projection_filled){
              h_TMstudies_2D_projection[icalo]->Fill(projeta,projphi);
          }
          // if a cluster is found at x,y=0 (inside the beampipe), then something went wrong
          if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_X==0 && (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Y==0) continue;
          // check deta/dphi or dx/dy difference
          float delta_1 = projeta-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta;
          float delta_2 = projphi-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi;
          h_TMstudies_2D_delta_projection[0][icalo][ialgo]->Fill(delta_1,delta_2);

          if( runSpecialCuts && ((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta)<0.5){
            h_TMstudies_2D_delta_projection_special[0][icalo][ialgo]->Fill(delta_1,delta_2);
          }

          h_TMstudies_dx_vsClsEta_projection[0][icalo][ialgo]->Fill(delta_1,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta);
          h_TMstudies_dx_vsTrkEta_projection[0][icalo][ialgo]->Fill(delta_1,projeta);

          if(h_TMstudies_dx_vsClsE_projection[0][icalo][ialgo][et]) h_TMstudies_dx_vsClsE_projection[0][icalo][ialgo][et]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E, delta_1);
          if(h_TMstudies_dy_vsClsE_projection[0][icalo][ialgo][et]) h_TMstudies_dy_vsClsE_projection[0][icalo][ialgo][et]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E, delta_2);

          float delta_1_lgst = 1000;
          float delta_2_lgst = 1000;
          if(fillcurrentlargest){
            delta_1_lgst =  projeta-(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_Eta;
            delta_2_lgst = projphi-(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_Phi;
            h_TMstudies_2D_delta_projection[1][icalo][ialgo]->Fill(delta_1_lgst,delta_2_lgst);
            h_TMstudies_dx_vsClsEta_projection[1][icalo][ialgo]->Fill(delta_1_lgst,(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_Eta);
            h_TMstudies_dx_vsTrkEta_projection[1][icalo][ialgo]->Fill(delta_1_lgst,projeta);
          }
          // if(TMath::Sqrt(TMath::Power(delta_1,2)+TMath::Power(_track_Proj_y[iproj]-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Y,2)) <15)matchfound=true;
          // h_TMstudies_2D_delta_projection[0][icalo][ialgo]->Fill(projeta-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta,projphi-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi);
          if(abs(delta_1) < matchingwdw){
            if(abs(delta_2) < matchingwdw){
              ismatched_projection = true;
            }
          }
          if(abs(delta_1_lgst) < matchingwdw){
            if(abs(delta_2_lgst) < matchingwdw){
              ismatched_largest_projection = true;
            }
          }
        }
        b_TMstudies_2D_projection_filled = true;

        if(ismatched_projection){
          h_TMstudies_clusSpec_Matched[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E);
        }

        for(Int_t itrk=0; itrk<m_nTracks; itrk++){
          // don't match backward tracks to forward calo and vise versa
          if(ReturnCaloFwdBarrelBckw(icalo)==0 && m_tr_pz[itrk]<0) continue;
          if(ReturnCaloFwdBarrelBckw(icalo)==2 && m_tr_pz[itrk]>0) continue;
          // check eta difference
          TVector3 trackvec(m_tr_px[itrk],m_tr_py[itrk],m_tr_pz[itrk]);

          float tracketa = trackvec.Eta();
          float trackphi = trackvec.Phi(); //(trackvec.Phi()<0 ? trackvec.Phi()+TMath::Pi() : trackvec.Phi()-TMath::Pi());
          if(!b_TMstudies_2D_track_filled[icalo]){
            h_TMstudies_2D_track[icalo]->Fill( tracketa, trackphi);
          }

          if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_X==0 && (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Y==0) continue;
          float delta_1 =  tracketa-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta;
          float delta_2 =trackphi-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi;
          h_TMstudies_2D_delta_track[0][icalo][ialgo]->Fill(delta_1,delta_2);
          h_TMstudies_dx_vsClsEta_track[0][icalo][ialgo]->Fill(delta_1,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta);
          h_TMstudies_dx_vsTrkEta_track[0][icalo][ialgo]->Fill(delta_1,tracketa);
          float delta_1_lgst = 1000;
          float delta_2_lgst = 1000;
          if(fillcurrentlargest){
            delta_1_lgst = tracketa-(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_Eta;
            delta_2_lgst =trackphi-(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_Phi;
            h_TMstudies_2D_delta_track[1][icalo][ialgo]->Fill(delta_1_lgst,delta_2_lgst);
            h_TMstudies_dx_vsClsEta_track[1][icalo][ialgo]->Fill(delta_1_lgst,(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_Eta);
            h_TMstudies_dx_vsTrkEta_track[1][icalo][ialgo]->Fill(delta_1_lgst,tracketa);
          }
          if(abs(delta_1) < matchingwdw){
            if(abs(delta_2) < matchingwdw){
              ismatched_track = true;
            }
          }
          if(abs(delta_1_lgst) < matchingwdw){
            if(abs(delta_2_lgst) < matchingwdw){
              ismatched_largest_track = true;
            }
          }
        }
        b_TMstudies_2D_track_filled[icalo] = true;

        for(Int_t ihit=0; ihit<_nHitsLayers; ihit++){
          if (_hits_layerID[ihit] == -1) continue;
          if( _hits_layerID[ihit]==ReturnProjectionIndexForCalorimeter(icalo, true) ){

            if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_X==0 && (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Y==0) continue;
            h_TMstudies_2D_delta_TTLhits[0][icalo][ialgo]->Fill(_hits_x[ihit]-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_X,_hits_y[ihit]-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Y);
            if(abs(_hits_x[ihit]-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_X) < matchingwdw){
              if(abs(_hits_y[ihit]-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Y) < matchingwdw){
                ismatched_hit = true;
              }
            }
          }
        }

        if(_combCalo[icalo] != -1  && IsHCALCalorimeter(icalo)){
          if (caloEnabled[_combCalo[icalo]]){
            for(Int_t iclus2=0; iclus2<(Int_t)_clusters_calo[ialgo][_combCalo[icalo]].size(); iclus2++){
              if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_X==0 && (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Y==0) continue;
              if((_clusters_calo[ialgo][_combCalo[icalo]].at(iclus2)).cluster_X==0 && (_clusters_calo[ialgo][_combCalo[icalo]].at(iclus2)).cluster_Y==0) continue;
              // if(mcID==(_clusters_calo[ialgo][_combCalo[icalo]].at(iclus2)).cluster_trueID){
              // if(_mcpart_PDG[mcID]==_mcpart_PDG[(_clusters_calo[ialgo][_combCalo[icalo]].at(iclus2)).cluster_trueID] && _mcpart_PDG[mcID]==211){
                if(abs((_clusters_calo[ialgo][_combCalo[icalo]].at(iclus2)).cluster_X-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_X) < 15){
                  if(abs((_clusters_calo[ialgo][_combCalo[icalo]].at(iclus2)).cluster_Y-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Y) < 15){
                    ismatched_calo = true;
                    h_TMstudies_2D_delta_calomatch[0][icalo][ialgo]->Fill((_clusters_calo[ialgo][_combCalo[icalo]].at(iclus2)).cluster_X-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_X,(_clusters_calo[ialgo][_combCalo[icalo]].at(iclus2)).cluster_Y-(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Y);
                    h_TMstudies_E1_E2_calomatch[0][icalo][ialgo]->Fill((_clusters_calo[ialgo][_combCalo[icalo]].at(iclus2)).cluster_E,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E);
                    // if(abs(_mcpart_PDG[mcID])==abs(_mcpart_PDG[(_clusters_calo[ialgo][_combCalo[icalo]].at(iclus2)).cluster_trueID]) && abs(_mcpart_PDG[mcID])==211)
                    h_TMstudies_clusSpec_MatchedCalo[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E);
                    if((_clusters_calo[ialgo][_combCalo[icalo]].at(iclus2)).cluster_E>(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E)
                      h_TMstudies_clusSpec_MatchedCalo_EDiff1[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E);
                    if((_clusters_calo[ialgo][_combCalo[icalo]].at(iclus2)).cluster_E<(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E)
                      h_TMstudies_clusSpec_MatchedCalo_EDiff2[icalo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E);
                  }
                }
              // }
            }
          }
        }
        // b_TMstudies_2D_track_filled[icalo] = true;

        if(ismatched_track){
          h_TMstudies_isMatched[0][icalo]->Fill(4*ialgo+1);
          h_TMstudies_isMatched_E[0][icalo]->Fill(4*ialgo+1,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E );
          h_TMstudies_isMatched_PDG[0][icalo]->Fill(4*ialgo+1,abs(_mcpart_PDG[mcID]) );
          if(isfindable){
            h_TMstudies_isMatched_findable[0][icalo]->Fill(4*ialgo+1);
            h_TMstudies_isMatched_findable_E[0][icalo]->Fill(4*ialgo+1,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E );
            h_TMstudies_isMatched_findable_PDG[0][icalo]->Fill(4*ialgo+1,abs(_mcpart_PDG[mcID]) );
          }
        }
        if(fillcurrentlargest && ismatched_largest_track){
          h_TMstudies_isMatched[1][icalo]->Fill(4*ialgo+1);
          h_TMstudies_isMatched_E[1][icalo]->Fill(4*ialgo+1,(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_E );
          h_TMstudies_isMatched_PDG[1][icalo]->Fill(4*ialgo+1,abs(_mcpart_PDG[(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_trueID]) );
          if(isfindable){
            h_TMstudies_isMatched_findable[1][icalo]->Fill(4*ialgo+1);
            h_TMstudies_isMatched_findable_E[1][icalo]->Fill(4*ialgo+1,(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_E );
            h_TMstudies_isMatched_findable_PDG[1][icalo]->Fill(4*ialgo+1,abs(_mcpart_PDG[(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_trueID]) );
          }
        }
        if(ismatched_projection){
          h_TMstudies_isMatched[0][icalo]->Fill(4*ialgo+2);
          h_TMstudies_isMatched_E[0][icalo]->Fill(4*ialgo+2,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E );
          h_TMstudies_isMatched_PDG[0][icalo]->Fill(4*ialgo+2,abs(_mcpart_PDG[mcID]) );
          if(isfindable){
            h_TMstudies_isMatched_findable[0][icalo]->Fill(4*ialgo+2);
            h_TMstudies_isMatched_findable_E[0][icalo]->Fill(4*ialgo+2,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E );
            h_TMstudies_isMatched_findable_PDG[0][icalo]->Fill(4*ialgo+2,abs(_mcpart_PDG[mcID]) );
          }
        }
        if(fillcurrentlargest && ismatched_largest_projection){
          h_TMstudies_isMatched[1][icalo]->Fill(4*ialgo+2);
          h_TMstudies_isMatched_E[1][icalo]->Fill(4*ialgo+2,(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_E );
          h_TMstudies_isMatched_PDG[1][icalo]->Fill(4*ialgo+2,abs(_mcpart_PDG[(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_trueID]) );
          if(isfindable){
            h_TMstudies_isMatched_findable[1][icalo]->Fill(4*ialgo+2);
            h_TMstudies_isMatched_findable_E[1][icalo]->Fill(4*ialgo+2,(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_E );
            h_TMstudies_isMatched_findable_PDG[1][icalo]->Fill(4*ialgo+2,abs(_mcpart_PDG[(_clusters_calo[ialgo][icalo].at(largestClusterID)).cluster_trueID]) );
          }
        }
        if(ismatched_hit){
          h_TMstudies_isMatched[0][icalo]->Fill(4*ialgo+3);
          h_TMstudies_isMatched_E[0][icalo]->Fill(4*ialgo+3,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E );
          h_TMstudies_isMatched_PDG[0][icalo]->Fill(4*ialgo+3,abs(_mcpart_PDG[mcID]) );
          if(isfindable){
            h_TMstudies_isMatched_findable[0][icalo]->Fill(4*ialgo+3);
            h_TMstudies_isMatched_findable_E[0][icalo]->Fill(4*ialgo+3,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E );
            h_TMstudies_isMatched_findable_PDG[0][icalo]->Fill(4*ialgo+3,abs(_mcpart_PDG[mcID]) );
          }
        }
        if(ismatched_calo){
          h_TMstudies_isMatched[0][icalo]->Fill(4*ialgo+4);
          h_TMstudies_isMatched_E[0][icalo]->Fill(4*ialgo+4,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E );
          h_TMstudies_isMatched_PDG[0][icalo]->Fill(4*ialgo+4,abs(_mcpart_PDG[mcID]) );
          if(isfindable){
            h_TMstudies_isMatched_findable[0][icalo]->Fill(4*ialgo+4);
            h_TMstudies_isMatched_findable_E[0][icalo]->Fill(4*ialgo+4,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E );
            h_TMstudies_isMatched_findable_PDG[0][icalo]->Fill(4*ialgo+4,abs(_mcpart_PDG[mcID]) );
          }
        }
      }
    }
  }
  //************************************************************************************************************
  // Fill histograms with MC variables
  //************************************************************************************************************
  float trueeta = -1000;
  float truephi = -1000;
  for(int imc=0; imc<_nMCPart; imc++){
    TVector3 mcpartvec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
    for(int icalo=0;icalo<_active_calo;icalo++){
      if (!caloEnabled[icalo]) continue;

        trueeta = _mcpart_Eta[imc];
        truephi = mcpartvec.Phi();


      if(ReturnCaloFwdBarrelBckw(icalo)==0 && _mcpart_pz[imc]>0){
        h_TMstudies_2D_true[icalo]->Fill(mcpartvec.X(),mcpartvec.Y());
        if(abs(_mcpart_PDG[imc])==11 && _mcpart_E[imc]>0.3)h_TMstudies_2D_true_electron[icalo]->Fill(mcpartvec.X(),mcpartvec.Y());
        if((abs(_mcpart_PDG[imc])==211 || abs(_mcpart_PDG[imc])==321 || abs(_mcpart_PDG[imc])==2212) && _mcpart_E[imc]>0.3)h_TMstudies_2D_true_hadron[icalo]->Fill(mcpartvec.X(),mcpartvec.Y());
      } else if(ReturnCaloFwdBarrelBckw(icalo)==1){
        h_TMstudies_2D_true[icalo]->Fill(trueeta,truephi);
        if(abs(_mcpart_PDG[imc])==11 && _mcpart_E[imc]>0.3)h_TMstudies_2D_true_electron[icalo]->Fill(trueeta,truephi);
        if((abs(_mcpart_PDG[imc])==211 || abs(_mcpart_PDG[imc])==321 || abs(_mcpart_PDG[imc])==2212) && _mcpart_E[imc]>0.3)h_TMstudies_2D_true_hadron[icalo]->Fill(trueeta,truephi);
      } else if(ReturnCaloFwdBarrelBckw(icalo)==2 && _mcpart_pz[imc]<0){
        h_TMstudies_2D_true[icalo]->Fill(mcpartvec.X(),mcpartvec.Y());
        if(abs(_mcpart_PDG[imc])==11 && _mcpart_E[imc]>0.3)h_TMstudies_2D_true_electron[icalo]->Fill(mcpartvec.X(),mcpartvec.Y());
        if((abs(_mcpart_PDG[imc])==211 || abs(_mcpart_PDG[imc])==321 || abs(_mcpart_PDG[imc])==2212) && _mcpart_E[imc]>0.3)h_TMstudies_2D_true_hadron[icalo]->Fill(mcpartvec.X(),mcpartvec.Y());
      }
    }
  }
  return false;
}

//************************************************************************************************************
//************************************************************************************************************
// Function to save the histograms for the track matching
//************************************************************************************************************
//************************************************************************************************************
void trackmatchingstudiesSave(){

  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_TMSTUD.root",outputDir.Data()),"RECREATE");
  // fileOutput->mkdir("etabins");
  for(int icalo=0;icalo<_active_calo;icalo++){

    if (!caloEnabled[icalo]) continue;
    Int_t caloDir           = GetCaloDirection(icalo);

    fileOutput->mkdir(Form("%s",str_calorimeter[icalo].Data()));
    fileOutput->cd(Form("%s",str_calorimeter[icalo].Data()));
    for(int icase=0;icase<diffCases_TMStudies;icase++){
      if(h_TMstudies_2D_track[icalo]&&icase==0)h_TMstudies_2D_track[icalo]->Write();
      if(h_TMstudies_2D_projection[icalo]&&icase==0)h_TMstudies_2D_projection[icalo]->Write();
      if(h_TMstudies_2D_true[icalo]&&icase==0)h_TMstudies_2D_true[icalo]->Write();
      if(h_TMstudies_2D_true_electron[icalo]&&icase==0)h_TMstudies_2D_true_electron[icalo]->Write();
      if(h_TMstudies_2D_true_hadron[icalo]&&icase==0)h_TMstudies_2D_true_hadron[icalo]->Write();
      if(h_TMstudies_clusSpec[icalo]&&icase==0)h_TMstudies_clusSpec[icalo]->Write();
      if(h_TMstudies_clusSpec_Matched[icalo]&&icase==0)h_TMstudies_clusSpec_Matched[icalo]->Write();
      if(h_TMstudies_clusSpec_MatchedCalo[icalo]&&icase==0)h_TMstudies_clusSpec_MatchedCalo[icalo]->Write();
      if(h_TMstudies_clusSpec_MatchedCalo_EDiff1[icalo]&&icase==0)h_TMstudies_clusSpec_MatchedCalo_EDiff1[icalo]->Write();
      if(h_TMstudies_clusSpec_MatchedCalo_EDiff2[icalo]&&icase==0)h_TMstudies_clusSpec_MatchedCalo_EDiff2[icalo]->Write();
      if(h_TMstudies_clusSpec_special[icalo]&&icase==0)h_TMstudies_clusSpec_special[icalo]->Write();
      if(h_TMstudies_clusSpec_special_Matched[icalo]&&icase==0)h_TMstudies_clusSpec_special_Matched[icalo]->Write();
      if(h_TMstudies_isMatched[icase][icalo])h_TMstudies_isMatched[icase][icalo]->Write();
      if(h_TMstudies_isMatched_E[icase][icalo])h_TMstudies_isMatched_E[icase][icalo]->Write();
      if(h_TMstudies_isMatched_PDG[icase][icalo])h_TMstudies_isMatched_PDG[icase][icalo]->Write();
      if(h_TMstudies_isMatched_findable[icase][icalo])h_TMstudies_isMatched_findable[icase][icalo]->Write();
      if(h_TMstudies_isMatched_findable_E[icase][icalo])h_TMstudies_isMatched_findable_E[icase][icalo]->Write();
      if(h_TMstudies_isMatched_findable_PDG[icase][icalo])h_TMstudies_isMatched_findable_PDG[icase][icalo]->Write();
      for(int ialgo=0;ialgo<_active_algo;ialgo++){
        if(h_TMstudies_2D_clus[icalo][ialgo]&&icase==0) h_TMstudies_2D_clus[icalo][ialgo]->Write();
        if(h_TMstudies_2D_delta_track[icase][icalo][ialgo]) h_TMstudies_2D_delta_track[icase][icalo][ialgo]->Write();
        if(h_TMstudies_2D_delta_projection[icase][icalo][ialgo]) h_TMstudies_2D_delta_projection[icase][icalo][ialgo]->Write();
        if(h_TMstudies_2D_delta_projection_special[icase][icalo][ialgo]) h_TMstudies_2D_delta_projection_special[icase][icalo][ialgo]->Write();
        if(h_TMstudies_2D_delta_TTLhits[icase][icalo][ialgo]) h_TMstudies_2D_delta_TTLhits[icase][icalo][ialgo]->Write();
        if(h_TMstudies_2D_delta_calomatch[icase][icalo][ialgo]) h_TMstudies_2D_delta_calomatch[icase][icalo][ialgo]->Write();
        if(h_TMstudies_E1_E2_calomatch[icase][icalo][ialgo]) h_TMstudies_E1_E2_calomatch[icase][icalo][ialgo]->Write();
        if(h_TMstudies_dx_vsClsEta_track[icase][icalo][ialgo]) h_TMstudies_dx_vsClsEta_track[icase][icalo][ialgo]->Write();
        if(h_TMstudies_dx_vsClsEta_projection[icase][icalo][ialgo]) h_TMstudies_dx_vsClsEta_projection[icase][icalo][ialgo]->Write();
        if(h_TMstudies_dx_vsTrkEta_track[icase][icalo][ialgo]) h_TMstudies_dx_vsTrkEta_track[icase][icalo][ialgo]->Write();
        if(h_TMstudies_dx_vsTrkEta_projection[icase][icalo][ialgo]) h_TMstudies_dx_vsTrkEta_projection[icase][icalo][ialgo]->Write();
        for (Int_t et = minEtaBinCalo[caloDir]; et<maxEtaBinCalo[caloDir]; et++){
          if(h_TMstudies_dx_vsClsE_projection[icase][icalo][ialgo][et]) h_TMstudies_dx_vsClsE_projection[icase][icalo][ialgo][et]->Write();
          if(h_TMstudies_dy_vsClsE_projection[icase][icalo][ialgo][et]) h_TMstudies_dy_vsClsE_projection[icase][icalo][ialgo][et]->Write();
        }
      }
      for (Int_t et = minEtaBinCalo[caloDir]; et<maxEtaBinCalo[caloDir]; et++){
        if(h_TMstudies_clusSpec_Eta[icalo][et]&&icase==0)h_TMstudies_clusSpec_Eta[icalo][et]->Write();
        if(h_TMstudies_clusSpec_Matched_Eta[icalo][et]&&icase==0)h_TMstudies_clusSpec_Matched_Eta[icalo][et]->Write();
      }
    }
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
