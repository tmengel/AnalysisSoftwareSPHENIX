#include <TString.h>

  const Int_t nBinsP              = 148;
  const Int_t nBinsPLow           = 99;
  Double_t binningP[149]          = { 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,  //10
                                      0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,         // 20
                                      1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,         // 30
                                      2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,         // 40
                                      3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0,         // 50
                                      4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0,         // 60
                                      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5,    // 70
                                      7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0,   // 80
                                      10.5, 11., 11.5, 12., 12.5, 13., 13.5, 14., 14.5, 15.0,   // 90
                                      15.5, 16., 16.5, 17., 17.5, 18., 18.5, 19., 19.5, 20.,    // 100
                                      21., 22., 23., 24., 25., 26., 27., 28., 29., 30.,         // 110
                                      31., 32., 33., 34., 35., 36., 37., 38., 39., 40.,         // 120
                                      42., 44., 46., 48., 50., 52., 54., 56., 58., 60.,         // 130
                                      65., 70., 75., 80., 85., 90., 95., 100., 110., 120.,      // 140
                                      130., 140., 150., 160., 170., 180., 190., 200., 250.};    // 149
                                       
  const Int_t nCuts               = 12;
  const Int_t nPID                = 6;
  const Int_t nResoSt             = 5;
  
  TString nameCuts[nCuts]            = {"N", "L3", "L3F", "LI", "LIF", "BE", "BEF", "AE", "AEF", "T", "LI2", "LI3" };
  TString nameResoAdd[5]          = {"All", "woT", "wT","LI", "LS"};
  TString nameAddBeta[4]          = {"", "AEMC", "LI3", "LI3AEMC"};
  
  TString labelCuts[12]           = { "no track cuts", "at least 3 hits", "#leq 3 hits & ( 1^{st} | 2^{nd} layer)", "only tracker", "only tracker & ( 1^{st} | 2^{nd} layer)", 
                                      "tracker & only LGAD before ECal", "tracker & only LGAD before ECal & ( 1^{st} | 2^{nd} layer)", "tracker & LGAD after ECal", "tracker & only LGAD after ECal & ( 1^{st} | 2^{nd} layer)", "any timing hit", 
                                      "#geq 2 tracker hits", "#geq 3 tracker hits"  };

  TString partName[6]             = {"All", "Electron", "Muon", "Pion", "Kaon", "Proton" };
  TString partLabel[6]            = {"(h/e)^{#pm}",  "e^{#pm}", "#mu^{#pm}", "#pi^{#pm}", "K^{#pm}", "p/#bar{p}"};

  TString partNameET2[6]            = {"all", "electron", "muon", "cpion", "ckaon", "proton" };
  TString partNameET[6]             = {"chargedpart", "electron", "muon", "cpion", "ckaon", "proton" };
  TString partLabelET[6]            = {"(h/e)^{#pm}",  "e^{#pm}", "#mu^{#pm}", "#pi^{#pm}", "K^{#pm}", "p/#bar{p}"};

  //************************** Read data **************************************************
  const Int_t nPt                  = 13;
  const static Double_t partPt[]     = {0., 0.5, 1.0,  2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 
                                      9.0, 10.0, 15.0, 20.0};
  const Int_t nP                   = 25;
  const static Double_t partP[]    = {0., 0.1, 0.2, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0,
                                      5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15., 20.0, 30.0,
                                      40.0, 50., 75., 100., 150, 200. };

  const Int_t nEta                = 15;                                        
  Double_t partEta[nEta+1]        = { -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.2, -0.4, 0.4, 1.2, 
                                       1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
                                      // 0: -4.0 < eta < -3.5 
                                      // 1: -3.5 < eta < -3.0 
                                      // 2: -3.0 < eta < -2.5 
                                      // 3: -2.5 < eta < -2.0 
                                      // 4: -2.0 < eta < -1.5 
                                      // 5: -1.5 < eta < -1.2 
                                      // 6: -1.2 < eta < -0.4 
                                      // 7: -0.4 < eta < 0.4 
                                      // 8: 0.4 < eta < 1.2 
                                      // 9: 1.2 < eta < 1.5 
                                      // 10: 1.5 < eta < 2 
                                      // 11: 2 < eta < 2.5 
                                      // 12: 2.5 < eta < 3.0 
                                      // 13: 3.0 < eta < 3.5 
                                      // 14: 3.5 < eta < 4.0
  Bool_t enablePlot[nEta+1]               = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                              1, 1, 1, 1, 0, 0};
  const Int_t rebinEta_PtResol[nEta+1]    = { 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 
                                              2, 2, 2, 2, 4, 1};
  const Int_t rebinEta_EtaResol[nEta+1]   = { 8, 8, 4, 1, 1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 4, 8,  1};
  const Int_t rebinEta_PhiResol[nEta+1]   = { 8, 8, 4, 1, 1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 4, 8,  1};
  const Float_t ptResolFit[nEta+1]        = { 0.3, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                              0.1, 0.1, 0.2, 0.3, 0.3, 0.2};                                    
  const Float_t etaResolFit[nEta+1]       = { 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 
                                              0.01, 0.01, 0.01, 0.02, 0.02, 0.01};                                    
  const Float_t phiResolFit[nEta+1]       = { 0.15, 0.15, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                                              0.05, 0.05, 0.1, 0.15, 0.2, 0.1};                                    
  const Int_t rebinEta_BetaResol[nEta+1]  = { 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 
                                              2, 2, 2, 2, 4,  1};
  const Float_t betaResolFit[nEta+1]      = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                                              0.05, 0.05, 0.05, 0.05, 0.05, 0.05};                                    

//   const Int_t rebinEta_PtResol[20]   = { 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 
//                                          2, 2, 2, 2, 2, 2, 2, 2, 4, 1};
//   const Int_t rebinEta_EtaResol[20]   = { 8, 8, 4, 1, 1, 1, 1, 1, 1, 1,
//                                           1, 1, 1, 1, 1, 1, 1, 4, 8,  1};
//   const Int_t rebinEta_PhiResol[20]   = { 8, 8, 4, 1, 1, 1, 1, 1, 1, 1,
//                                           1, 1, 1, 1, 1, 1, 1, 4, 8,  1};
//   Float_t ptResolFit[20]    = { 0.3, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
//                                 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3, 0.3, 0.2};                                    
//   Float_t etaResolFit[20]   = { 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 
//                                 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.01};                                    
//   Float_t phiResolFit[20]   = { 0.15, 0.15, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
//                                 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.15, 0.2, 0.1};                                    
                                
//   const Int_t rebinEta_BetaResol[20]   = { 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 
//                                          2, 2, 2, 2, 2, 2, 2, 2, 4,  1};
//   Float_t betaResolFit[20]    = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
//                                 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};                                    
  
  Color_t colorEta[nEta+1]         = {kBlue+1, kBlue-6, kViolet+2, kViolet-5, kMagenta-6, kPink-9, kRed+1,  kOrange+7, kOrange, kYellow-6, 
                                  kGreen+1, kGreen-5, kCyan+1, kCyan+3, kAzure+2, kBlack };
  Style_t markerStyleEta[nEta+1]   = {24, 25, 27, 28, 30, 42, 46, 24, 25, 27, 
                                 28, 30, 42, 46, 24, 20};
  Size_t markerSizeEta[nEta+1]     = {1.5, 1.4, 1.9, 1.5, 1.8, 1.8, 1.5, 1.5, 1.4, 1.9,
                                 1.5, 1.8, 1.8, 1.5, 1.5, 1.5 };

  const int nPhi = 16;
  Double_t partPhi[nPhi+1]          = {0.};
  void SetPhiBins(){
    for (Int_t iPhi = 0; iPhi < nPhi+1; iPhi++){
      partPhi[iPhi] = -TMath::Pi() +(2*TMath::Pi())/nPhi*iPhi;
      std::cout << partPhi[iPhi] << ", "; 
    }
    std::cout << std::endl;
  } 
  
  Color_t colorPhi[nPhi+1]         = {kBlue+1, kBlue-6, kViolet+2, kViolet-5, kMagenta-6, kPink-9, kRed+1,  kOrange+7, kOrange, kYellow-6, 
                                  kGreen+1, kGreen-5, kCyan+1, kCyan+3, kAzure+2,  kGray, kBlack };
  Style_t markerStylePhi[nPhi+1]   = {24, 25, 27, 28, 30, 42, 46, 24, 25, 27, 
                                 28, 30, 42, 46, 24, 25, 27};
  Size_t markerSizePhi[nPhi+1]     = {1.5, 1.4, 1.9, 1.5, 1.8, 1.8, 1.5, 1.5, 1.4, 1.9,
                                 1.5, 1.8, 1.8, 1.5, 1.5, 1.5, 1.5 };
                                 
//   Color_t colorEta[20]         = {kBlue+1, kBlue-6, kViolet+2, kViolet-5, kMagenta+1, kMagenta-6, kPink-5, kPink-9, kRed+1, kRed-7, 
//                                   kOrange+7, kOrange, kYellow-6, kSpring+5, kGreen+1, kGreen-5, kCyan+1, kCyan+3, kAzure+2, kBlack };
//   Style_t markerStyleEta[20]   = {24, 25, 27, 28, 30, 42, 46, 24, 25, 27, 
//                                  28, 30, 42, 46, 24, 25, 27, 28, 30, 20};
//   Size_t markerSizeEta[20]     = {1.5, 1.4, 1.9, 1.5, 1.8, 1.8, 1.5, 1.5, 1.4, 1.9,
//                                  1.5, 1.8, 1.8, 1.5, 1.5, 1.4, 1.9, 1.5, 1.8, 1.5 };
  Color_t colorPID[nPID]          = {kBlack, kRed+1, kGreen+2, kCyan+2, kBlue+1, kOrange};
  Style_t markerStylePID[nPID]    = {20, 24, 25, 27, 28, 30};
  Size_t markerSizePID[nPID]      = {1.5, 1.4, 1.5, 1.9, 1.6, 1.8};
  
  TString nameOutEtaRange[3]   = {"backward", "central", "forward"};
  TString labelEtaRange[3]     = {"     e-going", "     barrel", "     h-going"};
  Int_t minEtaBin[3]           = {1, 4, 10};
  Int_t maxEtaBin[3]           = {4, 8, 13};
  Int_t maxNEtaBins[3]         = {4, 5,  4};
    
  Int_t minEtaBinFull[3]       = {0, 4,  9};
  Int_t maxEtaBinFull[3]       = {5, 8, 14};
  Int_t maxNEtaBinsFull[3]     = {6, 5,  6};

  const Int_t nClus            = 7;
  TString nameClus[nClus]      = {"MA", "V1", "V3", "3x3", "5x5", "C3", "C5" };
  
  Bool_t enableClus[nClus]             = {1, 1, 0, 1, 0, 1, 1};
  Color_t colorClus[3][nClus]           = { {kMagenta+2, kBlack, kRed+1, kGreen+2, kCyan+2, kAzure,  kOrange }, 
                                            {kMagenta-1, kGray+1, kRed-1, kGreen-1, kCyan-1, kAzure+3,  kOrange-3 }, 
                                            {kMagenta-6, kGray+2, kRed-6, kGreen-6, kCyan-6, kAzure-5, kOrange+4}
                                          };
  Style_t markerStyleClus[nClus]    = {46, 20, 24, 25, 27, 28, 30 };
  Size_t markerSizeClus[nClus]      = {1.5, 1.4, 1.5, 1.9, 1.6, 1.8, 1.5};
