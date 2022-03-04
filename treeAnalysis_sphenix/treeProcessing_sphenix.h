#include <TROOT.h>
#include <TString.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TVector3.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>

TString outputDir;
// Constants are from the EventEvaluator class in coresoftware



const int _maxNHits = 10000;
// const int _maxNTowers = 50 * 50;
const int _maxNTowersCentral = 2000;
// const int _maxNTowersDR = 3000 * 3000;
const int _maxNTowersCalo = 5000000;
const int _maxNclusters = 100;
const int _maxNclustersCentral = 2000;
const int _maxNTracks = 200;
const int _maxNProjections = 2000;
const int _maxNMCPart = 100000;
const int _maxNHepmcp = 1000;
int verbosityBASE = 2;

float _nEventsTree;

enum calotype {
  kCEMC           = 0,
  kHCALIN         = 1,
  kHCALOUT        = 2,
};

typedef struct {
  float eta_Calo;
  float phi_Calo;
} projStrct;


bool caloEnabled[3]      = {0};
bool tracksEnabled        = 0;
bool vertexEnabled        = 0;
bool xSectionEnabled      = 0;

// Event level info
float _cross_section;
float _event_weight;
int _n_generator_accepted;




// track hits
int _nHitsLayers;
int* _hits_layerID              = new int[_maxNHits];
float* _hits_x             = new float[_maxNHits];
float* _hits_y             = new float[_maxNHits];
float* _hits_z             = new float[_maxNHits];
float* _hits_t             = new float[_maxNHits];
float* _hits_edep          = new float[_maxNHits];

// towers
int _nTowers_CEMC;
float* _tower_CEMC_E            = new float[_maxNTowersCentral];
int* _tower_CEMC_iEta         = new int[_maxNTowersCentral];
int* _tower_CEMC_iPhi         = new int[_maxNTowersCentral];
int* _tower_CEMC_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_HCALIN;
float* _tower_HCALIN_E            = new float[_maxNTowersCentral];
int* _tower_HCALIN_iEta           = new int[_maxNTowersCentral];
int* _tower_HCALIN_iPhi           = new int[_maxNTowersCentral];
int* _tower_HCALIN_trueID         = new int[_maxNTowersCentral];

// towers
int _nTowers_HCALOUT;
float* _tower_HCALOUT_E           = new float[_maxNTowersCentral];
int* _tower_HCALOUT_iEta          = new int[_maxNTowersCentral];
int* _tower_HCALOUT_iPhi          = new int[_maxNTowersCentral];
int* _tower_HCALOUT_trueID        = new int[_maxNTowersCentral];

// vertex
float _vertex_x;
float _vertex_y;
float _vertex_z;

// tracks
// int _nTracks;
// float* _track_ID                 = new float[_maxNTracks];
// float* _track_trueID             = new float[_maxNTracks];
// float* _track_px                 = new float[_maxNTracks];
// float* _track_py                 = new float[_maxNTracks];
// float* _track_pz                 = new float[_maxNTracks];
//unsigned short* _track_source             = new unsigned short[_maxNTracks];
std::array<std::vector<int>, _maxNTracks> _track_RefProjID;

int m_nTracks;
float* m_tr_px = new float[_maxNTracks];
float* m_tr_py = new float[_maxNTracks];
float* m_tr_pz = new float[_maxNTracks];
float* m_tr_p = new float[_maxNTracks];
float* m_tr_pt = new float[_maxNTracks];
float* m_tr_phi = new float[_maxNTracks];
float* m_tr_eta = new float[_maxNTracks];
int* m_charge = new int[_maxNTracks];
float* m_chisq = new float[_maxNTracks];
int* m_ndf = new int[_maxNTracks];
float* m_dca = new float[_maxNTracks];
float* m_tr_x = new float[_maxNTracks];
float* m_tr_y= new float[_maxNTracks];
float* m_tr_z= new float[_maxNTracks];
float* m_tr_ID = new float[_maxNTracks];

int* m_truth_is_primary= new int[_maxNTracks];
float* m_truthtrackpx = new float[_maxNTracks];
float* m_truthtrackpy = new float[_maxNTracks];
float* m_truthtrackpz = new float[_maxNTracks];
float* m_truthtrackp = new float[_maxNTracks];
float* m_truthtracke = new float[_maxNTracks];
float* m_truthtrackpt = new float[_maxNTracks];
float* m_truthtrackphi = new float[_maxNTracks];
float* m_truthtracketa = new float[_maxNTracks];
int* m_truthtrackpid = new int[_maxNTracks];
float* m_truthtrackID = new float[_maxNTracks];

// Initializes all elements to 0
std::array<bool, _maxNTracks> _track_hasTTL{{}};
std::array<int, _maxNTracks> _track_nTTL{{}};
std::array<int, _maxNTracks> _track_nTrL{{}};
std::array<bool, _maxNTracks> _track_hasIL{{}};
std::array<bool, _maxNTracks> _track_hasOL{{}};

int _nProjections;
float* _track_ProjTrackID        = new float[_maxNProjections];
int* _track_ProjLayer            = new int[_maxNProjections];
float* _track_Proj_x             = new float[_maxNProjections];
float* _track_Proj_y             = new float[_maxNProjections];
float* _track_Proj_z             = new float[_maxNProjections];
float* _track_Proj_t             = new float[_maxNProjections];
float* _track_Proj_true_x        = new float[_maxNProjections];
float* _track_Proj_true_y        = new float[_maxNProjections];
float* _track_Proj_true_z        = new float[_maxNProjections];
float* _track_Proj_true_t        = new float[_maxNProjections];
std::array<int, _maxNProjections> _track_Proj_Clas{{}};



// MC particles
int _nMCPart;
int* _mcpart_ID                   = new int[_maxNMCPart];
int* _mcpart_ID_parent            = new int[_maxNMCPart];
int* _mcpart_PDG                  = new int[_maxNMCPart];
float* _mcpart_E                  = new float[_maxNMCPart];
float* _mcpart_px                 = new float[_maxNMCPart];
float* _mcpart_py                 = new float[_maxNMCPart];
float* _mcpart_pz                 = new float[_maxNMCPart];
float* _mcpart_Eta                = new float[_maxNMCPart];
float* _mcpart_Phi                = new float[_maxNMCPart];
std::array<std::vector<int>, _maxNMCPart> _mcpart_RecTrackIDs;
std::array<std::vector<projStrct>, _maxNMCPart> _mcpart_EcalProjs;
std::array<std::vector<projStrct>, _maxNMCPart> _mcpart_HcalProjs;

TRandom3  _fRandom;                                  // random for effi generation


int _calogeom_ID;
int _calogeom_towers_N;
int*  _calogeom_towers_iEta = new int[_maxNTowersCalo];
int*  _calogeom_towers_iPhi = new int[_maxNTowersCalo];
//int*  _calogeom_towers_iL   = new int[_maxNTowersCalo];
float*  _calogeom_towers_Eta = new float[_maxNTowersCalo];
float*  _calogeom_towers_Phi = new float[_maxNTowersCalo];
float*  _calogeom_towers_x = new float[_maxNTowersCalo];
float*  _calogeom_towers_y = new float[_maxNTowersCalo];
float*  _calogeom_towers_z = new float[_maxNTowersCalo];

TTree * tt_geometry;

void SetBranchAddressesGeometryTree(TTree* inputTreeGeo){
    inputTreeGeo->SetBranchAddress("m_GEO_ID",              &_calogeom_ID);
    inputTreeGeo->SetBranchAddress("m_GEO_Ntowers",     &_calogeom_towers_N);
    inputTreeGeo->SetBranchAddress("m_GEO_towers_iEta",  _calogeom_towers_iEta);
    inputTreeGeo->SetBranchAddress("m_GEO_towers_iPhi",  _calogeom_towers_iPhi);
    //inputTreeGeo->SetBranchAddress("calo_towers_iL",    _calogeom_towers_iL);
    inputTreeGeo->SetBranchAddress("m_GEO_towers_Eta",   _calogeom_towers_Eta);
    inputTreeGeo->SetBranchAddress("m_GEO_towers_Phi",   _calogeom_towers_Phi);
    inputTreeGeo->SetBranchAddress("m_GEO_towers_x",     _calogeom_towers_x);
    inputTreeGeo->SetBranchAddress("m_GEO_towers_y",     _calogeom_towers_y);
    inputTreeGeo->SetBranchAddress("m_GEO_towers_z",     _calogeom_towers_z);
}
void SetBranchAddressesTree(TTree* inputTree){

    if (inputTree->GetBranchStatus("m_cross_section") ){
      xSectionEnabled = 1;
      inputTree->SetBranchAddress("m_cross_section", &_cross_section);
      inputTree->SetBranchAddress("m_event_weight", &_event_weight);
      inputTree->SetBranchAddress("m_n_generator_accepted", &_n_generator_accepted);
    }
    if (inputTree->GetBranchStatus("m_nHits") ){
      inputTree->SetBranchAddress("m_nHits",                        &_nHitsLayers);
      inputTree->SetBranchAddress("m_hits_layerID",                 _hits_layerID);
      inputTree->SetBranchAddress("m_hits_x",               _hits_x);
      inputTree->SetBranchAddress("m_hits_y",               _hits_y);
      inputTree->SetBranchAddress("m_hits_z",               _hits_z);
      inputTree->SetBranchAddress("m_hits_t",               _hits_t);
      //inputTree->SetBranchAddress("hits_edep",               _hits_edep);
    }

    if (inputTree->GetBranchStatus("m_nTracks") ){
      tracksEnabled = 1;
      inputTree->SetBranchAddress("m_nTracks", &m_nTracks );
      inputTree->SetBranchAddress("m_tr_px", m_tr_px);
      inputTree->SetBranchAddress("m_tr_py", m_tr_py);
      inputTree->SetBranchAddress("m_tr_pz", m_tr_pz);
      inputTree->SetBranchAddress("m_tr_p", m_tr_p);
      inputTree->SetBranchAddress("m_tr_pt", m_tr_pt);
      inputTree->SetBranchAddress("m_tr_phi", m_tr_phi);
      inputTree->SetBranchAddress("m_tr_eta", m_tr_eta);
      inputTree->SetBranchAddress("m_charge", m_charge);
      inputTree->SetBranchAddress("m_chisq", m_chisq);
      inputTree->SetBranchAddress("m_ndf", m_ndf);
      inputTree->SetBranchAddress("m_dca", m_dca);
      inputTree->SetBranchAddress("m_tr_x", m_tr_x);
      inputTree->SetBranchAddress("m_tr_y", m_tr_y);
      inputTree->SetBranchAddress("m_tr_z", m_tr_z);
      inputTree->SetBranchAddress("m_tr_ID", m_tr_ID);

      inputTree->SetBranchAddress("m_truth_is_primary", m_truth_is_primary);
      inputTree->SetBranchAddress("m_truthtrackpx", m_truthtrackpx);
      inputTree->SetBranchAddress("m_truthtrackpy", m_truthtrackpy);
      inputTree->SetBranchAddress("m_truthtrackpz", m_truthtrackpz);
      inputTree->SetBranchAddress("m_truthtrackp", m_truthtrackp);
      inputTree->SetBranchAddress("m_truthtracke", m_truthtracke);
      inputTree->SetBranchAddress("m_truthtrackpt", m_truthtrackpt);
      inputTree->SetBranchAddress("m_truthtrackphi", m_truthtrackphi);
      inputTree->SetBranchAddress("m_truthtracketa", m_truthtracketa);
      inputTree->SetBranchAddress("m_truthtrackpid", m_truthtrackpid);
      inputTree->SetBranchAddress("m_truthtrackID", m_truthtrackID);


      inputTree->SetBranchAddress("m_nProjections",         &_nProjections);
      inputTree->SetBranchAddress("m_track_ProjTrackID",    _track_ProjTrackID);
      inputTree->SetBranchAddress("m_track_ProjLayer",      _track_ProjLayer);

      inputTree->SetBranchAddress("m_track_TLP_x",           _track_Proj_x);
      inputTree->SetBranchAddress("m_track_TLP_y",           _track_Proj_y);
      inputTree->SetBranchAddress("m_track_TLP_z",           _track_Proj_z);
      inputTree->SetBranchAddress("m_track_TLP_t",           _track_Proj_t);
      inputTree->SetBranchAddress("m_track_TLP_true_x",      _track_Proj_true_x);
      inputTree->SetBranchAddress("m_track_TLP_true_y",      _track_Proj_true_y);
      inputTree->SetBranchAddress("m_track_TLP_true_z",      _track_Proj_true_z);
      inputTree->SetBranchAddress("m_track_TLP_true_t",      _track_Proj_true_t);
    }


      // towers HCALIN
    if( inputTree->GetBranchStatus("m_Ntowers_HCALIN") ){
      caloEnabled[kHCALIN] = 1;
      inputTree->SetBranchAddress("m_Ntowers_HCALIN",                &_nTowers_HCALIN);
      inputTree->SetBranchAddress("m_tower_E_HCALIN",                _tower_HCALIN_E);
      inputTree->SetBranchAddress("m_tower_iEta_HCALIN",             _tower_HCALIN_iEta);
      inputTree->SetBranchAddress("m_tower_iPhi_HCALIN",             _tower_HCALIN_iPhi);
      inputTree->SetBranchAddress("m_tower_trueID_HCALIN",           _tower_HCALIN_trueID);
    }

    // towers HCALOUT
    if( inputTree->GetBranchStatus("m_Ntowers_HCALOUT") ){
      caloEnabled[kHCALOUT] = 1;
      inputTree->SetBranchAddress("m_Ntowers_HCALOUT",                &_nTowers_HCALOUT);
      inputTree->SetBranchAddress("m_tower_E_HCALOUT",                _tower_HCALOUT_E);
      inputTree->SetBranchAddress("m_tower_iEta_HCALOUT",             _tower_HCALOUT_iEta);
      inputTree->SetBranchAddress("m_tower_iPhi_HCALOUT",             _tower_HCALOUT_iPhi);
      inputTree->SetBranchAddress("m_tower_trueID_HCALOUT",           _tower_HCALOUT_trueID);
    }

    // towers CEMC
    if( inputTree->GetBranchStatus("m_Ntowers_CEMC") ){
      caloEnabled[kCEMC] = 1;
      inputTree->SetBranchAddress("m_Ntowers_CEMC",                &_nTowers_CEMC);
      inputTree->SetBranchAddress("m_tower_E_CEMC",                _tower_CEMC_E);
      inputTree->SetBranchAddress("m_tower_iEta_CEMC",             _tower_CEMC_iEta);
      inputTree->SetBranchAddress("m_tower_iPhi_CEMC",             _tower_CEMC_iPhi);
      inputTree->SetBranchAddress("m_tower_trueID_CEMC",           _tower_CEMC_trueID);
    }




    if (inputTree->GetBranchStatus("m_vertex_x") ){
      vertexEnabled = 1;
      inputTree->SetBranchAddress("m_vertex_x",                     &_vertex_x);
      inputTree->SetBranchAddress("m_vertex_y",                     &_vertex_y);
      inputTree->SetBranchAddress("m_vertex_z",                     &_vertex_z);
    }
    // MC particles
    inputTree->SetBranchAddress("m_nMCparts",       &_nMCPart);
    inputTree->SetBranchAddress("m_MCpart_ID",     _mcpart_ID);
    inputTree->SetBranchAddress("m_MCpart_ID_parent",     _mcpart_ID_parent);
    inputTree->SetBranchAddress("m_MCpart_PDG",    _mcpart_PDG);
    inputTree->SetBranchAddress("m_MCpart_E",      _mcpart_E);
    inputTree->SetBranchAddress("m_MCpart_px",     _mcpart_px);
    inputTree->SetBranchAddress("m_MCpart_py",     _mcpart_py);
    inputTree->SetBranchAddress("m_MCpart_pz",     _mcpart_pz);
}



//__________________________________________________________________________________________________________
TString ReturnDateStr(){
    TDatime today;
    int iDate           = today.GetDate();
    int iYear           = iDate/10000;
    int iMonth          = (iDate%10000)/100;
    int iDay            = iDate%100;
    return Form("%i_%02d_%02d",iYear, iMonth, iDay);
}

bool _do_TimingStudies        = false;
bool _is_ALLSILICON           = false;
const int _maxProjectionLayers    = 60;
Int_t layerIndexHist[7]  =  {0, 1, 2, 3, 4, 5, 6};

// void ResetLayerIndexForward(){
//   if (!_is_ALLSILICON){
//     layerIndexHist[0]  = 50;
//     layerIndexHist[1]  = 51;
//     layerIndexHist[2]  = 52;
//     layerIndexHist[3]  = 53;
//     layerIndexHist[4]  = 54;
//   }
// }

TString GetProjectionNameFromIndex(int projindex)
{
  switch (projindex)
  {
    case 1:
      return "MVTX";
    case 2:
      return "INTT";
    case 3:
      return "TPC";
    case 4:
      return "MICROMEGAS";
    case 5:
      return "HCALIN";
    case 6:
      return "CEMC";
    case 7:
      return "HCALOUT";
    default:
      return "NOTHING";
  }
}

Int_t GetRegionFromIndex(int projindex)
{
  switch (projindex){
    // timing layers
    case 1:     return 1; // "CTTL_1";
    case 2:     return 1; //"FTTL_1";
    case 3:     return 1; //"FTTL_2";
    case 4:     return 1; // "ETTL_0";
    case 5:     return 1; // "ETTL_1";
    case 6:     return 1; // "ETTL_1";
    case 7:     return 1; // "CTTL_0";

    default:   return 1;
  }
}
//
Int_t ReturnIndexForwardLayer(Int_t layerID){
  for (Int_t i = 0; i < _maxProjectionLayers; i++){
    if (layerIndexHist[i] == layerID)
      return i;
  }
  return -1;
}
//

//Bool_t HasFirstTwoLayers(Int_t layerID){
//  if (_is_ALLSILICON){
//    switch (layerID){
//      case 20:
//      case 21:
//      case 30:
//      case 31:
//      case 10:
//      case 11:
//        return kTRUE;
//      default:
//        return kFALSE;
//    }
//  } else {
//    switch (layerID){
//      case 50:  // FST_0
//      case 51:  // FST_1
//      case 100: // EFST_0
//      case 101: // EFST_1
//      case 154: // SVTX_0
//      case 155: // SVTX_1
//      case 156: // SVTX_2
//      case 161: // SVTX
//      case 40:  // BARREL_0
//      case 41:  // BARREL_1
//        return kTRUE;
//      default:
//        return kFALSE;
//    }
//  }
//  return kFALSE;
//}

Bool_t IsTrackerLayer(Int_t layerID){
  if (_is_ALLSILICON){
    switch (layerID){
      case 1:
      case 2:
      case 3:
      case 4:
//      case 24:
//      case 30:
//      case 31:
//      case 32:
//      case 33:
//      case 34:
//      case 10:
//      case 11:
//      case 12:
//      case 13:
//      case 14:
//      case 15:
        return kTRUE;
      default:
        return kFALSE;
    }
  } else {
    switch (layerID){
        case 1:
        case 2:
        case 3:
        case 4:
//      // FST
//      case 50:
//      case 51:
//      case 52:
//      case 53:
//      case 54:
//      case 55:
//      // BARR by layer
//      case 150:
//      case 151:
//      case 152:
//      case 153:
//      case 154:
//      // SVTX by layer
//      case 155:
//      case 156:
//      case 157:
//      case 158:
//      case 159:
//      // BARR
//      case 160:
//      // SVTX
//      case 161:
//      // EFST
//      case 100:
//      case 101:
//      case 102:
//      case 103:
//      case 104:
//      case 105:
//      case 106:
//      // RWELL
//      case 110:
//      case 111:
//      case 112:
//      // FGEM
//      case 120:
//      case 121:
//      // EGEM
//      case 130:
//      case 131:
        return kTRUE;
      default:
        return kFALSE;
    }
  }
  return kFALSE;
}

Bool_t IsCaloProjection(Int_t layerID){
  switch (layerID){
//    case 60:
    case 5:
    case 6:
    case 7:
//    case 64:
//    case 65:
//    case 66:
//    case 67:
//    case 140:
//    case 141:
//    case 142:
//    case 143:
//    case 144:
//    case 145:
//    case 146:
//    case 147:
      return kTRUE;
    default:
      return kFALSE;
  }
  return kFALSE;
}

Bool_t IsECalProjection(Int_t layerID, bool alternate = false){
  if (!alternate){
    switch (layerID){

//      case 0:
      case 6:
//      case 65:
//      case 66:
        return kTRUE;
      default:
        return kFALSE;
    }
  } else {
    switch (layerID){
      case 6:
        return kTRUE;
      default:
        return kFALSE;
    }
  }
  return kFALSE;
}

Bool_t IsHCalProjection(Int_t layerID, bool alternate = false){
  if (!alternate){
    switch (layerID){
//      case 1:
//      case 2:
      case 5:
      case 7:
//      case 67:
//      case 140:
//      case 141:
//      case 142:
//      case 143:
//      case 144:
//      case 145:
//      case 146:
//      case 147:
          return kTRUE;
      default:
        return kFALSE;
    }
  }

  return kFALSE;
}

//Bool_t IsFarForwardProjection(Int_t layerID){
//switch (layerID){
//    case 70:
//    case 71:
//    case 72:
//    case 73:
//    case 74:
//    case 90:
//    case 91:
//    case 92:
//      return kTRUE;
//    default:
//      return kFALSE;
//  }
//  return kFALSE;
//}

Bool_t HasTimingLayer(Int_t layerID){
//  switch (layerID){
//    case 0:
//    case 1:
//    case 2:
//    case 3:
//    case 4:
//    case 5:
//    case 6:
//      return kTRUE;
//    default:
//      return kFALSE;
//  }
  return kFALSE;
}

// **********************************************************************************************
// ****************** resetting of MC references  ***********************************************
// **********************************************************************************************
int GetCorrectMCArrayEntry(float objectTrueID){
  for(Int_t imc=0; imc<_nMCPart; imc++){
    if(objectTrueID==_mcpart_ID[imc]){
      return imc;
    }
  }
  return -1;
}


// **********************************************************************************************
// ****************** create vectors for matching diff rec tracks to MC part ********************
// **********************************************************************************************
void prepareMCMatchInfo(){
  for(Int_t itrk=0; itrk<(Int_t)m_nTracks; itrk++){
    if (verbosityBASE > 3) std::cout << "processing track: " << itrk << std::endl;
      _mcpart_RecTrackIDs[(int)m_truthtrackID[itrk]].push_back(itrk);
  }
  Int_t nCurrProj = 0;
  for(Int_t itrk=0; itrk<m_nTracks; itrk++){
    if (verbosityBASE > 2) std::cout << "current track: " << itrk <<std::endl;
    for(Int_t iproj=nCurrProj; iproj<_nProjections; iproj++){
      if (itrk != _track_ProjTrackID[iproj])
        continue;
      double projectionR = TMath::Sqrt(_track_Proj_x[iproj]*_track_Proj_x[iproj]+_track_Proj_y[iproj]*_track_Proj_y[iproj]);       if(TMath::Abs(_track_Proj_t[iproj])< 2.e-20){
        if (verbosityBASE > 5) std::cout << iproj << "\t projection layer: "<< _track_ProjLayer[iproj] << "\t t: " << _track_Proj_t[iproj]
                                            << "\t x: " << _track_Proj_x[iproj] << "\t y: " << _track_Proj_y[iproj] << "\t z: " << _track_Proj_z[iproj] << "\t r: " << projectionR << std::endl;

        continue;
      }
      if (verbosityBASE > 2) std::cout << iproj << "\t projection layer: "<< _track_ProjLayer[iproj] << "\t t: " << _track_Proj_t[iproj]
                                            << "\t x: " << _track_Proj_x[iproj] << "\t y: " << _track_Proj_y[iproj] << "\t z: " << _track_Proj_z[iproj] << "\t r: " << projectionR << std::endl;
      _track_hasTTL[itrk]     = (_track_hasTTL[itrk] || HasTimingLayer(_track_ProjLayer[iproj]));
      if (HasTimingLayer(_track_ProjLayer[iproj])){
        _track_nTTL[itrk]++;
        _track_Proj_Clas[iproj] = 2;
      }
      if (verbosityBASE > 3) std::cout << "timing layer count: " << _track_nTTL[itrk] << std::endl;
      if (IsTrackerLayer(_track_ProjLayer[iproj])){
        _track_nTrL[itrk]++;
        _track_Proj_Clas[iproj] = 1;
      }
//      if (HasFirstTwoLayers(_track_ProjLayer[iproj])){
//        _track_Proj_Clas[iproj] = 4;
//      }
      if (IsCaloProjection(_track_ProjLayer[iproj])){
        _track_Proj_Clas[iproj] = 5;
      }
//      if (IsFarForwardProjection(_track_ProjLayer[iproj])){
//        _track_Proj_Clas[iproj] = 6;
//      }
      _track_RefProjID[itrk].push_back(iproj);

      projStrct tempProj;
      if ( (IsECalProjection(_track_ProjLayer[iproj],true) || IsHCalProjection(_track_ProjLayer[iproj],true))){
        TVector3 projvec(_track_Proj_true_x[iproj],_track_Proj_true_y[iproj],_track_Proj_true_z[iproj]);
        tempProj.eta_Calo = projvec.Eta();
        tempProj.phi_Calo = projvec.Phi();
        _mcpart_EcalProjs[(int)m_truthtrackID[itrk]].push_back(tempProj);
        _mcpart_HcalProjs[(int)m_truthtrackID[itrk]].push_back(tempProj);
      }
      nCurrProj = iproj;
    }


    if ((int)m_truthtrackID[itrk] < 0) continue;
      if (verbosityBASE > 2) std::cout << "found: " << _track_RefProjID[itrk].size() << "\t projections " << std::endl;

  }
}

// **********************************************************************************************
// ****************** cleanup vectors for matching diff rec tracks to MC part ********************
// **********************************************************************************************
void clearMCRecMatchVectors(){
  if (verbosityBASE > 2) std::cout << "clearing MC vectors" << std::endl;
  for(int imc=0; imc<_nMCPart+2 && imc < _maxNMCPart; imc++){
    if (_mcpart_RecTrackIDs[imc].size() > 0){
      _mcpart_RecTrackIDs[imc].clear();
      _mcpart_RecTrackIDs[imc].resize(0);
    }
    if (_mcpart_EcalProjs[imc].size() > 0){
       _mcpart_EcalProjs[imc].clear();
      _mcpart_EcalProjs[imc].resize(0);
    }
    if (_mcpart_HcalProjs[imc].size() > 0){
       _mcpart_HcalProjs[imc].clear();
      _mcpart_HcalProjs[imc].resize(0);
    }
  }

  for(Int_t itrk=0; itrk<m_nTracks; itrk++){
    _track_hasTTL[itrk]   = 0;
    _track_nTTL[itrk]     = 0;
    _track_nTrL[itrk]     = 0;
    _track_hasIL[itrk]    = 0;
    _track_hasOL[itrk]    = 0;
    if (_track_RefProjID[itrk].size() > 0){
      _track_RefProjID[itrk].clear();
      _track_RefProjID[itrk].resize(0);
    }
  }
  for (Int_t iproj=0; iproj<_nProjections; iproj++){
    _track_Proj_Clas[iproj] = 0;
  }
}
