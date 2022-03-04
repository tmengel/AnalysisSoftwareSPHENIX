#ifndef CALOHEADERSPHENIX_H
#define CALOHEADERSPHENIX_H

// ANCHOR basic struct for towers in clusterizer
typedef struct {
  float tower_E;
  int tower_iEta;
  int tower_iPhi;
  int tower_trueID;
} towersStrct;

typedef struct {
  float cluster_E;
  float cluster_seed;
  float cluster_Eta;
  float cluster_Phi;
  float cluster_Z;
  float cluster_X;
  float cluster_Y;
  float cluster_M02;
  float cluster_M20;
  bool cluster_isMatched;
  std::vector<int> cluster_matchedTrackIDs;
  int cluster_NTowers;
  int cluster_trueID;
  int cluster_NtrueID;
} clustersStrct;

typedef struct {
  int particle_ID;
  float highest_E;
  int nClusters;
} occuranceStrct;


const int maxcalo = 3;
int calogeomindex[maxcalo] = {0};

TString str_calorimeter[maxcalo] = {"CEMC", "HCALIN", "HCALOUT"};
int _combCalo[maxcalo]           = {      -1,  kHCALIN,   kHCALOUT};
int _combCalo2[maxcalo]          = {      -1,     -1,      -1};
int _caloTowersPhi[maxcalo]      = {     100,      64,       64};
const int _active_calo = 3;

void CheckBarrelCaloVersion (){
  if (caloEnabled[kCEMC]) {
    _combCalo[kHCALIN]  = kCEMC;
    _combCalo[kHCALOUT] = kCEMC;
  }
}

enum clusterizertype {
    kMA         = 0,
    kV3         = 1,
    kV1         = 2,
    k5x5        = 3,
    kC5         = 4,
    kC3         = 5,
    k3x3        = 6,
    kDummy      = 7
};
//                      CEMC  HCALIN  HCALOUT
float seedE[maxcalo]  = {0.5, 0.2,    0.5 };
float aggE[maxcalo]   = {0.1,  0.05,   0.1};

// TString str_clusterizer[7] = {"V1", "V3", "3x3", "5x5", "C3", "C5", "MA"};
TString str_clusterizer[7] = {"MA", "V3", "V1", "5x5", "C5", "C3", "3x3"};
// TString str_clusterizer[7] = {"V3", "V3", "V1", "5x5", "C5", "C3", "3x3"};
const int maxAlgo = 7;
const int _active_algo = 1;

float _ch_DRCALO_pos_z = 1;
float _ch_FHCAL_pos_z = 1;
float _ch_FEMC_pos_z = 1;
float _ch_EHCAL_pos_z = 1;
float _ch_EEMC_pos_z = 1;
float _ch_EEMCG_pos_z = 1;
float _ch_LFHCAL_pos_z = 1;
float _ch_FOCAL_pos_z = 1;

// sorting function for towers
bool acompare(towersStrct lhs, towersStrct rhs) { return lhs.tower_E > rhs.tower_E; }
bool acompareTrueID(towersStrct lhs, towersStrct rhs) { return lhs.tower_trueID > rhs.tower_trueID; }
bool acompareCl(clustersStrct lhs, clustersStrct rhs) { return lhs.cluster_E > rhs.cluster_E; }

void setINTClusterArrayToZero(int* &arrayinput){
  for(int iclus=0;iclus<_maxNclusters;iclus++){
    arrayinput[iclus]=0;
  }
}
void setFLOATClusterArrayToZero(float* &arrayinput){
  for(int iclus=0;iclus<_maxNclusters;iclus++){
    arrayinput[iclus]=0.;
  }
}
void setBOOLClusterArrayToZero(bool* &arrayinput){
  for(int iclus=0;iclus<_maxNclusters;iclus++){
    arrayinput[iclus]=false;
  }
}
// conversion functions for processing

float weightM02 = 4.5;
//float * CalculateM02andWeightedPosition(std::vector<towersStrct> cluster_towers, float w_0, float cluster_E_calc, int caloSelect, Bool_t debugOutput);


const int _maxNtowers1D = 200;
const int _maxNtowersL = 10;
TVector3 caloPositionArray[maxcalo][_maxNtowers1D/*iEta*/][_maxNtowers1D/*iPhi*/][_maxNtowersL/*iL*/];
void SetGeometryIndices(){
  Long64_t nEntriesTree                 = tt_geometry->GetEntries();
  for (int ireset=0; ireset<maxcalo;ireset++) {
    calogeomindex[ireset] = -1;
  }
  for (Long64_t i=0; i<nEntriesTree;i++) {
    tt_geometry->GetEntry(i);
    calogeomindex[_calogeom_ID] = i;

    for (Long64_t itow=0; itow<_calogeom_towers_N;itow++) {
    //  caloPositionArray[_calogeom_ID][_calogeom_towers_iEta[itow]][_calogeom_towers_iPhi[itow]][_calogeom_towers_iL[itow]] = TVector3(_calogeom_towers_x[itow],_calogeom_towers_y[itow],_calogeom_towers_z[itow]);
    }
  }

}

bool IsHCALCalorimeter(int caloID){
  switch (caloID){

    case kHCALIN: return true;
    case kHCALOUT: return true;
    case kCEMC: return false;

    default:
      std::cout << "IsHCALCalorimeter: caloID " << caloID << " not defined, returning false" << std::endl;
      return false;
  }
}


int GetCaloDirection(int caloID){
  switch (caloID){

    case kHCALIN: return 1;
    case kHCALOUT: return 1;
    case kCEMC: return 1;
    default:
      std::cout << "GetCaloDirection: caloID " << caloID << " not defined, returning -1" << std::endl;
      return 1;
  }

}




int ReturnCaloFwdBarrelBckw(int caloID){
  // 0 for forward calos
  // 1 for barrel calos
  // 2 for backwards calos
  switch (caloID){

    case kHCALIN: return 1;
    case kHCALOUT: return 1;
    case kCEMC: return 1;

    default:
      return 1;
  }

}


float* EtaPhiFromIndices(int ieta,int iphi,float energy = 0, int caloSelect = 0);


// ANCHOR function to return a TVector3 for the tower position based on iEta and iPhi indices

  // cout << caloSelect << "\t" << i_Eta<< "\t" << i_Phi<< "\t" << i_L << endl;
  // cout << "\tvec: " << twrPositionVec.x() << "\t" << twrPositionVec.y()<< "\t" << twrPositionVec.z() << endl;


int ReturnMaxTowerCalo(int caloID){
  switch (caloID){

    case kHCALIN: return _nTowers_HCALIN;
    case kHCALOUT: return _nTowers_HCALOUT;
    case kCEMC: return _nTowers_CEMC;

    default:
      std::cout << "ReturnMaxTowerCalo: caloID " << caloID << " not defined, returning -1" << std::endl;
      return -1;
  }
}



int ReturnProjectionIndexForCalorimeter(int caloID,  bool alternate){
  if (!alternate){
    switch (caloID){

      case kHCALIN: return 5;
      case kHCALOUT: return 7;
      case kCEMC: return 6;

      default:
        std::cout << "ReturnProjectionIndexForCalorimeter: caloID " << caloID << " not defined, returning -1" << std::endl;
        return -1;
    }
  }
  else {
    switch (caloID){

      case kHCALIN: return 5;
      case kHCALOUT: return 7;
      case kCEMC: return 6;

      default:
        std::cout << "ReturnProjectionIndexForCalorimeter: caloID " << caloID << " not defined, returning -1" << std::endl;
        return -1;
    }

  }

}

int ReturnCalorimeterFromProjectionIndex(int projID){
  switch (projID){

    case 5: return kHCALIN;
    case 7: return kHCALOUT;
    case 6: return kCEMC;
    default:
      // std::cout << "ReturnCalorimeterFromProjectionIndex: projID " << projID << " not defined, returning -1" << std::endl;
      return -1;
  }

}

float ReturnTrackMatchingWindowForCalo(int caloID){
  switch (caloID){

    case kHCALIN: return 0.1;
    case kHCALOUT: return 0.1;
    case kCEMC: return 0.05;

    default:
      std::cout << "ReturnTrackMatchingWindowForCalo: caloID " << caloID << " not defined, returning -1" << std::endl;
      return -1;
  }

}

//
//// ANCHOR function to determine shower shape
//float * CalculateM02andWeightedPosition(std::vector<towersStrct> cluster_towers, float w_0, float cluster_E_calc, int caloSelect, bool debugOutput){
//   static float returnVariables[4]; //0:M02, 1:M20, 2:eta, 3: phi
//   float w_tot = 0;
//   std::vector<float> w_i;
//
//
//   float zHC = 1;
//
//  std::vector<float> vecTwr = {0.,0.,0.};
//   //calculation of weights and weighted position vector
//   for(int cellI=0; cellI<(int)cluster_towers.size(); cellI++){
//       w_i.push_back(TMath::Max( (float)0, (float) (w_0 + TMath::Log(cluster_towers.at(cellI).tower_E/cluster_E_calc) )));
//       w_tot += w_i.at(cellI);
//
////         std::cout << caloSelect << "\t" << vecTwrTmp.X() << "\t" << vecTwrTmp.Y() << "\t" << vecTwrTmp.Z() << std::endl;
//       vecTwr += w_i.at(cellI)*vecTwr;
//
//   }
//   returnVariables[2]=vecTwr.Eta();
//   returnVariables[3]=vecTwr.Phi(); //(vecTwr.Phi()<0 ? vecTwr.Phi()+TMath::Pi() : vecTwr.Phi()-TMath::Pi());
////     std::cout << "X: "<< vecTwr.X() << "\t" << " Y: "<< vecTwr.Y() << "\t" << " Z: "<< vecTwr.Z() << "\t zHC: " <<  zHC << std::endl;
//
//   // returnVariables[4]=vecTwr.X();
//   // returnVariables[5]=vecTwr.Y();
//   // returnVariables[6]=vecTwr.Z();
//
//   //calculation of M02
//   float delta_phi_phi[4] = {0};
//   float delta_eta_eta[4] = {0};
//   float delta_eta_phi[4] = {0};
//
//   for(int cellI=0; cellI<(int)cluster_towers.size(); cellI++){
//       int iphi=cluster_towers.at(cellI).tower_iPhi;
//       int ieta=cluster_towers.at(cellI).tower_iEta;
//       delta_phi_phi[1] += (w_i.at(cellI)*iphi*iphi)/w_tot;
//       delta_phi_phi[2] += (w_i.at(cellI)*iphi)/w_tot;
//       delta_phi_phi[3] += (w_i.at(cellI)*iphi)/w_tot;
//
//       delta_eta_eta[1] += (w_i.at(cellI)*ieta*ieta)/w_tot;
//       delta_eta_eta[2] += (w_i.at(cellI)*ieta)/w_tot;
//       delta_eta_eta[3] += (w_i.at(cellI)*ieta)/w_tot;
//
//       delta_eta_phi[1] += (w_i.at(cellI)*ieta*iphi)/w_tot;
//       delta_eta_phi[2] += (w_i.at(cellI)*iphi)/w_tot;
//       delta_eta_phi[3] += (w_i.at(cellI)*ieta)/w_tot;
//   }
//   delta_phi_phi[0] = delta_phi_phi[1] - (delta_phi_phi[2] * delta_phi_phi[3]);
//   delta_eta_eta[0] = delta_eta_eta[1] - (delta_eta_eta[2] * delta_eta_eta[3]);
//   delta_eta_phi[0] = delta_eta_phi[1] - (delta_eta_phi[2] * delta_eta_phi[3]);
//
//   float calcM02 = 0.5 * ( delta_phi_phi[0] + delta_eta_eta[0] ) + TMath::Sqrt( 0.25 * TMath::Power( ( delta_phi_phi[0] - delta_eta_eta[0] ), 2 ) + TMath::Power( delta_eta_phi[0], 2 ) );
//   float calcM20 = 0.5 * ( delta_phi_phi[0] + delta_eta_eta[0] ) - TMath::Sqrt( 0.25 * TMath::Power( ( delta_phi_phi[0] - delta_eta_eta[0] ), 2 ) + TMath::Power( delta_eta_phi[0], 2 ) );
//   if(debugOutput) std::cout << "M02_calc: " << calcM02 << "\t\t = 0.5 * ( " << delta_phi_phi[0] <<" + "<<delta_eta_eta[0]<<" ) + TMath::Sqrt( 0.25 * TMath::Power( ( "<<delta_phi_phi[0]<<" - "<<delta_eta_eta[0]<<" ), 2 ) + TMath::Power( "<<delta_eta_phi[0]<<", 2 ) ) "<< std::endl;
//   returnVariables[0]=calcM02;
//   returnVariables[1]=calcM20;
//   return returnVariables;
//
//}


// 3x3 global cluster variables
std::vector<clustersStrct> _clusters_calo[maxAlgo][maxcalo];

bool loadClusterizerInput(
  int clusterizerEnum,
  int caloEnum
){
  // std::cout << clusterizerEnum << "\t" << caloEnum << std::std::endl;
  if (clusterizerEnum != kV1 ){
    if (!_clusters_calo[clusterizerEnum][caloEnum].empty() && caloEnabled[caloEnum])
      return true;
    else
      return false;
  } else {


      if (!_clusters_calo[clusterizerEnum][caloEnum].empty() && caloEnabled[caloEnum])
        return true;
      else
        return false;

  }

}

float getCalibrationValue(float clusterE, int caloEnum, int algoEnum){


  return 1.0;
}

float getEnergySmearing( int caloEnum, int algoEnum){
  // _fRandom.SetSeed(0);
  // if(caloEnum==kFHCAL){
  //   if(algoEnum==kMA){
  //     return _fRandom.Gaus(1,0.0985);
  //   } else if(algoEnum==kV1){
  //     return _fRandom.Gaus(1,0.10);
  //   } else if(algoEnum==kV3){
  //     return _fRandom.Gaus(1,0.04);
  //   } else {
  //     return 1.0;
  //   }
  // } else {
    return 1.0;
  // }
  // return 1.0;
}

void fillHCalClustersIntoJetFindingInputs( int caloEnum, int clusterizerEnum,
  std::vector<float> & jetf_hcal_E, std::vector<float> & jetf_hcal_px, std::vector<float> & jetf_hcal_py, std::vector<float> & jetf_hcal_pz,
  std::vector<float> & jetf_calo_E, std::vector<float> & jetf_calo_px, std::vector<float> & jetf_calo_py, std::vector<float> & jetf_calo_pz,
  std::vector<float> & jetf_all_E, std::vector<float> & jetf_all_px, std::vector<float> & jetf_all_py, std::vector<float> & jetf_all_pz
)
{
  for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[clusterizerEnum][caloEnum].size(); iclus++){
      if(!(_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_isMatched){
          double pt = (_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E / cosh((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_Eta);
          double px = pt * cos((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_Phi);
          double py = pt * sin((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_Phi);
          double pz = pt * sinh((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_Eta);
          jetf_hcal_px.push_back(px);
          jetf_hcal_py.push_back(py);
          jetf_hcal_pz.push_back(pz);
          jetf_hcal_E.push_back((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E);
          jetf_calo_px.push_back(px);
          jetf_calo_py.push_back(py);
          jetf_calo_pz.push_back(pz);
          jetf_calo_E.push_back((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E);
          jetf_all_px.push_back(px);
          jetf_all_py.push_back(py);
          jetf_all_pz.push_back(pz);
          jetf_all_E.push_back((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E);
      }
  }
}

void fillECalClustersIntoJetFindingInputs(
  int caloEnum, int clusterizerEnum,
  std::vector<float> & jetf_calo_E, std::vector<float> & jetf_calo_px, std::vector<float> & jetf_calo_py, std::vector<float> & jetf_calo_pz,
  std::vector<float> & jetf_all_E, std::vector<float> & jetf_all_px, std::vector<float> & jetf_all_py, std::vector<float> & jetf_all_pz,
  std::vector<float> & jetf_full_E, std::vector<float> & jetf_full_px, std::vector<float> & jetf_full_py, std::vector<float> & jetf_full_pz
)
{
    for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[clusterizerEnum][caloEnum].size(); iclus++){
      if(!(_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_isMatched){
          double pt = (_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E / cosh((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_Eta);
          double px = pt * cos((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_Phi);
          double py = pt * sin((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_Phi);
          double pz = pt * sinh((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_Eta);

          jetf_calo_px.push_back(px);
          jetf_calo_py.push_back(py);
          jetf_calo_pz.push_back(pz);
          jetf_calo_E.push_back((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E);
          jetf_all_px.push_back(px);
          jetf_all_py.push_back(py);
          jetf_all_pz.push_back(pz);
          jetf_all_E.push_back((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E);
          jetf_full_px.push_back(px);
          jetf_full_py.push_back(py);
          jetf_full_pz.push_back(pz);
          jetf_full_E.push_back((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E);
      }
  }
}

#endif
