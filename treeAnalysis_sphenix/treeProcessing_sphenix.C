
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TChain.h>
#include <TMath.h>
#include <TVector3.h>
#include <iostream>
#include <fstream>

#include "binningheader.h"
#include "plottingheader.h"
#include "treeProcessing_sphenix.h"
#include "event_utils.cxx"
#include "jet_finder.cxx"
#include "jet_observables_sphenix.cxx"
#include "caloheader_sphenix.h"
#include "clusterizer_sphenix.cxx"

#include "jetresolutionhistos_sphenix.cxx"
#include "caloresolutionhistos_sphenix.cxx"
#include "clusterstudies_sphenix.cxx"
#include "trackingefficiency.cxx"
#include "hitstudies_sphenix.cxx"
#include "trackmatchingstudies_sphenix.cxx"
#include "pi0studies_sphenix.cxx"




void tannerEventprocessing(
    TString inFile              = "",
    TString inFileGeometry      = "geometry.root",
    TString addOutputName       = "",
    Double_t maxNEvent          = -1,
    bool do_reclus              = true,
    bool do_jetfinding          = false,
    bool runCaloRes             = true,
    bool doCalibration          = false,
    Int_t verbosity             = 0,
    // Defaults to tracking from all layers.
    unsigned short primaryTrackSource = 0,
    std::string jetAlgorithm    = "anti-kt",
    double jetR                 = 0.5,
    double tracked_jet_max_pT   = 30,
    bool removeTracklets        = false
){
    // make output directory
    TString dateForOutput = ReturnDateStr();
    outputDir = Form("treeProcessing/%s",addOutputName.Data());
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec("mkdir -p "+outputDir + "/etaphi");

    if(do_jetfinding) _do_jetfinding = true;

    bool runPi0Reco = true;

    // load tree
    TChain *const tt_event = new TChain("event_tree");
    if (inFile.EndsWith(".root")) {                     // are we loading a single root tree?
        std::cout << "loading a single root file" << std::endl;
        tt_event->AddFile(inFile);
    }
    else {                                              // or are we loading a bunch?
        std::cout << "loading a list of files" << std::endl;
        std::ifstream files(inFile);
        std::string filePath;

        while (std::getline(files, filePath)) {
            tt_event->AddFile(filePath.c_str());
        }
        files.close();
    }
    if(!tt_event){ std::cout << "tree not found... returning!"<< std::endl; return;}

    // // load geometry tree
    tt_geometry =  (TTree *) (new TFile(inFileGeometry.Data(), "READ"))->Get("geometry_tree");
    if(!tt_geometry){ std::cout << "geometry tree not found... returning!"<< std::endl; return;}
    // load all branches (see header)
    SetBranchAddressesTree(tt_event);
    SetBranchAddressesGeometryTree(tt_geometry);
    SetGeometryIndices();

    for (Int_t c = 0; c < 12; c++){
      std::cout << str_calorimeter[c] << "\t" << caloEnabled[c] << std::endl;
    }

    Long64_t nEntriesTree                 = tt_event->GetEntries();
    std::cout << "Number of events in tree: " << nEntriesTree << std::endl;
    if(maxNEvent>0 && maxNEvent<nEntriesTree){
        nEntriesTree = maxNEvent;
        std::cout << "Will only analyze first " << maxNEvent << " events in the tree..." << std::endl;
    }
    _doClusterECalibration    = doCalibration;
    if(_doClusterECalibration){
        std::cout << "clusters will be energy-corrected and subsequently smeared to meet testbeam constant term!" << std::endl;
    }
    // Additional setup
    auto eventObservables = EventObservables();
    auto jetObservablesTrue = JetObservables{JetType_t::full, "true"};
    auto jetObservablesTrueCharged = JetObservables{JetType_t::charged, "true"};
    auto jetObservablesCharged = JetObservables{JetType_t::charged};
    auto jetObservablesCalo = JetObservables{JetType_t::calo};
    auto jetObservablesFull = JetObservables{JetType_t::full};

    if (_do_jetfinding) {
        // Event level
        eventObservables.Init();

        // TODO: Add fully to arguments
        std::vector<double> jetRParameters = {0.5};
        jetObservablesTrue.Init(jetRParameters);
        jetObservablesTrueCharged.Init(jetRParameters);
        jetObservablesCharged.Init(jetRParameters);
        jetObservablesCalo.Init(jetRParameters);
        jetObservablesFull.Init(jetRParameters);
    }

    _nEventsTree=0;
    // main event loop
    for (Long64_t i=0; i<nEntriesTree;i++) {
        // load current event
        tt_event->GetEntry(i);
        _nEventsTree++;

        // processing progress info
        if(i>0 && nEntriesTree>100 && i%(nEntriesTree/(50))==0) std::cout << "//processed " << 100*(i)/nEntriesTree << "%"  << std::endl;
        if(verbosity>0){
          std::cout << "***********************************************************************************************************" << std::endl;
          std::cout << "event " << i << std::endl;
          std::cout << "***********************************************************************************************************" << std::endl;
        }
        // calculate useful quantities
        for(Int_t imc=0; imc<_nMCPart; imc++){
            TVector3 truevec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
            _mcpart_Eta[imc]=truevec.Eta();
            _mcpart_Phi[imc]=truevec.Phi();
            int motherid = 0;
            for(Int_t imc2=0; imc2<_nMCPart; imc2++){
                if(_mcpart_ID[imc2]==_mcpart_ID_parent[imc]){
                    motherid=imc2;
                }
            }
            // if(_mcpart_PDG[imc]==111) cout << "pi0 found (" << _mcpart_ID[imc] << ") with mother: " << _mcpart_PDG[motherid] << " (" <<  _mcpart_ID_parent[imc] << ")" << endl;
        }
        if (verbosity > 0){
          for(int imc=0; imc<_nMCPart; imc++){
            std::cout << "MC part: \t" << imc << "\t E = " << _mcpart_E[imc] << "\tEta: " << _mcpart_Eta[imc]<< "\tPhi: " << _mcpart_Phi[imc]<<  std::endl;
          }
          for (Int_t icalo = 0; icalo < maxcalo; icalo++){
            if (caloEnabled[icalo])
              std::cout << "towers " <<   str_calorimeter[icalo].Data() << "\t" << ReturnMaxTowerCalo(icalo) << std::endl;
          }
        }

        Int_t nTowers[maxcalo] = {_nTowers_CEMC, _nTowers_HCALIN, _nTowers_HCALOUT};

        // run clusterizers normal calos
        for (int cal = 0; cal < maxcalo; cal++){
          if(do_reclus && nTowers[cal] > 0 && caloEnabled[cal]){
            for (int algo = 0; algo < _active_algo; algo++){
              if(verbosity>1) std::cout << "clusterizing " << str_clusterizer[algo].Data() <<  " for " << str_calorimeter[cal].Data() << std::endl;
              runclusterizer(algo, cal, seedE[cal], aggE[cal], primaryTrackSource);
            }
          }
        }


        for(Int_t itrk=0; itrk<_nTracks; itrk++){
            _track_trueID[itrk] = GetCorrectMCArrayEntry(_track_trueID[itrk]);
            // if(verbosity>1) std::cout << "\tTrack: track " << itrk << "\twith true ID " << _track_trueID[itrk] << "\tand X = " << _track_px[itrk] << " cm" << std::endl;
        }
        // obtain labels for different track sources
        prepareMCMatchInfo();

        // ANCHOR Hits loop variables:
        // float* _hits_x[ihit]
        // float* _hits_y[ihit]
        // float* _hits_z[ihit]
        // float* _hits_t[ihit]
        // float* _hits_edep[ihit]
        // for(Int_t ihit=0; ihit<_nHitsLayers; ihit++){
        //     if(verbosity>1) std::cout << "\tHIT: hit " << ihit << "\tin layer " << _hits_layerID[ihit] << "\twith X = " << _hits_x[ihit] << " cm" << std::endl;
        // }
        if(tracksEnabled) hitstudies(primaryTrackSource);


                // ANCHOR Track loop variables:
        // float* _track_ID[itrk]
        // float* _track_trueID[itrk]
        // float* _track_px[itrk]
        // float* _track_py[itrk]
        // float* _track_pz[itrk]
        std::vector<float> jetf_track_px;
        std::vector<float> jetf_track_py;
        std::vector<float> jetf_track_pz;
        std::vector<float> jetf_track_E;
        std::vector<float> jetf_full_px;
        std::vector<float> jetf_full_py;
        std::vector<float> jetf_full_pz;
        std::vector<float> jetf_full_E;
        std::vector<float> jetf_all_px;
        std::vector<float> jetf_all_py;
        std::vector<float> jetf_all_pz;
        std::vector<float> jetf_all_E;
        double massHypothesis = 0.139;
        if(_do_jetfinding){
            for(Int_t itrk=0; itrk<_nTracks; itrk++){
                // Skip tracks that aren't selected.
                if(_track_source[itrk] != primaryTrackSource) {
                    continue;
                }
                if (removeTracklets && (_track_source[itrk] == 0 && _track_hasIL[itrk])) {
                  continue;
                }

                TVector3 trackvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
                float Etrack = TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)+TMath::Power(_track_pz[itrk],2) + TMath::Power(massHypothesis, 2));
                // create track vector for jet finder
                float track_pT = Etrack / cosh(trackvec.Eta());
                if (track_pT < tracked_jet_max_pT) {
                    jetf_track_px.push_back(_track_px[itrk]);
                    jetf_track_py.push_back(_track_py[itrk]);
                    jetf_track_pz.push_back(_track_pz[itrk]);
                    jetf_track_E.push_back(Etrack);
                }
                jetf_full_px.push_back(_track_px[itrk]);
                jetf_full_py.push_back(_track_py[itrk]);
                jetf_full_pz.push_back(_track_pz[itrk]);
                jetf_full_E.push_back(Etrack);
                jetf_all_px.push_back(_track_px[itrk]);
                jetf_all_py.push_back(_track_py[itrk]);
                jetf_all_pz.push_back(_track_pz[itrk]);
                jetf_all_E.push_back(Etrack);
            }
        }



        // ANCHOR jet variables:
        std::vector<float> jetf_hcal_E;
        std::vector<float> jetf_hcal_px;
        std::vector<float> jetf_hcal_py;
        std::vector<float> jetf_hcal_pz;
        std::vector<float> jetf_calo_E;
        std::vector<float> jetf_calo_px;
        std::vector<float> jetf_calo_py;
        std::vector<float> jetf_calo_pz;



        for (Int_t icalo = 0; icalo < maxcalo; icalo++){
          for (Int_t ialgo = 0; ialgo < maxAlgo; ialgo++){
            if (loadClusterizerInput(ialgo, icalo) && verbosity>2) std::cout << str_calorimeter[icalo].Data() << "\t" << ialgo << std::endl;
          }
        }

        if (verbosity > 4){
        }

        // apply calibration if desired
        if(kV1<_active_algo){
          // TODO: Should kV1 in the calibration actually be _active_algo?

        }
        if(do_reclus && kMA<_active_algo && _do_jetfinding){
          // ANCHOR FEMC cluster loop variables:
          // for(Int_t iclus=0; iclus<_nclusters_FEMC; iclus++){
          //     if(verbosity>1) std::cout << "\tFEMC:  cluster " << iclus << "\twith E = " << (_clusters_calo[kV1][kFEMC].at(iclus)).cluster_E << " GeV" << std::endl;
          // }


          // CEMC
          if (caloEnabled[kCEMC]){
            fillECalClustersIntoJetFindingInputs(
                kCEMC, kMA,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz,
                jetf_full_E, jetf_full_px, jetf_full_py, jetf_full_pz
            );
          }
        // HCALIN
          if (caloEnabled[kHCALIN]){
            fillHCalClustersIntoJetFindingInputs(
                kHCALIN, kMA,
                jetf_hcal_E, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz
            );
          }
          // HCALOUT
          if (caloEnabled[kHCALOUT]){
            fillHCalClustersIntoJetFindingInputs(
                kHCALOUT, kMA,
                jetf_hcal_E, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz
            );
          }
        }

        // ANCHOR MC particle loop variables:
        // float * _mcpart_ID[imc]
        // float * _mcpart_ID_parent[imc]
        // float * _mcpart_PDG[imc]
        // float * _mcpart_E [imc]
        // float * _mcpart_px[imc]
        // float * _mcpart_py[imc]
        // float * _mcpart_pz[imc]
        if(_do_jetfinding){
            std::vector<float> jetf_truth_px;
            std::vector<float> jetf_truth_py;
            std::vector<float> jetf_truth_pz;
            std::vector<float> jetf_truth_E;
            std::vector<float> jetf_truthcharged_px;
            std::vector<float> jetf_truthcharged_py;
            std::vector<float> jetf_truthcharged_pz;
            std::vector<float> jetf_truthcharged_E;
            for(Int_t imc=0; imc<_nMCPart; imc++){
                TVector3 truevec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
//                 if(truevec.Eta()<1){
                    if(verbosity>1) std::cout << "\tMC:  particle " << imc << "\twith E = " << _mcpart_E[imc] << " GeV" << std::endl;
                    //  create truth vector for jet finder
                    jetf_truth_px.push_back(_mcpart_px[imc]);
                    jetf_truth_py.push_back(_mcpart_py[imc]);
                    jetf_truth_pz.push_back(_mcpart_pz[imc]);
                    jetf_truth_E.push_back(_mcpart_E[imc]);
//                 }
                // fill charged only inputs (reject neutral particles)
//                 std::cout << _mcpart_PDG[imc] << std::endl;
                if(abs(_mcpart_PDG[imc])==211 || abs(_mcpart_PDG[imc])==321 || abs(_mcpart_PDG[imc])==2212 || abs(_mcpart_PDG[imc])==11 || abs(_mcpart_PDG[imc])==13){
                    jetf_truthcharged_px.push_back(_mcpart_px[imc]);
                    jetf_truthcharged_py.push_back(_mcpart_py[imc]);
                    jetf_truthcharged_pz.push_back(_mcpart_pz[imc]);
                    jetf_truthcharged_E.push_back(_mcpart_E[imc]);
                }
            }

            std::vector<double> nocluster_px;
            std::vector<double> nocluster_py;
            std::vector<double> nocluster_pz;
            std::vector<double> nocluster_E;

            // ANCHOR JET FINDING
            // truth jets
            auto jetsTrue = findJets(jetR, jetAlgorithm, jetf_truth_px, jetf_truth_py, jetf_truth_pz, jetf_truth_E);
            if(verbosity>1) std::cout << "found " << std::get<1>(jetsTrue).size() << " true jets" << std::endl;        // printJets(std::get<1>(jetsTrue));

            auto jetsTrueCharged = findJets(jetR, jetAlgorithm, jetf_truthcharged_px, jetf_truthcharged_py, jetf_truthcharged_pz, jetf_truthcharged_E);
            if(verbosity>1) std::cout << "found " << std::get<1>(jetsTrueCharged).size() << " true charged jets" << std::endl;        // printJets(std::get<1>(jetsTrue));

            // track-based jets (rec)
            auto jetsTrackRec = findJets(jetR, jetAlgorithm, jetf_track_px, jetf_track_py, jetf_track_pz, jetf_track_E);
            if(verbosity>1) std::cout << "found " << std::get<1>(jetsTrackRec).size() << " rec track jets" << std::endl;        // printJets(std::get<1>(jetsTrackRec));

            // full jets (rec)
            auto jetsFullRec = findJets(jetR, jetAlgorithm, jetf_full_px, jetf_full_py, jetf_full_pz, jetf_full_E);
            if(verbosity>1) std::cout << "found " << std::get<1>(jetsFullRec).size() << " rec full jets" << std::endl;        // printJets(std::get<1>(jetsTrackRec));

            // hcal jets (rec)
            auto jetsHcalRec = findJets(jetR, jetAlgorithm, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz, jetf_hcal_E);
            if(verbosity>1) std::cout << "found " << std::get<1>(jetsHcalRec).size() << " rec hcal jets" << std::endl;        // printJets(std::get<1>(jetsTrackRec));

            // emcal jets (rec)
              // calo jets (rec)
            auto jetsCaloRec = findJets(jetR, jetAlgorithm, jetf_calo_px, jetf_calo_py, jetf_calo_pz, jetf_calo_E);
            if(verbosity>1) std::cout << "found " << std::get<1>(jetsCaloRec).size() << " rec calo jets" << std::endl;        // printJets(std::get<1>(jetsTrackRec));

            // all jets (rec)
            auto jetsAllRec = findJets(jetR, jetAlgorithm, jetf_all_px, jetf_all_py, jetf_all_pz, jetf_all_E);
            if(verbosity>1) std::cout << "found " << std::get<1>(jetsAllRec).size() << " rec all jets" << std::endl;        // printJets(std::get<1>(jetsTrackRec));

            // calo jets, no clustering
            auto jetsNoCluster = findJets(jetR, jetAlgorithm, nocluster_px, nocluster_py, nocluster_pz, nocluster_E);
            if(verbosity>1) std::cout << "found " << std::get<1>(jetsNoCluster).size() << " rec calo jets" << std::endl;        // printJets(std::get<1>(jetsTrackRec));


            // Jet observables
            fillEventObservables(eventObservables, primaryTrackSource);
            fillJetObservables(jetObservablesTrue, std::get<1>(jetsTrue), jetR);
            fillJetObservables(jetObservablesTrueCharged, std::get<1>(jetsTrueCharged), jetR);
            fillJetObservables(jetObservablesCharged, std::get<1>(jetsTrackRec), jetR);
            fillJetObservables(jetObservablesCalo, std::get<1>(jetsCaloRec), jetR);
            fillJetObservables(jetObservablesFull, std::get<1>(jetsFullRec), jetR);

            fillHadronObservables(jetObservablesTrue);

            jetresolutionhistos(jetsTrackRec, jetsTrueCharged, 0, jetR);
            jetresolutionhistos(jetsFullRec, jetsTrue, 1, jetR);
            jetresolutionhistos(jetsHcalRec, jetsTrue, 2, jetR);
            jetresolutionhistos(jetsCaloRec, jetsTrue, 3, jetR);
            jetresolutionhistos(jetsAllRec, jetsTrue, 4, jetR);
            jetresolutionhistos(jetsNoCluster,  jetsTrue,  5, jetR);
            // TString jettype[njettypes] = {"track", "full","hcal","calo","all"};
        }
        if(tracksEnabled){
          if(verbosity>1) std::cout << "running trackingefficiency" << std::endl;
          trackingefficiency();
          if(verbosity>1) std::cout << "running trackingresolution" << std::endl;
          trackingresolution();
          if(verbosity>1) std::cout << "running trackingcomparison" << std::endl;
          trackingcomparison();
          if(verbosity>1) std::cout << "finished tracking studies" << std::endl;
        }
        if (runPi0Reco){
          pi0studies();
        }
        if (runCaloRes){
          if(verbosity>1) std::cout << "running clusterstudies" << std::endl;
          clusterstudies_sphenix();
          if(verbosity>1) std::cout << "running  caloresolutionhistos" << std::endl;
          caloresolutionhistos_sphenix();
          if(verbosity>1) std::cout << "loop done ... next event" << std::endl;
        }
        if(tracksEnabled) trackmatchingstudies_sphenix();

        clearClusterVectors();
        clearMCRecMatchVectors();

    } // event loop end
    if (_do_jetfinding) {
        std::cout << "saving event level observables\n";
        eventObservables.Write(outputDir.Data());
        std::cout << "saving jet observables\n";
        jetObservablesTrue.Write(outputDir.Data(), "RECREATE");
        jetObservablesTrueCharged.Write(outputDir.Data());
        jetObservablesCharged.Write(outputDir.Data());
        jetObservablesCalo.Write(outputDir.Data());
        jetObservablesFull.Write(outputDir.Data());
    }
    std::cout << "running jetresolutionhistosSave" << std::endl;
    jetresolutionhistosSave();
    if (runCaloRes){
      std::cout << "running caloresolutionhistosSave" << std::endl;
      caloresolutionhistossave();
      std::cout << "running clusterstudiesSave" << std::endl;
      clusterstudiesSave();
    }
    if(tracksEnabled){
      std::cout << "running trackingefficiencyhistosSave" << std::endl;
      trackingefficiencyhistosSave();
    }
    if(tracksEnabled) {
      std::cout << "running trackingresolutionhistosSave" << std::endl;
      trackingresolutionhistosSave();
    }
    if(tracksEnabled) {
      std::cout << "running trackingcomparisonhistosSave" << std::endl;
      trackingcomparisonhistosSave();
    }
    if(tracksEnabled){
      std::cout << "running hitstudiesSave" << std::endl;
      hitstudiesSave();
    }
    if(runPi0Reco){
      pi0studiesSave();
    }
    if(tracksEnabled) {
      std::cout << "running trackmatchingstudiesSave" << std::endl;
      trackmatchingstudiesSave();
    }
    std::cout << "running clusterizerSave" << std::endl;
    clusterizerSave();
    std::cout << "all done :)" << std::endl;
}



// MAIN FUNCTION for non-ROOT compilation (prototype)
int main( int argc, char* argv[] )
{

    // Default arguments for treeProcessing
    TString inFile              = "";
    TString inFileGeometry      = "geometry.root";
    TString addOutputName       = "";
    Double_t maxNEvent          = -1;
    bool do_reclus              = true;
    bool do_jetfinding          = false;
    bool runCaloRes             = true;
    bool doCalibration          = false;
    Int_t verbosity             = 0;
    // Defaults to tracking from all layers.
    unsigned short primaryTrackSource = 0;
    std::string jetAlgorithm    = "anti-kt";
    double jetR                 = 0.5;
    double tracked_jet_max_pT   = 30;
    bool removeTracklets        = false;

    // Import main call arguments
        TString import;
        if( argc >  1 ) inFile                 = argv[1];
        if( argc >  2 ) inFileGeometry                  = argv[2];
        if( argc >  3 ) addOutputName          = argv[3];
        if( argc >  4 ) { // numberOfBins
            std::istringstream sstr(argv[4]);
            sstr >> maxNEvent;
        }
        if( argc >  5 ) { // UseTHnSparse
            import = argv[5];
            if( import.EqualTo("kTRUE") || import.EqualTo("true") )  do_reclus = kTRUE;
            if( import.EqualTo("kFALSE") || import.EqualTo("false") )  do_reclus = kFALSE;
        }
        if( argc >  6 ) { // UseTHnSparse
            import = argv[6];
            if( import.EqualTo("kTRUE") || import.EqualTo("true") )  do_jetfinding = kTRUE;
            if( import.EqualTo("kFALSE") || import.EqualTo("false") )  do_jetfinding = kFALSE;
        }
        if( argc >  7 ) { // UseTHnSparse
            import = argv[7];
            if( import.EqualTo("kTRUE") || import.EqualTo("true") )  runCaloRes = kTRUE;
            if( import.EqualTo("kFALSE") || import.EqualTo("false") )  runCaloRes = kFALSE;
        }
        if( argc >  8 ) { // UseTHnSparse
            import = argv[8];
            if( import.EqualTo("kTRUE") || import.EqualTo("true") )  doCalibration = kTRUE;
            if( import.EqualTo("kFALSE") || import.EqualTo("false") )  doCalibration = kFALSE;
        }
        if( argc >  9 ) {
          std::istringstream sstr(argv[9]);
          sstr >> verbosity;
        }
        if( argc > 10 ) {
          std::istringstream sstr(argv[10]);
          sstr >> primaryTrackSource;
        }
        if( argc > 11 ) jetAlgorithm = argv[11];
        if( argc > 12 ) { // numberOfBins
            std::stringstream sstr(argv[12]);
            sstr >> jetR;
        } if( argc > 13 ) { // numberOfBins
            std::stringstream sstr(argv[13]);
            sstr >> tracked_jet_max_pT;
        } if( argc > 14 ) { // UseTHnSparse
            import = argv[14];
            if( import.EqualTo("kTRUE") || import.EqualTo("true") )  removeTracklets = kTRUE;
            if( import.EqualTo("kFALSE") || import.EqualTo("false") )  removeTracklets = kFALSE;
        }

    // Function call ExtractSignalV2
        treeProcessing_sphenix(
            inFile,
            inFileGeometry,
            addOutputName,
            maxNEvent,
            do_reclus,
            do_jetfinding,
            runCaloRes,
            doCalibration,
            verbosity,
            primaryTrackSource,
            jetAlgorithm,
            jetR,
            tracked_jet_max_pT,
            removeTracklets
        );
    return 0;


}
