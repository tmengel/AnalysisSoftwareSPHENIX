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

#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TChain.h>
#include <TVector3.h>

#include <iostream>
#include <fstream>

void treeProcessing_simple_sphenix(
    TString inFile              = "",                 //
    TString inFileGeometry      = "geometry.root",    //
    TString addOutputName       = "",                 //
    Double_t maxNEvent          = -1,                 //
    bool do_reclus              = true,               //
    bool doCalibration          = false,              //
    Int_t verbosity             = 0,                  //
    // Defaults to tracking from all layers.
    unsigned short primaryTrackSource = 0             //
){
    // make output directory
    TString dateForOutput = ReturnDateStr();
    outputDir = Form("treeProcessingSimple/%s",addOutputName.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    // load tree
    TChain *const tt_event = new TChain("EventTree");
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
    tt_geometry =  (TTree *) (new TFile(inFileGeometry.Data(), "READ"))->Get("GeoTree");
    if(!tt_geometry){ cout << "geometry tree not found... returning!"<< endl; return;}
    // load all branches (see header)
    SetBranchAddressesTree(tt_event);
    SetBranchAddressesGeometryTree(tt_geometry);
    SetGeometryIndices();

    for (Int_t c = 0; c < 12; c++){
      cout << str_calorimeter[c] << "\t" << caloEnabled[c] << endl;
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

        Int_t nTowers[maxcalo] = {_nTowers_CEMC,
                                   _nTowers_HCALIN, _nTowers_HCALOUT  };

        // run clusterizers normal calos
        for (int cal = 0; cal < maxcalo; cal++){
          if(do_reclus && nTowers[cal] > 0 && caloEnabled[cal]){
            for (int algo = 0; algo < _active_algo; algo++){
              if(verbosity>1) cout << "clusterizing " << str_clusterizer[algo].Data() <<  " for " << str_calorimeter[cal].Data() << endl;
              runclusterizer(algo, cal, seedE[cal], aggE[cal], primaryTrackSource);
            }
          }
        }
        if((do_reclus) && verbosity>1) cout << "done with clusterization!" << endl;


        // Reset array entries for true ID
        for(Int_t itrk=0; itrk<m_nTracks; itrk++){
          m_truthtrackID[itrk] = GetCorrectMCArrayEntry(m_truthtrackID[itrk]);
          if(verbosity>2) std::cout << "\tTrack: track " << itrk << "\twith true ID " << m_truthtrackID[itrk] << "\tand X = " << m_tr_px[itrk] << " cm" << std::endl;
        }
        // obtain labels for different track sources
        prepareMCMatchInfo();


        for (Int_t icalo = 0; icalo < maxcalo; icalo++){
          for (Int_t ialgo = 0; ialgo < maxAlgo; ialgo++){
            if (loadClusterizerInput(ialgo, icalo) && verbosity>2) std::cout << str_calorimeter[icalo].Data() << "\t" << ialgo << endl;
          }
        }


        // simple TM validation
        if(tracksEnabled) trackmatchingstudies_sphenix();
        // simple cluster output QA
        clusterstudies_sphenix();

        // *****************************************************************************************************************
        // *****************************************************************************************************************
        // **********************   PUT your functions HERE ****************************************************************
        // *****************************************************************************************************************
        // *****************************************************************************************************************


        // cleanup of variables
        clearClusterVectors();
        clearMCRecMatchVectors();

    } // event loop end

    // Saving histos to output files
    if(tracksEnabled) {
      std::cout << "running trackmatchingstudiesSave" << std::endl;
      trackmatchingstudiesSave();
    }
    std::cout << "running clusterstudiesSave" << std::endl;
    clusterstudiesSave();

    // *****************************************************************************************************************
    // *****************************************************************************************************************
    // **********************   PUT your save functions HERE ***********************************************************
    // *****************************************************************************************************************
    // *****************************************************************************************************************


    std::cout << "running clusterizerSave" << std::endl;
    clusterizerSave();
    std::cout << "all done :)" << std::endl;
}
