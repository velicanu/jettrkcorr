//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 21 16:49:20 2014 by ROOT version 5.32/00
// from TTree HltTree/
// found on file: /mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet21_STARTHI53_LV1/merged3/HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0.root
//////////////////////////////////////////////////////////
#include "commonSetup.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

class Skims {
public :
   Skims(){};
   ~Skims(){};

   // Declaration of leaf types
   Int_t           photonStep;
   Int_t           extrapatstep;
   Int_t           ana_step;
   Int_t           phltJetHI;
   Int_t           pcollisionEventSelection;
   Int_t           pHBHENoiseFilter;
   Int_t           phfCoincFilter;
   Int_t           phfCoincFilter3;
   Int_t           pprimaryVertexFilter;
   Int_t           phltPixelClusterShapeFilter;
   Int_t           phiEcalRecHitSpikeFilter;
   Int_t           hltAna;
   Int_t           ppurityFractionFilter;

   // List of branches
   TBranch        *b_photonStep;   //!
   TBranch        *b_extrapatstep;   //!
   TBranch        *b_ana_step;   //!
   TBranch        *b_phltJetHI;   //!
   TBranch        *b_pcollisionEventSelection;   //!
   TBranch        *b_pHBHENoiseFilter;   //!
   TBranch        *b_phfCoincFilter;   //!
   TBranch        *b_phfCoincFilter3;   //!
   TBranch        *b_pprimaryVertexFilter;   //!
   TBranch        *b_phltPixelClusterShapeFilter;   //!
   TBranch        *b_phiEcalRecHitSpikeFilter;   //!
   TBranch        *b_hltAna;   //!
   TBranch        *b_ppurityFractionFilter;   //!

};


void setupSkimTree(TTree *t,Skims &tSkims,bool doCheck = 0)
{
   // Set branch addresses and branch pointers
   t->SetBranchAddress("photonStep", &tSkims.photonStep, &tSkims.b_photonStep);
   t->SetBranchAddress("extrapatstep", &tSkims.extrapatstep, &tSkims.b_extrapatstep);
   t->SetBranchAddress("ana_step", &tSkims.ana_step, &tSkims.b_ana_step);
   t->SetBranchAddress("phltJetHI", &tSkims.phltJetHI, &tSkims.b_phltJetHI);
   t->SetBranchAddress("pcollisionEventSelection", &tSkims.pcollisionEventSelection, &tSkims.b_pcollisionEventSelection);
   t->SetBranchAddress("pHBHENoiseFilter", &tSkims.pHBHENoiseFilter, &tSkims.b_pHBHENoiseFilter);
   t->SetBranchAddress("phfCoincFilter", &tSkims.phfCoincFilter, &tSkims.b_phfCoincFilter);
   t->SetBranchAddress("phfCoincFilter3", &tSkims.phfCoincFilter3, &tSkims.b_phfCoincFilter3);
   t->SetBranchAddress("pprimaryVertexFilter", &tSkims.pprimaryVertexFilter, &tSkims.b_pprimaryVertexFilter);
   t->SetBranchAddress("phltPixelClusterShapeFilter", &tSkims.phltPixelClusterShapeFilter, &tSkims.b_phltPixelClusterShapeFilter);
   t->SetBranchAddress("phiEcalRecHitSpikeFilter", &tSkims.phiEcalRecHitSpikeFilter, &tSkims.b_phiEcalRecHitSpikeFilter);
   t->SetBranchAddress("hltAna", &tSkims.hltAna, &tSkims.b_hltAna);
   t->SetBranchAddress("phfCoincFilter", &tSkims.phfCoincFilter, &tSkims.b_phfCoincFilter);
   t->SetBranchAddress("ppurityFractionFilter", &tSkims.ppurityFractionFilter, &tSkims.b_ppurityFractionFilter);
   if (doCheck) {
   }
}

