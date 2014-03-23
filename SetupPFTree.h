//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 21 16:49:28 2014 by ROOT version 5.32/00
// from TTree pfTree/dijet tree
// found on file: /mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet21_STARTHI53_LV1/merged3/HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0.root
//////////////////////////////////////////////////////////
#include "commonSetup.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

class PFs {
public :
   PFs(){};
   ~PFs(){};

   // Declaration of leaf types
   Int_t           nPFpart;
   Int_t           pfId[14010];   //[nPFpart]
   Float_t         pfPt[14010];   //[nPFpart]
   Float_t         pfVsPt[14010];   //[nPFpart]
   Float_t         pfVsPtInitial[14010];   //[nPFpart]
   Float_t         pfArea[14010];   //[nPFpart]
   Float_t         pfEta[14010];   //[nPFpart]
   Float_t         pfPhi[14010];   //[nPFpart]
   Float_t         vn[5][15];
   Float_t         psin[5][15];
   Float_t         sumpt[15];

   // List of branches
   TBranch        *b_nPFpart;   //!
   TBranch        *b_pfId;   //!
   TBranch        *b_pfPt;   //!
   TBranch        *b_pfVsPt;   //!
   TBranch        *b_pfVsPtInitial;   //!
   TBranch        *b_pfArea;   //!
   TBranch        *b_pfEta;   //!
   TBranch        *b_pfPhi;   //!
   TBranch        *b_vn;   //!
   TBranch        *b_vpsi;   //!
   TBranch        *b_sumpt;   //!

};


void setupPFTree(TTree *t,PFs &tPFs,bool doCheck = 0)
{
   // Set branch addresses and branch pointers
   t->SetBranchAddress("nPFpart", &tPFs.nPFpart, &tPFs.b_nPFpart);
   t->SetBranchAddress("pfId", tPFs.pfId, &tPFs.b_pfId);
   t->SetBranchAddress("pfPt", tPFs.pfPt, &tPFs.b_pfPt);
   t->SetBranchAddress("pfVsPt", tPFs.pfVsPt, &tPFs.b_pfVsPt);
   t->SetBranchAddress("pfVsPtInitial", tPFs.pfVsPtInitial, &tPFs.b_pfVsPtInitial);
   t->SetBranchAddress("pfArea", tPFs.pfArea, &tPFs.b_pfArea);
   t->SetBranchAddress("pfEta", tPFs.pfEta, &tPFs.b_pfEta);
   t->SetBranchAddress("pfPhi", tPFs.pfPhi, &tPFs.b_pfPhi);
   t->SetBranchAddress("vn", tPFs.vn, &tPFs.b_vn);
   t->SetBranchAddress("psin", tPFs.psin, &tPFs.b_vpsi);
   t->SetBranchAddress("sumpt", tPFs.sumpt, &tPFs.b_sumpt);
   if (doCheck) {
      if (t->GetMaximum("nPFpart")>14010) cout <<"FATAL ERROR: Arrary size of nPFpart too small!!!  "<<t->GetMaximum("nPFpart")<<endl;
   }
}

