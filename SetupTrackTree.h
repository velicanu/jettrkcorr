//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 21 16:49:25 2014 by ROOT version 5.32/00
// from TTree trackTree/v1
// found on file: /mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet21_STARTHI53_LV1/merged3/HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0.root
//////////////////////////////////////////////////////////
#include "commonSetup.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

class Tracks {
public :
   Tracks(){};
   ~Tracks(){};

   // Declaration of leaf types
   Int_t           nEv;
   Int_t           nLumi;
   Int_t           nBX;
   Int_t           nRun;
   Int_t           N;
   Int_t           nv;
   Float_t         vx[4];   //[nv]
   Float_t         vy[4];   //[nv]
   Float_t         vz[4];   //[nv]
   Float_t         vxErr[4];   //[nv]
   Float_t         vyErr[4];   //[nv]
   Float_t         vzErr[4];   //[nv]
   Int_t           nDaugher[4];   //[nv]
   Int_t           nVtx;
   Int_t           maxVtx;
   Int_t           nTrkVtx[2];   //[nVtx]
   Float_t         normChi2Vtx[2];   //[nVtx]
   Float_t         sumPtVtx[2];   //[nVtx]
   Float_t         xVtx[2];   //[nVtx]
   Float_t         yVtx[2];   //[nVtx]
   Float_t         zVtx[2];   //[nVtx]
   Float_t         xVtxErr[2];   //[nVtx]
   Float_t         yVtxErr[2];   //[nVtx]
   Float_t         zVtxErr[2];   //[nVtx]
   Float_t         vtxDist2D[2];   //[nVtx]
   Float_t         vtxDist2DErr[2];   //[nVtx]
   Float_t         vtxDist2DSig[2];   //[nVtx]
   Float_t         vtxDist3DErr[2];   //[nVtx]
   Float_t         vtxDist3DSig[2];   //[nVtx]
   Int_t           nTrk;
   Float_t         trkPt[12284];   //[nTrk]
   Float_t         trkPtError[12284];   //[nTrk]
   Int_t           trkNHit[12284];   //[nTrk]
   Int_t           trkNlayer[12284];   //[nTrk]
   Float_t         trkEta[12284];   //[nTrk]
   Float_t         trkPhi[12284];   //[nTrk]
   Int_t           trkCharge[12284];   //[nTrk]
   Int_t           trkVtxIndex[12284];   //[nTrk]
   Bool_t          highPurity[12284];   //[nTrk]
   Bool_t          highPuritySetWithPV[12284];   //[nTrk]
   Float_t         trkChi2[12284];   //[nTrk]
   Float_t         trkNdof[12284];   //[nTrk]
   Float_t         trkDxy1[12284];   //[nTrk]
   Float_t         trkDxyError1[12284];   //[nTrk]
   Float_t         trkDz1[12284];   //[nTrk]
   Float_t         trkDzError1[12284];   //[nTrk]
   Bool_t          trkFake[12284];   //[nTrk]
   Float_t         trkAlgo[12284];   //[nTrk]
   Float_t         trkStatus[12284];   //[nTrk]
   Float_t         trkPId[12284];   //[nTrk]
   Float_t         trkMPId[12284];   //[nTrk]
   Float_t         trkGMPId[12284];   //[nTrk]
   Int_t           nParticle;
   Float_t         pStatus[13496];   //[nParticle]
   Float_t         pPId[13496];   //[nParticle]
   Float_t         pEta[13496];   //[nParticle]
   Float_t         pPhi[13496];   //[nParticle]
   Float_t         pPt[13496];   //[nParticle]
   Float_t         pAcc[13496];   //[nParticle]
   Float_t         pAccPair[13496];   //[nParticle]
   Float_t         pNRec[13496];   //[nParticle]
   Int_t           pNHit[13496];   //[nParticle]
   Float_t         mtrkPt[13496];   //[nParticle]
   Float_t         mtrkPtError[13496];   //[nParticle]
   Int_t           mtrkNHit[13496];   //[nParticle]
   Int_t           mtrkNlayer[13496];   //[nParticle]
   Int_t           mtrkNlayer3D[13496];   //[nParticle]
   Int_t           mtrkQual[13496];   //[nParticle]
   Float_t         mtrkChi2[13496];   //[nParticle]
   Float_t         mtrkNdof[13496];   //[nParticle]
   Float_t         mtrkDz1[13496];   //[nParticle]
   Float_t         mtrkDzError1[13496];   //[nParticle]
   Float_t         mtrkDxy1[13496];   //[nParticle]
   Float_t         mtrkDxyError1[13496];   //[nParticle]
   Float_t         mtrkAlgo[13496];   //[nParticle]

   // List of branches
   TBranch        *b_nEv;   //!
   TBranch        *b_nLumi;   //!
   TBranch        *b_nBX;   //!
   TBranch        *b_nRun;   //!
   TBranch        *b_N;   //!
   TBranch        *b_nv;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_vxErr;   //!
   TBranch        *b_vyErr;   //!
   TBranch        *b_vzErr;   //!
   TBranch        *b_nDaugher;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_maxVtx;   //!
   TBranch        *b_nTrkVtx;   //!
   TBranch        *b_normChi2Vtx;   //!
   TBranch        *b_sumPtVtx;   //!
   TBranch        *b_xVtx;   //!
   TBranch        *b_yVtx;   //!
   TBranch        *b_zVtx;   //!
   TBranch        *b_xVtxErr;   //!
   TBranch        *b_yVtxErr;   //!
   TBranch        *b_zVtxErr;   //!
   TBranch        *b_vtxDist2D;   //!
   TBranch        *b_vtxDist2DErr;   //!
   TBranch        *b_vtxDist2DSig;   //!
   TBranch        *b_vtxDist3DErr;   //!
   TBranch        *b_vtxDist3DSig;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkPtError;   //!
   TBranch        *b_trkNHit;   //!
   TBranch        *b_trkNlayer;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_trkCharge;   //!
   TBranch        *b_trkVtxIndex;   //!
   TBranch        *b_highPurity;   //!
   TBranch        *b_highPuritySetWithPV;   //!
   TBranch        *b_trkChi2;   //!
   TBranch        *b_trkNdof;   //!
   TBranch        *b_trkDxy1;   //!
   TBranch        *b_trkDxyError1;   //!
   TBranch        *b_trkDz1;   //!
   TBranch        *b_trkDzError1;   //!
   TBranch        *b_trkFake;   //!
   TBranch        *b_trkAlgo;   //!
   TBranch        *b_trkStatus;   //!
   TBranch        *b_trkPId;   //!
   TBranch        *b_trkMPId;   //!
   TBranch        *b_trkGMPId;   //!
   TBranch        *b_nParticle;   //!
   TBranch        *b_pStatus;   //!
   TBranch        *b_pPId;   //!
   TBranch        *b_pEta;   //!
   TBranch        *b_pPhi;   //!
   TBranch        *b_pPt;   //!
   TBranch        *b_pAcc;   //!
   TBranch        *b_pAccPair;   //!
   TBranch        *b_pNRec;   //!
   TBranch        *b_pNHit;   //!
   TBranch        *b_mtrkPt;   //!
   TBranch        *b_mtrkPtError;   //!
   TBranch        *b_mtrkNHit;   //!
   TBranch        *b_mtrkNlayer;   //!
   TBranch        *b_mtrkNlayer3D;   //!
   TBranch        *b_mtrkQual;   //!
   TBranch        *b_mtrkChi2;   //!
   TBranch        *b_mtrkNdof;   //!
   TBranch        *b_mtrkDz1;   //!
   TBranch        *b_mtrkDzError1;   //!
   TBranch        *b_mtrkDxy1;   //!
   TBranch        *b_mtrkDxyError1;   //!
   TBranch        *b_mtrkAlgo;   //!

};


void setupTrackTree(TTree *t,Tracks &tTracks,bool doCheck = 0)
{
   // Set branch addresses and branch pointers
   t->SetBranchAddress("nEv", &tTracks.nEv, &tTracks.b_nEv);
   t->SetBranchAddress("nLumi", &tTracks.nLumi, &tTracks.b_nLumi);
   t->SetBranchAddress("nBX", &tTracks.nBX, &tTracks.b_nBX);
   t->SetBranchAddress("nRun", &tTracks.nRun, &tTracks.b_nRun);
   t->SetBranchAddress("N", &tTracks.N, &tTracks.b_N);
   t->SetBranchAddress("nv", &tTracks.nv, &tTracks.b_nv);
   t->SetBranchAddress("vx", tTracks.vx, &tTracks.b_vx);
   t->SetBranchAddress("vy", tTracks.vy, &tTracks.b_vy);
   t->SetBranchAddress("vz", tTracks.vz, &tTracks.b_vz);
   t->SetBranchAddress("vxErr", tTracks.vxErr, &tTracks.b_vxErr);
   t->SetBranchAddress("vyErr", tTracks.vyErr, &tTracks.b_vyErr);
   t->SetBranchAddress("vzErr", tTracks.vzErr, &tTracks.b_vzErr);
   t->SetBranchAddress("nDaugher", tTracks.nDaugher, &tTracks.b_nDaugher);
   t->SetBranchAddress("nVtx", &tTracks.nVtx, &tTracks.b_nVtx);
   t->SetBranchAddress("maxVtx", &tTracks.maxVtx, &tTracks.b_maxVtx);
   t->SetBranchAddress("nTrkVtx", tTracks.nTrkVtx, &tTracks.b_nTrkVtx);
   t->SetBranchAddress("normChi2Vtx", tTracks.normChi2Vtx, &tTracks.b_normChi2Vtx);
   t->SetBranchAddress("sumPtVtx", tTracks.sumPtVtx, &tTracks.b_sumPtVtx);
   t->SetBranchAddress("xVtx", tTracks.xVtx, &tTracks.b_xVtx);
   t->SetBranchAddress("yVtx", tTracks.yVtx, &tTracks.b_yVtx);
   t->SetBranchAddress("zVtx", tTracks.zVtx, &tTracks.b_zVtx);
   t->SetBranchAddress("xVtxErr", tTracks.xVtxErr, &tTracks.b_xVtxErr);
   t->SetBranchAddress("yVtxErr", tTracks.yVtxErr, &tTracks.b_yVtxErr);
   t->SetBranchAddress("zVtxErr", tTracks.zVtxErr, &tTracks.b_zVtxErr);
   t->SetBranchAddress("vtxDist2D", tTracks.vtxDist2D, &tTracks.b_vtxDist2D);
   t->SetBranchAddress("vtxDist2DErr", tTracks.vtxDist2DErr, &tTracks.b_vtxDist2DErr);
   t->SetBranchAddress("vtxDist2DSig", tTracks.vtxDist2DSig, &tTracks.b_vtxDist2DSig);
   t->SetBranchAddress("vtxDist2D", tTracks.vtxDist2D, &tTracks.b_vtxDist2D);
   t->SetBranchAddress("vtxDist3DErr", tTracks.vtxDist3DErr, &tTracks.b_vtxDist3DErr);
   t->SetBranchAddress("vtxDist3DSig", tTracks.vtxDist3DSig, &tTracks.b_vtxDist3DSig);
   t->SetBranchAddress("nTrk", &tTracks.nTrk, &tTracks.b_nTrk);
   t->SetBranchAddress("trkPt", tTracks.trkPt, &tTracks.b_trkPt);
   t->SetBranchAddress("trkPtError", tTracks.trkPtError, &tTracks.b_trkPtError);
   t->SetBranchAddress("trkNHit", tTracks.trkNHit, &tTracks.b_trkNHit);
   t->SetBranchAddress("trkNlayer", tTracks.trkNlayer, &tTracks.b_trkNlayer);
   t->SetBranchAddress("trkEta", tTracks.trkEta, &tTracks.b_trkEta);
   t->SetBranchAddress("trkPhi", tTracks.trkPhi, &tTracks.b_trkPhi);
   t->SetBranchAddress("trkCharge", tTracks.trkCharge, &tTracks.b_trkCharge);
   t->SetBranchAddress("trkVtxIndex", tTracks.trkVtxIndex, &tTracks.b_trkVtxIndex);
   t->SetBranchAddress("highPurity", tTracks.highPurity, &tTracks.b_highPurity);
   t->SetBranchAddress("highPuritySetWithPV", tTracks.highPuritySetWithPV, &tTracks.b_highPuritySetWithPV);
   t->SetBranchAddress("trkChi2", tTracks.trkChi2, &tTracks.b_trkChi2);
   t->SetBranchAddress("trkNdof", tTracks.trkNdof, &tTracks.b_trkNdof);
   t->SetBranchAddress("trkDxy1", tTracks.trkDxy1, &tTracks.b_trkDxy1);
   t->SetBranchAddress("trkDxyError1", tTracks.trkDxyError1, &tTracks.b_trkDxyError1);
   t->SetBranchAddress("trkDz1", tTracks.trkDz1, &tTracks.b_trkDz1);
   t->SetBranchAddress("trkDzError1", tTracks.trkDzError1, &tTracks.b_trkDzError1);
   t->SetBranchAddress("trkFake", tTracks.trkFake, &tTracks.b_trkFake);
   t->SetBranchAddress("trkAlgo", tTracks.trkAlgo, &tTracks.b_trkAlgo);
   t->SetBranchAddress("trkStatus", tTracks.trkStatus, &tTracks.b_trkStatus);
   t->SetBranchAddress("trkPId", tTracks.trkPId, &tTracks.b_trkPId);
   t->SetBranchAddress("trkMPId", tTracks.trkMPId, &tTracks.b_trkMPId);
   t->SetBranchAddress("trkGMPId", tTracks.trkGMPId, &tTracks.b_trkGMPId);
   t->SetBranchAddress("nParticle", &tTracks.nParticle, &tTracks.b_nParticle);
   t->SetBranchAddress("pStatus", tTracks.pStatus, &tTracks.b_pStatus);
   t->SetBranchAddress("pPId", tTracks.pPId, &tTracks.b_pPId);
   t->SetBranchAddress("pEta", tTracks.pEta, &tTracks.b_pEta);
   t->SetBranchAddress("pPhi", tTracks.pPhi, &tTracks.b_pPhi);
   t->SetBranchAddress("pPt", tTracks.pPt, &tTracks.b_pPt);
   t->SetBranchAddress("pAcc", tTracks.pAcc, &tTracks.b_pAcc);
   t->SetBranchAddress("pAccPair", tTracks.pAccPair, &tTracks.b_pAccPair);
   t->SetBranchAddress("pNRec", tTracks.pNRec, &tTracks.b_pNRec);
   t->SetBranchAddress("pNHit", tTracks.pNHit, &tTracks.b_pNHit);
   t->SetBranchAddress("mtrkPt", tTracks.mtrkPt, &tTracks.b_mtrkPt);
   t->SetBranchAddress("mtrkPtError", tTracks.mtrkPtError, &tTracks.b_mtrkPtError);
   t->SetBranchAddress("mtrkNHit", tTracks.mtrkNHit, &tTracks.b_mtrkNHit);
   t->SetBranchAddress("mtrkNlayer", tTracks.mtrkNlayer, &tTracks.b_mtrkNlayer);
   t->SetBranchAddress("mtrkNlayer3D", tTracks.mtrkNlayer3D, &tTracks.b_mtrkNlayer3D);
   t->SetBranchAddress("mtrkQual", tTracks.mtrkQual, &tTracks.b_mtrkQual);
   t->SetBranchAddress("mtrkChi2", tTracks.mtrkChi2, &tTracks.b_mtrkChi2);
   t->SetBranchAddress("mtrkNdof", tTracks.mtrkNdof, &tTracks.b_mtrkNdof);
   t->SetBranchAddress("mtrkDz1", tTracks.mtrkDz1, &tTracks.b_mtrkDz1);
   t->SetBranchAddress("mtrkDzError1", tTracks.mtrkDzError1, &tTracks.b_mtrkDzError1);
   t->SetBranchAddress("mtrkDxy1", tTracks.mtrkDxy1, &tTracks.b_mtrkDxy1);
   t->SetBranchAddress("mtrkDxyError1", tTracks.mtrkDxyError1, &tTracks.b_mtrkDxyError1);
   t->SetBranchAddress("mtrkAlgo", tTracks.mtrkAlgo, &tTracks.b_mtrkAlgo);
   if (doCheck) {
      if (t->GetMaximum("nv")>4) cout <<"FATAL ERROR: Arrary size of nv too small!!!  "<<t->GetMaximum("nv")<<endl;
      if (t->GetMaximum("nVtx")>2) cout <<"FATAL ERROR: Arrary size of nVtx too small!!!  "<<t->GetMaximum("nVtx")<<endl;
      if (t->GetMaximum("nTrk")>12284) cout <<"FATAL ERROR: Arrary size of nTrk too small!!!  "<<t->GetMaximum("nTrk")<<endl;
      if (t->GetMaximum("nParticle")>13496) cout <<"FATAL ERROR: Arrary size of nParticle too small!!!  "<<t->GetMaximum("nParticle")<<endl;
   }
}

