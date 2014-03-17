#include "HiForestAnalysis/hiForest.h"
#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <fstream>
#include <iostream>
#include <math.h>


void randomcone(string filename, int centmin = 0, int centmax = 4);
float getdR(float eta1, float eta2, float phi1, float phi2);

void randomcone(int condor_iter, string flist, string tag, int centmin, int centmax)
{
  using namespace std;
  
  string buffer;
  vector<string> listoffiles;
  int nlines = 0;
  ifstream infile(flist.data());
  if (!infile.is_open()) {
    cout << "Error opening file. Exiting." << endl;
    return;
  } else {
    while (!infile.eof()) {
      infile >> buffer;
      listoffiles.push_back(buffer);
      nlines++;
    }
  }
  
  HiForest * c = new HiForest(listoffiles[condor_iter].data(),"forest",cPPb,0);
  
  // TH2D * etaphi = new TH2D();
  int nevents = c->GetEntries();
  
  TRandom3 r;
  float rConeEta = 0;
  float rConePhi = 0;
  float dR = 0.3;
  float thisdR = 0;
  float pfPtSum;
  float pfVsPtSum;
  float pfVsPtInitialSum;
  
  TFile * outf = new TFile(Form("randcone_%s_%d_centmin%d_centmax%d.root",tag.data(),condor_iter,centmin,centmax),"recreate");
  TH1D * hpfPtSum = new TH1D("hpfPtSum","pfPtSum;#Sigma p_{T}",500,0,500);
  TH1D * hpfVsPtSum = new TH1D("hpfVsPtSum","pfVsPtSum;#Sigma p_{T}",500,0,50);
  TH1D * hpfVsPtInitialSum = new TH1D("hpfVsPtInitialSum","hpfVsPtInitialSum;#Sigma p_{T}",500,0,50);
  for(int i = 0 ; i < nevents ; ++i)
  {
    if(i%1000==0) cout<<i<<"/"<<nevents<<endl;
    c->GetEntry(i);
    if(c->evt.hiBin < centmin || c->evt.hiBin > centmax) continue;
    rConeEta = 2.4*(r.Rndm()-0.5); // uniform dist. +- 2.4 eta
    rConePhi = TMath::Pi()*(r.Rndm()-0.5); // uniform dist. +- pi
    pfPtSum = 0.0;
    pfVsPtSum = 0.0;
    pfVsPtInitialSum = 0.0;
    for(int j = 0 ; j < c->pf.nPFpart ; ++j)
    {
      thisdR = getdR(rConeEta,c->pf.pfEta[j],rConePhi,c->pf.pfPhi[j]);
      if(thisdR < dR)
      {
        pfPtSum+=c->pf.pfPt[j];
        pfVsPtSum+=c->pf.pfVsPt[j];
        pfVsPtInitialSum+=c->pf.pfVsPtInitial[j];
      }
    }
    hpfPtSum->Fill(pfPtSum);
    hpfVsPtSum->Fill(pfVsPtSum);
    hpfVsPtInitialSum->Fill(pfVsPtInitialSum);
  }
  // TCanvas c1;
  // hpfPtSum->Draw();
  // TCanvas c2;
  // hpfVsPtSum->Draw();
  // TCanvas c3;
  // hpfVsPtInitialSum->Draw();
  outf->Write();
  outf->Close();
}

float getdR(float eta1, float eta2, float phi1, float phi2)
{
  float detasq = (eta1-eta2)*(eta1-eta2);
  float dphi = fabs(phi1 - phi2);
  if(dphi > TMath::Pi())
    dphi = 2*TMath::Pi() - dphi;
  float dphisq = dphi*dphi;
  return sqrt( dphisq + detasq );
}

int main(int argc, char *argv[])
{
  if(argc != 6)
  {
    cout<<"usage ./randomcone.exe <filenum> <tag> <filelist> <centmin> <centmax>"<<endl;
    return 1;
  }
  randomcone(atoi(argv[1]),argv[2],argv[3],atof(argv[4]),atof(argv[5]));
  return 0;
}
