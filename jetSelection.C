#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TCut.h>

void jetSelection(char *infname = "skim_JetData.root",int weight = 0)
{
   TFile *inf = new TFile(infname);
   TTree *skimTree = (TTree*) inf->Get("skimanalysis/HltTree");
   TTree *evtTree  = (TTree*) inf->Get("hiEvtAnalyzer/HiTree");
   TTree *jetTree  = (TTree*) inf->Get("akVs3CaloJetAnalyzer/t");
   TTree *trackTree  = (TTree*) inf->Get("ppTrack/trackTree");
   if (trackTree==0) trackTree = (TTree*) inf->Get("anaTrack/trackTree");
   TCut centWeight = "(4.067+hiBin*(-0.06162)+hiBin*hiBin*(0.0003076)+hiBin*hiBin*hiBin*(-5.032e-7))";
   if (weight==0) centWeight="((pcollisionEventSelection&&pHBHENoiseFilter))";
   jetTree->AddFriend(trackTree);
   jetTree->AddFriend(skimTree);
   jetTree->AddFriend(evtTree);

   TCanvas *c = new TCanvas("c","Number of pixel hit",600,600);
   TH2D *h = new TH2D("h","",50,0,50,100,0,60000);
   h->SetXTitle("Number of Jets (|#eta|<2 && p_{T}>50GeV)");
   h->SetYTitle("Number of pixel hits (hiNpix)");
   h->SetTitleOffset(1.5,"Y");
   jetTree->Draw("hiNpix:Sum$(jtpt>50&&abs(jteta)<2)>>h",centWeight,"col");

   TCanvas *c2 = new TCanvas("c2","Number of pixel hit",600,600);
   TH2D *h2 = new TH2D("h2","",50,0,50,100,0,1000);
   h2->SetXTitle("Number of Jets (|#eta|<2 && p_{T}>50GeV)");
   h2->SetYTitle("Number of large dz tracks");
   h2->SetTitleOffset(1.5,"Y");

   jetTree->Draw("Sum$(highPurity&&abs(trkDxy1/trkDxyError1)<3&&abs(trkDz1/trkDzError1)>2&&trkPtError/trkPt<0.1):Sum$(jtpt>50&&abs(jteta)<2)>>h2",centWeight,"col");


   TCanvas *c3 = new TCanvas("c3","Number of pixel hit",600,600);
   TH2D *h3 = new TH2D("h3","",50,0,50,100,0,1);
   h3->SetXTitle("Number of Jets (|#eta|<2 && p_{T}>50GeV)");
   h3->SetYTitle("Fraction of large dz tracks");
   h3->SetTitleOffset(1.5,"Y");

   jetTree->Draw("Sum$(highPurity&&abs(trkDxy1/trkDxyError1)<3&&abs(trkDz1/trkDzError1)>2&&trkPtError/trkPt<0.1)/Sum$(highPurity&&abs(trkDxy1/trkDxyError1)<2&&abs(trkDz1/trkDzError1)<3&&trkPtError/trkPt<0.1):Sum$(jtpt>50&&abs(jteta)<2)>>h3",centWeight,"col");

//   TCanvas *c = new TCanvas("c","",600,600)
   
}
