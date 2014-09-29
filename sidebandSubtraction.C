#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <iostream>
#include <TH2D.h>
#include <TCanvas.h>

TH2D *sidebandSubtraction(TH2D *hInput, char *outputName="hTemp")
{
   
   int nBinX = hInput->GetNbinsX();
   int nBinY = hInput->GetNbinsY();
   
   TH2D *hOutput = (TH2D*)hInput->Clone(outputName);

   return hOutput;
}

void sidebandSubtraction()
{
   TFile *_file1 = new TFile("Data2011_Mixed_NewEff_trkPtCut1_Sept25.root");
   TH2D *hallieNew = (TH2D*)_file1->Get("hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2");
   TH2D *hallieMix = (TH2D*)_file1->Get("hists_hJetTrackMECent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2");
 
   TCanvas *cRatio = new TCanvas("cRatio","signal / ME", 600,600);
   TH2D *hRatio = hallieNew->Clone("hRatio");
   hRatio->Draw("surf1");

   TCanvas *cProjectionX = new TCanvas("cProjectionX","signal / ME", 600,600);
   TH2D *h = sidebandSubtraction(hallieNew);  
}
