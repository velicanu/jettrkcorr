#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <iostream>
#include <TH2D.h>
#include <TCanvas.h>

TH2D *sidebandSubtraction(TH2D *hInput, char *outputName="hTemp", double minEta =1.5, double maxEta = 3)
{
   
   int nBinX = hInput->GetNbinsX();
   int nBinY = hInput->GetNbinsY();
   
   TH1D *hX = (TH1D*)hInput->ProjectionX("hX");

   TH1D *hSidebandAvg = (TH1D*)hInput->ProjectionY("hY");

   TH2D *hOutput = (TH2D*)hInput->Clone(outputName);

   int negSidebandL = hX->FindBin(-maxEta+0.0000001);   
   int negSidebandH = hX->FindBin(-minEta-0.0000001);   
   int posSidebandL = hX->FindBin(minEta+0.0000001);   
   int posSidebandH = hX->FindBin(maxEta-0.0000001);   
   
   cout <<negSidebandL<<" "<<negSidebandH<<" "<<posSidebandL<<" "<<posSidebandH<<endl;

   // calculate the mean of sideband
   for (int j=1;j<=nBinY;j++){
      hSidebandAvg->SetBinContent(j,0);
      hSidebandAvg->SetBinError(j,0);
      double sumBinContent=0;
      double sumBinError=0;
      
      for (int i=negSidebandL;i<=negSidebandH;i++){
          sumBinContent=sumBinContent+hInput->GetBinContent(i,j);
          sumBinError=sumBinError+hInput->GetBinError(i,j)*hInput->GetBinError(i,j);
      }
        
      for (int i=posSidebandL;i<=posSidebandH;i++){
          sumBinContent=sumBinContent+hInput->GetBinContent(i,j);
          sumBinError=sumBinError+hInput->GetBinError(i,j)*hInput->GetBinError(i,j);
      }
      
      double nSidebandBin = negSidebandH-negSidebandL+1 +  posSidebandH-posSidebandL+1;
      cout <<sumBinContent<<endl;
      hSidebandAvg->SetBinContent(j,sumBinContent / nSidebandBin);
      hSidebandAvg->SetBinError(j,sqrt(sumBinError) / nSidebandBin);

   }

   // do subtraction
   for (int i=1;i<=nBinX;i++) {
      for (int j=1;j<=nBinY;j++) {
         double binContent = hInput->GetBinContent(i,j);
         double binError = hInput->GetBinError(i,j);
         //cout <<i<<" "<<j<<binContent<<" "<<hSidebandAvg->GetBinContent(j)<<endl;
         
         binContent-=hSidebandAvg->GetBinContent(j);
         //cout <<binContent<<" "<<hSidebandAvg->GetBinContent(j)<<endl;
         binError= sqrt(binError*binError+hSidebandAvg->GetBinError(j)*hSidebandAvg->GetBinError(j));
         hOutput->SetBinContent(i,j,binContent);
         hOutput->SetBinError(i,j,binError);
      }
   }
   return hOutput;
}

void sidebandSubtraction()
{
   // TFile *_file1 = new TFile("Data2011_Mixed_NewEff_trkPtCut1_Sept25.root");
   TFile *_file1 = new TFile("Data2011_Mixed_NewEff_trkPtCut1_Oct3.root");
   TH2D *hallieNew = (TH2D*)_file1->Get("hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2");
   TH2D *hallieMix = (TH2D*)_file1->Get("hists_hJetTrackMECent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2");
 
   TCanvas *cRatio = new TCanvas("cRatio","signal / ME", 600,600);
   TH2D *hRatio = (TH2D*)hallieNew->Clone("hRatio");
   hRatio->Divide(hallieMix);
   hRatio->SetAxisRange(-3,3,"x");
   hRatio->Draw("surf1");

   TCanvas *cSideSub = new TCanvas("cSideSub","sideband subtracted", 600,600);
   TH2D *hBckSubtracted = sidebandSubtraction(hRatio,"hBckSubtracted",1.5,3.0);  
   hBckSubtracted->Draw("surf1");
	 
	 TCanvas *cProjectionX = new TCanvas("cProjectionX","deta projection", 600,600);
   TH1D *hBckSubtractedx = hBckSubtracted->ProjectionX();  
   hBckSubtractedx->Draw();
	 
}
