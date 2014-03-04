#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TF2.h"
#include "TF1.h"
#include "TMath.h"
#include "TAttMarker.h"
#include "TH2D.h"
#include "TCanvas.h" 
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "/net/hisrv0001/home/dgulhan/run/CMSSW_4_4_4/src/UserCode/HiForest/V2/commonUtility.h"
#include "FitHist.C"

void jetplot(int trackqual = 1, string tag = "", string subleading = "", int nmin = 130, int nmax = 160, int centmin = 0, int centmax = 4, float ajmin = 0.0 , float ajmax = 1.0, float tptmin = 14, float tptmax = 16, float aptmin = 1, float aptmax = 2,int leadingjetptlow = 100 , int leadingjetpthigh = 300,  string draw = "corr", float zmax = -1 , float zmin = -1)
{
  int jetnum = 0;
  if(subleading=="") jetnum = 1;
	// HIHighPt_HIRun2011-15Apr2013-v1-HLT_HIFullTrack14_sort_trkqaul1_nmin185_nmax220_tptmin14_tptmax16_aptmin1_aptmax2.root
	// TFile * file = new TFile(Form(tag.data());
	// cout<<Form("merged/%s_trkqaul%d_nmin%d_nmax%d_tptmin%d_tptmax%d_aptmin%d_aptmax%d.root",tag.data(),trackqual,nmin,nmax,(int)tptmin,(int)tptmax,(int)aptmin,(int)aptmax)<<endl;
	// cout<<Form("merged/%s_trkqaul%d_nmin%d_nmax%d_cmin%d_cmax%d_aj%d_%d_tptmin%d_tptmax%d_aptmin%d_aptmax%d.root",tag.data(),trackqual,nmin,nmax,cmin,cmax,,(int)tptmin,(int)tptmax,(int)aptmin,(int)aptmax)<<endl;
	// system("Color 2B");
  cout<<Form("%s_trkqaul%d_nmin%d_nmax%d_cmin%d_cmax%d_aj%d_%d_tptmin%d_tptmax%d_aptmin%d_aptmax%d.root",tag.data(),trackqual,nmin,nmax,centmin,centmax,(int)(ajmin*10),(int)(ajmax*10),(int)tptmin,(int)tptmax,(int)aptmin,(int)aptmax)<<endl;
  // exit(1);
  TFile * file = new TFile(Form("merged/%s_trkqaul%d_nmin%d_nmax%d_cmin%d_cmax%d_aj%d_%d_tptmin%d_tptmax%d_aptmin%d_aptmax%d.root",tag.data(),trackqual,nmin,nmax,centmin,centmax,(int)(ajmin*10),(int)(ajmax*10),(int)tptmin,(int)tptmax,(int)aptmin,(int)aptmax));
  cout<<Form("signal_%sleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",subleading.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)aptmin,(int)aptmax,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100))<<endl;
	TH2D * Sig = (TH2D*)file->Get(Form("signal_%sleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",subleading.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)aptmin,(int)aptmax,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)));
	TH2D * Back = (TH2D*)file->Get(Form("background_%sleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",subleading.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)aptmin,(int)aptmax,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)));
	TH1D * Ntrig = (TH1D*)file->Get(Form("hjet%d_pt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetnum,leadingjetptlow,leadingjetpthigh,(int)aptmin,(int)aptmax,centmin,centmax));
	
	int x0 = Back->GetXaxis()->FindBin(0.0);
	int y0 = Back->GetYaxis()->FindBin(0.0);
	double B0 = Back->GetBinContent(x0,y0);
	// Back->Scale(1.0/B0);
	double bwx = Sig->GetXaxis()->GetBinWidth(2);
	double bwy = Sig->GetYaxis()->GetBinWidth(2);
	// Sig->Scale(1.0/(bwx*bwy*Ntrig->Integral(2,300)));
	
	TH2D * Corr = (TH2D*)Sig->Clone(Form("correlation_trg%d_%d_ass%d_%d_nmin%d_nmax%d",(int)tptmin,(int)tptmax,(int)aptmin,(int)aptmax,nmin,nmax));
	Corr->Divide(Back);
	Corr->Scale(B0/(bwx*bwy*Ntrig->Integral(2,300)));
	Corr->GetXaxis()->SetRangeUser(-4,4);
	Corr->GetYaxis()->SetRangeUser(-TMath::Pi()/2,3*TMath::Pi()/2);
  if(zmax > 0)
	{
		Corr->SetMaximum(zmax);
	}
	if(zmin > 0)
	{
		Corr->SetMinimum(zmin);
	}
	double max = Corr->GetMaximum();
	double min = Corr->GetMinimum();
	int xmn = Corr->GetXaxis()->FindBin(2);
	int xmx = Corr->GetXaxis()->FindBin(3.5);
	double nbins = xmx - xmn + 1;
	TH1D * Proj = Corr->ProjectionY(Form("projection_trg%d_%d_ass%d_%d_nmin%d_nmax%d",(int)tptmin,(int)tptmax,(int)aptmin,(int)aptmax,centmin,centmax),xmn,xmx,"e");
	Proj->Scale(1/nbins);
  Proj->SetAxisRange(0,TMath::Pi());
  Proj->SetMarkerColor(2);
	// double shift = Proj->GetMinimum();
	int rangedown = Proj->FindBin(0);
	int rangeup = Proj->FindBin(TMath::Pi());
	//!
	for(int i = rangedown; i < rangeup + 1; i++){
		// double value = Proj[imult][ipt]->GetBinContent(i);
		// double error = Proj->GetBinError(i)*0.001;
		
		// value = value - shift;
		// Proj[imult][ipt]->SetBinContent(i,value);		
		// Proj->SetBinError(i,error);		
	}
	TF1 * fourier = addFitToCanvas(Proj);
	float v2_low = sqrt(fourier->GetParameter(2));
	cout<<"v2(1-2 GeV) = "<<v2_low<<endl;
  string quality[] = {"","no selection","inverted"};
  
  TLatex * lmult;
  // if(tag.compare("PbPbUPC_pptracking_452p1_forest_sorted_v3")==0)
  lmult = new TLatex(0.66, 0.95, Form("CMS PbPb %2.0f-%2.0f%%", centmin*2.5, centmax*2.5));
  if(tag.compare("HIRun2013-PromptReco-v1-HLT_PAPixelTrackMultiplicity130_FullTrack12-sort")==0)
    lmult = new TLatex(0.01, 0.95, Form("CMS pPb %d #leq N < %d", nmin, nmax));
	if(tag.compare("HIRun2013-PromptReco-v1-HLT_PAPixelTracks_Multiplicity190_v1")==0)
    lmult = new TLatex(0.01, 0.95, Form("CMS pPb' %d #leq N < %d", nmin, nmax));
	
  
  lmult->SetNDC(1);
  lmult->SetTextSize(lmult->GetTextSize()*1.5);

  TLatex * ltrackqual = new TLatex(0.8, 0.95, Form("%s", quality[0].data()));
  ltrackqual->SetNDC(1);
  ltrackqual->SetTextSize(ltrackqual->GetTextSize()*1.5);

  TLatex * lpt = new TLatex(0.01, 0.95, Form("%d #leq %sleading-jet p_{T}^{trig} < %d GeV/c", (int)leadingjetptlow, subleading.data(), (int)leadingjetpthigh));
  TLatex * lpta = new TLatex(0.01, 0.88, Form("%d #leq p_{T}^{ass} < %d GeV/c", (int)aptmin, (int)aptmax));
  TLatex * laj = new TLatex(0.66, 0.88, Form("%2.2f #leq A_{j} < %2.2f ", ajmin, ajmax));
  lpt->SetNDC(1);
  lpta->SetNDC(1);
  laj->SetNDC(1);
  lpt->SetTextSize(lpt->GetTextSize()*1.5);
  lpta->SetTextSize(lpta->GetTextSize()*1.5);
  laj->SetTextSize(laj->GetTextSize()*1.5);
    
  
	
	TCanvas * c1 = new TCanvas(Form("signal_%sleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",subleading.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)aptmin,(int)aptmax,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),"",700,700);
	if(draw.compare("corr")==0)
		Corr->Draw("Surf1");
	if(draw.compare("proj")==0)
	{
		Proj->Draw();
		fourier->Draw("same");
		// TLatex * lv2 = new TLatex(0.5, 0.77, Form("v_{2} = %2.2f", v2));
		// lv2->SetNDC(1);
		// lv2->SetTextSize(lv2->GetTextSize()*1.5);
	}
  // if(draw.compare("diff")==0)
    // Diff->Draw("surf1");
  // if(draw.compare("ratio")==0)
    // Diff->Draw("surf1");
  
  cout<<"5"<<endl;
	lmult->Draw("Same");
	lpt->Draw("Same");
	lpta->Draw("Same");
	laj->Draw("Same");
	// lv2->Draw("Same");
  if(draw.compare("corr")==0)
	{
    c1->SaveAs(Form("signal_%sleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d.png",subleading.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)aptmin,(int)aptmax,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)));
  }
	else if(draw.compare("proj")==0)
	{
    c1->SaveAs(Form("proj_%s_trg%d_%d_ass%d_%d_nmin%d_nmax%d.png",tag.data(),(int)tptmin,(int)tptmax,(int)aptmin,(int)aptmax,nmin,nmax));
	}

}

