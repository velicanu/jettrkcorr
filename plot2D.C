#include <string>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TString.h>

void plot2D(string path , string tag , int jetptmin, int jetptmax)
{
  TFile *_file60 = TFile::Open(Form("%s/%s_%d_%d_akPu3Calo_trkqaul-1_nmin-1_nmax-1_cmin0_cmax0_aj0_10_tptmin-1_tptmax-1_aptmin1_aptmax2.root",path.data(),tag.data(),jetptmin,jetptmax));
  TH2D * sig60 = (TH2D*)_file60->Get("signal_akPu3Calo_leadingjet60_80_ass1_2_cmin0_cmax0_ajmin0_ajmax100");
  TH2D * bak60 = (TH2D*)_file60->Get("background_akPu3Calo_leadingjet60_80_ass1_2_cmin0_cmax0_ajmin0_ajmax100");
  sig60->Divide(bak60);
  TH1D * sigx60 = (TH1D*) sig60->ProjectionX();
  TH1D * sigy60 = (TH1D*) sig60->ProjectionY();
  TLatex * l60pt1 = new TLatex(0.01, 0.94, "60<jet p_{T}<80 GeV/c");
  TLatex * l60pt2 = new TLatex(0.01, 0.86, "1<trkPt<2 GeV/c");
  TLatex * l60pt3 = new TLatex(0.65, 0.94, "#sqrt{s}=5.02 TeV");
  l60pt1->SetNDC(1);
  l60pt1->SetTextSize(l60pt1->GetTextSize()*1.5);
  l60pt2->SetNDC(1);
  l60pt2->SetTextSize(l60pt2->GetTextSize()*1.5);
  l60pt3->SetNDC(1);
  l60pt3->SetTextSize(l60pt3->GetTextSize()*1.5);
  
  TCanvas * c1 = new TCanvas("c1");
  
  sig60->Draw("surf1");
  l60pt1->Draw();
  l60pt2->Draw();
  l60pt3->Draw();
  c1->SaveAs("pPb502TeVinit60_80jets_akPu3Calo_CF.png");
  sigx60->Draw();
  l60pt1->Draw();
  l60pt2->Draw();
  l60pt3->Draw();
  c1->SaveAs("pPb502TeVinit60_80jets_akPu3Calo_deta.png");
  sigy60->Draw();
  l60pt1->Draw();
  l60pt2->Draw();
  l60pt3->Draw();
  c1->SaveAs("pPb502TeVinit60_80jets_akPu3Calo_dphi.png");
}