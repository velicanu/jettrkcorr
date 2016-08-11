#include <string>
#include <iostream>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>

void multiplot(string path, string tag, string calo, int jetptmin, int jetptmax, int aptmin, int aptmax, int detaordphi = 0, string tag2 = "", string calo2 = "", string tag3 = "", string calo3 = "")
{
  std::string ppb5 = "pPb_5TeV";
  std::string pp2 = "pp_276TeV";
  std::string pp7 = "pp_7TeV";
  bool is_ppb5 = (tag.find(ppb5) != std::string::npos);
  bool is_pp2 = (tag.find(pp2) != std::string::npos);
  bool is_pp7 = (tag.find(pp7) != std::string::npos);
  string sqrts = "";
  if(is_pp2) sqrts = "2.76";
  if(is_pp7) sqrts = "7";
  if(is_ppb5) sqrts = "5.02";
  
  TFile *_file = TFile::Open(Form("%s/%s_%d_%d_%s_trkqaul-1_nmin-1_nmax-1_cmin0_cmax0_aj0_10_tptmin-1_tptmax-1_aptmin%d_aptmax%d.root",path.data(),tag.data(),jetptmin,jetptmax,calo.data(),aptmin,aptmax));
  TH2D * sig = (TH2D*)_file->Get(Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin0_cmax0_ajmin0_ajmax100", calo.data(), jetptmin, jetptmax,aptmin,aptmax));
  TH1D * jetpt = (TH1D*)_file->Get(Form("hjet0_pt_trg%d_%d_ass%d_%d_nmin0_nmax0", jetptmin, jetptmax,aptmin,aptmax));
  double njets = jetpt->Integral();
  sig->Scale(1/njets);
  sig->SetAxisRange(-0.3,0.3,"X");
  sig->SetAxisRange(-0.3,0.3,"Y");
  TH1D * sigx = (TH1D*) sig->ProjectionX();
  TH1D * sigy = (TH1D*) sig->ProjectionY();
  TLatex * lpt1 = new TLatex(0.01, 0.94, Form("%d<jet p_{T}<%d GeV/c", jetptmin, jetptmax));
  TLatex * lpt2 = new TLatex(0.01, 0.86, Form("%d<trkPt<%d GeV/c",aptmin,aptmax));
  TLatex * lpt3 = new TLatex(0.65, 0.94, Form("#sqrt{s}=%s TeV", sqrts.data()));
  lpt1->SetNDC(1);
  lpt1->SetTextSize(lpt1->GetTextSize()*1.5);
  lpt2->SetNDC(1);
  lpt2->SetTextSize(lpt2->GetTextSize()*1.5);
  lpt3->SetNDC(1);
  lpt3->SetTextSize(lpt3->GetTextSize()*1.5);
  // TCanvas * c1 = new TCanvas("c1");
  // sig->Draw("surf1");
  // lpt1->Draw();
  // lpt2->Draw();
  // lpt3->Draw();
  // c1->SetLogz();
  // c1->SaveAs(Form("%s_init%d_%d_jets_%s_apt%d_%d_CF.png", tag.data(), jetptmin, jetptmax, calo.data(),aptmin,aptmax));
  TH1D * sig2x;
  TH1D * sig2y;
  TH1D * sig3x;
  TH1D * sig3y;
  if(!tag2.empty() && !calo2.empty())
  {
    cout<<Form("%s/%s_%d_%d_%s_trkqaul-1_nmin-1_nmax-1_cmin0_cmax0_aj0_10_tptmin-1_tptmax-1_aptmin%d_aptmax%d.root",path.data(),tag2.data(),jetptmin,jetptmax,calo2.data(),aptmin,aptmax)<<endl;
    TFile *_file2 = TFile::Open(Form("%s/%s_%d_%d_%s_trkqaul-1_nmin-1_nmax-1_cmin0_cmax0_aj0_10_tptmin-1_tptmax-1_aptmin%d_aptmax%d.root",path.data(),tag2.data(),jetptmin,jetptmax,calo2.data(),aptmin,aptmax));
    cout<<Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin0_cmax0_ajmin0_ajmax100", calo2.data(), jetptmin, jetptmax,aptmin,aptmax)<<endl;
    TH2D * sig2 = (TH2D*)_file2->Get(Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin0_cmax0_ajmin0_ajmax100", calo2.data(), jetptmin, jetptmax,aptmin,aptmax));
    TH1D * jetpt2 = (TH1D*)_file2->Get(Form("hjet0_pt_trg%d_%d_ass%d_%d_nmin0_nmax0", jetptmin, jetptmax,aptmin,aptmax));
    double njets2 = jetpt2->Integral();
    sig2->Scale(1/njets2);
    sig2->SetAxisRange(-0.3,0.3,"X");
    sig2->SetAxisRange(-0.3,0.3,"Y");
    sig2x = (TH1D*) sig2->ProjectionX();
    sig2y = (TH1D*) sig2->ProjectionY();
    sig2x->SetMarkerColor(kRed);
    sig2y->SetMarkerColor(kRed);
    sig2x->SetMarkerStyle(21);
    sig2y->SetMarkerStyle(21);
  }
  if(!tag3.empty() && !calo3.empty())
  {
    cout<<Form("%s/%s_%d_%d_%s_trkqaul-1_nmin-1_nmax-1_cmin0_cmax0_aj0_10_tptmin-1_tptmax-1_aptmin%d_aptmax%d.root",path.data(),tag3.data(),jetptmin,jetptmax,calo3.data(),aptmin,aptmax)<<endl;
    TFile *_file3 = TFile::Open(Form("%s/%s_%d_%d_%s_trkqaul-1_nmin-1_nmax-1_cmin0_cmax0_aj0_10_tptmin-1_tptmax-1_aptmin%d_aptmax%d.root",path.data(),tag3.data(),jetptmin,jetptmax,calo3.data(),aptmin,aptmax));
    TH2D * sig3 = (TH2D*)_file3->Get(Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin0_cmax0_ajmin0_ajmax100", calo3.data(), jetptmin, jetptmax,aptmin,aptmax));
    TH1D * jetpt3 = (TH1D*)_file3->Get(Form("hjet0_pt_trg%d_%d_ass%d_%d_nmin0_nmax0", jetptmin, jetptmax,aptmin,aptmax));
    double njets3 = jetpt3->Integral();
    sig3->Scale(1/njets3);
    sig3->SetAxisRange(-0.3,0.3,"X");
    sig3->SetAxisRange(-0.3,0.3,"Y");
    sig3x = (TH1D*) sig3->ProjectionX();
    sig3y = (TH1D*) sig3->ProjectionY();
    sig3x->SetMarkerColor(kBlue);
    sig3y->SetMarkerColor(kBlue);
    sig3x->SetMarkerStyle(24);
    sig3y->SetMarkerStyle(24);
  }

  
  TF1 * fa1 = new TF1("fa1","gaus(0)",-0.045,0.045);
  fa1->SetParameters(2.0e-02,0.0,2.0e-02);
  // fa1->FixParameter(1,0);
  // sigx->Fit(fa1,"RNO");
  if(detaordphi==0)
  {
    sigx->Fit(fa1,"RNO");
    sigx->Draw();
  }
  else
  {
    sigy->Fit(fa1,"RNO");
    sigy->Draw();
  }
  
  fa1->Draw("same");
  lpt1->Draw();
  lpt2->Draw();
  TLegend *leg = new TLegend(0.62,0.75,0.96,0.90,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.039);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(19);
  leg->SetFillStyle(0);
  
  leg->AddEntry(sigx,Form("%s",tag.data()),"lpfe");
  if(!tag2.empty() && !calo2.empty())
  {
    sig2x->Draw("same");
    leg->AddEntry(sig2x,Form("%s",tag2.data()),"lpfe");
  }
  if(!tag3.empty() && !calo3.empty())
  {
    sig3x->Draw("same");
    leg->AddEntry(sig3x,Form("%s",tag3.data()),"lpfe");
  }
  
  leg->Draw();
  lpt3->Draw();
  // c1->SetLogy();
  // c1->SaveAs(Form("%s_init%d_%d_jets_%s_apt%d_%d_deta_%s_%s.png", tag.data(), jetptmin, jetptmax, calo.data(),aptmin,aptmax, tag2.data(), tag3.data()));
  /*
  sigy->Draw();
  TLegend *leg2 = new TLegend(0.62,0.75,0.96,0.90,NULL,"brNDC");
  leg2->SetBorderSize(0);
  leg2->SetTextFont(62);
  leg2->SetTextSize(0.039);
  leg2->SetLineColor(1);
  leg2->SetLineStyle(1);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(19);
  leg2->SetFillStyle(0);
  
  leg2->AddEntry(sigy,Form("%s",tag.data()),"lpfe");
  if(!tag2.empty() && !calo2.empty())
  {
    sig2y->Draw("same");
    leg2->AddEntry(sig2y,Form("%s",tag2.data()),"lpfe");
  }
  if(!tag3.empty() && !calo3.empty())
  {
    sig3y->Draw("same");
    leg2->AddEntry(sig3y,Form("%s",tag3.data()),"lpfe");
  }
  lpt1->Draw();
  lpt2->Draw();
  leg2->Draw();
  lpt3->Draw();
  // c1->SaveAs(Form("%s_init%d_%d_jets_%s_apt%d_%d_dphi_%s_%s.png", tag.data(), jetptmin, jetptmax, calo.data(),aptmin,aptmax, tag2.data(), tag3.data()));
  */
}
