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

void sigmaplot(string path, string tag, string calo, int jetptmin, int jetptmax, int ietaphi = 0, int nptbins = 3)
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
  int aaptmin [3] = {60,73,86};
  int aaptmax [3] = {73,86,103};
  string tags [3] = {"pp_276TeV","pPb_5TeV","pp_7TeV"};
  // string tags [3] = {"mcpp_276TeV","mcpPb_5TeV","mcpp_7TeV"};
  string calos[3] = {"ak3PF","akPu3PF","ak3PF"};
  // string calos[3] = {"akPu3PF","akPu3PF","ak3PF"};
  float roots [3] = { 2.76 , 5.02 , 7 };
  int colors [3] = {kRed,kBlue,kBlack};
  int styles [3] = {20,21,24};
  
  TFile *_file;
  TH2D * sig;
  TH1D * jetpt;

  TH1D * sigx;
  TH1D * sigy;
  TLatex * lpt1 = new TLatex(0.01, 0.94, Form("%d<jet p_{T}<%d GeV/c", jetptmin, jetptmax));
  TLatex * lpt3 = new TLatex(0.65, 0.94, Form("#sqrt{s}=%s TeV", sqrts.data()));
  TCanvas * c1 = new TCanvas("c1");
  lpt1->SetNDC(1);
  lpt1->SetTextSize(lpt1->GetTextSize()*1.5);
  lpt3->SetNDC(1);
  lpt3->SetTextSize(lpt3->GetTextSize()*1.5);
  TH1D * sigmavsroots [3];
  TLegend *leg = new TLegend(0.54,0.24,0.88,0.39,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.039);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(19);
  leg->SetFillStyle(0);
  
  string detadphi [2] = {"#eta","#phi"};
  
  for(int i = 0 ; i < 3 ; ++i) //trk pt loop
  {
    sigmavsroots[i] = new TH1D(Form("sigmavsroots_%d",i),Form(";#sqrt{s} (TeV);#sigma#Delta%s",detadphi[ietaphi].data()),50,1,8);
    for(int j = 0 ; j < 3 ; ++j) //roots loop
    {          
      _file = TFile::Open(Form("%s/%s_%d_%d_%s_trkqaul-1_nmin-1_nmax-1_cmin0_cmax0_aj0_10_tptmin-1_tptmax-1_aptmin%d_aptmax%d.root",path.data(),tags[j].data(),jetptmin,jetptmax,calos[j].data(),aaptmin[i],aaptmax[i]));
      sig = (TH2D*)_file->Get(Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin0_cmax0_ajmin0_ajmax100", calos[j].data(), jetptmin, jetptmax,aaptmin[i],aaptmax[i]));
      jetpt = (TH1D*)_file->Get(Form("hjet0_pt_trg%d_%d_ass%d_%d_nmin0_nmax0", jetptmin, jetptmax,aaptmin[i],aaptmax[i]));
      double njets = jetpt->Integral();
      sig->Scale(1/njets);
      sig->SetAxisRange(-0.3,0.3,"X");
      sig->SetAxisRange(-0.3,0.3,"Y");
      sigx = (TH1D*) sig->ProjectionX();
      sigy = (TH1D*) sig->ProjectionY();
      // lpt1->Draw();
      // lpt3->Draw();
      
      TF1 * fa1 = new TF1("fa1","gaus(0)",-0.045,0.045);
      fa1->SetParameters(2.0e-02,0.0,2.0e-02);
      // fa1->FixParameter(1,0);
      if( ietaphi == 0) sigx->Fit(fa1,"RNO");
      else             sigy->Fit(fa1,"RNO");
      // sigx->Draw();  
      // fa1->Draw("same");
      // lpt1->Draw();
      // cout<<fa1->GetParameter(2)<<endl;
      sigmavsroots[i]->SetBinContent(sigmavsroots[i]->FindBin(roots[j])+i,fabs(fa1->GetParameter(2)));
      sigmavsroots[i]->SetBinError(sigmavsroots[i]->FindBin(roots[j])+i,fa1->GetParError(2));
    }
    sigmavsroots[i]->SetMarkerColor(colors[i]);
    sigmavsroots[i]->SetMarkerStyle(styles[i]);
    sigmavsroots[i]->SetAxisRange(0,0.035,"Y");
    sigmavsroots[i]->GetXaxis()->CenterTitle();
    sigmavsroots[i]->GetYaxis()->CenterTitle();
    sigmavsroots[i]->GetYaxis()->SetTitleOffset(1.5);
    if(i==0) sigmavsroots[i]->Draw("pe");
    else sigmavsroots[i]->Draw("pesame");
    
    // leg->AddEntry(sigmavsroots[i],Form("%d < trkPt < %d MC",aaptmin[i],aaptmax[i]),"lpfe");
    leg->AddEntry(sigmavsroots[i],Form("%d < trkPt < %d",aaptmin[i],aaptmax[i]),"lpfe");
  }
  lpt1->Draw();
  leg->Draw();
  string ddetadphi [2] = {"eta","phi"};
  
  c1->SaveAs(Form("d%s_%s_widths.png",ddetadphi[ietaphi].data(),tag.data()));
  c1->SaveAs(Form("d%s_%s_widths.pdf",ddetadphi[ietaphi].data(),tag.data()));
  
  // sig->Draw("surf1");
  

}
