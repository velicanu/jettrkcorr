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

// void sigmavpt(string path, string tag, string calo, int jetptmin, int jetptmax, int ietaphi = 0, int nptbins = 3)
void sigmavpt(int jetptmin, int jetptmax, int ietaphi = 0)
{
  int aaptmin [9] = {19,24,29,35,41,48,61,73,86};
  int aaptmax [9] = {24,29,35,41,48,61,73,86,103};
  string datatags [3] = {"genpp_276TeV","genpPb_5TeV","genpp_7TeV"};
  string mctags [3] = {"mcpp_276TeV","mcpPb_5TeV","mcpp_7TeV"};
  // string datacalos[3] = {"ak3PF","akPu3PF","ak3PF"};
  string datacalos[3] = {"ak3PF","ak3PF","ak3PF"};
  string mccalos[3] = {"ak3PF","ak3PF","ak3PF"};
  float roots [3] = { 2.76 , 5.02 , 7 };
  int colors [3] = {kRed,kBlue,kBlack};
  int datastyles [3] = {20,21,33};
  int mcstyles [3] = {24,25,27};
  
  TFile *_file;
  TH2D * sig;
  TH1D * jetpt;

  TH1D * sigx;
  TH1D * sigy;
  TLatex * lpt1 = new TLatex(0.01, 0.94, Form("%d<jet p_{T}<%d GeV/c", jetptmin, jetptmax));
  // TLatex * lpt3 = new TLatex(0.65, 0.94, Form("#sqrt{s}=%s TeV", sqrts.data()));
  TCanvas * c1 = new TCanvas("c1");
  lpt1->SetNDC(1);
  lpt1->SetTextSize(lpt1->GetTextSize()*1.5);
  // lpt3->SetNDC(1);
  // lpt3->SetTextSize(lpt3->GetTextSize()*1.5);
  TH1D * sigmavspt [2][3]; // [data or mc][roots]
  TLegend *leg = new TLegend(0.25,0.20,0.59,0.45,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.039);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(19);
  leg->SetFillStyle(0);
  
  string detadphi [2] = {"#eta","#phi"};
  
  for(int j = 0 ; j < 3 ; ++j) //roots loop
  {
    sigmavspt[0][j] = new TH1D(Form("sigmavspt_%d",j),Form(";p_{T} (GeV);#sigma#Delta%s",detadphi[ietaphi].data()),109,0.5,109.5);
    for(int i = 0 ; i < 9 ; ++i) //trk pt loop
    {
      // cout<<i<<" "<<j<<" "<<Form("merged/%s_%d_%d_%s_trkqaul-1_nmin-1_nmax-1_cmin0_cmax0_aj0_10_tptmin-1_tptmax-1_aptmin%d_aptmax%d.root",datatags[j].data(),jetptmin,jetptmax,datacalos[j].data(),aaptmin[i],aaptmax[i])<<endl;
      _file = TFile::Open(Form("merged/%s_%d_%d_%s_trkqaul-1_nmin-1_nmax-1_cmin0_cmax0_aj0_10_tptmin-1_tptmax-1_aptmin%d_aptmax%d.root",datatags[j].data(),jetptmin,jetptmax,datacalos[j].data(),aaptmin[i],aaptmax[i]));
      cout<<Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin0_cmax0_ajmin0_ajmax100", datacalos[j].data(), jetptmin, jetptmax,aaptmin[i],aaptmax[i])<<endl;
      sig = (TH2D*)_file->Get(Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin0_cmax0_ajmin0_ajmax100", datacalos[j].data(), jetptmin, jetptmax,aaptmin[i],aaptmax[i]));
      jetpt = (TH1D*)_file->Get(Form("hjet0_pt_trg%d_%d_ass%d_%d_nmin0_nmax0", jetptmin, jetptmax,aaptmin[i],aaptmax[i]));
      double njets = jetpt->Integral();
      sig->Scale(1/njets);
      sig->SetAxisRange(-0.3,0.3,"X");
      sig->SetAxisRange(-0.3,0.3,"Y");
      sigx = (TH1D*) sig->ProjectionX();
      sigy = (TH1D*) sig->ProjectionY();
      TF1 * fa1 = new TF1("fa1","gaus(0)",-0.045,0.045);
      fa1->SetParameters(2.0e-02,0.0,2.0e-02);
      if( ietaphi == 0) sigx->Fit(fa1,"RNO");
      else             sigy->Fit(fa1,"RNO");
      sigmavspt[0][j]->SetBinContent(sigmavspt[0][j]->FindBin(aaptmin[i]),fabs(fa1->GetParameter(2)));
      sigmavspt[0][j]->SetBinError(sigmavspt[0][j]->FindBin(aaptmin[i]),fa1->GetParError(2));
      // cout<<"end "<<i<<" "<<j<<endl;
    }
    sigmavspt[0][j]->SetMarkerColor(colors[j]);
    sigmavspt[0][j]->SetMarkerStyle(datastyles[j]);
    sigmavspt[0][j]->SetAxisRange(0,0.035,"Y");
    sigmavspt[0][j]->GetXaxis()->CenterTitle();
    sigmavspt[0][j]->GetYaxis()->CenterTitle();
    sigmavspt[0][j]->GetYaxis()->SetTitleOffset(1.5);
    if(j==0) sigmavspt[0][j]->Draw("pe");
    else sigmavspt[0][j]->Draw("pesame");
    // leg->AddEntry(sigmavspt[j],Form("%d < trkPt < %d",aaptmin[j],aaptmax[j]),"lpfe");
    leg->AddEntry(sigmavspt[0][j],Form("%s",datatags[j].data()),"lpfe");
  }
  for(int j = 0 ; j < 3 ; ++j) //roots loop
  {
    sigmavspt[1][j] = new TH1D(Form("sigmavspt_%d",j),Form(";p_{T} (GeV);#sigma#Delta%s",detadphi[ietaphi].data()),109,0.5,109.5);
    for(int i = 0 ; i < 9 ; ++i) //trk pt loop
    {
      // cout<<i<<" "<<j<<" "<<Form("merged/%s_%d_%d_%s_trkqaul-1_nmin-1_nmax-1_cmin0_cmax0_aj0_10_tptmin-1_tptmax-1_aptmin%d_aptmax%d.root",mctags[j].data(),jetptmin,jetptmax,mccalos[j].data(),aaptmin[i],aaptmax[i])<<endl;
      _file = TFile::Open(Form("merged/%s_%d_%d_%s_trkqaul-1_nmin-1_nmax-1_cmin0_cmax0_aj0_10_tptmin-1_tptmax-1_aptmin%d_aptmax%d.root",mctags[j].data(),jetptmin,jetptmax,mccalos[j].data(),aaptmin[i],aaptmax[i]));
      cout<<Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin0_cmax0_ajmin0_ajmax100", mccalos[j].data(), jetptmin, jetptmax,aaptmin[i],aaptmax[i])<<endl;
      sig = (TH2D*)_file->Get(Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin0_cmax0_ajmin0_ajmax100", mccalos[j].data(), jetptmin, jetptmax,aaptmin[i],aaptmax[i]));
      jetpt = (TH1D*)_file->Get(Form("hjet0_pt_trg%d_%d_ass%d_%d_nmin0_nmax0", jetptmin, jetptmax,aaptmin[i],aaptmax[i]));
      double njets = jetpt->Integral();
      sig->Scale(1/njets);
      sig->SetAxisRange(-0.3,0.3,"X");
      sig->SetAxisRange(-0.3,0.3,"Y");
      sigx = (TH1D*) sig->ProjectionX();
      sigy = (TH1D*) sig->ProjectionY();
      TF1 * fa1 = new TF1("fa1","gaus(0)",-0.045,0.045);
      fa1->SetParameters(2.0e-02,0.0,2.0e-02);
      if( ietaphi == 0) sigx->Fit(fa1,"RNO");
      else             sigy->Fit(fa1,"RNO");
      sigmavspt[1][j]->SetBinContent(sigmavspt[1][j]->FindBin(aaptmin[i]),fabs(fa1->GetParameter(2)));
      sigmavspt[1][j]->SetBinError(sigmavspt[1][j]->FindBin(aaptmin[i]),fa1->GetParError(2));
      // cout<<"end "<<i<<" "<<j<<endl;
    }
    sigmavspt[1][j]->SetMarkerColor(colors[j]);
    sigmavspt[1][j]->SetMarkerStyle(mcstyles[j]);
    sigmavspt[1][j]->SetAxisRange(0,0.035,"Y");
    sigmavspt[1][j]->GetXaxis()->CenterTitle();
    sigmavspt[1][j]->GetYaxis()->CenterTitle();
    sigmavspt[1][j]->GetYaxis()->SetTitleOffset(1.5);
    // if(j==0) sigmavspt[1][j]->Draw("pe");
    sigmavspt[1][j]->Draw("pesame");
    // leg->AddEntry(sigmavspt[j],Form("%d < trkPt < %d",aaptmin[j],aaptmax[j]),"lpfe");
    leg->AddEntry(sigmavspt[1][j],Form("%s",mctags[j].data()),"lpfe");
  }    
    
    // sigmavspt[i] = new TH1D(Form("sigmavspt_%d",i),Form(";#sqrt{s} (TeV);#sigma#Delta%s",detadphi[ietaphi].data()),50,1,8);
      // lpt1->Draw();
      // lpt3->Draw();
      
      // sigx->Draw();  
      // fa1->Draw("same");
      // lpt1->Draw();
      // cout<<fa1->GetParameter(2)<<endl;
    //}
  //}
  lpt1->Draw();
  leg->Draw();
  string ddetadphi [2] = {"eta","phi"};
  
  c1->SaveAs(Form("d%s_jet%d_%d_widths.png",ddetadphi[ietaphi].data(),jetptmin,jetptmax));
  c1->SaveAs(Form("d%s_jet%d_%d_widths.pdf",ddetadphi[ietaphi].data(),jetptmin,jetptmax));
  // c1->SaveAs(Form("mcd%s_jet%d_%d_widths.png",ddetadphi[ietaphi].data(),jetptmin,jetptmax));
  // c1->SaveAs(Form("mcd%s_jet%d_%d_widths.pdf",ddetadphi[ietaphi].data(),jetptmin,jetptmax));
  
  // sig->Draw("surf1");
  

}
