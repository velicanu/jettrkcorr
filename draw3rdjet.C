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
#include <TMath.h>


void draw3rdjet(int filenum = 0, int whichplot = 0, bool save = 0)
{
  string filename;
  string savetag;
  TLatex * l_njets;
  if(filenum == 0) {
    savetag = "only2";
    filename = "merged/only2jets04_ak3PF_trkqaul-1_nmin-1_nmax-1_cmin0_cmax0_aj0_10_tptmin-1_tptmax-1_aptmin1_aptmax2.root";
    l_njets = new TLatex(.03,.95,"nJets = 2");
  } else {
    filename = "merged/morethan2jets04_ak3PF_trkqaul-1_nmin-1_nmax-1_cmin0_cmax0_aj0_10_tptmin-1_tptmax-1_aptmin1_aptmax2.root";
    savetag = "morethan2";
    l_njets = new TLatex(.03,.95,"nJets > 2");
  }
  TLatex * l_jetcut = new TLatex(0.2,.95,"  gen jet p_{T} > 25     gen trk 1 > p_{T} > 2");;
  TLatex * l_leading = new TLatex(0.03,.88,"leading jet");;
  TLatex * l_subleading = new TLatex(0.03,.88,"subleading jet");;
  l_jetcut->SetNDC(1);
  l_njets->SetNDC(1);
  l_leading->SetNDC(1);
  l_subleading->SetNDC(1);
  float pi = TMath::Pi();
  
  TFile *_file0 = TFile::Open(filename.data());
  TH2D * signal_jet0 = (TH2D*) _file0->Get("signal_ak3PF_leadingjet120_300_ass1_2_cmin0_cmax0_ajmin0_ajmax100");
  TH2D * signal_jet1 = (TH2D*) _file0->Get("signal_ak3PF_subleadingjet120_300_ass1_2_cmin0_cmax0_ajmin0_ajmax100");
  // TH2D * shiftsignal_jet1 = (TH2D*) _file0->Get("signal_ak3PF_shifted_subleadingjet120_300_ass1_2_cmin0_cmax0_ajmin0_ajmax100");
  TH2D * signal_jet2 = (TH2D*) _file0->Get("signal_ak3PF_3rdjet120_300_ass1_2_cmin0_cmax0_ajmin0_ajmax100");
  TH2D * background_jet0 = (TH2D*) _file0->Get("background_ak3PF_leadingjet120_300_ass1_2_cmin0_cmax0_ajmin0_ajmax100");
  TH2D * background_jet1 = (TH2D*) _file0->Get("background_ak3PF_subleadingjet120_300_ass1_2_cmin0_cmax0_ajmin0_ajmax100");
  TH2D * background_jet2 = (TH2D*) _file0->Get("background_ak3PF_3rdjet120_300_ass1_2_cmin0_cmax0_ajmin0_ajmax100");
  TH2D * hjet0_phi = (TH2D*) _file0->Get("hjet0_phi_trg120_300_ass1_2_nmin0_nmax0");
  TH2D * hjet1_phi = (TH2D*) _file0->Get("hjet1_phi_trg120_300_ass1_2_nmin0_nmax0");
  TH2D * hjet2_phi = (TH2D*) _file0->Get("hjet2_phi_trg120_300_ass1_2_nmin0_nmax0");
  TH1D * hj3pt    = (TH1D*) _file0->Get("hj2_hj3pt_trg120_300_ass1_2_nmin0_nmax0");
  TH1D * hj13dphi = (TH1D*) _file0->Get("hj2_hj13dphi_trg120_300_ass1_2_nmin0_nmax0");
  TH1D * hj23dphi = (TH1D*) _file0->Get("hj2_hj23dphi_trg120_300_ass1_2_nmin0_nmax0");
  TH1D * hj13deta = (TH1D*) _file0->Get("hj2_hj13deta_trg120_300_ass1_2_nmin0_nmax0");
  TH1D * hj23deta = (TH1D*) _file0->Get("hj2_hj23deta_trg120_300_ass1_2_nmin0_nmax0");
  TH1D * hj123dpt = (TH1D*) _file0->Get("hj2_hj123dpt_trg120_300_ass1_2_nmin0_nmax0");
  
  int nbinsx = signal_jet0->GetNbinsX();
  int nbinsy = signal_jet0->GetNbinsY();

  // signal_jet0->Scale(-1);
  // signal_jet0->Add(signal_jet1);
  // signal_jet0->Draw("surf1");
  // return;
  
  TCanvas * c1 = new TCanvas("c1");
  
  // signal_jet0->Draw("surf1");
  // c1->SaveAs("mcsignal_leadingjet.png");
  // signal_jet1->Draw("surf1");
  // c1->SaveAs("mcsignal_subleadingjet.png");
  // signal_jet2->Draw("surf1");
  // c1->SaveAs("mcsignal_3rdjet.png");
  // background_jet0->Draw("surf1");
  // c1->SaveAs("mcbackground_leadingjet.png");
  // background_jet1->Draw("surf1");
  // c1->SaveAs("mcbackground_subleadingjet.png");
  // background_jet2->Draw("surf1");
  // c1->SaveAs("mcbackground_3rdjet.png");
  
  // TLatex * label_jet0 = new TLatex(.03,.94,"leading");
  // label_jet0->SetNDC(1);
  double d_background_jet0 = background_jet0->GetBinContent(background_jet0->FindBin(0,0));
  double d_binarea_jet0 = signal_jet0->GetXaxis()->GetBinWidth(1)*signal_jet0->GetYaxis()->GetBinWidth(1);
  double n_jet0 = hjet0_phi->Integral();
  signal_jet0->Divide(background_jet0);
  signal_jet0->Scale(d_background_jet0);
  signal_jet0->Scale(1/d_binarea_jet0);
  signal_jet0->Scale(1/n_jet0);
  
  // TLatex * label_jet1 = new TLatex(.03,.94,"leading");
  // label_jet1->SetNDC(1);
  double d_background_jet1 = background_jet1->GetBinContent(background_jet1->FindBin(0,0));
  double d_binarea_jet1 = signal_jet1->GetXaxis()->GetBinWidth(1)*signal_jet1->GetYaxis()->GetBinWidth(1);
  double n_jet1 = hjet1_phi->Integral();
  signal_jet1->Divide(background_jet1);
  signal_jet1->Scale(d_background_jet1);
  signal_jet1->Scale(1/d_binarea_jet1);
  signal_jet1->Scale(1/n_jet1);
  
  
  
  signal_jet0->SetAxisRange(-3,3,"X");
  signal_jet0->SetAxisRange(0,2,"Z");
  signal_jet1->SetAxisRange(-3,3,"X");
  signal_jet1->SetAxisRange(0,2,"Z");

  if( whichplot == 0 ) {
    signal_jet0->Draw("surf1");
    // cout<<"here"<<endl;
    l_njets->Draw();
    l_jetcut->Draw();
    l_leading->Draw();
    if(save) { c1->SaveAs(Form("correlation_leading_%s.png",savetag.data())); }
    return;
  }
  if( whichplot == -1 ) {
    signal_jet1->Draw("surf1");
    l_njets->Draw();
    l_jetcut->Draw();
    l_subleading->Draw();
    if(save) { c1->SaveAs(Form("correlation_subleading_%s.png",savetag.data())); }
    return;
  }
  // return;
  double nbins_jet0 = signal_jet0->GetXaxis()->FindBin(3)-signal_jet0->GetXaxis()->FindBin(1.5);
  TH1D * dphi_jet0 = signal_jet0->ProjectionY("dphi_jet0",signal_jet0->GetXaxis()->FindBin(1.5),signal_jet0->GetXaxis()->FindBin(3));
  TH1D * dphi_jet1 = signal_jet1->ProjectionY("dphi_jet1",signal_jet1->GetXaxis()->FindBin(1.5),signal_jet1->GetXaxis()->FindBin(3));
  dphi_jet0->Scale(1/nbins_jet0);
  dphi_jet1->Scale(1/nbins_jet0);
  for(int k = 1; k<nbinsx/2+1; k++) {
    float diff = signal_jet0->GetBinContent(k) - dphi_jet1->GetBinContent(k);
    cout << diff << endl;
    // cout<<signal_jet0<<" "<<signal_jet1<< endl;

    dphi_jet0->SetBinContent(k+nbinsx/2,dphi_jet1->GetBinContent(k));
    dphi_jet0->SetBinError(k+nbinsx/2,dphi_jet1->GetBinError(k));
  }
  
  
  
  string formula_jet0 = "[0]*(1+2*[1]*cos(x)+2*[2]*cos(2*x)+[3]*exp(-0.5*((x-pi)/[4])**2))";
  TF1 * sideband_jet0 = new TF1("sideband_jet0",formula_jet0.data(),-pi/2.0,3*pi/2.0);
  sideband_jet0->SetParameters(0.3,0.1,0.1,0.001,0.1);
  dphi_jet0->Fit("sideband_jet0","RNO");
  if( whichplot == 1 ) {
    dphi_jet0->Draw();
    sideband_jet0->Draw("same");
    // label_jet0->Draw("same");
    l_njets->Draw();
    l_jetcut->Draw();
    // l_subleading->Draw();
    if(save) { c1->SaveAs(Form("sidebandfit_leading_%s.png",savetag.data())); }
    return;
  }

  // return;
  
  TH2D * subtract_jet0 = new TH2D("subtracted_jet0","",40,-5,5,40,-0.5*TMath::Pi(),1.5*TMath::Pi());
  for(int i = 0 ; i < 40 ; ++i)
  {
    double val = sideband_jet0->Eval(dphi_jet0->GetBinCenter(i));
    for(int j = 0 ; j < 40 ; ++j)
    {
      subtract_jet0->SetBinContent(j,i,val);
    }
  }
  // dphi_jet0->Draw();
  // sideband_jet0->Draw("same");
  signal_jet0->Add(subtract_jet0,-1);
  signal_jet0->SetAxisRange(-1.5,1.5,"X");
  signal_jet0->SetAxisRange(-1.5,1.5,"Y");
  if( whichplot == 2 ) {
    signal_jet0->Draw("surf1");
    l_njets->Draw();
    l_jetcut->Draw();
    if(save) { c1->SaveAs(Form("correlation_sidebandsub_leading_%s.png",savetag.data())); }
    return;
  }
  
  double nbinsub_jet0 = signal_jet0->GetXaxis()->FindBin(1.5)-signal_jet0->GetXaxis()->FindBin(-1.5);
  TH1D * dphisub_jet0 = signal_jet0->ProjectionY("dphisub_jet0",signal_jet0->GetXaxis()->FindBin(-1.5),signal_jet0->GetXaxis()->FindBin(1.5));
  dphisub_jet0->Scale(1/nbinsub_jet0);
  if( whichplot == 3 ) {
    dphisub_jet0->Draw();
    l_njets->Draw();
    l_jetcut->Draw();
    if(save) { c1->SaveAs(Form("dphi_postsub_leading_%s.png",savetag.data())); }
    return;
  }
  double nbinsub_jet0 = signal_jet0->GetYaxis()->FindBin(1.5)-signal_jet0->GetYaxis()->FindBin(-1.5);
  TH1D * detasub_jet0 = signal_jet0->ProjectionX("detasub_jet0",signal_jet0->GetYaxis()->FindBin(-1.5),signal_jet0->GetYaxis()->FindBin(1.5));
  detasub_jet0->Scale(1/nbinsub_jet0);
  if( whichplot == 4 ) {
    detasub_jet0->Draw();
    l_njets->Draw();
    l_jetcut->Draw();
    if(save) { c1->SaveAs(Form("deta_postsub_leading_%s.png",savetag.data())); }
    return;
  }
}
