// #include "/net/hisrv0001/home/dav2105/run/CMSSW_4_4_4/src/CmsHi/JetAnalysis/macros/forest/hiForest_charge.h"
#include "HiForestAnalysis/hiForest.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>

using namespace std;

const double pi = 3.1415926;
const double pi23 = 2.0 * pi / 3.0;
double maxetatrg = 2.4;
double maxetaass = 2.4;
int ntottrig = 0;
int cent_index_start[41];
int mytrackquality = 0;

TH1D * hmult;
TH1D * hpttrg;
TH1D * hetatrg;
TH1D * hphitrg;
TH1D * hmulttrg;
TH1D * hptass;
TH1D * hetaass;
TH1D * hphiass;

TH2D * hetaphi;
TH1D * hdz;
TH1D * hdxy;
TH1D * hdz_dzerr;
TH1D * hdxy_dxyerr;
TH1D * hpterr_pt;
TH1D * hnhits;
TH1D * hchi2ndof;

TH1D * hjeteta;
TH1D * hjetphi;
TH1D * hjetpt;

TH1D * hjetdphi;

TH1D * hjtrketa;
TH1D * hjtrkphi;
TH1D * hjtrkdeta;
TH1D * hjtrkdphi;

HiForest *c;
TH2D * JetTrackSignal(int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int doptweight , int mccommand , double jetamin , double jetamax , double vzrange = 15.0, double dijetdphicut = pi23);
TH2D * JetTrackBackground(int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int doptweight , int mccommand , double jetamin , double jetamax , double vzrange = 15.0, double dijetdphicut = pi23);

int GetNTotTrig();

void jettrkana(const char * infname = "/mnt/hadoop/cms/store/user/velicanu/mergedv1_sortedforest/mergesortv1_0.root", int trackquality = 0)
{
  mytrackquality = trackquality;
  cout<<"initializing hiforest and building centrality index"<<endl;
  cout<<"running on: "<<infname<<endl;
  // c = new HiForest(infname,"forest",cPbPb);
  c = new HiForest(infname,"forest",cPPb);
  c->LoadNoTrees();
  c->hasTrackTree = true;
  c->hasEvtTree = true;
}

bool skipevent(double vzrange, int runboundary)
{
  bool doskip = false;
  if(fabs(c->evt.vz)>vzrange) doskip = true;
  if(c->evt.run>runboundary) doskip = true;
  return doskip;
}

TH2D * JetTrackSignal(int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int doptweight , int mccommand , double jetamin , double jetamax , double vzrange, double dijetdphicut)
{
  Long64_t nentries = c->GetEntries();
  if(mccommand==2)  doptweight=0;
  TH2D::SetDefaultSumw2(true);
  c->LoadNoTrees();
  c->hasEvtTree = true;
  c->hasHltTree = true;
  c->hasTrackTree = true;
  c->hasAkPu3JetTree = true;
  c->hasSkimTree = true;
  
  hjeteta    = new TH1D(Form("hjet%d_eta_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",80,-4,4);
  hjetphi    = new TH1D(Form("hjet%d_phi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#phi",80,-TMath::Pi(),TMath::Pi());
  hjetpt    = new TH1D(Form("hjet%d_pt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#phi",80,0,leadingjetpthigh);
  
  
  hjetdphi    = new TH1D(Form("hjet%d_dphi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#Delta#phi",40,0,TMath::Pi());

  hjtrketa    = new TH1D(Form("hj%d_trketa_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",80,-4,4);
  hjtrkphi    = new TH1D(Form("hj%d_trkphi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#phi",80,-TMath::Pi(),TMath::Pi());

  hjtrkdeta    = new TH1D(Form("hj%d_trkdeta_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",100,-5,5);
  hjtrkdphi    = new TH1D(Form("hj%d_trkdphi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#phi",80,0,TMath::Pi());


  // Long64_t nentries = c->GetEntries();
  // Long64_t nbytes = 0, nb = 0;
  TH2D * hJetTrackSignal;
  cout<<"before first exit"<<endl;
  if(jetindex==0) //leading jet
  {
    hJetTrackSignal = new TH2D(Form("signal_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,75,-pi,2*pi);
  }
  else if(jetindex==1) //subleading jet
  {
    hJetTrackSignal = new TH2D(Form("signal_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,75,-pi,2*pi);
  }
  else //leading jet
  {
    cout<<"Error: only jetindex 0 and 1 currently supported (leading + subleading jet)"<<endl;
    exit(1);
  }
  
  double ntottrig = 0;
  // int n_entries_in_cent_range = cent_index_start[centmax] - cent_index_start[centmin];


  // for (Long64_t jentry=0; jentry<10000;jentry++) {
  // for (Long64_t jentry=cent_index_start[centmin]; jentry<cent_index_start[centmax];jentry++) {
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    // if(jentry%1000==0) cout<<jentry-cent_index_start[centmin]<<"/"<<n_entries_in_cent_range<<endl;
    if(jentry%1000==0) cout<<jentry<<"/"<<nentries<<endl;
    c->GetEntry(jentry);
    // cout<<"gotentry "<<n_entries_in_cent_range<<endl;
    if(c->evt.hiBin < centmin || c->evt.hiBin > centmax) continue;
    //! jet cuts, vz range, more than 2 jettss, dijetdphi 7/8 pi, both jet |eta| < 2, pt range cuts, aj cuts
    // cout<<mccommand<<endl;
    // cout<<"selectEvent"<<endl;
    if(mccommand==0)
    {
      bool select = (c->skim.pHBHENoiseFilter || c->mc) && c->skim.pcollisionEventSelection;
      
      if(!select) continue; 
      // if(!c->hlt.HLT_HIJet80_v1) continue; 
    }
    // cout<<"event selected"<<endl;
    if( fabs(c->evt.vz) > vzrange ) continue;
    // if( c->akPu3PF.nref < 2 )         continue;
    double dijetdphi = fabs(c->akPu3PF.jtphi[0] - c->akPu3PF.jtphi[1]);
    if( dijetdphi > pi ) dijetdphi = 2*pi - dijetdphi;
    // cout<< dijetdphi << " vs " << dijetdphicut << endl;
    if( dijetdphi < dijetdphicut ) continue;
    // cout<<"dphicut "<<c->akPu3PF.nref<<endl;
    // int jetindex = -1;
    // while(true)
    // {
    // cout<<"here"<<endl;
    // jetindex++;
    hjeteta->Fill(c->akPu3PF.jteta[jetindex]);
    hjetphi->Fill(c->akPu3PF.jtphi[jetindex]);
    hjetdphi->Fill(dijetdphi);
    hjetpt->Fill(c->akPu3PF.jtpt[jetindex]);
    
    if( c->akPu3PF.nref < 2 ) continue;
    if( c->akPu3PF.jtpt[0] > leadingjetpthigh ) continue;
    if( c->akPu3PF.jtpt[0] < leadingjetptlow ) continue;
    if( c->akPu3PF.jtpt[1] > subleadingjetpthigh ) continue;
    if( c->akPu3PF.jtpt[1] < subleadingjetptlow ) continue;
    
    if( fabs(c->akPu3PF.jteta[0]) > jetamax || fabs(c->akPu3PF.jteta[0]) < jetamin ) continue;
    if( fabs(c->akPu3PF.jteta[1]) > jetamax || fabs(c->akPu3PF.jteta[1]) < jetamin ) continue;
    if ((c->akPu3PF.trackMax[0]/c->akPu3PF.jtpt[0])<0.01) continue;
    if ((c->akPu3PF.trackMax[1]/c->akPu3PF.jtpt[1])<0.01) continue;
    // if(jentry%1000==0) { cout<<"here"<<endl; break;}
    
    double jeteta = c->akPu3PF.jteta[jetindex];
    double jetphi = c->akPu3PF.jtphi[jetindex];
    // cout<<jeteta<<endl;
    int ntrig = 1 , nass = 0;
    
    // vector<double> asstrkEta(c->track.nTrk);
    // vector<double> asstrkPhi(c->track.nTrk);
    // vector<double> asstrkPt(c->track.nTrk);
    // vector<int> asstrkIndex(c->track.nTrk);
    if(mccommand<2)
    {
      for(int j = 0 ; j < c->track.nTrk ; ++j)
      {
        
        // if( c->track.trkPt[j]<ptasslow || c->track.trkPt[j]>ptasshigh || fabs(c->track.trkEta[j])>maxetaass) continue;
        // if( c->track.trkAlgo[j]<4 || c->track.highPurity[j])
        // {
        if(mytrackquality==0&&fabs(c->track.trkEta[j])<2.4&&c->track.highPurity[j]&&fabs(c->track.trkDz1[j]/c->track.trkDzError1[j])<3&&fabs(c->track.trkDxy1[j]/c->track.trkDxyError1[j])<3&&c->track.trkPtError[j]/c->track.trkPt[j]<0.1)
        {
          ++nass;
          // cout<<c->track.trkEta[j]<<endl;
          // asstrkEta.push_back(c->track.trkEta[j]);
          // asstrkPhi.push_back(c->track.trkPhi[j]);
          // asstrkPt.push_back(c->track.trkPt[j]);
          // asstrkIndex.push_back(j);
          // cout<<j<<" : "<<c->track.trkEta[j]<<" : "<<nass<<endl;
        }
      }
    }

    for(int i = 0 ; i < ntrig ; ++i)
    {
      ntottrig += 1;

      // cout<<nass<<endl;
      for(int j = 0 ; j < c->track.nTrk ; ++j)
      {
        // cout<<mytrackquality<<endl; //mytrackquality==0&&
        if(fabs(c->track.trkEta[j])<2.4&&c->track.highPurity[j]&&fabs(c->track.trkDz1[j]/c->track.trkDzError1[j])<3&&fabs(c->track.trkDxy1[j]/c->track.trkDxyError1[j])<3&&c->track.trkPtError[j]/c->track.trkPt[j]<0.1)
        {
          // cout<<asstrkIndex[k]<<" : "<<asstrkEta[asstrkIndex[k]]<<endl;
          double deta = fabs(jeteta-c->track.trkEta[j]);
          double dphi = fabs(jetphi-c->track.trkPhi[j]);
          // double ptw = asstrkPt[k];
          if( dphi > pi ) dphi = 2*pi - dphi;
          
          hjtrketa->Fill(c->track.trkEta[j]);
          hjtrkphi->Fill(c->track.trkPhi[j]);
          hjtrkdeta->Fill(deta);
          hjtrkdphi->Fill(dphi);
          

          double effweight = c->track.trkPt[j];
          // effweight = 1.0;
          // if(doptweight!=0) effweight = c->getTrackCorrection(j);
          hJetTrackSignal->Fill(deta,dphi,effweight);
          hJetTrackSignal->Fill(-deta,dphi,effweight);
          hJetTrackSignal->Fill(deta,-dphi,effweight);
          hJetTrackSignal->Fill(-deta,-dphi,effweight);
          hJetTrackSignal->Fill(deta,(2*pi)-dphi,effweight);
          hJetTrackSignal->Fill(-deta,(2*pi)-dphi,effweight);
          
        }
      }
    }
    // cout<<"end"<<endl;
    // }

    // if(jentry>100000) break;
  }
  // if( ntottrig != 0 ) hJetTrackSignal->Scale(1/ntottrig);
  c->ResetBooleans();
  return hJetTrackSignal;
}


TH2D * JetTrackBackground(int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int doptweight , int mccommand , double jetamin , double jetamax , double vzrange, double dijetdphicut)
{
  Long64_t nentries = c->GetEntries();
  if(mccommand==2)  doptweight=0;
  TH2D::SetDefaultSumw2(true);
  c->LoadNoTrees();
  c->hasEvtTree = true;
  c->hasHltTree = true;
  c->hasTrackTree = true;
  c->hasAkPu3JetTree = true;
  c->hasSkimTree = true;
  

  // Long64_t nentries = c->GetEntries();
  // Long64_t nbytes = 0, nb = 0;
  TH2D * hJetTrackBackground;
  cout<<"before first exit"<<endl;
  if(jetindex==0) //leading jet
  {
    hJetTrackBackground = new TH2D(Form("background_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,75,-pi,2*pi);
  }
  else if(jetindex==1) //subleading jet
  {
    hJetTrackBackground = new TH2D(Form("background_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,75,-pi,2*pi);
  }
  else //leading jet
  {
    cout<<"Error: only jetindex 0 and 1 currently supported (leading + subleading jet)"<<endl;
    exit(1);
  }
  
  double ntottrig = 0;
  int n_entries_in_cent_range = cent_index_start[centmax] - cent_index_start[centmin];


  cout<<"before event for loop "<<n_entries_in_cent_range<<endl;
  // for (Long64_t jentry=0; jentry<10000;jentry++) {
  // for (Long64_t jentry=cent_index_start[centmin]; jentry<cent_index_start[centmax];jentry++) {
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    // if(jentry%1000==0) cout<<jentry-cent_index_start[centmin]<<"/"<<n_entries_in_cent_range<<endl;
    if(jentry%1000==0) cout<<jentry<<"/"<<nentries<<endl;
    c->GetEntry(jentry);
    if(c->evt.hiBin < centmin || c->evt.hiBin > centmax) continue;
    
    //! jet cuts, vz range, more than 2 jettss, dijetdphi 7/8 pi, both jet |eta| < 2, pt range cuts, aj cuts
    // cout<<mccommand<<endl;
    if(mccommand==0)
    {
      // if(!c->selectEvent()) continue; 
      // if(!c->hlt.HLT_HIJet80_v1) continue; 
      bool select = (c->skim.pHBHENoiseFilter || c->mc) && c->skim.pcollisionEventSelection;

      if(!select) continue; 
    }
    if( fabs(c->evt.vz) > vzrange ) continue;
    // if( c->akPu3PF.nref < 2 )         continue;
    double dijetdphi = fabs(c->akPu3PF.jtphi[0] - c->akPu3PF.jtphi[1]);
    if( dijetdphi > pi ) dijetdphi = 2*pi - dijetdphi;
    if( dijetdphi < dijetdphicut ) continue;
    // while(true)
    // {
    // jetindex++;
    if( c->akPu3PF.nref < 2 ) continue;
    if( c->akPu3PF.jtpt[0] > leadingjetpthigh ) continue;
    if( c->akPu3PF.jtpt[0] < leadingjetptlow ) continue;
    if( c->akPu3PF.jtpt[1] > subleadingjetpthigh ) continue;
    if( c->akPu3PF.jtpt[1] < subleadingjetptlow ) continue;
    
    if( fabs(c->akPu3PF.jteta[0]) > jetamax || fabs(c->akPu3PF.jteta[0]) < jetamin ) continue;
    if( fabs(c->akPu3PF.jteta[1]) > jetamax || fabs(c->akPu3PF.jteta[1]) < jetamin ) continue;
    if ((c->akPu3PF.trackMax[0]/c->akPu3PF.jtpt[0])<0.01) continue;
    if ((c->akPu3PF.trackMax[1]/c->akPu3PF.jtpt[1])<0.01) continue;
    // cout<<"hhere"<<endl;
    // cout<<"ihere"<<endl;
    // cout<<"ghere "<<fabs(c->akPu3PF.jteta[1])<<" > "<<jetamax<<" ... "<<fabs(c->akPu3PF.jteta[1])<<" < "<<jetamin<<endl;
    // cout<<"jhere"<<endl;
    // if(jentry%1000==0) { cout<<"here"<<endl; break;}
    
    double jeteta = c->akPu3PF.jteta[jetindex];
    double jetphi = c->akPu3PF.jtphi[jetindex];
    // cout<<jeteta<<endl;
    int ntrig = 1;
    int mixentry = jentry;
    int k = 0;
    int statfactor = 10;
    // for(int k = 0 ; k < statfactor ; ++k)
    while(k<statfactor)
    {
      mixentry++;
      if(mixentry > nentries) break;
      k++;
      c->GetEntry(mixentry);
      for(int i = 0 ; i < ntrig ; ++i)
      {
        ntottrig += 1;
        // cout<<nass<<endl;
        for(int j = 0 ; j < c->track.nTrk ; ++j)
        {
          // mytrackquality==0&&
          // cout<<"here front"<<endl;
          if(fabs(c->track.trkEta[j])<2.4&&c->track.highPurity[j]&&fabs(c->track.trkDz1[j]/c->track.trkDzError1[j])<3&&fabs(c->track.trkDxy1[j]/c->track.trkDxyError1[j])<3&&c->track.trkPtError[j]/c->track.trkPt[j]<0.1)
          {
            // cout<<"here back"<<endl;
            // cout<<asstrkIndex[k]<<" : "<<asstrkEta[asstrkIndex[k]]<<endl;
            double deta = fabs(jeteta-c->track.trkEta[j]);
            double dphi = fabs(jetphi-c->track.trkPhi[j]);
            // double ptw = asstrkPt[k];
            if( dphi > pi ) dphi = 2*pi - dphi;
            
            hjtrketa->Fill(c->track.trkEta[j]);
            hjtrkphi->Fill(c->track.trkPhi[j]);
            hjtrkdeta->Fill(deta);
            hjtrkdphi->Fill(dphi);
            
            double effweight = c->track.trkPt[j];
            // effweight = 1;
            // if(doptweight!=0) effweight = c->getTrackCorrection(j);
            hJetTrackBackground->Fill(deta,dphi,effweight);
            hJetTrackBackground->Fill(-deta,dphi,effweight);
            hJetTrackBackground->Fill(deta,-dphi,effweight);
            hJetTrackBackground->Fill(-deta,-dphi,effweight);
            hJetTrackBackground->Fill(deta,(2*pi)-dphi,effweight);
            hJetTrackBackground->Fill(-deta,(2*pi)-dphi,effweight);
          }
        }
      }
    }
    // }

    // if(jentry>100000) break;
  }
  // if( ntottrig != 0 ) hJetTrackBackground->Scale(1/ntottrig);
  c->ResetBooleans();
  return hJetTrackBackground;
}


int GetNTotTrig()
{
  int tmp = ntottrig;
  ntottrig = 0;
  return tmp;
}
