// #include "/net/hisrv0001/home/dav2105/run/CMSSW_4_4_4/src/CmsHi/JetAnalysis/macros/forest/hiForest_charge.h"
#include "HiForestAnalysis/hiForest.h"
// #include "hiForest.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>
#include "factorizedPtCorr.h"


using namespace std;

const double pi = 3.1415926;
const double pi23 = 2.0 * pi / 3.0;
double maxetatrg = 2.4;
double maxetaass = 2.4;
int ntottrig = 0;
int cent_index_start[41];
int mytrackquality = 0;
int sigcentarray[100000];
int bakcentarray[100000];
float sigvzarray[100000];
float bakvzarray[100000];
bool bakpcollarray[100000];
bool bakhbhearray[100000];
bool bakhltjet80[100000];
sampleType sType = kHIDATA;

TH1D * hcent;
TH1D * hcentpostcut;

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
TH1D * hjtrkpt;
TH1D * hjtrkdeta;
TH1D * hjtrkdphi;

HiForest *c;
HiForest *bk;
TH2D * JetTrackSignal(int condor_iter, int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int doptweight , int mccommand , double jetamin , double jetamax , double vzrange = 15.0, double dijetdphicut = pi23, string whichjet = "");
TH2D * JetTrackBackground(int condor_iter, int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int doptweight , int mccommand , double jetamin , double jetamax , double vzrange = 15.0, double dijetdphicut = pi23, string whichjet = "");

int GetNTotTrig();

void jettrkana(const char * infname = "/mnt/hadoop/cms/store/user/velicanu/mergedv1_sortedforest/mergesortv1_0.root", const char * bkfname = "/mnt/hadoop/cms/store/user/velicanu/mergedv1_sortedforest/mergesortv1_0.root", int trackquality = 0, string whichjet = "")
{
  mytrackquality = trackquality;
  cout<<"initializing hiforest and building centrality index"<<endl;
  cout<<"running on: "<<infname<<endl;
  // c = new HiForest(infname,"forest",cPbPb);

  // bk= new HiForest(bkfname,"forest2",cPPb,0,whichjet);
  c = new HiForest(infname,"forest1",cPPb,0,whichjet);
  bk = c;
  
  Long64_t nentries = c->GetEntries();
  // Long64_t nbkentries = bk->GetEntries();

  
  InitCorrFiles(sType);   // Yen-Jie: assuming that we are working on Heavy Ion Data!!!!!!!!
  cout <<"Loading HEAVY ION DATA correction tables !!"<<endl;
  InitCorrHists(sType);

  c->LoadNoTrees();
  c->hasSkimTree = true;
  c->hasEvtTree = true;
  cout<<"consturcting cent, vz, pcoll, hbhe arrays... ";
  for(int i = 0 ; i < nentries ; ++i)
  {
    c->GetEntry(i);
    sigcentarray[i]=c->evt.hiBin;
    sigvzarray[i]=c->evt.vz;
  // }
  // for(int i = 0 ; i < nbkentries ; ++i)
  // {
    // bk->GetEntry(i);
    bakcentarray[i]=bk->evt.hiBin;
    bakvzarray[i]=bk->evt.vz;
    bakpcollarray[i]=c->skim.pcollisionEventSelection;
    bakhbhearray[i]=c->skim.pHBHENoiseFilter;
    // bakhltjet80[i]=c->hlt.HLT_HIJet80_v1;
  }
  cout<<"done."<<endl;
  c->ResetBooleans();
}

bool skipevent(double vzrange)
{
  bool doskip = false;
  if(fabs(c->evt.vz)>vzrange) doskip = true;
  return doskip;
}

TH2D * JetTrackSignal(int condor_iter, int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int doptweight , int mccommand , double jetamin , double jetamax , double vzrange, double dijetdphicut, string whichjet)
{
  int parallelization = 20;
  cout<<"whichjet: "<<whichjet<<endl;
  Long64_t nentries = c->GetEntries();
  if(mccommand==2)  doptweight=0;
  TH2D::SetDefaultSumw2(true);
  c->LoadNoTrees();
  c->hasEvtTree = true;
  // c->hasHltTree = true;
  // c->hasTrackTree = true;
  // c->hasAkPu3JetTree = true;
  c->hasSkimTree = true;
  c->hasAkVs3CaloJetTree = true;
  
  hjeteta    = new TH1D(Form("hjet%d_eta_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",80,-4,4);
  hjetphi    = new TH1D(Form("hjet%d_phi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#phi",80,-TMath::Pi(),TMath::Pi());
  hjetpt    = new TH1D(Form("hjet%d_pt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";p_{T}",80,0,leadingjetpthigh);
  
  hcent    = new TH1D(Form("hcent%d_pt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";p_{T}",80,0,leadingjetpthigh);
  hcentpostcut    = new TH1D(Form("hcentpostcut%d_pt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";p_{T}",80,0,leadingjetpthigh);
  
  
  hjetdphi    = new TH1D(Form("hjet%d_dphi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#Delta#phi",40,0,TMath::Pi());

  hjtrketa    = new TH1D(Form("hj%d_trketa_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",80,-4,4);
  hjtrkphi    = new TH1D(Form("hj%d_trkphi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#phi",80,-TMath::Pi(),TMath::Pi());
  hjtrkpt    = new TH1D(Form("hj%d_trkpt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";p_{T}",80,0,ptasshigh*2);

  hjtrkdeta    = new TH1D(Form("hj%d_trkdeta_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",100,-5,5);
  hjtrkdphi    = new TH1D(Form("hj%d_trkdphi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#phi",80,0,TMath::Pi());


  // Long64_t nentries = c->GetEntries();
  // Long64_t nbytes = 0, nb = 0;
  TH2D * hJetTrackSignal;
  // cout<<"before first exit"<<endl;
  if(jetindex==0) //leading jet
  {
    //Hallie Binning
    // hJetTrackSignal = new TH2D(Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());
    // hJetTrackSignal = new TH2D(Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",50,-5,5,50,-2*pi,2*pi);
    
    //Pelin Binning
    hJetTrackSignal = new TH2D(Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,75,-1.0*TMath::Pi(),2.0*TMath::Pi());
    // hJetTrackSignal = new TH2D(Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",50,-5,5,50,-2*pi,2*pi);
  
  }
  else if(jetindex==1) //subleading jet
  {
    //Hallie Binning
    // hJetTrackSignal = new TH2D(Form("signal_%s_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());
    // hJetTrackSignal = new TH2D(Form("signal_%s_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",50,-5,5,50,-2*pi,2*pi);
    
    //Pelin Binning
    hJetTrackSignal = new TH2D(Form("signal_%s_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,75,-1.0*TMath::Pi(),2.0*TMath::Pi());
    // hJetTrackSignal = new TH2D(Form("signal_%s_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",50,-5,5,50,-2*pi,2*pi);
  
  }
  else //leading jet
  {
    cout<<"Error: only jetindex 0 and 1 currently supported (leading + subleading jet)"<<endl;
    exit(1);
  }
  // return hJetTrackSignal;
  
  double ntottrig = 0;
  // int n_entries_in_cent_range = cent_index_start[centmax] - cent_index_start[centmin];
// cout<<parallelization<<endl;
  int increment = nentries/parallelization;
  // int increment = nentries;
  // for (Long64_t jentry=0; jentry<10000;jentry++) {
  // for (Long64_t jentry=cent_index_start[centmin]; jentry<cent_index_start[centmax];jentry++) {
  for (Long64_t jentry=condor_iter*increment; jentry<(condor_iter+1)*increment;jentry++) {
    // if(jentry%1000==0) cout<<jentry-cent_index_start[centmin]<<"/"<<n_entries_in_cent_range<<endl;
    if(jentry%1000==0) cout<<jentry<<"/"<<nentries<<endl;
    c->GetEntry(jentry);
    hcent->Fill(c->evt.hiBin);
    

    // ### START Event Selection and quality cuts ###
    if(c->evt.hiBin < centmin || c->evt.hiBin >= centmax) continue;
    if(c->hasGenParticleTree)
    {
      //if is PbPb MC
      if(!c->skim.pcollisionEventSelection) continue;
    }
    else
    {
      //if is PbPb Data
      if(!(c->skim.pcollisionEventSelection && c->skim.pHBHENoiseFilter )) continue;
      // if(!(c->skim.pcollisionEventSelection && c->skim.pHBHENoiseFilter && c->hlt.HLT_HIJet80_v1 )) continue;
    }
    if( fabs(c->evt.vz) > vzrange ) continue;
    // ### END Event Selection and quality cuts ###
    

    



    
    // ### Find leading and subleading jets in |eta| < 2
    int leadindex = -1 , subleadindex = -1;
    bool foundleading = false;
    double oldjetpt = 999999999;
    for(int i = 0 ; i < c->myjet.nref ; ++i)
    {
      double deltajetpt = oldjetpt - c->myjet.jtpt[i];
      if(deltajetpt < 0) { cout<<"FATAL ERROR: unsorted jets"<<endl; exit(1); }
      oldjetpt = c->myjet.jtpt[i];
      if(!foundleading)
      {
        if(fabs(c->myjet.jteta[i])<2.0 )
        {
          leadindex = i;
          foundleading = true;
        }
      }
      else
      {
        if(fabs(c->myjet.jteta[i])<2.0 ) 
        {
          subleadindex = i;
          break;
        }
      }
    }
    if(leadindex < 0 || subleadindex < 0 ) continue; //Didn't find jet pair in |eta|<2
    
    // ### START et cuts ###
    double dijetdphi = fabs(c->myjet.jtphi[leadindex] - c->myjet.jtphi[subleadindex]);
    if( dijetdphi > pi ) dijetdphi = 2*pi - dijetdphi;
    if( fabs(dijetdphi) < dijetdphicut ) continue;  // dphi cut
    if( c->myjet.jtpt[leadindex] > leadingjetpthigh ) continue;       //lead pt cut
    if( c->myjet.jtpt[leadindex] < leadingjetptlow ) continue;        //lead pt cut
    if( c->myjet.jtpt[subleadindex] > subleadingjetpthigh ) continue; //sublead pt cut
    if( c->myjet.jtpt[subleadindex] < subleadingjetptlow ) continue;  //sublead pt cut
    if( fabs(c->myjet.jteta[leadindex]) > jetamax || fabs(c->myjet.jteta[leadindex]) < jetamin ) continue;        //jet eta cut
    if( fabs(c->myjet.jteta[subleadindex]) > jetamax || fabs(c->myjet.jteta[subleadindex]) < jetamin ) continue;  //jet eta cut
    if ((c->myjet.trackMax[leadindex]/c->myjet.jtpt[leadindex])<0.01) continue;         //trackmax cut
    if ((c->myjet.trackMax[subleadindex]/c->myjet.jtpt[subleadindex])<0.01) continue;   //trackmax cut
    float aj = (c->myjet.jtpt[leadindex] - c->myjet.jtpt[subleadindex]) / (c->myjet.jtpt[leadindex] + c->myjet.jtpt[subleadindex]);
    if(aj < ajmin || aj > ajmax) continue;  //aj cut
    // ### END jet cuts ###


    
    // ### Set which jet to correlate with , 0 = leading , 1 = subleadings
    int dojet = -1;
    if(jetindex==0) dojet = leadindex;
    else if(jetindex==1) dojet = subleadindex;
    else { cout<<"Error: in JetTrackSignal, only jetindex 0 and 1 supported"<<endl; exit(1); }
    
    
    
    
    double jeteta = c->myjet.jteta[dojet];
    double jetphi = c->myjet.jtphi[dojet];
    
    hjeteta->Fill(c->myjet.jteta[dojet]);
    hjetphi->Fill(c->myjet.jtphi[dojet]);
    hjetdphi->Fill(dijetdphi);
    hjetpt->Fill(c->myjet.jtpt[dojet]);
    ntottrig += 1;
    c->hasTrackTree = true;
    c->GetEntry(jentry);
    InitPosArrPbPb(c->evt.hiBin);
    for(int j = 0 ; j < c->track.nTrk ; ++j)
    {
      if(fabs(c->track.trkEta[j])<2.4&&c->track.highPurity[j]&&fabs(c->track.trkDz1[j]/c->track.trkDzError1[j])<3&&fabs(c->track.trkDxy1[j]/c->track.trkDxyError1[j])<3&&c->track.trkPtError[j]/c->track.trkPt[j]<0.1)
      {
        if(c->track.trkPt[j]>ptasshigh || c->track.trkPt[j]<ptasslow) continue;

        float correction = 1.0;
        float trkRMin = getTrkRMin(c->track.trkPhi[j], c->track.trkEta[j],c->myjet.nref,c->myjet.jtphi,c->myjet.jteta);
        float effweight = factorizedPtCorr(getPtBin(c->track.trkPt[j], sType), c->evt.hiBin, c->track.trkPt[j], c->track.trkPhi[j], c->track.trkEta[j], trkRMin, sType);
        
        // effweight = 1.0;
        double deta = -99, dphi = -99;
        
        if(false) //! Dragos:
        {
          deta = fabs(jeteta-c->track.trkEta[j]);
          dphi = fabs(jetphi-c->track.trkPhi[j]);
          if( dphi > pi ) dphi = 2*pi - dphi;
          // double effweight = c->track.trkPt[j]*correction;
          hJetTrackSignal->Fill(deta,dphi,effweight);
          hJetTrackSignal->Fill(-deta,dphi,effweight);
          hJetTrackSignal->Fill(deta,-dphi,effweight);
          hJetTrackSignal->Fill(-deta,-dphi,effweight);
          hJetTrackSignal->Fill(deta,(2*pi)-dphi,effweight);
          hJetTrackSignal->Fill(-deta,(2*pi)-dphi,effweight);
        }
        else //! Hallie
        {
          deta = jeteta - c->track.trkEta[j];
          dphi = jetphi - c->track.trkPhi[j];
          if(dphi>(1.5*TMath::Pi()))  { dphi+= -2*TMath::Pi(); }
          if(dphi<(-0.5*TMath::Pi())) { dphi+=  2*TMath::Pi(); }
          
          // hJetTrackSignal->Fill(deta,dphi, effweight);
          hJetTrackSignal->Fill(deta,dphi, effweight);
        }
        hjtrketa->Fill(c->track.trkEta[j],correction);
        hjtrkphi->Fill(c->track.trkPhi[j],correction);
        hjtrkpt->Fill(c->track.trkPt[j],correction);
        hjtrkdeta->Fill(deta);
        hjtrkdphi->Fill(dphi);

      }
    }
    c->hasTrackTree = false;
    hcentpostcut->Fill(c->evt.hiBin);
  }
  cout<<ntottrig<<endl;
  return hJetTrackSignal;
}


TH2D * JetTrackBackground(int condor_iter, int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int doptweight , int mccommand , double jetamin , double jetamax , double vzrange, double dijetdphicut, string whichjet)
{
  Long64_t nentries = c->GetEntries();
  Long64_t nbkentries = bk->GetEntries();
  
  int parallelization = 20;

  c->ResetBooleans();
  
  if(mccommand==2)  doptweight=0;
  TH2D::SetDefaultSumw2(true);
  bk->LoadNoTrees();
  // bk->hasTrackTree = true;
  bk->hasEvtTree = true;
  bk->hasSkimTree = true;

  c->hasAkVs3CaloJetTree = true;

  // bk->hasGenParticleTree = (bk->genParticleTree != 0);
  // c->hasEvtTree = true;
  // c->hasHltTree = true;
  // c->hasTrackTree = true;
  // c->hasAkPu3JetTree = true;
  // c->hasSkimTree = true;
  

  // Long64_t nentries = c->GetEntries();
  // Long64_t nbytes = 0, nb = 0;
  TH2D * hJetTrackBackground;
  cout<<"===================================================="<<endl;
  if(jetindex==0) //leading jet
  {
    //Hallie Binning
    // hJetTrackBackground = new TH2D(Form("background_%s_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());

    //Pelin Binning
    hJetTrackBackground = new TH2D(Form("background_%s_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,75,-1.0*TMath::Pi(),2.0*TMath::Pi());
  }
  else if(jetindex==1) //subleading jet
  {
    //Hallie Binning
    // hJetTrackBackground = new TH2D(Form("background_%s_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());

    //Pelin Binning
    hJetTrackBackground = new TH2D(Form("background_%s_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,75,-1.0*TMath::Pi(),2.0*TMath::Pi());


  }
  else //leading jet
  {
    cout<<"Error: only jetindex 0 and 1 currently supported (leading + subleading jet)"<<endl;
    exit(1);
  }

  // return hJetTrackBackground;
  
  double ntottrig = 0;
  int n_entries_in_cent_range = cent_index_start[centmax] - cent_index_start[centmin];


  cout<<"before event for loop "<<n_entries_in_cent_range<<endl;
  int sigindex = 0;
  int bakindex = 0;
  float ntries = 0.0;
  float njetsevts = 0.0;
  int old_bakindex = 0;
  // for (Long64_t jentry=0; jentry<10000;jentry++) {
  // for (Long64_t jentry=cent_index_start[centmin]; jentry<cent_index_start[centmax];jentry++) {
cout<<parallelization<<endl;


  int increment = nentries/parallelization;
  // int increment = nentries;

  for (Long64_t jentry=condor_iter*increment; jentry<(condor_iter+1)*increment;jentry++) {
    // if(jentry%1000==0) cout<<jentry-cent_index_start[centmin]<<"/"<<n_entries_in_cent_range<<endl;
    if(jentry%1000==0) cout<<jentry<<"/"<<nentries<<endl;
    // if(jentry>100) break;
    c->GetEntry(jentry);
    sigindex++;
    if(c->evt.hiBin < centmin || c->evt.hiBin > centmax) continue;
    
    //! jet cuts, vz range, more than 2 jettss, dijetdphi 7/8 pi, both jet |eta| < 2, pt range cuts, aj cuts
    // cout<<mccommand<<endl;
    if(c->hasGenParticleTree)
    {
      //if is PbPb MC
      if(!c->skim.pcollisionEventSelection) continue;
    }
    else
    {
      //if is PbPb Data
      if(!(c->skim.pcollisionEventSelection && c->skim.pHBHENoiseFilter )) continue;
      // if(!(c->skim.pcollisionEventSelection && c->skim.pHBHENoiseFilter && c->hlt.HLT_HIJet80_v1)) continue;
    }
    
    if( fabs(c->evt.vz) > vzrange ) continue;
    if( c->myjet.nref < 2 ) continue;
    int leadindex = -1 , subleadindex = -1;
    bool foundleading = false;
    for(int i = 0 ; i < c->myjet.nref ; ++i)
    {
      if(!foundleading)
      {
        if(fabs(c->myjet.jteta[i])<2.0 )
        {
          leadindex = i;
          foundleading = true;
        }
      }
      else
      {
        if(fabs(c->myjet.jteta[i])<2.0 ) 
        {
          subleadindex = i;
          break;
        }
      }
    }
    if(leadindex < 0 || subleadindex < 0 ) continue; //Didn't find jet pair in |eta|<2
    
    int dojet = -1;
    if(jetindex==0) dojet = leadindex;
    else if(jetindex==1) dojet = subleadindex;
    else { cout<<"Error: in JetTrackBackground, only jetindex 0 and 1 supported"<<endl; exit(1); }
    
    
    double dijetdphi = fabs(c->myjet.jtphi[0] - c->myjet.jtphi[1]);
    if( dijetdphi > pi ) dijetdphi = 2*pi - dijetdphi;
    if( fabs(dijetdphi) < dijetdphicut ) continue;
    // while(true)
    // {
    // jetindex++;
    if( c->myjet.jtpt[leadindex] > leadingjetpthigh ) continue;
    if( c->myjet.jtpt[leadindex] < leadingjetptlow ) continue;
    if( c->myjet.jtpt[subleadindex] > subleadingjetpthigh ) continue;
    if( c->myjet.jtpt[subleadindex] < subleadingjetptlow ) continue;
    
    if( fabs(c->myjet.jteta[leadindex]) > jetamax || fabs(c->myjet.jteta[leadindex]) < jetamin ) continue;
    if( fabs(c->myjet.jteta[subleadindex]) > jetamax || fabs(c->myjet.jteta[subleadindex]) < jetamin ) continue;
    if ((c->myjet.trackMax[leadindex]/c->myjet.jtpt[leadindex])<0.01) continue;
    if ((c->myjet.trackMax[subleadindex]/c->myjet.jtpt[subleadindex])<0.01) continue;
    // cout<<"hhere"<<endl;
    // cout<<"ihere"<<endl;
    // cout<<"ghere "<<fabs(c->myjet.jteta[1])<<" > "<<jetamax<<" ... "<<fabs(c->myjet.jteta[1])<<" < "<<jetamin<<endl;
    // cout<<"jhere"<<endl;
    // if(jentry%1000==0) { cout<<"here"<<endl; break;}
    
    double jeteta = c->myjet.jteta[dojet];
    double jetphi = c->myjet.jtphi[dojet];
    cout<<"jetetaphi: "<<jeteta<<" "<<jetphi<<endl;
    // int ntrig = 1;
    int k = 0;
    int statfactor = 10;
    int thisntries = 0;
    // for(int k = 0 ; k < statfactor ; ++k)
    // bool first = true;
    // bakindex = old_bakindex;
    bakindex = jentry;
    old_bakindex = jentry-1;
    njetsevts+=1.0;
    while(k<statfactor)
    {
      ntries+=1.0;
      thisntries++;
      bakindex++;
      if(bakindex > nbkentries) bakindex = 0;
      if(bakindex==jentry) continue;
      if(old_bakindex==bakindex)
      {
        cout<<"wraparound, got "<<k<<" : cent,vz = "<<sigcentarray[jentry]<<","<<sigvzarray[jentry]<<endl;
        break; //we wrapped around
      }
      // cout<<"a"<<endl;
      if(abs(sigcentarray[jentry]-bakcentarray[bakindex])/(double)sigcentarray[jentry]>0.1) continue;
      if(fabs(sigvzarray[jentry] - bakvzarray[bakindex])>5) continue;
      // cout<<abs(sigcentarray[jentry]-bakcentarray[bakindex])<<"/"<<(double)sigcentarray[jentry]<<" = "<<
      // abs(sigcentarray[jentry]-bakcentarray[bakindex])/(double)sigcentarray[jentry]<<" and |"<<
      // sigvzarray[jentry]<<" - "<<bakvzarray[bakindex]<<"| = "<<fabs(sigvzarray[jentry] - bakvzarray[bakindex])<<" < "<<5<<endl;
      
      // cout<<"b"<<endl;
      //
      if((bk->genParticleTree != 0))
      {
        //if is PbPb MC
        if(!bakpcollarray[bakindex]) continue;
      }
      else
      {
        //if is PbPb Data
        if(!(bakpcollarray[bakindex] && bakhbhearray[bakindex] )) continue;
        // if(fabs(bakvzarray[]))
      }
      // cout<<"c"<<endl; 
      bk->hasTrackTree = true;
      bk->GetEntry(bakindex);
      // cout<<"a"<<endl;
      k++;

      // for(int i = 0 ; i < ntrig ; ++i)
      // {
      ntottrig += 1;
      InitPosArrPbPb(c->evt.hiBin);
      for(int j = 0 ; j < bk->track.nTrk ; ++j)
      {
        if(fabs(bk->track.trkEta[j])<2.4&&bk->track.highPurity[j]&&fabs(bk->track.trkDz1[j]/bk->track.trkDzError1[j])<3&&fabs(bk->track.trkDxy1[j]/bk->track.trkDxyError1[j])<3&&bk->track.trkPtError[j]/bk->track.trkPt[j]<0.1)
        {
          if(bk->track.trkPt[j]>ptasshigh || bk->track.trkPt[j]<ptasslow) continue;
          double deta = fabs(jeteta-bk->track.trkEta[j]);
          // cout<<"aq"<<endl;
          double dphi = fabs(jetphi-bk->track.trkPhi[j]);
          if( dphi > pi ) dphi = 2*pi - dphi;
          
          // float correction = c->getTrackCorrection(j);
          // correction = 1.0;
          float trkRMin = getTrkRMin(c->track.trkPhi[j], c->track.trkEta[j],c->myjet.nref,c->myjet.jtphi,c->myjet.jteta);
          float effweight = factorizedPtCorr(getPtBin(c->track.trkPt[j], sType), c->evt.hiBin, c->track.trkPt[j], c->track.trkPhi[j], c->track.trkEta[j], trkRMin, sType);
          // effweight = 1.0;
          
          hJetTrackBackground->Fill(deta,dphi,effweight);
          hJetTrackBackground->Fill(-deta,dphi,effweight);
          hJetTrackBackground->Fill(deta,-dphi,effweight);
          hJetTrackBackground->Fill(-deta,-dphi,effweight);
          hJetTrackBackground->Fill(deta,(2*pi)-dphi,effweight);
          hJetTrackBackground->Fill(-deta,(2*pi)-dphi,effweight);
        }
      }
      bk->hasTrackTree = false;
      // }
      // cout<<"done loop"<<endl;
    }
    // }

    // if(jentry>100000) break;
  }
  cout<<"done back: "<<ntottrig<<endl;
  cout<<"njetsevts: "<<njetsevts<<" , ntries: "<<ntries<<" , efficiency: "<<njetsevts/ntries<<endl;
  // if( ntottrig != 0 ) hJetTrackBackground->Scale(1/ntottrig);
  bk->ResetBooleans();
  cout<<"-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"<<endl;
  return hJetTrackBackground;
}


int GetNTotTrig()
{
  int tmp = ntottrig;
  ntottrig = 0;
  return tmp;
}
