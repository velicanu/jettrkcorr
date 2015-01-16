// #include "/net/hisrv0001/home/dav2105/run/CMSSW_4_4_4/src/CmsHi/JetAnalysis/macros/forest/hiForest_charge.h"
#include "HiForestAnalysis/hiForest.h"
// #include "hiForest.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>
#include "factorizedPtCorr.h"
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <iomanip>  
#include "array2D.C"
// #include <limits>


using namespace std;

// const double pi = 3.141592653589793238462643;
const double pi = TMath::Pi();
const double pi56 = 5.0 * pi / 6.0;
double maxetatrg = 2.4;
double maxetaass = 2.4;
int ntottrig = 0;
int cent_index_start[41];
int mytrackquality = 0;
int sigcentarray[100000];
int bakcentarray[987198];
float sigvzarray[100000];
float bakvzarray[100000];
bool bakpcollarray[987198];
bool bakhbhearray[987198];
bool bakhltjet80[987198];
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
TH1D * hjtrketanocorr;
TH1D * hjtrkphi;
TH1D * hjtrkpt;
TH1D * hjtrkdeta;
TH1D * hjtrkdphi;

TH1D * hsigvz;
TH1D * hsigcent;
TH1D * hsigvz2;
TH1D * hmixvz;
TH1D * hmixdvz;
TH1D * hmixtrkpt;
TH1D * hmixtrketa;
TH1D * hmixtrketanocorr;
TH1D * hmixtrkphi;
TH1D * hmixcent;
TH1D * hmixdcent; 

TH1D * hdvz;
TH1D * hnmatched;
int njet;
std::unordered_set<long long int> visitedevents;
std::vector<std::unordered_set<long long int> > pairs;
// std::unordered_set< std::pair<int,long long int> > mixevents;
std::vector<int> vjetindices;
std::vector<long long int> backgroundids;

HiForest *c;
HiForest *bk;
TH2D * JetTrackSignal(int condor_iter, int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int dotrkcorr , int mccommand , double jetamin , double jetamax , double vzrange = 15.0, double dijetdphicut = pi56, string whichjet = "");
TH2D * JetTrackBackground(int condor_iter, int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int dotrkcorr , int mccommand , double jetamin , double jetamax , double vzrange = 15.0, double dijetdphicut = pi56, string whichjet = "", double dvz = 0.2);

int GetNTotTrig();

bool rejectevent(int centmin, int centmax, double vzrange, double dijetdphicut, double leadingjetpthigh, double leadingjetptlow, double subleadingjetpthigh, double subleadingjetptlow, double jetamax, double jetamin, double ajmin, double ajmax) //all cuts go here
{
  bool willreject = true;
  for (int pseudoindex = 0; pseudoindex < 1; ++pseudoindex)
  {
    
    // ### START Event Selection and quality cuts ###
    // if(c->evt.hiBin < centmin || c->evt.hiBin >= centmax) continue;
    if(c->hasGenParticleTree)
    {
      //if is PbPb MC
      // if(!c->skim.pcollisionEventSelection) continue;
    }
    else
    {
      //if is PbPb Data
      // if(!(c->skim.pcollisionEventSelection && c->skim.pHBHENoiseFilter )) continue;
      // if(!(c->skim.pcollisionEventSelection && c->skim.pHBHENoiseFilter && c->hlt.HLT_HIJet80_v1 )) continue;
    }
    // if( fabs(c->evt.vz) > vzrange ) continue;
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
    
    // ### START jet cuts ###
    double dijetdphi = fabs(c->myjet.jtphi[leadindex] - c->myjet.jtphi[subleadindex]);
    if( dijetdphi > pi ) dijetdphi = 2*pi - dijetdphi;
    if( fabs(dijetdphi) < dijetdphicut ) continue;  // dphi cut
    if( c->myjet.jtpt[leadindex] > leadingjetpthigh ) continue;       //lead pt cut
    if( c->myjet.jtpt[leadindex] < leadingjetptlow ) continue;        //lead pt cut
    if( c->myjet.jtpt[subleadindex] > subleadingjetpthigh ) continue; //sublead pt cut
    if( c->myjet.jtpt[subleadindex] < subleadingjetptlow ) continue;  //sublead pt cut
    if( fabs(c->myjet.jteta[leadindex]) > jetamax || fabs(c->myjet.jteta[leadindex]) < jetamin ) continue;        //jet eta cut
    if( fabs(c->myjet.jteta[subleadindex]) > jetamax || fabs(c->myjet.jteta[subleadindex]) < jetamin ) continue;  //jet eta cut
    float aj = (c->myjet.jtpt[leadindex] - c->myjet.jtpt[subleadindex]) / (c->myjet.jtpt[leadindex] + c->myjet.jtpt[subleadindex]);
    if(aj < ajmin || aj > ajmax) continue;  //aj cut
    // ### END jet cuts ###
    willreject = false;
  }
  return willreject;
}

void jettrkana(const char * infname = "/mnt/hadoop/cms/store/user/velicanu/mergedv1_sortedforest/mergesortv1_0.root", const char * bkfname = "/mnt/hadoop/cms/store/user/velicanu/mergedv1_sortedforest/mergesortv1_0.root", int trackquality = 0, string whichjet = "")
{
  njet = 0;
  mytrackquality = trackquality;
  cout<<"initializing hiforest and building centrality index"<<endl;
  cout<<"running on: "<<infname<<endl;
  // c = new HiForest(infname,"forest",cPbPb);

  // bk= new HiForest(bkfname,"forest2",cPPb,0,whichjet);
  c = new HiForest(infname,"forest1",cPP,0,whichjet);
  // bk = c;
  
  Long64_t nentries = c->GetEntries();

  
  //InitCorrFiles(sType);   // Yen-Jie: assuming that we are working on Heavy Ion Data!!!!!!!!
  // cout <<"Loading HEAVY ION DATA correction tables !!"<<endl;
  //InitCorrHists(sType);

  // ### Some jet cuts ###
  double leadingjetpthigh = 300;
  double leadingjetptlow = 120;
  double subleadingjetpthigh = 300;
  double subleadingjetptlow = 50;
  double jetamax = 1.6;
  double jetamin = 0.0;
  int centmin = 0;
  int centmax = 20;
  int vzrange = 15;
  double dijetdphicut = 5 * TMath::Pi() / 6.0;
  double ajmin = 0.0;
  double ajmax = 1.0;

  // ### End jet cuts  ###

  // c->LoadNoTrees();
  // c->hasSkimTree = true;
  // c->hasEvtTree = true;
  // c->hasAk3CaloJetTree = true;
  cout<<"consturcting cent, vz, pcoll, hbhe arrays... ";
  for(int i = 0 ; i < nentries ; ++i)
  {
  break;
    c->GetEntry(i);
    if(i%1000==0) cout<<i<<"/"<<nentries<<endl;
    sigcentarray[i] = c->evt.hiBin;
    sigvzarray[i] = c->evt.vz;
    if(rejectevent(centmin,centmax,vzrange,dijetdphicut,leadingjetpthigh,leadingjetptlow,subleadingjetpthigh,subleadingjetptlow,jetamax,jetamin,ajmin,ajmax)) continue;
    vjetindices.push_back(i);
    njet++;
  }
  // cout<<"back init loop "<<endl;
  cout<<"njet = "<<njet<<endl;
  bk = new HiForest(bkfname,"forest2",cPPb,0,whichjet);
  Long64_t nbkentries = bk->GetEntries();
  // bk->LoadNoTrees();
  // bk->hasEvtTree = true;
  // bk->hasSkimTree = true;
  // bk->hasTrackTree = true;
  for(int i = 0 ; i < nbkentries ; ++i)
  {
  break;
		if(i%10000==0) cout<<i<<"/"<<nbkentries<<endl;
    bk->GetEntry(i);
    bakcentarray[i]=bk->evt.hiBin;
    bakvzarray[i]=bk->evt.vz;
    bakpcollarray[i]=bk->skim.pcollisionEventSelection;
    bakhbhearray[i]=bk->skim.pHBHENoiseFilter;
    long long int big = 1000000000;
    long long int thisid = (big*bk->evt.run) + bk->evt.evt; // format is run|0|evt 
    backgroundids.push_back(thisid);
    
    // cout<<bk->track.nTrk<<endl;
    // bakhltjet80[i]=c->hlt.HLT_HIJet80_v1;
  }
  // cout<<"done."<<endl;
  c->ResetBooleans();
  bk->ResetBooleans();




}

TH2D * JetTrackSignal(int condor_iter, int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int dotrkcorr , int mccommand , double jetamin , double jetamax , double vzrange, double dijetdphicut, string whichjet)
{
  // ofstream myfile;
  // myfile.open ("jetinfo.csv");
  
  
  c->ResetBooleans();
  // int parallelization = 40;
  int parallelization = 1;
  cout<<"whichjet: "<<whichjet<<endl;
  Long64_t nentries = c->GetEntries();
  if(mccommand==2)  dotrkcorr=0;
  TH2D::SetDefaultSumw2(true);
  // c->LoadNoTrees();
  // c->hasEvtTree = true;
  // c->hasHltTree = true;
  // c->hasTrackTree = true;
  // c->hasAkPu3JetTree = true;
  // c->hasSkimTree = true;
  // c->hasAk3CaloJetTree = true;
  
  hjeteta    = new TH1D(Form("hjet%d_eta_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",80,-4,4);
  hjetphi    = new TH1D(Form("hjet%d_phi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#phi",80,-TMath::Pi(),TMath::Pi());
  hjetpt    = new TH1D(Form("hjet%d_pt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";p_{T}",80,0,leadingjetpthigh);
  
  hcent    = new TH1D(Form("hcent%d_pt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";p_{T}",80,0,leadingjetpthigh);
  hcentpostcut    = new TH1D(Form("hcentpostcut%d_pt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";p_{T}",80,0,leadingjetpthigh);
  
  
  hjetdphi    = new TH1D(Form("hjet%d_dphi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#Delta#phi",40,0,TMath::Pi());

  hjtrketa    = new TH1D(Form("hj%d_trketa_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",100,-2.5,2.5);
  hjtrketanocorr    = new TH1D(Form("hj%d_trketanocorr_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",100,-2.5,2.5);
  hjtrkphi    = new TH1D(Form("hj%d_trkphi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#phi",80,-TMath::Pi(),TMath::Pi());
  hjtrkpt    = new TH1D(Form("hj%d_trkpt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";p_{T}",80,0,ptasshigh*2);

  hjtrkdeta    = new TH1D(Form("hj%d_trkdeta_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",100,-5,5);
  hjtrkdphi    = new TH1D(Form("hj%d_trkdphi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#phi",80,0,TMath::Pi());

  hsigvz = new TH1D("hsigvz",";vz;",100,-16,16);
  hsigcent = new TH1D("hsigcent",";vz;",200,0,200);

  // Long64_t nentries = c->GetEntries();
  // Long64_t nbytes = 0, nb = 0;
  TH2D * hJetTrackSignal;
  // cout<<"before first exit"<<endl;
  if(jetindex==0) //leading jet
  {
    //Hallie Binning
    hJetTrackSignal = new TH2D(Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());
    
    //Pelin Binning
    // hJetTrackSignal = new TH2D(Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,75,-1.0*TMath::Pi(),2.0*TMath::Pi());

  
  }
  else if(jetindex==1) //subleading jet
  {
    //Hallie Binning
    hJetTrackSignal = new TH2D(Form("signal_%s_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());
    
    //Pelin Binning
    // hJetTrackSignal = new TH2D(Form("signal_%s_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,75,-1.0*TMath::Pi(),2.0*TMath::Pi());
  
  }
  else //leading jet
  {
    cout<<"Error: only jetindex 0 and 1 currently supported (leading + subleading jet)"<<endl;
    exit(1);
  }
  // return hJetTrackSignal;
  
  double ntottrig = 0;
  // nentries = (int)vjetindices.size();
  int increment = nentries/parallelization;
  int stopentry = (condor_iter+1)*increment;
  if(stopentry > nentries)
  {
    cout<<"FATAL ERROR: trying to run on more events than we have"<<endl; exit(1);
  }
  if((condor_iter+2)*increment > nentries) stopentry = nentries;

  cout<<increment<<endl;
  cout<<"start: "<<condor_iter*increment<<" to end "<<stopentry<<" of "<<nentries<<endl;
  // myfile << std::setprecision(8)<< "Event Number,Jet pT,Jet phi,Jet eta,centrality,vz"<<endl;

  for (Long64_t jentry=condor_iter*increment; jentry<stopentry;jentry++) {
    c->GetEntry(jentry);
    if(jentry>10) break;
    if(jentry%1000==0) cout<<jentry<<"/"<<nentries<<endl;
    // long long int big = 1000000000;
    // long long int thisid = (big*c->evt.run) + c->evt.evt; // format is run|0|evt 
    // auto search = visitedevents.find(thisid);
    // if(search != visitedevents.end()) {
      // cout<<"this data sample has duplicate events :( , but we're not analyzing them :) "<<endl;
      // continue; // found duplicate
    // }
    // else // no duplicate found, add this to visited events
    // {
      // visitedevents.insert(thisid);
    // }

    // if(rejectevent(centmin,centmax,vzrange,dijetdphicut,leadingjetpthigh,leadingjetptlow,subleadingjetpthigh,subleadingjetptlow,jetamax,jetamin,ajmin,ajmax)) continue;


    // ### Find leading and subleading jets in |eta| < 2
    int leadindex = -1 , subleadindex = -1;
    bool foundleading = false;
    double oldjetpt = 999999999;
    cout<<c->myjet.nref<<endl;
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
    cout<<"here"<<endl;
    // ### Set which jet to correlate with , 0 = leading , 1 = subleadings
    int dojet = -1;
    if(jetindex==0) dojet = leadindex;
    else if(jetindex==1) dojet = subleadindex;
    else { cout<<"Error: in JetTrackSignal, only jetindex 0 and 1 supported"<<endl; exit(1); }
     
    double jeteta = c->myjet.jteta[dojet];
    double jetphi = c->myjet.jtphi[dojet];
    
    hjeteta->Fill(c->myjet.jteta[dojet]);
    hjetphi->Fill(c->myjet.jtphi[dojet]);
    double dijetdphi = fabs(c->myjet.jtphi[leadindex] - c->myjet.jtphi[subleadindex]);
    if( dijetdphi > pi ) dijetdphi = 2*pi - dijetdphi;
    hjetdphi->Fill(dijetdphi);
    hjetpt->Fill(c->myjet.jtpt[dojet]);
    ntottrig += 1;
    c->hasTrackTree = true;
    c->GetEntry(vjetindices[jentry]);
    InitPosArrPbPb(c->evt.hiBin);
    hsigvz->Fill(c->evt.vz);
    hsigcent->Fill(c->evt.hiBin);
    for(int j = 0 ; j < c->track.nTrk ; ++j)
    {
      if(fabs(c->track.trkEta[j])<2.4&&c->track.highPurity[j]&&fabs(c->track.trkDz1[j]/c->track.trkDzError1[j])<3&&fabs(c->track.trkDxy1[j]/c->track.trkDxyError1[j])<3&&c->track.trkPtError[j]/c->track.trkPt[j]<0.1)
      {
        if(c->track.trkPt[j]>ptasshigh || c->track.trkPt[j]<ptasslow) continue;

        float trkRMin = getTrkRMin(c->track.trkPhi[j], c->track.trkEta[j],c->myjet.nref,c->myjet.jtphi,c->myjet.jteta,c->myjet.jtpt);
        float effweight = factorizedPtCorr(getPtBin(c->track.trkPt[j], sType), c->evt.hiBin, c->track.trkPt[j], c->track.trkPhi[j], c->track.trkEta[j], trkRMin, sType);
        if(dotrkcorr==0) effweight = 1.0; //turns off tracking correction

        double deta = -99, dphi = -99;
        // cout<<c->track.trkPt[j]<<","<<c->track.trkEta[j]<<","<<c->track.trkPhi[j]<<","<<trkRMin<<","<<effweight<<endl;
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
        hjtrketa->Fill(c->track.trkEta[j],effweight);
        hjtrketanocorr->Fill(c->track.trkEta[j]);
        hjtrkphi->Fill(c->track.trkPhi[j],effweight);
        hjtrkpt->Fill(c->track.trkPt[j],effweight);
        hjtrkdeta->Fill(deta);
        hjtrkdphi->Fill(dphi);

      }
    }
    // exit(1);
    c->hasTrackTree = false;
    hcentpostcut->Fill(c->evt.hiBin);
  }
  cout<<ntottrig<<endl;
  // myfile.close();
  return hJetTrackSignal;
}


TH2D * JetTrackBackground(int condor_iter, int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int dotrkcorr , int mccommand , double jetamin , double jetamax , double vzrange, double dijetdphicut, string whichjet, double dvz)
{

  //### Start Histogram Definitions ###
  hmixvz = new TH1D("hmixvz",";vz;",100,-16,16);
  hsigvz2 = new TH1D("hsigvz2",";vz;",100,-16,16);
  hmixdvz = new TH1D("hmixdvz",";#Deltavz;",100,-1,1);
  hmixtrkpt = new TH1D("hmixtrkpt",";trkPt;",80,0,ptasshigh*2);
  hmixtrketa = new TH1D("hmixtrketa",";trkEta;",100,-2.5,2.5);
  hmixtrketanocorr = new TH1D("hmixtrketanocorr",";trkEta;",100,-2.5,2.5);
  hmixtrkphi = new TH1D("hmixtrkphi",";trkPhi;",80,-TMath::Pi(),TMath::Pi());
  hmixcent = new TH1D("hmixcent",";hiBin;",200,0,200);
  hmixdcent = new TH1D("hmixdcent",";#DeltahiBin;",5,0,5);
  //### End Histogram Definitions   ###
  
  Long64_t nentries = c->GetEntries();
  Long64_t nbkentries = bk->GetEntries();
  TString csvname;
  // if(condor_iter<10)
    // csvname = Form("mixevents_0%d.csv",condor_iter);
  // else
    // csvname = Form("mixevents_%d.csv",condor_iter);
  std::vector< std::vector< long long int > > mixevents = csvtoarraylong(csvname.Data()); 
  // cout<<(int)mixevents.size()<<" x "<<(int)mixevents[0].size()<<endl;
  

  int parallelization = 40;
  // int parallelization = 1;

  c->LoadNoTrees();
  bk->LoadNoTrees();
  
  if(mccommand==2)  dotrkcorr=0;
  TH2D::SetDefaultSumw2(true);
  // bk->hasTrackTree = true;
  c->hasEvtTree = true;
  c->hasSkimTree = true;
  c->hasAk3CaloJetTree = true;
  bk->hasEvtTree = true;
  bk->hasSkimTree = true;
  bk->hasAk3CaloJetTree = true;
  bk->hasTrackTree = true;

  // bk->hasGenParticleTree = (bk->genParticleTree != 0);
  // bk->hasGenParticleTree = (bk->genParticleTree != 0);
  // c->hasEvtTree = true;
  // c->hasHltTree = true;
  // c->hasTrackTree = true;
  // c->hasAkPu3JetTree = true;
  // c->hasSkimTree = true;
  hdvz    = new TH1D(Form("hdvz%d_eta_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",80,-5,5);
  hnmatched    = new TH1D(Form("hnmatched%d_eta_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",21,0,21);
  

  // Long64_t nentries = c->GetEntries();
  // Long64_t nbytes = 0, nb = 0;
  TH2D * hJetTrackBackground;
  cout<<"===================================================="<<endl;
  if(jetindex==0) //leading jet
  {
    //Hallie Binning
    hJetTrackBackground = new TH2D(Form("background_%s_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());

    //Pelin Binning
    // hJetTrackBackground = new TH2D(Form("background_%s_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,75,-1.0*TMath::Pi(),2.0*TMath::Pi());
  }
  else if(jetindex==1) //subleading jet
  {
    //Hallie Binning
    hJetTrackBackground = new TH2D(Form("background_%s_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());

    //Pelin Binning
    // hJetTrackBackground = new TH2D(Form("background_%s_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,75,-1.0*TMath::Pi(),2.0*TMath::Pi());


  }
  else //leading jet
  {
    cout<<"Error: only jetindex 0 and 1 currently supported (leading + subleading jet)"<<endl;
    exit(1);
  }

  return hJetTrackBackground;
  
  double ntottrig = 0;
  int n_entries_in_cent_range = cent_index_start[centmax] - cent_index_start[centmin];


  cout<<"before event for loop "<<n_entries_in_cent_range<<endl;
  int sigindex = 0;
  float ntries = 0.0;
  float njetsevts = 0.0;
  nentries = (int)vjetindices.size();
  int increment = nentries/parallelization;
  int stopentry = (condor_iter+1)*increment;
  if(stopentry > nentries)
  {
    cout<<"FATAL ERROR: trying to run on more events than we have"<<endl; exit(1);
  }
  if((condor_iter+2)*increment > nentries) stopentry = nentries;

  cout<<increment<<endl;
  cout<<"start: "<<condor_iter*increment<<" to end "<<stopentry<<" of "<<nentries<<endl;

  // for (int i = 0; i < (stopentry-(condor_iter*increment)); ++i)
  // {
  //   pairs.push_back(std::unordered_set<long long int> tempset);
  // }
  // int nextindex = 0;
  // for (int i = 0; i < (int)mixevents.size(); ++i)
  // {
  //   if(previndex!=mixevents[i][0])
  //   {
  //     previndex = mixevents[i][0];
  //     whichset++;
  //   }
  //   pairs[whichset].insert(mixevents[i][1]);
  // }
  // exit(1);
	int currentbackindex = 0;
  for (Long64_t jentry=condor_iter*increment; jentry<stopentry;jentry++) 
  {
    // if(jentry%100==0) 
		cout<<jentry<<"/"<<nentries<<endl;
    c->GetEntry(vjetindices[jentry]);
    // cout<<mixevents[nextindex][0]<<" "<<vjetindices[jentry]<<endl;
    std::unordered_set<long long int> tempset;
    // int tempindex = nextindex;
    // while(true)
    // {
      // if(tempindex == (int)mixevents.size()) break;
      // if(vjetindices[jentry]!=mixevents[tempindex][0])
      // {
      //   nextindex = tempindex;
      //   break;
      // }
      // cout<<"inserting: mixevents["<<tempindex<<"][1] = "<<mixevents[tempindex][1]<<endl;
      // tempset.insert(mixevents[tempindex][1]);
      // tempindex++;
    // }
    // cout<<vjetindices[jentry]<<" "<<mixevents[tempindex-1][0]<<" "<<mixevents[tempindex-1][1]<<endl;
    // continue;
    if(rejectevent(centmin,centmax,vzrange,dijetdphicut,leadingjetpthigh,leadingjetptlow,subleadingjetpthigh,subleadingjetptlow,jetamax,jetamin,ajmin,ajmax)) continue;
    sigindex++;

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
    
    
    double jeteta = c->myjet.jteta[dojet];
    double jetphi = c->myjet.jtphi[dojet];
    int k = 0;
    int statfactor = 20;
    njetsevts+=1.0;
    hsigvz2->Fill(c->evt.vz);
    // for(int backindex = 0 ; backindex < nbkentries ; ++backindex)
		
    for(int backindex = 0 ; backindex < nbkentries ; ++backindex)
    {
      if(k>=statfactor) break;
      currentbackindex = ( currentbackindex + 1 ) % nbkentries;
			ntries+=1.0;
      if(abs((sigcentarray[vjetindices[jentry]]-bakcentarray[currentbackindex]))>5) continue;
      // if(abs((sigcentarray[vjetindices[jentry]]-bakcentarray[currentbackindex]))>0) continue;
      if(fabs(sigvzarray[vjetindices[jentry]] - bakvzarray[currentbackindex])>dvz) continue;
      // cout<<"searching for: "<<backgroundids[backindex]<<endl;
      
			//! fancy part turned off
      // auto search = tempset.find(backgroundids[backindex]);
			// if(search == tempset.end()) continue;
      bk->GetEntry(currentbackindex);

      if(c->evt.run == bk->evt.run && c->evt.evt == bk->evt.evt) continue; //same event
      if((bk->genParticleTree != 0))
      {
        //if is PbPb MC
        if(!bakpcollarray[currentbackindex]) continue;
      }
      else
      {
        //if is PbPb Data
        if(!(bakpcollarray[currentbackindex] && bakhbhearray[currentbackindex] )) continue;
      }
      k++;
      ntottrig += 1;
      InitPosArrPbPb(bk->evt.hiBin);
      hmixvz->Fill(bk->evt.vz);
      hmixdvz->Fill(c->evt.vz - bk->evt.vz);
      hmixcent->Fill(bk->evt.hiBin);
      hmixdcent->Fill(c->evt.hiBin - bk->evt.hiBin); 
      
      // cout<<bk->evt.run<<","<<bk->evt.evt<<","<<bk->evt.hiBin<<","<<bk->evt.vz<<endl;
      int thisntracks = 0;
      // cout<<"=============="<<endl;
      // cout<<"trkPt[j]"<<","<<"trkPhi[j]"<<","<<"trkEta[j]"<<","<<"effweight"<<endl;
      for(int j = 0 ; j < bk->track.nTrk ; ++j)
      {
        if(fabs(bk->track.trkEta[j])<2.4&&bk->track.highPurity[j]&&fabs(bk->track.trkDz1[j]/bk->track.trkDzError1[j])<3&&fabs(bk->track.trkDxy1[j]/bk->track.trkDxyError1[j])<3&&bk->track.trkPtError[j]/bk->track.trkPt[j]<0.1)
        {
          if(bk->track.trkPt[j]>ptasshigh || bk->track.trkPt[j]<ptasslow) continue;
		      double deta = jeteta - bk->track.trkEta[j];
          double dphi = jetphi - bk->track.trkPhi[j];
          if(dphi>(1.5*TMath::Pi()))  { dphi+= -2*TMath::Pi(); }
          if(dphi<(-0.5*TMath::Pi())) { dphi+=  2*TMath::Pi(); }
          float trkRMin = getTrkRMin(bk->track.trkPhi[j], bk->track.trkEta[j],bk->myjet.nref,bk->myjet.jtphi,bk->myjet.jteta,bk->myjet.jtpt);
          float effweight = factorizedPtCorr(getPtBin(bk->track.trkPt[j], sType), bk->evt.hiBin, bk->track.trkPt[j], bk->track.trkPhi[j], bk->track.trkEta[j], trkRMin, sType);
          if(dotrkcorr==0) effweight = 1.0; //turns off tracking correction
          hmixtrkpt->Fill(bk->track.trkPt[j],effweight);
          hmixtrketa->Fill(bk->track.trkEta[j],effweight);
          hmixtrketanocorr->Fill(bk->track.trkEta[j]);
          hmixtrkphi->Fill(bk->track.trkPhi[j],effweight);
          hJetTrackBackground->Fill(deta,dphi,effweight);
          // cout<<bk->track.trkPt[j]<<","<<bk->track.trkPhi[j]<<","<<bk->track.trkEta[j]<<","<<effweight<<endl;
          thisntracks++;
        }
      }
      // cout<<"number of tracks = " <<thisntracks<<endl;
      // exit(1);

      // bk->hasTrackTree = false;
    }
    // cout<<k<<endl;
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
