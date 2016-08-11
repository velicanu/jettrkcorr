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
const unsigned int JETSIZE = 2;
const double pi = TMath::Pi();
const double pi56 = 5.0 * pi / 6.0;
double maxetatrg = 2.4;
double maxetaass = 2.4;
int ntottrig = 0;
int cent_index_start[41];
int mytrackquality = 0;
int sigcentarray[444444];
int bakcentarray[444444];
float sigvzarray[444444];
float bakvzarray[444444];
bool bakpcollarray[444444];
bool bakhbhearray[444444];
bool bakhltjet80[444444];
// sampleType sType = kHIDATA;
sampleType sType = kPPDATA;
bool hasEVT = true;

TH1D * hcuts;
TH1D * hnjets;
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
TH2D * hjettrk;
TH2D * JetTrackSignalShifted;
TH2D * JetTrackBackgroundShifted;
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
TH1D * halljetpt;

TH1D * hjetdphi;

TH1D * hjtrketa;
TH1D * hjtrketanocorr;
TH1D * hjtrkphi;
TH1D * hjtrkpt;
TH1D * halltrkpt;
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
TH2D * hgenrefeta;
TH2D * hgenrefphi;
TH2D * hgenrefpt;
TH1D * hgenrefdeta;
TH1D * hgenrefdphi;
TH1D * hgenrefdpt;

TH1D * hj3pt;
TH1D * hj13dphi;
TH1D * hj23dphi;
TH1D * hj13deta;
TH1D * hj23deta;
TH1D * hj123dpt;

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
TH2D * JetTrackBackground1(int condor_iter, int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int dotrkcorr , int mccommand , double jetamin , double jetamax , double vzrange = 15.0, double dijetdphicut = pi56, string whichjet = "", double dvz = 0.2);

int GetNTotTrig();
double getRmin(HiForest * h, double trkEta, double trkPhi);
double getdR2(double jet_eta, double jet_phi, double track_eta, double track_phi);

//calculates dr^2 to avoid the slow TMath() Sqrt function
double getdR2(double jet_eta, double jet_phi, double track_eta, double track_phi)
{
  return TMath::Power(jet_eta-track_eta,2)+TMath::Power(acos(cos(jet_phi-track_phi)),2);
}

double getRmin(HiForest * h, double trkEta, double trkPhi)
{
  double r_min = 9;
  for(int j2 = 0; j2<h->ak3PF.nref; j2++)
  {
    if(TMath::Abs(h->ak3PF.jteta[j2])>2 || TMath::Abs(h->ak3PF.jtpt[j2]) < 50) continue;
    double r_min_temp = TMath::Power(getdR2(h->ak3PF.jteta[j2],h->ak3PF.jtphi[j2],trkEta,trkPhi),0.5);
    if(r_min_temp < r_min) r_min = r_min_temp;
  }
  return r_min;
}

bool jet80()
{
  if(c->hlt.HLT_Jet80_v1 || c->hlt.HLT_Jet80_v2 || c->hlt.HLT_Jet80_v3 || c->hlt.HLT_Jet80_v4 || c->hlt.HLT_Jet80_v5 || c->hlt.HLT_Jet80_v6 || c->hlt.HLT_PAJet80_NoJetID_v1 ) return true;
  return false;
}
          
bool rejectevent(int centmin, int centmax, double vzrange, double dijetdphicut, double leadingjetpthigh, double leadingjetptlow, double subleadingjetpthigh, double subleadingjetptlow, double jetamax, double jetamin, double ajmin, double ajmax, bool writecuts = false) //all cuts go here
{
  bool debug = false;
  if(c->hasGenParticleTree)
  {
    if(!c->skim.pcollisionEventSelection) return true;
  }
  else
  {
    if( c->hlt.Run >= 210498 && c->hlt.Run < 211313 )
    {
      if( c->hlt.Run >= 210676 ) // austin run cut
      {
        if(writecuts) hcuts->Fill(1);
      }
      else
      {
        return true;
      }
    }
    else // no run cut outside pPb range
    {
      if(writecuts) hcuts->Fill(1);
    }  
    if(debug ) 
    {
      cout<<"skim "<<c->hlt.Run<<" "<<c->hlt.Event<<": "<<c->skim.pPAcollisionEventSelectionPA << c->skim.pHBHENoiseFilter<<endl;
    }
    if(writecuts) hcuts->Fill(2);
    if(!c->skim.pPAcollisionEventSelectionPA ) return true;
    if(writecuts) hcuts->Fill(3);
    if(!c->skim.pHBHENoiseFilter ) return true;
    if(writecuts) hcuts->Fill(4);
  }
  return false;
}

void jettrkana(const char * infname = "/mnt/hadoop/cms/store/user/velicanu/mergedv1_sortedforest/mergesortv1_0.root", const char * bkfname = "/mnt/hadoop/cms/store/user/velicanu/mergedv1_sortedforest/mergesortv1_0.root", int trackquality = 0, string whichjet = "")
{
  bool debug = false;
  njet = 0;
  mytrackquality = trackquality;
  cout<<"initializing hiforest and building centrality index"<<endl;
  cout<<"running on: "<<infname<<endl;
  c = new HiForest(infname,"forest1",cPP,0,whichjet);
  bk = c;
  Long64_t nentries = c->GetEntries();

  
  InitCorrFiles(sType);   // Yen-Jie: assuming that we are working on Heavy Ion Data!!!!!!!!
  cout <<"Loading PP correction tables !!"<<endl;
  InitCorrHists(sType);

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

  cout<<"consturcting cent, vz, pcoll, hbhe arrays... ";
  for(int i = 0 ; i < nentries ; ++i)
  {
  // break;
    c->GetEntry(i);
    if(i%1000==0) cout<<i<<"/"<<nentries<<endl;
    if(debug) cout<<centmin<<" "<<centmax<<" "<<vzrange<<" "<<dijetdphicut<<" "<<leadingjetpthigh<<" "<<leadingjetptlow<<" "<<subleadingjetpthigh<<" "<<subleadingjetptlow<<" "<<jetamax<<" "<<jetamin<<" "<<ajmin<<" "<<ajmax<<endl;
    if(rejectevent(centmin,centmax,vzrange,dijetdphicut,leadingjetpthigh,leadingjetptlow,subleadingjetpthigh,subleadingjetptlow,jetamax,jetamin,ajmin,ajmax)) continue;
    if(debug) cout<<"some event pass"<<endl;
    vjetindices.push_back(i);
    njet++;
  }
  cout<<"njet = "<<njet<<endl;
  cout<<"done."<<endl;




}

TH2D * JetTrackSignal(int condor_iter, int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int dotrkcorr , int mccommand , double jetamin , double jetamax , double vzrange, double dijetdphicut, string whichjet)
{
  cout<<"Signal "<<jetindex<<": "<<endl;
  c->ResetBooleans();
  cout<<"whichjet: "<<whichjet<<endl;
  Long64_t nentries = c->GetEntries();
  if(mccommand==2)  dotrkcorr=0;
  TH2D::SetDefaultSumw2(true);
  
  hcuts    = new TH1D(Form("hcuts%d_eta_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),"",10,0.5,10.5);
  hcuts->GetXaxis()->SetBinLabel(1,"total");
  hcuts->GetXaxis()->SetBinLabel(2,"run");
  hcuts->GetXaxis()->SetBinLabel(3,"pcoll");
  hcuts->GetXaxis()->SetBinLabel(4,"hbhe");
  hcuts->GetXaxis()->SetBinLabel(5,"jet80");
  
  float yAxis[23] = {0.5, 0.63, 0.77,  1.03,1.38, 1.84, 2.46, 3.29,  4.40, 5.88,  7.87,  10.52, 14.06,  18.8, 25.13,  33.58,  44.89,  60, 80, 100, 120, 140, 200};
  float xAxis[121];
  for (int i = 0; i < 121; ++i)
  {
    xAxis[i] = 2.5*i;
  }
  
  hjettrk = new TH2D(Form("hjettrk_%s_%djet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",120,xAxis,22,yAxis);
  
  hnjets    = new TH1D(Form("hnjets%d_eta_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),"",11,-0.5,10.5);
  hgenrefeta    = new TH2D(Form("hgenref%d_eta_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",80,-4,4,80,-4,4);
  hgenrefphi    = new TH2D(Form("hgenref%d_phi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#phi",80,-TMath::Pi(),TMath::Pi(),80,-TMath::Pi(),TMath::Pi());
  hgenrefpt    = new TH2D(Form("hgenref%d_pt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";p_{T}",80,0,leadingjetpthigh,80,0,leadingjetpthigh);
  hgenrefdeta    = new TH1D(Form("hgenref%d_deta_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#Delta#eta",400,-1,1);
  hgenrefdphi    = new TH1D(Form("hgenref%d_dphi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#Delta#phi",400,-0.3,0.3);
  hgenrefdpt    = new TH1D(Form("hgenref%d_dpt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#Delta p_{T}",800,-20,20);
  hjeteta    = new TH1D(Form("hjet%d_eta_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",80,-4,4);
  hjetphi    = new TH1D(Form("hjet%d_phi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#phi",80,-TMath::Pi(),TMath::Pi());
  hjetpt    = new TH1D(Form("hjet%d_pt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";p_{T}",80,0,leadingjetpthigh);
  halljetpt    = new TH1D(Form("halljet%d_pt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";p_{T}",80,0,leadingjetpthigh);
  
  hcent    = new TH1D(Form("hcent%d_pt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";p_{T}",80,0,leadingjetpthigh);
  hcentpostcut    = new TH1D(Form("hcentpostcut%d_pt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";p_{T}",80,0,leadingjetpthigh);
  
  
  hjetdphi    = new TH1D(Form("hjet%d_dphi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#Delta#phi",40,0,TMath::Pi());

  hjtrketa    = new TH1D(Form("hj%d_trketa_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",100,-2.5,2.5);
  hjtrketanocorr    = new TH1D(Form("hj%d_trketanocorr_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",100,-2.5,2.5);
  hjtrkphi    = new TH1D(Form("hj%d_trkphi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#phi",80,-TMath::Pi(),TMath::Pi());
  hjtrkpt    = new TH1D(Form("hj%d_trkpt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";p_{T}",80,0,ptasshigh*2);
  halltrkpt    = new TH1D(Form("hj%d_alltrkpt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";p_{T}",120,0.5,120.5);

  hjtrkdeta    = new TH1D(Form("hj%d_trkdeta_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#eta",100,-5,5);
  hjtrkdphi    = new TH1D(Form("hj%d_trkdphi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#phi",80,0,TMath::Pi());

  if(jetindex==2)
  {
    hj3pt = new TH1D(Form("hj%d_hj3pt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";3^{rd} jet p_{T}",100,0,200);
    hj13dphi = new TH1D(Form("hj%d_hj13dphi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#Delta#Phi Jet_{1-3}",50,-pi,pi);
    hj23dphi = new TH1D(Form("hj%d_hj23dphi_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#Delta#Phi Jet_{2-3}",50,-pi,pi);
    hj13deta = new TH1D(Form("hj%d_hj13deta_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#Delta#Eta Jet_{1-3}",50,-5,5);
    hj23deta = new TH1D(Form("hj%d_hj23deta_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#Delta#Eta Jet_{2-3}",50,-5,5);
    hj123dpt = new TH1D(Form("hj%d_hj123dpt_trg%d_%d_ass%d_%d_nmin%d_nmax%d",jetindex,(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax),";#Delta p_{T} Jet_{1-2-3}",50,-200,200);
  }
  
  
  
  // hsigvz = new TH1D("hsigvz",";vz;",100,-16,16);
  // hsigcent = new TH1D("hsigcent",";vz;",200,0,200);

  // Long64_t nentries = c->GetEntries();
  // Long64_t nbytes = 0, nb = 0;
  TH2D * hJetTrackSignal;
  // cout<<"before first exit"<<endl;
  if(jetindex==0) //leading jet
  {
    //Hallie Binning
    hJetTrackSignal = new TH2D(Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",40,-5,5,40,-0.5*TMath::Pi(),1.5*TMath::Pi());
    
    //Pelin Binning
    // hJetTrackSignal = new TH2D(Form("signal_%s_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,75,-1.0*TMath::Pi(),2.0*TMath::Pi());

  
  }
  else if(jetindex==1) //subleading jet
  {
    //Hallie Binning
    hJetTrackSignal = new TH2D(Form("signal_%s_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",40,-5,5,40,-0.5*TMath::Pi(),1.5*TMath::Pi());
    JetTrackSignalShifted = new TH2D(Form("signal_%s_shifted_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",40,-5,5,40,-0.5*TMath::Pi(),1.5*TMath::Pi());
    
    //Pelin Binning
    // hJetTrackSignal = new TH2D(Form("signal_%s_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,75,-1.0*TMath::Pi(),2.0*TMath::Pi());
  
  }
  else if(jetindex==2) //3rd jet
  {
    //Hallie Binning
    hJetTrackSignal = new TH2D(Form("signal_%s_3rdjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",40,-5,5,40,-0.5*TMath::Pi(),1.5*TMath::Pi());
    
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
  bool writecuts = true; 
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    c->GetEntry(jentry);
    if(jentry%1000==0) cout<<jentry<<"/"<<nentries<<endl;

    float cmetashift = 0.0;
    if( c->hlt.Run >= 210498 && c->hlt.Run < 211313 ) cmetashift = 0.4654;  //pPb range
    if( c->hlt.Run >= 211313 && c->hlt.Run <= 211631 ) cmetashift = -0.4654; //Pbp range

    
    if(rejectevent(centmin,centmax,vzrange,dijetdphicut,leadingjetpthigh,leadingjetptlow,subleadingjetpthigh,subleadingjetptlow,jetamax,jetamin,ajmin,ajmax,writecuts)) continue;
    bool debug = false;
    if(debug && !jet80()) cout<<c->hlt.Run<<" "<<c->hlt.Event<<" "<<c->hlt.HLT_PAJet80_NoJetID_v1<<endl;
    if(!jet80() && mccommand == 0) continue;
    c->hasTrackTree = true;
    c->GetEntry(jentry);
    if(debug) cout<<"ngen: "<<c->myjet.ngen<<endl;
    for (int i = 0; i < c->myjet.ngen; ++i)
    {
      if(fabs(c->myjet.geneta[i] + cmetashift) > jetamax) continue;
      if(c->myjet.genpt[i] < 15) continue;
      for (int j = 0; j < c->track.nTrk; ++j)
      {
        if(fabs(c->track.trkEta[j])<2.4&&c->track.highPurity[j]&&fabs(c->track.trkDz1[j]/c->track.trkDzError1[j])<3&&fabs(c->track.trkDxy1[j]/c->track.trkDxyError1[j])<3&&c->track.trkPtError[j]/c->track.trkPt[j]<0.1)
        {
          if(c->track.trkPt[j]<=0.5) continue;
          // cout<<c->track.trkPt[j]<<endl;
          double trkRMin = getRmin(c, c->track.trkEta[j], c->track.trkPhi[j]);
          float effweight = factorizedPtCorr(getPtBin(c->track.trkPt[j], sType), 1, c->track.trkPt[j], c->track.trkPhi[j], c->track.trkEta[j], trkRMin, sType);
          if(dotrkcorr==0) effweight = 1.0; //turns off tracking correction
          double deta = -99, dphi = -99;
          deta = fabs(c->myjet.geneta[i] - c->track.trkEta[j]);
          dphi = fabs(c->myjet.genphi[i] - c->track.trkPhi[j]);
          if( dphi > pi ) dphi = 2*pi - dphi;
          double dr2 = (deta*deta)+(dphi*dphi);
          double dr = sqrt(dr2);
          if(dr>0.3)
          {
            hjettrk->Fill(c->myjet.genpt[i],c->track.trkPt[j],effweight);
          }
        }
      }
    }



    std::vector<int> vjets ; 
    for(int i = 0 ; i < c->myjet.ngen ; ++i)
    {
      if(fabs(c->myjet.geneta[i] + cmetashift) > jetamax) continue;
      if(fabs(c->myjet.genpt[i]) < 25) continue;
      vjets.push_back(i);
    }
    hnjets->Fill(vjets.size());
    if(vjets.size() <= 1) continue;
    float dijetdphi = fabs(c->myjet.genphi[vjets[0]] - c->myjet.genphi[vjets[1]]);
    if(dijetdphi > pi) dijetdphi = 2 * pi - dijetdphi;
    if(dijetdphi < dijetdphicut) continue;
    //! HARDCODED, PLZ FIX
    if(c->myjet.genpt[vjets[0]] < 120 || c->myjet.genpt[vjets[0]] > 300) continue;
    if(c->myjet.genpt[vjets[1]] < 50 || c->myjet.genpt[vjets[1]] > 300) continue;
    if((int)vjets.size() <= jetindex) continue;
    if(jetindex==2)
    {
      if(c->myjet.genpt[vjets[jetindex]] < 15) continue;
    }
    // if(vjets.size() != JETSIZE) continue;
    if(vjets.size() < 3) continue;
    // cout<<vjets.size()<<endl;
    
    // for(int i = 0 ; i < (int)vjets.size() ; ++i)
    // {
      hgenrefeta->Fill(c->myjet.geneta[vjets[jetindex]],c->myjet.refeta[vjets[jetindex]]);
      hgenrefphi->Fill(c->myjet.genphi[vjets[jetindex]],c->myjet.refphi[vjets[jetindex]]);
      hgenrefpt->Fill(c->myjet.genpt[vjets[jetindex]],c->myjet.refpt[vjets[jetindex]]);
      hgenrefdeta->Fill(c->myjet.geneta[vjets[jetindex]]-c->myjet.refeta[vjets[jetindex]]);
      hgenrefdphi->Fill(c->myjet.genphi[vjets[jetindex]]-c->myjet.refphi[vjets[jetindex]]);
      hgenrefdpt->Fill(c->myjet.genpt[vjets[jetindex]]-c->myjet.refpt[vjets[jetindex]]);
      double jeteta = c->myjet.geneta[vjets[jetindex]];
      double jetphi = c->myjet.genphi[vjets[jetindex]];
      if(jetindex==2)
      {
        //! Fill third jet histograms
        hj3pt->Fill(c->myjet.jtpt[vjets[2]]);
        float j13dphi, j23dphi, j13deta, j23deta, j123dpt;
        j13dphi = c->myjet.jtphi[vjets[0]] - c->myjet.jtphi[vjets[2]];
        j23dphi = c->myjet.jtphi[vjets[1]] - c->myjet.jtphi[vjets[2]];
        j13deta = c->myjet.jteta[vjets[0]] - c->myjet.jteta[vjets[2]];
        j23deta = c->myjet.jteta[vjets[1]] - c->myjet.jteta[vjets[2]];
        j123dpt = c->myjet.jtpt[vjets[0]] - c->myjet.jtpt[vjets[1]] - c->myjet.jtpt[vjets[2]];
        hj13dphi->Fill(j13dphi);
        hj23dphi->Fill(j23dphi);
        hj13deta->Fill(j13deta);
        hj23deta->Fill(j23deta);
        hj123dpt->Fill(j123dpt);
      }
      
      hjeteta->Fill(c->myjet.geneta[vjets[jetindex]]);
      hjetphi->Fill(c->myjet.genphi[vjets[jetindex]]);
      // double dijetdphi = fabs(c->myjet.genphi[leadindex] - c->myjet.genphi[subleadindex]);
      // if( dijetdphi > pi ) dijetdphi = 2*pi - dijetdphi;
      // hjetdphi->Fill(dijetdphi);
      hjetpt->Fill(c->myjet.genpt[vjets[jetindex]]);
      ntottrig += 1;
      c->hasTrackTree = true;
      c->GetEntry(jentry);
      // InitPosArrPbPb(c->evt.hiBin);
      // hsigvz->Fill(c->evt.vz);
      // hsigcent->Fill(c->evt.hiBin);
      for(int j = 0 ; j < c->track.nTrk ; ++j)
      {
        if(fabs(c->track.trkEta[j])<2.4&&c->track.highPurity[j]&&fabs(c->track.trkDz1[j]/c->track.trkDzError1[j])<3&&fabs(c->track.trkDxy1[j]/c->track.trkDxyError1[j])<3&&c->track.trkPtError[j]/c->track.trkPt[j]<0.1)
        {
          halltrkpt->Fill(c->track.trkPt[j]);
          if(c->track.trkPt[j]>ptasshigh || c->track.trkPt[j]<ptasslow) continue;

          double trkRMin = getRmin(c, c->track.trkEta[j], c->track.trkPhi[j]);
          float effweight = factorizedPtCorr(getPtBin(c->track.trkPt[j], sType), 1, c->track.trkPt[j], c->track.trkPhi[j], c->track.trkEta[j], trkRMin, sType);
          if(dotrkcorr==0) effweight = 1.0; //turns off tracking correction

          double deta = -99, dphi = -99, dphi_shifted = -99;
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
            if(jetindex==1)
            {
              dphi_shifted = dphi + pi;
              if( dphi_shifted > pi ) dphi_shifted = 2*pi - dphi_shifted;
              JetTrackSignalShifted->Fill(deta,dphi_shifted,effweight);
            }
          }
          hjtrketa->Fill(c->track.trkEta[j],effweight);
          hjtrketanocorr->Fill(c->track.trkEta[j]);
          hjtrkphi->Fill(c->track.trkPhi[j],effweight);
          hjtrkpt->Fill(c->track.trkPt[j],effweight);
          hjtrkdeta->Fill(deta);
          hjtrkdphi->Fill(dphi);

        }
      }
    // }

    // exit(1);
    c->hasTrackTree = false;
    if(hasEVT)
    {
      hcentpostcut->Fill(c->evt.hiBin);
    }
    
  }
  cout<<ntottrig<<endl;
  // myfile.close();
  return hJetTrackSignal;
}

// TH2D * JetTrackBackground1(int condor_iter, int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int dotrkcorr , int mccommand , double jetamin , double jetamax , double vzrange, double dijetdphicut, string whichjet, double dvz)
TH2D * JetTrackBackground(int condor_iter, int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int dotrkcorr , int mccommand , double jetamin , double jetamax , double vzrange, double dijetdphicut, string whichjet, double dvz)
{
  cout<<"Background"<<jetindex<<": "<<endl;
  c->ResetBooleans();
  cout<<"whichjet: "<<whichjet<<endl;
  Long64_t nentries = c->GetEntries();
  if(mccommand==2)  dotrkcorr=0;
  TH2D::SetDefaultSumw2(true);
  
  

  // Long64_t nentries = c->GetEntries();
  // Long64_t nbytes = 0, nb = 0;
  TH2D * hJetTrackBackground;
  // cout<<"before first exit"<<endl;
  if(jetindex==0) //leading jet
  {
    //Hallie Binning
    hJetTrackBackground = new TH2D(Form("background_%s_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",40,-5,5,40,-0.5*TMath::Pi(),1.5*TMath::Pi());
  
  }
  else if(jetindex==1) //subleading jet
  {
    //Hallie Binning
    JetTrackBackgroundShifted = new TH2D(Form("background_%s_shifted_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",40,-5,5,40,-0.5*TMath::Pi(),1.5*TMath::Pi());
    hJetTrackBackground = new TH2D(Form("background_%s_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",40,-5,5,40,-0.5*TMath::Pi(),1.5*TMath::Pi());
  }
  else if(jetindex==2) //3rd jet
  {
    //Hallie Binning
    hJetTrackBackground = new TH2D(Form("background_%s_3rdjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",40,-5,5,40,-0.5*TMath::Pi(),1.5*TMath::Pi());
  
  }
  else //leading jet
  {
    cout<<"Error: only jetindex 0 and 1 currently supported (leading + subleading jet)"<<endl;
    exit(1);
  }
  // return hJetTrackBackground;
  
  double ntottrig = 0;
  bool writecuts = true; 
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    c->GetEntry(jentry);
    if(jentry%1000==0) cout<<jentry<<"/"<<nentries<<endl;

    float cmetashift = 0.0;
    if( c->hlt.Run >= 210498 && c->hlt.Run < 211313 ) cmetashift = 0.4654;  //pPb range
    if( c->hlt.Run >= 211313 && c->hlt.Run <= 211631 ) cmetashift = -0.4654; //Pbp range

    
    if(rejectevent(centmin,centmax,vzrange,dijetdphicut,leadingjetpthigh,leadingjetptlow,subleadingjetpthigh,subleadingjetptlow,jetamax,jetamin,ajmin,ajmax,writecuts)) continue;
    bool debug = false;
    if(debug && !jet80()) cout<<c->hlt.Run<<" "<<c->hlt.Event<<" "<<c->hlt.HLT_PAJet80_NoJetID_v1<<endl;
    if(!jet80() && mccommand == 0) continue;
    c->hasTrackTree = true;
    c->GetEntry(jentry);
    // if(debug) cout<<"ngen: "<<c->myjet.ngen<<endl;


    std::vector<int> vjets ; 
    //! note this isn't exactly dijet selection, it's top 3 jets in |eta| < 1.5
    for(int i = 0 ; i < c->myjet.ngen ; ++i)
    {
      // if(c->myjet.genpt[i] > leadingjetpthigh) continue;
      // if(c->myjet.genpt[i] < leadingjetptlow) continue;
      if(fabs(c->myjet.geneta[i] + cmetashift) > jetamax) continue;
      if(fabs(c->myjet.genpt[i]) < 25) continue;
      vjets.push_back(i);
    }
    // hnjets->Fill(vjets.size());
    // for(int i = 0 ; i < (int)vjets.size() ; ++i)
    // {
    if(vjets.size() <= 1) continue;
    // cout<<vjets.size()<<endl;
    float dijetdphi = fabs(c->myjet.genphi[vjets[0]] - c->myjet.genphi[vjets[1]]);
    if(dijetdphi > pi) dijetdphi = 2 * pi - dijetdphi;
    if(dijetdphi < dijetdphicut) continue;
    //! HARDCODED, PLZ FIX
    if(c->myjet.genpt[vjets[0]] < 120 || c->myjet.genpt[vjets[0]] > 300) continue;
    if(c->myjet.genpt[vjets[1]] < 50 || c->myjet.genpt[vjets[1]] > 300) continue;
    if((int)vjets.size() <= jetindex) continue;
    // if(vjets.size() != JETSIZE) continue;
    if(vjets.size() < 3) continue;
    // cout<<vjets.size()<<endl;
    double jeteta = c->myjet.geneta[vjets[jetindex]];
    double jetphi = c->myjet.genphi[vjets[jetindex]];
    
    ntottrig += 1;
    // c->hasTrackTree = true;
    // c->GetEntry(jentry);
    for(int imix = jentry+1; imix < jentry + 11 ; ++imix)
    {
      c->GetEntry(imix);
      for(int j = 0 ; j < c->track.nTrk ; ++j)
      {
        if(fabs(c->track.trkEta[j])<2.4&&c->track.highPurity[j]&&fabs(c->track.trkDz1[j]/c->track.trkDzError1[j])<3&&fabs(c->track.trkDxy1[j]/c->track.trkDxyError1[j])<3&&c->track.trkPtError[j]/c->track.trkPt[j]<0.1)
        {
          halltrkpt->Fill(c->track.trkPt[j]);
          if(c->track.trkPt[j]>ptasshigh || c->track.trkPt[j]<ptasslow) continue;

          double trkRMin = getRmin(c, c->track.trkEta[j], c->track.trkPhi[j]);
          float effweight = factorizedPtCorr(getPtBin(c->track.trkPt[j], sType), 1, c->track.trkPt[j], c->track.trkPhi[j], c->track.trkEta[j], trkRMin, sType);
          if(dotrkcorr==0) effweight = 1.0; //turns off tracking correction

          double deta = -99, dphi = -99, dphi_shifted = -99;
          // cout<<c->track.trkPt[j]<<","<<c->track.trkEta[j]<<","<<c->track.trkPhi[j]<<","<<trkRMin<<","<<effweight<<endl;
          if(false) //! Dragos:
          {
            deta = fabs(jeteta-c->track.trkEta[j]);
            dphi = fabs(jetphi-c->track.trkPhi[j]);
            if( dphi > pi ) dphi = 2*pi - dphi;
            // double effweight = c->track.trkPt[j]*correction;
            hJetTrackBackground->Fill(deta,dphi,effweight);
            hJetTrackBackground->Fill(-deta,dphi,effweight);
            hJetTrackBackground->Fill(deta,-dphi,effweight);
            hJetTrackBackground->Fill(-deta,-dphi,effweight);
            hJetTrackBackground->Fill(deta,(2*pi)-dphi,effweight);
            hJetTrackBackground->Fill(-deta,(2*pi)-dphi,effweight);
          }
          else //! Hallie
          {
            deta = jeteta - c->track.trkEta[j];
            dphi = jetphi - c->track.trkPhi[j];
            if(dphi>(1.5*TMath::Pi()))  { dphi+= -2*TMath::Pi(); }
            if(dphi<(-0.5*TMath::Pi())) { dphi+=  2*TMath::Pi(); }
            
            // hJetTrackBackground->Fill(deta,dphi, effweight);
            hJetTrackBackground->Fill(deta,dphi, effweight);
            if(jetindex==1)
            {
              dphi_shifted = dphi + pi;
              if(dphi_shifted>(1.5*TMath::Pi()))  { dphi_shifted+= -2*TMath::Pi(); }
              if(dphi_shifted<(-0.5*TMath::Pi())) { dphi_shifted+=  2*TMath::Pi(); }
              JetTrackBackgroundShifted->Fill(deta,dphi_shifted,effweight);
            }
          }
          hjtrketa->Fill(c->track.trkEta[j],effweight);
          hjtrketanocorr->Fill(c->track.trkEta[j]);
          hjtrkphi->Fill(c->track.trkPhi[j],effweight);
          hjtrkpt->Fill(c->track.trkPt[j],effweight);
          hjtrkdeta->Fill(deta);
          hjtrkdphi->Fill(dphi);

        }
      }
    }

    c->hasTrackTree = false;
    if(hasEVT)
    {
      hcentpostcut->Fill(c->evt.hiBin);
    }
    
  }
  cout<<ntottrig<<endl;
  return hJetTrackBackground;
}


TH2D * JetTrackBackground1(int condor_iter, int jetindex, double leadingjetptlow , double leadingjetpthigh , double subleadingjetptlow , double subleadingjetpthigh , double ptasslow , double ptasshigh, int centmin, int centmax, float ajmin , float ajmax , int dotrkcorr , int mccommand , double jetamin , double jetamax , double vzrange, double dijetdphicut, string whichjet, double dvz)
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
  // TString csvname;
  // if(condor_iter<10)
    // csvname = Form("mixevents_0%d.csv",condor_iter);
  // else
    // csvname = Form("mixevents_%d.csv",condor_iter);
  // std::vector< std::vector< long long int > > mixevents = csvtoarraylong(csvname.Data()); 
  // cout<<(int)mixevents.size()<<" x "<<(int)mixevents[0].size()<<endl;
  

  // int parallelization = 40;
  // int parallelization = 1;

  // c->LoadNoTrees();
  // bk->LoadNoTrees();
  
  if(mccommand==2)  dotrkcorr=0;
  TH2D::SetDefaultSumw2(true);
  // bk->hasTrackTree = true;
  // c->hasEvtTree = true;
  // c->hasSkimTree = true;
  // c->hasAk3CaloJetTree = true;
  // bk->hasEvtTree = true;
  // bk->hasSkimTree = true;
  // bk->hasAk3CaloJetTree = true;
  // bk->hasTrackTree = true;

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
    hJetTrackBackground = new TH2D(Form("background_%s_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",40,-5,5,40,-0.5*TMath::Pi(),1.5*TMath::Pi());

    //Pelin Binning
    // hJetTrackBackground = new TH2D(Form("background_%s_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,75,-1.0*TMath::Pi(),2.0*TMath::Pi());
  }
  else if(jetindex==1) //subleading jet
  {
    //Hallie Binning
    hJetTrackBackground = new TH2D(Form("background_%s_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",40,-5,5,40,-0.5*TMath::Pi(),1.5*TMath::Pi());

    //Pelin Binning
    // hJetTrackBackground = new TH2D(Form("background_%s_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptasslow,(int)ptasshigh,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)),";#Delta#eta;#Delta#phi",100,-5,5,75,-1.0*TMath::Pi(),2.0*TMath::Pi());


  }
  else //leading jet
  {
    cout<<"Error: only jetindex 0 and 1 currently supported (leading + subleading jet)"<<endl;
    exit(1);
  }

  // return hJetTrackBackground;
  
  double ntottrig = 0;
  // int n_entries_in_cent_range = cent_index_start[centmax] - cent_index_start[centmin];


  // cout<<"before event for loop "<<n_entries_in_cent_range<<endl;
  int sigindex = 0;
  float ntries = 0.0;
  float njetsevts = 0.0;
  // nentries = (int)vjetindices.size();
  // int increment = nentries/parallelization;
  // int stopentry = (condor_iter+1)*increment;
  // if(stopentry > nentries)
  // {
    // cout<<"FATAL ERROR: trying to run on more events than we have"<<endl; exit(1);
  // }
  // if((condor_iter+2)*increment > nentries) stopentry = nentries;

  // cout<<increment<<endl;
  // cout<<"start: "<<condor_iter*increment<<" to end "<<stopentry<<" of "<<nentries<<endl;

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
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
    if(jentry%1000==0) 
		cout<<jentry<<"/"<<nentries<<endl;
    c->GetEntry(jentry);
    // cout<<mixevents[nextindex][0]<<" "<<vjetindices[jentry]<<endl;
    // std::unordered_set<long long int> tempset;
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
    // cout<<"here back"<<endl;
    // if(jentry>100) break;
    sigindex++;

    if( c->myjet.ngen < 2 ) continue;
    int leadindex = -1 , subleadindex = -1;
    bool foundleading = false;
    for(int i = 0 ; i < c->myjet.ngen ; ++i)
    {
      if(!foundleading)
      {
        if(fabs(c->myjet.geneta[i])<2.0 )
        {
          leadindex = i;
          foundleading = true;
        }
      }
      else
      {
        if(fabs(c->myjet.geneta[i])<2.0 ) 
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
    
    
    double jeteta = c->myjet.geneta[dojet];
    double jetphi = c->myjet.genphi[dojet];
    int k = 0;
    int statfactor = 20;
    njetsevts+=1.0;
    if(hasEVT)
    {
      hsigvz2->Fill(c->evt.vz);
    }
    // for(int backindex = 0 ; backindex < nbkentries ; ++backindex)
    for(int backindex = 0 ; backindex < nbkentries ; ++backindex)
    {
      if(k>=statfactor) break;
      currentbackindex = ( currentbackindex + 1 ) % nbkentries;
			ntries+=1.0;
      // if(abs((sigcentarray[vjetindices[jentry]]-bakcentarray[currentbackindex]))>5) continue;
      // if(abs((sigcentarray[vjetindices[jentry]]-bakcentarray[currentbackindex]))>0) continue;
      // if(fabs(sigvzarray[vjetindices[jentry]] - bakvzarray[currentbackindex])>dvz) continue;
      // cout<<"searching for: "<<backgroundids[backindex]<<endl;
      
			//! fancy part turned off
      // auto search = tempset.find(backgroundids[backindex]);
			// if(search == tempset.end()) continue;
      bk->GetEntry(currentbackindex);
      // cout<<"really here"<<endl;

      // if(c->hlt.Run == bk->evt.run && c->evt.evt == bk->evt.evt) continue; //same event
      if((bk->genParticleTree != 0))
      {
        //if is PbPb MC
        if(!bakpcollarray[currentbackindex]) continue;
      }
      else
      {
        //if is PbPb Data
        // if(!(bakpcollarray[currentbackindex] && bakhbhearray[currentbackindex] )) continue;
      }
      k++;
      ntottrig += 1;
      // InitPosArrPbPb(bk->evt.hiBin);
      if(hasEVT)
      {
        hmixvz->Fill(bk->evt.vz);
        hmixdvz->Fill(c->evt.vz - bk->evt.vz);
        hmixcent->Fill(bk->evt.hiBin);
        hmixdcent->Fill(c->evt.hiBin - bk->evt.hiBin); 
      }
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
          // float trkRMin = getTrkRMin(bk->track.trkPhi[j], bk->track.trkEta[j],bk->myjet.ngen,bk->myjet.genphi,bk->myjet.geneta,bk->myjet.genpt);
          // float effweight = factorizedPtCorr(getPtBin(bk->track.trkPt[j], sType), bk->evt.hiBin, bk->track.trkPt[j], bk->track.trkPhi[j], bk->track.trkEta[j], trkRMin, sType);
          // if(dotrkcorr==0) effweight = 1.0; //turns off tracking correction
          float effweight = 1.0; //turns off tracking correction
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
  // bk->ResetBooleans();
  cout<<"-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"<<endl;
  return hJetTrackBackground;
}


int GetNTotTrig()
{
  int tmp = ntottrig;
  ntottrig = 0;
  return tmp;
}
