#include "HiForestAnalysis/hiForest.h"
#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <iostream>

// Example of forest skim

void skimTree(char *infname = "../JetSample/hiForest_Jet80or95_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21_0.root", char *outfname = "skim_Jet.root")
{
  // Define the input file and HiForest
  HiForest *c = new HiForest(infname);
  c->hasHltTree=0;
  c->hasPFTree=0;
  c->hasPhotonTree=0;
  c->hasTowerTree=0;
  c->hasHbheTree=0;
  c->hasEbTree=0;
  c->hasGenpTree=0;
  c->hasGenParticleTree=0;   
  c->hasAk5CaloJetTree=0;
  c->hasAkPu2CaloJetTree=0;
  c->hasAkPu3CaloJetTree=0;
  c->hasAkPu4CaloJetTree=0;
  c->hasAkPu5CaloJetTree=0;
  c->hasAkPu2JetTree=0;
  c->hasAkPu3JetTree=0;
  c->hasAkPu4JetTree=0;
  c->hasAkPu5JetTree=0;
  c->hasAkVs2PFJetTree=0;
  c->hasAkVs3PFJetTree=0;
  c->hasAkVs4PFJetTree=0;
  c->hasAkVs5PFJetTree=0;
  c->SetOutputFile(outfname);

  int filtered=0;
  
  // ### START Define CUTS
  float dijetdphicut = 5.0 * TMath::Pi() / 6.0 ; 
  float vzrange = 15;
  float leadingjetpthigh = 300;
  float leadingjetptlow = 120;
  float subleadingjetpthigh = 300;
  float subleadingjetptlow = 50;
  float jetamax = 1.6;
  float jetamin = 0;
  float ajmax = 1.0;
  float ajmin = 0;
float pi = TMath::Pi();

  // ### END Define Cuts

  // Main loop
  for (int i=0;i<c->GetEntries();i++)
  {
    c->GetEntry(i);
    if (i%1000==0) cout <<filtered<<" "<<i<<" / "<<c->GetEntries()<<endl;
    
    int flag=0;    // # of jets with |eta|<2 and pt > 120
    int flag2=0;   // # of jets with |eta|<2 and pt > 50
    for (int j=0;j<c->akVs3Calo.nref;j++) {
      if (fabs(c->akVs3Calo.jteta[j])>2) continue;
      if (c->akVs3Calo.jtpt[j]>120) flag=1;
      if (c->akVs3Calo.jtpt[j]>50) flag2++;
      if (flag>=1&&flag2>=2) break;
    }
    
    
    // ### START CUTS ###
    
    if(!(c->skim.pcollisionEventSelection && c->skim.pHBHENoiseFilter )) continue;
    if( fabs(c->evt.vz) > vzrange ) continue;
    
    
    // ### Find leading and subleading jets in |eta| < 2
    int leadindex = -1 , subleadindex = -1;
    bool foundleading = false;
    double oldjetpt = 999999999;
    for(int i = 0 ; i < c->akVs3Calo.nref ; ++i)
    {
      double deltajetpt = oldjetpt - c->akVs3Calo.jtpt[i];
      if(deltajetpt < 0) { cout<<"FATAL ERROR: unsorted jets"<<endl; exit(1); }
      oldjetpt = c->akVs3Calo.jtpt[i];
      if(!foundleading)
      {
        if(fabs(c->akVs3Calo.jteta[i])<2.0 )
        {
          leadindex = i;
          foundleading = true;
        }
      }
      else
      {
        if(fabs(c->akVs3Calo.jteta[i])<2.0 ) 
        {
          subleadindex = i;
          break;
        }
      }
    }
    if(leadindex < 0 || subleadindex < 0 ) continue; //Didn't find jet pair in |eta|<2
    
    // ### START jet cuts ###
    double dijetdphi = fabs(c->akVs3Calo.jtphi[leadindex] - c->akVs3Calo.jtphi[subleadindex]);
    if( dijetdphi > pi ) dijetdphi = 2*pi - dijetdphi;
    if( fabs(dijetdphi) < dijetdphicut ) continue;  // dphi cut
    if( c->akVs3Calo.jtpt[leadindex] > leadingjetpthigh ) continue;       //lead pt cut
    if( c->akVs3Calo.jtpt[leadindex] < leadingjetptlow ) continue;        //lead pt cut
    if( c->akVs3Calo.jtpt[subleadindex] > subleadingjetpthigh ) continue; //sublead pt cut
    if( c->akVs3Calo.jtpt[subleadindex] < subleadingjetptlow ) continue;  //sublead pt cut
    if( fabs(c->akVs3Calo.jteta[leadindex]) > jetamax || fabs(c->akVs3Calo.jteta[leadindex]) < jetamin ) continue;        //jet eta cut
    if( fabs(c->akVs3Calo.jteta[subleadindex]) > jetamax || fabs(c->akVs3Calo.jteta[subleadindex]) < jetamin ) continue;  //jet eta cut
    // if ((c->akVs3Calo.trackMax[leadindex]/c->akVs3Calo.jtpt[leadindex])<0.01) continue;         //trackmax cut
    // if ((c->akVs3Calo.trackMax[subleadindex]/c->akVs3Calo.jtpt[subleadindex])<0.01) continue;   //trackmax cut
    float aj = (c->akVs3Calo.jtpt[leadindex] - c->akVs3Calo.jtpt[subleadindex]) / (c->akVs3Calo.jtpt[leadindex] + c->akVs3Calo.jtpt[subleadindex]);
    if(aj < ajmin || aj > ajmax) continue;  //aj cut
      
// ### END CUTS ###
      
      
    if (flag>=1&&flag2>=2) {
      c->FillOutput(); // Write output forest
      filtered++;
    }  
  }

  delete c;
}

int main(int argc, char *argv[])
{
  if(argc != 3)
  {
    std::cout << "Usage: ./skimTree.exe <input> <output>" << std::endl;
    return 1;
  }
  skimTree(argv[1], argv[2]);

  return 0;
}

