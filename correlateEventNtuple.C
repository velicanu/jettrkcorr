#include "HiForestAnalysis/hiForest.h"
#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <iostream>

// Example of forest skim

void turnOffBranches(HiForest *c)
{
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
}

void correlateEventNtuple(char *infname = "0.root", char *mbname = "mb.root", char *outfname = "matched_Jet.root",int startEntry=0,int endEntry=0)
{
   // Define the input file and HiForest

   HiForest *d = new HiForest(mbname);
   turnOffBranches(d);
   d->SetOutputFile(outfname);
   cout <<"good good"<<endl;

   HiForest *c = new HiForest(infname);
   turnOffBranches(c);
   int filtered=0;
   int idx=0;

   TFile *outf = new TFile(outfname,"recreate");
   int m[20];
   TTree *nt = new TTree("nt","");  
   nt->Branch("m",m,"m[20]/I");
   // Main loop
   if (endEntry==0) endEntry=c->GetEntries();

   for (int i=startEntry;i<=endEntry;i++)
   {
      c->evtTree->GetEntry(i);
      c->skimTree->GetEntry(i);
      if (i%10==0) cout <<filtered<<" "<<i<<" / "<<c->GetEntries()<<endl;
      int nMatched=0;
      while (nMatched<20)
      {
         idx++;
         cout <<i<<" / " << c->GetEntries() << " / "<<idx<<" / "<<nMatched<<"\r";
         if (idx>=d->GetEntries()) idx=0;
         if (idx==i) continue;
         d->evtTree->GetEntry(idx);
         if (fabs(c->evt.hiBin-d->evt.hiBin)>0) continue;
         if ((fabs(c->evt.vz-d->evt.vz)>0.5)) continue;

         d->skimTree->GetEntry(idx);
         if (!(d->skim.pcollisionEventSelection && d->skim.pHBHENoiseFilter )) continue;        
         m[nMatched]=idx;
         nMatched++;
         //d->GetEntry(idx);
         //d->FillOutput(); // Write output forest 
      }   
         cout <<std::endl;
         nt->Fill();
   }

   delete c;
   delete d;
   outf->Write();
   outf->Close();
}
