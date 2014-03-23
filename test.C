{
TFile a("output-801.root");
t->Draw("trkPt>>h(100,0,10)","trkWt");
t->Draw("pPt>>h2(100,0,10)","");
h->Divide(h2)
h->Draw();
}

TFile a("full_ntuple_HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0_pt0_10_akVs3Calo_doGenJet1.root")
TFile b("../../mpt/output-801.root")
TH1D *h2 = (TH1D*)a->Get("hPtUncor")
t->Draw("trkPt>>h(100,0,10)","")
h2->Divide(h)
h2->Draw()

TFile a("full_ntuple_HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0_pt0_10_akVs3Calo_doGenJet1.root")
TFile b("../../mpt/output-801.root")
TH1D *h1 = (TH1D*)a->Get("hRmin");
TH1D *h2 = (TH1D*)b->Get("hRmin");
h1->Divide(h2);
h1->Draw();
