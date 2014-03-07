#include "jettrkana.C"
#include <fstream>

// void runcorr(int condor_iter, int trackqual, string flist = "", string tag = "", int nmin = 220, int nmax = 1000, float pttrigmin = 1, float pttrigmax = 2, float ptassmin = 1, float ptassmax = 2)
void runcorr(int condor_iter = 0, int trackqual = 0, string flist = "", string tag = "", int nmin = 220, int nmax = 1000, float pttrigmin = 1, float pttrigmax = 2, float ptassmin = 1, float ptassmax = 2, float leadingjetptlow = 100, float leadingjetpthigh = 300 , float subleadingjetptlow = 50, float subleadingjetpthigh = 100, int centmin = 0, int centmax = 40, float ajmin = 0.0, float ajmax = 1.0 , int doptweight = 0, int mccommand = 0, float jetamin = 0.0, float jetamax = 2.0, string whichjet = "")
{
  float pi78 = 7.0*TMath::Pi()/8.0;
  cout<<"runjettrk.C whichjet: "<<whichjet<<endl;
  string buffer;
  vector<string> listoffiles;
  int nlines = 0;
  ifstream infile(flist.data());
  if (!infile.is_open()) {
    cout << "Error opening file. Exiting." << endl;
    return;
  } else {
    while (!infile.eof()) {
      infile >> buffer;
      listoffiles.push_back(buffer);
      nlines++;
    }
  }
  jettrkana(listoffiles[condor_iter].data(),trackqual,whichjet);
  cout<<Form("%s_%s_trkqaul%d_nmin%d_nmax%d_tptmin%d_tptmax%d_aptmin%d_aptmax%d_%d.root",tag.data(),whichjet.data(),trackqual,nmin,nmax,(int)pttrigmin,(int)pttrigmax,(int)ptassmin,(int)ptassmax,condor_iter)<<endl;
  TFile * outf = new TFile(Form("%s_%s_trkqaul%d_nmin%d_nmax%d_cmin%d_cmax%d_aj%d_%d_tptmin%d_tptmax%d_aptmin%d_aptmax%d_%d.root",tag.data(),whichjet.data(),trackqual,nmin,nmax,centmin,centmax,(int)(ajmin*10),(int)(ajmax*10),(int)pttrigmin,(int)pttrigmax,(int)ptassmin,(int)ptassmax,condor_iter),"recreate");
  
	TH2D * ljtlead = JetTrackSignal    (0, leadingjetptlow, leadingjetpthigh , subleadingjetptlow , subleadingjetpthigh , ptassmin , ptassmax, centmin, centmax,ajmin,ajmax,doptweight,mccommand,jetamin,jetamax,15,pi78,whichjet);
	TH2D * bakjtlead = JetTrackBackground    (0, leadingjetptlow, leadingjetpthigh , subleadingjetptlow , subleadingjetpthigh , ptassmin , ptassmax, centmin, centmax,ajmin,ajmax,doptweight,mccommand,jetamin,jetamax,15,pi78,whichjet);
  
	TH2D * ljtsubl = JetTrackSignal    (1, leadingjetptlow, leadingjetpthigh , subleadingjetptlow , subleadingjetpthigh , ptassmin , ptassmax, centmin, centmax,ajmin,ajmax,doptweight,mccommand,jetamin,jetamax,15,pi78,whichjet);
	TH2D * bakjtsubl = JetTrackBackground    (1, leadingjetptlow, leadingjetpthigh , subleadingjetptlow , subleadingjetpthigh , ptassmin , ptassmax, centmin, centmax,ajmin,ajmax,doptweight,mccommand,jetamin,jetamax,15,pi78,whichjet);
  
  TH2D * ljtleadcorr = (TH2D*)ljtlead->Clone(Form("corr_%s_leadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptassmin,(int)ptassmax,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)));
  ljtleadcorr->Divide(bakjtlead);
  ljtleadcorr->Scale(bakjtlead->GetBinContent(bakjtlead->FindBin(0,0)));
  ljtleadcorr->GetXaxis()->SetRange(ljtleadcorr->GetXaxis()->FindBin(-4.0),ljtleadcorr->GetXaxis()->FindBin(4.0));
  ljtleadcorr->GetYaxis()->SetRange(ljtleadcorr->GetYaxis()->FindBin(-3.1415926/2.0),ljtleadcorr->GetYaxis()->FindBin(3*3.1415926/2.0));

	TH2D * ljtsublcorr = (TH2D*)ljtsubl->Clone(Form("corr_%s_subleadingjet%d_%d_ass%d_%d_cmin%d_cmax%d_ajmin%d_ajmax%d",whichjet.data(),(int)leadingjetptlow,(int)leadingjetpthigh,(int)ptassmin,(int)ptassmax,centmin,centmax,(int)(ajmin*100),(int)(ajmax*100)));
  ljtsublcorr->Divide(bakjtsubl);
  ljtsublcorr->Scale(bakjtsubl->GetBinContent(bakjtsubl->FindBin(0,0)));
  ljtsublcorr->GetXaxis()->SetRange(ljtsublcorr->GetXaxis()->FindBin(-4.0),ljtsublcorr->GetXaxis()->FindBin(4.0));
  ljtsublcorr->GetYaxis()->SetRange(ljtsublcorr->GetYaxis()->FindBin(-3.1415926/2.0),ljtsublcorr->GetYaxis()->FindBin(3*3.1415926/2.0));

	
	// TH2D * ttsig = TrackTrackSignal(pttrigmin,pttrigmax,ptassmin,ptassmax,nmin,nmax,nmin,nmax);
  // TH2D * ttbak = TrackTrackBackground(pttrigmin,pttrigmax,ptassmin,ptassmax,nmin,nmax,nmin,nmax);

  TH1I * hntottrig = new TH1I(Form("nttottrig_%s_trg%d_%d_ass%d_%d_nmin%d_nmax%d",whichjet.data(),(int)pttrigmin,(int)pttrigmax,(int)ptassmin,(int)ptassmax,nmin,nmax),"",1,0.5,1.5);
  int myntottrig = GetNTotTrig();
  hntottrig->Fill(1,myntottrig);
  cout<<"ntottrig: "<<myntottrig<<endl;

  // TH2D * ttcorr = (TH2D*)ttsig->Clone(Form("corr_trg%d_%d_ass%d_%d_nmin%d_nmax%d",(int)pttrigmin,(int)pttrigmax,(int)ptassmin,(int)ptassmax,nmin,nmax));
  // ttcorr->Divide(ttbak);
  // ttcorr->Scale(ttbak->GetBinContent(ttbak->FindBin(0,0)));
  // ttcorr->Scale(1/0.0594998609); //! bin width
  // ttcorr->GetXaxis()->SetRange(ttcorr->GetXaxis()->FindBin(-4.0),ttcorr->GetXaxis()->FindBin(4.0));
  // ttcorr->GetYaxis()->SetRange(ttcorr->GetYaxis()->FindBin(-3.1415926/2.0),ttcorr->GetYaxis()->FindBin(3*3.1415926/2.0));
  
  outf->Write();
  outf->Close();
}


int main(int argc, char *argv[])
{
  if(argc != 24)
  {
    std::cout << "Usage: runcorr <condor_iter> <trackqual> <file-list> <tag> <nmin> <nmax> <pttrigmin> <pttrigmax> <ptassmin> <ptassmax> <leadingjetptlow> <leadingjetpthigh> <subleadingjetptlow> <subleadingjetpthigh> <centmin> <centmax> <ajmin> <ajmax> <doptweight> <mccommand> <jetetamin> <jetetamax> <whichjet>" << std::endl;
    return 1;
  }
  runcorr(std::atoi(argv[1]), std::atoi(argv[2]), argv[3], argv[4], std::atoi(argv[5]), std::atoi(argv[6]), std::atof(argv[7]), std::atof(argv[8]), std::atof(argv[9]), std::atof(argv[10]), std::atof(argv[11]), std::atof(argv[12]), std::atof(argv[13]), std::atof(argv[14]), std::atoi(argv[15]), std::atoi(argv[16]), std::atof(argv[17]), std::atof(argv[18]), std::atoi(argv[19]), std::atoi(argv[20]), std::atof(argv[21]), std::atof(argv[22]), argv[23]);

  return 0;
}

