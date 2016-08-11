#include "multiplot.C"

void drawmultiplots(int ietaphi = 0)
{
  int aaptmin [3] = {60,73,86};
  int aaptmax [3] = {73,86,103};
  string tags [3] = {"pp_276TeV","pPb_5TeV","pp_7TeV"};
  string calos[3] = {"ak3PF","akPu3PF","ak3PF"};
  string detadphi [2] = {"eta","phi"};
  for(int i = 0 ; i < 3 ; ++i)
  {
    for(int j = 0 ; j < 3 ; ++j)
    {
      TCanvas * c1 = new TCanvas(Form("c1_%s_%d_%d",tags[i].data(),aaptmin[j],aaptmax[j]));
      multiplot("merged",tags[i],calos[i],140,160,aaptmin[j],aaptmax[j],ietaphi);
      c1->SaveAs(Form("d%s_fits_%s_%d_%d.png",detadphi[ietaphi].data(),tags[i].data(),aaptmin[j],aaptmax[j]));
    }
  }
}
