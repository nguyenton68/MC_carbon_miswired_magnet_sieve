/*
July 2017
Nguyen Ton
- Purpose: get cross section for sieve slit run from experimental data
  The cut is either a rectangular or a CUTG. You need to use graph cut, then save it.
  Additional cuts are: PID, number of track and W2 cut.
  Solid angle is determined from simulation.
  + Need to check:
  target thickness, density.
- Input: 2 sieve slit IN root files.
- Output: cross section in [ub]
 */

#include <TF1.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TCut.h>
const int nbin =200;
double hth_min(0.01),hth_max(0.065);
double hph_min(-0.03),hph_max(0.04);
void sieve_xs()
{
  //gStyle->SetOptStat(0);
  // Normalization
  double Heff_vdc,Heff_pid,Heff_sci;
  double Ceff_vdc,Ceff_pid,Ceff_sci;
  double Cps1,Hps1;
  double Cq,Hq;
  double Clt,Hlt;
  double runnum=1816;
  double P0,Ebeam;
  Heff_vdc = 1; Heff_pid = 1; Heff_sci = 1;
  Ceff_vdc = 1;  Ceff_pid = 1; Ceff_sci = 1;
  Clt = 1; Hlt = 1;
  if(runnum==1866)
    {
      Cps1 = 10; Hps1 = 3;
      Cq = 2041.0; Hq = 963.5;
      P0 = 1148.85;
      Ebeam = 1.14954;
    }
  if(runnum==1805)
    {
      Cps1 = 2; Hps1 = 2;
      Cq = 1927; Hq = 1485;
      P0 = 1148.85;
      Ebeam = 1.14887;
    }
  if(runnum==1809)
    {
      Cps1 = 1; Hps1 = 1;
      Cq = 924.3; Hq = 973.2;
      P0 = 1171.82;
      Ebeam = 1.14884;
    }
  if(runnum==1816)
    {
      Cps1 = 2; Hps1 = 2;
      Cq = 2317.7; Hq = 785;
      P0 = 1125.88;
      Ebeam = 1.14927;
    }
  if(runnum==1863)
    {
      Cps1 = 4; Hps1 = 1;
      Cq = 1796; Hq = 1172;
      P0 = 1183.29;
      Ebeam = 1.14961;
    }
  if(runnum==1869)
    {
      Cps1 = 11; Hps1 = 4;
      Cq = 2077; Hq = 1058;
      P0 = 1114.39;
      Ebeam = 1.14967;
    }
  double scl = (Hps1*Clt*Cq*Ceff_vdc*Ceff_pid*Ceff_sci)/(Cps1*Hlt*Hq*Heff_vdc*Heff_pid*Heff_sci);
  cout<<"Scale factor = "<<scl<<endl;


  TFile *fCarb = new TFile("gdh_1816.root"); 
  TFile *fHe4  = new TFile("gdh_1813.root");
  TTree *tdata=(TTree*)fCarb->Get("T");
  TTree *the=(TTree*)fHe4->Get("T");

  TH1F *hC = new TH1F("hC","",nbin,hth_min,hth_max);
  TH1F *hH = new TH1F("hH","",nbin,hth_min,hth_max);
  TH1F *hCyt = new TH1F("hCyt","",nbin,-0.03,0.03);
  TH1F *hHyt = new TH1F("hHyt","",nbin,-0.03,0.03);
  TH2F *hThPh = new TH2F("hThPh","",nbin,hph_min,hph_max,nbin,hth_min,hth_max);

  TString cutg ="./cut_1816/cut_1816_22f";
  gROOT->Macro(cutg);

  TCut cutPID =Form("CUTG && R.tr.n==1&&R.cer.asum_c>300&&(DR.evtypebits&(1<<1))>0&&(R.ps.e/(%f*(1+R.gold.dp)))>0.105&&((R.ps.e+R.sh.e)/(%f*(1+R.gold.dp)))>0.78",P0,P0);//P0 MeV

  TCut cutw2=Form("(sqrt(11.178**2+2*11.178*(%f-(%f*(1+R.gold.dp)))-4*%f*%f*(1+R.gold.dp)*sin(0.5*(TMath::ACos((0.9945+0.1045*R.gold.ph)/(1+R.gold.ph**2+R.gold.th**2))))**2)-11.178)<0.01",Ebeam,P0/1000.,Ebeam,P0/1000.);

  tdata->Draw("R.gold.th>>hC",cutPID&&cutw2);
  tdata->Draw("R.gold.y>>hCyt",cutPID&&cutw2);
  the->Draw("R.gold.th>>hH",cutPID&&cutw2);
  the->Draw("R.gold.y>>hHyt",cutPID&&cutw2);
  tdata->Draw("R.gold.th:R.gold.ph>>hThPh",cutPID&&cutw2);

  double nhe= hH->GetEntries();
  cout<<"Number of He4= "<<nhe<<endl;
  double nc = hC->GetEntries();

  double SA=2.08*1e-6;
  double scl_carb=Cps1*1.6*1e-19/(Ceff_vdc*Ceff_sci*Ceff_pid*Clt*Cq*1e-6);
  double scl_he  =Hps1*1.6*1e-19/(Heff_vdc*Heff_sci*Heff_pid*Hlt*Hq*1e-6);
  double dens = 1.933*0.0254*6.02*1e23;

  double xsection = (nc*scl_carb-nhe*scl_he)*1e30*12/(dens*SA);
  cout<<" cross section = [ub]"<<xsection<<endl;
  cout<<"Number of clean carbon= "<<(nc-scl*nhe)<<endl;
  TCanvas * cc0= new TCanvas("cc0","",800,600);
  cc0->Clear();
  hHyt->Scale(scl);
  hCyt->Add(hHyt,-1);
  hCyt->Draw();
  TCanvas * cc2 = new TCanvas("cc2","",800,600);
  cc2->Clear();
  hThPh->SetTitle("Single foil, dp = 0%;#phi_{tg}(rad);#theta_{tg}(rad)");
  hThPh->Draw("colz");
}
