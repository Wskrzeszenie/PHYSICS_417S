
#include <iostream>
#include <fstream>
#include <vector>

Double_t Pol2(Double_t *x, Double_t *par)
//Polynomial test
{ 
  Double_t arg1= 0; 
  arg1=par[0]+ par[1]*x[0]+ par[2]*x[0]*x[0] ; 
  Double_t fitval = arg1 ; 
  /* 
  cout <<par[0]<<" "<<par[1]<<" "<<par[2]<<"\n"; 
  cout <<par[3]<<" "<<par[4]<<" "<<par[5]<<"\n"; 
  cout <<x[0]<<" "<<arg1<<" "<<arg2<<" "<<arg3<<" "<<fitval<<" "<<"\n"; 
  */ test
  return fitval; 
}

Double_t DLorentz(Double_t *x, Double_t *par) 
//Derivative of Lorentz
{
  Double_t amp=par[0];
  Double_t avg=par[1];
  Double_t halfs=(par[2]/2.)*(par[2]/2.);
  Double_t fitval=2.0*amp*halfs/( pow( (x[0]-avg)*(x[0]-avg)+halfs, 2.0));
  fitval=(x[0]-avg)*fitval;
  return fitval;
}

Double_t DGauss(Double_t *x, Double_t *par)  
//Gauss fitting 
{
    Double_t xnew=(x[0]-par[1]);
    Double_t fitval=par[0]*exp(-xnew*xnew/2.0/par[2]/par[2]);
    return fitval;
}

Double_t GaussPoly(Double_t *x, Double_t *par)  
//Gauss+Pol3 fitting 
{
    Double_t xnew=(x[0]-par[1]);
    Double_t fitval=par[0]*exp(-xnew*xnew/2.0/par[2]/par[2]);
    fitval=fitval+par[3]+par[4]*xnew+par[5]*xnew*xnew+par[6]*xnew*xnew*xnew;
    return fitval;
}

Double_t AplusSin(Double_t *x, Double_t *par)  
// A + sin(wt+phi) fitting 
{
    Double_t fitval=par[0] + par[1]*sin( par[2]*x[0]+par[3]);
    return fitval;
}


void MOSplotF(){
  // read in data
  vector<Double_t> vx;
  vector<Double_t> vy;
  Double_t xdat,ydat;

  // Read in Data
  fstream infile;
  infile.open("FE2O3.TKA", ios_base::in);
  while (infile>>ydat){
    //cout <<ydat<<"\n"; 
    vy.push_back(ydat) ;
  }
  vy[1]=0.0;
  vy[2]=0.0;
  vy[3]=0.0;
  infile.close();
  Int_t vsize = vy.size();
  //cout <<" Here I am "<<vsize<<"\n";

 
  // histogram parameters
  Int_t ntbin = 16384;       // number of total bins in the data
  Int_t nbin = 800;       // number of  bins yin the final plot
  double xmax=1.4;
  Double_t convert=0.00007;  // change channle to time

  // book histogram
  TH1F *hist1 = new TH1F("hist1","Raw Data",                    ntbin,0.0,16384.0);
  TH1F *hist2 = new TH1F("hist2","Rebinned  & Converted to T",  nbin,0.0,xmax);
  TH1F *hist3 = new TH1F("hist3","Background ",  nbin,0.0,xmax);
  TH1F *hist4 = new TH1F("hist4","Data - background ",  nbin,0.0,xmax);

  hist3 ->Sumw2();  //take care of error properly
  hist4 ->Sumw2();  //take care of error properly
  
  //plot histogram
  for (Int_t i=2; i!=16384; i++)
    { 
  Double_t xx=i+0.5;
  hist1->Fill(xx,vy[i]);
  hist1 ->SetBinError(i,sqrt(vy[i]));
  hist2->Fill(xx*convert,vy[i]);
    }

  for (Int_t i=1; i<nbin; i++)
    { 
    double yyy=hist2->GetBinContent(i);
    hist2 ->SetBinError(i,sqrt(yyy));
    }

  TCanvas *myc1 =new TCanvas("myc1","Sig1");
  //hist1->SetStats(kFALSE);
  hist1->GetYaxis()->SetTitle("number of counts");
  hist1->GetXaxis()->SetTitle(" channel ");
  hist1->Draw("*");

  TCanvas *myc2 =new TCanvas("myc2","Sig2");
  //hist2->SetStats(kFALSE);
  hist2->GetYaxis()->SetTitle("number of counts");
  hist2->GetXaxis()->SetTitle(" t (sec)");
  hist2->Draw("e");

  //fitting. Otaine w and phase angle for velocity convertion
  //Define parameters
  Double_t par[3];  //par[0]= Constant, par[1]=Amp, par[2]=omega, par[3]= phase angle

  //define fit
  TF1 *fit1 = new TF1("fit1",AplusSin, 0.0,2.0,4); // range & number of parameters 
  fit1 ->SetParameters(60000., 50000., 6.283, 1.5);  //initial parameter A+Bsin(wt+phi)
  //fit1 ->FixParameter(1, 0.0); 
  //fit1 ->FixParameter(2, 6.283);   //fix to 1 Hz 
  //fit1 ->FixParameter(3, 0.0); 
  Double_t xl1=0.015;   //fit range
  Double_t xh1=0.95;
  //Do the fit
  hist2  ->Fit("fit1","R"," ",xl1,xh1);
  hist2->Draw("same e");
  
  //Get  signal minus background
   Double_t a0 =fit1 ->GetParameter(0);
   Double_t a1 =fit1 ->GetParameter(1);
   Double_t a2 =fit1 ->GetParameter(2);
   Double_t a3 =fit1 ->GetParameter(3);
   
  for (Int_t i=0; i!=nbin; i++)
    { 
  Double_t xx=(i+0.5)*xmax/nbin;
  //Double_t background=par[0]+par[1]*sin(par[2]*xx+par[3]);
  Double_t background=a0+a1*sin(a2*xx+a3);
  hist3 ->Fill(xx,background*1.01);
  hist3->GetYaxis()->SetTitle("number of counts");
  hist3->GetXaxis()->SetTitle(" t (sec)");
  hist3 ->SetBinError(i,1.0);
    }
  TCanvas *myc3 =new TCanvas("myc3","sig3");
  hist3->Draw();

  TCanvas *myc4 =new TCanvas("myc4","Final");
  hist4 ->Add(hist2,hist3,-1.0,1.0);
  hist4->GetYaxis()->SetTitle("number of counts");
  hist4->GetXaxis()->SetTitle(" t (sec)");
  hist4->Draw("h");
  
}

